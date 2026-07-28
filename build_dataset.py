import os, urllib.request
import pandas as pd, numpy as np

DATA = "data/"
OUT = "out/"
MO = 30.4375  # days per month

# --- ensure folders exist ---
os.makedirs(DATA, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

# --- download raw files from the cBioPortal datahub if missing ---
BASE = ("https://media.githubusercontent.com/media/cBioPortal/datahub/"
        "master/public/msk_chord_2024")
RAW_FILES = ["data_clinical_patient.txt", "data_clinical_sample.txt",
             "data_timeline_diagnosis.txt", "data_timeline_specimen.txt"]
for f in RAW_FILES:
    dest = DATA + f
    if not os.path.exists(dest):
        print("downloading", f, "...")
        urllib.request.urlretrieve(BASE + "/" + f, dest)

def read_clinical(path):
    # 4 comment lines (#...) then real header line, then data
    return pd.read_csv(path, sep="\t", skiprows=4)

pat = read_clinical(DATA+"data_clinical_patient.txt")
samp = read_clinical(DATA+"data_clinical_sample.txt")
dx = pd.read_csv(DATA+"data_timeline_diagnosis.txt", sep="\t", dtype=str)
spec = pd.read_csv(DATA+"data_timeline_specimen.txt", sep="\t", dtype=str)

print("patients:", len(pat), "| samples:", len(samp))
print("pat cols:", list(pat.columns))
print("dx rows:", len(dx), "| spec rows:", len(spec))

# --- specimen: first sequencing per patient; verify day 0 anchor ---
spec["START_DATE"] = pd.to_numeric(spec["START_DATE"], errors="coerce")
spec_min = spec.groupby("PATIENT_ID")["START_DATE"].min()
print("\n[CHECK] specimen min START_DATE per patient: describe")
print(spec_min.describe())
print("fraction of patients whose first specimen == day 0:",
      float((spec_min == 0).mean()))

# --- diagnosis: earliest diagnosis per patient (most negative offset) ---
dx["START_DATE"] = pd.to_numeric(dx["START_DATE"], errors="coerce")
dx_min = dx.groupby("PATIENT_ID")["START_DATE"].min().rename("dx_offset_days")
print("\n[CHECK] earliest diagnosis START_DATE (days, should be <=0):")
print(dx_min.describe())
print("fraction diagnosis offset <= 0:", float((dx_min <= 0).mean()))
print("fraction diagnosis offset > 0 (dx AFTER seq, anomalous):", float((dx_min > 0).mean()))

# delay diagnosis -> first sequencing (immortal time), in months
delay_days = (-dx_min).clip(lower=None)  # positive = days from dx to seq
delay = (delay_days / MO).rename("delay_months")

# --- OS from patient file (anchored at first sequencing) ---
pat["OS_MONTHS"] = pd.to_numeric(pat["OS_MONTHS"], errors="coerce")
pat["deceased"] = pat["OS_STATUS"].astype(str).str.startswith("1").astype(int)
os_raw_na = pat["OS_MONTHS"].isna().sum()
os_zero = ((pat["OS_MONTHS"] == 0)).sum()
os_zero_dead = ((pat["OS_MONTHS"] == 0) & (pat["deceased"] == 1)).sum()
print(f"\n[OS] NA OS_MONTHS: {os_raw_na} ({os_raw_na/len(pat)*100:.1f}%) | "
      f"OS==0: {os_zero} | OS==0 & deceased: {os_zero_dead}")

# --- cancer type (first sample per patient) ---
ctype = (samp.sort_values("SAMPLE_ID")
         .groupby("PATIENT_ID")["CANCER_TYPE"].first().rename("cancer_type"))

# --- assemble ---
df = pat[["PATIENT_ID","OS_MONTHS","deceased","OS_STATUS",
          "STAGE_HIGHEST_RECORDED","CURRENT_AGE_DEID","GENDER"]].copy()
df = df.merge(delay, left_on="PATIENT_ID", right_index=True, how="left")
df = df.merge(ctype, left_on="PATIENT_ID", right_index=True, how="left")
df["CURRENT_AGE_DEID"] = pd.to_numeric(df["CURRENT_AGE_DEID"], errors="coerce")

df.to_csv("out/merged_raw.csv", index=False)
print("\nmerged rows:", len(df))
print("with delay:", df["delay_months"].notna().sum(),
      "| with OS:", df["OS_MONTHS"].notna().sum(),
      "| with both:", (df["delay_months"].notna() & df["OS_MONTHS"].notna()).sum())
print("\ndelay_months describe (all):")
print(df["delay_months"].describe(percentiles=[.25,.5,.75,.9]))
print("\nstage counts:\n", df["STAGE_HIGHEST_RECORDED"].value_counts(dropna=False))
print("\ntop cancer types:\n", df["cancer_type"].value_counts().head(12))
