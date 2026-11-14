"""Build the patient-level analytic dataset from the raw MSK-CHORD files.

Input : data/data_clinical_patient.txt, data/data_clinical_sample.txt,
        data/data_timeline_diagnosis.txt (cBioPortal datahub, msk_chord_2024)
Output: out/analytic_dataset.csv (one row per patient, n = 24,950)

Run `python build_dataset.py` first, or place the four raw files in data/.
"""
import os
import numpy as np
import pandas as pd

DATA, OUT, MO = "data/", "out/", 30.4375
os.makedirs(OUT, exist_ok=True)

pat = pd.read_csv(DATA + "data_clinical_patient.txt", sep="\t", skiprows=4)
samp = pd.read_csv(DATA + "data_clinical_sample.txt", sep="\t", skiprows=4)
dx = pd.read_csv(DATA + "data_timeline_diagnosis.txt", sep="\t", dtype=str)
dx["START_DATE"] = pd.to_numeric(dx["START_DATE"], errors="coerce")

# index diagnosis = earliest diagnosis event (most negative offset to first sequencing)
first = dx.sort_values(["PATIENT_ID", "START_DATE"]).groupby("PATIENT_ID").first()
ndx = dx.groupby("PATIENT_ID").size()
ctype = samp.sort_values("SAMPLE_ID").groupby("PATIENT_ID")["CANCER_TYPE"].first()

stage_map = {"Stage 1-3": "I-III", "Stage 4": "IV"}
f = pd.DataFrame({"patient_id": pat["PATIENT_ID"]})
f["cancer_type"] = f.patient_id.map(ctype)
f["stage_highest_recorded"] = pat["STAGE_HIGHEST_RECORDED"].map(stage_map)
f["stage_at_index_diagnosis"] = f.patient_id.map(first["STAGE_CDM_DERIVED"]).map(stage_map)
f["sex"] = pat["GENDER"]
f["current_age_deid"] = pd.to_numeric(pat["CURRENT_AGE_DEID"], errors="coerce")
f["os_status"] = pat["OS_STATUS"].astype(str).str.startswith("1").astype(int)
f["os_months_from_sequencing"] = pd.to_numeric(pat["OS_MONTHS"], errors="coerce")
f["dx_offset_days"] = f.patient_id.map(first["START_DATE"])
f["n_diagnosis_events"] = f.patient_id.map(ndx).fillna(0).astype(int)

f["exclusion_reason"] = np.select(
    [f.dx_offset_days.isna(), f.dx_offset_days > 0],
    ["No diagnosis event in timeline", "Diagnosis recorded after sequencing"], "")
f["included"] = (f.exclusion_reason == "").astype(int)
inc = f.included == 1
f["dx_to_seq_months"] = np.where(inc, -f.dx_offset_days / MO, np.nan)
f["os_months_from_diagnosis"] = f.dx_to_seq_months + f.os_months_from_sequencing
med = f.loc[inc, "dx_to_seq_months"].median()
f["seq_timing_group"] = ""
f.loc[inc, "seq_timing_group"] = np.where(f.loc[inc, "dx_to_seq_months"] >= med, "Late", "Early")

f.to_csv(OUT + "analytic_dataset.csv", index=False)
print(f"patients {len(f)} | included {int(inc.sum())} | deaths {int(f.loc[inc, 'os_status'].sum())}"
      f" | median interval {med:.2f} months")
