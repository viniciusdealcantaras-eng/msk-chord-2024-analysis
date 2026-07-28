# Immortal time bias in survival reanalysis of MSK-CHORD — reproducibility code

Code to reproduce every figure and table of the manuscript
*"Viés de tempo imortal e truncamento à esquerda na reanálise de sobrevida a
partir de coortes de sequenciamento clínico"* (MSK-CHORD, cBioPortal).

The analysis shows that estimating overall survival **from the date of
diagnosis** in a cohort anchored at **first clinical sequencing** introduces
immortal time bias, and demonstrates the correction via left truncation
(entry time = diagnosis-to-sequencing interval).

## Data

Public cBioPortal study `msk_chord_2024`. The four raw files are downloaded
automatically by `build_dataset.py` (and, independently, by
`verificar_CHORD.R`) from the cBioPortal datahub, whose Git LFS content is
served at `media.githubusercontent.com`:

- `data_clinical_patient.txt`
- `data_clinical_sample.txt`
- `data_timeline_diagnosis.txt`
- `data_timeline_specimen.txt`

No raw data is redistributed in this repository. The data are the property of
Memorial Sloan Kettering Cancer Center / cBioPortal and are subject to their
terms of use; please cite the original MSK-CHORD study when using them.

The cohort is exactly **five cancer types**: non-small cell lung (NSCLC),
colorectal, breast, prostate, and pancreatic.

## Environment

Python >= 3.10 with the packages in `requirements.txt`
(`pip install -r requirements.txt`, i.e. `pandas numpy scipy matplotlib
lifelines pillow`). The independent check in R needs only base R (>= 4.0)
plus the `survival` package (already bundled with base R).

## Pipeline (run in order, from this folder)

Downloaded data land in `data/`; all outputs (CSVs, JSON results, and figure
PNGs) are written to `out/`. Both folders are created automatically and are
git-ignored.

1. `python build_dataset.py` — downloads the raw files (if missing), merges
   them, computes the diagnosis→sequencing interval (immortal time), and
   writes `out/merged_raw.csv`.
2. `python analysis.py` — delay distribution, immortal person-time, naive vs
   left-truncated vs from-sequencing Kaplan-Meier, spurious-group
   demonstration, and Cox models; writes `out/results.json`.
3. `python analysis2.py` — multiple-diagnosis frequency and single-diagnosis
   sensitivity analysis; writes `out/results2.json`.
4. `python analysis3.py` — Cox for stage IV **stratified by cancer type**
   (separate baseline hazard per histology); writes `out/results3_strat.json`.
5. `python per_type.py` — per-cancer-type KM (Figure 5) and the survival
   overestimation heatmap (Figure 6); writes `out/per_type.csv`,
   `out/fig5_km_by_type.png`, `out/fig6_heatmap_bias.png`.
6. `python figures.py` — Figures 1–4 (fig4 reads `out/results3_strat.json`);
   writes `out/fig1_immortal_time.png` … `out/fig4_cox_forest.png`.
7. `python verify.py` — internal check: a manual risk-set Kaplan-Meier with
   left truncation is compared against `lifelines` (medians and survival
   probabilities must match).

## Independent verification in R

`Rscript verificar_CHORD.R` downloads the data again and recomputes the key
numbers with the `survival` package (left truncation via
`Surv(entry, exit, event)`), printing a claimed-vs-obtained table.

## Notes

- **Time origin:** MSK-CHORD anchors OS at the date of first sequencing. The
  diagnosis-to-sequencing interval is used as the left-truncation entry time.
- Left-truncated KM/Cox use lifelines' `entry` / `entry_col` arguments.
- Patients with a diagnosis recorded *after* first sequencing (negative
  interval) are dropped as anomalous before analysis.
