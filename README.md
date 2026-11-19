# Immortal time bias in survival reanalysis of MSK-CHORD: reproducibility code

Code to reproduce every number, figure and table of the manuscript
*"Reanchoring overall survival at diagnosis without accounting for delayed entry
overstates survival in clinical genomic cohorts: immortal time bias in MSK-CHORD"*.

The analysis shows that estimating overall survival **from the date of
diagnosis** in a cohort anchored at **first clinical sequencing** introduces
immortal time bias, and demonstrates the correction by left truncation
(entry time = diagnosis-to-sequencing interval).

## Data and license

Public cBioPortal study `msk_chord_2024` (Jee et al., Nature 2024). The four raw
files are downloaded automatically by `build_dataset.py` (and, independently, by
`verificar_CHORD.R`) from the cBioPortal datahub, whose Git LFS content is
served at `media.githubusercontent.com`:

- `data_clinical_patient.txt`
- `data_clinical_sample.txt`
- `data_timeline_diagnosis.txt`
- `data_timeline_specimen.txt`

MSK-CHORD is distributed by Memorial Sloan Kettering Cancer Center under the
[Creative Commons BY-NC-ND 4.0 license](https://creativecommons.org/licenses/by-nc-nd/4.0/).
For this reason **no raw or derived patient-level data are redistributed in this
repository**. `make_analytic_dataset.py` rebuilds the analytic dataset used in
the paper from the original files in one step. Please cite the original
MSK-CHORD study when using the data.

The cohort comprises exactly **five cancer types**: non-small cell lung,
colorectal, breast, prostate and pancreatic.

## Environment

- Python >= 3.10 with the packages in `requirements.txt`
  (`pip install -r requirements.txt`).
- R >= 4.0 with the `survival` package (bundled with base R) for the figures
  and the independent verification.

## Pipeline (run in order, from this folder)

Downloaded data land in `data/` and all outputs are written to `out/`. Both
folders are created automatically and are git-ignored.

1. `python build_dataset.py`: downloads the raw files (if missing) and writes
   `out/merged_raw.csv`.
2. `python make_analytic_dataset.py`: builds the patient-level analytic dataset
   `out/analytic_dataset.csv` (24,950 patients, 24,890 analysed), with the
   diagnosis-to-sequencing interval, survival from diagnosis and from
   sequencing, stage (highest recorded and at the index diagnosis), exclusion
   flags and the early/late sequencing groups.
3. `python analysis.py`: interval distribution, immortal person-time, naive vs
   left-truncated vs from-sequencing Kaplan-Meier, and the early/late
   demonstration. Writes `out/results.json`.
4. `python analysis2.py`: multiple-diagnosis frequency and single-diagnosis
   sensitivity analysis. Writes `out/results2.json`.
5. `python analysis3.py`: Cox models for stage IV stratified by cancer type.
   Writes `out/results3_strat.json`.
6. `python analysis4_landmark_tau.py`: landmark restricted to patients already
   sequenced by the landmark (reported for comparison in the note to
   Supplementary Table S3), subsample tau test and Schoenfeld test. Writes `out/results4_landmark_tau.json`.
7. `python analysis5_sensitivity.py`: sensitivity analyses: stage at the index
   diagnosis, exclusion of intervals longer than 10 and 20 years, early risk
   sets, medians by cancer type and the full-cohort conditional Kendall's tau.
   Writes `out/results5_sensitivity.json`.
8. `python analysis6_landmark_tau_by_type.py`: landmark at 6 and 12 months
   with delayed entry for patients sequenced after the landmark (Supplementary
   Table S3), conditional Kendall's tau with U-statistic standard error overall,
   by cancer type and by cancer type and stage (Supplementary Table S5), and
   characteristics, follow-up and stage IV hazard ratio by cancer type
   (Supplementary Table S4). Writes `out/results6_landmark_tau_by_type.json`.
9. `Rscript make_figures.R`: Figures 1 to 4 and Supplementary Figures S1 and
   S2, generated in R (base graphics and `survival`) from
   `out/analytic_dataset.csv`. Writes TIFF (600 dpi), PNG and PDF files to
   `out/figures/`.
10. `python verify.py`: internal check, a manual risk-set Kaplan-Meier with left
   truncation compared against `lifelines` (medians and survival
   probabilities must match).

## Independent verification in R

`Rscript verificar_CHORD.R` downloads the data again and recomputes the key
numbers with the `survival` package (left truncation via
`Surv(entry, exit, event)`), printing a claimed-vs-obtained table.

## Stata

`landmark_hr_by_type.do` reads `out/analytic_dataset.csv` and fits the
left-truncated Cox model, the stage IV hazard ratio by cancer type and the
landmark analysis with delayed entry. It does not save or overwrite any dataset.

## Notes

- **Time origin:** MSK-CHORD anchors OS at the date of first sequencing. The
  diagnosis-to-sequencing interval is used as the left-truncation entry time.
- Left-truncated KM and Cox models use `entry` / `entry_col` in lifelines and
  `Surv(entry, exit, event)` in R.
- Patients with no diagnosis in the timeline (n = 10) or with a diagnosis
  recorded after first sequencing (n = 50) are excluded.
- Entry times are kept at full numerical precision.
