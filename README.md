# Covariate-conditional reference change values for HbA1c from routine laboratory data

Analysis code for:

> Akbay HI, Çubukcu HC, Theodorsson E, Alp HH. *Covariate-conditional reference change values for HbA1c derived from routine laboratory data: type 2 diabetes with and without diabetic nephropathy.* Submitted to Clinical Chemistry, 2026.

Zenodo DOI: [10.5281/zenodo.23449051](https://doi.org/10.5281/zenodo.23449051)

## What is here

| Folder | Content |
|---|---|
| `data/` | `Synthetic_HbA1c_v3_CVI1.6_CVA0.8.xlsx` — synthetic validation cohort (no patient data), usable as a runnable example. `INPUT_SCHEMA.md` — the column layout your own LIS export must have. `DATA_DICTIONARY.md` — how the synthetic cohort was generated. |
| `scripts/` | Python 3.11 analysis pipeline (MIT licence). |
| `results/` | Frozen outputs used for the submitted manuscript (`resultsD.json` = all numbers in Tables 2–3 and the Results text; `ademp_grid.json` = Supplemental Table 9; `coverage_D.json`; flow-chart counts; `SIM_v3_truth.json` = injected synthetic truth). |
| `run_all.sh` | Runs the whole pipeline in order. |

Patient data are not shared. The repository contains the complete analysis code and a synthetic cohort so that the pipeline can be run end to end and applied to other laboratories' own serial HbA1c data (see `data/INPUT_SCHEMA.md`).

## Pipeline

```
build_cohort.py             LIS-schema workbook -> consecutive pairs (1-180 d and all) + flow-chart counts; run it on your
                            own export (data/INPUT_SCHEMA.md) and on the synthetic workbook
build_D.py                  Student-t central scale per phenotype, CVTOT, CVI, asymmetric RCV, patient-cluster bootstrap
                            (B = 1000), Gaussian-mixture and truncation benchmarks (Table 2, Supplemental Tables 2 and 7)
diff_D.py                   bootstrap CIs of the between-phenotype differences
fit_hier_D.py / run_hierD.sh hierarchical Student-t dispersion model in PyMC (NUTS, 4 chains x 1000/1000):
                            full model, sex-stratified, without baseline HbA1c, HbA1c only, pair-mean HbA1c
                            (Table 3, Supplemental Tables 3 and 8)
make_sim_v3.py              synthetic cohort generator (mechanism-based, Gaussian noise; seed 20261007)
ademp_grid.py               ADEMP simulation grid, 10 replicates per scenario (Supplemental Table 9)
coverage_D.py               stable-pair coverage of +/-1.96 sigma_t in the synthetic cohort
collect_D.py                gathers all outputs into work/resultsD.json
make_main_figs_D.py, make_supp_figs_D.py, make_flow_D.py   Figures 1-3 and Supplemental Figures 1-5 (600 dpi TIFF)
tfit.py                     ECME maximum-likelihood fit of the location-scale Student-t (nu, mu, sigma)
six_month.py, table1_D.py, supp_numbers_D.py   descriptive numbers (Table 1, interval strata) from the measurement-level export
```

## Running

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
bash run_all.sh                      # synthetic cohort only (demo); edit the first lines to point to your own export
```

`run_all.sh` builds `work/pairs_D.pkl` from the workbook named in its first line; place your own export there (or at
`work/real.xlsx` as the script expects) to reproduce the full real-data analysis. The synthetic cohort reproduces the
synthetic rows of Tables 2–3 and Supplemental Table 9 of the article (CVI,Phenotype 1.74% / 1.55% against an injected
1.6%; γ1 compatible with zero). The bootstrap (B = 1000 with refitting) takes minutes; the PyMC models take 1–3 hours in
total on a laptop. Results reproduce to the reported precision with the seeds in the scripts; small differences in the
third decimal may arise across PyMC/NumPy versions.

Analytical variation (CVA = 1.1%, root mean square of 0.9% and 1.2% at two IQC levels over 12 months) is entered as a
constant in `build_D.py`; the daily IQC records are not part of this repository.

## Ethics

Non-Interventional Clinical Research Ethics Committee of Van Yüzüncü Yıl University, decision no. 2026/07-39
(14 August 2026); the requirement for individual informed consent was waived because only anonymized routine
laboratory data were used.

## Licences and citation

Code and synthetic data: MIT (`LICENSE`). Please cite the Zenodo record (`CITATION.cff`) and the article.
