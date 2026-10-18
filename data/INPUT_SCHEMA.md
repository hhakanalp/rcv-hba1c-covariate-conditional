# Input data schema

Patient data are not distributed with this repository. To run the pipeline on your own laboratory data, export one
workbook (.xlsx) or pickle with one row per HbA1c measurement and the columns below (names must match exactly);
`scripts/build_cohort.py` then builds the consecutive within-subject pairs and applies the eligibility, exclusion and
phenotype rules of the article.

| Column | Type | Description |
|---|---|---|
| Patient ID | string | Any patient identifier (used only for grouping; never written to outputs) |
| Date | date | Sampling date |
| HbA1c percent | float | HbA1c, NGSP % |
| Gender | string | `Female` / `Male` |
| Age | int | Age at sampling, years |
| eGFR | float | CKD-EPI 2021 creatinine eGFR, mL/min/1.73 m² (may be missing) |
| ACR mg/g | float | Urinary albumin-to-creatinine ratio, mg/g (missing when no urine sample) |
| Hemoglobin | float | g/dL (may be missing) |
| MCV_fL | float | fL (may be missing) |
| Nephropathy_patient | string | Patient-level label from your own albuminuria rule: must start with `None` (no confirmed albuminuria), `No ACR` (no ACR result at all) or anything else for candidate nephropathy; `build_cohort.py` confirms nephropathy only when ≥ 2 visits have ACR ≥ 30 mg/g |
| Exclusion | string | Empty/NaN for interpretable results; any text flags the measurement as uninterpretable (hemoglobin < 10 g/dL, MCV < 70 fL, eGFR < 15 mL/min/1.73 m² in the article) |

T2DM eligibility (≥ 1 HbA1c ≥ 6.5% during follow-up), the ≥ 2-measurement rule, pairing of consecutive results and the
1–180-day interval restriction are applied inside `build_cohort.py`. The synthetic workbook in this folder uses exactly
this schema and serves as a runnable example; `data/DATA_DICTIONARY.md` describes how it was generated.
Analytical variation (CVA) is a constant in `scripts/build_D.py` (1.1% for the article's HPLC method; 0.8% for the
synthetic cohort) — replace it with your own long-term IQC estimate.
