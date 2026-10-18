# -*- coding: utf-8 -*-
"""Synthetic validation cohort v3 (ADEMP; Morris et al. 2019). Mechanism-based data-generating process at the
measurement level; NO Student-t (or any heavy-tailed law) is used to generate the data, so the estimator's own
distributional assumption cannot act as its goodness-of-fit criterion.
  Aims      : (1) recover the stable within-subject variation (CVI) from routine-like serial HbA1c with unknown
              clinical-change contamination; (2) return no spurious phenotype or covariate effect on dispersion when
              none exists (null) despite a realistic case-mix difference between phenotypes.
  DGP       : ln X_ij = ln S_i(t_ij) + e_bio + e_ana,  e_bio ~ N(0, sdln(CVI)),  e_ana ~ N(0, sdln(CVA)),
              sdln(CV) = sqrt(ln(1 + CV^2)); CVI = 1.6% (EFLM database median), CVA = 0.8% (desirable APS, 0.5 x CVI),
              identical in both phenotypes.  Treated set point S_i(t):
                baseline lognormal (uncomplicated median 7.4%, albuminuric DKD 8.1%; floor 6.6%),
                persistent treatment step at a follow-up visit with probability 0.15, size ~ N(-0.03, 0.08) ln-units,
                transient clinical episode with probability 0.05, amplitude ~ U(0.10, 0.30) ln-units, sign + with
                probability 0.6, decaying with time constant 60 days (erythrocyte kinetics),
                progression slope ~ N(0.01, 0.015) ln-units/year.
  Estimands : CVI (1.6), RCV = (exp(+/-1.96*sqrt(2)*sdln(CVTOT)) - 1)*100, gamma1 = 0, all beta = 0.
  Schema    : identical to the LIS workbook so the real pipeline runs unchanged (eligibility, exclusions, phenotype rule).
  usage: python3 -I make_sim_v3.py <out.pkl> [--cvi_unc 1.6 --cvi_neph 1.6 --cva 0.8 --p_step 0.15 --p_epi 0.05 --seed 20261007 --n_unc 1100 --n_neph 900 --noise gauss]"""
import numpy as np, pandas as pd, json, sys, argparse
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--cvi_unc', type=float, default=1.6); ap.add_argument('--cvi_neph', type=float, default=1.6)
ap.add_argument('--cva', type=float, default=0.8); ap.add_argument('--p_step', type=float, default=0.15); ap.add_argument('--p_epi', type=float, default=0.05)
ap.add_argument('--step_sd', type=float, default=0.08); ap.add_argument('--seed', type=int, default=20261007); ap.add_argument('--n_unc', type=int, default=1100); ap.add_argument('--n_neph', type=int, default=900)
ap.add_argument('--noise', default='gauss', help='gauss|laplace|lognormal|t4|mixture  (distribution of e_bio; same SD)')
ap.add_argument('--quiet', action='store_true')
A = ap.parse_args(); rng = np.random.default_rng(A.seed)
sdl = lambda cv: np.sqrt(np.log(1 + (cv / 100) ** 2))
def noise(sd, n):
    if A.noise == 'gauss': return rng.normal(0, sd, n)
    if A.noise == 'laplace': return rng.laplace(0, sd / np.sqrt(2), n)
    if A.noise == 'lognormal':                               # right-skewed, centred, unit-variance scaled
        s = 0.5; x = rng.lognormal(0, s, n); x = (x - np.exp(s**2/2)) / np.sqrt((np.exp(s**2) - 1) * np.exp(s**2)); return x * sd
    if A.noise == 't4': return rng.standard_t(4, n) * sd / np.sqrt(2)      # var of t4 = 2
    if A.noise == 'mixture':                                 # 90% N(0,0.8sd') + 10% N(0,2.2sd'), scaled to sd
        k = rng.random(n) < 0.1; s2 = 0.9 * 0.8**2 + 0.1 * 2.2**2; x = np.where(k, rng.normal(0, 2.2, n), rng.normal(0, 0.8, n)) / np.sqrt(s2); return x * sd
    raise ValueError(A.noise)
def ckd_epi_2021(scr, age, female):
    k, a = (0.7, -0.241) if female else (0.9, -0.302)
    return 142 * min(scr / k, 1) ** a * max(scr / k, 1) ** -1.200 * 0.9938 ** age * (1.012 if female else 1)
def lab_patient(acr_vals):
    v = [a for a in acr_vals if not np.isnan(a)]
    if not v: return 'No ACR (unassigned)', np.nan
    m = float(np.median(v))
    return ('None (ACR<30)' if m < 30 else 'Microalbuminuria (30-300)' if m <= 300 else 'Macroalbuminuria (>300)'), m
rows = []; pid = 0; nvis_p = {2: .15, 3: .25, 4: .25, 5: .17, 6: .10, 7: .05, 8: .03}
for neph, n, cvi in [(0, A.n_unc, A.cvi_unc), (1, A.n_neph, A.cvi_neph)]:
    for _ in range(n):
        pid += 1; female = rng.random() < 0.52
        age0 = int(np.clip(rng.normal(62 if neph else 60, 11), 30, 89))
        s0 = max(np.exp(rng.normal(np.log(8.1 if neph else 7.4), 0.15)), 6.6 + rng.random() * 0.4)
        scr_base = np.exp(rng.normal(np.log(1.15 if neph else 0.85), 0.25 if neph else 0.18))
        hb_base = rng.normal(13.2 if female else 14.6, 1.0) - (0.4 if neph else 0); mcv_base = rng.normal(86, 4.5)
        has_acr = rng.random() > (0.03 if neph else 0.06)
        acr_base = np.exp(rng.normal(np.log(110), 0.6)) if neph else np.exp(rng.normal(np.log(6), 0.75))
        slope = rng.normal(0.01, 0.015)
        nvis = rng.choice(list(nvis_p), p=list(nvis_p.values()))
        date0 = pd.Timestamp('2020-03-01') + pd.Timedelta(days=int(rng.integers(0, 1800)))
        t = 0.0; lnS = np.log(s0); episodes = []; steps = 0; epis = 0; acrs = []
        for v in range(1, nvis + 1):
            if v > 1:
                dt = float(np.clip(rng.gamma(5, 22), 21, 400)); t += dt
                if rng.random() < A.p_step: lnS += rng.normal(-0.03, A.step_sd); steps += 1
                if rng.random() < A.p_epi: episodes.append((t, rng.uniform(0.10, 0.30) * (1 if rng.random() < 0.6 else -1))); epis += 1
            tr = sum(amp * np.exp(-(t - t0) / 60) for t0, amp in episodes if t >= t0)
            x = np.exp(lnS + slope * t / 365.25 + tr + noise(sdl(cvi), 1)[0] + rng.normal(0, sdl(A.cva)))
            date = date0 + pd.Timedelta(days=int(round(t))); age = age0 + int(t // 365.25)
            scr = scr_base * np.exp(rng.normal(0, 0.06)); egfr = ckd_epi_2021(scr, age, female)
            if has_acr and rng.random() < (0.85 if neph else 0.7):
                acr = acr_base * np.exp(rng.normal(0, 0.35 if neph else 0.35))
            else: acr = np.nan
            acrs.append(acr)
            hb = hb_base + rng.normal(0, 0.4); mcv = mcv_base + rng.normal(0, 1.0)
            if rng.random() < (0.04 if neph else 0.012): hb = rng.uniform(8.0, 9.9)      # anaemia -> clinical exclusion
            if rng.random() < 0.01: mcv = rng.uniform(60, 69.9)
            if neph and rng.random() < 0.02: egfr = rng.uniform(8, 14.9)
            excl = '; '.join([s for s, c in [('eGFR<15 (advanced CKD)', egfr < 15), ('Hb<10 (anaemia)', hb < 10), ('MCV<70 (possible haemoglobinopathy)', mcv < 70)] if c]) or np.nan
            rows.append(dict(**{'Patient ID': f'S{pid:05d}', 'Measurement_No': v, 'Date': date.normalize(), 'HbA1c percent': round(x, 1),
                                'Gender': 'Female' if female else 'Male', 'Age': age, 'T2DM eligible': 'Yes',
                                'Glucose mg/dL': round((28.7 * x - 46.7) * np.exp(rng.normal(0, 0.12))), 'Creatinine serum mg/dL': round(scr, 2),
                                'eGFR': round(egfr, 2), 'Microalbumin mg/L': np.nan, 'Urine Creatinine mg/dL': np.nan, 'ACR mg/g': (round(acr, 1) if not np.isnan(acr) else np.nan),
                                'Hemoglobin': round(hb, 1), 'MCV_fL': round(mcv, 1), 'Nephropathy_measurement': '', 'ACR_median_patient': np.nan, 'Nephropathy_patient': '',
                                'Exclusion': excl, 'RCV_eligible': 'Yes', 'True_phenotype': 'neph' if neph else 'unc', 'True_steps': steps, 'True_episodes': epis, 'HbA1c_exact': x}))
df = pd.DataFrame(rows)
for p_, g in df.groupby('Patient ID'):
    lab, med = lab_patient(g['ACR mg/g'].values); df.loc[g.index, 'Nephropathy_patient'] = lab; df.loc[g.index, 'ACR_median_patient'] = med
df['Nephropathy_measurement'] = np.where(df['ACR mg/g'].isna(), 'No ACR (unassigned)', np.where(df['ACR mg/g'] < 30, 'None (ACR<30)', np.where(df['ACR mg/g'] <= 300, 'Microalbuminuria (30-300)', 'Macroalbuminuria (>300)')))
df['T2DM eligible'] = np.where(df.groupby('Patient ID')['HbA1c percent'].transform('max') >= 6.5, 'Yes', 'No')
df.to_pickle(A.out)
truth = dict(cvi_unc=A.cvi_unc, cvi_neph=A.cvi_neph, cva=A.cva, p_step=A.p_step, p_epi=A.p_epi, noise=A.noise, seed=A.seed,
             cvtot_unc=float(np.sqrt(A.cvi_unc**2 + A.cva**2)), cvtot_neph=float(np.sqrt(A.cvi_neph**2 + A.cva**2)))
for k in ['unc', 'neph']:
    sp = np.sqrt(2) * sdl(truth[f'cvtot_{k}']); truth[f'rcv_up_{k}'] = float((np.exp(1.96 * sp) - 1) * 100); truth[f'rcv_dn_{k}'] = float((np.exp(-1.96 * sp) - 1) * 100)
json.dump(truth, open(A.out.replace('.pkl', '_truth.json'), 'w'), indent=1)
if not A.quiet:
    print('patients', df['Patient ID'].nunique(), 'measurements', len(df), 'exclusions', df['Exclusion'].notna().sum())
    print(df.groupby('True_phenotype').agg(pts=('Patient ID', 'nunique'), meas=('Patient ID', 'size'), hba1c=('HbA1c percent', 'median'), egfr=('eGFR', 'median'), acr=('ACR mg/g', 'median'), steps=('True_steps', 'max')))
    print(pd.crosstab(df.drop_duplicates('Patient ID')['True_phenotype'], df.drop_duplicates('Patient ID')['Nephropathy_patient']))
    print('truth', truth)
