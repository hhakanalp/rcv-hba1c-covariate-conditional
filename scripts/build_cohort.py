# -*- coding: utf-8 -*-
"""Pipeline D cohort builder: identical to pipeline C (T2DM eligibility, clinical exclusions, >=2 visits, phenotype
from Nephropathy_patient with >=2 ACR >= 30 mg/g for confirmed albuminuria) + restriction of consecutive pairs to an
inter-measurement interval of 1-180 days. Works on the real LIS workbook and on synthetic workbooks with the same schema.
usage: python3 -I build_cohort.py <input.pkl|xlsx> <out_pairs.pkl> [cohort_json]"""
import pandas as pd, numpy as np, sys, json
src, out = sys.argv[1], sys.argv[2]
d = pd.read_excel(src) if src.endswith('.xlsx') else pd.read_pickle(src)
d = d.copy(); d['Date'] = pd.to_datetime(d['Date'])
d.loc[d['HbA1c percent'] == 56.0, 'HbA1c percent'] = 5.6          # transcription error (confirmed: 5.4 / 5.6 patient)
flow = {'records': int(len(d)), 'patients': int(d['Patient ID'].nunique())}
elig = d.groupby('Patient ID')['HbA1c percent'].transform('max') >= 6.5
e = d[elig]; flow['t2dm_meas'] = int(len(e)); flow['t2dm_pts'] = int(e['Patient ID'].nunique())
r = e[e['Exclusion'].isna()].copy(); flow['excluded_meas'] = int(len(e) - len(r)); flow['after_excl_meas'] = int(len(r)); flow['after_excl_pts'] = int(r['Patient ID'].nunique())
r = r[r.groupby('Patient ID')['Patient ID'].transform('size') >= 2].sort_values(['Patient ID', 'Date'])
flow['ge2_meas'] = int(len(r)); flow['ge2_pts'] = int(r['Patient ID'].nunique())
# phenotype: patient label from workbook; nephropathy must be confirmed by >=2 ACR >= 30 mg/g over all visits
n30 = d.groupby('Patient ID')['ACR mg/g'].apply(lambda a: (a >= 30).sum())
lab = r.groupby('Patient ID')['Nephropathy_patient'].first()
ph = pd.Series(np.where(lab.str.startswith('None'), 'unc', np.where(lab.str.startswith('No ACR'), 'noACR', np.where(n30.reindex(lab.index).fillna(0) >= 2, 'neph', 'unconf'))), index=lab.index)
g = r.groupby('Patient ID'); r['X2'] = g['HbA1c percent'].shift(-1); r['D2'] = g['Date'].shift(-1)
p = r.dropna(subset=['X2']).copy(); p['pheno'] = p['Patient ID'].map(ph)
p['delta'] = np.log(p['X2'] / p['HbA1c percent']); p['days'] = (p['D2'] - p['Date']).dt.days
flow['pairs_all'] = int(len(p)); flow['pairs_all_pts'] = int(p['Patient ID'].nunique())
for k in ['noACR', 'unconf']: flow[f'{k}_pairs'] = int((p.pheno == k).sum()); flow[f'{k}_pts'] = int(p[p.pheno == k]['Patient ID'].nunique())
p = p[p.pheno.isin(['unc', 'neph'])].copy()
flow['phenotyped_pairs'] = int(len(p)); flow['phenotyped_pts'] = int(p['Patient ID'].nunique())
flow['phenotyped_meas'] = int(r[r['Patient ID'].isin(p['Patient ID'])].shape[0])
for k in ['unc', 'neph']: flow[f'{k}_pairs_nolimit'] = int((p.pheno == k).sum()); flow[f'{k}_pts_nolimit'] = int(p[p.pheno == k]['Patient ID'].nunique())
p180 = p[(p.days >= 1) & (p.days <= 180)].copy()
flow['pairs_gt180'] = int((p.days > 180).sum()); flow['pairs_lt1'] = int((p.days < 1).sum())
flow['D_pairs'] = int(len(p180)); flow['D_pts'] = int(p180['Patient ID'].nunique())
for k in ['unc', 'neph']: flow[f'{k}_pairs'] = int((p180.pheno == k).sum()); flow[f'{k}_pts'] = int(p180[p180.pheno == k]['Patient ID'].nunique())
flow['D_meas'] = int(len(set(zip(p180['Patient ID'], p180['Date'])) | set(zip(p180['Patient ID'], p180['D2']))))
p180.to_pickle(out); p.to_pickle(out.replace('.pkl', '_nolimit.pkl'))
if len(sys.argv) > 3: json.dump(flow, open(sys.argv[3], 'w'), indent=1)
print(json.dumps(flow, indent=1))
