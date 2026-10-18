# -*- coding: utf-8 -*-
"""Coverage of the +/-1.96*sigma_t band among synthetic pairs with no new treatment step or clinical episode between the two
measurements (True_steps / True_episodes are cumulative counters per measurement) and among all pairs; also the P(|T|<1.96) of the
fitted Student-t. Writes work/coverage_D.json.  usage: python3 -I coverage_D.py"""
import pandas as pd, numpy as np, json
from scipy import stats
p = pd.read_pickle('work/pairs_D_sim.pkl'); s = pd.read_pickle('work/SIM_v3.pkl')
s = s.sort_values(['Patient ID', 'Measurement_No']); T2 = pd.read_pickle('work/table2_D.pkl')
key = s.set_index(['Patient ID', 'Measurement_No'])[['True_steps', 'True_episodes']]
# pair rows carry the first measurement; D2 is the date of the second; find the second measurement's counters
nxt = s.groupby('Patient ID')[['True_steps', 'True_episodes']].shift(-1); s['steps2'] = nxt['True_steps']; s['epis2'] = nxt['True_episodes']
m = p.merge(s[['Patient ID', 'Measurement_No', 'steps2', 'epis2']], on=['Patient ID', 'Measurement_No'], how='left')
m['stable'] = (m.steps2 == m.True_steps) & (m.epis2 == m.True_episodes)
out = {}
for ph in ['unc', 'neph']:
    r = T2[('sim', ph)]; sig = r['sig_t']; nu = r['nu']; mu = r.get('mu', 0.0)
    q = m[m.pheno == ph]; st = q[q.stable]
    cov_all = float((np.abs(q.delta - mu) <= 1.96 * sig).mean()); cov_st = float((np.abs(st.delta - mu) <= 1.96 * sig).mean())
    out[ph] = dict(n=int(len(q)), n_stable=int(len(st)), frac_stable=float(len(st) / len(q)), cov_all=cov_all * 100, cov_stable=cov_st * 100,
                   sd_stable_cv=float(st.delta.std(ddof=1) / np.sqrt(2) * 100), sigt_cv=float(sig / np.sqrt(2) * 100), p_t_1_96=float(stats.t.cdf(1.96, nu) - stats.t.cdf(-1.96, nu)) * 100, nu=float(nu))
for k in [('real', 'unc'), ('real', 'neph')]:
    r = T2[k]; out['|'.join(k)] = dict(p_t_1_96=float(stats.t.cdf(1.96, r['nu']) - stats.t.cdf(-1.96, r['nu'])) * 100, nu=float(r['nu']))
json.dump(out, open('work/coverage_D.json', 'w'), indent=1); print(json.dumps(out, indent=1))
