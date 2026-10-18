# -*- coding: utf-8 -*-
"""Patient-cluster bootstrap CI for the nephropathy-minus-uncomplicated difference in CVI and RCV limits (pipeline D)."""
import pandas as pd, numpy as np, json, sys
sys.path.insert(0, 'scripts')
from tfit import t_fit
q = pd.read_pickle('work/pairs_D.pkl'); sp = pd.read_pickle('work/pairs_D_sim.pkl')
CVA2 = {'real': 1.1 ** 2, 'sim': 0.8 ** 2}
def est(x, cva2):
    nu, mu, sg = t_fit(x); cvt = sg / np.sqrt(2) * 100
    return np.array([np.sqrt(max(cvt ** 2 - cva2, 0)), (np.exp(1.96 * sg) - 1) * 100, (np.exp(-1.96 * sg) - 1) * 100, sg])
def prep(s):
    x = s['delta'].values; ids, inv = np.unique(s['Patient ID'].values, return_inverse=True)
    return x, [np.where(inv == i)[0] for i in range(len(ids))]
out = {}
rng = np.random.default_rng(7)
for key, df in [('real', q), ('sim', sp)]:
    for sx in [None, 'Female', 'Male']:
        d = df if sx is None else df[df.Gender == sx]
        xu, gu = prep(d[d.pheno == 'unc']); xn, gn = prep(d[d.pheno == 'neph'])
        pt = est(xn, CVA2[key]) - est(xu, CVA2[key])
        Bs = []
        for b in range(1000):
            iu = np.concatenate([gu[i] for i in rng.integers(0, len(gu), len(gu))]); inn = np.concatenate([gn[i] for i in rng.integers(0, len(gn), len(gn))])
            Bs.append(est(xn[inn], CVA2[key]) - est(xu[iu], CVA2[key]))
        Bs = np.array(Bs); lo, hi = np.percentile(Bs, 2.5, axis=0), np.percentile(Bs, 97.5, axis=0)
        ratio = np.percentile(np.array([1] * 0 + list(Bs[:, 3])), [2.5, 97.5])
        tag = f"{key}_{sx or 'all'}"
        out[tag] = dict(dcvi=(pt[0], lo[0], hi[0]), dup=(pt[1], lo[1], hi[1]), ddn=(pt[2], lo[2], hi[2]), dsig=(pt[3], lo[3], hi[3]), p_le0=float((Bs[:, 0] <= 0).mean()))
        print(tag, 'ΔCVI %.2f (%.2f to %.2f)  Δup %.1f (%.1f to %.1f)  Δdn %.1f (%.1f to %.1f)  P(Δ<=0)=%.3f' % (pt[0], lo[0], hi[0], pt[1], lo[1], hi[1], pt[2], lo[2], hi[2], out[tag]['p_le0']), flush=True)
json.dump(out, open('work/diff_D.json', 'w'), indent=1)
