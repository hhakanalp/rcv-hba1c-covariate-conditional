# -*- coding: utf-8 -*-
"""Pipeline D: pipeline C (no |z|<=2 truncation) on consecutive pairs 1-180 days apart. Central scale = ML scale of a location-scale Student-t fitted to ALL
consecutive log-differences of each stratum (robust analogue of Jones's central Gaussian). CVTOT = sigma_t/sqrt2,
CVI = sqrt(CVTOT^2 - CVA^2), RCV+/- = (exp(+/-1.96 sigma_t) - 1)*100. Patient-cluster bootstrap CIs (B=1000).
Benchmarks: 2-component Gaussian mixture central component; classical symmetric RCV."""
import pandas as pd, numpy as np, pickle, json, sys
sys.path.insert(0, 'scripts')
from tfit import t_fit
from sklearn.mixture import GaussianMixture

q = pd.read_pickle('work/pairs_D.pkl'); sp = pd.read_pickle('work/pairs_D_sim.pkl')
CVA2 = {'real': 1.1 ** 2, 'sim': 0.8 ** 2}
B = 1000

def est(x, cva2):
    nu, mu, sg = t_fit(x)
    cvt = sg / np.sqrt(2) * 100
    cvi = np.sqrt(max(cvt ** 2 - cva2, 0))
    return cvi, (np.exp(1.96 * sg) - 1) * 100, (np.exp(-1.96 * sg) - 1) * 100, cvt, sg, nu, mu

def gmm_central(x, seed=0):
    g = GaussianMixture(2, n_init=3, random_state=seed).fit(x.reshape(-1, 1))
    k = int(np.argmin(g.covariances_.ravel()))
    return float(np.sqrt(g.covariances_.ravel()[k])), float(g.weights_[k])

def jones_iter(x, z=1.96, corr=0.8711, tol=1e-6, maxit=200):
    # iterative symmetric truncation: keep |x-median| <= z*sd_corrected until stable; SD corrected for truncation (0.8711 at +/-1.96)
    keep = np.ones(len(x), bool)
    for _ in range(maxit):
        m = np.median(x[keep]); sd = x[keep].std(ddof=1) / corr; new = np.abs(x - m) <= z * sd
        if new.sum() == keep.sum() and np.all(new == keep): break
        keep = new
    return float(sd), float(keep.mean())

def block(s, cva2, seed=42):
    rng = np.random.default_rng(seed)
    x = s['delta'].values; pid = s['Patient ID'].values
    pt = est(x, cva2)
    ids, inv = np.unique(pid, return_inverse=True); groups = [np.where(inv == i)[0] for i in range(len(ids))]
    Bs = []
    for b in range(B):
        idx = np.concatenate([groups[i] for i in rng.integers(0, len(ids), len(ids))])
        Bs.append(est(x[idx], cva2)[:3])
    Bs = np.array(Bs); lo, hi = np.percentile(Bs, 2.5, axis=0), np.percentile(Bs, 97.5, axis=0)
    nu, mu, sg = pt[5], pt[6], pt[4]
    w = (nu + 1) / (nu + ((x - mu) / sg) ** 2)
    gs, gw = gmm_central(x); js, jk = jones_iter(x)
    return dict(n=len(s), n_pat=int(len(ids)), nu=nu, mu=mu, sig_t=sg, sd=sg, cvtot=pt[3], cvi=pt[0], cvi_ci=(lo[0], hi[0]),
                up=pt[1], up_ci=(lo[1], hi[1]), dn=pt[2], dn_ci=(lo[2], hi[2]),
                classical=1.96 * np.sqrt(2) * np.sqrt(pt[0] ** 2 + cva2),
                w_lt05=float((w < 0.5).mean() * 100), raw_sd=float(x.std(ddof=1)),
                gmm_sig=gs, gmm_w=gw, gmm_cvtot=gs / np.sqrt(2) * 100, gmm_cvi=float(np.sqrt(max((gs / np.sqrt(2) * 100) ** 2 - cva2, 0))),
                jones_sig=js, jones_keep=jk, jones_cvi=float(np.sqrt(max((js / np.sqrt(2) * 100) ** 2 - cva2, 0))),
                median_days=float(s['days'].median()), iqr_days=(float(s['days'].quantile(.25)), float(s['days'].quantile(.75))))

def show(tag, r):
    print(f"{tag:28s} n={r['n']:6d} pts={r['n_pat']:5d} nu={r['nu']:.2f} sig_t={r['sig_t']:.4f} CVTOT={r['cvtot']:.2f} CVI={r['cvi']:.2f} ({r['cvi_ci'][0]:.2f}-{r['cvi_ci'][1]:.2f}) "
          f"RCV {r['up']:+.1f} ({r['up_ci'][0]:.1f}-{r['up_ci'][1]:.1f}) / {r['dn']:.1f} ({r['dn_ci'][0]:.1f} to {r['dn_ci'][1]:.1f}) classical {r['classical']:.1f} | w<0.5 {r['w_lt05']:.1f}% rawSD {r['raw_sd']:.3f} | GMM2 CVI {r['gmm_cvi']:.2f} (w {r['gmm_w']:.2f})", flush=True)

T2 = {}
for key, df in [('real', q), ('sim', sp)]:
    for ph in ['All', 'unc', 'neph']:
        s = df if ph == 'All' else df[df.pheno == ph]
        T2[(key, ph)] = block(s, CVA2[key]); show(f'{key} {ph}', T2[(key, ph)])
pickle.dump(T2, open('work/table2_D.pkl', 'wb'))
S = {}
for key, df in [('real', q), ('sim', sp)]:
    for sx in ['Female', 'Male']:
        for grp in ['All', 'unc', 'neph']:
            s = df[df.Gender == sx] if grp == 'All' else df[(df.Gender == sx) & (df.pheno == grp)]
            S[(key, sx, grp)] = block(s, CVA2[key]); r = S[(key, sx, grp)]; r['n_raw'] = r['n']
            show(f'{key} {sx} {grp}', r)
pickle.dump(S, open('work/tableS_sex_D.pkl', 'wb'))
print('done')
