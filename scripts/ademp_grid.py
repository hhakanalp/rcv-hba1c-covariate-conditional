# -*- coding: utf-8 -*-
"""ADEMP Monte Carlo grid for the pair-level estimators on synthetic v3 cohorts (pairs 1-180 d), 10 replicates each.
Estimators: A Student-t scale (pipeline D); B Student-t fit then |z|<=2 retained-pair SD (previous pipeline);
C adaptive central 95% Student-t band; D two-component Gaussian mixture central component; E Jones-type iterative truncation.
usage: python3 -I ademp_grid.py <out.json> [reps]"""
import subprocess, sys, json, numpy as np, pandas as pd, os, tempfile
from scipy import stats
sys.path.insert(0, 'scripts'); from tfit import t_fit
from sklearn.mixture import GaussianMixture
OUT = sys.argv[1]; REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 10
SC = {
 'G1 null, no clinical change':            dict(p_step=0.0, p_epi=0.0),
 'G2 null, main contamination (v3)':       dict(),
 'G3 null, heavy contamination':           dict(p_step=0.30, p_epi=0.10),
 'G4 injected difference 1.6 vs 2.24':     dict(cvi_neph=2.24),
 'G5 real-like 5.5 vs 6.0, CVA 1.1':       dict(cvi_unc=5.5, cvi_neph=6.0, cva=1.1, step_sd=0.12),
 'G6 null, Laplace biological noise':      dict(noise='laplace'),
 'G7 null, t4 biological noise':           dict(noise='t4'),
 'G8 null, mixture biological noise':      dict(noise='mixture'),
 'G9 null, lognormal (skewed) noise':      dict(noise='lognormal'),
}
def est_all(x, cva):
    nu, mu, sg = t_fit(x); z = (x - mu) / sg
    cvi = lambda s: float(np.sqrt(max((s / np.sqrt(2) * 100) ** 2 - cva ** 2, 0)))
    A = cvi(sg); B = cvi(x[np.abs(z) <= 2].std(ddof=1)); C = cvi(x[np.abs(z) <= stats.t.ppf(0.975, nu)].std(ddof=1))
    g = GaussianMixture(2, n_init=3, random_state=0).fit(x.reshape(-1, 1)); k = int(np.argmin(g.covariances_.ravel())); D = cvi(float(np.sqrt(g.covariances_.ravel()[k])))
    keep = np.ones(len(x), bool)
    for _ in range(200):
        m = np.median(x[keep]); sd = x[keep].std(ddof=1) / 0.8711; new = np.abs(x - m) <= 1.96 * sd
        if np.all(new == keep): break
        keep = new
    E = cvi(sd)
    return dict(A=A, B=B, C=C, D=D, E=E, nu=float(nu), kept_B=float((np.abs(z) <= 2).mean()), kept_E=float(keep.mean()))
res = {}
tmp = tempfile.mkdtemp()
for name, kw in SC.items():
    rows = []
    for rep in range(REPS):
        args = ['python3', '-I', 'scripts/make_sim_v3.py', f'{tmp}/s.pkl', '--seed', str(1000 + rep), '--quiet'] + sum([[f'--{k}', str(v)] for k, v in kw.items()], [])
        subprocess.run(args, check=True)
        subprocess.run(['python3', '-I', 'scripts/build_cohort.py', f'{tmp}/s.pkl', f'{tmp}/p.pkl'], check=True, stdout=subprocess.DEVNULL)
        p = pd.read_pickle(f'{tmp}/p.pkl'); truth = json.load(open(f'{tmp}/s_truth.json'))
        for ph in ['unc', 'neph']:
            e = est_all(p[p.pheno == ph]['delta'].values, truth['cva']); e.update(rep=rep, pheno=ph, truth=truth[f'cvi_{ph}'], n=int((p.pheno == ph).sum())); rows.append(e)
    df = pd.DataFrame(rows); summ = {}
    for ph in ['unc', 'neph']:
        s = df[df.pheno == ph]; t = float(s.truth.iloc[0])
        summ[ph] = dict(truth=t, n_pairs=float(s.n.mean()), nu=float(s.nu.mean()), **{k: dict(mean=float(s[k].mean()), sd=float(s[k].std()), bias_pct=float((s[k].mean() - t) / t * 100)) for k in 'ABCDE'}, kept_B=float(s.kept_B.mean()), kept_E=float(s.kept_E.mean()))
    r = df.pivot(index='rep', columns='pheno', values='A'); summ['ratio_A'] = dict(mean=float((r['neph'] / r['unc']).mean()), truth=float(summ['neph']['truth'] / summ['unc']['truth']))
    rb = df.pivot(index='rep', columns='pheno', values='B'); summ['ratio_B'] = float((rb['neph'] / rb['unc']).mean())
    res[name] = summ
    print(f"{name:40s} " + ' | '.join(f"{ph}: truth {summ[ph]['truth']:.2f} A {summ[ph]['A']['mean']:.2f} ({summ[ph]['A']['bias_pct']:+.0f}%) B {summ[ph]['B']['mean']:.2f} ({summ[ph]['B']['bias_pct']:+.0f}%) C {summ[ph]['C']['mean']:.2f} D {summ[ph]['D']['mean']:.2f} E {summ[ph]['E']['mean']:.2f} nu {summ[ph]['nu']:.1f}" for ph in ['unc', 'neph']) + f" | ratio A {summ['ratio_A']['mean']:.2f} B {summ['ratio_B']:.2f} (truth {summ['ratio_A']['truth']:.2f})", flush=True)
    json.dump(res, open(OUT, 'w'), indent=1)
