# -*- coding: utf-8 -*-
"""Fast ML fit of a location-scale Student-t (nu, mu, sigma) by ECME (Liu & Rubin 1995), vectorized."""
import numpy as np
from scipy.special import digamma
from scipy.optimize import brentq

def t_fit(x, tol=1e-8, maxit=500, nu0=3.0):
    x = np.asarray(x, float); n = len(x)
    mu = np.median(x); sg = 1.4826 * np.median(np.abs(x - mu)); nu = nu0
    for it in range(maxit):
        z2 = ((x - mu) / sg) ** 2
        w = (nu + 1.0) / (nu + z2)
        mu_new = np.sum(w * x) / np.sum(w)
        sg_new = np.sqrt(np.sum(w * (x - mu_new) ** 2) / n)
        # CML step for nu given (mu, sigma): solve  -digamma(nu/2)+log(nu/2)+1+mean(log w - w)+digamma((nu+1)/2)-log((nu+1)/2)=0
        z2 = ((x - mu_new) / sg_new) ** 2
        w = (nu + 1.0) / (nu + z2)
        c = np.mean(np.log(w) - w) + 1.0
        f = lambda v: -digamma(v / 2) + np.log(v / 2) + c + digamma((v + 1) / 2) - np.log((v + 1) / 2)
        try:
            nu_new = brentq(f, 0.1, 1000.0)
        except ValueError:
            nu_new = 1000.0 if f(1000.0) > 0 else 0.1
        done = abs(mu_new - mu) < tol * sg and abs(sg_new - sg) < tol * sg and abs(nu_new - nu) < 1e-6 * max(1, nu)
        mu, sg, nu = mu_new, sg_new, nu_new
        if done: break
    return nu, mu, sg

if __name__ == '__main__':
    from scipy import stats
    import pandas as pd, time
    q = pd.read_pickle('work/pairs_D.pkl')
    for ph in ['unc', 'neph']:
        g = q[q.pheno == ph]['delta'].values
        t0 = time.time(); a = t_fit(g); t1 = time.time(); b = stats.t.fit(g); t2 = time.time()
        print(ph, 'ECME', np.round(a, 5), round(t1 - t0, 3), 's | scipy', np.round(b, 5), round(t2 - t1, 3), 's')
