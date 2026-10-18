# -*- coding: utf-8 -*-
"""Pipeline D (pairs <=180 d): hierarchical Student-t dispersion model on ALL pairs (no |z|<=2 trimming).
Baseline (Table 3) and extended (+age, +sex) fits. NUTS 2 chains x (700 tune + 700 draws)."""
import pandas as pd, numpy as np, json, sys, time
import pymc as pm, arviz as az
from scipy import stats
OUT='work/hierD'
import os; os.makedirs(OUT,exist_ok=True)
DATA=sys.argv[2] if len(sys.argv)>2 else 'real'
if DATA=='real':
    p=pd.read_pickle('work/pairs_D.pkl')
else:
    p=pd.read_pickle('work/pairs_D_sim.pkl')
SEX=sys.argv[3] if len(sys.argv)>3 else None
if SEX: p=p[p['Gender']==SEX].copy()
OUT=OUT+'/'+DATA+('_'+SEX if SEX else ''); os.makedirs(OUT,exist_ok=True)
tr=p.copy()
tr=tr.dropna(subset=['HbA1c percent','MCV_fL','Hemoglobin','eGFR','Age','Gender']).copy()
print('model pairs',len(tr),'patients',tr['Patient ID'].nunique(), tr.pheno.value_counts().to_dict(), flush=True)
y=tr['delta'].values
neph=(tr.pheno=='neph').astype(float).values
def z(v): v=v.astype(float).values; return (v-v.mean())/v.std()
hb_cov=(tr['HbA1c percent']+tr['X2'])/2 if os.environ.get('HBMEAN')=='1' else tr['HbA1c percent']
X_base=np.column_stack([z(hb_cov),z(tr['MCV_fL']),z(tr['Hemoglobin']),z(tr['eGFR'])])
fem=(tr['Gender']=='Female').astype(float).values
X_ext=np.column_stack([X_base,z(tr['Age'])]) if SEX else np.column_stack([X_base,z(tr['Age']),fem-fem.mean()])
pid=pd.factorize(tr['Patient ID'])[0]; npat=int(pid.max()+1)
def fit(X,names,tag):
    t0=time.time()
    with pm.Model() as m:
        g0=pm.Normal('g0',0,1); g1=pm.Normal('g1',0,1)
        beta=pm.Normal('beta',0,1,shape=X.shape[1])
        tau=pm.HalfNormal('tau',0.5); u_raw=pm.Normal('u_raw',0,1,shape=npat); u=pm.Deterministic('u',u_raw*tau)
        mu=pm.Normal('mu',0,0.1); nu=pm.Gamma('nu',2,0.1)
        lnsig=g0+g1*neph+pm.math.dot(X,beta)+u[pid]
        pm.StudentT('y',nu=nu,mu=mu,sigma=pm.math.exp(lnsig),observed=y)
        idata=pm.sample(draws=int(os.environ.get('DRAWS','1000')),tune=int(os.environ.get('TUNE','1000')),chains=int(os.environ.get('CHAINS','4')),cores=2,target_accept=float(os.environ.get('TA','0.9')),random_seed=42,progressbar=False,nuts_sampler=os.environ.get('SAMPLER','pymc'))
    post=idata.posterior
    summ=az.summary(idata,var_names=['g0','g1','beta','tau','mu','nu'],ci_prob=0.95)
    print(tag,'\n',summ.to_string(),flush=True)
    res={'n_pairs':int(len(y)),'n_patients':int(npat),'time_min':round((time.time()-t0)/60,1),'names':names}
    for v in ['g0','g1','tau','mu','nu']:
        s=post[v].values.ravel(); res[v]={'mean':float(s.mean()),'sd':float(s.std()),'lo':float(np.percentile(s,2.5)),'hi':float(np.percentile(s,97.5))}
    b=post['beta'].values.reshape(-1,X.shape[1])
    res['beta']={n:{'mean':float(b[:,i].mean()),'sd':float(b[:,i].std()),'lo':float(np.percentile(b[:,i],2.5)),'hi':float(np.percentile(b[:,i],97.5))} for i,n in enumerate(names)}
    # RCVcond at mean profile (posterior of sigma_c), CVA=1.5
    g0s=post['g0'].values.ravel(); g1s=post['g1'].values.ravel()
    for lab,ind in [('unc',0),('neph',1)]:
        sc=np.exp(g0s+g1s*ind)
        up=(np.exp(1.96*sc)-1)*100; dn=(np.exp(-1.96*sc)-1)*100; cvi=np.sqrt(np.maximum((sc*100/np.sqrt(2))**2-(0.64 if DATA=='sim' else 1.21),0))
        res[f'rcv_{lab}']={'up':float(np.median(up)),'up_lo':float(np.percentile(up,2.5)),'up_hi':float(np.percentile(up,97.5)),
                          'dn':float(np.median(dn)),'dn_lo':float(np.percentile(dn,2.5)),'dn_hi':float(np.percentile(dn,97.5)),
                          'cvi':float(np.median(cvi)),'cvi_lo':float(np.percentile(cvi,2.5)),'cvi_hi':float(np.percentile(cvi,97.5)),'sigma_c':float(np.median(sc))}
    res['diag']={'max_rhat':float(summ['r_hat'].max()),'min_ess':float(summ['ess_bulk'].min()),'divergences':int(idata.sample_stats['diverging'].values.sum())}
    json.dump(res,open(f'{OUT}/{tag}.json','w'),indent=1)
    az.to_netcdf(idata.sel(draw=slice(None,None,1)),f'{OUT}/{tag}.nc') if False else None
    print(tag,'RCVcond',res['rcv_unc'],res['rcv_neph'],res['diag'],flush=True)
    return res
which=sys.argv[1] if len(sys.argv)>1 else 'both'
if which=='sel':
    names_all=['hba1c','mcv','hb','egfr','age','female']; X_all=np.column_stack([X_base,z(tr['Age']),fem-fem.mean()])
    sel=[c for c in os.environ.get('COVS','').split(',') if c]
    idx=[names_all.index(c) for c in sel]
    Xs=X_all[:,idx] if idx else np.zeros((len(y),0))
    fit(Xs,sel,os.environ.get('TAG','sel'))
if which in ('base','both'): fit(X_base,['hba1c','mcv','hb','egfr'],'base')
if which in ('ext','both'): fit(X_ext,['hba1c','mcv','hb','egfr','age'] if SEX else ['hba1c','mcv','hb','egfr','age','female'],os.environ.get('TAG','ext'))
