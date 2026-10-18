import pandas as pd, numpy as np, sys, json
sys.path.insert(0,'scripts'); from tfit import t_fit
d=pd.read_pickle('work/sheet_0.pkl'); d.loc[d['HbA1c percent']==56.0,'HbA1c percent']=5.6
elig=d.groupby('Patient ID')['HbA1c percent'].transform('max')>=6.5
r=d[elig & d['Exclusion'].isna()].copy(); r=r[r.groupby('Patient ID')['Patient ID'].transform('size')>=2].sort_values(['Patient ID','Date'])
n30=d.groupby('Patient ID')['ACR mg/g'].apply(lambda a:(a>=30).sum()); np_=r.groupby('Patient ID')['Nephropathy_patient'].first()
ph=pd.Series(np.where(np_.str.startswith('None'),'unc',np.where(np_.str.startswith('No ACR'),'noACR',np.where(n30.reindex(np_.index).fillna(0)>=2,'neph','unconf'))),index=np_.index)
g=r.groupby('Patient ID'); r['X2']=g['HbA1c percent'].shift(-1); r['D2']=g['Date'].shift(-1); r['first']=g['Date'].transform('min')
p=r.dropna(subset=['X2']).copy(); p['pheno']=p['Patient ID'].map(ph); p=p[p.pheno.isin(['unc','neph'])]
p['delta']=np.log(p.X2/p['HbA1c percent']); p['days']=(p.D2-p.Date).dt.days
CVA=1.1
def est(x):
    nu,mu,sg=t_fit(x); cvt=sg/np.sqrt(2)*100; return np.sqrt(max(cvt**2-CVA**2,0)),nu,sg
def boot(s,B=1000,seed=1):
    rng=np.random.default_rng(seed); x=s.delta.values; ids,inv=np.unique(s['Patient ID'].values,return_inverse=True)
    grp=[np.where(inv==i)[0] for i in range(len(ids))]; out=[]
    for b in range(B): out.append(est(x[np.concatenate([grp[i] for i in rng.integers(0,len(ids),len(ids))])])[0])
    return np.percentile(out,[2.5,97.5])
sets={
 'A) interval 1-180 d (pair gap)': p[(p.days>=1)&(p.days<=180)],
 'B) interval <= calendar 6 months': p[(p.days>=1)&(p.D2<=p.Date+pd.DateOffset(months=6))],
 'C) first 6 months of follow-up (both measurements within 180 d of first visit)': p[(p.D2-p['first']).dt.days<=180],
 'D) no interval limit': p,
}
res={}
for k,s in sets.items():
    print(k)
    for grp in ['All','unc','neph']:
        q=s if grp=='All' else s[s.pheno==grp]
        cvi,nu,sg=est(q.delta.values); ci=boot(q,B=1000 if k.startswith('A') else 300)
        up,dn=(np.exp(1.96*sg)-1)*100,(np.exp(-1.96*sg)-1)*100
        print(f"   {grp:5s} pts {q['Patient ID'].nunique():5d} pairs {len(q):6d} nu {nu:5.2f} CVI {cvi:.2f} ({ci[0]:.2f}-{ci[1]:.2f}) RCV {up:+.1f}/{dn:.1f}  median gap {q.days.median():.0f} d",flush=True)
        res[f'{k}|{grp}']=dict(pts=int(q['Patient ID'].nunique()),pairs=len(q),nu=nu,cvi=cvi,ci=ci.tolist(),rcv=[up,dn])
json.dump(res,open('work/six_month.json','w'),indent=1)
