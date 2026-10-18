import pandas as pd, numpy as np, json, sys
sys.path.insert(0,'scripts'); from tfit import t_fit
W='work/'
d=pd.read_pickle(W+'sheet_revize.pkl'); d['Date']=pd.to_datetime(d['Date']); d.loc[d['HbA1c percent']==56.0,'HbA1c percent']=5.6
p=pd.read_pickle(W+'pairs_D.pkl'); pn=pd.read_pickle(W+'pairs_D_nolimit.pkl')
keys=set(zip(p['Patient ID'],p['Date']))|set(zip(p['Patient ID'],p['D2']))
m=d[[(a,b) in keys for a,b in zip(d['Patient ID'],d['Date'])]].copy()
ph=p.drop_duplicates('Patient ID').set_index('Patient ID')['pheno']; m['pheno']=m['Patient ID'].map(ph)
out={'sex':{}}
for sx in ['Female','Male']:
    mm=m[m.Gender==sx]; first=mm.sort_values('Date').groupby('Patient ID').first(); q=p[p.Gender==sx]
    def miqr(s): return f"{s.median():.1f} ({s.quantile(.25):.1f}–{s.quantile(.75):.1f})"
    out['sex'][sx]=dict(patients=int(first.shape[0]),age_mean=float(first.Age.mean()),age_sd=float(first.Age.std()),age_med=float(first.Age.median()),age_q1=float(first.Age.quantile(.25)),age_q3=float(first.Age.quantile(.75)),
        neph_n=int((first.pheno=='neph').sum()),neph_pct=float((first.pheno=='neph').mean()*100),measurements=int(len(mm)),meas_per_pt=float(len(mm)/first.shape[0]),pairs=int(len(q)),
        hba1c=miqr(mm['HbA1c percent']),egfr=f"{mm.eGFR.median():.0f} ({mm.eGFR.quantile(.25):.0f}–{mm.eGFR.quantile(.75):.0f})",acr=miqr(mm['ACR mg/g'].dropna()),hb=miqr(mm.Hemoglobin),mcv=miqr(mm.MCV_fL))
# interval sensitivity
CVA2=1.21
def est(x):
    nu,mu,sg=t_fit(x); cvt=sg/np.sqrt(2)*100; return dict(n=int(len(x)),nu=float(nu),cvi=float(np.sqrt(max(cvt**2-CVA2,0))),up=float((np.exp(1.96*sg)-1)*100),dn=float((np.exp(-1.96*sg)-1)*100))
out['interval']={}
for lab,lo,hi in [('1–90',1,90),('91–180',91,180),('181–365',181,365),('>365',366,100000),('1–180 (main)',1,180),('all intervals',1,100000)]:
    for g in ['unc','neph']:
        s=pn[(pn.pheno==g)&(pn.days>=lo)&(pn.days<=hi)]; r=est(s.delta.values); r['pts']=int(s['Patient ID'].nunique()); r['days_med']=float(s.days.median()); out['interval'][f'{lab}|{g}']=r
        print(f"{lab:14s} {g:5s} n={r['n']:6d} pts {r['pts']:5d} med {r['days_med']:4.0f} d nu {r['nu']:.2f} CVI {r['cvi']:.2f} RCV {r['up']:+.1f}/{r['dn']:.1f}")
json.dump(out,open(W+'supp_numbers_D.json','w'),indent=1)
print(json.dumps(out['sex'],indent=1)[:1500])
