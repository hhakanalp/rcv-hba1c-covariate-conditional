import pandas as pd, numpy as np, json
out={}
for key,src,pp in [('real','work/sheet_revize.pkl','work/pairs_D.pkl'),('sim','work/SIM_v3.pkl','work/pairs_D_sim.pkl')]:
    d=pd.read_pickle(src); d['Date']=pd.to_datetime(d['Date']); d.loc[d['HbA1c percent']==56.0,'HbA1c percent']=5.6
    p=pd.read_pickle(pp)
    keys=set(zip(p['Patient ID'],p['Date']))|set(zip(p['Patient ID'],p['D2']))
    m=d[[ (a,b) in keys for a,b in zip(d['Patient ID'],d['Date'])]].copy()
    ph=p.drop_duplicates('Patient ID').set_index('Patient ID')['pheno']; m['pheno']=m['Patient ID'].map(ph)
    res={}
    for grp in ['All','unc','neph']:
        mm=m if grp=='All' else m[m.pheno==grp]; q=p if grp=='All' else p[p.pheno==grp]
        pt=mm.drop_duplicates('Patient ID'); first=mm.sort_values('Date').groupby('Patient ID').first()
        npp=q.groupby('Patient ID').size()
        res[grp]=dict(patients=int(pt['Patient ID'].nunique()),age_mean=float(first['Age'].mean()),age_sd=float(first['Age'].std()),
            female=int((first['Gender']=='Female').sum()),female_pct=float((first['Gender']=='Female').mean()*100),
            measurements=int(len(mm)),meas_per_pt=float(len(mm)/pt['Patient ID'].nunique()),pairs=int(len(q)),
            pairs_per_pt_median=float(npp.median()),pairs_per_pt_iqr=(float(npp.quantile(.25)),float(npp.quantile(.75))),pct_one_pair=float((npp==1).mean()*100),
            egfr_med=float(mm['eGFR'].median()),acr_med=float(mm['ACR mg/g'].median()),hb_med=float(mm['Hemoglobin'].median()),mcv_med=float(mm['MCV_fL'].median()),hba1c_med=float(mm['HbA1c percent'].median()),
            days_med=float(q.days.median()),days_iqr=(float(q.days.quantile(.25)),float(q.days.quantile(.75))),acr_visits_med=float(mm.groupby('Patient ID')['ACR mg/g'].count().median()))
    out[key]=res
    print(key); print(pd.DataFrame(res).T.to_string())
json.dump(out,open('work/table1_D.json','w'),indent=1)
