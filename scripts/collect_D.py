import json, pickle, numpy as np, pandas as pd
W='work/'
T2=pickle.load(open(W+'table2_D.pkl','rb')); S=pickle.load(open(W+'tableS_sex_D.pkl','rb'))
H={k:json.load(open(W+f'hierD/{k}/ext.json')) for k in ['real','sim','real_Female','real_Male','sim_Female','sim_Male']}
H['nohba1c']=json.load(open(W+'hierD/real/nohba1c.json')); H['hba1c']=json.load(open(W+'hierD/real/hba1c.json')); H['pairmean']=json.load(open(W+'hierD/real/ext_pairmean.json'))
R=dict(T2={f'{k[0]}|{k[1]}':v for k,v in T2.items()},S={f'{k[0]}|{k[1]}|{k[2]}':v for k,v in S.items()},diff=json.load(open(W+'diff_D.json')),reclass=json.load(open(W+'reclass_D.json')),H=H,
       flow=dict(real=json.load(open(W+'flow_real.json')),sim=json.load(open(W+'flow_sim.json'))),table1=json.load(open(W+'table1_D.json')),ademp=json.load(open(W+'ademp_grid.json')),truth=json.load(open(W+'SIM_v3_truth.json')),six=json.load(open(W+'six_month.json')))
# steady: measurements per patient & intervals
p=pd.read_pickle(W+'pairs_D.pkl'); R['extra']=dict(days_med=float(p.days.median()),days_iqr=[float(p.days.quantile(.25)),float(p.days.quantile(.75))],days_unc=float(p[p.pheno=='unc'].days.median()),days_neph=float(p[p.pheno=='neph'].days.median()))
json.dump(R,open(W+'resultsD.json','w'),indent=1,default=float)
for k,v in H.items(): print(f"{k:12s} pairs {v['n_pairs']} pts {v['n_patients']} nu {v['nu']['mean']:.1f} tau {v['tau']['mean']:.2f} g1 {v['g1']['mean']:+.3f} ({v['g1']['lo']:+.3f},{v['g1']['hi']:+.3f}) | "+' '.join(f"{n} {b['mean']:+.3f}({b['lo']:+.3f},{b['hi']:+.3f})" for n,b in v['beta'].items())+f" | RCVc unc {v['rcv_unc']['up']:+.1f}/{v['rcv_unc']['dn']:.1f} neph {v['rcv_neph']['up']:+.1f}/{v['rcv_neph']['dn']:.1f} | rhat {v['diag']['max_rhat']:.3f} ess {v['diag']['min_ess']:.0f} div {v['diag']['divergences']}")
for k,v in R['T2'].items(): print(k, f"n {v['n']} pts {v['n_pat']} nu {v['nu']:.2f} cvtot {v['cvtot']:.2f} cvi {v['cvi']:.2f} ({v['cvi_ci'][0]:.2f}-{v['cvi_ci'][1]:.2f}) up {v['up']:+.1f} ({v['up_ci'][0]:.1f}-{v['up_ci'][1]:.1f}) dn {v['dn']:.1f} ({v['dn_ci'][0]:.1f} to {v['dn_ci'][1]:.1f}) classical {v['classical']:.1f} w<0.5 {v['w_lt05']:.1f} gmm {v['gmm_cvi']:.2f} jones {v['jones_cvi']:.2f} keep {v['jones_keep']*100:.0f}% days {v['median_days']:.0f}")
print(R['extra']); print(R['reclass'])
