# -*- coding: utf-8 -*-
"""Regenerate Figures 2, 3, 5, 6, 7 (600 dpi) from reproduced data + hierarchical model results."""
import pandas as pd, numpy as np, json, pickle, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy import stats
import os
MODE=os.environ.get('MODE','A')
NOTITLE=os.environ.get('NOTITLE','0')=='1'
if NOTITLE:
    import matplotlib.figure as _mf; _mf.Figure.suptitle=lambda self,*a,**k: None
OUT='work/figs_D'; os.makedirs(OUT,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titleweight':'bold'})
BLUE='#1f77b4'; ORANGE='#d95f02'; GREY='#6e6e6e'; LBLUE='#9ecae1'; LORANGE='#fdae6b'; DORANGE='#a63603'

T2=pickle.load(open('work/table2_D.pkl','rb'))
Rr=json.load(open('work/hierD/real/ext.json')); Rs=json.load(open('work/hierD/sim/ext.json'))
p=pd.read_pickle('work/pairs_D.pkl'); sp=pd.read_pickle('work/pairs_D_sim.pkl')
def tscale(g): nu,mu,sg=stats.t.fit(g); return sg
cond={'real':{'unc':(Rr['rcv_unc']['up'],Rr['rcv_unc']['dn']),'neph':(Rr['rcv_neph']['up'],Rr['rcv_neph']['dn'])},
      'sim':{'unc':(Rs['rcv_unc']['up'],Rs['rcv_unc']['dn']),'neph':(Rs['rcv_neph']['up'],Rs['rcv_neph']['dn'])}}
marg={k:{g:(T2[(k,g)]['up'],T2[(k,g)]['dn']) for g in ['unc','neph','All']} for k in ['real','sim']}

# ================= Figure 2 =================
fig,axes=plt.subplots(1,2,figsize=(15,7.6))
fig.suptitle('Figure 2. Log-difference (Δ) distributions by phenotype with marginal and conditional RCV limits',fontsize=15,fontweight='bold',y=0.985)
for ax,df,key,title in zip(axes,[p,sp],['real','sim'],['A. Real data','B. Synthetic data']):
    bins=np.linspace(-0.6,0.6,121)
    for ph,col in [('unc',BLUE),('neph',ORANGE)]:
        g=df[df.pheno==ph]['delta'].values
        ax.hist(g,bins=bins,density=True,alpha=0.5,color=col)
        up,dn=marg[key][ph]
        for v in (np.log(1+up/100),np.log(1+dn/100)): ax.axvline(v,color=col,ls='-',lw=2)
        cu,cd=cond[key][ph]
        for v in (np.log(1+cu/100),np.log(1+cd/100)): ax.axvline(v,color=col,ls='-.',lw=1.6)
        r_=T2[(key,ph)]; xx=np.linspace(-0.6,0.6,600)
        ax.plot(xx,stats.t.pdf(xx,r_['nu'],r_['mu'],r_['sig_t']),color=col,ls='--',lw=1.4)
    ax.set_xlim(-0.62,0.62); ax.set_xlabel('Log-difference Δ (ln HbA1c$_{t+1}$ / HbA1c$_t$)',fontsize=13); ax.set_ylabel('Density',fontsize=13)
    ax.set_title(title,loc='left',fontsize=14); ax.spines[['top','right']].set_visible(False)
    nu_,nn=T2[(key,'unc')],T2[(key,'neph')]; cu,cn=cond[key]['unc'],cond[key]['neph']
    nunc=int((df.pheno=='unc').sum()); nneph=int((df.pheno=='neph').sum())
    ax.text(0.5,-0.25,f"Uncomplicated (n = {nunc:,}):  RCV$_{{Phenotype}}$ +{nu_['up']:.1f}/−{abs(nu_['dn']):.1f}%   •   RCV$_{{Cond}}$ +{cu[0]:.1f}/−{abs(cu[1]):.1f}%",transform=ax.transAxes,ha='center',fontsize=12,color=BLUE)
    ax.text(0.5,-0.33,f"Nephropathy (n = {nneph:,}):  RCV$_{{Phenotype}}$ +{nn['up']:.1f}/−{abs(nn['dn']):.1f}%   •   RCV$_{{Cond}}$ +{cn[0]:.1f}/−{abs(cn[1]):.1f}%",transform=ax.transAxes,ha='center',fontsize=12,color=ORANGE)
handles=[Patch(fc=BLUE,alpha=0.5,label='Uncomplicated T2DM'),Patch(fc=ORANGE,alpha=0.5,label='T2DM + nephropathy'),
         Line2D([],[],color='#444',ls='-',lw=2,label='RCV$_{Phenotype}$ limits (marginal)'),Line2D([],[],color='#444',ls='-.',lw=1.6,label='RCV$_{Cond}$ limits (conditional)'),
         Line2D([],[],color='#444',ls='--',lw=1.4,label='Fitted Student-t density')]
fig.legend(handles=handles,loc='lower center',ncol=5,fontsize=11.5,frameon=False,bbox_to_anchor=(0.5,0.012))
plt.subplots_adjust(top=0.90,bottom=0.30,left=0.06,right=0.98,wspace=0.2)
fig.savefig(f'{OUT}/Figure2_delta_600dpi.png',dpi=600,facecolor='white'); plt.close(fig)

# ================= Figure 3 =================
fig,axes=plt.subplots(2,2,figsize=(13,8.2))
fig.suptitle('Figure 3. Central-scale CV$_I$ and asymmetric RCV with 95% bootstrap CIs',fontsize=15,fontweight='bold',y=0.985)
labels=['All','Uncomplicated\nT2DM','T2DM +\nnephropathy']; cols=[GREY,BLUE,ORANGE]; ys=[2,1,0]
for row,key,ptitle in [(0,'real','A. Real data'),(1,'sim','B. Synthetic data')]:
    axL,axR=axes[row]
    for g,y,c in zip(['All','unc','neph'],ys,cols):
        r=T2[(key,g)]
        axL.errorbar(r['cvi'],y,xerr=[[r['cvi']-r['cvi_ci'][0]],[r['cvi_ci'][1]-r['cvi']]],fmt='o',color=c,ms=9,capsize=4,lw=2)
        axL.text(r['cvi_ci'][1]+(0.04 if key=='real' else 0.02),y+0.13,f"{r['cvi']:.2f}",color=c,fontsize=11,va='bottom')
        axR.errorbar(r['up'],y+0.12,xerr=[[r['up']-r['up_ci'][0]],[r['up_ci'][1]-r['up']]],fmt='s',color=c,ms=8,capsize=3,lw=1.8)
        axR.errorbar(r['dn'],y-0.12,xerr=[[r['dn']-r['dn_ci'][0]],[r['dn_ci'][1]-r['dn']]],fmt='D',color=c,ms=8,capsize=3,lw=1.8)
        axR.text(r['up_ci'][1]+(0.4 if key=='real' else 0.12),y+0.12,f"+{r['up']:.1f}",color=c,fontsize=11,va='center')
        axR.text(r['dn_ci'][0]-(0.4 if key=='real' else 0.12),y-0.12,f"{r['dn']:.1f}",color=c,fontsize=11,va='center',ha='right')
    axL.set_yticks(ys); axL.set_yticklabels(labels,fontsize=11); axL.set_ylim(-0.7,2.7); axL.set_xlabel('CV$_I$ (%)',fontsize=12)
    axL.set_title('Within-subject CV$_I$',fontsize=12.5); axL.spines[['top','right']].set_visible(False)
    axR.axvline(0,color='#999',lw=1); axR.set_yticks(ys); axR.set_yticklabels([]); axR.set_ylim(-0.7,2.7); axR.set_xlabel('RCV (%)',fontsize=12)
    axR.set_title('Asymmetric RCV (■ upper / ◆ lower)',fontsize=12.5); axR.spines[['top','right']].set_visible(False)
    if key=='real': axL.set_xlim(5.0,7.0); axR.set_xlim(-20,24)
    else: axL.set_xlim(1.2,2.1); axR.set_xlim(-6.5,7.0)
    if key=='sim':
        axL.axvline(1.6,color='#444',ls=':',lw=1.4)
        axL.text(1.6,2.55,'injected CV$_I$ = 1.6 (both phenotypes)',color='#444',fontsize=9.5,ha='center')
    axL.text(-0.32,1.08,ptitle,transform=axL.transAxes,fontsize=13.5,fontweight='bold')
plt.subplots_adjust(top=0.90,bottom=0.08,left=0.13,right=0.98,hspace=0.55,wspace=0.12)
fig.savefig(f'{OUT}/Figure3_cvi_rcv_600dpi.png',dpi=600,facecolor='white'); plt.close(fig)

# ================= Figure 5 =================
# panel A: |delta| vs eGFR on real central pairs (per-phenotype trimming)
tr=p.dropna(subset=['eGFR'])
fig,axes=plt.subplots(1,2,figsize=(15,6.6),gridspec_kw={'width_ratios':[1,1.1]})
fig.suptitle('Figure 5. Hierarchical Student-t dispersion model (real data, tail-trimmed central scale)',fontsize=15,fontweight='bold',y=0.985)
ax=axes[0]; ax.scatter(tr.eGFR,tr.delta.abs(),s=3,alpha=0.12,color=BLUE,rasterized=True)
edges=np.arange(15,175,10); mids=[]; meds=[]
for lo,hi in zip(edges[:-1],edges[1:]):
    m=(tr.eGFR>=lo)&(tr.eGFR<hi)
    if m.sum()>=30: mids.append((lo+hi)/2); meds.append(tr.delta[m].abs().median())
ax.plot(mids,meds,'-o',color=ORANGE,lw=2.2,ms=6,label='Binned median |Δ| (n ≥ 30 per bin)')
ax.set_xlabel('eGFR (mL/min/1.73 m²)',fontsize=13); ax.set_ylabel('|Δ| (log-difference)',fontsize=13); ax.set_ylim(0,0.6); ax.set_xlim(10,175)
ax.set_title('A. |Δ| versus eGFR (all pairs)',loc='left',fontsize=13.5); ax.legend(frameon=False,fontsize=11,loc='upper left'); ax.spines[['top','right']].set_visible(False)
# panel B: forest
ax=axes[1]
terms=[('Baseline HbA1c (β₁)','hba1c'),('MCV (β₂)','mcv'),('Hemoglobin (β₃)','hb'),('eGFR (β₄)','egfr'),('Age (β₅)','age'),('Female sex (β₆)','female')]
rows=[(lab,Rr['beta'][k]['mean'],Rr['beta'][k]['lo'],Rr['beta'][k]['hi']) for lab,k in terms]+[('Nephropathy (γ₁)',Rr['g1']['mean'],Rr['g1']['lo'],Rr['g1']['hi'])]
n=len(rows); ypos=np.arange(n)[::-1]
ax.axvline(0,color='#999',lw=1)
for (lab,m,lo,hi),y in zip(rows,ypos):
    sig=(lo>0)|(hi<0); col=ORANGE if 'DKD' in lab else (BLUE if sig else GREY)
    ax.errorbar(m,y,xerr=[[m-lo],[hi-m]],fmt='o',color=col,ms=9,capsize=4,lw=2)
    if m>0.2: ax.text(lo-0.008,y+0.02,f"{m:+.3f} [{lo:+.3f}, {hi:+.3f}]",color=col,fontsize=10.5,va='center',ha='right')
    else: ax.text(hi+0.008,y+0.02,f"{m:+.3f} [{lo:+.3f}, {hi:+.3f}]",color=col,fontsize=10.5,va='center')
ax.set_yticks(ypos); ax.set_yticklabels([r[0] for r in rows],fontsize=11.5); ax.set_xlim(-0.2,0.42)
ax.set_xlabel('Posterior coefficient on ln σ (per SD; 95% CrI)',fontsize=12.5); ax.set_title('B. Dispersion model coefficients',loc='left',fontsize=13.5)
ax.set_ylim(-0.7,n-0.3)

fig.text(0.5,0.025,f"Residual ν = {Rr['nu']['mean']:.1f}; {Rr['n_pairs']:,} pairs from {Rr['n_patients']:,} patients; continuous covariates standardized (per SD), sex centered",ha='center',fontsize=10.5,color='#555',style='italic')
ax.spines[['top','right']].set_visible(False)
plt.subplots_adjust(top=0.88,bottom=0.16,left=0.06,right=0.99,wspace=0.42)
fig.savefig(f'{OUT}/Figure5_model_600dpi.png',dpi=600,facecolor='white'); plt.close(fig)

# ================= Figure 6 =================
fig,axes=plt.subplots(1,2,figsize=(14,6.4))
fig.suptitle('Figure 6. Conditional phenotype-specific RCV (RCV$_{Cond}$) versus marginal RCV$_{Phenotype}$ (central scale)',fontsize=14.5,fontweight='bold',y=0.985)
for ax,key,title in zip(axes,['real','sim'],['A. Real data','B. Synthetic data']):
    x=np.array([0,1]); w=0.19
    for i,ph in enumerate(['unc','neph']):
        mu_,md=marg[key][ph]; cu,cd=cond[key][ph]
        bars=[(x[i]-1.5*w,mu_,BLUE),(x[i]-0.5*w,md,LBLUE),(x[i]+0.5*w,cu,DORANGE),(x[i]+1.5*w,cd,LORANGE)]
        for xx,v,c in bars:
            ax.bar(xx,v,w*0.95,color=c,edgecolor='none')
            ax.text(xx,v+(0.5 if v>0 else -0.5)*(1 if key=='real' else 0.3),f"{v:+.1f}",ha='center',va='bottom' if v>0 else 'top',fontsize=11)
    ax.axhline(0,color='#333',lw=1); ax.set_xticks(x); ax.set_xticklabels(['Uncomplicated T2DM','T2DM + nephropathy'],fontsize=12)
    ax.set_ylabel('RCV (%)',fontsize=13); ax.set_title(title,loc='left',fontsize=13.5); ax.spines[['top','right']].set_visible(False)
    ax.set_ylim((-22,25) if key=='real' else (-8,9))
handles=[Patch(fc=BLUE,label='RCV$_{Phenotype}$ upper'),Patch(fc=LBLUE,label='RCV$_{Phenotype}$ lower'),Patch(fc=DORANGE,label='RCV$_{Cond}$ upper'),Patch(fc=LORANGE,label='RCV$_{Cond}$ lower')]
fig.legend(handles=handles,loc='lower center',ncol=4,fontsize=12,frameon=False,bbox_to_anchor=(0.5,0.01))
plt.subplots_adjust(top=0.88,bottom=0.2,left=0.07,right=0.98,wspace=0.22)
fig.savefig(f'{OUT}/Figure6_cond_vs_marginal_600dpi.png',dpi=600,facecolor='white'); plt.close(fig)

# ================= Figure 7 =================
recl={}
for key,df in [('real',p),('sim',sp)]:
    n=df[df.pheno=='neph']['delta'].values; N=len(n)
    su=T2[(key,'unc')]['sd']; sn=T2[(key,'neph')]['sd']
    cu,cd=cond[key]['neph']; lu,ld=np.log(1+cu/100),np.log(1+cd/100)
    fixed=int(((n>1.96*su)|(n<-1.96*su)).sum()); pheno=int(((n>1.96*sn)|(n<-1.96*sn)).sum()); condn=int(((n>lu)|(n<ld)).sum())
    recl[key]=dict(N=N,fixed=fixed,pheno=pheno,cond=condn)
json.dump(recl,open('work/reclass_D.json','w'),indent=1); print('reclass',recl)
fig,axes=plt.subplots(1,2,figsize=(14,6.4))
fig.suptitle('Figure 7. Flagging of nephropathy pairs under fixed, phenotype-specific and conditional RCV limits',fontsize=14.5,fontweight='bold',y=0.985)
for ax,key,title in zip(axes,['real','sim'],['A. Real data (nephropathy pairs)','B. Synthetic data (nephropathy pairs)']):
    r=recl[key]; vals=[r['fixed']/r['N']*100,r['pheno']/r['N']*100,r['cond']/r['N']*100]; ns=[r['fixed'],r['pheno'],r['cond']]
    bars=ax.bar([0,1,2],vals,0.62,color=[GREY,ORANGE,DORANGE],edgecolor='#333',lw=1)
    for b,v,nn in zip(bars,vals,ns): ax.text(b.get_x()+b.get_width()/2,v+0.6,f"{v:.1f}%\n(n = {nn:,})",ha='center',va='bottom',fontsize=11,fontweight='bold')
    ax.set_xticks([0,1,2]); ax.set_xticklabels(['Fixed uncomplicated\nRCV$_{Phenotype}$ limits','Phenotype-specific\nRCV$_{Phenotype}$ limits','Conditional\nRCV$_{Cond}$ limits'],fontsize=11.5)
    ax.set_ylabel('Nephropathy pairs flagged (%)',fontsize=13); ax.set_ylim(0,32); ax.set_title(title,loc='left',fontsize=13.5); ax.spines[['top','right']].set_visible(False)
fig.text(0.5,0.095,f"Discordant pairs (flagged only under the fixed uncomplicated limit): {recl['real']['fixed']-recl['real']['pheno']} real / {recl['sim']['fixed']-recl['sim']['pheno']} simulated.\n"
         f"Every pair flagged by the wider phenotype-specific limits is also flagged by the fixed limit (structural property). RCV$_{{Cond}}$ is evaluated at the mean covariate profile.",ha='center',fontsize=10.5,style='italic',color='#444')
plt.subplots_adjust(top=0.88,bottom=0.25,left=0.07,right=0.98,wspace=0.22)
fig.savefig(f'{OUT}/Figure7_reclass_600dpi.png',dpi=600,facecolor='white'); plt.close(fig)
print('figures done')
