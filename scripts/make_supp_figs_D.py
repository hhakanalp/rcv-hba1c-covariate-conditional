# -*- coding: utf-8 -*-
import pandas as pd, numpy as np, json, pickle, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os
MODE=os.environ.get('MODE','A')
NOTITLE=os.environ.get('NOTITLE','0')=='1'
if NOTITLE:
    import matplotlib.figure as _mf; _mf.Figure.suptitle=lambda self,*a,**k: None
OUT='work/figs_D'; os.makedirs(OUT,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titleweight':'bold'})
BLUE='#1f77b4'; ORANGE='#d95f02'; GREY='#6e6e6e'; LBLUE='#9ecae1'; LORANGE='#fdae6b'; DORANGE='#a63603'
S=pickle.load(open('work/tableS_sex_D.pkl','rb'))
M={(ds,sx):json.load(open(f"work/hierD/{ds}_{sx}/ext.json")) for ds in ['real','sim'] for sx in ['Female','Male']}

# ---------- Figure S1: CVI and RCV by sex x phenotype ----------
fig,axes=plt.subplots(2,2,figsize=(13,9))
fig.suptitle('Figure S1. Central-scale CV$_I$ and asymmetric RCV$_{Phenotype}$ by sex and phenotype (95% bootstrap CIs)',fontsize=14.5,fontweight='bold',y=0.985)
rows=[('Female','unc','Women — uncomplicated',BLUE,'o'),('Male','unc','Men — uncomplicated',BLUE,'s'),('Female','neph','Women — nephropathy',ORANGE,'o'),('Male','neph','Men — nephropathy',ORANGE,'s')]
for r_,key,ptitle in [(0,'real','A. Real data'),(1,'sim','B. Simulated data')]:
    axL,axR=axes[r_]; ys=[3,2,1,0]
    for (sx,ph,lab,c,mk),y in zip(rows,ys):
        r=S[(key,sx,ph)]
        axL.errorbar(r['cvi'],y,xerr=[[r['cvi']-r['cvi_ci'][0]],[r['cvi_ci'][1]-r['cvi']]],fmt=mk,color=c,ms=9,capsize=4,lw=2)
        axL.text(r['cvi_ci'][1]+(0.05 if key=='real' else 0.03),y+0.12,f"{r['cvi']:.2f}",color=c,fontsize=11,va='bottom')
        axR.errorbar(r['up'],y+0.12,xerr=[[r['up']-r['up_ci'][0]],[r['up_ci'][1]-r['up']]],fmt=mk,color=c,ms=8,capsize=3,lw=1.8)
        axR.errorbar(r['dn'],y-0.12,xerr=[[r['dn']-r['dn_ci'][0]],[r['dn_ci'][1]-r['dn']]],fmt=mk,color=c,ms=8,capsize=3,lw=1.8,mfc='white')
        axR.text(r['up_ci'][1]+(0.5 if key=='real' else 0.15),y+0.12,f"+{r['up']:.1f}",color=c,fontsize=11,va='center')
        axR.text(r['dn_ci'][0]-(0.5 if key=='real' else 0.15),y-0.12,f"{r['dn']:.1f}",color=c,fontsize=11,va='center',ha='right')
    axL.set_yticks(ys); axL.set_yticklabels([r[2] for r in rows],fontsize=11); axL.set_ylim(-0.7,3.7); axL.set_xlabel('CV$_I$ (%)',fontsize=12); axL.set_title('Within-subject CV$_I$',fontsize=12.5)
    axR.axvline(0,color='#999',lw=1); axR.set_yticks(ys); axR.set_yticklabels([]); axR.set_ylim(-0.7,3.7); axR.set_xlabel('RCV (%)',fontsize=12); axR.set_title('Asymmetric RCV (filled: upper / open: lower)',fontsize=12.5)
    for a in (axL,axR): a.spines[['top','right']].set_visible(False)
    if key=='real': axL.set_xlim(4.6,7.4); axR.set_xlim(-21,25)
    else: axL.set_xlim(1.1,2.3); axR.set_xlim(-7,8)
    axL.text(-0.42,1.06,ptitle,transform=axL.transAxes,fontsize=13.5,fontweight='bold')
plt.subplots_adjust(top=0.90,bottom=0.07,left=0.17,right=0.98,hspace=0.5,wspace=0.12)
fig.savefig(f'{OUT}/FigureS1_sex_cvi_rcv_600dpi.png',dpi=600,facecolor='white'); plt.close(fig)

# ---------- Figure S2: dispersion model by sex (coefficients + RCVcond) ----------
fig,axes=plt.subplots(1,2,figsize=(15,6.4),gridspec_kw={'width_ratios':[1.25,1]})
fig.suptitle('Figure S2. Sex-stratified hierarchical Student-t dispersion models (real data)',fontsize=14.5,fontweight='bold',y=0.985)
ax=axes[0]
terms=[('Baseline HbA1c (β₁)','hba1c'),('MCV (β₂)','mcv'),('Hemoglobin (β₃)','hb'),('eGFR (β₄)','egfr'),('Age (β₅)','age')]
labels=[t[0] for t in terms]+['Nephropathy (γ₁)']; n=len(labels); ypos=np.arange(n)[::-1]
ax.axvline(0,color='#999',lw=1)
for sx,off,c,mk in [('Female',+0.17,'#7b3294','o'),('Male',-0.17,'#008837','s')]:
    R=M[('real',sx)]
    vals=[(R['beta'][k]['mean'],R['beta'][k]['lo'],R['beta'][k]['hi']) for _,k in terms]+[(R['g1']['mean'],R['g1']['lo'],R['g1']['hi'])]
    for (m,lo,hi),y in zip(vals,ypos):
        ax.errorbar(m,y+off,xerr=[[m-lo],[hi-m]],fmt=mk,color=c,ms=8,capsize=3,lw=1.8)
        if m>0.2: ax.text(lo-0.006,y+off,f"{m:+.3f} [{lo:+.3f}, {hi:+.3f}]",color=c,fontsize=9.5,va='center',ha='right')
        else: ax.text(hi+0.006,y+off,f"{m:+.3f} [{lo:+.3f}, {hi:+.3f}]",color=c,fontsize=9.5,va='center')
ax.set_yticks(ypos); ax.set_yticklabels(labels,fontsize=11.5); ax.set_xlim(-0.15,0.37); ax.set_ylim(-0.7,n+0.5)
ax.set_xlabel('Posterior coefficient on ln σ (per SD; 95% CrI)',fontsize=12.5); ax.set_title('A. Model coefficients by sex',loc='left',fontsize=13.5); ax.spines[['top','right']].set_visible(False)
ax.legend(handles=[plt.Line2D([],[],marker='o',color='#7b3294',ls='',ms=8,label=f"Women (n = {M[('real','Female')]['n_pairs']:,} pairs, {M[('real','Female')]['n_patients']:,} patients)"),
                   plt.Line2D([],[],marker='s',color='#008837',ls='',ms=8,label=f"Men (n = {M[('real','Male')]['n_pairs']:,} pairs, {M[('real','Male')]['n_patients']:,} patients)")],loc='upper left',frameon=False,fontsize=10.5)
ax=axes[1]; x=np.array([0,1]); w=0.19
for i,sx in enumerate(['Female','Male']):
    R=M[('real',sx)]; mu_,md=S[('real',sx,'unc')]['up'],S[('real',sx,'unc')]['dn']; mn,mdn=S[('real',sx,'neph')]['up'],S[('real',sx,'neph')]['dn']
    bars=[(x[i]-1.5*w,mu_,BLUE),(x[i]-0.5*w,R['rcv_unc']['up'],LBLUE),(x[i]+0.5*w,mn,ORANGE),(x[i]+1.5*w,R['rcv_neph']['up'],LORANGE)]
    lows=[(x[i]-1.5*w,md,BLUE),(x[i]-0.5*w,R['rcv_unc']['dn'],LBLUE),(x[i]+0.5*w,mdn,ORANGE),(x[i]+1.5*w,R['rcv_neph']['dn'],LORANGE)]
    for xx,v,c in bars+lows:
        ax.bar(xx,v,w*0.95,color=c,edgecolor='none'); ax.text(xx,v+(0.5 if v>0 else -0.5),f"{v:+.1f}",ha='center',va='bottom' if v>0 else 'top',fontsize=9.5)
ax.axhline(0,color='#333',lw=1); ax.set_xticks(x); ax.set_xticklabels(['Women','Men'],fontsize=12); ax.set_ylabel('RCV (%)',fontsize=13); ax.set_ylim(-22,26)
ax.set_title('B. Marginal vs conditional RCV by sex',loc='left',fontsize=13.5); ax.spines[['top','right']].set_visible(False)
ax.legend(handles=[Patch(fc=BLUE,label='Uncomplicated RCV$_{Phenotype}$'),Patch(fc=LBLUE,label='Uncomplicated RCV$_{Cond}$'),Patch(fc=ORANGE,label='Nephropathy RCV$_{Phenotype}$'),Patch(fc=LORANGE,label='Nephropathy RCV$_{Cond}$')],loc='lower center',ncol=2,frameon=False,fontsize=10,bbox_to_anchor=(0.5,-0.32))
plt.subplots_adjust(top=0.88,bottom=0.2,left=0.11,right=0.99,wspace=0.3)
fig.savefig(f'{OUT}/FigureS2_sex_model_600dpi.png',dpi=600,facecolor='white'); plt.close(fig)
print('supp figs done')
