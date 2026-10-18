# -*- coding: utf-8 -*-
"""Regenerate Figures 1, 2, 4 at 600 dpi with corrected cohort counts (56.0 -> 5.6)."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy import stats
import os
MODE='B'
NOTITLE=os.environ.get('NOTITLE','0')=='1'
if NOTITLE:
    import matplotlib.figure as _mf; _mf.Figure.suptitle=lambda self,*a,**k: None
OUT='work/figs_D'; os.makedirs(OUT,exist_ok=True)
plt.rcParams['font.family']='DejaVu Sans'
BLUE='#1f77b4'; ORANGE='#d95f02'

# ---------------- Figure 1 ----------------
def box(ax,x,y,w,h,text,fill='#eef3f7',edge='#2f3e46',lw=1.8,fs=10.8,pad=0.0):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle=f'round,pad={pad},rounding_size=0.02',fc=fill,ec=edge,lw=lw))
    ax.text(x,y,text,ha='center',va='center',fontsize=fs,family='DejaVu Serif')
def arrow(ax,x1,y1,x2,y2,color='#2f3e46'):
    ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle='-|>',lw=1.8,color=color,mutation_scale=16))
fig,axes=plt.subplots(1,2,figsize=(14,12))
fig.suptitle('Figure 1. Study data isolation flow (real and simulated datasets)',fontsize=16,fontweight='bold',family='DejaVu Serif',y=0.97)
import json
fr=json.load(open('work/flow_real.json')); fs=json.load(open('work/flow_sim.json'))
def steps_for(f,sim=False):
    return [(f"All {'synthetic' if sim else 'LIS HbA1c'} records\n{f['records']:,} measurements / {f['patients']:,} patients",'n'),
      (f"T2DM eligible (≥1 HbA1c ≥6.5%)\n{f['t2dm_meas']:,} measurements / {f['t2dm_pts']:,} patients",'n'),
      (f"Exclude uninterpretable measurements\n(Hb<10 g/dL, MCV<70 fL, eGFR<15): −{f['excluded_meas']:,}\n{f['after_excl_meas']:,} measurements / {f['after_excl_pts']:,} patients",'x'),
      (f"≥2 retained visits → consecutive pairs\n{f['ge2_meas']:,} measurements / {f['ge2_pts']:,} patients / {f['pairs_all']:,} pairs",'n'),
      (f"Drop no ACR result ({f['noACR_pts']} patients, −{f['noACR_pairs']} pairs) and\nunconfirmed albuminuria, single ACR ≥30 mg/g\n({f['unconf_pts']} patients, −{f['unconf_pairs']:,} pairs)",'x'),
      (f"Phenotyped cohort, all intervals: {f['phenotyped_pts']:,} patients / {f['phenotyped_pairs']:,} pairs\nUncomplicated {f['unc_pairs_nolimit']:,}  |  Nephropathy {f['neph_pairs_nolimit']:,}",'n'),
      (f"Restrict to inter-measurement interval 1–180 days: −{f['pairs_gt180']:,} pairs\nAnalysis cohort: {f['D_pts']:,} patients / {f['D_meas']:,} measurements / {f['D_pairs']:,} pairs",'x')]
real=steps_for(fr); sim=steps_for(fs,True)
bottoms=[(f"Uncomplicated T2DM\n{fr['unc_pairs']:,} pairs / {fr['unc_pts']:,} patients",f"T2DM with nephropathy\n{fr['neph_pairs']:,} pairs / {fr['neph_pts']:,} patients"),
         (f"Uncomplicated T2DM\n{fs['unc_pairs']:,} pairs / {fs['unc_pts']:,} patients",f"T2DM with nephropathy\n{fs['neph_pairs']:,} pairs / {fs['neph_pts']:,} patients")]
for ax,steps,title,bottom in zip(axes,[real,sim],['A. Real-data (LIS) cohort','B. Synthetic validation cohort (v3)'],bottoms):
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
    ax.text(0.5,0.955,title,ha='center',fontsize=15,fontweight='bold',family='DejaVu Serif')
    ys=np.linspace(0.905,0.225,7); w=0.94; h=0.078
    for (t,k),y in zip(steps,ys):
        box(ax,0.5,y,w,h,t,fill='#fdeee2' if k=='x' else '#eef3f7')
    for y1,y2 in zip(ys[:-1],ys[1:]): arrow(ax,0.5,y1-h/2-0.005,0.5,y2+h/2+0.008)
    yb=0.085
    box(ax,0.27,yb,0.40,0.072,bottom[0],fill='white',edge=BLUE,lw=2.2)
    box(ax,0.73,yb,0.40,0.072,bottom[1],fill='white',edge=ORANGE,lw=2.2)
    arrow(ax,0.42,ys[-1]-h/2-0.003,0.3,yb+0.04,color=BLUE); arrow(ax,0.58,ys[-1]-h/2-0.003,0.7,yb+0.04,color=ORANGE)
plt.subplots_adjust(top=0.93,bottom=0.02,left=0.03,right=0.97,wspace=0.08)
fig.savefig(f'{OUT}/FigS1_flow_600dpi.png',dpi=600,facecolor='white'); plt.close(fig)
print('flow done')
