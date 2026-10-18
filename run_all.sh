#!/bin/bash
# Runs the pipeline on the synthetic cohort (demo) or on your own export. Run from the repository root. Python 3.11 + requirements.txt.
set -e
# Real data (not distributed): export your LIS data in the layout of data/INPUT_SCHEMA.md to work/real.xlsx, then uncomment:
# python3 -I scripts/build_cohort.py work/real.xlsx work/pairs_D.pkl work/flow_real.json
python3 -I scripts/build_cohort.py data/Synthetic_HbA1c_v3_CVI1.6_CVA0.8.xlsx work/pairs_D_sim.pkl work/flow_sim.json   # synthetic cohort (demo; used as the main set when no real export is given)
[ -f work/pairs_D.pkl ] || cp work/pairs_D_sim.pkl work/pairs_D.pkl
python3 -I scripts/build_D.py                                                   # Table 2 + benchmarks  -> work/table2_D.pkl, tableS_sex_D.pkl
python3 -I scripts/diff_D.py                                                    # bootstrap CIs of between-phenotype differences -> work/diff_D.json
bash scripts/run_hierD.sh                                                       # hierarchical Student-t dispersion models (PyMC) -> work/hierD/*/ext.json
python3 -I scripts/ademp_grid.py work/ademp_grid.json 10                        # ADEMP simulation grid (Supplemental Table 9)
python3 -I scripts/collect_D.py                                                 # collects everything -> work/resultsD.json
python3 -I scripts/make_main_figs_D.py                                          # Figures 1-3 -> work/figs_D
python3 -I scripts/make_supp_figs_D.py                                          # Supplemental Figures 2-5
python3 -I scripts/make_flow_D.py                                               # Supplemental Figure 1 (flow chart)
