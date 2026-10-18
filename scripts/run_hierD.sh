#!/bin/bash
mkdir -p work/hierD
run(){ echo "=== $* $(date)"; env CHAINS=4 DRAWS=1000 TUNE=1000 "$@" 2>&1 | grep -v -i warn | grep -E "model pairs|took|RCVcond|Error|rror" ; }
run python3 -I fit_hier_D.py ext real
run TA=0.97 python3 -I fit_hier_D.py ext sim
run python3 -I fit_hier_D.py ext real Female
run python3 -I fit_hier_D.py ext real Male
run COVS=mcv,hb,egfr,age,female TAG=nohba1c python3 -I fit_hier_D.py sel real
run TA=0.97 python3 -I fit_hier_D.py ext sim Female
run TA=0.97 python3 -I fit_hier_D.py ext sim Male
run COVS=hba1c TAG=hba1c python3 -I fit_hier_D.py sel real
run HBMEAN=1 TAG=ext_pairmean python3 -I fit_hier_D.py ext real
echo "ALLDONE $(date)"
