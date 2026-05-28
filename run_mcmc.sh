#!/bin/bash
rm -f results/outputs/tep_cobaya_sne.*
export PYTHONUNBUFFERED=1
mpirun -n 2 python scripts/steps/step_03_04_cobaya_mcmc.py 2>&1 | perl -pe 's/\x1b\[[0-9;]*[mGK]//g' | tr -d '\000' > logs/step_03_04_cobaya_mcmc.log &
