#!/bin/bash
cdir=${PWD}
ron=4.0
roff=8.0
rnuc=0.5 #0.1 #1.0
rGTP=0.5
rGDP=0.5
roff=10.0
ron=6.0
brAng=12
seed=2
rbr=1.0 #0.4
roff=10.0
numPhis=15
numPhisM=14
Ndiff=2100
sigma_d=6.0
r_a=6.0
mass_d=150.0
xbox=70.0 #30.0
yzbox=50.0
eps_dd=3.0
lamda=0.5
eps_dd=0.5
r_a=6.0
eps_df=5.0
lamda=0.2
OvaR1=0.7
OvaR2=3.5
temp=1.0
factor=10
volfract=0.8
Pdist=35.0
for mass_d in 1.0 #5.0 #10.0 50.0 #5.0 10.0 50.0 #10.0 #10.0 #50.0 
do
for eps_dd in 1.0 #1.0 3.0 #2.0 3.0
do
for volfract in 0.7
do
for seed in 1
do
filepattern='runs_restart/runsReactions_ron_'${ron}'_roff_'${roff}'_rnuc_'${rnuc}'_rGTP_'${rGTP}'_rGDP_'${rGDP}'_rbr_'${rbr}'_brAng_'${brAng}'_numPhis_'${numPhis}'_N_'${Ndiff}'_sigD_'${sigma_d}'_Ra_'${r_a}'_epsDD_'${eps_dd}'_epsDF_'${eps_df}'_lamda_'${lamda}'_OvaR1_'${OvaR1}'_OvaR2_'${OvaR2}'_mass_'${mass_d}'_factor_'${factor}'_Pdist_'${Pdist}'_vfrac_'${volfract}'_seed_'${seed}'/'
mkdir ${filepattern}
cp -r 'Inputs/Reactions/' ${filepattern}
#cp -r 'Inputs/restartsPostSpindleNoReactepsDD0.5/spindle.restart.11500000' ${filepattern}'/restart.postwall'
cp -r 'Inputs/restartsPostSpindleReactionsepsDD1.0/spindle.restart.9500000' ${filepattern}'/restart.postwall'
cd ${filepattern}
python /nfs/scistore15/saricgrp/bmeadowc/Scratch/CytoDiv/Reaction_model/Asters_Molecules/build_spindle_diffusors_deletemidzoneNoReactionsLoopsRS.py ${ron} ${roff} ${rnuc} ${rGTP} ${rGDP} ${rbr} ${brAng} ${numPhis} ${sigma_d} ${r_a} ${eps_dd} ${eps_df} ${lamda} ${OvaR1} ${OvaR2} ${factor} ${Pdist} ${mass_d} ${seed}
for i in $(seq 0 $numPhisM)
do
filefile='post_branch1_'${i}'.txt'
mv -f ${filefile} 'Reactions/Branch'
done
sbatch runscriptDel.sh
cd ..
cd ..
done
done
done
done
