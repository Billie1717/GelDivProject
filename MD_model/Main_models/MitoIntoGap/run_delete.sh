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
rbr=0.4
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
factor=100
volfract=0.7
Pdist=35.0
for mass_d in 1.0 5.0 #10.0 #1.0 5.0 10.0 50.0 #10.0 #10.0 #50.0 
do
for rstep in 1 10 #100 #3.0 6.0
do
for eps_dd in 0.5 1.0 3.0 #0.5 1.0 2.0 3.0
do
for seed in 2
do
filepattern='runs_reactionsLoops/runsReactions_ron_'${ron}'_roff_'${roff}'_rnuc_'${rnuc}'_rGTP_'${rGTP}'_rGDP_'${rGDP}'_rbr_'${rbr}'_brAng_'${brAng}'_numPhis_'${numPhis}'_N_'${Ndiff}'_sigD_'${sigma_d}'_Ra_'${r_a}'_epsDD_'${eps_dd}'_epsDF_'${eps_df}'_lamda_'${lamda}'_OvaR1_'${OvaR1}'_OvaR2_'${OvaR2}'_mass_'${mass_d}'_rstep_'${rstep}'_Pdist_'${Pdist}'_vfrac_'${volfract}'_seed_'${seed}'/'
mkdir ${filepattern}
cp -r 'Inputs/Reactions/' ${filepattern}
#cp -r 'Inputs/restart.postasterDense' ${filepattern}
cp -r 'Inputs/restartsPostAsterPbr0.4/spindle.restart.7500000' ${filepattern}'/restart.postaster'
cd ${filepattern}
#echo "1" ${ron} "2" ${roff} "3" ${rnuc} "4" ${rGTP} "5" ${rGDP} "6" ${rbr} "7" ${brAng} ${numPhis} ${sigma_d} ${r_a} ${eps_dd} ${eps_df} ${lamda} ${OvaR1} ${OvaR2} "F" ${factor} "P" ${Pdist} "S" ${seed}
python /nfs/scistore15/saricgrp/bmeadowc/Scratch/CytoDiv/Reaction_model/Asters_Molecules/build_spindle_diffusors_deletemidzone_ReactionsLoops.py ${ron} ${roff} ${rnuc} ${rGTP} ${rGDP} ${rbr} ${brAng} ${numPhis} ${sigma_d} ${r_a} ${eps_dd} ${eps_df} ${lamda} ${OvaR1} ${OvaR2} ${rstep} ${Pdist} ${mass_d} ${seed}
#python /nfs/scistore15/saricgrp/bmeadowc/Scratch/CytoDiv/Reaction_model/Asters_Molecules/build_spindle_diffusors_deletemidzoneNoReactions.py ${ron} ${roff} ${rnuc} ${rGTP} ${rGDP} ${rbr} ${brAng} ${numPhis} ${sigma_d} ${r_a} ${eps_dd} ${eps_df} ${lamda} ${OvaR1} ${OvaR2} ${factor} ${Pdist} ${mass_d} ${seed}
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
