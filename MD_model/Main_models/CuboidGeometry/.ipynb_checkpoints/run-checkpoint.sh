#!/bin/bash
cdir=${PWD}
rnuc=0.03 #Nucleation rate at poles
rbr=0.5 #Branching rate (on other MTs)
Ndiff=400 #Number of mitochondria (re-defined below)
sigma_d=6.0 #Diameter of mitocondria (for volume exclusion)
eps_dd=3.0 #Interaction strength between two mitochondria (Cosine-squared potential to ensure fluidity)
eps_df=1.0 #Interaction strength between MT and mitochondria (Turns out not to matter too much)
factor=100 #The MT reaction rate modulation, e.g for changing the ratio between reaction rates and diffusion. Rates ~1/Factor. Check in build code
mass_d=1.0 #The mass of the mitochondria, alter to change diffusion speed
numPhisM=14 #The number of possible branching angles, 14 is enough
seed=1 #Simulation seed
for eps_dd in 1.0 #0.1 1.0 2.0 3.0 4.0 5.0
do
for eps_df in 1.0 #0.1 1.0 2.0 3.0 4.0 5.0
do
for seed in 1 #2 3 
do
for Ndiff in 400 #200 300
do 
filepattern='runs/runsRc1.5_rnuc_'${rnuc}'_rbr_'${rbr}'_N_'${Ndiff}'_sigD_'${sigma_d}'_epsDD_'${eps_dd}'_epsDF_'${eps_df}'_mass_'${mass_d}'_factor_'${factor}'_seed_'${seed}'/'
mkdir ${filepattern}
cp -r 'Inputs/Reactions/' ${filepattern}
cd ${filepattern}
python /nfs/scistore15/saricgrp/bmeadowc/Scratch/CytoDiv/Reaction_model/Experiment1/build_cuboid_mito.py ${rnuc} ${rbr} ${sigma_d} ${eps_dd} ${eps_df} ${factor} ${seed} #Change to your own directory
python /nfs/scistore15/saricgrp/bmeadowc/Scratch/CytoDiv/Reaction_model/Experiment1/build_config.py ${Ndiff} ${sigma_d} ${mass_d} #Change to your own directory
for i in $(seq 0 $numPhisM)
do
filefile='post_branch1_'${i}'.txt'
mv -f ${filefile} 'Reactions/Branch/'
done
sbatch runscript.sh
cd ..
cd ..
done
done
done
done