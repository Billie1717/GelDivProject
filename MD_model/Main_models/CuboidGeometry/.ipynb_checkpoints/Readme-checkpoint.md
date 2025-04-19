# Intro
This folder contains the scripts required to run simulations where MT networks grow from two 'poles' which are the square ends of a cuboid. The figure that these scripts produced is 'MitoAsBarrier.pdf'. The mitochondria is initialised asymmetrically so that the role of mitochondria barrier, and not distance from the pole, is properly explored. 

## The location of some simulation runs on my cluster space:

- From heatplot of mito-mito vs MT-mito interaction strength showing we need above a threshold to see good excl. zone maintenance,
    eps_mito-mito = 1, eps_mito-MT=1 : Scratch/CytoDiv/Reaction_model/Experiment1/runsExp1/runs_rnuc_0.03_rbr_0.5_N_400_sigD_6.0_epsDD_1.0_epsDF_1.0_mass_1.0_factor_100_seed_1/output.xyz
    eps_mito-mito = 3, eps_mito-MT=1 : Scratch/CytoDiv/Reaction_model/Experiment1/runsExp1/runs_rnuc_0.03_rbr_0.5_N_400_sigD_6.0_epsDD_3.0_epsDF_1.0_mass_1.0_factor_100_seed_1/output.xyz
    
The plot was created with the parameters currently specified in the run.sh file, with the sweep over each epsilon: 0.1 1.0 2.0 3.0 4.0 5.0

## build_spindle_diffusors.py

Syntax : $ python build_cuboid_mito.py ${rnuc} ${rbr} ${sigma_d} ${eps_dd} ${eps_df} ${factor} ${seed}
(parameters are explained in the run.sh script)

This python file does 3 things. 1) it creates the all-important lammps input script including all the specified parameters (specified in, e.g. run.sh which calls this python script). 2) it creates a runscript.sh which is a slurm script for use on the cluster. 3) it creates some extra reaction templates for all the branching angles you want (we set 14 branching angles so this makes 14 extra reaction templates). If you only had one branching angle then your MTs would end up very ordered. Our MTs randomly branch at any of the 14 angles. If you only care about the input script, you can find an example commented input script example.in  

## build_config.py

Syntax : $ python build_config.py ${Ndiff} ${sigma_d} ${mass_d}

Creates the lammps configuration file which consists of the nucleation 'walls' and the mitochondria. Takes the arguments : {# mitochondria} {mito diameter} {mito mass} (as can be seen when it is called in run.sh). An example output with 400 mitochondria is in the file configuration.txt