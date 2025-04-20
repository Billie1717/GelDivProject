# Explaining the model

- The model is run in lammps and uses the REACTION package. The reactions allow the addition and removal of beads as well as the change in 'states' of the beads. This allows us to simulate the dynamic growth and shrinkage of microtubules and altering these dynamics with other chemical species. 
- Some parameters, particularly the parameters governing MT dynamics, were explored previously and we subsequently keep fixed as they are broadly mapped to experimental values. These are listed below. 
- It is best when understanding the model to first look at "CuboidGeometry" then "AsterGeometry" then "MitoIntoGap" as these are increasing in complexity and the simpler parameters introduced in CuboidGeometry are not re-commented in the more complex setups. 
- Unfortunately, the nucleation pattern of the microtubules is coded in the lammps source code. I have set up a handful of useful nucleation patterns which can be called in the fix bond/react nucleation line. Only if one wants to nucleate in a different way than those already coded, you must go into the source code and create a new nucleation 'rule'. Nucleation is explained below.
- Please first download and combile the correct Lammps version with the altered fix_bond_react source codes, as described in Lammps/LammpsReadme.md


## Parameters kept fixed:


## Coding the MT nucleation:

- Nucleation of an MT is written in the input file with in the bond/react fix as :
$ react Nucleation all ${rstep} 0.900000 1.100000 mPreNucleation mPostNucleation Reactions/Maps/map_Nucleation.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} modify_create overlap ${OvaR1} ${OvaR2} modify_create nuc xwalls -2 &
- Particularly important is the flag that comes after 'nuc'. Here xwalls means we nucleation microtubules with a random y- and z- coordinate, but at each x-wall and they will be nucleated with the polarity pointing inwards (plus-end pointing towards the middle of the simulation box). Here, the flag '-2' has no physical meaning but is just the number assigned to this geometry. 
- OvaR1 should always be ~0.9 and this is the minimum distance between a newly added MT bead and any other MT bead
- OvaR2 specifies the minimum distance between a newly added MT and a bead of type 7. This should be set to 0.5*(sigma + sigma_d) where sigma = MT diameter and sigma_d = mito diameter
- Please ask me how more overlap distances are specified for 'size sorting' experiments. 
- See NucleationInstructions.txt for the documentary on other nucleation geometries.



## (note to self) Currently a working simulation showing pressure-driven flow:
- Pushing : Scratch/CytoDiv/Reaction_model/Asters_Molecules/runs_reactionsLoops/runsReactions_ron_6.0_roff_10.0_rnuc_0.5_rGTP_0.5_rGDP_0.5_rbr_0.4_brAng_12_numPhis_15_N_2100_sigD_6.0_Ra_6.0_epsDD_1.0_epsDF_5.0_lamda_0.2_OvaR1_0.7_OvaR2_3.5_mass_1.0_rstep_1_Pdist_35.0_vfrac_0.7_seed_2/output2.xyz
- No pushing : Scratch/CytoDiv/Reaction_model/Asters_Molecules/runs_noreactionsLoops/runsNoReactions_ron_6.0_roff_10.0_rnuc_0.5_rGTP_0.5_rGDP_0.5_rbr_1.0_brAng_12_numPhis_15_N_2100_sigD_6.0_Ra_6.0_epsDD_1.0_epsDF_5.0_lamda_0.2_OvaR1_0.7_OvaR2_3.5_mass_1.0_factor_10_Pdist_35.0_vfrac_0.7_seed_1
