# Explaining the model

These simulation set-ups result in asters growing towards walls of mitochondria, a gap forming between two growing asters and mitochondria flowing into this gap (for 'Pushing') or not flowing into this gap (for 'NoPushing'). In the 'No pushing' example, the MTs are frozen when the gap appears such that they are no longer growing towards the mitochondria boundaries. In the 'Pushing' simulations, the MTs grow as normal throughout. The additional detail from the 'aster' setup in the adjacent folder is a man-made gap which is enforced some way through the simulation. This gap is maintained by periodically deleting all MT beads in a fixed region at the centre of the simulation box. A snapshot of the simulation setup post-gap can be found in XXX

## Currently a working simulation showing pressure-driven flow:
- Pushing : Scratch/CytoDiv/Reaction_model/Asters_Molecules/runs_reactionsLoops/runsReactions_ron_6.0_roff_10.0_rnuc_0.5_rGTP_0.5_rGDP_0.5_rbr_0.4_brAng_12_numPhis_15_N_2100_sigD_6.0_Ra_6.0_epsDD_1.0_epsDF_5.0_lamda_0.2_OvaR1_0.7_OvaR2_3.5_mass_1.0_rstep_1_Pdist_35.0_vfrac_0.7_seed_2/output2.xyz
- No pushing : Scratch/CytoDiv/Reaction_model/Asters_Molecules/runs_noreactionsLoops/runsNoReactions_ron_6.0_roff_10.0_rnuc_0.5_rGTP_0.5_rGDP_0.5_rbr_1.0_brAng_12_numPhis_15_N_2100_sigD_6.0_Ra_6.0_epsDD_1.0_epsDF_5.0_lamda_0.2_OvaR1_0.7_OvaR2_3.5_mass_1.0_factor_10_Pdist_35.0_vfrac_0.7_seed_1







## initial conditions


## inital data file


## input file

the file input_example.in is not complete for running a simulation but contains all the ingredients of the input file included comments, it skips the many 'cycles' of runs explained below.

the file input.in is complete and should be able to be run with the correct lammps version + src code edits (Note : it is a very long lammps input file!)

Because the command fix create and fix break for bonds in lammps shouldn't be used simultaneously, the main point to note is that the code requires bonds to be make and broken separately, with type-IDs of atoms being updated in between. These cycles have to be short enough so that atoms types can be updated sufficiently frequently without a single atom bonding more than once. This makes the lammps input script very long as it may have 100s of cycles of breaking and making bonds. 

