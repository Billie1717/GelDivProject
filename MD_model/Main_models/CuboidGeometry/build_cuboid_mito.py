from __future__ import division, print_function
import numpy as np
#import matplotlib.pyplot as plt
import math
import time
import random
import sys, getopt
import os

#'''
#
#'''
    
def main():
    ################ Individual filament parameters ################
    ron = 6.0
    roff = 10.0
    rnuc = float(sys.argv[1])
    rGTP = 0.5
    rGDP = 0.5
    rbr = float(sys.argv[2])
    BrAng= 12
    numPhis = 15
    sigma_d = float(sys.argv[3])
    R_a = sigma_d
    eps_dd = float(sys.argv[4])
    eps_df = float(sys.argv[5])
    lamda = 0.2
    OvaR1 = 0.9
    OvaR2 = 3.5
    factor = float(sys.argv[6])
    seed = float(sys.argv[7])
    
    Epsilon = 0.5 
    EpsilonB = 0.0
    kbond = 1000
    kbend = 1000
    tstep = 0.0005
    run_time = 4000.0 
    run_time1 = 1000.0
    frame_rate = 20.0

    ##Calculating branch placement
    sigma = 1.0
    sigma_df = (sigma +sigma_d)/2
    # VOLUME EXCLUSION ONLY
    rc = sigma
    rc_df = sigma_df*1.5
    rc_dd = sigma_d+R_a
    
    theta = BrAng
    phi = np.random.rand()*(2*np.pi)
    print('Phi: '+str(phi)+'\n Theta: '+str(theta)+'\n')
    x1 = 2+sigma*np.cos(theta*(np.pi/180))
    y1 = 0 + sigma*np.sin(theta*(np.pi/180))*np.cos(phi)
    z1 = 0 + sigma*np.sin(theta*(np.pi/180))*np.sin(phi)
    x2 = 2+2*sigma*np.cos(theta*(np.pi/180))
    y2 = 0 + 2*sigma*np.sin(theta*(np.pi/180))*np.cos(phi)
    z2 = 0 + 2*sigma*np.sin(theta*(np.pi/180))*np.sin(phi)
    
    print('ron = %.1f,roff = %.1f, rnuc = %.1f, rGTP =  %.1f, rGDP = %.1f seed = %d'%(ron,roff,rnuc,rGTP,rGDP,seed))
    ################  Write in.local ################
    
    print('Writing input file\n')
    infile = 'in.local'
    fin=open(infile,'w')
    fin.write(
        '''#####################################################
log                 log.txt
units               lj
dimension           3
atom_style          molecular
boundary            p p p
read_data           configuration.txt extra/bond/per/atom 5  extra/special/per/atom 20  extra/angle/per/atom 3\n''')
    fin.write('variable            ron equal '+str(ron)+'                              # growth rate [monomers/s]\n')
    fin.write('variable            roff equal '+str(roff)+'                                # shrink time [monomers/s]\n')
    fin.write('variable            rnuc equal '+str(rnuc)+'                                # nucleation rate [filaments/s]\n')
    fin.write('variable            rbr equal '+str(rbr/numPhis)+'                                # branching rate [filaments/s]\n')
    fin.write('variable            rGTP equal '+str(rGTP)+'                               # switch rate from GTP->GDP [/s]\n')
    fin.write('variable            rGDP equal '+str(rGDP)+'                               # switch rate from GTP->GDP [/s]\n\n')
    fin.write('variable            Kbond equal '+str(kbond)+'                               # bond constant [kT/sigma2]\n')
    fin.write('variable            Kbend equal '+str(kbend)+'                               # bend constant [kT/sigma2]\n')
    fin.write('variable            BranchAng equal '+str(BrAng)+'                               # Total branching angle [degrees]\n')
    fin.write('variable            eps equal '+str(Epsilon)+'                               # interaction strength between MTs \n')
    fin.write('variable            epsB equal '+str(EpsilonB)+'                               # interaction strength of branch junction bead \n')
    fin.write('variable            eps_dd equal '+str(eps_dd)+'                               # interaction strength of mito-mito \n')
    fin.write('variable            eps_df equal '+str(eps_df)+'                               # interaction strength of mito-MT \n')
    fin.write('variable            lamda equal '+str(lamda)+'                               # how soft is the repulsion of mito-MT (interesting to interrogate pushing) \n')
    fin.write('variable            sigma equal '+str(sigma)+'                               # diameter of filament beads \n')
    fin.write('variable            sigma_dd equal '+str(sigma_d)+'                               # diameter of mito \n')
    fin.write('variable            sigma_df equal '+str(sigma_df)+'                               # resulting combo MT-mito diameter\n')
    fin.write('variable            rc_dd equal '+str(rc_dd)+'                               # radius cutoff  mito-mito \n')
    fin.write('variable            rc_df equal '+str(rc_df)+'                               # radius cutoff  mito-MT  \n')
    fin.write('variable            rc equal '+str(rc)+'                               #  radius cutoff  MT-MT \n')
    fin.write('variable            tstep equal  '+str(tstep)+'                                 # simulation timestep size [s]\n')
    fin.write('variable            run_time equal '+str(run_time)+'                            # simulation run time post MT nucleation [s]\n')
    fin.write('variable            run_time1 equal '+str(run_time1)+'                            # simulation run time pre MT nucleation (only mito) [s]\n')
    fin.write('variable            frame_rate equal '+str(frame_rate)+'                          # dumping interval [s]\n')
    fin.write('variable            OvaR1 equal '+str(OvaR1)+'                          # How much gap do I need before nucleating/polymerising a new MT bead (eg = sigma_df for perfect ratchet)\n')
    fin.write('variable            OvaR2 equal '+str(OvaR2)+'                          # An allowed second value for this (when introducing different sized mito) \n')
    fin.write('variable            factor equal '+str(factor)+'                          # Modulation of rates \n')
    fin.write('variable            seed equal '+str(seed)+'                                  # random number generator seed\n')
    
    fin.write(
        '''variable            angleImp equal ${BranchAng}/10.0                          # growth probability
variable            kon atom ${ron}/${factor}                       # growth probability
variable            knuc atom ${rnuc}/1.0                        # nucleation probability
variable            koff atom ${roff}/${factor}                        # shrink probability 
variable            sGTP atom ${rGTP}/${factor}                        # switch rate from GTP->GDP
variable            sGDP atom ${rGDP}/${factor}                        # switch rate from GTP->GDP
variable            kbr atom ${rbr}/${factor}                         # branching probability
variable            rstep equal 0.1/${tstep}                      # reaction interval [simulation steps]
variable            run_steps equal ${run_time}/${tstep}          # simulation run time [simulation steps]
variable            run_steps1 equal ${run_time1}/${tstep}          # simulation run time [simulation steps]
variable            dump_time equal ${frame_rate}/${tstep}        # dumping interval [simulation steps]
variable 		    yCoord atom y
variable 		    yEdge atom yhi
variable 		    zCoord atom z
variable 		    zEdge atom zhi
variable 		    xCoord atom x
variable 		    xEdge atom xhi
variable 		    ID atom id
variable            stab_steps equal 1

################################################################################
###--------------------- GROUPS --------------------- ###
################################################################################

group               ghosts type 4
group               NucEnd type 2
group               Mito type 7
group               MTs type 1 2 3 5 6

################################################################################
###--------------------- BONDS AND ANGLES --------------------- ###
################################################################################

special_bonds       lj 1.0 1.0 1.0
bond_style          harmonic
bond_coeff          1 ${Kbond} 1.0

angle_style         harmonic
angle_coeff         1 ${Kbend} 180.0
angle_coeff	        2 ${Kbend} ${BranchAng}

################################################################################
###--------------------- PAIR STYLES AND COEFFICIENTS --------------------- ###
################################################################################


pair_style          hybrid/overlay lj/cut/soft 1 0.5 ${rc_df} cosine/squared ${rc_df}
pair_coeff          * * cosine/squared 0.00 2.00 ${rc_df}
pair_coeff          * * lj/cut/soft ${eps} ${sigma} ${lamda} ${rc_df}
pair_coeff          1 1 cosine/squared ${eps} ${sigma} ${rc} wca
pair_coeff          1 2 cosine/squared ${eps} ${sigma} ${rc} wca
pair_coeff          1 3 cosine/squared ${eps} ${sigma} ${rc} wca
pair_coeff          2 2 cosine/squared ${eps} ${sigma} ${rc} wca
pair_coeff          2 3 cosine/squared ${eps} ${sigma} ${rc} wca
pair_coeff          3 3 cosine/squared ${eps} ${sigma} ${rc} wca
pair_coeff          5 1 cosine/squared ${eps} ${sigma} ${rc} wca
pair_coeff          5 2 cosine/squared ${eps} ${sigma} ${rc} wca
pair_coeff          5 3 cosine/squared ${eps} ${sigma} ${rc} wca
pair_coeff          5 5 cosine/squared ${eps} ${sigma} ${rc} wca
pair_coeff          6 1 cosine/squared ${epsB} ${sigma} ${rc} wca
pair_coeff          6 2 cosine/squared ${epsB} ${sigma} ${rc} wca
pair_coeff          6 3 cosine/squared ${epsB} ${sigma} ${rc} wca
pair_coeff          6 5 cosine/squared ${epsB} ${sigma} ${rc} wca
pair_coeff          6 6 cosine/squared ${epsB} ${sigma} ${rc} wca
pair_coeff          7 6 cosine/squared ${eps_df} ${sigma_df} ${rc_df} wca
#pair_coeff          7 1 cosine/squared ${eps_df} ${sigma_df} ${rc_df} wca
#pair_coeff          7 2 cosine/squared ${eps_df} ${sigma_df} ${rc_df} wca
#pair_coeff          7 3 cosine/squared ${eps_df} ${sigma_df} ${rc_df} wca
#pair_coeff          7 5 cosine/squared ${eps_df} ${sigma_df} ${rc_df} wca
pair_coeff          7 1 lj/cut/soft ${eps_df} ${sigma_df} ${lamda} ${rc_df} #Mitochondria have soft core repulsion with the MTs
pair_coeff          7 2 lj/cut/soft ${eps_df} ${sigma_df} ${lamda} ${rc_df}
pair_coeff          7 3 lj/cut/soft ${eps_df} ${sigma_df} ${lamda} ${rc_df}
pair_coeff          7 5 lj/cut/soft ${eps_df} ${sigma_df} ${lamda} ${rc_df}
pair_coeff          7 7 cosine/squared ${eps_dd} ${sigma_dd} ${rc_dd} wca

################################################################################
###--------------------- REACTION MAPS --------------------- ###
################################################################################

# Nucleation reactions molecular templates
molecule            mPreNucleation Reactions/Nuc/pre_Nucleation.txt
molecule            mPostNucleation Reactions/Nuc/post_Nucleation.txt

# Death reactions molecular templates
molecule            mPreDeath Reactions/Nuc/pre_Death.txt
molecule            mPostDeath Reactions/Nuc/post_Death.txt

# Poly reactions molecular templates
molecule            mPreDimerOn Reactions/Poly/pre_DimerOn.txt
molecule            mPostDimerOn Reactions/Poly/post_DimerOn.txt
molecule            mPreTrimerOn Reactions/Poly/pre_TrimerOn.txt
molecule            mPostTrimerOn Reactions/Poly/post_TrimerOn.txt
molecule            mPreOligomerOn Reactions/Poly/pre_OligomerOn.txt
molecule            mPostOligomerOn Reactions/Poly/post_OligomerOn.txt

# Depol reactions molecular templates
molecule            mPreDimerOff Reactions/Depoly/pre_DimerOff.txt
molecule            mPostDimerOff Reactions/Depoly/post_DimerOff.txt
molecule            mPreTrimerOff Reactions/Depoly/pre_TrimerOff.txt
molecule            mPostTrimerOff Reactions/Depoly/post_TrimerOff.txt
molecule            mPreOligomerOff Reactions/Depoly/pre_OligomerOff.txt
molecule            mPostOligomerOff Reactions/Depoly/post_OligomerOff.txt

# Tip switch reactions molecular templates
molecule            mPreDimerGDP Reactions/Tip/pre_DimerTipGDP.txt
molecule            mPreDimerGTP Reactions/Tip/pre_DimerTipGTP.txt
molecule            mPreTrimerGDP Reactions/Tip/pre_TrimerTipGDP.txt
molecule            mPreTrimerGTP Reactions/Tip/pre_TrimerTipGTP.txt
molecule            mPreOligomerGDP Reactions/Tip/pre_OligomerTipGDP.txt
molecule            mPreOligomerGTP Reactions/Tip/pre_OligomerTipGTP.txt

# Branching reactions molecular templates
''')
    for i in range(numPhis):
        fin.write('molecule            mPostBranch1_'+str(i)+' Reactions/Branch/post_branch1_'+str(i)+'.txt\n')
        
    fin.write('''
molecule            mPreBranch1 Reactions/Branch/pre_branch1.txt
molecule            mPreBranch3 Reactions/Branch/pre_branch3.txt
molecule            mPostBranch3 Reactions/Branch/post_branch3.txt
molecule            mPreBranch4 Reactions/Branch/pre_branch4.txt
molecule            mPostBranch4 Reactions/Branch/post_branch4.txt
molecule            mPreBranch3Off_v2 Reactions/Branch/pre_branchDepol3_v2.txt
molecule            mPostBranch3Off_v2 Reactions/Branch/post_branchDepol3_v2.txt
molecule            mPreBranch3Off Reactions/Branch/pre_branchDepol3.txt
molecule            mPostBranch3Off Reactions/Branch/post_branchDepol3.txt
molecule            mPreBranch4Off Reactions/Branch/pre_branchDepol4.txt
molecule            mPostBranch4Off Reactions/Branch/post_branchDepol4.txt
molecule            mPreBranch4Off_v2 Reactions/Branch/pre_branchDepol4_v2.txt
molecule            mPostBranch4Off_v2 Reactions/Branch/post_branchDepol4_v2.txt

#Cleaning Up Branch Junctions
molecule            mPreBranchAng Reactions/Depoly/pre_BranchAngle.txt
molecule            mPostBranchAng Reactions/Depoly/post_BranchAngle.txt

################################################################################
###------------ CLEANING BRANCHING JUNCTIONS WHEN STARTING FROM RESTART ------------ ###
################################################################################

fix             freact all bond/react stabilization yes AllAtoms 0.1 reset_mol_ids no &
                react BranchTidy all ${rstep} 0.900000 1.100000 mPreBranchAng mPostBranchAng Reactions/Maps/map_BranchAngle.txt prob 1.0 ${seed} stabilize_steps ${stab_steps}
''')
        
    fin.write('''


################################################################################
###------------ INTEGRATION ------------ ###
################################################################################

fix                 fLang all langevin 1.0 1.0 0.1 ${seed}
fix                 fNVE AllAtoms_REACT nve

compute             displacement Mito displace/atom
variable            realtime equal step*${tstep}

################################################################################
###------------ OUTPUTS ------------ ###
################################################################################


dump                1 all custom ${dump_time} output.xyz id mol type x y z 
dump_modify         1 format line "%d %d %d %.2f %.2f %.2f"

#dump                3 Mito custom ${dump_time} outputMito.xyz id mol type x y z vx vy vz c_displacement[1] c_displacement[2] c_displacement[3] c_displacement[4]
#dump_modify         3 format line "%d %d %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f"

thermo              ${dump_time}
compute_modify      thermo_temp dynamic/dof yes
thermo_style        custom step v_realtime temp pe ke etotal epair ebond eangle press vol density atoms

compute             cBonds all property/local batom1 batom2 btype

dump                2 all local ${dump_time} bonds.dump c_cBonds[*]
dump_modify         2 format line "%f %f %f"

timestep            ${tstep}

################################################################################
###------------ INITIALISATION, ONLY MITOCHONDRIA ------------ ###
################################################################################

run                 ${run_steps1}

write_restart	restart.postwall

unfix           fLang
unfix           fNVE
unfix           freact

################################################################################
###------------ ALLOW NUCLEATION AND GROWTH OF MICROTUBULES ------------ ###
################################################################################


fix             freact all bond/react stabilization yes AllAtoms 0.1 reset_mol_ids no &
                react BranchTidy all ${rstep} 0.900000 1.100000 mPreBranchAng mPostBranchAng Reactions/Maps/map_BranchAngle.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} &
                react Nucleation all ${rstep} 0.900000 1.100000 mPreNucleation mPostNucleation Reactions/Maps/map_Nucleation.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} modify_create overlap ${OvaR1} ${OvaR2} modify_create nuc xwalls -2 &
                react Death all ${rstep} 0.900000 1.100000 mPreDeath mPostDeath Reactions/Maps/map_Death.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} &
                react DimerOn all ${rstep} 0.900000 1.100000 mPreDimerOn mPostDimerOn Reactions/Maps/map_DimerOn.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} modify_create fit 1 modify_create overlap ${OvaR1} ${OvaR2}  &
                react TrimerOn all ${rstep} 0.900000 1.100000 mPreTrimerOn mPostTrimerOn Reactions/Maps/map_TrimerOn.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} modify_create fit 1 modify_create overlap ${OvaR1} ${OvaR2} &
                react OligomerOn all ${rstep} 0.900000 1.100000 mPreOligomerOn mPostOligomerOn Reactions/Maps/map_OligomerOn.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} modify_create fit 1 modify_create overlap ${OvaR1} ${OvaR2} &
                react OligomerOff all ${rstep} 0.900000 1.100000 mPreOligomerOff mPostOligomerOff Reactions/Maps/map_OligomerOff.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} &
                react TrimerOff all ${rstep} 0.900000 1.100000 mPreTrimerOff mPostTrimerOff Reactions/Maps/map_TrimerOff.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} &
                react DimerOff all ${rstep} 0.900000 1.100000 mPreDimerOff mPostDimerOff Reactions/Maps/map_DimerOff.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} &
                react tipGTP all ${rstep} 0.900000 1.100000 mPreOligomerGDP mPreOligomerGTP Reactions/Maps/map_TipGTP.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} &
                react tipGDP all ${rstep} 0.900000 1.100000 mPreOligomerGTP mPreOligomerGDP Reactions/Maps/map_TipGDP.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} &
                react tipTrimGTP all ${rstep} 0.900000 1.100000 mPreTrimerGDP mPreTrimerGTP Reactions/Maps/map_TipGTPtrim.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} &
                react tipTrimGDP all ${rstep} 0.900000 1.100000 mPreTrimerGTP mPreTrimerGDP Reactions/Maps/map_TipGDPtrim.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} &
                react tipDimGTP all ${rstep} 0.900000 1.100000 mPreDimerGDP mPreDimerGTP Reactions/Maps/map_TipGTP.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} &
                react tipDimGDP all ${rstep} 0.900000 1.100000 mPreDimerGTP mPreDimerGDP Reactions/Maps/map_TipGDP.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} &
                react Branch3 all ${rstep} 0.900000 1.100000 mPreBranch3 mPostBranch3 Reactions/Maps/map_Branch3.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} modify_create fit 3 modify_create overlap ${OvaR1} ${OvaR2} &
                react Branch4 all ${rstep} 0.900000 1.100000 mPreBranch4 mPostBranch4 Reactions/Maps/map_OligomerOn.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} modify_create fit 3 modify_create overlap ${OvaR1} ${OvaR2} &
                react DepBranch3 all ${rstep} 0.900000 1.100000 mPreBranch3Off mPostBranch3Off Reactions/Maps/map_DepolBranch3.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} &
                react DepBranch3_v2 all ${rstep} 0.900000 1.100000 mPreBranch3Off_v2 mPostBranch3Off_v2 Reactions/Maps/map_DepolBranch3.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} &
                react DepBranch4 all ${rstep} 0.900000 1.100000 mPreBranch4Off mPostBranch4Off Reactions/Maps/map_DepolBranch4.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} &
                react DepBranch4_v2 all ${rstep} 0.900000 1.100000 mPreBranch4Off_v2 mPostBranch4Off_v2 Reactions/Maps/map_DepolBranch4.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} &
''')
    for i in range(numPhis):
        fin.write('\t\t\t\treact Branch_'+str(i)+' all ${rstep} 0.900000 1.100000 mPreBranch1 mPostBranch1_'+str(i)+' Reactions/Maps/map_Branch1.txt prob 0.1 ${seed} stabilize_steps ${stab_steps} modify_create fit 3 modify_create overlap ${OvaR1} ${OvaR2} &\n')

    fin.write('''

################################################################################
###------------ GROUPS AND COUNTING ------------ ###
################################################################################

variable            vMaskGTPType atom "type==3"
group               FilGTP dynamic all var vMaskGTPType every 1
variable            vNumGTP equal count(FilGTP)
variable            vMaskGDPType atom "type==5"
group               FilGDP dynamic all var vMaskGDPType every 1
variable            vNumGDP equal count(FilGDP)

variable            vMaskAlive atom "type==1 || type==2 || type==3 || type==5"
group               alive dynamic all var vMaskAlive every 1
variable            vNumAlive equal count(alive)

variable            vNumNuc equal f_freact[1]
variable            vDeath equal f_freact[2]
variable            vNumPoly equal f_freact[5]
variable            vNumDepoly equal f_freact[6]


fix                 fLang all langevin 1.0 1.0 0.1 ${seed}
fix                 fNVE AllAtoms_REACT nve

fix                 print_thermo all ave/time ${dump_time} 1 ${dump_time} v_vNumGTP v_vNumGDP v_vNumAlive v_vNumPoly v_vNumDepoly v_vNumNuc v_vDeath file thermo.txt

run ${run_steps}


''')
    
    
    fin.close
    
    print('Writing submission file\n')
    runscriptname = 'runscript.sh'
    fsub = open(runscriptname, "w")
    fsub.write(
        '''#!/bin/bash
#SBATCH --job-name=Job_ReactionModel
#SBATCH --output=OUTFILE_reaction.dat
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --mail-user=bmeadowc@ist.ac.at
#SBATCH --mail-type=NONE
#SBATCH --no-requeue
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV
export OMP_NUM_THREADS=1

module load gcc
module load openmpi
module load lammps/20220623b\n


mpirun -np 1 /nfs/scistore15/saricgrp/bmeadowc/Scratch/lammps-15Jun2023/src/lmp_mpi -in in.local
\n''')   
    fsub.close()
    for i in range(numPhis):
        theta = 90 #BrAng
        phi = (2*np.pi)*(i/numPhis)
        #print('Phi: '+str(phi)+'\n Theta: '+str(theta)+'\n')
        x1 = 2+sigma*np.cos(theta*(np.pi/180))
        y1 = 0 + sigma*np.sin(theta*(np.pi/180))*np.cos(phi)
        z1 = 0 + sigma*np.sin(theta*(np.pi/180))*np.sin(phi)
        x2 = 2+2*sigma*np.cos(theta*(np.pi/180))
        y2 = 0 + 2*sigma*np.sin(theta*(np.pi/180))*np.cos(phi)
        z2 = 0 + 2*sigma*np.sin(theta*(np.pi/180))*np.sin(phi)
        topologyname = 'post_branch1_'+str(i)+'.txt'
        ftop = open(topologyname, "w")
        ftop.write(
            '''molecule template: end of chain plus polymerized styrene
    
     5 atoms
     4 bonds
     3 angles
     3 fragments

Coords

    1 1.0 0.0 0.0
    2 2.0 0.0 0.0
    3 3.0 0.0 0.0
''')
    
        ftop.write('\t4 '+str(np.round(x1,2))+' '+str(np.round(y1,2))+' '+str(np.round(z1,2))+'\n')
        ftop.write('\t5 '+str(np.round(x2,2))+' '+str(np.round(y2,2))+' '+str(np.round(z2,2))+'\n')
        ftop.write(
        '''
Types

    1 1
    2 1 
    3 1
    4 6
    5 3

Molecules

    1 0
    2 0 
    3 0 
    4 0
    5 0

Bonds

    1 1 1 2
    2 1 2 3
    3 1 2 4
    4 1 4 5

Angles

    1 1 1 2 3
    2 2 3 2 4
    3 1 2 4 5

Fragments

    1 1
    2 3
    3 2 3 4''') 
    
        ftop.close()
    
    
if __name__ == "__main__":
    main()
#module load openmpi/4.1.1
#module load pmix
#module load lammps