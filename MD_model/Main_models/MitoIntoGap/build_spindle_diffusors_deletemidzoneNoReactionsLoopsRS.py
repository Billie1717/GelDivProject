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
    ron = float(sys.argv[1])
    roff = float(sys.argv[2])
    rnuc = float(sys.argv[3])
    rGTP = float(sys.argv[4])
    rGDP = float(sys.argv[5])
    rbr = float(sys.argv[6])
    BrAng= float(sys.argv[7])
    numPhis = int(sys.argv[8])
    sigma_d = float(sys.argv[9])
    R_a = float(sys.argv[10])
    eps_dd = float(sys.argv[11])
    eps_df = float(sys.argv[12])
    lamda = float(sys.argv[13])
    OvaR1 = float(sys.argv[14])
    OvaR2 = float(sys.argv[15])
    factor = float(sys.argv[16])
    Pdist = float(sys.argv[17])
    mass_d = float(sys.argv[18])
    seed = float(sys.argv[19])
    
    Epsilon = 0.5 
    EpsilonB = 0.0
    kbond = 700
    kbend = 1000
    tstep = 0.0005 #1/2000
    run_time = 2000.0 #4000.0
    frame_rate = 10.0 #1.0 #0 #
    frame_rate2 = 50.0 #0 #
    TstepsofLoop =frame_rate/tstep
    numloops = int(run_time/frame_rate)

    ##Calculating branch placement
    sigma = 1.0
    sigma_df = (sigma +sigma_d)/2
    # VOLUME EXCLUSION ONLY
    rc = sigma
    rc_df = sigma_df #*1.2
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
    infile = 'inDEL.local'
    fin=open(infile,'w')
    fin.write(
        '''#####################################################
log                 logDEL.txt
units               lj
dimension           3
atom_style          molecular
boundary            p p p
#read_restart        spindle.restart.6600000
read_restart        restart.postwall\n
''')
    fin.write('mass                7 '+str(mass_d)+' #mass of mitos \n')
    fin.write('variable            ron equal '+str(ron)+'                              # growth rate [monomers/s]\n')
    fin.write('variable            roff equal '+str(roff)+'                                # shrink time [monomers/s]\n')
    fin.write('variable            rnuc equal '+str(rnuc)+'                                # nucleation rate [filaments/s]\n')
    fin.write('variable            rbr equal '+str(rbr/numPhis)+'                                # branching rate [filaments/s]\n')
    fin.write('variable            rGTP equal '+str(rGTP)+'                               # switch rate from GTP->GDP [/s]\n')
    fin.write('variable            rGDP equal '+str(rGDP)+'                               # switch rate from GTP->GDP [/s]\n\n')
    fin.write('variable            Kbond equal '+str(kbond)+'                               # bond constant [kT/sigma2]\n')
    fin.write('variable            Kbend equal '+str(kbend)+'                               # bend constant [kT/sigma2]\n')
    fin.write('variable            BranchAng equal '+str(BrAng)+'                               # Total branching angle (will be split into two) [degrees]\n')
    fin.write('variable            eps equal '+str(Epsilon)+'                               # interaction strength \n')
    fin.write('variable            epsB equal '+str(EpsilonB)+'                               # interaction strength of branch junction bead \n')
    fin.write('variable            eps_dd equal '+str(eps_dd)+'                               # interaction strength of diffusors wt dif \n')
    fin.write('variable            eps_df equal '+str(eps_df)+'                               # interaction strength of diffusors wt fil \n')
    fin.write('variable            lamda equal '+str(lamda)+'                               # interaction strength of diffusors wt fil \n')
    fin.write('variable            sigma equal '+str(sigma)+'                               # diameter of filament beads \n')
    fin.write('variable            sigma_dd equal '+str(sigma_d)+'                               # diameter of diffusors \n')
    fin.write('variable            sigma_df equal '+str(sigma_df)+'                               #  \n')
    fin.write('variable            rc_dd equal '+str(rc_dd)+'                               # radius cutoff \n')
    fin.write('variable            rc_df equal '+str(rc_df)+'                               #  \n')
    fin.write('variable            rc equal '+str(rc)+'                               #  \n')
    fin.write('variable            tstep equal  '+str(tstep)+'                                 # simulation timestep size [seconds]\n')
    fin.write('variable            run_time equal '+str(run_time)+'                            # simulation run time [seconds]\n')
    fin.write('variable            tsteploop equal '+str(TstepsofLoop)+'                            # simulation run time [seconds]\n')
    fin.write('variable            frame_rate equal '+str(frame_rate)+'                          # dumping interval [seconds]\n')
    fin.write('variable            frame_rate2 equal '+str(frame_rate2)+'                          # dumping interval [seconds]\n')
    fin.write('variable            OvaR1 equal '+str(OvaR1)+'                          # \n')
    fin.write('variable            OvaR2 equal '+str(OvaR2)+'                          # \n')
    fin.write('variable            pole_distance equal '+str(Pdist)+'                          # \n')
    fin.write('variable            factor equal '+str(factor)+'                          # \n')
    fin.write('variable            seed equal '+str(seed)+'                                  # random number generator seed\n')
    
    fin.write(
        '''variable            angleImp equal ${BranchAng}/10.0                          # growth probability
variable            kon atom ${ron}/${factor}                       # growth probability
variable            knuc atom ${rnuc}/1.0                        # nucleation probability
variable            koff atom ${roff}/${factor}                        # shrink probability 
variable            sGTP atom ${rGTP}/${factor}                        # switch rate from GTP->GDP
variable            sGDP atom ${rGDP}/${factor}                        # switch rate from GTP->GDP
variable            kbr atom ${rbr}/${factor}                         # branching probability
variable            rstep equal 100/${tstep}                      # reaction interval [simulation steps]
variable            run_steps equal ${run_time}/${tstep}          # simulation run time [simulation steps]
variable            dump_time equal ${frame_rate}/${tstep}        # dumping interval [simulation steps]
variable            restart_time equal ${frame_rate2}/${tstep}        # dumping interval [simulation steps]
variable 		    yCoord atom y
variable 		    yEdge atom yhi
variable 		    zCoord atom z
variable 		    zEdge atom zhi
variable 		    xCoord atom x
variable 		    xEdge atom xhi
variable 		    ID atom id
variable            stab_steps equal 1

group               ghosts type 4
group               NucEnd type 2
group               MTs type 1 2 3 5 6
group               Mito type 7

special_bonds       lj 1.0 1.0 1.0
bond_style          harmonic
bond_coeff          1 ${Kbond} 1.0

angle_style         harmonic
angle_coeff         1 ${Kbend} 180.0
angle_coeff	        2 ${Kbend} ${BranchAng}

pair_style          hybrid/overlay lj/cut/soft 1 0.5 ${rc_df} cosine/squared ${rc_df}
pair_coeff          * * cosine/squared 0.00 2.00 ${rc_df}
pair_coeff          * * lj/cut/soft 0.0 ${sigma} ${lamda} ${rc_df}
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
pair_coeff          7 6 cosine/squared ${epsB} ${sigma_df} ${rc_df} wca
pair_coeff          7 1 lj/cut/soft ${eps_df} ${sigma_df} ${lamda} ${rc_df} #Mitochondria have soft core repulsion with the MTs
pair_coeff          7 2 lj/cut/soft ${eps_df} ${sigma_df} ${lamda} ${rc_df}
pair_coeff          7 3 lj/cut/soft ${eps_df} ${sigma_df} ${lamda} ${rc_df}
pair_coeff          7 5 lj/cut/soft ${eps_df} ${sigma_df} ${lamda} ${rc_df}
pair_coeff          7 7 cosine/squared ${eps_dd} ${sigma_dd} ${rc_dd} wca

special_bonds  	lj 0.0 0.0 0.5 

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


variable            realtime equal step*${tstep}


############# GROUPS AND COUNTING ###############

variable            vMaskGTPType atom "type==3"
group               FilGTP dynamic all var vMaskGTPType every 1
variable            vNumGTP equal count(FilGTP)
variable            vMaskGDPType atom "type==5"
group               FilGDP dynamic all var vMaskGDPType every 1
variable            vNumGDP equal count(FilGDP)

variable            vMaskAlive atom "type==1 || type==2 || type==3 || type==5"
group               alive dynamic all var vMaskAlive every 1
variable            vNumAlive equal count(alive)

fix                 fLang all langevin 1.0 1.0 1.0 ${seed}
fix                 fNVE AllAtoms_REACT nve

#-------HISTORICAL DELETING MIDDLE ATOMS --------------#
#fix                 fBreak MtsInMidzoneO bond/break ${run_steps} 1 0.9 prob 1 1234
#delete_bonds        MtsInMidzoneO multi remove
#---------------------------------------------------#


dump                1 all custom ${dump_time} outputDEL.xyz id mol type x y z 
dump_modify         1 format line "%d %d %d %.2f %.2f %.2f"

thermo              ${dump_time}
compute_modify      thermo_temp dynamic/dof yes
thermo_style        custom step v_realtime temp pe ke etotal epair ebond eangle press vol density atoms

compute             cBonds all property/local batom1 batom2 btype

dump                2 all local ${dump_time} bondsDEL.dump c_cBonds[*]
dump_modify         2 format line "%f %f %f"

restart             ${restart_time} spindle.restart

fix                 print_thermo all ave/time ${dump_time} 1 ${dump_time} v_vNumGTP v_vNumGDP v_vNumAlive file thermo3.txt
timestep      ${tstep}

#-------HISTORICAL DELETING MIDDLE ATOMS --------------#
#run 100
#unfix fBreak
#reset_atoms mol all
#delete_atoms        group MtsInMidzoneI mol no bond yes
#run                 ${run_steps}
#fix                 fBreak MtsInMidzone bond/break ${run_steps} 1 0.9 prob 1 1234
#delete_bonds        MtsInMidzone multi remove
#run 100
#unfix fBreak
#reset_atoms mol all
#delete_atoms        group MtsInMidzoneI mol no bond yes 
#reset_atoms mol MTs
#---------------------------------------------------#
#region regionboxO block -9.0 9.0 -50.0 50.0 -50.0 50.0
#group InMidzoneO region regionboxO
#group MtsInMidzoneO intersect MTs InMidzoneO
#group MtsInMidzoneO dynamic MTs region regionboxO every 100
#set group MtsInMidzoneO type 4
run 100
#set group MTs type 4
#fix                 fBreak MTs bond/break ${run_steps} 1 0.9 prob 1 1234
#delete_bonds        MTs multi remove
reset_atoms mol all
delete_atoms group MTs mol yes bond yes

run                 ${run_steps}
write_restart	restart.postrelax

''')
    
    
    fin.close
    
    print('Writing submission file\n')
    runscriptname = 'runscriptDel.sh'
    fsub = open(runscriptname, "w")
    fsub.write(
        '''#!/bin/bash
#SBATCH --job-name=Job_ReactionModel
#SBATCH --output=OUTFILE_reaction.dat
#SBATCH --time=96:00:00
#SBATCH --ntasks=4
#SBATCH --mem=1G
#SBATCH --mail-user=bmeadowc@ist.ac.at
#SBATCH --mail-type=NONE
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
export OMP_NUM_THREADS=1

module load openmpi/4.1.1
module load pmix
module load lammps
module load gcc

mpirun -np 1 /nfs/scistore15/saricgrp/bmeadowc/Scratch/lammps-15Jun2023/src/lmp_mpi -in inDEL.local
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
    4 1
    5 1

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
