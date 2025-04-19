from __future__ import division, print_function
import numpy as np
import math
import time
import random
import sys, getopt
import os

def placeIn(xbox1,ybox1,zbox1,sigma_d,xarray,yarray,zarray):
    c=1
    while c==1:
        x = np.random.rand()*(xbox1*2)-xbox1+20 #Add some amount to make off-centre
        y = np.random.rand()*(ybox1*2)-ybox1
        z = np.random.rand()*(zbox1*2)-zbox1
        c=0
        for i in range(len(xarray)):
            Dist = np.sqrt((x-xarray[i])**2 + (y-yarray[i])**2+(z-zarray[i])**2)
            if Dist <= sigma_d:
                c = 1
        if c==0:
            Dist2 = (x-0)**2+(y-0.5)**2+(z+2)**2
            Dist3 = (x-0)**2+(y+0.5)**2+(z+2)**2
            if Dist2 <= (sigma_d+1)/2 or Dist3 <= (sigma_d+1)/2:
                c=1
    return x,y,z
        
def main():
    ################ Individual filament parameters ################
    Ndiff = int(sys.argv[1])
    sigma_d = float(sys.argv[2])
    mass_d = float(sys.argv[3])
    xbox = 100 #100
    ybox = 25 #20
    zbox = ybox 
    xbox1,ybox1,zbox1 = 20*(Ndiff/225)-sigma_d*0.6,ybox-sigma_d*0.1,zbox-sigma_d*0.1
    #Nside = int(np.cbrt(Ndiff))
    molid = 0
    Natoms = 2+Ndiff
    xarray = []
    yarray = []
    zarray = []
    ################  Write in.local ################
    
    print('Writing configuration file\n')
    infile = 'configuration.txt'
    fin=open(infile,'w')
    fin.write('Configuration file for spindle and diffusors\n')
    fin.write(str(Natoms)+' atoms\n')
    fin.write('''1 bonds
0 angles
7 atom types
1 bond types
2 angle types

''')
    fin.write(str(-xbox)+' '+str(xbox)+' xlo xhi\n')
    fin.write(str(-ybox)+' '+str(ybox)+' ylo yhi\n')
    fin.write(str(-zbox)+' '+str(zbox)+' zlo zhi\n')
    fin.write('''
Masses

1 1.0
2 100000
3 1.0
4 1.0
5 1.0
6 1.0
''')
    fin.write('7 '+str(mass_d)+'\n')
    fin.write('''

Atoms

1 0 4 0.0 -0.5 -2.0
2 0 4 0.0 0.5 -2.0
''')
    for n in range(Ndiff):
        molid +=1
        #print("here", n)
        x,y,z = placeIn(xbox1,ybox1,zbox1,sigma_d,xarray,yarray,zarray) #-xbox1+nx*xbox1*2/Nside,-ybox1+ny*ybox1*2/Nside,-zbox1+nz*zbox1*2/Nside
        xarray.append(x)
        yarray.append(y)
        zarray.append(z)

        fin.write(str(molid+2)+' '+str(molid)+' 7 '+str(x)+' '+str(y)+' '+str(z)+'\n')

    fin.write('''
    
Bonds

1 1 1 2''')
    
    
    fin.close
    
    
if __name__ == "__main__":
    main()
