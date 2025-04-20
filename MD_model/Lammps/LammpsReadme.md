From the lammps downloads:

https://download.lammps.org/tars/

The version used in this work is: lammps-15Jun2023

With some packages edited, namely the files in the REACTION package (see below)

all simulations were run on Lammps built with the following packages :

make yes-dipole
make yes-extra-pair
make yes-rigid
make yes-molecule
make yes-reaction

# Files edited in folder /src

For specifying the geometry of newly nucleated microtubule beads as well as the different sizes for microtubule beads not to grow against

REACTION/fix_bond_react.cpp
REACTION/fix_bond_react.h

The model only works with these edited files, they have to replace the original files in the lammps src folder and lammps has to be built with these files.
