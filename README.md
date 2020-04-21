# Solvent_sites_clustering
This repo contains the information on how to obtain solvent sites according to https://watclust.wordpress.com/methodology/
using as input molecular dynamics in AMBER format (https://ambermd.org/).

It requieres: A molecular dynamics trajectory (*.netcdf) and a topology file (*.prmtop)

To obtain the solvent sites the command is:

python scipy_ss2.py -nc traj.nc -p top.prmtop -bs :1-20,23,25,28-40 -dist 0.23 -outdir eta_C1_sites -r ref.pdb -solv ETA -probe C1 -ff 0 -lf 999 -watNmin 80

-bs is the binding site using AMBER selection sintax. 

-dist is the distance to join the probe atoms.

-ff is the first frame of the molecular dynamic we want to analyze.

-lf is the last frame.

-watNmin is the minimum number of waters that a cluster can have.

-r is the reference PDB that will be used to align the trajectory.

The default density corresponds to the ethanol density. If you want to get the water sites you should use -dens 0.0334 as argument.
