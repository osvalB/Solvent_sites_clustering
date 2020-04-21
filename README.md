# Solvent_sites_clustering
This repo contains a python script to obtain solvent sites according to https://watclust.wordpress.com/methodology/
using as input molecular dynamics in AMBER format (https://ambermd.org/).

It requieres: A molecular dynamics trajectory file (*.netcdf) and a topology file (*.prmtop).
The python libraries pytraj, scipy and networkx.

To obtain the solvent sites run the following command:

python scipy_ss2.py -nc traj.nc -p top.prmtop -bs :1-20,23,25,28-40 -dist 0.23 -outdir eta_C1_sites -r ref.pdb -solv ETA -probe C1 -ff 0 -lf 999 -watNmin 80

-bs is the binding site using AMBER selection sintax. 

-dist is the distance to join the probe atoms.

-ff is the first frame of the molecular dynamic we want to analyze.

-lf is the last frame.

-watNmin is the minimum number of waters that a cluster can have.

-r is the reference PDB that will be used to align the trajectory.

The default density corresponds to the ethanol density. If you want to get the water sites you should use -dens 0.0334 as argument. The unit of density are #molecules / (A^3).

For more details about co-solvent sites please read:

 Arcon, J. P., Defelipe, L. A., Lopez, E. D., Burastero, O., Modenutti, C. P., Barril, X., ... & Turjanski, A. G. (2019). Cosolvent-based protein pharmacophore for ligand enrichment in virtual screening. Journal of chemical information and modeling, 59(8), 3572-3583.
 
 Arcon, J. P., Defelipe, L. A., Modenutti, C. P., López, E. D., Alvarez-Garcia, D., Barril, X., ... & Martí, M. A. (2017). Molecular dynamics in mixed solvents reveals protein–ligand interactions, improves docking, and allows accurate binding free energy predictions. Journal of chemical information and modeling, 57(4), 846-863.
