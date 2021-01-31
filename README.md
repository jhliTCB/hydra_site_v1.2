# hydra_site_v1.2

This is an update from a previous project in a summer school: https://gitlab.com/jhli/PDC_summer_school

This update is just for fun.

More complicated functions could be achived by using mobywat (http://www.mobywat.com/index.php)

Changes to the summer school version:

    1). added support for Amber netcdf trajectories
    
    2). smaller pool (only grid around the protein)
    
    3). threshold for counting can be input from argv (default 20%)
    
    4). remove omp and cuda parallelizations
    
    5). add a flag to the structure ATOM

To compile the code, libnetcdf is require, can use the one from amber package, e.g.

export LIBRARY_PATH=/data/soft/amber18/lib

make

USAGE:

./hydra_site.x pdb_file trajectory(*.xtc/*.nc)

./hydra_site.x pdb_file trajectory 0.3


Requirement for the xtc/nc file: PBC removed and structure aligned.

Gromacs:

gmx_mpi make_ndx -f md.pdb -o protein.ndx >& /dev/null <<_EOF

a 1-3000

q

_EOF

gmx_mpi trjconv -s md.tpr -f md.xtc -pbc mol -ur compact -n index.ndx -center -o md_np1.xtc # Choose the group a1-3000 for center, system for output

gmx_mpi trjconv -s md.tpr -f md_np1.xtc -o md_np1_fit.xtc -ur compact -fit rot+trans



Amber:

cpptraj -p p_sol.prmtop < remove_pbc.in

remove_pbc.in:

trajin md.nc

autoimage # or image familiar

reference referece.restrt

trajout md_fit.nc

go





