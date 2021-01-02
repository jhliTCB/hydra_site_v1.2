# hydra_site_v1.2

This is an update from a previous project in a summer school: https://gitlab.com/jhli/PDC_summer_school

This update is just for fun.

More complicated functions could be achived by using mobywat (http://www.mobywat.com/index.php)

Changes to the summer school version:

    1). added support for Amber netcdf trajectories
    
    2). smaller pool (only grid around the protein)
    
    3). threshold for counting can be input from argv
    
    4). remove omp and cuda parallelizations
    
    5). add a flag to the structure ATOM

To compile the code, libnetcdf is require, can use the one from amber package, e.g.

export LIBRARY_PATH=/data/soft/amber18/lib

make
