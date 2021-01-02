
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "main.h"

int main(int argc, char *argv[])
{
   if (argc < 3) {
      printf("Error, please input the topology(*.pdb) and trajectory(*.nc/*.xtc)\n");
      exit(1);
   }
   char * trj_ext;
   int atomnum, watnum, frames, les_num, init_count, flag_nc, flag_xtc;
   int tot_wat;
   float t_star, t_end1, t_end2;
   float init_grid[5000][5] = {0};
   float rangeWAT = 0.2; // default threshold for keeping water molecules
   vt  num_grid; 
   vt  min_grid_edge;
   ATOM * atom;
   awt  * mywat;
   size_t atom_size  = sizeof(ATOM) * MAXATOM;
   size_t mywat_size = sizeof(awt) * MAXWAT;
   atom   = (ATOM *) malloc(atom_size);
   mywat  = (awt *) malloc(mywat_size);

   proc_pdb(argv[1], &atomnum, &watnum, atom);
   printf("mywat_size is %u, atom_size is %u\n", mywat_size, atom_size);

   if (DEBUG == 1) {
      printf("\natoms %d; water % d\n", atomnum, watnum);
      printf("omp_get_num_threads() is %d\n", omp_get_num_threads());
   }

   t_star = omp_get_wtime();

   trj_ext  = strrchr(argv[2], '.');
   flag_nc  = strcmp(trj_ext+1, "nc");
   flag_xtc = strcmp(trj_ext+1, "xtc");
   if (flag_nc == 0) {
      printf("Reading %s as amber netcdf trajectory\n", argv[2]);
      struct AmberNetcdf * A;
      double * X;
      A = (struct AmberNetcdf*) malloc( MAXFRAME * sizeof(struct AmberNetcdf) );
      netcdfLoad(A, argv[2]);
      frames = A->ncframe;
      X      = (double *) malloc(sizeof(double) * atom_size);            
      proc_nc(atom, mywat, A, X, &num_grid, &min_grid_edge, atomnum, &tot_wat);
   }
   else if (flag_xtc == 0) {
      printf("Reading %s as gromacs xdr trajectory\n", argv[2]);
      proc_xtc(argv[2], atomnum, watnum, atom, mywat, &num_grid, &min_grid_edge, &frames, &tot_wat);
   }
   else {
      printf("Only supporting Amber netcdf (*.nc) and Gromacs xdrfile (*.xtc) formats of trajectories\n");
      exit(1);
   }
   
   t_end1 = omp_get_wtime();
   printf("\nTime duration for input_reading: %f s\n", t_end1 - t_star);

   // Here the griding box is defined! Default is within the range of 0.3 nm of all protein atoms

   if (argv[3]) {rangeWAT =  strtof(argv[3],NULL);}
   //int findwat(int tot_wat, awt * mywat, int frames, float (*init_grid)[5], int * init_count, float box_center[3], float box_size[3], int flag_nc, float rangeWAT)
   findwat(tot_wat, mywat, frames, init_grid, &init_count, num_grid, min_grid_edge, flag_nc, rangeWAT);
   t_end2 = omp_get_wtime();
   printf("\nTime duration for findwat(): %f s\n", t_end2 - t_end1);
   
   iteration(mywat, tot_wat, init_grid, init_count, frames, flag_nc, rangeWAT);
   t_end1 = omp_get_wtime();
   printf("\nTime duration for iteration(): %f s\n", t_end1 - t_end2);
   printf("\nTOTAL running Time duration: %f s (%f min)\n", t_end1 - t_star, (t_end1 - t_star)/60.0);
   
   free(atom);
   free(mywat);

   atom = NULL;
   mywat = NULL;

   return 0;
}