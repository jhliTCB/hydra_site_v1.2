#include "netcdf.h"
#include "AmberNetcdf.h"
#include "global.h"
/*
int proc_pdb(char * filename, int * atomnum, int * watnum, ATOM * atom);
int proc_xtc(char * filename, int atomnum, int watnum, ATOM * atom, awt * mywat, vt * p_grid_ct, vt * max_wat_edge, int * frames, int * len_awt);
int proc_nc(ATOM * atom, awt * mywat, struct AmberNetcdf * A, double * X, vt * p_grid_ct, vt * max_wat_edge, int atomnum, int * tot_wat);
int checkout(int len_awt, awt * mywat, vt * boxct);
int findwat(int tot_wat, awt * mywat, int frames, float (*init_grid)[5], int * init_count, float box_center[3], float box_size[3], int flag_nc, float rangeWAT)
int merge(float (*init_grid)[5], int init_count, float (*merge_grid)[5], int * merge_count);
int iteration(awt * mywat, int tot_wat, float (*init_grid)[5], int init_count, int frames, int flag_nc, float rangeWAT);
*/

extern int proc_pdb(char *, int *, int *, ATOM *);
extern int proc_xtc(char *, int, int, ATOM *, awt *, vt *, vt *, int *, int *);
extern int proc_nc(ATOM *, awt *, struct AmberNetcdf *, double *, vt *, vt *, int, int *);
extern int findwat(int, awt *, int, float (*)[5], int *, vt, vt, int, float);
extern int merge(float (*)[5], int, float (*)[5], int *);
extern int iteration(awt *, int, float (*)[5], int, int, int, float);
