// A try for integration
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 
#include <memory.h>
#include "global.h"
#include "utils.h"

int findwat(int tot_wat, awt * mywat, int frames, float (*init_grid)[5], int * init_count, vt num_grid, vt min_grid_edge, int flag_nc, float rangeWAT)
{
   int i, j,k, x, y, z;
   int ***count;
   float occ; // Occupancy
   int max_count=0;
   int valid_count = 0;
   float conv = 10.0;
   float fframe = (float) frames;
   int tmp1     = (int) num_grid.v1;
   int tmp2     = (int) num_grid.v2;
   int tmp3     = (int) num_grid.v3;
   int grid_num[3] = {tmp1, tmp2, tmp3};
   FILE * InitWatSites = fopen("WatSites_inits.pdb", "w");
   
   if (flag_nc == 0) {
      conv = 1.0;
   }

   printf("number of grids on xyz: %d %d %d\n", grid_num[0], grid_num[1], grid_num[2]);
   printf("grid_edge_min on xyz: %f %f %f\n", min_grid_edge.v1, min_grid_edge.v2, min_grid_edge.v3);

   // count is 3D array to store the grid
   count = (int***)calloc(grid_num[0],sizeof(int**));
   for(i=0;i<grid_num[0];i++)
   {
     count[i]=(int**)calloc(grid_num[1],sizeof(int*));
     for(j=0;j<grid_num[1];j++)
     {
       count[i][j]=(int*)calloc(grid_num[2],sizeof(int));
       for(k=0;k<grid_num[2];k++)
       {
         count[i][j][k]=0;
       }
     }
   }

   // Get the initial countings:
   for (i = 0; i < tot_wat; i++) {
      x = (int)((mywat[i].x - min_grid_edge.v1)/(GRID_SIZE/conv));
      y = (int)((mywat[i].y - min_grid_edge.v2)/(GRID_SIZE/conv));
      z = (int)((mywat[i].z - min_grid_edge.v3)/(GRID_SIZE/conv));
      if (x < grid_num[0] && y < grid_num[1] && z < grid_num[2] && x>=0 && y>=0 && z>=0) {
         count[x][y][z] += 1;
      }
      else if (DEBUG == 2) {
         printf("aborted count x, y, z: %d %d %d; at i %d, mywat[i]'s xyz: %f %f %f; mywat[i].fid %d\n", x,y,z,i,mywat[i].x,mywat[i].y,mywat[i].z,mywat[i].fid);
      }
   }
   
   for (x = 0; x < grid_num[0]; x++) {
      for (y = 0; y < grid_num[1]; y++) {
         for (z = 0; z < grid_num[2]; z++) {
            if (count[x][y][z]>max_count) {
               max_count=count[x][y][z];
               // To recover the center's xyz from edge, the shift should be added back
               printf("max count: %d in grid center: %f %f %f\n", max_count, x * GRID_SIZE/conv  + min_grid_edge.v1 + GRID_SIZE/conv / 2.0, y * GRID_SIZE/conv + min_grid_edge.v2 + GRID_SIZE/conv / 2.0, z * GRID_SIZE/conv + min_grid_edge.v3 + GRID_SIZE/conv / 2.0);
            }
            float fcount = (float) count[x][y][z];
	    
            if (fcount > fframe * rangeWAT) {
               init_grid[valid_count][0] = min_grid_edge.v1 + x * GRID_SIZE/conv + GRID_SIZE/conv / 2.0;
               init_grid[valid_count][1] = min_grid_edge.v2 + y * GRID_SIZE/conv + GRID_SIZE/conv / 2.0;
               init_grid[valid_count][2] = min_grid_edge.v3 + z * GRID_SIZE/conv + GRID_SIZE/conv / 2.0;
               init_grid[valid_count][3] = fcount;
               init_grid[valid_count][4] = 0.0; //The tag for merging (shoul be structure, for storing float and char)
               valid_count +=1;
            }
         }
      }
   }

   *init_count = valid_count;
   //printf("valid_count grid_num: %d %d\n", valid_count, grid_num);
   for (i = 0; i < valid_count; i++) {
      occ = init_grid[i][3] / fframe * 100.0;
      printf("Initial hydration site %2d: %f %f %f; Frames: %f, Occupancies: %f%\n", i, init_grid[i][0], init_grid[i][1], init_grid[i][2], init_grid[i][3], occ);
      fprintf(InitWatSites, "HETATM%5d O    H2O%6d      %5.3f  %5.3f  %5.3f %4.2f %4.2f           O-1\n", i, i, init_grid[i][0] * conv, init_grid[i][1] * conv, init_grid[i][2] * conv, 1.0, occ);
    }
   fclose(InitWatSites);

   for(i=0;i<grid_num[0];i++) {
      for(j=0;j<grid_num[1];j++) {
         free(count[i][j]);
      }
   }
   for(i=0;i<grid_num[0];i++) {
       free(count[i]);
   }
   free(count);

   return 0;
}

int merge(float (*init_grid)[5], int init_count, float (*merge_grid)[5], int * merge_count) {
   int i, j, k, temp_count;
   float xx, yy, zz, dist2;
   k = 0;
   for (i=0;i<init_count;i++) {
      if (DEBUG == 1) {
         printf("Checking if label changed: init_count ID %d and its label: %f\n", i, init_grid[i][4]);
      }
      if (init_grid[i][4] == 0.0) {
         temp_count = 1;
         xx = init_grid[i][0];
         yy = init_grid[i][1];
         zz = init_grid[i][2];
         init_grid[i][4] = 1.0;

         for (j=0;j<init_count;j++) {
            dist2 = calc_dist_squr(init_grid[j][0], init_grid[j][1], init_grid[j][2], xx, yy, zz); 
            if (dist2 < pow(0.1, 2) && dist2 > 0.0) { // max 6 adjacent sites, and excludes itself
               xx += init_grid[j][0];
               yy += init_grid[j][1];
               zz += init_grid[j][2];
               temp_count += 1;
               init_grid[j][4] = 1.0;
               if (DEBUG == 1) {
                  printf("Found %d adjacent site for the %d th Site; adding site ID %d to the %d th merged site\n", temp_count-1, i, j, k);
               }
            }
         }

         float tmp3 = (float) temp_count;
         xx /= tmp3;
         yy /= tmp3;
         zz /= tmp3;
         merge_grid[k][0] = xx;
         merge_grid[k][1] = yy;
         merge_grid[k][2] = zz;
	 merge_grid[k][3] = temp_count;
	 merge_grid[k][4] = 0.0;
         k += 1;
      }
   }
   *merge_count = k;
   if (DEBUG == 1) {
      printf("Initial count is %d; Merged count is %d\n", init_count, *merge_count);
   }
   return 0;
}

int iteration(awt * mywat, int tot_wat, float (*init_grid)[5], int init_count, int frames, int flag_nc, float rangeWAT) {
   int h, i, j, k, count, dcount;
   int temp_count, merge_count;
   float xnew, ynew, znew, xs_min, ys_min, zs_min, xs_max, ys_max, zs_max, occ;
   double tmpde, deviation;
   float merge_grid[500][5] = {0};
   float new_grid[500][5] = {0};
   float fframe = (float) frames;
   float (*ingr)[5];
   float conv = 10.0;
   FILE * printsite = fopen("WatSites_merged.pdb", "w");

   if (flag_nc == 0) {
      conv   = 1.0;
   }
   
   //merge_count = init_count;
   ingr = init_grid;
   
   for (h=0;h<ITER_TIMES;h++) {
      deviation = 0.000;
      tmpde     = 0.000;
      if (DEBUG == 1) {
         printf("\n\n\n\nThis is a seperator for printing merge info. at %d th round\n", h);
      }

      if (h == 0) {
         merge(ingr, init_count, merge_grid, &merge_count);
      }
      else {
         merge(ingr, temp_count, merge_grid, &merge_count);
      }
      if (DEBUG == 1) {
         printf("MERGE_COUNT after merge() is %d\n", merge_count);
      }

      dcount = 0;
// find water aroud each merged site
      for (i=0;i<merge_count;i++) {
         // The diameter here is 0.13, plus or minus CUTOFF/2
         xs_min = merge_grid[i][0] - CUTOFF / conv / 2.0;
         ys_min = merge_grid[i][1] - CUTOFF / conv / 2.0;
         zs_min = merge_grid[i][2] - CUTOFF / conv / 2.0;
         xs_max = merge_grid[i][0] + CUTOFF / conv / 2.0;
         ys_max = merge_grid[i][1] + CUTOFF / conv / 2.0;
         zs_max = merge_grid[i][2] + CUTOFF / conv / 2.0;
         xnew = 0.0;
         ynew = 0.0;
         znew = 0.0;
         count = 0;
//#pragma omp parallel for reduction(+:xnew, ynew, znew, count) private(j,k)
         for (k=0;k<tot_wat;k++) {
            if (mywat[k].x > xs_min && mywat[k].x < xs_max \
                  && mywat[k].y > ys_min && mywat[k].y < ys_max \
                  && mywat[k].z > zs_min && mywat[k].z < zs_max) {
               xnew += mywat[k].x;
               ynew += mywat[k].y;
               znew += mywat[k].z;
               count += 1;
            }
         }

         float fcount = (float) count;
         if (fcount > fframe * rangeWAT) {
            xnew /= count;
            ynew /= count;
            znew /= count;
            tmpde = calc_dist(merge_grid[i][0], merge_grid[i][1], merge_grid[i][2], xnew, ynew, znew);
            
            if (tmpde > deviation) {
               deviation = tmpde;
            }
            
            new_grid[i][0] = xnew;
            new_grid[i][1] = ynew;
            new_grid[i][2] = znew;
	         new_grid[i][3] = fcount;
            new_grid[i][4] = 0.0;
            dcount += 1;
            //printf("TEMP hydration site %2d: %f %f %f; Frames: %f at %d round %f%\n", i, new_grid[i][0], new_grid[i][1], new_grid[i][2], h, new_grid[i][3]);
         }
      }
      
      if (dcount == 0) {
         printf("Warning: There are no desired water sites found in the %d th round's iteration\n", h);
      }
      // Here the deviation should not be the average value. It should be compared to the largest distance:
      else if (deviation < 0.01/conv) {
         printf(" %d sites were found at round %d, max_deviation is %lf\n", temp_count, h, deviation);
         printf("The max deviation reached the drifting creterion of 0.01 Angstrom after %d rounds of iteration\n", h);
         break;
      }
      
      ingr = new_grid;
      temp_count = i;
      if (DEBUG == 1) {
         printf(" %d sites were found at round %d, max_deviation is %lf\n", temp_count, h, deviation);
      }
   } // Ending of the iteration
   
   for (i = 0; i < dcount; i ++)
   {
      occ = new_grid[i][3] / fframe * 100.0;
      printf("Final hydration site %2d: %f %f %f; Frames: %f, Occupancies: %f%\n", i, new_grid[i][0], new_grid[i][1], new_grid[i][2], new_grid[i][3], occ);
      fprintf(printsite, "HETATM%5d O    H2O%6d      %5.3f  %5.3f  %5.3f %4.2f %4.2f           O-1\n", i, i, new_grid[i][0] * conv, new_grid[i][1] * conv, new_grid[i][2] * conv, 1.0, occ);
   }
   fclose(printsite);

   return 0;
}
