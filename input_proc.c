#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<memory.h>
#include<limits.h>
#include "global.h"
#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "AmberNetcdf.h"
#include "netcdf.h"
#include "utils.h"

int proc_pdb(char * filename, int * atomnum, int * watnum, ATOM * atom)
{
   int id            = -1;
   int terindex      = 0;
   int numatom       = 0;
   int numwat        = 0;
   int resno         = -1;
   char atomname[10] = {"INITIAL"};
   char resname[20]  = {"INITIAL"};
   char line[LINELEN_MAX];
   double x, y, z;
   int tmpint, next;
   FILE *fpin;

   if ((fpin = fopen(filename, "r")) == NULL) {
      fprintf(stdout, "Error: Cannot open the input file: %s in rpdb(), exit\n", filename);
      exit(1);
   }
   //initial(50, atom);
   for (;;) {
      if (fgets(line, LINELEN_MAX, fpin) == NULL) {
         break;
      }

      if (strncmp("TER", line, 3) == 0) {
         terindex = 1;
         continue;
      }

      if (strncmp("ATOM", line, 4) == 0 || strncmp("HETATM", line, 6) == 0) {
         if (line[12] != ' ') {
            atomname[0] = line[12];
            atomname[1] = line[13];
            atomname[2] = line[14];
            atomname[3] = line[15];
            next        = 4;
         }
         else if (line[13] != ' ') {
            atomname[0] = line[13];
            atomname[1] = line[14];
            atomname[2] = line[15];
            next        = 3;
         }
         else {
            atomname[0] = line[14];
            atomname[1] = line[15];
            next        = 2;
         }
         
         if (line[16] == ' ') {
            atomname[next] = '\0';          
         }
         else {
            atomname[next] = line[16];
            atomname[next+1] = '\0';
         }
         
         if (strchr("0123456789", line[5]) != NULL) {
            sscanf(line + 5, "%d", &id); /* 'ATOM\s\d+' */
         }
         else {
            sscanf(line + 6, "%d", &id); /* 'HETATM' or 'ATOM\s\s' */
         }
         resname[0] = line[17];
         resname[1] = line[18];
         resname[2] = line[19];

         if (line[20] != ' ') {
            resname[3] = line[20];
            resname[4] = '\0';
         }
         else {
            resname[3] = '\0'; 
         }

         sscanf(&line[22], "%d%lf%lf%lf", &tmpint, &x, &y, &z);
         resno = tmpint;

         atom[numatom].id = id;
         strcpy(atom[numatom].name, atomname);
         strcpy(atom[numatom].aa, resname);
         atom[numatom].resno = resno;
         if (pro_res_comp(resname) == 1) {
            atom[numatom].flag = 0;
         }
         else if (wat_res_comp(resname) == 1 && judge_hydrogen(atomname) == 0) {
            numwat ++;
            atom[numatom].flag = 1;
         }
         else {
            atom[numatom].flag =2;
         }
         atom[numatom].x = x;
         atom[numatom].y = y;
         atom[numatom].z = z;

         numatom++;
      }
   }
   //sscanf(&numatom, "%d", atomnum);
   *atomnum = numatom;
   *watnum  = numwat;
   fclose(fpin);
   //return numatom;
   return 0;
}  


int proc_xtc(char * filename, int atomnum, int watnum, ATOM * atom, awt * mywat, vt * num_grid, vt * min_grid_edge, int * frames, int * tot_wat)
{
   int natoms, natom, step, read_return, numwat, i, j;
   float time, p;
   float * proX, * proY, * proZ, * f_ctX, * f_ctY, * f_ctZ, * f_edX, * f_edY, * f_edZ;
   float a = 0, b = 0, c = 0, aa = 0, bb = 0, cc = 0;
   vt * boxct;
   //float * distance;
   FILE *fpin;
   int    frame = 0;
   int len_wat  = 1;
   double diag1 = 0.0;
   double diag2 = 0.0;
   double diag3 = 0.0;

   if ((fpin = fopen(filename, "r")) == NULL) {
      fprintf(stdout, "Cannot open the input file: %s in rxtc(), exit\n",
                     filename);
      exit(1);
   }
   
   proX  = (float*) malloc(sizeof(float) * atomnum);
   proY  = (float*) malloc(sizeof(float) * atomnum);
   proZ  = (float*) malloc(sizeof(float) * atomnum);
   f_ctX = (float*) malloc(sizeof(float) * atomnum);
   f_ctY = (float*) malloc(sizeof(float) * atomnum);
   f_ctZ = (float*) malloc(sizeof(float) * atomnum);
   f_edX = (float*) malloc(sizeof(float) * atomnum);
   f_edY = (float*) malloc(sizeof(float) * atomnum);
   f_edZ = (float*) malloc(sizeof(float) * atomnum);
  
   matrix box; // an empty matrix just for read_xtc()
   rvec *x;    // coordinates
   XDRFILE *xtc;
   xtc = xdrfile_open(filename, "r");
   read_return = read_xtc_natoms(filename, &natoms);
   x = (rvec *)calloc(natoms, sizeof(x[0]));

   if (natoms != atomnum) {
      printf("ERROR: Mismatch between trajectory (%d atoms) and topology (%d atoms)\n", natoms,atomnum);
      exit(1);
   }

   while (1) {
      read_return = read_xtc(xtc,natoms,&step,&time,box,x,&p);
      if (read_return != 0) {
         boxct->v1         = 0.5 * diag1 / frame; // This is the basis vector, it needs to multiply the real box size!
         boxct->v2         = 0.5 * diag2 / frame;
         boxct->v3         = 0.5 * diag3 / frame;
         num_grid->v1      = ((sum(f_edX, frame) - sum(f_ctX, frame)) / frame + 0.8) / (GRID_SIZE/10) + 1;
         num_grid->v2      = ((sum(f_edY, frame) - sum(f_ctY, frame)) / frame + 0.8) / (GRID_SIZE/10) + 1;
         num_grid->v3      = ((sum(f_edZ, frame) - sum(f_ctZ, frame)) / frame + 0.8) / (GRID_SIZE/10) + 1;
         printf("check grid_num: %f %f %f\n", num_grid->v1, num_grid->v2, num_grid->v3);
         min_grid_edge->v1 = sum(f_ctX, frame) / frame - 0.4;
         min_grid_edge->v2 = sum(f_ctY, frame) / frame - 0.4;
         min_grid_edge->v3 = sum(f_ctZ, frame) / frame - 0.4;
         
         free(proX);   free(proY);     free(proZ);
         free(f_ctX);  free(f_ctY);    free(f_ctZ);
         free(f_edX);  free(f_edY);    free(f_edZ);
         proX  = NULL;
         proY  = NULL;
         proZ  = NULL;
         f_ctX = NULL;
         f_ctY = NULL;
         f_ctZ = NULL;
         f_edX = NULL;
         f_edY = NULL;
         f_edZ = NULL;
         break;
      }
      
      if (frame == 0) { //abort the first frame, which is the topology with the x of (0, 0, 0)
         frame = 1;
         continue;
      }
      
      diag1 += box[0][0]; // The diagnose for rectangle simulation box
      diag2 += box[1][1];
      diag3 += box[2][2];

      numwat = 1;
      i      = 0;
      // Process all pdb atoms frame by frame
      for (natom=0;natom<=natoms;natom++) {
         // assume that protein is always in front of water, so "continue" is used here
         if (atom[natom].flag == 0) {
            proX[i] = x[natom][0];
            proY[i] = x[natom][1];
            proZ[i] = x[natom][2];
            i += 1;
            continue;
         }
         else if (natom == i + 1) {
            a  = max_fbl(proX, i);
            b  = max_fbl(proY, i);
            c  = max_fbl(proZ, i);
            aa = min_fbl(proX, i);
            bb = min_fbl(proY, i);
            cc = min_fbl(proZ, i);
            f_ctX[frame] = aa;
            f_ctY[frame] = bb;
            f_ctZ[frame] = cc;
            f_edX[frame] = a;
            f_edY[frame] = b;
            f_edZ[frame] = c;
         }

         //if (wat_res_comp(atom[natom].aa) == 1 && judge_hydrogen(atom[natom].name) == 0) {
         if (atom[natom].flag == 1) {
            if (x[natom][0] > aa-0.4 && x[natom][0] < a+0.4 && \
                x[natom][1] > bb-0.4 && x[natom][1] < b+0.4 && \
                x[natom][2] > cc-0.4 && x[natom][2] < c+0.4  ) {
               mywat[len_wat].fid = frame; 
               mywat[len_wat].nid = numwat;
               mywat[len_wat].aid = atom[natom].id;
               mywat[len_wat].x   = x[natom][0];
               mywat[len_wat].y   = x[natom][1];
               mywat[len_wat].z   = x[natom][2];
               numwat ++;
               len_wat ++;
            }
         }
      }
      
      frame ++;
   }

   *frames  = frame;
   *tot_wat = len_wat;
   i        = len_wat - 1;
   
   printf("*tot_wat is %d, tot_wat is %p\n", len_wat, tot_wat);
   printf("check mywat:\ni=%d\nmywat[i].fid=%d\nmywat[i].nid=%d\nmywat[i].aid=%d\nmywat[i].x=%f\nmywat[i].y=%f\nmywat[i].z=%f\n", i, mywat[i].fid, mywat[i].nid, mywat[i].aid, mywat[i].x, mywat[i].y, mywat[i].z);
   printf("check mywat again: mywat[1].fid %d, mywat[1].nid %d, mywat[1].aid %d, mywat[1].x %f, mywat[1].y %f, mywat[1].z %f, box %f\n", mywat[1].fid, mywat[1].nid, mywat[1].aid, mywat[1].x, mywat[1].y, mywat[1].z, boxct->v1);

   xdrfile_close(xtc);
   return 0;
}

int proc_nc(ATOM * atom, awt * mywat, struct AmberNetcdf * A, double * X, vt * num_grid, vt * min_grid_edge, int atomnum, int * tot_wat) {
   netcdfInfo(A); // print the basic information of the loaded nc file
   int i = 0;
   int h, j, numwat, natom; //h number of protein atoms
   int frames = A->ncframe;
   //float * distance;
   float * proX, * proY, * proZ, * f_ctX, * f_ctY, * f_ctZ, * f_edX, * f_edY, * f_edZ;
   float a = 0, b = 0, c = 0, aa = 0, bb = 0, cc = 0;
   double box[9]; // just for reading the nc frames
   int r_return = 1;
   int cycle = 0;
   A->currentFrame = 1;
   
   proX  = (float*) malloc(sizeof(float) * atomnum);
   proY  = (float*) malloc(sizeof(float) * atomnum);
   proZ  = (float*) malloc(sizeof(float) * atomnum);
   f_ctX = (float*) malloc(sizeof(float) * frames);
   f_ctY = (float*) malloc(sizeof(float) * frames);
   f_ctZ = (float*) malloc(sizeof(float) * frames);
   f_edX = (float*) malloc(sizeof(float) * frames);
   f_edY = (float*) malloc(sizeof(float) * frames);
   f_edZ = (float*) malloc(sizeof(float) * frames);
   //ARRAY_INIT(proX, atomnum); // no need to initial here

   //for (int i=0;i<A->ncframe;i++) {
   while (1) {
      cycle += 1;
      h      = 0;
      numwat = 0;
      natom  = 0;
      ARRAY_INIT(box, 9);
      netcdfGetFrame(A, A->currentFrame, X, box);
      for (natom=0; natom<=atomnum; natom++) {
         /* assume that protein is always in front of water, so "continue" is used here */
         if (atom[natom].flag == 0) {
            proX[natom] = X[natom*3];
            proY[natom] = X[natom*3+1];
            proZ[natom] = X[natom*3+2];
            h += 1;
            continue;
         }
         else if (natom == h + 1) {  // make sure that all protein atoms were placed at the begining!
            a  = max_fbl(proX, h);
            b  = max_fbl(proY, h);
            c  = max_fbl(proZ, h);
            aa = min_fbl(proX, h);
            bb = min_fbl(proY, h);
            cc = min_fbl(proZ, h);
            f_ctX[A->currentFrame-1] = aa;
            f_ctY[A->currentFrame-1] = bb;
            f_ctZ[A->currentFrame-1] = cc;
            f_edX[A->currentFrame-1] = a;
            f_edY[A->currentFrame-1] = b;
            f_edZ[A->currentFrame-1] = c;
         }
         
         if (atom[natom].flag == 1) {
            if (X[natom*3]   > aa-4 && X[natom*3]   < a+4 &&  \
                X[natom*3+1] > bb-4 && X[natom*3+1] < b+4 &&  \
                X[natom*3+2] > cc-4 && X[natom*3+2] < c+4   ) {
               memcpy(&(mywat[i].fid), &(A->currentFrame), sizeof(int));
               memcpy(&(mywat[i].nid), &numwat, sizeof(int));
               memcpy(&(mywat[i].aid), &(atom[natom].id), sizeof(int));
               memcpy(&(mywat[i].x), &(X[natom*3]), sizeof(double));
               memcpy(&(mywat[i].y), &(X[natom*3+1]), sizeof(double));
               memcpy(&(mywat[i].z), &(X[natom*3+2]), sizeof(double));
               i += 1;
               numwat += 1;
            }
         }
      }
      
      r_return = netcdfGetNextFrame(A, X, box);
      if (r_return == 0) {
         *tot_wat = i - 1;
         i = i - 1;
         frames = frames - 1;
         // average the protein min and max, assuming there is no large flutuation
         num_grid->v1      = ((sum(f_edX, frames) - sum(f_ctX, frames)) / frames + 8)/ GRID_SIZE + 1;
         num_grid->v2      = ((sum(f_edY, frames) - sum(f_ctY, frames)) / frames + 8)/ GRID_SIZE + 1;
         num_grid->v3      = ((sum(f_edZ, frames) - sum(f_ctZ, frames)) / frames + 8)/ GRID_SIZE + 1;
         printf("check grid_num: %f %f %f\n", num_grid->v1, num_grid->v2, num_grid->v3);
         min_grid_edge->v1 = sum(f_ctX, frames) / frames - 4;
         min_grid_edge->v2 = sum(f_ctY, frames) / frames - 4;
         min_grid_edge->v3 = sum(f_ctZ, frames) / frames - 4;
            
         printf("i-1 is %d, *tot_wat is %d, tot_wat is %p, cycle is %d\n", i, *tot_wat, tot_wat, cycle);
         printf("check mywat: mywat[i].fid %d, mywat[i].nid %d, mywat[i].aid %d, mywat[i].x %f, mywat[i].y %f, mywat[i].z %f\n", mywat[i].fid, mywat[i].nid, mywat[i].aid, mywat[i].x, mywat[i].y, mywat[i].z);
         printf("check mywat again: mywat[1].fid %d, mywat[1].nid %d, mywat[1].aid %d, mywat[1].x %f, mywat[1].y %f, mywat[1].z %f, box %f\n", mywat[1].fid, mywat[1].nid, mywat[1].aid, mywat[1].x, mywat[1].y, mywat[1].z, box);
         
         free(proX);  free(proY);  free(proZ);
         free(f_ctX); free(f_ctY); free(f_ctZ);
         free(f_edX); free(f_edY); free(f_edZ);
         proX  = NULL;
         proY  = NULL;
         proZ  = NULL;
         f_ctX = NULL;
         f_ctY = NULL;
         f_ctZ = NULL;
         f_edX = NULL;
         f_edY = NULL;
         f_edZ = NULL;
         break;
      }
   }
   return 0;
}
