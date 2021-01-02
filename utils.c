#include <math.h>
#include <string.h>

int wat_res_comp(char sample[]) {
   static char *wat_res_name[] = {"WAT", "SOL", "H2O", "HOH", "TIP3"};
   int indicator = 0;
   int i;
   for (i = 0; i < 4; i++) {if ( strncmp(wat_res_name[i],sample,3) == 0 ) {indicator = 1; break;}}
   if ( strncmp(wat_res_name[4],sample,4) == 0 ) {indicator = 1;}
   return (indicator); 
}
   
int pro_res_comp(char sample[]) { // LIG flag should be added here for those systems with ligand(s) on the surface
   static char *pro_res_name[] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN",
                                  "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
                                  "MET", "PHE", "PRO", "SER", "THR", "TRP",
                                  "TYR", "VAL", "HID", "HIE", "HIP", "GLH", "CYI", "HEM" 
                                  "ASH", "HSD", "HSE", "HSP", "GLUP", "ASPP"};
   int indicator = 0;
   int i;
   for (i=0;i<27;i++) {if ( strncmp(pro_res_name[i], sample, 3) == 0 ) {indicator = 1; break;}}
   if ( strncmp(pro_res_name[28], sample, 4) == 0 ) {indicator = 1;}
   if ( strncmp(pro_res_name[29], sample, 4) == 0 ) {indicator = 1;}
   return(indicator);
}

int judge_hydrogen(char sample[]) {
   static char *gmx_num_H[] = {"1H", "2H", "3H", "4H", "5H", "6H"};
   int indicator = 0;
   //int len_atm = strlen(sample);
   int i;
   for (i = 0; i < 6; i++) {if (strncmp(gmx_num_H[i],sample,2) == 0) {indicator = 1; break;}}
   if (strncmp("H", sample, 1) == 0) {indicator = 1;}
   return indicator;
}
   
float calc_dist(float x1, float y1, float z1, float x2, float y2, float z2) {
   return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}

float calc_dist_squr(float x1, float y1, float z1, float x2, float y2, float z2) {
   return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
}

float min_fbl(float * a, int n) {
   int i;
   float temp;
   temp = *a;
   for (i=0;i<n;i++)
      if (temp > *(a+i))
         temp = *(a+i);
   
   return temp;
}

float max_fbl(float * a, int n) {
   int i;
   float temp;
   temp = *a;
   for (i=0;i<n;i++)
      if (temp < *(a+i))
         temp = *(a+i);
   
   return temp;
}

float sum(float * a, int n) {
   int i;
   float sum = 0;
   for (i=0;i<n;i++)
      sum += *(a+i);
   
   return(sum);
}
