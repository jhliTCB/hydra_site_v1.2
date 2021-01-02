#include <string.h>
#include <omp.h>

#define LINELEN_MAX 120
#define MAXCHAR 200
#define MAXATOM 200000
#define MAXWAT 1000000000
#define MAXFRAME 100000
// SOME fixed parameters are defined here, no need to pass them
#define GRID_SIZE 1.20000 /* 1.2 nm? */
#define ITER_TIMES 100
#define CUTOFF 1.300000
#define DEBUG 1 // debug or verbose mode
//#define BOX_SIZE 3.200000

#define ARRAY_INIT(array, size)           \
{                                        \
   for (int i=0;i<size;i++)              \
   {                                     \
      array[i] = 0.0;                      \
   }                                     \
}

struct VECTOR {
   double v1;
   double v2;
   double v3;
};
typedef struct VECTOR vt;

typedef struct {
   int id;        // PDB atom id
   char name[10]; // PDB atom name
   char aa[20];   // PDB residue name
   int resno;     // PDB residue number
   int flag;      // protein: flag==0; water(oxygen): flag==1; orther: flag==2
   double x;
   double y;
   double z;
} ATOM;


typedef struct {
   int fid; //frameID
   int nid; //Wat counting ID at each frame!
   int aid; //original atom ID in topology
   double x;
   double y;
   double z;
} awt; // big structure of all water molecules among all frames
