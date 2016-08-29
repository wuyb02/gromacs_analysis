/**
 ** do_xxx.h
 ** By Yanbin, 01/22/2008
 **/

//#define __MG_X_DEBUG__ 1

#include "Mgromacs.h"
#include "Mmatrix.h"

#define PBCSEP 0.6
#define PI 3.1415926
#define maxList 200
#define MASSO 15.99940
#define MASSH 1.00800

float periodicity(float seem,float box);

int MolInGRP(t_topology top, int size, atom_id *index, char *grpname);
void PrintAtomType(t_topology top, int size, atom_id *index, char *grpname);
double LJSigmaByAtomName(char *atomname);

int RemovePBCSep(t_trxframe fr, atom_id *index, int per, int start);
int RemovePBCSep2(t_trxframe fr, t_topology top, 
                 int size, atom_id *index, char *grpname);

real GetDisOSOW(MMatrix_t posOW, int ii,
                MMatrix_t posOS, int nOS, int start, float L[3]);
int GetDisOSOW2(MMatrix_t posOW, int ii,
                MMatrix_t posOS, int nOS, 
		real *dis, int *iOS, int start, float L[3]);

int ReadRBin(MMatrix_t *rbin, char *file);
int fillDonorInformation(float donor[3], float hydrogen1[3], float hydrogen2[3],
                         MMatrix_t posO, MMatrix_t posH1, MMatrix_t posH2, 
                         int i);
int WhichBin(MMatrix_t rbin, int nbin, real dis);
int fillAcceptorInformation(float acceptor[3], MMatrix_t posO, int i);
float isHbond(float d[3],float h[3],float a[3], float rcut,float acut,float L[3]);
double i_prod(float a[3],float b[3]);
double dist1(float a[3], float b[3], float L[3]);
float cos_angle_new(float a[3],float b[3]);
double normNEW(float a[3]);
int NormVec(float a[3]);
float dist2(MMatrix_t posOS,int i,MMatrix_t posOW,int j,float L[3]);
int neighbourList(int nWt,MMatrixInt_t List,float rcut2,MMatrix_t posO,float L[3]);
int StorePosition(t_trxframe fr, t_topology top, MMatrix_t posOW, atom_id *OWIndex, int size, 
		  int nWt, int start, bool bMW);
void CombinePosition(MMatrix_t pos, MMatrix_t half1, MMatrix_t half2);
void ShiftPosition(MMatrix_t posOW, int half);
int neighbourList2(MMatrixInt_t List,MMatrix_t posOS,MMatrix_t posO,int nOS,int nWt,float rcut2,float L[3]);
void GetOrientation(MMatrix_t A, MMatrix_t aa);
