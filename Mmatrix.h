#ifndef __MMATRIX_H__
#define __MMATRIX_H__

#include "Mgromacs.h"

typedef struct
{
  double **m;
  int rows, cols;
}MMatrix, *MMatrix_t;

typedef struct
{
  double ***m;
  int nx, ny, nz;
}MMatrix2, *MMatrix2_t;

typedef struct
{
  int **m;
  int rows, cols;
}MMatrixInt, *MMatrixInt_t;

MMatrix_t CreateMatrix(int rows,int cols);
void FreeMatrix(MMatrix_t A);
MMatrixInt_t CreateMatrixInt(int rows,int cols);
void FreeMatrixInt(MMatrixInt_t A);
double *CreateVector(int rows);
void FreeVector(double *v);

MMatrix2_t CreateMatrix2(int nx,int ny, int nz);
void FreeMatrix2(MMatrix2_t A);

void ResetMatrix(MMatrix_t A);
void PrintMatrix(MMatrix_t A, char *);
void MatrixNorm2(MMatrix_t A);
void MatrixNorm(MMatrix_t A);
void MatrixCopy(MMatrix_t B, MMatrix_t A);

MMatrix_t MatrixProduct(MMatrix_t A, MMatrix_t B);
MMatrix_t MatrixATA(MMatrix_t A);
double MatrixDet(MMatrix_t A);
MMatrix_t MatrixTranspose(MMatrix_t A);

void EqnSolver(MMatrix_t A, MMatrix_t b, MMatrix_t aa);//AX=b, get X

#endif
