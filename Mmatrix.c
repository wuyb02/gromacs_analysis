#include "Mmatrix.h"

MMatrix_t CreateMatrix(int rows,int cols)
{
  int i, j;
  MMatrix_t A = (MMatrix_t)calloc(1,sizeof(MMatrix));

  A->rows=rows;
  A->cols=cols;

  if(A->rows==0)
  {
    A->m = NULL;
    return A;
  }
  A->m = calloc((unsigned int)A->rows,sizeof(double *));
  for(i=0; i<A->rows; i++) 
  {
    A->m[i] = calloc((unsigned int)A->cols,sizeof(double));
    for(j=0; j<A->cols; j++)
    {
      A->m[i][j] = 0.0;
    }
  }

  return A;
}

void FreeMatrix(MMatrix_t A)
{
  int i;

  for(i=0; i<A->rows; i++)
  {
    free(A->m[i]);
  }
  free(A->m);
  free(A);
}

MMatrixInt_t CreateMatrixInt(int rows,int cols)
{
  int i, j, k;
  MMatrixInt_t A=(MMatrixInt_t)calloc(1,sizeof(MMatrixInt));

  A->rows = rows;
  A->m = (int **)calloc((unsigned int)A->rows,sizeof(int *));
  A->cols = cols;
  for(i=0; i<A->rows; i++)
  {
     A->m[i] = (int *)calloc((unsigned int)A->cols,sizeof(int ));
     for(j=0; j<A->cols; j++)
     {
       A->m[i][j] = 0;
     }
  }

  return A;
}

void FreeMatrixInt(MMatrixInt_t A)
{
  int i;

  for(i=0; i<A->rows; i++)
  {
    free(A->m[i]);
  }
  free(A->m);
  free(A);
}

double *CreateVector(int rows) 
{ 
  return calloc((unsigned int)rows, sizeof(double)); 
}

void FreeVector(double *v) 
{ 
  free(v);
}

void ResetMatrix(MMatrix_t A)
{
  int i, j, k;

  for(i=0; i<A->rows; i++)
  {
    for(j=0; j<A->cols; j++)
    {
      A->m[i][j] = 0.0;
    }
  }
}

void PrintMatrix(MMatrix_t A, char *name)
{
  int i,j,k;

  printf("###Matrix %s ###\n", name);
  for(i=0;i<A->rows;i++)
  {
    for(j=0;j<A->cols;j++)
    {
      printf("%11.3f", A->m[i][j]);
    }
    printf("\n");
  }
}

MMatrix2_t CreateMatrix2(int nx,int ny, int nz)
{
  int i, j, k;
  MMatrix2_t A = (MMatrix2_t)calloc(1,sizeof(MMatrix2));
  A->nx = nx;
  A->ny = ny;
  A->nz = nz;

  if(A->nx==0)
  {
    A->m = NULL;
    return A;
  }
  A->m = calloc((unsigned int)A->nx,sizeof(double **));

  for(i=0; i<A->nx; i++) 
  {
    A->m[i] = calloc((unsigned int)A->ny, sizeof(double *));
    for(j=0; j<A->ny; j++)
    {
      A->m[i][j] = calloc((unsigned int)A->nz, sizeof(double));
      for(k=0; k<A->nz; k++)
      {
	A->m[i][j][k] = 0.0;
      }
    }
  }

  return A;
}

void FreeMatrix2(MMatrix2_t A)
{
  int i, j, k;

  for(i=0; i<A->nx; i++)
  {
    for(j=0; j<A->ny; j++)
    {
      free(A->m[i][j]);
    }
    free(A->m[i]);
  }
  free(A->m);
  free(A);
}


void MatrixNorm(MMatrix_t A)
{
  int i,j,k;
  double norm;

  for(i=0;i<A->rows;i++)
  {
    for(norm=0.0,j=0;j<A->cols;j++)
    {
      norm += A->m[i][j]*A->m[i][j];
    }
    norm=sqrt(norm);

    printf("\nnorm=%11.3f\n", norm);

    if(fabs(norm)<1e-10)
    {
      continue;
    }
    for(j=0;j<A->cols;j++)
    {
      A->m[i][j]=A->m[i][j]/norm;
    }
  }
}

void MatrixNorm2(MMatrix_t A)
{
  int i,j,k;
  double norm;

  for(i=0;i<A->cols;i++)
  {
    for(norm=0.0,j=0;j<A->rows;j++)
    {
      norm += A->m[j][i]*A->m[j][i];
    }
    norm=sqrt(norm);

//    printf("\nnorm=%11.3f\n", norm);

    if(fabs(norm)<1e-10)
    {
      continue;
    }
    for(j=0;j<A->rows;j++)
    {
      A->m[j][i]=A->m[j][i]/norm;
    }
  }
}

void MatrixCopy(MMatrix_t B, MMatrix_t A)
{
  int i,j,k;

  if(B->rows!=A->rows)
  {
    printf("MMatrixCopy, A->rows!=B->rows, %d, %d\n", A->rows, B->rows);
    exit(0);
  }
  if(B->cols!=A->cols)
  {
    printf("MMatrixCopy, A->cols!=B->cols, %d, %d\n", A->cols, B->cols);
    exit(0);
  }

  for(i=0;i<A->rows;i++)
  {
    for(j=0;j<A->cols;j++)
    {
      B->m[i][j]=A->m[i][j];
    }
  }
}

MMatrix_t MatrixTranspose(MMatrix_t A)
{
  int i,j,k;
  MMatrix_t B=CreateMatrix(A->cols,A->rows);

  for(i=0;i<A->rows;i++)
  {
    for(j=0;j<A->cols;j++)
    {
      B->m[j][i]=A->m[i][j];
    }
  }

  return B;
}

MMatrix_t MatrixATA(MMatrix_t A)
{
  int i, j, k;

  MMatrix_t AA = CreateMatrix(A->cols, A->cols);

  for(i=0; i<A->cols; i++)
  {
    for(j=0; j<A->cols; j++)
    {
      AA->m[i][j]=0.0;
      for(k=0; k<A->rows; k++)
      {
	AA->m[i][j]+=A->m[i][k]*A->m[j][k];
      }
    }
  }

  return AA;
}

MMatrix_t MatrixProduct(MMatrix_t A, MMatrix_t B)
{
  int i, j, k;

  if(A->cols!=B->rows)
  {
    printf("MatrixProduct, A->cols!=B->rows\n");
    exit(0);
  }

  MMatrix_t AB = CreateMatrix(A->rows, B->cols);
  for(i=0; i<A->rows; i++)
  {
    for(j=0; j<B->cols; j++)
    {
      AB->m[i][j]=0.0;
      for(k=0; k<A->cols; k++)
      {
	AB->m[i][j]+=A->m[i][k]*B->m[k][j];
      }
    }
  }

  return AB;
}

double MatrixDet(MMatrix_t A)
{
  int i, j, k;
  int n=A->rows;
  double result=0.0;

  if(A->rows!=A->cols)
  {
    printf("MatrixDet, A->rows!=A->cols, %d, %d\n", A->rows, A->cols);
    exit(0);
  }

  if(n==1)
  {
    result=A->m[0][0];
    return result;
  }

  for(i=0; i<n; i++)
  {
    MMatrix_t B=CreateMatrix(n-1,n-1);
    for(j=1; j<n; j++)
    {
      for(k=0; k<i; k++)
      {
	B->m[j-1][k]=A->m[j][k];
      }
      for(k=i+1; k<n; k++)
      {
	B->m[j-1][k-1]=A->m[j][k];
      }
    }
    if(i%2)
    {
      result+=A->m[0][i]*MatrixDet(B);
    }else
    {
      result-=A->m[0][i]*MatrixDet(B);
    }
    FreeMatrix(B);
  }
  return result;
}

void EqnSolver(MMatrix_t A, MMatrix_t b, MMatrix_t aa)
{
  int i,j,k;
  int rows=A->rows;

  if(A->rows!=b->rows)
  {
    printf("A->rows!=b->rows\n");
    exit(0);
  }

//  PrintMatrix(A, "A");
  double D=MatrixDet(A);
  if(D<=1e-5&&D>=-1e-5)
  {
    printf("D=%11.f, too small\n", D);
  }
  MMatrix_t TempA=CreateMatrix(A->rows, A->cols);
  for(i=0;i<rows;i++)
  {
    MatrixCopy(TempA, A);
    for(j=0;j<rows;j++)
    {
      TempA->m[j][i]=b->m[j][0];
    }
//    PrintMatrix(TempA, "TempA");
    aa->m[i][0]=MatrixDet(TempA)/D;
  }
  FreeMatrix(TempA);
}
