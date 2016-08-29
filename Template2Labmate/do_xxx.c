# include "do_xxx.h"

float periodicity(float seem,float box)
{
  float actual;
  if(seem>box/2)
    actual = seem-box;
  else if(seem<-box/2)
    actual = box+seem;
  else
    actual = seem;
  return actual;
}

int MolInGRP(t_topology top,
             int size, atom_id *index, char *grpname)
{
  int i, j, k;
  int num;

  if(size<=0)
  {
    printf("MolInGRP, size<0\n");
    exit(0);
  }

  for(i=1, num=1; i<size; i++)
  {
    if(top.atoms.atom[index[i]].resnr!=top.atoms.atom[index[i-1]].resnr)
    {
      num++;
    }
    top.atoms.atom[index[i]].type;
    *(top.atoms.atomtype[index[i]]);
  }
//  printf("\nMolInGRP returned\n");

  return num;
}

void PrintAtomType(t_topology top, int size, atom_id *index, char *grpname)
{
  int i, j, k;

  if(size<=0)
  {
    printf("MolInGRP, size<0\n");
    exit(0);
  }

  for(i=0; i<size; i++)
  {
    printf("i=%10d, type=%10d, atomtypeName=%10s, sigma=%10.4f\n", 
	index[i]+1, top.atoms.atom[index[i]].type, 
	*(top.atoms.atomtype[index[i]]), LJSigmaByAtomName(*(top.atoms.atomtype[index[i]])));
  }
}

double LJSigmaByAtomName(char *atomname)
{
  if(strcmp(atomname, "OP")==0)
  {
    return 0.2850;
  }

  if(strcmp(atomname, "CP")==0)
  {
    return 0.3440;
  }

  if(strcmp(atomname, "HP")==0)
  {
    return 0.3000;
  }

  if(strcmp(atomname, "CAR")==0)
  {
    return 0.3550;
  }

  if(strcmp(atomname, "CA1R")==0)
  {
    return 0.3550;
  }

  if(strcmp(atomname, "CT1R")==0)
  {
    return 0.4054;
  }

  if(strcmp(atomname, "CT2R")==0)
  {
    return 0.3875;
  }

  if(strcmp(atomname, "CT3R")==0)
  {
    return 0.3671;
  }

  if(strcmp(atomname, "CDR")==0)
  {
    return 0.3564;
  }

  if(strcmp(atomname, "OSR")==0)
  {
    return 0.3154;
  }

  if(strcmp(atomname, "OBR")==0)
  {
    return 0.3029;
  }

  if(strcmp(atomname, "NC2R")==0)
  {
    return 0.3296;
  }

  if(strcmp(atomname, "HPR")==0)
  {
    return 0.2420;
  }

  if(strcmp(atomname, "HCR")==0)
  {
    return 0.0400;
  }

  if(strcmp(atomname, "HAR")==0)
  {
    return 0.2352;
  }

  if(strcmp(atomname, "CL")==0)
  {
    return 0.4045;
  }

  if(strcmp(atomname, "NA")==0)
  {
    return 0.2430;
  }

  if(strcmp(atomname, "OW")==0)
  {
    return 0.3166;
  }

  if(strcmp(atomname, "HW")==0)
  {
    return 0.0000;
  }

  printf("do_xxx.c, LJSigmaByAtomName, No match, %10s\n", atomname);
  exit(0);
}

int RemovePBCSep(t_trxframe fr, atom_id *index, int per, int start)
{
  int i, j, k, ll;
  int num;

//  int isize = size;
//  for(j=0; j<isize; j++)
//  {
//    printf("%d  ", index[j]);
//    if(!((j+1)%10))
//      printf("\n");
//  }
//  printf("\n\n");

  for(j=1,num=0; j<per; j++)
  {
    ll=j+start;
    for(k=0; k<DIM; k++)
    {
//      printf("ll=%d,j=%d,k=%d,%15.4f  ", ll, j, k, fr.x[index[j]][k]);
      if(fr.x[index[ll]][k]-fr.x[index[ll-1]][k] > PBCSEP*fr.box[k][k])
      {
	fr.x[index[ll]][k] -= fr.box[k][k];
	num++;
//	printf("ll=%d,j=%d,k=%d,%f-->%f\n", ll, j, k, fr.x[index[j]][k]+fr.box[k][k], fr.x[index[j]][k]);
      }
      if(fr.x[index[ll]][k]-fr.x[index[ll-1]][k] < -PBCSEP*fr.box[k][k])
      {
	fr.x[index[ll]][k] += fr.box[k][k];
	num++;
//	printf("ll=%d,j=%d,k=%d,%f-->%f\n", ll, j, k, fr.x[index[j]][k]-fr.box[k][k], fr.x[index[j]][k]);
      }
    }
//    printf("\n");
  }

  return num;
}

int RemovePBCSep2(t_trxframe fr, t_topology top, 
                 int size, atom_id *index, char *grpname)
{
  int i, j, k;
  int num;

//  int isize = size;
//  for(j=0; j<isize; j++)
//  {
//    printf("%d  ", index[j]);
//    if(!((j+1)%10))
//      printf("\n");
//  }
//  printf("\n\n");

  int resnr = MolInGRP(top, size, index, grpname);
  int per=size/resnr;
//  printf("resnr=%d, per=%d\n", resnr, per);

  if(size%resnr)
  {
    printf("RemovePBCSep, size%%resnr!=0\n");
    exit(0);
  }

  for(i=0, num=0; i<resnr; i++)
  {
    for(j=1+i*per; j<per; j++)
    {
      for(k=0; k<DIM; k++)
      {
	if(fr.x[index[j]][k]-fr.x[index[j-1]][k] > PBCSEP*fr.box[k][k])
	{
	  fr.x[index[j]][k] -= fr.box[k][k];
	  num++;
//	  printf("%f-->%f\n", fr.x[index[j]][k]+fr.box[k][k], fr.x[index[j]][k]);
	}
	if(fr.x[index[j]][k]-fr.x[index[j-1]][k] < -PBCSEP*fr.box[k][k])
	{
	  fr.x[index[j]][k] += fr.box[k][k];
	  num++;
//	  printf("%f-->%f\n", fr.x[index[j]][k]-fr.box[k][k], fr.x[index[j]][k]);
	}
      }
    }
  }

  return num;
}

real GetDisOSOW(MMatrix_t posOW, int ii,
                MMatrix_t posOS, int nOS, int start, float L[3])
{
  int i, j, k, ll, mm, nn;
  real dis, tempDis;

  if(ii>=posOW->rows)
  {
    printf("do_msd.c, GetDisOSOW, ii>=pos->rows, %d>=%d", ii, posOW->rows);
    exit(0);
  }

  for(i=0, dis=100.0; i<nOS; i++)
  {
    ll = start+i;
    for(tempDis=0.0, k=0; k<DIM; k++)
    {
#ifdef __MG_DEBUG__
      tempDis += pow(periodicity(2.0-posOW->m[ii][k], L[k]), 2);
      printf("define __MG_DEBUG\n");
#elif __MG_X_DEBUG__
      tempDis += pow(periodicity(2.0-posOW->m[ii][k], L[k]), 2);
      printf("define __MG_X_DEBUG\n");
      break;
#else
      tempDis += pow(periodicity(posOS->m[ll][k]-posOW->m[ii][k], L[k]), 2);
#endif
    }
    tempDis = sqrt(tempDis);
    if(tempDis<dis)
    {
      dis = tempDis;
    }
//    printf("%8.3f,%8.3f,%8.3f, %8.3f,%8.3f,%8.3f, tempDis=%11.3f, dis=%11.3f, ii=%d, ll=%d\n", 
//	   posOW->m[ii][0], posOW->m[ii][1], posOW->m[ii][2], 
//	   posOS->m[ll][0], posOS->m[ll][1], posOS->m[ll][2],
//	   tempDis, dis, ii, ll);
  }

  return dis;
}

int GetDisOSOW2(MMatrix_t posOW, int ii,
                MMatrix_t posOS, int nOS, 
		real *dis, int *iOS, int start, float L[3])
{
  int i, j, k, ll, mm, nn;
  real tempDis;

  if(ii>=posOW->rows)
  {
    printf("do_msd.c, GetDisOSOW, ii>=pos->rows, %d>=%d", ii, posOW->rows);
    exit(0);
  }

  for(i=0, *dis=100.0; i<nOS; i++)
  {
    ll = start+i;
    for(tempDis=0.0, k=0; k<DIM; k++)
    {
#ifdef __MG_DEBUG__
      tempDis += pow(periodicity(2.0-posOW->m[ii][k], L[k]), 2);
      printf("define __MG_DEBUG\n");
#elif __MG_X_DEBUG__
      tempDis += pow(periodicity(2.0-posOW->m[ii][k], L[k]), 2);
      printf("define __MG_X_DEBUG\n");
      break;
#else
      tempDis += pow(periodicity(posOS->m[ll][k]-posOW->m[ii][k], L[k]), 2);
#endif
    }
    tempDis = sqrt(tempDis);
    if(tempDis<*dis)
    {
      *dis = tempDis;
      *iOS = i;
    }
//    printf("%8.3f,%8.3f,%8.3f, %8.3f,%8.3f,%8.3f, tempDis=%11.3f, dis=%11.3f, ii=%d, ll=%d\n", 
//	   posOW->m[ii][0], posOW->m[ii][1], posOW->m[ii][2], 
//	   posOS->m[ll][0], posOS->m[ll][1], posOS->m[ll][2],
//	   tempDis, *dis, ii, ll);
  }

  return 1;
}

int ReadRBin(MMatrix_t *rbin, char *file)
{
  int i, j, k;
  int nrbin;
  float temp;

  FILE *RBIN = fopen(file, "r");
  if(RBIN==NULL)
  {
    printf("ReadRBIN, %s doesnot exist\n", file);
    exit(0);
  }

  fscanf(RBIN, "%d", &nrbin);
  printf("nrbin=%d\n", nrbin);
  *rbin = CreateMatrix(nrbin+1, 1);
  MMatrix_t rbin2 = *rbin;

  rbin2->m[0][0] = 0.0;
  for(i=0; i<nrbin; i++)
  {
    fscanf(RBIN, "%f", &temp);
    rbin2->m[i+1][0] = temp;
  }
//  for(i=0; i<nrbin+1; i++)
//  {
//    printf("%10.5f", rbin2->m[i][0]);
//  }
//  printf("\n");

  return nrbin;
}

int fillDonorInformation(float donor[3], float hydrogen1[3], float hydrogen2[3],
                         MMatrix_t posO, MMatrix_t posH1, MMatrix_t posH2, 
                         int i)
{
  donor[0] = posO->m[i][0];
  donor[1] = posO->m[i][1];
  donor[2] = posO->m[i][2];
 
  hydrogen1[0] = posH1->m[i][0];
  hydrogen1[1] = posH1->m[i][1];
  hydrogen1[2] = posH1->m[i][2];
 
  hydrogen2[0] = posH2->m[i][0];
  hydrogen2[1] = posH2->m[i][1];
  hydrogen2[2] = posH2->m[i][2];
 
  return 0;
}

int fillAcceptorInformation(float acceptor[3], MMatrix_t posO, int i)
{
  acceptor[0] = posO->m[i][0];
  acceptor[1] = posO->m[i][1];
  acceptor[2] = posO->m[i][2];
 
  return 0;
}

float isHbond(float d[3],float h[3],float a[3],
              float rcut,float acut,float L[3])
{    
  int i,j;
  float r_da[3],r_dh[3];
  double CosAngle;

  acut = cos(acut/180.0*3.1415926); // convert from angle to cos(angle)
  
  // compute vector: donor->acceptor and hydrogen->donor
  for(i=0;i<3;i++)
  {
    r_da[i] = periodicity(a[i]-d[i],L[i]); // r_da = a - d
    r_dh[i] = periodicity(h[i]-d[i],L[i]); // r_dh = h - d
    
    if(r_da[i]>rcut || r_da[i]< -rcut)
      return 2003;
  }
  
  if(i_prod(r_da,r_da)>rcut*rcut)
    return 2003;
  else
  {
    CosAngle = cos_angle_new(r_da,r_dh);
    
    if(CosAngle < acut) // small angle gives big cos(angle)
      return 2003;
    else
      return acos(CosAngle)/3.1415926*180;  
  }
}

double i_prod(float a[3],float b[3])
{
  double aa=a[0];
  double aaa=a[1];
  double aaaa=a[2];

  double bb=b[0];
  double bbb=b[1];
  double bbbb=b[2];

  double res = aa*bb+aaa*bbb+aaaa*bbbb;

  return res;
}

double dist1(float a[3], float b[3], float L[3])
{
  float r2=0, dr;
  int k;

  for(k=0; k<3; k++)
  {
    dr = periodicity(a[k]-b[k], L[k]);
    r2+= dr*dr;
  }

  return sqrt(r2);
}

float cos_angle_new(float a[3],float b[3])
{
  /* 
   *                  ax*bx + ay*by + az*bz
   * cos-vec (a,b) =  ---------------------
   *                      ||a|| * ||b||
   */
  double r_ab_2,r_a,r_b;
  double cosAng;
  
  r_ab_2 = i_prod(a,b);
  r_a  = normNEW(a);
  r_b  = normNEW(b);
  
  cosAng = r_ab_2/(r_a*r_b);
  
  if(cosAng<-(1+1e-8) || cosAng > (1.0+1e-8))
    {
      printf("Houston, we got a problem: |cos(alpha)| > 1 !, cosAng=%20.15f, %10.4f, %10.4f, %10.4f, %10.4f, %10.4f, %10.4f\n", 
             cosAng, a[0], a[1], a[2], b[0], b[1], b[2]);
      exit(0);
    }  
  else
    return cosAng;
}

int WhichBin(MMatrix_t rbin, int nbin, real dis)
{
  int i, j, k;

  for(i=0; i<nbin && dis>rbin->m[i+1][0]; i++);

  if(i>=nbin)
  {
    i=nbin-1;
  }

  return i;
}

int neighbourList(int nWt,MMatrixInt_t List,float rcut2,MMatrix_t posO,float L[3])
{
  int i,j;
  for(i=0;i<nWt;i++)
  {
    List->m[i][0] = 0;
    for(j=0;j<nWt;j++)
    {
      if(i!=j && dist2(posO,i,posO,j,L)<rcut2)
      {
        List->m[i][0]++;
        List->m[i][List->m[i][0]] = j;
      }
      if(List->m[i][0]>maxList)
      {
        printf("Number of atoms in neighbour exceeds %d\n",maxList);
        exit(0);
      }
    }
  }
  return 0;
}

int neighbourList2(MMatrixInt_t List,MMatrix_t posOS,MMatrix_t posO,
                   int nOS,int nWt,float rcut2,float L[3])
{
  int i,j;

  if(List->rows<nOS)
  {
    printf("do_xxx.c, neighbourList2, List->rows<nWt\n");
    exit(0);
  }

  for(i=0;i<nOS;i++)
  {
    List->m[i][0] = 0;
    for(j=0;j<nWt;j++)
    {
      if(dist2(posOS,i,posO,j,L)<rcut2)
      {
        List->m[i][++List->m[i][0]] = j;
      }
      if(List->m[i][0]>maxList)
      {
        printf("Number of atoms in neighbour exceeds %d\n",maxList);
        exit(0);
      }
    }
  }
  return 0;
}

float dist2(MMatrix_t posOS,int i,MMatrix_t posOW,int j,float L[3])
{
  float r2 = 0,dr;
  int k;
  for(k=0;k<3;k++)
  {
    dr = periodicity(posOS->m[i][k]-posOW->m[j][k],L[k]);
    r2+= dr*dr;
  }
  return r2;
}

double normNEW(float a[3])
{
  double aa=a[0];
  double aaa=a[1];
  double aaaa=a[2];

  double res = aa*aa+aaa*aaa+aaaa*aaaa;

  return sqrt(res);
}

int StorePosition(t_trxframe fr, t_topology top,
                  MMatrix_t posOW, atom_id *OWIndex, int size, 
		  int nWt, int start, bool bMW)
{
  int i, j, k, ll, mm, nn;
  int per = size/nWt;
  int num;

  if(nWt==0)
  {
    return 0;
  }
  if(posOW->rows<start+nWt)
  {
    printf("do_hbond.c, posOW->rows<start+nWt\n");
    exit(0);
  }
  if(posOW->cols<3)
  {
    printf("do_hbond.c, posOW->cols<3\n");
    exit(0);
  }
  if(size%nWt)
  {
    printf("do_hbond.c, size%%nWt!=0\n");
    exit(0);
  }

  if(per>1)
  {
    for(i=0; i<nWt; i++)
    {
      num=RemovePBCSep(fr, OWIndex, per, i*per);
      printf("i=%d, num=%d\n", i, num);
    }
  }

  if(bMW==FALSE)
  {
    for(i=0;i<nWt;i++)
    {
      for(k=0; k<DIM; k++)
      {
	posOW->m[i+start][k] = fr.x[OWIndex[i]][k];
      }
    }
    return 0;
  }

  real mass;
//  printf("per=%d, size=%d, nWt=%d\n", per, size, nWt);
  for(i=0; i<nWt; i++)
  {
    for(k=0; k<DIM; k++)
    {
      posOW->m[i+start][k] = 0.0;
    }
    for(j=0, mass=0.0; j<per; j++)
    {
      ll=i*per+j;
      mass += top.atoms.atom[OWIndex[ll]].m;
      for(k=0; k<DIM; k++)
      {
	posOW->m[i+start][k] += fr.x[OWIndex[ll]][k]*top.atoms.atom[OWIndex[ll]].m;
      }
    }
    for(k=0; k<DIM; k++)
    {
      posOW->m[i+start][k] /= mass;
      printf("i=%d, k=%d, %15.4f ", i, k, posOW->m[i+start][k]);
    }
//    printf("mass=%15.4f\n", mass);
  }

//  printf("size=%d, start=%d, per=%d, DIM=%d\n", size, start, per, DIM);
//  printf("mass weighted, size=%d, start=%d, per=%d, DIM=%d\n", size, start, per, DIM);

  return 0;
}

void CombinePosition(MMatrix_t pos, MMatrix_t half1, MMatrix_t half2)
{
  int i, j, k;

  if(half1->cols!=half2->cols)
  {
    printf("do_msd2.c, half1->cols!=half2->cols, %d, %d\n", half1->cols, half2->cols);
    exit(0);
  }

  if(pos->rows!=half1->rows+half2->rows)
  {
    printf("do_msd2.c, pos->rows!=half1->rows+half2->rows, %d, %d, %d\n", pos->rows, half1->rows, half2->rows);
    exit(0);
  }

  for(i=0; i<pos->rows/2; i++)
  {
    for(j=0; j<pos->cols; j++)
    {
      pos->m[i][j] = half1->m[i][j];
      pos->m[i+pos->rows/2][j] = half2->m[i][j];
    }
  }
}

void ShiftPosition(MMatrix_t posOW, int half)
{
  int i, j, k;

  if(posOW->rows!=2*half)
  {
    printf("do_msd2, ShiftPosition, posOW->rows!=half*2, %d, %d\n", posOW->rows, half);
    exit(0);
  }

  for(i=0; i<half; i++)
  {
    for(j=0; j<posOW->cols; j++)
    {
      posOW->m[i][j] = posOW->m[i+half][j];
    }
  }
}

void GetOrientation(MMatrix_t A, MMatrix_t aa)
{
  int i, j, k;

  MMatrix_t R=CreateMatrix(A->rows,1);
  for(i=0;i<R->rows;i++)
  {
    R->m[i][0]=1.0;
  }
//  PrintMatrix(R, "R");
  MMatrix_t AT=MatrixTranspose(A);
  MMatrix_t AA=MatrixProduct(AT,A);
  MMatrix_t AR=MatrixProduct(AT,R);
//  PrintMatrix(A, "A");
//  PrintMatrix(AA, "AA");
//  PrintMatrix(AR, "AR");

  EqnSolver(AA,AR,aa);
//  PrintMatrix(aa, "aa");
  MatrixNorm2(aa);
//  PrintMatrix(aa, "aa");

  FreeMatrix(AT);
  FreeMatrix(AA);
  FreeMatrix(AR);
  FreeMatrix(R);
}

void GetEEOrientation(MMatrix_t A, MMatrix_t aa)
{
  int i, j, k, ll, mm, nn;

  for(i=0; i<aa->rows; i++)
  {
    ll = A->rows-1;
    aa->m[i][0] = A->m[ll][i]-A->m[0][i];
  }
//  PrintMatrix(aa, "aa");
  MatrixNorm2(aa);
//  PrintMatrix(aa, "aa");
}

int NormVec(float a[3])
{
  int i, j, k;

  float r = normNEW(a);

  if(r>1e-6)
  {
    for(k=0; k<3; k++)
    {
      a[k] /= r;
    }
  }

  return 1;
}
