#include "Mgromacs.h"
#include "Mmatrix.h"

int storePosition(t_trxframe, double **, int *, int);

void do_rdf(char *fnTPS, char *fnTRX, char *fnRDF, 
            int ng, int *isize, char **grpname, atom_id **index,
	    real bin, real Lx, real Ly, real Lz, int natom, int snp)
{
  FILE       *RDF;
  int        status;
  int        g,natoms,nbin,j0,j1,n,nframes;
  int        **count;
  int        isize_cm=0,nrdf=0,max_i;
  atom_id    *index_cm=NULL;
  unsigned long int *sum;
  t_topology top;
  t_block    *excl;
  t_trxframe fr;
  excl=NULL;
  int i, j, k;
  
  real L[3], minL=0.0;
  L[0] = Lx; 
  L[1] = Ly; 
  L[2] = Lz;
  for(i=1, minL=L[0]; i<3; i++)
  {
    if(minL>L[i])
    {
      minL = L[i];
    }
  }
  nbin = (int)(minL/bin/2);
  real rSpacing = bin;
  
  double  **posOS, **posOW;
  double  *rdfCount;
  int *OSIndex, *OWIndex, nWa, nWt;
  rdfCount = CreateVector(nbin+1);
  posOS = CreateMatrix(isize[0]+1,3);
  posOW = CreateMatrix(isize[1]+1,3);
  for(i=0;i<nbin;i++)
  {
    rdfCount[i] = 0.0;
  }
  OSIndex = index[0];
  OWIndex = index[1];
  nWa = isize[0];
  nWt = isize[1];

  double r, diff_r[3];

  int flags = TRX_READ_X;
  read_first_frame(&status, fnTRX, &fr, flags);
  int step, me;
  for(step=0; step<snp; step++)
  {
    read_next_frame(status,&fr);
    storePosition(fr, posOS, OSIndex, nWa);
    storePosition(fr, posOW, OWIndex, nWt);

    for(i=0; i<nWa; i++)
    {
      for(j=0; j<nWt; j++)
      {
	if(i/natom==j/natom)
	  continue;
	r = 0.0;
	//Remove PBC
	for(k=0; k<3; k++)
	{
	  diff_r[k] = posOS[i][k] - posOW[j][k];
	  diff_r[k] -= L[k]*nint(diff_r[k]/L[k]); 
	  r += diff_r[k]*diff_r[k];
	}
        r = sqrt(r);
        me = (int)(r/rSpacing);
        if(me>=0 && me <nbin) rdfCount[me] += 1.0;
      }
    }
  }

  //report
  real volume, density;
  for(i=0;i<nbin;i++)
  {
    volume = rSpacing*4*3.1415926*(i+0.5)*rSpacing*(i+0.5)*rSpacing;
//    volume = rSpacing*2*3.1415926*(i+0.5)*rSpacing*Lz;
    density = nWt/Lx/Ly/Lz;
//    rdfCount[i] /= (nWa*(snp-1)*density*volume);
    rdfCount[i] /= (nWa*(snp)*density*volume);
  }
  RDF = fopen(fnRDF, "w");
  for(i=0;i<nbin;i++)
  {
    fprintf(RDF, "%8.6f %8.6f\n",(i+0.5)*rSpacing, rdfCount[i]);
    printf("%6.2f %6.2f\n",(i+0.5)*rSpacing, rdfCount[i]);
  }
  fclose(RDF);
}

int storePosition(t_trxframe fr, double **positionOS, int *OSIndex, int nWa)
{
  int i, j, k;

  for(i=0;i<nWa;i++)
  {
    positionOS[i][0] = fr.x[OSIndex[i]][0];
    positionOS[i][1] = fr.x[OSIndex[i]][1];
    positionOS[i][2] = fr.x[OSIndex[i]][2];
  }
  return 0;
}
