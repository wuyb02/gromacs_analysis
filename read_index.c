#include "read_index.h"

void read_index(char *fnNDX, int *ng, int *isize, atom_id **index, char **grpname)
{
  int i, j, k;
  
  rd_index(fnNDX,*ng,isize,index,grpname);

  /*
  for(i=0; i<*ng; i++)
  {
    printf("[ %s ]\n", grpname[i]);
    for(j=0; j<isize[i]; j++)
    {
      printf("%d  ", index[i][j]);
      if(!((j+1)%10))
        printf("\n");
    }
    printf("\n\n");
  }
  */
}
