/*
permall.c 

creates all symmetric and 
asymmetric permutations 
of the matrix size 'size'
with 'connections' conn


*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int count = 0;

/*  Free memory used for data  */

void FreeData(void) {
  return;
}

int find_n(int n, int idata[], int size)
{
  int j;
  int cnt = 0;
  for(j=0; j<=size; j++)
    {
      if (idata[j] == 1)
        cnt++;
      if(cnt == n)
        return j;
    }

  return j;
}

int find_shift_elmt(int connections, int idata[], int size)
{
  int x;
  for(x=connections-1; x>=1; x--)
    {
      int curr_n = find_n(x, idata, size);
      if(idata[curr_n+1] == 0)
        return curr_n;
    }
  return 0;
}

/*  Insert data 0 at position 'pos'  */

void insert_zero(unsigned pos, int size, int connections, int idata[], int recentdata[]) {
  int selmt_pos;
  if(idata[pos+1] != 0 || pos+1>=size){
    selmt_pos = find_shift_elmt(connections, idata, size);
    pos = selmt_pos;
  }
  int j,m,k,p;
  int bitcount = 0;
  int n = 0;
  /* create tmp for shifting */
  for(j=0; j<size; j++){
    recentdata[j] = idata[j];
  }
  int i = (signed) size;

  /* set position of interest */
  idata[pos] = n;

  /* set those to left which should 
     be static but keep track of them*/

  for(k=0; k<pos; k++){
    idata[k] = recentdata[k];
    if(idata[k] == 1){ bitcount++;}
  }

  /* now remaining active bits not set on left
     should be set immediately to the right */

  for(m=pos+1; m<=pos+(connections-bitcount); m++)
    {
      idata[m] = 1;
    }
  for(p=pos+connections-bitcount+1;p<size; p++)
    {
      idata[p] = 0;
    }
}

int get_total(int size, int connections)
{
  int top = (size-connections)+1;
  int statop = top;
  int total=0;
  int j;

  for(j=1;j<=statop;j++)
    {
      total += top*j;
      top--;
    }
  printf("total is %i\n", total);
  return total;
}

int find_last(int size, int idata[])
  {
    int j;
    int last;
  int count = 0;
  for(j=0; j<size; j++)
    {
       if (idata[j] == 1)
	 last = j;
    }
  return last;
  }

void Perm(int size, int connections, int idata[], int recentdata[], int rowcol, float weight) {

  int i,j,c,m,n,p,q,scount;
  int k = 0;
  int symm = 0;
  int asymm = 0;
  int total = get_total(size, connections);
  int mxdat[rowcol][rowcol];
  for(i=0; i<total-1; i++)
    {
      int symm = 0;
      int asymm = 0;
      insert_zero(find_last(size,idata), size, connections, idata, recentdata);

      /* load into grid */
      /* check for symmetry */
      while(k<size)
	{
	  for(m=0; m<rowcol; m++)
	    {

	      for(n=0; n<rowcol; n++)
		{  
		  mxdat[m][n] = idata[k];
		  if((m == n) && (mxdat[m][n] == 1))
		    symm = 1;
		  if((m != n) && (mxdat[m][n] == 1))
		    asymm = 1;
		  k++;
		}
	    }
	}
      if((symm == 1) && (asymm == 0)){
	  FILE *datfile;
	  char fname [10];
	  const char ext[] = ".sdat";
	  sprintf(fname, "%d%s", count,ext);
	  datfile = fopen(fname, "wt");
	  int x,y;
	  for (x=0; x<rowcol; x++)
	    {
	      for(y=0; y<rowcol; y++)
		{
		  fprintf(datfile, "%i ",mxdat[x][y]);
		}
	      fputc('\n',datfile);
	    }
	  fputc('\n', datfile);
	  fclose(datfile);
      }
      if((symm == 0) && (asymm == 1)){
	  FILE *datfile;
	  char fname [10];
	  const char ext[] = ".adat";
	  sprintf(fname, "%d%s", count,ext);
	  datfile = fopen(fname, "wt");
	  int x,y;
	  for (x=0; x<rowcol; x++)
	    {
	      for(y=0; y<rowcol; y++)
		{
		  fprintf(datfile, "%i ",mxdat[x][y]);
		}
	      fputc('\n',datfile);
	    }
	  fputc('\n', datfile);
	  fclose(datfile);
      }
	  count++;
	
	  k=0;
    }
}

int run_perm(int *con, int *sz, int *nc, float *wt) {
  int numcol = nc[0];
  int connections = con[0];
  int size = sz[0];
  int idata[size];
  float weight = wt[0];
  int recentdata[size];
  /* will use later when symm conn's not 
     on the diag are permitted */
  int Ssize = ((size - numcol)/2) - numcol;

  int i,j;
  for(i=0; i<size; i++)
    {
      if(i<connections)
	idata[i] = 1;
      else
	idata[i] = 0;
    }
  for(j=0; j<size; j++)
    {
      recentdata[j] = idata[j];
    }

  printf("connections: %i\n",connections);
  printf("size: %i\n",size);
 
  printf("Permutations of %u items:\n\n", connections);
  Perm(size, connections, idata, recentdata, numcol, weight);

  printf("\n%d permutations in all.\n", count);
  FreeData();
  return EXIT_SUCCESS;
    }

int main(int argc, char *argv[])
{
  printf("in main\n");
  return 0;
}
