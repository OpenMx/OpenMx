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

static int count = 1;

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

void Perm(int size, int connections, int idata[], int recentdata[], int rowcol, int permnum) {

  int i,j,c,m,n,p,q,scount;
  int k = 0;
  int symm = 0;
  int asymm = 0;
  int total = get_total(size, connections);
  int mxdat_a[rowcol][rowcol];
  int mxdat_s[rowcol][rowcol];

  /* get initial matrix */
  if(permnum == 0)
    {
   FILE *sdatfile, *adatfile;
   char sname [10], aname [10];
   const char sext[] = ".sdat";
   const char aext[] = ".adat";
   sprintf(sname, "%d%s", permnum,sext);
   sprintf(aname, "%d%s", permnum,aext);
   sdatfile = fopen(sname, "wt");
   adatfile = fopen(aname, "wt");
   int x,y;
   for (x=0; x<rowcol; x++)
     {
       for(y=0; y<rowcol; y++)
	 {
	   fprintf(sdatfile, "%i ",mxdat_s[x][y]);
	   fprintf(adatfile, "%i ",mxdat_a[x][y]);
	 }
       fputc('\n',sdatfile);
       fputc('\n',adatfile);
     }
   fputc('\n', sdatfile);
   fputc('\n', adatfile);
   fclose(sdatfile);
   fclose(adatfile);
   count++;
    }

   /* begin permutations */

  int last_shift_pos = size - connections;

  while(find_n(1,idata,size) != last_shift_pos)
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
		
		  if((m == n) && (idata[k] == 1))
		    {
		      mxdat_s[m][n] = idata[k];
		      mxdat_a[m][n] = 0;
		    }
		  else{
		    mxdat_a[m][n] = idata[k];
		    mxdat_s[m][n] = 0;
		  }
		    k++;
		}
	    }
	}
      if(permnum == count)
	{
	  FILE *sdatfile, *adatfile;
	  char sname [10], aname [10];
	  const char sext[] = ".sdat";
	  const char aext[] = ".adat";
	  sprintf(sname, "%d%s", count,sext);
	  sprintf(aname, "%d%s", count,aext);
	  sdatfile = fopen(sname, "wt");
	  adatfile = fopen(aname, "wt");
	  int x,y;
	  for (x=0; x<rowcol; x++)
	    {
	      for(y=0; y<rowcol; y++)
		{
		  fprintf(sdatfile, "%i ",mxdat_s[x][y]);
		  fprintf(adatfile, "%i ",mxdat_a[x][y]);
		}
	      fputc('\n',sdatfile);
	      fputc('\n',adatfile);
	    }
	  fputc('\n', sdatfile);
	  fputc('\n', adatfile);
	  fclose(sdatfile);
	  fclose(adatfile);
	}
   	  count++;
	  k=0;
    }
}   

int run_perm(int *con, int *sz, int *nc, int *total, int *pnum) {
  int numcol = nc[0];
  int connections = con[0];
  int size = sz[0];
  int permnum = pnum[0];
  int idata[size];
  int recentdata[size];

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
  Perm(size, connections, idata, recentdata, numcol, permnum);

  printf("\n%d permutations in all.\n", count);
  total[0] = count;
  FreeData();
  return 1;
    }

int main(int argc, char *argv[])
{
  printf("in main\n");
  return 0;
}
