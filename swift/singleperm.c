/*
singleperm.c 

generates S matrices for residuals and A matrices for 
one-way symmetric connections for all possible models 
of the given size with the given number of connections. 
this can be called on it's own or by an R script.

skenny@uchicago.edu 2009


*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int PERMCOUNT = 0;
static int ASSYMCOUNT = 0;
static int STATCOUNT = 0;
static int data[4] = {0,0,0,0};


/*  Free memory used for data  */

void FreeData(void) {
  return;
}


int factorial( int n )
{
  if ( n <= 1 )
    return 1;
  else
    return  n * factorial( n-1 );
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

int get_stats(int connections, int idata[], int size, int recentdata[])
{
  int r = 0;
  int ch = 0;
  printf("[%3i] ", STATCOUNT);
  while(r<size)
    {
      printf("%2i", idata[r]);
      r++;
    }
  printf("\n");
  STATCOUNT++;
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

int find_last(int size, int idata[])
  {
    int j;
    int last;
  int ct = 0;
  for(j=0; j<size; j++)
    {
       if (idata[j] == 1)
	 last = j;
    }
  return last;
  }

void Perm(int size, int connections, int idata[], int recentdata[], int rowcol, int permnum) {

  int i,j,c,m,n,p,q,r,s,scount;
  int k = 0;

  int asymm = 0;
  int mxdat_a[rowcol][rowcol];
  int mxdat_s[rowcol][rowcol];
  int zeros = size - connections;
  printf("getting total, size is %i, connections is %i, zeros is %i\n", size, connections, zeros);
  /*  int total = factorial(size)/(factorial(connections)*factorial(size-connections));*/
  int last_shift_pos = size - connections;
  /*  printf("total is %i:\n",total);*/

   /* begin permutations */

  while(find_n(1,idata,size) <= last_shift_pos)
    {
      /*      get_stats(connections, idata, size, recentdata);*/
      /* load into grid */

      while(k<size)
	{
	  for(m=0; m<rowcol; m++)
	    {

	      for(n=0; n<rowcol; n++)
		{  
		  /* check for diag symmetry */		
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
      /* check for remaining symmetry */
      int symm = 0;
      for(r=0; r<rowcol; r++)
	{
	  for(s=0; s<rowcol; s++)
	    {
	      if(mxdat_a[r][s] == 1 && mxdat_a[s][r] == 1)
		symm = 1;
	    }
	}
      if(!symm)
	{
	  if(permnum == ASSYMCOUNT)
	    {
	      printf("permnum %i PERMCOUNT %i\n",permnum,PERMCOUNT);
	      FILE *sdatfile, *adatfile;
	      char sname [50], aname [50];
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
	    }
	  ASSYMCOUNT++;

	}
      k=0;
      PERMCOUNT++;
      if(find_n(1,idata,size) == last_shift_pos)
	{
	  printf("ASSYMCOUNT %i\n", ASSYMCOUNT);
	  printf("PERMCOUNT %i\n", PERMCOUNT);
	  break;
	}
	  insert_zero(find_last(size,idata), size, connections, idata, recentdata);
    }

}

int run_perm(int *sz, int *nc, int *total, int *pnum, int *connections) {

  int i,j;
  int numcol = nc[0];
  int size = sz[0];
  int permnum = pnum[0];
  int conn = connections[0];
  int idata[size];
  int recentdata[size];
  /*  int max_con = floor(numcol*(numcol+1)/2);*/

  for(i=0; i<size; i++)
    {
      if(i<conn)
	idata[i] = 1;
      else
	idata[i] = 0;
    }
  for(j=0; j<size; j++)
    {
      recentdata[j] = idata[j];
    }
  Perm(size, conn, idata, recentdata, numcol, permnum);
  return 0;


}

int main(int argc, char *argv[])
{

  int pnum, ncol, conn;
  printf("enter permnum, number of columns and conn:\n");
  scanf("%i %i %i", &pnum, &ncol, &conn);
  int size = ncol*ncol;
  int numcol = ncol;
  int idata[size];
  int recentdata[size];
  /*  int max_con = floor(numcol*(numcol+1)/2);*/

  int i,j;

  for(i=0; i<size; i++)
    {
      if(i<conn)
	idata[i] = 1;
      else
	idata[i] = 0;
    }
  for(j=0; j<size; j++)
    {
      recentdata[j] = idata[j];
    }
  Perm(size, conn, idata, recentdata, ncol, pnum);
  return 0;
}
