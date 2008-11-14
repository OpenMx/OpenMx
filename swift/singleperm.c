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

static int PERMCOUNT = 0;

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

  int i,j,c,m,n,p,q,scount;
  int k = 0;
  int symm = 0;
  int asymm = 0;
  int mxdat_a[rowcol][rowcol];
  int mxdat_s[rowcol][rowcol];
  int batch_counter = 0;
  /*  int total = factorial(size)/(factorial(connections)*factorial(size-connections));*/

  printf(".....BEGINNING PERMS FOR %i connections.....\n",connections);
  /*  printf("total is %i:\n",total);*/

   /* begin permutations */

  int last_shift_pos = size - connections;

  while(find_n(1,idata,size) < last_shift_pos)
    {
      int symm = 0;
      int asymm = 0;


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

      if(permnum == PERMCOUNT)
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
  else{
      insert_zero(find_last(size,idata), size, connections, idata, recentdata);
  }
      PERMCOUNT++;
      k=0;
    }
  if(find_n(1,idata,size) == last_shift_pos)
    {

      int symm = 0;
      int asymm = 0;


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
      if(permnum == PERMCOUNT)
	{
	      printf("permnum %i PERMCOUNT %i\n",permnum,PERMCOUNT);
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
	}
      PERMCOUNT++;
    }
}

int run_perm(int *sz, int *nc, int *total, int *pnum) {

  int conn,i,j;
  int numcol = nc[0];
  int size = sz[0];
  int permnum = pnum[0];
  int idata[size];
  int recentdata[size];

  for(conn=1;conn<=size;conn++)
    {
      /* init matrices*/

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
    }
      printf("count is %i\n",PERMCOUNT);
      return 0;
}

int main(int argc, char *argv[])
{
  int size, pnum, conn, ncol, i, j;
  printf("enter size, permnum and number of columns:\n");
  scanf("%i %i %i", &size, &pnum, &ncol);

  int numcol = ncol;
  int idata[size];
  int recentdata[size];

  for(conn=1;conn<=size;conn++)
    {
      /* init matrices*/

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
      printf("count is %i\n",PERMCOUNT);

    }
      return 0;
}
