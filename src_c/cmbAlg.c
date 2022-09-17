int mainj (int i, char *pfname, int *ifs, char* jobname,int *ialgOut);
#include "stdlib.h"
#include <stdio.h>
#include <string.h>
void do_chunk_(int *i,int *one, int *idir);
void do_chunkx_(int *i,int *one, int *idir);
void dealloc_chunk_(int *i);
void radarretsub2_(int *nmu2,  int *nmfreq2,   int *icL, float *tbRgrid,  
		  float *dprrain, int *ichunk, int *orbNumb, int *ialg, int *idir);
void radarretsub3_(int *nmu2,  int *nmfreq2,   int *icL, float *tbRgrid,  
		  float *dprrain, int *ichunk, int *orbNumb, int *ialg, int *idir);
void radarretsub4_(int *nmu2,  int *nmfreq2,   int *icL, float *tbRgrid,  
		   float *dprrain, int *ichunk, int *orbNumb, int *ialg, int *idir, int *nscans_c);
void radarretsub4_fs_(int *nmu2,  int *nmfreq2,   int *icL, float *tbRgrid,  
		      float *dprrain, int *ichunk, int *orbNumb, int *ialg, int *idir,
		      int *nscans_c);
void dealloc_struct_(int *i);
void closefiles_(int *i);
void rewindc_(int *ic);
void writescan_fs_300_(int *isc);
void writescan_300_(int *isc);
void writescant_300_(int *isc);
//begin 2/9/22 WSO declare status variable for GMI/TMI
#define iConus __globaldata_MOD_iconus
int st_GMITMI2_c;
//end 2/9/22 WSO
int main(int argc, char *argv[])
{
  char  fname[100];
  char jobname[255];
  int ifs;
//begin 2/9/22 WSO declare external variable in main
  extern int st_GMITMI2_c;
  extern int __globaldata_MOD_iconus;
//  end 2/9/22
  if(argc != 3)
    {fprintf(stderr, 
	     "\nCommand Line ERROR-should be 2 arguments (jobname, parameterFile)\n");
      exit(1);}
  
  strcpy(jobname, argv[1]);
  
  strcpy(&fname[0],argv[2]);
  printf("%s \n",&fname[0]);
  int ialg;
  int ndpr=mainj(1,fname,&ifs,&jobname[0],&ialg);
  printf("Back from mainj() %i %i\n",ndpr,ialg);
  printf("value of st_GMITMI2_c: %i\n", st_GMITMI2_c);
  
  int ny=49;
  int nx=300;
  int nz=88;
  printf("%i %i\n",ndpr, ifs);
  if(ndpr<0)
    exit(0);
  int i,one=1;
  int nmu=5, nmfreq=8, orbNumb=0;
  float *dprrain, *tbRgrid;
  
  tbRgrid=(float*) malloc(sizeof(float)*9300*49*14);
  dprrain=(float*) malloc(sizeof(float)*49*300);
 
  int idir;
  int icL;
  int nchunk=ndpr/300;
  //nchunk=6;
  //printf("nchunk = %d\n",nchunk);
  //printf("here \n");
  if(ndpr>0)
    {
      for(i=0;i<=nchunk;i++)
	{
	  if(ifs==1)
	    {
	      printf("Calling do_chunkx()\n");
	      do_chunkx_(&i,&one,&idir);
	    }
	  else
	    do_chunk_(&i,&one,&idir);
	  icL=i*300;
	  printf("iConus=%i\n",iConus);
	  iConus=0;
	  if(iConus==1)
	    {
	      printf("GMIretsub\n");
	      //gmiretsub_(&icL, &i, &orbNumb, &ialg, &idir, &st_GMITMI2_c);
	      printf("Calling radarretsub2\n");
	      radarretsub2_(&nmu,  &nmfreq,   &icL, tbRgrid,  
			    dprrain, &i, &orbNumb, &ialg, &idir);
	      printf("Calling radarretsub3\n");
	      radarretsub3_(&nmu,  &nmfreq,   &icL, tbRgrid,  
			    dprrain, &i, &orbNumb, &ialg, &idir);
	      if(ifs==1) 
		{
		  int nscans_c;
		  printf("Calling radarretsub4_fs() \n");
		  radarretsub4_fs_(&nmu,  &nmfreq,   &icL, tbRgrid,  
				   dprrain, &i, &orbNumb, &ialg, &idir, &nscans_c);
		  int j=0;
		  /*for(j=0;j<nscans_c;j++)
		    writescan_fs_300_(&j);*/
		}
	      else
		{
		  int nscans_c;
		  printf("Calling radarretsub4()\n");
		  radarretsub4_(&nmu,  &nmfreq,   &icL, tbRgrid,  
				dprrain, &i, &orbNumb, &ialg, &idir, &nscans_c);
		  int j=0;
		  /*if(ialg==1)
		    for(j=0;j<-nscans_c;j++)
		    writescan_300_(&j);
		    else
		    for(j=0;j<-nscans_c;j++)
		    writescant_300_(&j);
		  */
		  
		}
	      dealloc_struct_(&i);
	    }
	  
	  //rewindc_(&icL);
	  dealloc_chunk_(&i);
	}
      printf("Closing Files \n");
      closefiles_(&one);
    }
}

