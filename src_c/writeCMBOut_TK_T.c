//  SFM 04/06/2013  Code module added in merge from M.Grecu's code
//  SFM 05/06/2013  Modifications from LW to facilitate using job names
//  SFM 06/27/2013  Parameter name changes from W.Olson; reduce unused code
//  SFM 07/19/2013  Large volume of code added for M.Grecu
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "TKheaders.h"
#include "TK_2BCMB.h"
#ifdef GFOR 
extern int __nbinmod_MOD_imemb;
#define nbins __nbinmod_MOD_imemb
//begin  WSO 9/15/13 
extern float __missingmod_MOD_missing_r4;
#define missing_r4c __missingmod_MOD_missing_r4
extern short __missingmod_MOD_missing_i2;
#define missing_i2c __missingmod_MOD_missing_i2
extern long __missingmod_MOD_missing_i4;
#define missing_i4c __missingmod_MOD_missing_i4
extern int __nbinmod_MOD_ntransition;
#define ntransitions __nbinmod_MOD_ntransition
//end    WSO 9/15/13
#endif

#ifdef IFORT 
extern int nbinmod_mp_nbin_;
//begin  WSO 8/8/13
extern int nbinmod_mp_ntransition_;
//end    WSO 8/8/13
#define nbins nbinmod_mp_nbin_
//begin  WSO 8/8/13
#define ntransitions nbinmod_mp_ntransition_
//end    WSO 8/8/13
//begin  WSO 9/15/13
extern float missingmod_mp_missing_r4_;
#define missing_r4c missingmod_mp_missing_r4_
extern short missingmod_mp_missing_i2_;
#define missing_i2c missingmod_mp_missing_i2_
extern long  missingmod_mp_missing_i4_;
#define missing_i4c missingmod_mp_missing_i4_
//end    WSO 9/15/13
#endif

//begin WSO 04/07/2013
//Note that there were many structure/variable name changes in this
//version to be compatible with TKIO 3.50.8
//All S1 and S2 were changed to NS and MS, respectively
//The variable ending "Out" was removed because a separate Input structure
//was created
//end WSO 04/07/2013

extern TKINFO dprtkfile;
TKINFO ctkfile;
TKINFO ctkfileIn;

L2BCMB_SWATHS swath;
L2ADPR_SWATHS dprswath;
L2ADPR_SWATHS dprxswath;
L2BCMB_SWATHS swath1;
L2AKu_FS     L2AKuDataX;



void setlatlon_t_(float *lat, float *lon, float *sfcPrecip, 
                  float *sfcPrecipStd, float *piaOut)
{
  int i;
  extern L2BCMBT_SWATHS swath_t;

  for(i=0;i<49;i++)
    {
      swath_t.KuTMI.Latitude[i]=lat[i];
      swath_t.KuTMI.Longitude[i]=lon[i];
      if(swath_t.KuTMI.Longitude[i]>180)
	swath_t.KuTMI.Longitude[i]-=360;
      swath_t.KuTMI.nearSurfPrecipTotRate[i]=sfcPrecip[i];
      swath_t.KuTMI.nearSurfPrecipTotRateSigma[i]=sfcPrecipStd[i];
      swath_t.KuTMI.pia[i]=piaOut[i];
    }
}


void copy_oe_missing_staff_t_(int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
swath_t.KuTMI.OptEst.OEestimSurfPrecipLiqRate[*i]=missing_r4c;
swath_t.KuTMI.OptEst.OEestimSurfPrecipTotRate[*i]=missing_r4c;
swath_t.KuTMI.OptEst.OEestimSurfPrecipTotRateSigma[*i]=missing_r4c;
}


void copy_oecloudicepath_kutmi_(float *OECloudLiqPath, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
swath_t.KuTMI.OptEst.OEcolumnCloudIceWater[*i]=missing_r4c;
}



void copyrrate_t_(float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  for(k=0;k<nbins;k++)
    {
      swath_t.KuTMI.precipTotRate[*i][k]=rrate[k];
      swath_t.KuTMI.precipTotRateSigma[*i][k]=rratestd[k];
    }
}

void copysflfract_t_(float *lfract, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  swath_t.KuTMI.nearSurfPrecipLiqRate[*i]=*lfract;
  
}


//begin  WSO 8/30/13
void copyenvsfqv_t_(float *envQv, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  swath_t.KuTMI.surfaceVaporDensity[*i]=envQv[nbins-1];
}



void copyenvqv_t_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  for(k=0;k<10;k++)
      swath_t.KuTMI.vaporDensity[*i][k]=envQv[envnodes[k]-1];
}


void copyenvpress_t_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  for(k=0;k<10;k++)
    swath_t.KuTMI.airPressure[*i][k]=envQv[envnodes[k]-1];
}

void copyenvtemp_t_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  for(k=0;k<10;k++)
     {
      swath_t.KuTMI.envParamNode[*i][k]=envnodes[k]-1;
      swath_t.KuTMI.airTemperature[*i][k]=envQv[envnodes[k]-1];
     }
}



void copyenvsftemp_t_(float *envQv, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  swath_t.KuTMI.surfaceAirTemperature[*i]=envQv[nbins-1];
}



//end    WSO 8/30/13

void copypwc_t_(float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  for(k=0;k<nbins;k++)
    {
      swath_t.KuTMI.precipTotWaterCont[*i][k]=rrate[k];
      swath_t.KuTMI.precipTotWaterContSigma[*i][k]=rratestd[k];
    }
}

//begin  WSO 8/7/13
void copylwcfrac_t_(float *mlwc_frac, float *mrate_frac, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  for(k=0;k<ntransitions;k++)
    {
      //swath_t.KuTMI.liqMassFracTrans[*i][k]=mlwc_frac[k];
      //swath_t.KuTMI.liqRateFracTrans[*i][k]=mrate_frac[k];
    }
}

void copysfcrainliqfrac_t_(float *sfcrainliq_frac, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  swath_t.KuTMI.nearSurfPrecipLiqRate[*i]=*sfcrainliq_frac;

}
//end    WSO 8/7/13

void copyd0_t_(float *dm, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  for(k=0;k<nbins;k++)
    {
      swath_t.KuTMI.precipTotDm[*i][k]=dm[k];
    }
}

void copyd0_a_t_(float *dm, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  for(k=0;k<nbins;k++)
    {
      swath_t.KuTMI.aPriori.initPrecipTotDm[*i][k]=dm[k];
    }
}


void copynw_t_(float *Nw, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  for(k=0;k<nbins;k++)
    {
      swath_t.KuTMI.precipTotLogNw[*i][k]=Nw[k];
      swath_t.KuTMI.precipTotMu[*i][k]=2.0;
    }
}

void copynw_a_t_(float *Nw, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  for(k=0;k<nbins;k++)
    {
      swath_t.KuTMI.aPriori.initPrecipTotLogNw[*i][k]=Nw[k];
    }
}


void copyzcku_t_(float *zc, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

   for(k=0;k<nbins;k++)
    {
      if(zc[k] > -90.)
        swath_t.KuTMI.correctedReflectFactor[*i][k] = zc[k];
      else
//begin  WSO 9/17/13 standardized missing flags
        swath_t.KuTMI.correctedReflectFactor[*i][k] = missing_r4c;
//end    WSO 9/17/13
    }
}

void copynodes_t_(int *node, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  for(k=0;k<5;k++)
    {
      swath_t.KuTMI.phaseBinNodes[*i][k]=node[k];
    }
}



void copy_tot_to_liqwatercont_ku_t_(int *ij, float *precipTot, int *node)
{
  int i;
  extern L2BCMBT_SWATHS swath_t;
  for(i=node[0];i<=node[1];i++)
    swath_t.KuTMI.precipLiqWaterCont[*ij][i]=0.0;
  int n3=node[3];
  if(n3>node[4])
    n3=node[4];
  for(i=node[1]+1;i<=n3;i++)
    if(precipTot[i]>=-1e-9)
      {
	float f=(i-node[1])/(node[3]-node[1]+1e-3);
	swath_t.KuTMI.precipLiqWaterCont[*ij][i]=f*precipTot[i];
      }
    else
      swath_t.KuTMI.precipLiqWaterCont[*ij][i]=missing_r4c;
  int node4=node[4]+1;
  if(node4>87)
    node4=87;
  for(i=node[3]+1;i<=node4;i++)
    if(precipTot[i]>=-1e-9)
      swath_t.KuTMI.precipLiqWaterCont[*ij][i]=precipTot[i];
    else
      swath_t.KuTMI.precipLiqWaterCont[*ij][i]=missing_r4c;
  for(i=node[4]+1;i<88;i++)
    swath_t.KuTMI.precipLiqWaterCont[*ij][i]=missing_r4c;
}

void copy_tot_to_liqrate_ku_t_(int *ij, float *precipTot, int *node)
{
  int i;
  extern L2BCMBT_SWATHS swath_t;
  for(i=node[0];i<=node[1];i++)
    swath_t.KuTMI.precipLiqRate[*ij][i]=0.0;
  int n3=node[3];
  if(n3>node[4])
    n3=node[4];
  for(i=node[1]+1;i<=n3;i++)
    if(precipTot[i]>=-1e-9)
      {
	float f=(i-node[1])/(node[3]-node[1]+1e-3);
	swath_t.KuTMI.precipLiqRate[*ij][i]=f*precipTot[i];
      }
    else
      swath_t.KuTMI.precipLiqRate[*ij][i]=missing_r4c;
  int node4=node[4]+1;
  if(node4>87)
    node4=87;
  for(i=node[3]+1;i<=node4;i++)
    if(precipTot[i]>=-1e-9)
      	swath_t.KuTMI.precipLiqRate[*ij][i]=precipTot[i];
    else
      swath_t.KuTMI.precipLiqRate[*ij][i]=missing_r4c;
  for(i=node[4]+1;i<88;i++)
    swath_t.KuTMI.precipLiqRate[*ij][i]=missing_r4c;
}


void rewind_t_(int *ic)
{
  extern TKINFO       granuleHandle2AKu;
  int status = TKseek(&granuleHandle2AKu, *ic, TK_ABS_SCAN_OFF); 
}


void frominput_t_(long *st_2adpr, int *flagScanPattern)
{
//  SFM  begin  12/13/2013
  extern TKINFO       granuleHandle2AKu;
  extern L2AKu_FS        L2AKuDataX;
  extern L2BCMBT_SWATHS swath_t;
//begin  WSO 9/1/13
  extern L2ADPR_SWATHS dprswath;
  extern L2ADPR_SWATHS dprxswath;
//end    WSO 9/1/13
  int j;
  int status, status_dpr ;
//for diagnostic
  float dummyPIA[49];
  int k, printPIA[49];
//end for diagnostic
//

//  SFM  begin  12/13/2013; add conditional to dpr read
  status=TKreadScan(&granuleHandle2AKu,&L2AKuDataX);
 
//  SFM  begin  12/13/2013
  
  
  *flagScanPattern=0;
  for( j=0; j<49; j++)
    {
      //swath_t.KuTMI.Input.piaEffective[j]=L2AKuData.SRT.pathAtten[j];
      swath_t.KuTMI.Input.piaEffective[j]=L2AKuDataX.SRT.PIAhybrid[j];  //MG  7/31/18, use hybrid PIA
//begin  WSO 9/5/13 remove flag assignment
//       swath_t.KuTMI.Input.piaEffectiveSigma[j]=-99;
//end    WSO 9/5/13
   //   swath_t.KuTMI.Input.piaEffectiveReliabFlag[j]=
	// L2AKuData.SRT.reliabFlag[j];
      swath_t.KuTMI.Input.piaEffectiveReliabFlag[j]=
	L2AKuDataX.SRT.reliabFlagHY[j];                              //WSO  8/2/18 use hybrid flag
      swath_t.KuTMI.Input.precipitationType[j]=
	L2AKuDataX.CSF.typePrecip[j];
      swath_t.KuTMI.Input.precipTypeQualityFlag[j]=
	L2AKuDataX.CSF.qualityTypePrecip[j];
      swath_t.KuTMI.Input.surfaceElevation[j]=L2AKuDataX.PRE.elevation[j];
      swath_t.KuTMI.Input.localZenithAngle[j]=L2AKuDataX.PRE.localZenithAngle[j];
      swath_t.KuTMI.Input.surfaceType[j]=L2AKuDataX.PRE.landSurfaceType[j];
//begin  WSO 9/28/13 use alternate rain flag that includes missing for bad scans
//      swath_t.KuTMI.Input.precipitationFlag[j]=L2AKuData.PRE.flagPrecip[j];
//end    WSO 9/28/13
      swath_t.KuTMI.Input.surfaceRangeBin[j]=(L2AKuDataX.PRE.binRealSurface[j]-1)/2; // MG 04/11/2014
      swath_t.KuTMI.Input.stormTopBin[j]=(L2AKuDataX.PRE.binStormTop[j]-1)/2; // MG 04/11/2014
      if(swath_t.KuTMI.Input.stormTopBin[j]<0)
	swath_t.KuTMI.Input.stormTopBin[j]=missing_i2c;
      swath_t.KuTMI.Input.stormTopAltitude[j]=L2AKuDataX.PRE.heightStormTop[j];
//begin  WSO 09/30/15 add one bin to the binClutterFreeBottom to temporarily compensate for
//the subtraction of one bin by the radar team
//      swath_t.KuTMI.Input.lowestClutterFreeBin[j]=
//	(L2AKuData.PRE.binClutterFreeBottom[j]-1)/2; // MG 04/11/2014
//begin  WSO 10/19/15 subtract one 125 m bin from binClutterFreeBottom, and 
//restore V3 definition of lowestClutterFreeBin as a test
//      swath_t.KuTMI.Input.lowestClutterFreeBin[j]=
//	(L2AKuData.PRE.binClutterFreeBottom[j])/2;
      swath_t.KuTMI.Input.lowestClutterFreeBin[j]=
	(L2AKuDataX.PRE.binClutterFreeBottom[j] - 2)/2;
//end    WSO 10/15/15
//end    WSO 09/30/15
//begin  WSO 9/17/13 correction for two bin average location in combined
      swath_t.KuTMI.Input.ellipsoidBinOffset[j]=
	    L2AKuDataX.PRE.ellipsoidBinOffset[j] + 0.125/2.;
//end    WSO 9/17/13
//begin  WSO 8/19/13
      swath_t.KuTMI.Input.zeroDegAltitude[j] = L2AKuDataX.VER.heightZeroDeg[j];
      swath_t.KuTMI.Input.zeroDegBin[j] = (L2AKuDataX.VER.binZeroDeg[j]-1)/2; // MG 04/11/2014
//end    WSO 8/19/13

    }

//diagnostic
//       if(L2AKuData.Latitude[24] > 30. &&  L2AKuData.Latitude[24] < 40. && L2AKuData.Longitude[24] > -165. && L2AKuData.Longitude[24] <-155.)
//         {
//           for(k=0;k<49;k++)
//             if(dummyPIA[k] < -99.)
//               {
//                 printPIA[k] = 99;
//               }
//             else
//               printPIA[k] = dummyPIA[k]*10.;
//           printf("lon: %10.2f,  ", L2AKuData.Longitude[24]);
//           for(k=0;k<49;k++)
//             printf("%2i", printPIA[k]);
//           printf("\n");
//         }
//end diagnostic

//begin  WSO 9/1/13 scanStatus variables copied from 2AKu
    swath_t.KuTMI.scanStatus.FractionalGranuleNumber =    
     L2AKuDataX.scanStatus.FractionalGranuleNumber;
    swath_t.KuTMI.scanStatus.SCorientation =
     L2AKuDataX.scanStatus.SCorientation;
    swath_t.KuTMI.scanStatus.acsModeMidScan =
     L2AKuDataX.scanStatus.acsModeMidScan;
    swath_t.KuTMI.scanStatus.dataQuality =
     L2AKuDataX.scanStatus.dataQuality;
    swath_t.KuTMI.scanStatus.dataWarning =
     L2AKuDataX.scanStatus.dataWarning;
    swath_t.KuTMI.scanStatus.geoError =
     L2AKuDataX.scanStatus.geoError;
    swath_t.KuTMI.scanStatus.geoWarning =
     L2AKuDataX.scanStatus.geoWarning;
    swath_t.KuTMI.scanStatus.limitErrorFlag =
     L2AKuDataX.scanStatus.limitErrorFlag;
    swath_t.KuTMI.scanStatus.missing =
     L2AKuDataX.scanStatus.missing;
    swath_t.KuTMI.scanStatus.modeStatus =
     L2AKuDataX.scanStatus.modeStatus;
    swath_t.KuTMI.scanStatus.operationalMode =
     L2AKuDataX.scanStatus.operationalMode;
    swath_t.KuTMI.scanStatus.pointingStatus =
     L2AKuDataX.scanStatus.pointingStatus;
    swath_t.KuTMI.scanStatus.targetSelectionMidScan =
     L2AKuDataX.scanStatus.targetSelectionMidScan;



}

void copyscantime_t_(int *i)
{
  extern L2BCMBT_SWATHS swath_t;
  extern int DayOfMonth[300], DayOfYear[300], Hour[300], MilliSecond[300],
    Minute[300], Month[300], Second[300], Year[300], SecondOfDay[300];
  extern NAVIGATION navigation[300];

 swath_t.KuTMI.ScanTime.DayOfMonth=DayOfMonth[*i];
 swath_t.KuTMI.ScanTime.DayOfYear=DayOfYear[*i];
 swath_t.KuTMI.ScanTime.Hour=Hour[*i];
 swath_t.KuTMI.ScanTime.MilliSecond=MilliSecond[*i];
 swath_t.KuTMI.ScanTime.Minute=Minute[*i];
 swath_t.KuTMI.ScanTime.Month=Month[*i];
 swath_t.KuTMI.ScanTime.Second=Second[*i];
 swath_t.KuTMI.ScanTime.Year=Year[*i];
 swath_t.KuTMI.ScanTime.SecondOfDay=SecondOfDay[*i];
 memcpy(&swath_t.KuTMI.navigation, &navigation[*i], sizeof(NAVIGATION));

//end WSO 04/07/2013
}

void copypreciptype_t_(int *ptype, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  //swath.S1.precipitationType[*i]=*ptype;
}

void copyw10_t_(float *w10, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  swath_t.KuTMI.tenMeterWindSpeed[*i]=*w10;
}

void copyw10sigma_t_(float *w10s, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  swath_t.KuTMI.tenMeterWindSigma[*i]=*w10s;
}





//  begin  SFM  12/26/2013
void write_empty_t_(void)

//    brief utility to put empty keyword into output file header
//    when needed
{
  char emptygranuletext[100];

  strcpy(emptygranuletext,"EMPTY") ;
  TKsetMetaString(&ctkfile, "FileHeader", "EmptyGranule", 
                  emptygranuletext);	  
}
//  end    SFM  12/26/2013

//  begin  SFM  11/27/2013
void writescan_t_(void)
{
  int ret;
  char emptygranuletext[100];
  extern L2BCMBT_SWATHS swath_t;
  // TKgetMetaString(&ctkfile, "FileHeader", "EmptyGranule", 
  //                emptygranuletext);	  
  // if (strncmp(emptygranuletext,"NOT_EMPTY",9) == 0)
  ret= TKwriteScan(&ctkfile,&swath_t);
}
//  end    SFM  11/27/2013

void copysfcairtemp_t_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  swath_t.KuTMI.surfaceAirTemperature[*i]=*sfcVar;
}


void copysfcairpress_t_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  swath_t.KuTMI.surfaceAirPressure[*i]=*sfcVar;
}


void copyskintemp_t_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  swath_t.KuTMI.skinTemperature[*i]=*sfcVar;
}

//write skin temperature estimate uncertainty
void copyskintempsigma_t_(float *skinsigma, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  swath_t.KuTMI.skinTempSigma[*i] = *skinsigma;
}


//write column vapor estimate uncertainty
void copycolumnvaporsigma_t_(float *colvaporsigma, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  //  swath_t.KuTMI.columnVaporSigma[*i] = *colvaporsigma;
}



//write column cloud liquid estimate uncerainty
void copycolumncloudliqsigma_t_(float *colcldsigma, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  
  //  swath_t.KuTMI.columnCloudLiqSigma[*i] = *colcldsigma;
}


//write algorithm type flag
void copyalgotype_t_(int *algotype, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  //  swath_t.KuTMI.FLG.algoType[*i] = *algotype;
}


//write error of non-raining data fit
void copyerrorofdatafit_t_(float *erroroffit, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  
  //  swath_t.KuTMI.errorOfDataFit[*i] = *erroroffit;
}

void copysfcemissout_t_(float *tbout, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swath_t.KuTMI.surfEmissivity[*i][k]=tbout[k];
    else
      swath_t.KuTMI.surfEmissivity[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swath_t.KuTMI.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissoutsigma_t_(float *tbout, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
    {
//      printf("%g ",tbout[k]);
      if(tbout[k] > 0.)
	swath_t.KuTMI.surfEmissSigma[*i][k]=tbout[k];
      else
	swath_t.KuTMI.surfEmissSigma[*i][k]=missing_r4c;
    }
//  printf(" from c\n");
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swath_t.KuTMI.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copytbout_t_(float *tbout, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
//begin  WSO 9/16/13
  for(k=0;k<9;k++)
    if(tbout[k] > -90.)
      swath_t.KuTMI.simulatedBrightTemp[*i][k]=tbout[k];
    else
      swath_t.KuTMI.simulatedBrightTemp[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels  
//  for(k=13;k<15;k++)
//    swath_t.KuTMI.simulatedBrightTemp[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
}



void copyrainflag_t_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  swath_t.KuTMI.Input.precipitationFlag[*i]=*sfcVar;
}

//begin  WSO 8/20/14 write new ioquality flags
void copyioquality_t_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  swath_t.KuTMI.FLG.ioQuality[*i]=*sfcVar;
}


void copysnowice_t_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  swath_t.KuTMI.Input.snowIceCover[*i]=*sfcVar;
}


void copysfcliqfract_t_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;
  swath_t.KuTMI.nearSurfPrecipLiqRate[*i]=*sfcVar;
}


void copycldwater_t_(float *var1d, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  for(k=0;k<nbins;k++)
    {
        swath_t.KuTMI.cloudLiqWaterCont[*i][k]=var1d[k];
    }
}

void copycldice_t_(float *var1d, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  for(k=0;k<nbins;k++)
    {
      swath_t.KuTMI.cloudIceWaterCont[*i][k]=var1d[k];
    }
}

//begin  WSO 9/5/13 new copy routine for SRT and DSRT pia effective sigma's
void copysigmapia_t_(float *sigmapia, int *i)
{
  extern L2BCMBT_SWATHS swath_t;
//diagnostic
//  printf("in writeCMBOut i: %5i,  sigmapia: %10.4f\n", *sigmapia, *i);
//end diagnostic
  swath_t.KuTMI.Input.piaEffectiveSigma[*i] = *sigmapia;
}

void copyprincomp_t_(float *princomp, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  for(k=0;k<5;k++)
    {
      //.. swath_t.KuTMI.aPriori.prinComp[*i][k] = princomp[k];
    }
}
//

//write profile class
void copyprofclass_t_(int *profclass, int *i)
{
  extern L2BCMBT_SWATHS swath_t;

  swath_t.KuTMI.aPriori.profClass[*i] = *profclass;
}

//write surface precip bias ratio
void copysurfprecipbiasratio_t_(float *biasratio, int *i)
{ 
  extern L2BCMBT_SWATHS swath_t;

  //  swath_t.KuTMI.aPriori.surfPrecipBiasRatio[*i] = *biasratio;
}


//write initial log10 of the PSD intercept
void copyinitnw_t_(float *initlogNw, int *n9, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  for(k=0;k<9;k++)
    {
     if(n9[k]>1 && n9[k]<=88)
       {
         //swath_t.KuTMI.aPriori.initNw[*i][k] = initlogNw[n9[k]-1];
       }
     else
       {
         //swath_t.KuTMI.aPriori.initNw[*i][k] = missing_r4c;
       }
    }
}


//write sub-footprint variability parameter
void copysubfootvariability_t_(float *subfoot, int *i)
{ 
  extern L2BCMBT_SWATHS swath_t;

  swath_t.KuTMI.nubfPIAfactor[*i] = *subfoot;
} 
  
//write multiple scattering flag
void copymultiscatcalc_t_(int *multiscat, int *i)
{ 
  extern L2BCMBT_SWATHS swath_t;
  
  swath_t.KuTMI.FLG.multiScatCalc[*i] = *multiscat;
}


//write multiple scattering surface parameter
void copymultiscatsurface_t_(float *multisfc, int *i)
{
  extern L2BCMBT_SWATHS swath_t;

  swath_t.KuTMI.multiScatMaxContrib[*i] = *multisfc;
}


void copysigmazero_t_(float *sigmazeroku, int *i)
{
  extern L2BCMBT_SWATHS swath_t;
  swath_t.KuTMI.Input.sigmaZeroMeasured[*i] = *sigmazeroku;
}

void copylognw_t_(float *logNw, int *n9, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  for(k=0;k<9;k++)
    {
//begin WSO 9/7/14 added upper limit for n9 in next line
      if(n9[k]>1 && n9[k]<=88)
        {
	  //swath_t.KuTMI.PSDparamLowNode[*i][k] = n9[k]-1;
	  //swath_t.KuTMI.precipTotPSDparamLow[*i][k][0]=logNw[n9[k]-1];
        }
      else
        {
	  //swath_t.KuTMI.PSDparamLowNode[*i][k] = missing_i2c;
	  //swath_t.KuTMI.precipTotPSDparamLow[*i][k][0]= missing_r4c;
        }
    }
}

void copymu_t_(float *mu_prof, int *n9,int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath_t;

  for(k=0;k<9;k++)
    {
      if(n9[k]>1)
        {
	  //         swath_t.KuTMI.precipTotPSDparamLow[*i][k][1]=mu_prof[n9[k]-1];
        }
      else
        {
	  //         swath_t.KuTMI.precipTotPSDparamLow[*i][k][1]=missing_r4c;
        }
    }
}



