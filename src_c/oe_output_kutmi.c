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

void copy_oeqv_kutmi_(float *OEQv,short *n10,int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
for(k=0;k<10;k++)
{
swath_t.KuTMI.OptEst.OEvaporDensity[*i][k]=OEQv[n10[k]-1];
}
}

void copy_oetemp_kutmi_(float *OETemp,short *n10,int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
for(k=0;k<10;k++)
{
swath_t.KuTMI.OptEst.OEairTemperature[*i][k]=OETemp[n10[k]-1];
}
}

void copy_oecloud_kutmi_(float *OECloud, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
for(k=0;k<nbins;k++)
{
swath_t.KuTMI.OptEst.OEcloudLiqWaterCont[*i][k]=OECloud[k];
}
}

void copy_oecloudliqpath_kutmi_(float *OECloudLiqPath, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
swath_t.KuTMI.OptEst.OEcolumnCloudLiqWater[*i]=*OECloudLiqPath;
}

void copy_oecloudliqsigma_kutmi_(float *OEcloudLiqSigma, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
swath_t.KuTMI.OptEst.OEcolumnCloudLiqSigma[*i]=*OEcloudLiqSigma;
}

void copy_oesfcwind_kutmi_(float *OESfcWind, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
swath_t.KuTMI.OptEst.OEtenMeterWindSpeed[*i]=*OESfcWind;
}

void copy_oesfcwindsigma_kutmi_(float *OESfcWindSigma, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
swath_t.KuTMI.OptEst.OEtenMeterWindSigma[*i]=*OESfcWindSigma;
}

void copy_oeskntemp_kutmi_(float *OESknTemp, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
swath_t.KuTMI.OptEst.OEskinTemperature[*i]=*OESknTemp;
}

void copy_oeskintempsigma_kutmi_(float *OEskinTempSigma, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
swath_t.KuTMI.OptEst.OEskinTempSigma[*i]=*OEskinTempSigma;
}

void copy_oesfctemp_kutmi_(float *OESfcTemp, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
swath_t.KuTMI.OptEst.OEsurfaceAirTemperature[*i]=*OESfcTemp;
}

void copy_oesfcqv_kutmi_(float *OESfcQv, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
swath_t.KuTMI.OptEst.OEsurfaceVaporDensity[*i]=*OESfcQv;
}

void copy_oetpw_kutmi_(float *OEtpw, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
swath_t.KuTMI.OptEst.OEcolumnWaterVapor[*i]=*OEtpw;
}

void copy_oetpwsigma_kutmi_(float *OEtpwSigma, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
swath_t.KuTMI.OptEst.OEcolumnVaporSigma[*i]=*OEtpwSigma;
}

void copy_oechisq_kutmi_(float *OEchiSq, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
swath_t.KuTMI.OptEst.OEerrorOfDataFit[*i]=*OEchiSq;
}

void copy_oepianonrain_kutmi_(float *OEpiaNonRain, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
for(k=0;k<2;k++)
{
swath_t.KuTMI.OptEst.OEpiaNoPrecip[*i]=OEpiaNonRain[0];
}
}

void copy_oesimtbnonrain_kutmi_(float *OEsimTbNonRain, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
for(k=0;k<13;k++)
{
swath_t.KuTMI.OptEst.OEsimulatedBrightTemp[*i][k]=OEsimTbNonRain[k];
}
}

void copy_oeemis_kutmi_(float *OEemis, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
for(k=0;k<13;k++)
{
swath_t.KuTMI.OptEst.OEsurfEmissivity[*i][k]=OEemis[k];
}
}

void copy_oeemissigma_kutmi_(float *OEemisSigma, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
for(k=0;k<13;k++)
{
swath_t.KuTMI.OptEst.OEsurfEmissSigma[*i][k]=OEemisSigma[k];
}
}

void copy_oeemisa_kutmi_(float *OEemisA, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
for(k=0;k<13;k++)
{
swath_t.KuTMI.OptEst.OEemissivityAvgKernel[*i][k]=OEemisA[k];
}
}

void copy_oestype_kutmi_(int *OEstype, int *i)
{int k;
extern L2BCMBT_SWATHS swath_t;
swath_t.KuTMI.OptEst.OEsurfaceClass[*i]=*OEstype;
}

