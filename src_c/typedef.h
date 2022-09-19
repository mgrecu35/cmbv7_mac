#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef GFOR 
#define nb 5
#endif
#ifdef IFORT 
#define nb 8
#endif

typedef struct stormStructType
{
  int     nodes[5];     // 5 nodes that define the storm structure 
                      // for stratiform, with BB, nodes[0] is the storm top, 
                      // nodes[1] is the BB top,
                      // nodes[2] is the BB peak, nodes[3] is the BB bottom, 
                      // and nodes[4] is the lowest
                      // clutter free gate (0<=nodes[0]<nodes[1]<nodes[2]
                      // <nodes[3]<nodes[4]<ngates)
                      // for convective, nodes[3] is not used
                      // the C convention is used the gates are 
                      // numbered from 0 to ngates-1
  int     iSurf;      // surface gate number
  float   freezH;     // freezing height - not used
  int     rainType;   // rainType -1 stratiform
                      //          -2 convective
} stormStructType;

typedef struct radarDataType
{
  int     ngates;     // number of gates
  float   z13obs[88];    // observed Ku-reflectivities
  float   z35obs[88];    // observed Ka-reflectivities
  float   xlong;      // longitude -not used
  float   xlat;       // latitude -not used
  float   pia13srt;   // Ku-band SRT PIA 
  float   relPia13srt;   // Ku-band SRT PIA 
  float   pia35srt;   // Ka-band SRT PIA
  float   relPia35srt;   // Ka-band SRT PIA
  float   dr;         // gate size
  float   hh[88];        // hh[i] is the height of gate [i]
  float   hfreez;
  float   sigmaZeroKu; //SJM 12/3/2014
  float   sigmaZeroKa; //SJM 12/3/2014
} radarDataType;

#define nbins 88
#define nEnsM  50
#define nc 50
#define nmfreq_ 8
typedef struct radarRetType
{

  int     ngates;     // number of gates
  int     nMemb;
  int     nmfreq;     // # of simulated passive microwave frequencies
  float   z13c[nbins*nEnsM];      // effective reflectivity (attenuation corrected)
                      // at Ku-band
  float   z35mod0[nbins*nEnsM];   // simulated observations at Ka-band
  float   dz35ms[nbins*nEnsM];   // simulated observations at Ka-band
  float   z35[nbins*nEnsM];   // simulated observations at Ka-band
  float   pwc[nbins*nEnsM];       // precipitation water content (g/m3) 
  float   rrate[nbins*nEnsM];     // rain rate (mm/h)
  float   d0[nbins*nEnsM];        // median diameter (mm)
  float   log10dNw[nbins*nEnsM];  // 
  float   tb[nmfreq_*2*nEnsM];        // simulated brightness temperatures
  float   emTb[nmfreq_*2*nEnsM];      // simulated emission brightness temperatures
  float   emis[nmfreq_*2*nEnsM];      // simulated emissivity
  int     imu[nbins*nEnsM];       // mu index of the look up table
  int     iwc[nc];       // index of surface wind profile (from 1 to nc)
  int     icc[nc];       // index of RH profile (from 1 to nc)
                      // nc is the number of possible RH profiles see cloud.f90
  int     jcc[nc];       // index of cloud profile (from 1 to nc)
                      // nc is the number of possible RH profiles see cloud.f90
  float  sfc_wind[nEnsM];   // surface wind speed
  float  sfc_windU[nEnsM];   // surface wind speed U
  float  sfc_windV[nEnsM];   // surface wind speed V
  float  pia13[nEnsM];   // surface wind speed
  float  pia35[nEnsM];   // surface wind speed
  float  simSigmaZeroKu[nEnsM];   // simulated backscatter
  float  simSigmaZeroKa[nEnsM];   // simulated backscatter
  float  z35mMean[nbins];   // surface wind speed
  float  z35mSig[nbins];    // surface wind speed
  float  pia13mMean, pia35mMean, pia13mSig, pia35mSig;    // surface wind speed
  int    ntpw;
  float  tpw[300];
  float  tpwCldMod[nc];
  float  logdNw[9*nEnsM];
 
} radarRetType;


typedef struct retParamType
{
  /*
    F=wz*SUM(zsim,ka-zobs,ka)**2+w13*(pia13-pia13srt)**2+w35*(pia35,sim-pia35srt)**2
    wz the weight of the reflectivity term
    w13 the weight of the pia squared difference at ku-band
    w35 the weight of the pia squared difference at ka-band
    z13thresh - threshold (dBZ) to determine the Ku-band observations used in the algorithm
    z35thresh - threshold (dBZ) to determine the Ka-band observations used in the algorithm
  */
  
  float wz, w13,  w35,  z13thresh,  z35thresh;
} retParamType;
