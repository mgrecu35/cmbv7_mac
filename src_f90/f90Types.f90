
module f90Types
  use iso_c_binding
  integer :: nbins,nEnsM,nc
  parameter(nbins=88)
  parameter(nEnsM=50)
  parameter(nmfreq_=8)
  parameter(nc=50)
 
  type,bind(c) :: stormStructType
     integer(c_int) :: nodes (5)    
     ! 5 nodes that define the storm structure 
     ! for stratiform, with BB, nodes[0] is the storm top, 
     ! nodes[1] is the BB top,
     ! nodes[2] is the BB peak, nodes[3] is the BB bottom, 
     ! and nodes[4] is the lowest
     ! clutter free gate (0<=nodes[0]<nodes[1]<nodes[2]
     ! <nodes[3]<nodes[4]<ngates)
     ! for convective, nodes[3] is not used
     ! the C convention is used the gates are 
     ! numbered from 0 to ngates-1
     integer(c_int)   ::  iSurf      ! surface gate number
     real(c_float)    ::  freezH     ! freezing height - not used
     integer(c_int)   ::  rainType;  ! rainType -1 stratiform
  end type stormStructType

  type,bind(c):: radarDataType
     integer(c_int)                    :: ngates     ! number of gates
     real(c_float)  :: z13obs(nbins)     ! observed Ku-reflectivities
     real(c_float)  :: z35obs(nbins)     ! observed Ka-reflectivities
     real(c_float)                       :: xlong      ! longitude -not used
     real(c_float)                       :: xlat       ! latitude -not used
     real(c_float)                       :: pia13srt   ! Ku-band SRT PIA 
     real(c_float)                       :: relPia13srt   ! Ku-band SRT PIA 
     real(c_float)                       :: pia35srt   ! Ka-band SRT PIA
     real(c_float)                       :: relPia35srt   ! Ka-band SRT PIA
     real(c_float)                       :: dr         ! gate size
     real(c_float)  :: hh(nbins);        ! hh[i] is the height of gate [i]
     real(c_float)                       :: hfreez
     real(c_float)			:: sigmaZeroKu
     real(c_float)			:: sigmaZeroKa
  end type radarDataType

  type,bind(c)::  radarRetType
    
     integer(c_int)   ::  ngates     ! number of gates
     integer(c_int)   ::  nMemb      ! number of members
     integer(c_int)   ::  nmfreq     ! # of simulated passive microwave frequencies
     real(c_float) :: z13c(nbins*nEnsM)        ! effective reflectivity 
                                              ! (attenuation corrected)
                                              ! at Ku-band
     real(c_float) :: z35mod0(nbins*nEnsM)     ! simulated observations at Ka-band
!  SFM  begin  06/16/2014; for M. Grecu, multiple scattering
     real(c_float) :: dz35ms(nbins*nEnsM)      !  multiple scattering effect
!  SFM  end    06/16/2014
     real(c_float) :: z35(nbins*nEnsM)         ! simulated observations at Ka-band
     real(c_float) :: pwc(nbins*nEnsM)         ! precipitation water content 
                                              ! (g/m3) 
     real(c_float) :: rrate(nbins*nEnsM)       ! rain rate (mm/h)
     real(c_float) :: d0(nbins*nEnsM)          ! median diameter (mm)
     real(c_float) :: log10dNw(nbins*nEnsM)    !
     real(c_float) :: tb(nmfreq_*2*nEnsM)          ! simulated brightness temperatures
     real(c_float) :: emTb(nmfreq_*2*nEnsM)          ! simulated brightness temperatures
     real(c_float) :: emis(nmfreq_*2*nEnsM)          ! emissivity
     integer(c_int)  ::  imu(nbins*nEnsM)    ! mu index of the look up table
     integer(c_int)  ::  iwc(nc)    
     integer(c_int)  ::  icc(nc)        
                                              ! index of RH profile (from 1 to nc)
                                              ! nc is the number of possible 
                                              ! RH profiles see cloud.f90
     integer(c_int)   ::  jcc(nc)       
                                              ! index of cloud profile (from 1 to nc)
     real(c_float)      ::  sfc_wind(nEnsM), sfc_windU(nEnsM), sfc_windV(nEnsM)   
  
     real(c_float)     ::  pia13(nEnsM)
     real(c_float)     ::  pia35(nEnsM)
     real(c_float)     ::  simSigmaZeroKu(nEnsM)
     real(c_float)     ::  simSigmaZeroKa(nEnsM)
     real(c_float)     ::  z35mMean(nbins)
     real(c_float)     ::  z35mSig(nbins)
     real(c_float)                           ::  pia13mMean(nEnsM), pia35mMean(nEnsM), pia13mSig(nEnsM), pia35mSig(nEnsM)
     integer(c_int)    ::  ntpw
     real(c_float)     ::  tpw(300)
     real(c_float)     ::  tpwCldMod(nc)
     real(c_float)     ::  logdNw(9*nEnsM)
                                             ! the last four variables are not retrieved
                                             ! they are randomly set and 
                                             ! specify the conditions 
     ! in which the retrievals and 
     ! associated brightness temperatures are derived
  end type radarRetType

  type  retParamType
     
  !
  !  F=wz*SUM(zsim,ka-zobs,ka)**2+w13*(pia13-pia13srt)**2+
  !  w35*(pia35,sim-pia35srt)**2
  !  wz the weight of the reflectivity term
  !  w13 the weight of the pia squared difference at ku-band
  !  w35 the weight of the pia squared difference at ka-band
  !  z13thresh - threshold (dBZ) to determine 
  !  the Ku-band observations used in the algorithm
  !  z35thresh - threshold (dBZ) to determine the Ka-band observations 
  !  used in the algorithm
  
  
     real :: wz, w13,  w35,  z13thresh,  z35thresh
  end type retParamType


end module f90Types
