module local_RD_var
  use f90Types
  real :: covTb(49,300,15,15), tbMax(49,300,15), tbMin(49,300,15), &
       tbMean(49,300,15)
  real :: invCovTb(49,300,15,15)
  real :: tbout2D(49,300,15), tb(49,300,15), tbNoOcean(49,300,15), &
       tbout2DNoOcean(49,300,15), tbObs(49,300,15)
  real :: dfdtb(49,300,15), rerr(15), tb0(49,300,15), fem(15) , &
       tb0MS(49,300,15), tbNoOceanMS(49,300,15), tbout2DNoOceanMS(49,300,15),&
       tbout2DMS(49,300,15)
  integer :: ipol(15), ifreq(15), iobs(15), ifreq1
  integer :: ifreqG(15), ipolG(15) 
  integer,parameter :: nbinL=88
  real :: sfcRain(49,300),sfcRainStd(49,300)
  real :: rRate3D(nbinL,49,300),  rRate3Dstd(nbinL,49,300)
  real :: pwc3D(nbinL,49,300),  pwc3Dstd(nbinL,49,300)
  real :: zcKu3D(nbinL,49,300), d03D(nbinL,49,300), piaOut(49,300)

  real :: sfcRainMS(49,300),sfcRainStdMS(49,300),pia35m(49,300)
  real :: rRate3DMS(nbinL,49,300),  rRate3DstdMS(nbinL,49,300)
  real :: pwc3DMS(nbinL,49,300),  pwc3DstdMS(nbinL,49,300)
  real :: zcKu3DMS(nbinL,49,300), zcKa3DMS(nbinL,49,300), d03DMS(nbinL,49,300), &
          piaOutKuMS(49,300), piaOutKaMS(49,300)
  real, allocatable :: emissoutL(:,:,:), emis_out_NS(:,:,:), emis_out_MS(:,:,:) !sjm 8/10/15
  type (stormStructType) :: stormStruct
  real, allocatable ::  Yens(:,:), Xens(:,:), Yobs(:), Xup(:)
  real, allocatable :: geoloc(:), hFreqTbs(:,:), PRgeoloc(:), hFreqPRg(:,:,:)
  real, allocatable :: ndn(:), ndnp(:,:), xscalev(:), logdNwf(:), randemiss(:), dwind(:)
  real, allocatable  :: rhPCij(:,:), cldwPCij(:,:)
  type (radarRetType)    :: radarRet
  type (radarDataType)   :: radarData
  real :: emis_eofs(100,12) !SJM 7/9/2015
  real :: scLatPR(49,300),scLonPR(49,300),wfmap(49,300), fpmap(49,300,15), fpmapN(49,300,15)
  real :: w10(49,300), w10_out_NS(49,300), w10_out_MS(49,300), w10_min, w10_max, emis, relAz
  real :: w10_rms_NS(49,300), emis_rms_NS(49,300,13), w10_rms_MS(49,300), emis_rms_MS(49,300,13)
  real :: S1eiaPR(49,300), S2eiaPR(49,300)
  real :: sigmaZeroVarKu(49,300), sigmaZeroVarKa(49,300), sigmaZeroCov(49,300)
  integer :: iiad
  real :: dZms(49,300) !! MS addition Feb 10, 2017
  integer :: msFlag(49, 300) !!WSO addition Feb 11, 2017

end module local_RD_var

subroutine get_rain_type(raintype, nscans)
  use globalData
  integer :: raintype(49,300)
  integer :: nscans
  raintype=dPRData%rainType
  nscans=dPRData%n1c21
  print*, nscans
end subroutine get_rain_type

subroutine testr(it)
  integer :: it
  it=136
end subroutine testr

subroutine radarRetSub2(nmu2,  nmfreq2,   icL, tbRgrid,               &
      dprrain,ichunk,orbNumb,ialg,idir)
!  SFM  end    12/13/2013
  use local_RD_var
  use globalData
  use f90DataTypes
  use f90Types
  use cldclass
  use ran_mod
  use geophysEns
  use nbinMod
  !use tables2
  use weight
  Use BMCVparameters
  use emissMod
!begin  MG 10/29/15 add gEnv module
  use gEnv
!end    MG 10/29/15
!begin  WSO 9/14/13 incorporate missing flags
  use missingMod
!end    WSO 9/14/13
!begin  WSO 6/5/18 add limits to output variables
  use outputminmax
!end    WSO 6/5/18
  use LUT_def !SJM 7/9/2015
  implicit none
!  type (geoDataType) :: geoData
!  type (gridEnvDataType) :: gridEnvData
!  type (dPRDataType)     :: dPRData
!  type (dPRRetType)      :: dPRRet
!  type (gmi2GridType)    :: gmi2Grid
!  type (cgMIDataType)     :: gmiData
  integer :: nmu2, nmfreq2
  integer*4 :: ichunk
!  integer :: st_2adpr              ! file open status for 2adpr file
  real :: tbRgrid(14,49,9300), dprrain(49,300)
 ! integer :: tbRgridIn(9,49,9300)

  real                   :: meansfcRain,stddevSfcRain, tbout(14)
 
  !type (stormStructType) :: stormStruct
  type (retParamType)    :: retParam
  real :: xin363(363)
 
  integer :: ii, jj, iGMI, jGMI
  integer :: di(8), dj(8)  
  integer :: i, j, ig, jg, ntpw, nmemb1, itop, irand
  real    :: pia13m, rms1, rms2,  unSortedRR(200), corrcoef, sfcRain2, tpw_ij
  integer :: iy(200), kflag, it
  real    :: sysdNl, pia13mean
  integer :: iLandSea(5,5), i1, j1, igetlandsea, ic2, nobs, iit
  real :: a0(2,8), emiss(2), tb19v, tb19h, tb37v, tb37h, tb22v, tb85v, tb85h
  real :: meanTb85
  real :: stdTb85,kgain,kgain37, ymean(3)
!...Ensemble parameters

  integer :: ibatch
  real    :: stddev, srtpiaf
  real    :: FWHMx, FWHMy, tbconv(2), tbconvEns(2,100)
  integer :: dnx,dny, ik
  
  real :: cldw(nlayer), rh(nlayer), pia13s
  integer :: nx,ny, icount, imin
  real ::  xm1,xs,rmsmin, prob, probtot, rmstot
  real :: piaR(100), fPIA, z13m
  integer :: ntbpix, ntbpix2
  real :: emtbm(9)
  real :: zminsc
  real :: realOut(49)
 
  !real :: dZms(49,300) !! MS addition Feb 10, 2017
  !integer :: msFlag(49, 300) !!WSO addition Feb 11, 2017
!begin  WSO 2/8/17 new variables
  integer :: multiscatcalc_NS(49, 300), multiscatcalc_MS(49, 300)
  integer :: algotype_NS(49, 300), algotype_MS(49, 300)
  integer :: profclass_NS(49, 300), profclass_MS(49, 300)
  real :: subfootvariability_NS(49, 300), subfootvariability_MS(49, 300)
  real :: multiscatsurface_NS(49, 300), multiscatsurface_MS(49, 300)
  real :: skintempsigma_NS(49, 300), skintempsigma_MS(49, 300)
  real :: columnvaporsigma_NS(49, 300), columnvaporsigma_MS(49, 300)
  real :: columncloudliqsigma_NS(49, 300), columncloudliqsigma_MS(49, 300)
  real :: errorofdatafit_NS(49, 300), errorofdatafit_MS(49, 300)
  real :: initnw_NS(nbin, 49, 300), initnw_MS(nbin, 49, 300) 
  real :: princomp_NS(5, 49, 300), princomp_MS(5, 49, 300)
  real :: surfprecipbiasratio_NS(49, 300), surfprecipbiasratio_MS(49, 300)
!end    WSO 2/8/17 
  integer :: l, ipias(2)
  character*3 :: ifdpr, iftest
  character*90 :: outfile
  integer :: ink
  !integer :: ifreqG(15), ipolG(15) 
  DOUBLE PRECISION input(6)
  DOUBLE PRECISION output(2)
  real :: wfract(5,5), wfractm, wfractsd
  real                    emissv(n_chan)
  real                    emissh(n_chan)
  real                    emissv_std(n_chan)
  real                    emissh_std(n_chan)
  integer :: stype!SJM 7/9/2015
  real  :: vLand(18,18), vOcean(10,10)
  real  :: pMLand(18), pMOcean(10)
  real  :: mTbLand(9), mTbOcean(9)
  real  :: stTbLand(9), stTbOcean(9)
  double precision  :: xin(18), xpred(18), yout(9)

!begin  WSO 8/19/13 change Nw variable name (not dN) and add mu
  real :: cldwprof(88), cldiprof(88), log10NwMean(88), mu_mean_prof(88)
  integer *2 :: env_nodes(10, 49)
  real :: env_levs(10), ray_angle, pi
!end    WSO 8/19/13
  real :: lFract(49,300), sprobs, probs(100), rmsS(100)
  real :: covar, xf
  integer  :: orbNumb
  !begin SJM 7/25/14
  real :: s0Ku, s0Ka, s0stdKu, s0stdKa, s0corr, ds0Ku, ds0Ka
 
  !end SJM 7/25/2014
!begin WSO 8/8/13
  real :: gatelength
  real :: depthBB, depthML, depth
  real :: mu_mean(49, 300)
  real :: mu_meanMS(49, 300)

  real :: mlwc_frac(10, 49, 300)
  real :: mrate_frac(10, 49, 300)
  real :: mlwc_fracMS(10, 49, 300)
  real :: mrate_fracMS(10, 49, 300)
  real :: sfcRainLiqFrac(49, 300)
  real :: sfcRainLiqFracMS(49, 300)
  real :: tbMax1(15), tbMin1(15)

!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
!  SFM  begin  06/22/2014
  real :: wfractPix, windPert(100), windPertU(100), windPertV(100), qvPert(100), dnqv
!  SFM  end    06/22/2014
!  SFM  end    07/29/2014
!end   WSO 8/8/13
 

  integer :: actOb(49,300), iactOb
  integer :: jk, nf
  integer :: dig               ! SFM  04/16/2014  for M.Grecu
  real   :: cl(9,25), xin25(25),dtb(9)
  real   :: ebar, minl
  !real, allocatable :: geoloc(:), hFreqTbs(:,:), PRgeoloc(:), hFreqPRg(:,:,:)
  !integer*4 :: istart, iend
  integer :: iconv, ialg, icL
  real :: nubfc, stdpia35

  integer,parameter :: nscans=300, npixs=25, nlev=88, nchans=13
  integer :: nfreq, idir
  integer :: pType(nscans,npixs)
  real :: sfcTemp(nscans,npixs), cldw3d(nscans,npixs,nlev)
  integer :: clutFree(nscans,npixs)
  real :: pRate(nscans,npixs,nlev), swc3d(nscans,npixs,nlev), tbobsT(nscans,npixs,nchans)
  real :: z13(nscans,npixs,nlev),emiss2d(nscans,npixs,nchans)
  real :: nw3d(nscans,npixs,nlev), press3d(nscans,npixs,nlev), &
       airTemp3d(nscans,npixs,nlev),qv3d(nscans,npixs,nlev)
  integer :: binNodes(nscans,npixs,5)
  integer :: envNode(nscans,npixs,10)
  real    :: pRateOut(nscans,npixs,nlev), swcOut(nscans,npixs,nlev), nwOut(nscans,npixs,nlev)
  integer :: sfcBin(nscans,npixs)
  real    :: tbsim(nscans,npixs,nchans)
!begin  WSO 8/30/13 prescribe levels for environmental parameters
  data env_levs/18., 14., 10., 8., 6., 4., 2., 1., 0.5, 0./
  data pi/3.14159265/
  tbMax1(1:9)=(/ 285.96795654,284.71334839,286.23388672,&
       284.92977905,286.37451172,&
       287.90933228,286.43515015,284.8597412,284.51959229/)
  tbMin1(1:9)=(/175.59274292,95.61299896,&
       215.97111511,161.71897888,205.35340881,&
       143.96331787,143.96331787,87.26193237,87.26193237/)
  tbMin1(1:9)=(/167.43344116,   86.27742767, &
       200.01838684,  134.10232544, &
       233.02462769,&
       221.44563293,  159.45098877,  267.19430542,  242.64634705/)
  print*, orbNumb,ichunk
  !call openascii(orbNumb,ichunk)
  
  print*, dPRData%n1c21
  dPRRet%cldwcoeff=0
  ifreqG(1:13)=(/1,1,2,2,3,4,4,5,5,6,6,7,8/)
  ipolG(1:13)=(/1,2,1,2,1,1,2,1,2,1,2,1,1/)
! 1 is V
! 2 is H
  !call readclust()
!begin  MG Sept 15 2015
  !call openenkffilensl(orbNumb,ichunk)
  allocate(hFreqPRg(49,dPRData%n1c21,4))
  hFreqPRg=missing_r4
!end  MG
  dZms=missing_r4 !! MS addition Feb 10, 2017
  msFlag = missing_i2 !! WSO addition Feb 11, 2017
  print*, iEnd, iStart, ichunk
  if(iEnd-iStart>3) then
     allocate(geoloc(2*(iEnd+1-iStart)*81))
     allocate(hFreqTbs((iEnd+1-iStart)*81,4))
     allocate(prgeoloc(dPRData%n1c21*49*2))
     ic2=1
     
     do i=iStart,iEnd
        do j=70,150
           geoloc(ic2)=gmiData%S2lat(j,i)
           ic2=ic2+1
           geoloc(ic2)=gmiData%S2lon(j,i)
           ic2=ic2+1
        enddo
     enddo
     print*, minval(gmiData%S2lat(70:150,iStart:iEnd))
     print*, maxval(gmiData%S2lat(70:150,iStart:iEnd))
     print*, minval(gmiData%S2lon(70:150,iStart:iEnd))
     print*, maxval(gmiData%S2lon(70:150,iStart:iEnd))
     minl=minval(gmiData%S2lon(70:150,iStart:iEnd))
     do k=1,4
        ic2=1
        do i=iStart,iEnd
           do j=70,150
              hFreqTbs(ic2,k)=gmiData%gmiS2(k,j,i)
              ic2=ic2+1
           enddo
        enddo
     enddo
     
     
     ic2=1
     do i=1,dPRData%n1c21
        do j=1,49
           prgeoloc(ic2)=dPRData%xlat(j,i)
           ic2=ic2+1
           if(dPRData%xlon(j,i)+180>minl) then
              prgeoloc(ic2)=dPRData%xlon(j,i)
           else
              prgeoloc(ic2)=dPRData%xlon(j,i)!+360
           endif
           ic2=ic2+1
        enddo
     enddo
     print*, minval(dPRData%xlat)
     print*, maxval(dPRData%xlat)
     print*, minval(dPRData%xlon)
     print*, maxval(dPRData%xlon)


     do k=1,4
        call flannint(geoloc,  prgeoloc, hFreqTbs(:,k), hFreqPRg(:,:,k), &
             (iEnd+1-iStart)*81, 2, &
             dPRData%n1c21*49)
     enddo

     deallocate(geoloc)
     deallocate(prgeoloc)
     deallocate(hFreqTbs)
  end if

  di=(/0,  0, -1, 1, -2, 2, 0, 0/)
  dj=(/-1, 1,  0, 0, 0, 0, 2, -2/)

  a0(1,1:8)=(/0.5098,4.4664E-3,-6.0427E-6,&
       -2.5285E-3,-2.3725E-3,9.8163E-4,-2.2269E-3,-1.3193E-3/)
  a0(2,1:8)=(/0.3186,-1.5225E-3,1.7213E-3,&
       -3.7164E-4,6.5607E-3,8.1213E-4,-1.7678E-3,-1.7250E-3/)
  
  nmemb1=nmemb
  nobs=2
  allocate(ndn(nmemb1),xscalev(nmemb1),&
       logdNwf(9*nmemb1), randemiss(nmfreq2*nmemb1*2)) 
  allocate(ndnp(10,nmemb1)) 
  allocate(emissoutL(49,dPRData%n1c21,13))

  ntbpix=0
  ntbpix2=0
!  begin  SFM  07/29/2014; for M.Grecu,  eliminate NANs
  do k=0,nmemb1-1
     windPert(k+1) = normal2(0.,.15) !SJM 2/4/15
     windPertU(k+1) = normal2(0.,.15) !SJM 2/4/15
     windPertV(k+1) = normal2(0.,.15)
     qvPert(k+1)=ran1()
  end do
!  end    SFM  07/29/2014
  iactOb=0
  do k=0,nmemb1-1
     ndn(k+1)=normal2(0.,1.0)
     xscalev(k+1)=ran1()
     do i=1,9
        ndnp(i,k+1)=normal2(0.,1.0)
     enddo
  enddo
  do k=1,2*nmfreq2*nmemb1
     randemiss(k)=.5*ran1()
  enddo
  !set random emissivity/sigm0 PCs
  do k=0,nmemb1-1
     do i=1,12
       emis_eofs(k+1,i) = normal2(0.,1.)
     end do
  end do
 

  rRate3D=0.
  rRate3Dstd=0.
  pwc3D=0.
  pwc3Dstd=0.
!begin  WSO 9/15/13 set to flag instead of 0
  d03D = missing_r4
  zcKu3D = missing_r4
  piaOut = missing_r4
!end    WSO 9/15/13 
  sfcRain=0.
  sfcRainStd=0.
  dPRRet%z35mod0=0.
!begin WSO 04/07/13
  rRate3DMS=0.
  rRate3DstdMS=0.
  pwc3DMS=0.
  pwc3DstdMS=0.
!begin  WSO 9/15/13 set to flag instead of 0
  d03DMS = missing_r4
  zcKu3DMS = missing_r4
  zcKa3DMS = missing_r4
  piaOutKuMS = missing_r4
  piaOutKaMS = missing_r4
!end    WSO 9/15/13 
  sfcRainMS=0.
  sfcRainStdMS=0.
!begin  WSO 9/15/13 set to flag instead of 0
  w10 = missing_r4
  w10_rms_NS = missing_r4
  emis_rms_NS = missing_r4
  w10_rms_MS = missing_r4
  emis_rms_MS = missing_r4
!end    WSO 9/15/13
!end WSO 04/07/13

!end    WSO 8/8/13
!begin  WSO 8/19/13
  mu_mean = missing_r4
  mu_meanMS = missing_r4
!end    WSO 8/19/13
  !end    WSO 9/15/13
  
  print*, nmfreq2
  print*, nmemb1, ngates
  call allocateDPRProfRet(radarRet,nmfreq2,nmemb1,ngates, 9)   ! allocates memory for the 1D 
  print*, 'after allocate DPR'
  !...retrieval structures
  radarRet%rrate = 0.0  
  print*, 'radar assign'
  radarRet%tb=-99

  call allocateDPRProfData(radarData, ngates)                 ! allocates memory for the 
  print*, 'affter DPRPRofData'                               ! 1-D DPR observations
  call allocateStormStructData(stormStruct)                   ! allocates memory for the 5-node 
  print*, 'C allocation'                                                            ! storm structure
  call setRetParam(retParam)
  print*, 'setRetParam'
  call setrandclass(radarRet, nmu2)
  print*, 'setRand'
  print*, 'after allocate'
  dPRRet%sfc_wind(1:nmemb)=radarRet%sfc_wind(1:nmemb)
  
  dPRRet%sfcRainEns=0
  stormStruct%iSurf=ngates
  radarData%ngates=ngates
  radarData%dr=0.25
  dPRRet%tb=-99
  dPRRet%emtb=-99
  dPRRet%emis=-99
  dPRRet%n9=0
  call allocGeophys(6,61,9,nmemb1,nmfreq2*nmemb1*2)
  
  call setdNwIcJcL(sysdN,nmemb1)
  nx=30+nbin*7
!begin  MG 9/18/13 changed 110 to 130
  ny=130
!end    MG 9/18/
! SFM  begin  03/27/2014; execution protections

! SFM  end    03/27/2014

  print*, nx, nmemb1
  allocate(Xens(nx,nmemb1),Yens(ny,nmemb1),Yobs(ny), Xup(nx))
  
  allocate(rhPCij(nmemb1,nRhEofs), cldwPCij(nmemb1,nCldwEofs))
!00000000000000000000000000000000000000000000000000000000000000000000000000000
  lFract=-99.

  tb19v=0
  tb19h=0
  tb22v=0
  tb37v=0
  tb37h=0
  tb85v=0
  tb85h=0
  ic2=0
!begin  WSO 9/16/13
  emissoutL=missing_r4
!end    WSO 9/16/13
!begin WSO 04/07/13
   ifdpr(1:1) = 'N'
   iftest(1:1) = 'N'


  sfcRain=0
  !begin SJM 12/9/2014
  w10=dPRData%envSfcWind
  w10_out_NS=dPRData%envSfcWind
  w10_out_MS=dPRData%envSfcWind
  
  sigmaZeroVarKu = 0.
  sigmaZeroVarKa = 0.
  sigmaZeroCov = 0.
  S1eiaPR = 52.7
  S2eiaPR = 49.1
  print*,'max_gmiData%tpw3=', maxval(gmiData%tpw3)
!   write(outfile,'(A64,I6.6,A4)') '/PANFS/user/home/smunchak/data/combAlg-tests/ensdata/InitialEns.',orbNumb,'.bin'
!   if(ichunk .eq. 0) then
!     open(31,file=outfile,form='unformatted') 
!   else
!     open(31,file=outfile,access='append',form='unformatted') 
!   endif
!   write(outfile,'(A64,I6.6,A4)') '/PANFS/user/home/smunchak/data/combAlg-tests/ensdata/NSFiltered.',orbNumb,'.bin'
!   if(ichunk .eq. 0) then
!     open(32,file=outfile,form='unformatted') 
!   else
!     open(32,file=outfile,access='append',form='unformatted') 
!   endif
!   write(outfile,'(A64,I6.6,A4)') '/PANFS/user/home/smunchak/data/combAlg-tests/ensdata/MSFiltered.',orbNumb,'.bin'
!   if(ichunk .eq. 0) then
!     open(33,file=outfile,form='unformatted') 
!   else
!     open(33,file=outfile,access='append',form='unformatted') 
!   endif
!   if(ichunk .eq. 0) then
!     open(34,file='deconvTb.bin',form='unformatted')
!   else
!     open(34,file='deconvTb.bin',access='append',form='unformatted') 
!   endif
  !end SJM 12/9/2014
  w10=dPRData%envSfcWind
  actOb=0

  dPRRet%convtb=-99
  scLonPR=-99
  scLatPR=-99
  tbMean=-99
  tbMax=-99
  tbMin=-99
  !call startprofs()
  rrate3D=0!-99
  tbNoOcean=-99
  pia35m=0.
  print*, gmi2Grid%xmin, gmi2Grid%ymin
  !return 
  tbRgrid=-9999
  do j=1,dPRData%n1c21
  !do j=108,108
     do i=1,49
     !do i=25,25      
        eLon=dPRData%xlon(i,j)
        eLat=dPRData%xlat(i,j)
        call getwfraction(eLat,&
             eLon,wfmap(i,j))
        wfmap(i,j)=wfmap(i,j)/100.
        !print*, 'here in the loop'      
        !begin  WSO 8/21/14 initialize ioquality flag to bad
        dPRData%ioqualityflagku(i, j) = dPRData%ioqualityflagku(i,j) + 900000
        if(i > 12 .and. i < 38) then
           dPRData%ioqualityflagdpr(i, j) = dPRData%ioqualityflagdpr(i, j) + 900000
        endif
        !end    WSO 8/21/14
        !  SFM  begin  04/16/2014; for M.Grecu, revision of nodes processing
        if(dPRData%xlon(i,j)>-998) then  !4/15/14 MG begin
           iLandSea=-2
           iGMI=-99  !4/22/14 MG
           jGMI=-99  !4/22/14 MG
           !begin  WSO 3/16/17 initialize ig, jg
           ig = -99
           jg = -99
           !end    WSO 3/16/17
           if(gmi2Grid%xmin>-998 .and. gmi2Grid%ymin>-998) then
              ig=(dPRData%xlon(i,j)-gmi2Grid%xmin)/gmi2Grid%dx+1
              jg=(dPRData%xlat(i,j)-gmi2Grid%ymin)/gmi2Grid%dx+1
              if(ig>0 .and. ig<gmi2Grid%nx .and. jg>0 .and. jg<gmi2Grid%ny) then
                 iGMI=gmi2Grid%ig(ig,jg)
                 jGMI=gmi2Grid%jg(ig,jg)
                 if(gmi2Grid%actOb(ig,jg)==1) then
                    actOb(i,j)=1
                    iactob=iactob+1
                 endif
                 
              else
                 if(jg>0 .and. jg<gmi2Grid%ny) then
                    dig=int(360/gmi2Grid%dx)
                    if(ig+dig>0 .and.  ig+dig<gmi2Grid%nx) then
                       ig=ig+dig
                       iGMI=gmi2Grid%ig(ig,jg)
                       jGMI=gmi2Grid%jg(ig,jg)
                    else
                       if(ig+2*dig>0 .and.  ig+2*dig<gmi2Grid%nx) then
                          ig=ig+2*dig
                          iGMI=gmi2Grid%ig(ig,jg)
                          jGMI=gmi2Grid%jg(ig,jg)
                       else
                          iGMI=-99
                          jGMI=-99
                       endif
                    endif
                 endif
              endif
           else
              iGMI=-99
              jGMI=-99
           endif   ! 4/15/14 MG End
           !print*, iGMI, jGMI, i, j
           !  SFM  end    04/16/2104
           k=1
           if(iGMI<0) then
              do while(iGMI<0 .and. k<8)
                 if( ig+di(k)>0 .and. ig+di(k)<=gmi2Grid%nx .and.                 &
                      jg+dj(k)>0 .and. jg+dj(k)<=gmi2Grid%ny) then
                    iGMI=gmi2Grid%ig(ig+di(k),jg+dj(k))
                    jGMI=gmi2Grid%jg(ig+di(k),jg+dj(k))
                    
                 endif
                 k=k+1
              end do
           endif
           dPRData%ig(i,j)=iGMI
           dPRData%jg(i,j)=jGMI   
           !print*, iGMI,jGMI, i, j, gmiData%tpw3(iGMI,jGMI), gmi2Grid%xmin, icL  
           if(jGMI>-99) then
              scLonPR(i,j)=gmiData%SCLon3(jGMI)
              scLatPR(i,j)=gmiData%SCLat3(jGMI)
              S1eiaPR(i,j)=gmidata%S1eia3(iGMI,jGMI)
              S2eiaPR(i,j)=gmidata%S2eia3(iGMI,jGMI)
              !print*,i,j,S1eiaPR(i,j), S2eiaPR(i,j)
           endif
           
           if(iGMI>0 .and. jGMI>0 .and. gmi2Grid%xmin>-998 ) then !4/22/14 MG
              if(gmiData%tpw3(iGMI,jGMI)>0) then
                 if(gmiData%sfc_wind3(iGMI,jGMI)>0) then
                    radarRet%sfc_wind=gmiData%sfc_wind3(iGMI,jGMI)
                    tbRgrid(14,i,j+icL)=gmiData%tpw3(iGMI,jGMI) 
                 else
                    radarRet%sfc_wind=dPRData%envSfcWind(i,j)
                    tbRgrid(14,i,j+icL)=-99
                 endif
                 !print*, ifdpr(1_, iftest(1), icL
                 if(ifdpr(1:1).ne.'Y'.and.iftest(1:1).ne.'Y') then
                    tbRgrid(1:9,i,j+icL)=gmiData%gmiS13(1:9,iGMI,jGMI)
                    
                    !begin  WSO 8/21/2014 for at least some good data, re-set quality to nominal
                    
                    if(maxval(tbRgrid(1:9, i, j + icL)) > 0) then
                       dPRData%ioqualityflagku(i, j) = dPRData%ioqualityflagku(i, j) - 900000
                       if(i > 12 .and. i < 38) then
                          dPRData%ioqualityflagdpr(i, j) = dPRData%ioqualityflagdpr(i, j) - 900000
                       endif
                    endif
                    
                 endif
              else
                 if(ifdpr(1:1).ne.'Y'.and.iftest(1:1).ne.'Y') then
                    tbRgrid(1:9,i,j+icL)=gmiData%gmiS13(1:9,iGMI,jGMI)
                    
                    if(maxval(tbRgrid(1:9, i, j + icL)) > 0) then
                       dPRData%ioqualityflagku(i, j) = &
                            dPRData%ioqualityflagku(i, j) - 900000
                       if(i > 12 .and. i < 38) then
                          dPRData%ioqualityflagdpr(i, j) = &
                               dPRData%ioqualityflagdpr(i, j) - 900000
                       endif
                    endif
                    
                    
                 endif
                 radarRet%sfc_wind=dPRData%envSfcWind(i,j)
              endif
           endif
           
           !begin  WSO 9/5/13 rename SRT PIA
           if(dPRData%rainType(i,j)>=100 ) then
              !end    WSO 9/5/13 
              !print*, dPRData%node(:,i,j)
              dPRData%node(5,i,j)=dPRData%node(5,i,j)
              
              !end    WSO 9/5/13 
              if (dPRData%rainType(i,j)/100==2) then
                 dPRData%node(2,i,j)=dPRData%node(2,i,j)-1
              endif
              if (dPRData%rainType(i,j)/100==1) then
                 dPRData%node(2,i,j)=dPRData%node(2,i,j)-1
                 dPRData%node(4,i,j)=dPRData%node(4,i,j)+1
              endif
              stormStruct%nodes  = dPRData%node(:,i,j)
              radarData%z13obs   = dPRData%zku1c21(:,i,j)
              radarData%z35obs   = dPRData%zka1c21(:,i,j)
              !radarData%z35obs   = -99
              !begin  WSO 9/5/13 rename SRT and DSRT PIA's and reliability factor
              radarData%pia13srt = dPRData%srtPIAku(i,j)
              radarData%relpia13srt = dPRData%srtrelPIAku(i,j)
              radarData%pia35srt = -99 
              radarData%pia35srt = dPRData%dsrtPIAka(i,j)
              !end    WSO 9/5/13
              !begin SJM 7/10/14 add sigma-zero
              radarData%sigmaZeroKu = dPRData%sigmaZeroKu(i,j)
              radarData%sigmaZeroKa = dPRData%sigmaZeroKa(i,j)
              !print*, radarData%sigmaZeroKu, radarData%sigmaZeroKa
              !end SJM 7/10/14
              radarData%hfreez   = dPRData%freezH(i,j) /1000. 
              stormStruct%iSurf = dPRData%binRealSurface(i,j)+1
              
              do k=1,nbin
                 if(radarData%z13obs(k)<-99) radarData%z13obs(k)=-99
              enddo
              
              !if(iGMI>0 .and. jGMI>0 ) then
              !   if(gmiData%tpw3(iGMI,jGMI)<0) then
              !      call landEmiss(emissout(i,j,1:9),dPRData,i,j,&
              !           dPRData%envSfcWind(i,j),dPRData%envSfcTemp(i,j))
              !   endif
              !endif
              
              stormStruct%rainType=dPRData%rainType(i,j)
              stormStruct%rainType=stormStruct%rainType/100
              
              itop=1
              radarRet%sfc_wind(1)=dPRData%envSfcWind(i,j)
              radarRet%sfc_windU(1)=dPRData%envSfcWindU(i,j)
              radarRet%sfc_windV(1)=dPRData%envSfcWindV(i,j)
              w10(i,j)=radarRet%sfc_wind(1)
              !print*, i,j,radarRet%sfc_wind(1),radarRet%sfc_windU(1),radarRet%sfc_windV(1)
              !w10(i,j)=0.5*(radarRet%sfc_wind(1)+dPRData%envSfcWind(i,j))
              
              if(dPRData%rainType(i,j)/100>=1) then

                 do ibatch=1,1
                    eLon=dPRData%xlon(i,j)
                    eLat=dPRData%xlat(i,j)
                    call getwfraction(eLat,&
                         eLon,wfractPix)
                    call interpoldNw(i,j, logdNwf)
                    
                    if(stormStruct%rainType == 1) then
                       do k=0,nmemb1-1 
                          radarRet%logdNw(k*9+1:(k+1)*9)=sysdn-0.1+ & !Sept 17, 2015 MG
                               logdNwf(k*9+1:(k+1)*9)
                       enddo
                    else
                       if(dPRData%rainType(i,j)/100==2) then
                          do k=0,nmemb1-1
                             !  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
                             !  SFM  begin  06/22/2014; for M.Grecu  (unknown justification)
                             radarRet%logdNw(k*9+1:(k+1)*9)=sysdn+0.0+            &
                                  logdNwf(k*9+1:(k+1)*9)
                             !  SFM  end    06/22/2014
                             !  SFM  end    07/29/2014
                          enddo
                       else
                          do k=0,nmemb1-1
                             radarRet%logdNw(k*9+1:(k+1)*9)=sysdn+0.0+            &
                                  logdNwf(k*9+1:(k+1)*9)
                          enddo
                       endif
                    endif
                    if(stormStruct%rainType .eq. 1 .and. dPRData%BBbin(i,j)<=0)   &
                         then
                       stormStruct%rainType = 3
                    endif
                    iRad=i
                    jRad=j
                    if(wfractPix<90 .and.  stormStruct%rainType ==1) then
                       do k=0,nmemb1-1 
                          radarRet%logdNw(k*9+1:(k+1)*9)=&
                               radarRet%logdNw(k*9+1:(k+1)*9)-0.1
                       enddo
                    endif
                    do k=0,nmemb1-1
                       xscalev(k+1)=1
                    enddo
                    
                    
                    reliabFlag=dPRData%NSRelibFlag(i,j)
                    !print*, 'before ens'
                    call  ensRadRetStCvKu(radarData,stormStruct,                  &
                         retParam, nmu2,radarRet, itop, rms1, rms2, sysdN, iit, &
                         xscalev, randemiss, dPRData%localZenithAngle(i,j), &
                         wfractPix, ichunk, i, j, dZms(i,j), msFlag(i, j)) 
                    !print*, 'but not here'
                    !begin WSO 2/11/17 assign dZms to output variable
                    multiscatsurface_MS(i, j) = dZms(i, j)
                    multiscatcalc_MS(i, j) = msFlag(i, j)
                    !end   WSO 2/11/17
                    
                    !begin SJM 10/16/15
                    !generate simulated sigma_zero over water
                    !print*, wfractPix
                    

                    dPRRet%n9(:,i,j)=n9+1
                    
                    dPRRet%cldwcoeff(i,j,:,:)=cldwcoeff(1:10,1:nmemb)
                    
                   
                    if(radarRet%tb(2)>0) ntbpix=ntbpix+1
                    
                    do k=0,nmemb1-1
                       iy(k+1)=k+1
                    enddo
                    do k=0,nmemb1-1
                       dPRRet%tb(i,j,2,:,k+1+(ibatch-1)*nmemb1)=                  &
                            radarRet%tb((iy(k+1)-1)*2*radarRet%nmfreq+1:          &
                            2*(iy(k+1)-1)*radarRet%nmfreq+radarRet%nmfreq)
                       dPRRet%tb(i,j,1,:,k+1+(ibatch-1)*nmemb1)=                  &
                            radarRet%tb(((iy(k+1)-1)*2+1)*radarRet%nmfreq+1:      &
                            2*((iy(k+1)-1)+1)*radarRet%nmfreq)
                       dPRRet%emtb(i,j,2,:,k+1+(ibatch-1)*nmemb1)=           &
                            radarRet%emtb((iy(k+1)-1)*2*radarRet%nmfreq+1:    &
                            2*(iy(k+1)-1)*radarRet%nmfreq+radarRet%nmfreq)
                       dPRRet%emtb(i,j,1,:,k+1+(ibatch-1)*nmemb1)=             &
                            radarRet%emtb(((iy(k+1)-1)*2+1)*radarRet%nmfreq+1:   &
                            2*((iy(k+1)-1)+1)*radarRet%nmfreq)
                       !begin SJM 10/16/15
                       dPRRet%emis(i,j,1,:,k+1+(ibatch-1)*nmemb1)=             &
                            radarRet%emis(((iy(k+1)-1)*2+1)*radarRet%nmfreq+1:   &
                            2*((iy(k+1)-1)+1)*radarRet%nmfreq)
                       dPRRet%emis(i,j,2,:,k+1+(ibatch-1)*nmemb1)=           &
                            radarRet%emis((iy(k+1)-1)*2*radarRet%nmfreq+1:    &
                            2*(iy(k+1)-1)*radarRet%nmfreq+radarRet%nmfreq)
                       !end SJM 10/16/15
                       dPRRet%log10dNw (k+1+(ibatch-1)*nmemb1,:,i,j)=-99
                       dPRRet%d0 (k+1+(ibatch-1)*nmemb1,:,i,j)=-99
                    enddo
                    !print*, dPRRet%emis(i,j,1,5,:)
                    !                 print*, dPRRet%tb(i,j,1,:,1)
                    !  ifreqG(1:9)=(/1,1,2,2,3,4,4,5,5/)
                    !  ipolG(1:9)=(/1,2,1,2,1,1,2,1,2/)
                    do ik=1,9
                       tbMean(i,j,ik)=&
                            sum(dPRRet%tb(i,j,ipolG(ik),ifreqG(ik),1:nmemb1))/&
                            nmemb1
                       tbMin(i,j,ik)=&
                            minval(dPRRet%tb(i,j,ipolG(ik),ifreqG(ik),1:nmemb1))
                       tbMax(i,j,ik)=&
                            maxval(dPRRet%tb(i,j,ipolG(ik),ifreqG(ik),1:nmemb1))
                       tbNoOcean(i,j,ik)=tbMean(i,j,ik)
                    enddo
                    
                    do ik=1,9
                       do jk=1,9
                          covTb(i,j,ik,jk)=&
                               covar(dPRRet%tb(i,j,ipolG(ik),ifreqG(ik),1:nmemb1),&
                               dPRRet%tb(i,j,ipolG(jk),ifreqG(jk),1:nmemb1),nmemb1)
                       enddo
                       
                    enddo
121                 format(81(F8.2,1x))
                   
                   
                    
                    do k=0,nmemb1-1
                       dPRRet%log10dNw (k+1+(ibatch-1)*nmemb1,1:ngates,i,j)=      &
                            radarRet%log10dNw((iy(k+1)-1)*ngates+1:(iy(k+1))*ngates)
                       dPRRet%d0 (k+1+(ibatch-1)*nmemb1,1:ngates,i,j)=            &
                            radarRet%d0((iy(k+1)-1)*ngates+1:(iy(k+1))*ngates)
                       dPRRet%rrate (k+1+(ibatch-1)*nmemb1,1:ngates,i,j)=         &
                            radarRet%rrate((iy(k+1)-1)*ngates+1:(iy(k+1))*ngates)
                       dPRRet%pwc (k+1+(ibatch-1)*nmemb1,1:ngates,i,j)=           &
                            (radarRet%pwc((iy(k+1)-1)*ngates+1:(iy(k+1))*ngates))
                       dPRRet%z13c(k+1+(ibatch-1)*nmemb1,1:ngates,i,j)=           &
                            radarRet%z13c((iy(k+1)-1)*ngates+1:(iy(k+1))*ngates)
                       dPRRet%z35mod0(k+1+(ibatch-1)*nmemb1,1:ngates,i,j)=        &
                            radarRet%z35mod0((iy(k+1)-1)*ngates+1:(iy(k+1))*ngates)
                       dPRRet%z35(k+1+(ibatch-1)*nmemb1,1:ngates,i,j)=            &
                            radarRet%z35((iy(k+1)-1)*ngates+1:(iy(k+1))*ngates)
                    enddo
                    
                    meansfcRain=0
                    !print*, dPRData%node(5,i,j)
                    do k=0,nmemb1-1
                       dPRRet%sfcRainEns(i,j,k+1+(ibatch-1)*nmemb1)=              &
                            radarRet%rrate((iy(k+1)-1)*ngates+dPRData%node(5,i,j))
                       dPRRet%sfcd0Ens(i,j,k+1+(ibatch-1)*nmemb1)=                &
                            radarRet%d0((iy(k+1)-1)*ngates+1+dPRData%node(5,i,j))
                       meansfcRain=meansfcRain+ dPRRet%sfcRainEns(i,j,k+1)
                       dPRRet%sfcNwEns(i,j,k+1+(ibatch-1)*nmemb1)=                &
                            radarRet%log10dNw((iy(k+1)-1)*ngates+1                &
                            +dPRData%node(5,i,j))
                       dPRRet%sfcWindEns(i,j,k+1+(ibatch-1)*nmemb1)=radarRet%sfc_wind(k+1) !SJM 12/4/2014
                       dPRRet%pia13mod(i,j,k+1+(ibatch-1)*nmemb1)=                &
                            radarRet%pia13(iy(k+1))
                       dPRRet%pia35mod(i,j,k+1+(ibatch-1)*nmemb1)=&
                            radarRet%pia35(iy(k+1))
                       dPRRet%simSigmaZeroKu(i,j,k+1+(ibatch-1)*nmemb1)=                & !SJM 2/4/2015
                            radarRet%simSigmaZeroKu(iy(k+1))
                       dPRRet%simSigmaZeroKa(i,j,k+1+(ibatch-1)*nmemb1)=                & !SJM 2/4/2015
                            radarRet%simSigmaZeroKa(iy(k+1))
                    enddo
                    pia13m=sum(dPRRet%pia13mod(i,j,1:nmemb1))/nmemb1
                    pia35m(i,j)=sum(dPRRet%pia35mod(i,j,1:nmemb1))/nmemb1
                    pia13s=(sum((dPRRet%pia13mod(i,j,1:nmemb1)-pia13m)**2)/&
                         nmemb1)**.5
                    sfcRain(i,j)=meansfcRain/nmemb1
                    
                    sfcRainStd(i,j)=sqrt(sum((dPRRet%sfcRainEns(i,j,1:nmemb1)- &
                         sfcRain(i,j))**2)/(nmemb1-1))
                    piaOut(i,j)=pia13m
                 enddo
              endif
           else
              do i1=1,5
                 do j1=1,5
                    eLon=dPRData%xlon(i,j)
                    eLat=dPRData%xlat(i,j)
                    call getwfraction(eLat+(i1-2)*0.15,&
                         eLon+(j1-2)*0.15,wfract(i1,j1))
                 enddo
              enddo
              
              do ik=1,9
                 tbNoOcean(i,j,ik)=tbRgrid(ik,i,j+icL)
              enddo
              do ik=1,9
                 if(wfract(3,3)<5) then
                    tbNoOcean(i,j,ik)=tbRgrid(ik,i,j+icL)
                 endif
              enddo
           endif
        endif
        !if(maxval(emissout(i,j,1:13)) .gt. 2.) print*, i,j, emissout(i,j,1:13)
        !if(maxval(dPRRet%emis(i,j,:,:,:)) .gt. 2.) print*, dPRRet%emis(i,j,:,:,:)
     enddo
  enddo
end subroutine radarRetSub2

subroutine radarRetSub3(nmu2,  nmfreq2,   icL, tbRgrid,               &
      dprrain,ichunk,orbNumb,ialg,idir)
!  SFM  end    12/13/2013
  use local_RD_var
  use globalData
  use f90DataTypes
  use f90Types
  use cldclass
  use ran_mod
  use geophysEns
  use nbinMod
  !use tables2
  use weight
  Use BMCVparameters
  use emissMod
!begin  MG 10/29/15 add gEnv module
  use gEnv
!end    MG 10/29/15
!begin  WSO 9/14/13 incorporate missing flags
  use missingMod
!end    WSO 9/14/13
!begin  WSO 6/5/18 add limits to output variables
  use outputminmax
!end    WSO 6/5/18
  use LUT_def !SJM 7/9/2015
  implicit none
!  type (geoDataType) :: geoData
!  type (gridEnvDataType) :: gridEnvData
!  type (dPRDataType)     :: dPRData
!  type (dPRRetType)      :: dPRRet
!  type (gmi2GridType)    :: gmi2Grid
!  type (cgMIDataType)     :: gmiData
  integer :: nmu2, nmfreq2
  integer*4 :: ichunk
!  integer :: st_2adpr              ! file open status for 2adpr file
  real :: tbRgrid(14,49,9300), dprrain(49,300)
 ! integer :: tbRgridIn(9,49,9300)

  real                   :: meansfcRain,stddevSfcRain, tbout(14)
  !type (radarRetType)    :: radarRet
  !type (radarDataType)   :: radarData
  type (retParamType)    :: retParam
  real :: xin363(363)
 
  integer :: ii, jj, iGMI, jGMI
  integer :: di(8), dj(8)  
  integer :: i, j, ig, jg, ntpw, nmemb1, itop, irand
  real    :: pia13m, rms1, rms2,  unSortedRR(200), corrcoef, sfcRain2, tpw_ij
  integer :: iy(200), kflag, it
  real    :: sysdNl, pia13mean
  integer :: iLandSea(5,5), i1, j1, igetlandsea, ic2, nobs, iit
  real :: a0(2,8), emiss(2), tb19v, tb19h, tb37v, tb37h, tb22v, tb85v, tb85h
  real :: meanTb85
  real :: stdTb85,kgain,kgain37, ymean(3)
!...Ensemble parameters
  !real, allocatable ::  Yens(:,:), Xens(:,:), Yobs(:), Xup(:)
  integer :: ibatch
  real    :: stddev, srtpiaf
  real    :: FWHMx, FWHMy, tbconv(2), tbconvEns(2,100)
  integer :: dnx,dny, ik
  !real, allocatable :: ndn(:), ndnp(:,:), xscalev(:), logdNwf(:), randemiss(:), dwind(:)
  !real, allocatable  :: rhPCij(:,:), cldwPCij(:,:)
  real :: cldw(nlayer), rh(nlayer), pia13s
  integer :: nx,ny, icount, imin
  real ::  xm1,xs,rmsmin, prob, probtot, rmstot
  real :: piaR(100), fPIA, z13m
  integer :: ntbpix, ntbpix2
  real :: emtbm(9)
  real :: zminsc
  real :: realOut(49)
  !real :: w10(49,300), w10_out_NS(49,300), w10_out_MS(49,300), w10_min, w10_max, emis, relAz
  !real :: w10_rms_NS(49,300), emis_rms_NS(49,300,13), w10_rms_MS(49,300), emis_rms_MS(49,300,13)
  !real :: dZms(49,300) !! MS addition Feb 10, 2017
  !integer :: msFlag(49, 300) !!WSO addition Feb 11, 2017
!begin  WSO 2/8/17 new variables
  integer :: multiscatcalc_NS(49, 300), multiscatcalc_MS(49, 300)
  integer :: algotype_NS(49, 300), algotype_MS(49, 300)
  integer :: profclass_NS(49, 300), profclass_MS(49, 300)
  real :: subfootvariability_NS(49, 300), subfootvariability_MS(49, 300)
  real :: multiscatsurface_NS(49, 300), multiscatsurface_MS(49, 300)
  real :: skintempsigma_NS(49, 300), skintempsigma_MS(49, 300)
  real :: columnvaporsigma_NS(49, 300), columnvaporsigma_MS(49, 300)
  real :: columncloudliqsigma_NS(49, 300), columncloudliqsigma_MS(49, 300)
  real :: errorofdatafit_NS(49, 300), errorofdatafit_MS(49, 300)
  real :: initnw_NS(nbin, 49, 300), initnw_MS(nbin, 49, 300) 
  real :: princomp_NS(5, 49, 300), princomp_MS(5, 49, 300)
  real :: surfprecipbiasratio_NS(49, 300), surfprecipbiasratio_MS(49, 300)
!end    WSO 2/8/17 
  integer :: l, ipias(2)
  character*3 :: ifdpr, iftest
  character*90 :: outfile
  integer :: ink
  DOUBLE PRECISION input(6)
  DOUBLE PRECISION output(2)
  real :: wfract(5,5), wfractm, wfractsd
  real                    emissv(n_chan)
  real                    emissh(n_chan)
  real                    emissv_std(n_chan)
  real                    emissh_std(n_chan)
  integer :: stype!SJM 7/9/2015
  real  :: vLand(18,18), vOcean(10,10)
  real  :: pMLand(18), pMOcean(10)
  real  :: mTbLand(9), mTbOcean(9)
  real  :: stTbLand(9), stTbOcean(9)
  double precision  :: xin(18), xpred(18), yout(9)
  !real, allocatable :: emissoutL(:,:,:), emis_out_NS(:,:,:), emis_out_MS(:,:,:) !sjm 8/10/15
!begin  WSO 8/19/13 change Nw variable name (not dN) and add mu
  real :: cldwprof(88), cldiprof(88), log10NwMean(88), mu_mean_prof(88)
  integer *2 :: env_nodes(10, 49)
  real :: env_levs(10), ray_angle, pi
!end    WSO 8/19/13
  real :: lFract(49,300), sprobs, probs(100), rmsS(100)
  real :: covar, xf
  integer  :: orbNumb
  !begin SJM 7/25/14
  real :: s0Ku, s0Ka, s0stdKu, s0stdKa, s0corr, ds0Ku, ds0Ka
  !real :: sigmaZeroVarKu(49,300), sigmaZeroVarKa(49,300), sigmaZeroCov(49,300)
  !end SJM 7/25/2014
!begin WSO 8/8/13
  real :: gatelength
  real :: depthBB, depthML, depth
  real :: mu_mean(49, 300)
  real :: mu_meanMS(49, 300)
!  real :: scLatPR(49,300),scLonPR(49,300),wfmap(49,300), fpmap(49,300,15), fpmapN(49,300,15)
!  real :: S1eiaPR(49,300), S2eiaPR(49,300)
  real :: mlwc_frac(10, 49, 300)
  real :: mrate_frac(10, 49, 300)
  real :: mlwc_fracMS(10, 49, 300)
  real :: mrate_fracMS(10, 49, 300)
  real :: sfcRainLiqFrac(49, 300)
  real :: sfcRainLiqFracMS(49, 300)
  real :: tbMax1(15), tbMin1(15)

!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
!  SFM  begin  06/22/2014
  real :: wfractPix, windPert(100), windPertU(100), windPertV(100), qvPert(100), dnqv
!  SFM  end    06/22/2014
!  SFM  end    07/29/2014
!end   WSO 8/8/13
  

  integer :: actOb(49,300), iactOb
  integer :: jk, nf
  integer :: dig               ! SFM  04/16/2014  for M.Grecu
  real   :: cl(9,25), xin25(25),dtb(9)
  real   :: ebar, minl

  !integer*4 :: istart, iend
  integer :: iconv, ialg, icL
  real :: nubfc, stdpia35

  integer,parameter :: nscans=300, npixs=25, nlev=88, nchans=13
  integer :: nfreq, idir
  integer :: pType(nscans,npixs)
  real :: sfcTemp(nscans,npixs), cldw3d(nscans,npixs,nlev)
  integer :: clutFree(nscans,npixs)
  real :: pRate(nscans,npixs,nlev), swc3d(nscans,npixs,nlev), tbobsT(nscans,npixs,nchans)
  real :: z13(nscans,npixs,nlev),emiss2d(nscans,npixs,nchans)
  real :: nw3d(nscans,npixs,nlev), press3d(nscans,npixs,nlev), &
       airTemp3d(nscans,npixs,nlev),qv3d(nscans,npixs,nlev)
  integer :: binNodes(nscans,npixs,5)
  integer :: envNode(nscans,npixs,10)
  real    :: pRateOut(nscans,npixs,nlev), swcOut(nscans,npixs,nlev), nwOut(nscans,npixs,nlev)
  integer :: sfcBin(nscans,npixs)
  real    :: tbsim(nscans,npixs,nchans)


end subroutine radarRetSub3



subroutine dealloc_struct(i)
  use local_RD_var
  use geophysEns
  print*,i
  IF (ALLOCATED(emissoutL)) deallocate(emissoutL)
  if(allocated(hFreqPRg)) deallocate(hFreqPRg)

  IF (ALLOCATED(Yobs)) deallocate(Yobs)
  IF (ALLOCATED(Xup)) deallocate(Xup)
  IF (ALLOCATED(randemiss)) deallocate(randemiss)
  IF (ALLOCATED(Xens)) deallocate(Xens)
  IF (ALLOCATED(Yens)) deallocate(Yens)

  call deallocGeophys()
  call deallocateStormStructData(stormStruct)
  call deallocateDPRProfRet(radarRet)
  call deallocateDPRProfData(radarData)
  
  IF (ALLOCATED(ndn)) deallocate(ndn) 
  IF (ALLOCATED(xscalev)) deallocate(xscalev) 
  IF (ALLOCATED(logdnwf)) deallocate(logdnwf) 
  IF (ALLOCATED(ndnp)) deallocate(ndnp) 
  IF (ALLOCATED(rhPCij)) deallocate(rhPCij)
  IF (ALLOCATED(cldwPCij)) deallocate(cldwPCij)

end subroutine dealloc_struct

subroutine dealloc_chunk(i)
  use globalData
  integer :: i
  print*, i
  call deallocateDPRRetSpace(dPRRet)
  call deallocateHRescGMI(gmiData,gmi2Grid)
end subroutine dealloc_chunk
