! test_keras.f90

! TO RUN
! ./test_keras $NF_PATH/ExampleModels/simple_model.txt

! this file is used in $NF_PATH/KerasWeightsProcessing/examples/test_network.py
! load a specified network from cmd line arg
! pass simple input through it
! print result


subroutine init_keras_chase()
  use keras_def
  call net_ch % load("FKB/chase_model.txt")

  input_ch = [10.0,1.0,-5.0]

  ! run test input through network
  result1 = net_ch % output(input_ch)
  !print *, result1

end subroutine init_keras_chase

subroutine init_keras2()
  use keras_def
  call net % load("FKB/nw_dm_Kuo_model.txt")
  print*, 'here'
  
  input2 = [10.0,1.0,-5.0]

  ! run test input through network
  result2 = net % output(input2)
  print *, result2

end subroutine init_keras2


subroutine call_keras(input_data,output_data,n)
  use keras_def
  integer :: n
  real :: input_data(n,3)
  real, intent(out) :: output_data(n,3)
  integer :: i
  real :: xm(3),xs(3),ym(3),ys(3), dwr, zku,zka
  real :: yhat(3)
  real :: IWC,Nw,Dm,pi
  xm=(/11.61196423,   1.62957636, -19.86298585/)
  xs=(/13.30597739,  2.15601794, 11.72233551/)
  ym=(/8.00513416, -0.26086958,  0.11921838/)
  ys=(/0.78638963, 0.21327072, 0.291193/)
  pi=3.141592653589793
  do i=1,n
     zku=input_data(i,1)
     zka=input_data(i,2)
     if (zku<0) zku=0
     if (zka<0) zka=0
     dwr=zku-zka
     if(dwr<-1) dwr=-1
     if(dwr>11) dwr=11
     input(1)=(zku-xm(1))/xs(1)
     input(2)=(dwr-xm(2))/xs(2)
     input(3)=(input_data(i,3)-xm(3))/xs(3)
     result1 = net % output(input)
     do k=1,3
        yhat(k)=result1(k)*ys(k)+ym(k)
     enddo
     yhat(2) = 10.0**yhat(2) 
     Nw=yhat(1)
     Dm=yhat(2)
     Nw = 10**Nw !undo log, should be in m^-4
     Dm = Dm/1000. !convert to m ^4
     IWC = (Nw*(Dm)**4*1000*pi)/4.0**(4) ! the 1000 is density of water (kg/m^3)
     IWC = IWC*1000 !convert to g/m^3 
     output_data(i,1)=Nw
     output_data(i,2)=IWC
     output_data(i,3)=Dm
  end do
end subroutine call_keras


subroutine call_keras2(input_data,output_data,n)
  use keras_def
  integer :: n
  real :: input_data(n,3)
  real, intent(out) :: output_data(n,3)
  integer :: i
  real :: xm(3),xs(3),ym(3),ys(3), dwr, zku,zka
  real :: yhat(3)
  real :: IWC,Nw,Dm,pi
  xm=(/11.61196423,   1.62957636, -19.86298585/)
  xs=(/13.30597739,  2.15601794, 11.72233551/)
  ym=(/8.00513416, -0.26086958,  0.11921838/)
  ys=(/0.78638963, 0.21327072, 0.291193/)
  !----------Kuo's Tables-----------!
  xm=(/ 13.73949658,   1.58678511, -11.68680432/)
  xs=(/12.34509855,  3.46974159,  7.77404693/)
  ym(1:2)=(/ 0.34320764, -0.28632255/)
  ys(1:2)=(/1.11526479, 0.27130908/)
  pi=3.141592653589793
  do i=1,n
     zku=input_data(i,1)
     zka=input_data(i,2)
     if (zku<0) zku=0
     if (zka<0) zka=0
     dwr=zku-zka
     if(dwr<-1) dwr=-1
     if(dwr>11) dwr=11
     input2(1)=(zku-xm(1))/xs(1)
     input2(2)=(dwr-xm(2))/xs(2)
     input2(3)=(input_data(i,3)-xm(3))/xs(3)
     result2 = net % output(input2)
     do k=1,2
        yhat(k)=result2(k)*ys(k)+ym(k)
     enddo
     yhat(2) = 10.0**(yhat(2))
     Nw=yhat(1)
     Dm=yhat(2)
     Nw = 10**Nw !undo log, should be in m^-4
     Dm = Dm/1000. !convert to m ^4
     IWC = (Nw*(Dm)**4*1000*pi)/4.0**(4) ! the 1000 is density of water (kg/m^3)
     IWC = IWC*1000 !convert to g/m^3 
     output_data(i,1)=Nw
     output_data(i,2)=IWC
     output_data(i,3)=Dm
  end do
end subroutine call_keras2

subroutine call_keras21(zku,zka,nw,iwc,dm)
  use keras_def
  integer :: n
  integer :: i
  !real :: zku, zka
  real :: xm(3),xs(3),ym(3),ys(3), dwr, zku,zka
  real :: yhat(3)
  real :: IWC,Nw,Dm,pi
  xm=(/11.61196423,   1.62957636, -19.86298585/)
  xs=(/13.30597739,  2.15601794, 11.72233551/)
  ym=(/8.00513416, -0.26086958,  0.11921838/)
  ys=(/0.78638963, 0.21327072, 0.291193/)
  !----------Kuo's Tables-----------!
  xm=(/ 13.73949658,   1.58678511, -11.68680432/)
  xs=(/12.34509855,  3.46974159,  7.77404693/)
  ym(1:2)=(/ 0.34320764, -0.28632255/)
  ys(1:2)=(/1.11526479, 0.27130908/)
  pi=3.141592653589793
  !print*, 'call_keras21',zku,zka
  if (zku<0) zku=0
  if (zka<0) zka=0
  dwr=zku-zka
  if(dwr<-1) dwr=-1
  if(dwr>11) dwr=11
  input2(1)=(zku-xm(1))/xs(1)
  input2(2)=(dwr-xm(2))/xs(2)
  input2(3)=(-5-xm(3))/xs(3)
  result2 = net % output(input2)
  do k=1,2
     yhat(k)=result2(k)*ys(k)+ym(k)
  enddo
  yhat(2) = 10.0**(yhat(2))
  Nw=yhat(1)
  Dm=yhat(2)
  Nw = 10**Nw !undo log, should be in m^-4
  Dm = Dm/1000. !convert to m ^4
  IWC = (Nw*(Dm)**4*1000*pi)/4.0**(4) ! the 1000 is density of water (kg/m^3)
  IWC = IWC*1000 !convert to g/m^3 
  
end subroutine call_keras21



subroutine call_keras_chase(input_data,output_data,n)
  use keras_def
  integer :: n
  real :: input_data(n,3)
  real, intent(out) :: output_data(n,3)
  integer :: i
  real :: xm(3),xs(3),ym(3),ys(3), dwr, zku,zka
  real :: yhat(3)
  real :: IWC,Nw,Dm,pi
  xm=(/11.61196423,   1.62957636, -19.86298585/)
  xs=(/13.30597739,  2.15601794, 11.72233551/)
  ym=(/8.00513416, -0.26086958,  0.11921838/)
  ys=(/0.78638963, 0.21327072, 0.291193/)

  do i=1,n
     zku=input_data(i,1)
     zka=input_data(i,2)
     if (zku<0) zku=0
     if (zka<0) zka=0
     dwr=zku-zka
     if(dwr<-1) dwr=-1
     if(dwr>11) dwr=11
     input_ch(1)=(zku-xm(1))/xs(1)
     input_ch(2)=(dwr-xm(2))/xs(2)
     input_ch(3)=(input_data(i,3)-xm(3))/xs(3)
     result1 = net_ch % output(input_ch)
     do k=1,2
        yhat(k)=result1(k)*ys(k)+ym(k)
     enddo
     yhat(2) = 10.0**(yhat(2))
     Nw=yhat(1)
     Dm=yhat(2)
     Nw = 10**Nw !undo log, should be in m^-4
     Dm = Dm/1000. !convert to m ^4
     IWC = (Nw*(Dm)**4*1000*pi)/4.0**(4) ! the 1000 is density of water (kg/m^3)
     IWC = IWC*1000 !convert to g/m^3 
     output_data(i,1)=Nw
     output_data(i,2)=IWC
     output_data(i,3)=Dm
  end do
end subroutine call_keras_chase

