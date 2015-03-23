!*** ****** PLL Module ****** ***!
module PLL

!*** Vertical Mode ***!

Logical :: sigFF_z,oscFF_z,XOR_z
logical :: dSignalOld=.false.
logical :: dOscillOld=.false.

logical :: sigFDReset=.false.
logical :: oscFDSet=.false.

! double precision :: p1_z=0.0d0,p1o_z=0.0d0      !Phase Detector Output Signal (present and previous)
! double precision :: tau1_z=55.0E-6,exp1_z       !First Low Pass Filter time constant and Exp term
! 
! double precision :: p2_z=0.0d0,p2o_z=0.0d0      !Filtered signal (present and previous)
! double precision :: tau2_z=160.0E-6,exp2_z      !Second Low Pass Filter time constant and Exp term

double precision :: RC_z,pref_z,prefo_z

double precision :: ut_z=0.0d0,uto_z=0.0d0      !Filters output voltage (present and previous)
double precision :: uGain_z=6.0d0               !VCO voltage gain
double precision :: K0_z=220.0d0                !Volt/Hertz conversion factor

double precision :: DeltaF_z=0.0d0              !Frequency detune
double precision :: TotalF_z=0.0d0              !Total signal frequency
double precision :: IPhi_z=0.0d0                !Signal phase

double precision :: SinPLL_z=0.0d0,SinPLLo_z=0.0d0                  !Excitation reference (present and previous)
double precision :: refPhi_z,cosRefPhi_z,sinRefPhi_z                !Constant phase lag with its cos and sin
double precision :: xRefCos_z=0.0d0, xRefSin_z=0.0d0                !Excitation reference cos and sin parts
!*** ************* ***!

!*** Diethering Mode ***!
double precision :: pllyos
double precision :: phd_z                       !Phase detector output
double precision :: p1_y=0.0d0,p1o_y=0.0d0      !Phase Detector Output after 2 RC filters (present and previous)
double precision :: tau1_y=55.0E-6,exp1_y       !First Low Pass Filter time constant and Exp term

double precision :: p2_y=0.0d0,p2o_y=0.0d0      !Filtered signal (present and previous)
double precision :: tau2_y=160.0E-6,exp2_y      !Second Low Pass Filter time constant and Exp term

double precision :: RC_y
double precision :: phd_y,pref_y,prefo_y
double precision :: pref2_y, pref2o_y

double precision :: ut_y=0.0d0,uto_y=0.0d0      !Filters output voltage (present and previous)
double precision :: uGain_y=0.0d0               !VCO voltage gain
double precision :: K0_y=0.0d0                !Volt/Hertz conversion factor

double precision :: DeltaF_y=0.0d0              !Frequency detune
double precision :: TotalF_y=0.0d0              !Total signal frequency
double precision :: IPhi_y=0.0d0                !Signal phase

double precision :: SinPLL_y=0.0d0,SinPLLo_y=0.0d0      !Excitation reference (present and previous)
double precision :: refPhi_y,cosRefPhi_y,sinRefPhi_y    !Constant phase lag with its cos and sin
double precision :: xRefCos_y=0.0d0, xRefSin_y=0.0d0    !Excitation reference cos and sin parts
!*** ************* ***!

double precision :: bp1_z,bp1o_z,bp2_z,bp2o_z,but_z,buto_z,bDF_z,bTF_z,bSinPLL_z,bSinPLLo_z,bIphi_z
double precision :: bp1_y,bp1o_y,bp2_y,bp2o_y,but_y,buto_y,bDF_y,bTF_y,bSinPLL_y,bSinPLLo_y,bIphi_y


double precision :: RealF_z=0.0d0,RealFo_z=0.0d0,StartT_z=0.0d0
double precision :: RealF_y=0.0d0,RealFo_y=0.0d0,StartT_y=0.0d0
logical :: EXC_V=.true.
double precision, Allocatable :: Filter_z(:),Filtero_z(:)
double precision, Allocatable :: Filter_y(:),Filtero_y(:)
integer :: FilterOrder_z = 5
integer :: FilterOrder_y = 5

contains


!*** PLL Initializer ***
!   Input variables:
!   (DBL)dt         timestep
subroutine PLLInit(dt)
  use consts
  implicit none
  double precision, intent(in) :: dt

  !allocate filters
  Allocate(Filter_z(0:FilterOrder_z))
  Allocate(Filtero_z(0:FilterOrder_z))
  Filter_z=0.0d0
  Filtero_z=0.0d0
  Allocate(Filter_y(0:FilterOrder_y))
  Allocate(Filtero_y(0:FilterOrder_y))
  Filter_y=0.0d0
  Filtero_y=0.0d0
  
  write(*,'(A)') 'Frequency Detector circuits initialized.'
  write(*,'(A)') 'FD parameters:'
  write(*,'(A F12.3 A I3.1 A)')'Filter vertical : ',1.0d0/(RC_z*dPi),' Hz (order ',FilterOrder_z,')'
  write(*,'(A F12.3 A I3.1 A)')'Filter lateral  : ',1.0d0/(RC_y*dPi),' Hz (order ',FilterOrder_y,')'
  write(*,*)
  
  
  !write(*,*) refPhi_z,SinRefPhi_z,CosRefPhi_z
  !stop

end subroutine PLLInit





!*** ****** Frequency Detector ****** ***!
!   Input variables:
!   (DBL) dt          timestep
!   Description:
!   These 2 routines are great and pwn a lot!
!   They measure the time delay between 2 signal peaks
!   and get the frequency. Since this would lead to a
!   fuzzy noisy signal, the lowpass filter is applied twice.
!
subroutine FM_UpdateZ(dt)
  use consts
  use cantilever
  use timer
  implicit none
  double precision,intent(in) :: dt
  logical :: tik
  integer :: i
  
  tik=(Xo(3)>Xoo(3)).and.(Xo(3)>X(3))   !detect the peak
  
  if(tik)then     !if there is a tick,
    Filter_z(0)=(1.0d0/(t-StartT_z))-F1    !compute the frequency from the time lag
    StartT_z=t                        !reset the time lag reference cursor
  endif

  do i=1,FilterOrder_z
    Filter_z(i)=dt/(RC_z+dt)*Filter_z(i-1) + (1.0-dt/(RC_z+dt))*Filtero_z(i)
    Filtero_z(i)=Filter_z(i)    !back collector
  enddo
  
  DeltaF_z=Filter_z(FilterOrder_z)-RealFo_z     !set the frequency shift
  TotalF_z=F1+DeltaF_z 
  
  
end subroutine FM_UpdateZ

subroutine FM_UpdateY(dt)
  use consts
  use cantilever
  use timer
  implicit none
  double precision,intent(in) :: dt
  logical :: tik
  integer :: i
  
  tik=(Xo(2)>Xoo(2)).and.(Xo(2)>X(2))   !detect the peak
  
  if(tik)then     !if there is a tick,
    Filter_y(0)=(1.0d0/(t-StartT_y))-F2    !compute the frequency from the time lag
    StartT_y=t                        !reset the time lag reference cursor
  endif

  do i=1,FilterOrder_y
    Filter_y(i)=dt/(RC_y+dt)*Filter_y(i-1) + (1.0-dt/(RC_y+dt))*Filtero_y(i)
    Filtero_y(i)=Filter_y(i)    !back collector
  enddo
  
  DeltaF_y=Filter_y(FilterOrder_y)-RealFo_y     !set the frequency shift
  TotalF_y=F2+DeltaF_y                          !set the total frequency (for the AGC)
  
end subroutine FM_UpdateY



!*** ****** Excitation signal Generator ****** ***!
!   Description:
!   Generates the excitation signals SinPLL for Z and Y,
!   If EXC_V is .TRUE. the excitation is a digital signal
!   locked to the velocity of the cantilever.
!   An experimental guy told me itz common practice in real cases,
!   since probably they too have problems with that fucking PLLs.
! 
subroutine Excite_Update()
  use cantilever
  implicit none
  
  !set the excitation signal locked to the velocity
  if(EXC_V) then
    
    SinPLL_z=-1.0d0
    if(V(3)>0.0d0) SinPLL_z=1.0d0
    
    SinPLL_y=-1.0d0
    if(V(2)>0.0d0) SinPLL_y=1.0d0
  
  endif

end subroutine Excite_Update







end module PLL