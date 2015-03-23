!*** ****** AGC Module ****** ***!
module AGC

  !*** vertical variables ***!
  double precision :: Aset_z=10.0e-9    !Amplitude Setpoint
  double precision :: Aex0_z=1.0d0			!excitation at relaxed state
  double precision :: FB_z=0.0d0        !Driving force intensity
  double precision :: FBo_z=0.010d0     !Feedback offset
  double precision :: FBp_z=0.000d0     !Feedback Proportional part
  double precision :: FBi_z=0.0d0       !Feedback Integral part
  double precision :: KP_z=0.10d0       !Proportional feedback gain
  double precision :: KI_z=200000.000d0      !Integral feedback gain
  
  double precision :: tauD_z=2.00E-3              !Amplitude detector time constant
  double precision :: expD_z=0.0d0                !Exp(-tstep/TauD_z)
  double precision :: At_Z=1.0d0, Ato_z=1.0d0     !measured amplitude present and previous
  double precision :: wv_z=0.0d0, wvo_z=0.0d0     !Module of the Signal present and previous
  logical :: AGCON_z=.true.
  logical :: ExciteON_z=.true.
  
  !convergence parameters ***!
  double precision :: RelaxTol_z=2.00e-12         !tolerance for amplitude relaxation (m) (vertical)
  double precision :: RelaxTol_y=0.01e-12         !tolerance for amplitude relaxation (m) (lateral)
  double precision :: RelaxTimeSpan_z=0.01        !Time span during which the Amplitude has to stay within the tolerance (s)
  double precision :: RelaxTimer=0.0d0            !the relaxation timer
  
  !*** ****************** ***!
  
  !*** Diethering variables ***!
  double precision :: Aset_y=1.0e-10    !Amplitude Setpoint
  double precision :: Aex0_y=0.0d0			!excitation at relaxed state
  double precision :: FB_y=0.0d0        !Driving force intensity
  double precision :: FBo_y=0.00d0     !Feedback offset
  double precision :: FBp_y=0.00d0     !Feedback Proportional part
  double precision :: FBi_y=0.0d0       !Feedback Integral part
  double precision :: KP_y=0.25d0       !Proportional feedback gain
  double precision :: KI_y=0.000d0      !Integral feedback gain
  
  double precision :: tauD_y=1.0E-5               !Amplitude detector time constant
  double precision :: expD_y=0.0d0                !Exp(-tstep/TauD_z)F1=150000.0d0
  double precision :: At_y=0.0d0, Ato_y=0.0d0     !measured amplitude present and previous
  double precision :: wv_y=0.0d0, wvo_y=0.0d0     !Module of the Signal present and previous
  
  !Amplitude detector parameters
  integer :: AFilterOrder_z = 3						!number of RC filters (vertical)
  integer :: AFilterOrder_y = 3						!number of RC filters (lateral)
  double precision :: ARC_z, ARC_y				!RC value for vertical and lateral filters
  double precision :: Aerr_z,AErr_y				!error correction offset for the filters
  
  double precision, Allocatable :: AFilter_z(:),AFiltero_z(:)			!the filters (vertical)
  double precision, Allocatable :: AFilter_y(:),AFiltero_y(:)			!the filters (lateral)
  
  logical :: AGCON_y=.false.
  logical :: ExciteON_y=.true.
  logical :: ExciteSimple=.false.
  !*** ******************** ***!
  
  !backup variables
  double precision :: bFB_z, bAt_z, bAto_z, bwv_z, bwvo_z
  double precision :: bFB_y, bAt_y, bAto_y, bwv_y, bwvo_y
  
double precision :: Dissipation_z=0.0d0,DissiFilter_z(0:2)=0.0d0,Fzo=0.0d0,zo=0.0d0,DissiFiltero_z(0:2)=0.0d0
double precision :: DissipationI_z=0.0d0
double precision :: DFilterRC_z=0.001590d0,tmpDissi_z=0.0d0
double precision :: Dissipation_y=0.0d0,DissiFilter_y(0:2)=0.0d0,DissiFiltero_y(0:2)=0.0d0
double precision :: DissipationI_y=0.0d0
double precision :: DFilterRC_y=0.001590d0,tmpDissi_y=0.0d0
double precision :: E0_z, E0_y

contains


!*** AGC Initializer ***
!   Input variables:
!   (DBL)dt         timestep
subroutine AGCInit(dt)
  use Consts
  use cantilever
  implicit none
  double precision,intent(in) :: dt
  
  Aerr_z=1.0d0
  Aerr_y=1.0d0
  allocate(AFilter_z(0:AFilterOrder_z),AFilter_y(0:AFilterOrder_y))
  allocate(AFiltero_z(0:AFilterOrder_z),AFiltero_y(0:AFilterOrder_y))
  AFiltero_z=0.0d0
  AFiltero_y=0.0d0
  
  if(FBo_z<0.0d0) FBo_z=k_z*Aset_z/Q1       !if a negative offset is given, autocalculate the right one
  if(FBo_y<0.0d0) FBo_y=k_y*Aset_y/Q2
	
	Aex0_z=FBo_z
	Aex0_y=FBo_y
	
  write(*,'(A)') 'AGC circuits initialized.'
  write(*,'(A F8.3 A I3.1 A)') 'Amplitude.Z detection bandwidth ',1.0d0/(dPi*ARC_z), ' Hz (order ',AFilterOrder_z,')'
  write(*,'(A F8.3 A I3.1 A)') 'Amplitude.Y detection bandwidth ',1.0d0/(dPi*ARC_y), ' Hz (order ',AFilterOrder_y,')'
  write(*,*)
  write(*,'(A)')'AGC Feedback parameters (vertical and lateral):'
  write(*,'(A ES14.8 A ES14.8)')'Offset: ',FBo_z,'  ',FBo_y
  write(*,'(A ES14.8 A ES14.8)')'P Gain: ',KP_z,'  ',KP_y
  write(*,'(A ES14.8 A ES14.8)')'I Gain: ',KI_z,'  ',KI_y
  write(*,*)
  
	E0_z = Pi*K_z*Aset_z**2/Q1
	E0_y = Pi*K_y*Aset_y**2/Q2
  
end subroutine AGCInit

! subroutine AGCCalcTau(dt)
!   implicit none
!   double precision, intent(in) :: dt
!   
!   expD_z=exp(-dt/tauD_z)
!   expD_y=exp(-dt/tauD_y)
!   
!   skipNext_z=.true.
!   skipNext_y=.true.
!   
! end subroutine AGCCalcTau

subroutine AGCy_Tune()
  implicit none
  character*80 :: cmd
  
  write(*,'(A)')'ADC Tuning'
  write(*,'(A)')'Type done to save and exit.'
  do
    write(*,'(A)')'What parameter do you wish to tune?'
    read(*,'(A)')cmd
    
    cmd=trim(adjustl(cmd))
    
    select case (cmd)
      case('help')
	write(*,*)'You can tune: enabled, pgain, igain.'
      case('done')
	write(*,*)'End tuning AGC y.'
	exit
      case('enabled')
	write(*,'(A)',advance='no')'AGC y circuit is now '
	AGCON_y=.not.AGCON_y
	if(AGCON_y) then
	  write(*,*)'ON'
	else
	  write(*,*)'OFF'
	endif
      case('pgain')
        write(*,'(A,ES16.8)')'The proportional gain is ',KP_y
	write(*,'(A)')'new value:'
        read(*,*)KP_y
        write(*,'(A,ES16.8)')'proportional gain is now ',KP_y
      case('igain')
	write(*,'(A,ES16.8)')'The integral gain is ',KI_y
	write(*,'(A)')'new value:'
	read(*,*)KI_y
	write(*,'(A,ES16.8)')'integral gain is now ',KI_y
      case DEFAULT
	write(*,'(A)')'Unknown command!'
      end select
    
  enddo
  
  
end subroutine AGCy_Tune


subroutine GetDissi(dt)
	use cantilever
	use consts
	implicit none
	double precision, intent(in) :: dt
	integer :: i
	
	!Aex_z=FB_z*Q1/(K_z*Sqrt(1.0d0-0.250d0*Q1**2))
	!Dissipation_z=J2eV*E0_z*(FB_z-Aex0_z)/Aex0_z
	!Dissipation_y=J2eV*E0_y*(FB_y-Aex0_y)/Aex0_y
	Dissipation_z=J2eV*E0_z*(FB_z-Aex0_z)/Aex0_z
	
! 	do i=1,2					!process the filters
! 		DissiFilter(i)=DissiFilter(i-1)*(dt/(DFilterRC+dt)) + (1.0d0-dt/(DFilterRC+dt))*DissiFiltero(i)
! 		DissiFiltero(i)=DissiFilter(i)    !back collector
! 	enddo
! 	Dissipation_z=DissiFilter(2)
! 		if(phase>=1) then
! 			write(*,*)DissiFilter(0)*(dt/(DFilterRC+dt)),(1.0d0-dt/(DFilterRC+dt))
! 			stop
! 		endif
		
! 	write(*,*)DissiFilter(:)
	
end subroutine GetDissi


subroutine DissiCalc(dt,Freq_z,Freq_y,zTick,yTick)
	use cantilever
	use consts
	implicit none
	double precision, intent(in) :: dt,Freq_z,Freq_y
	logical, intent(in) :: zTick,yTick
	double precision :: dtmp_z,tfq_z
	double precision :: dtmp_y,tfq_y
	double precision :: zFmatrix(3,3), zFb(3),ipiv(3)
	integer :: i
	
	
	Dissipation_z=J2eV*E0_z*(FB_z-Aex0_z)/Aex0_z				!compute the dissipation from the feedback signal
	Dissipation_y=J2eV*E0_y*(FB_y-Aex0_y)/Aex0_y
	
	!dissipation from Integrate[F.ds]

	dtmp_z = 0.50d0*(Ftipsurf(3)+Ftipsurfo(3))*(x(3)-xo(3))*J2eV
	tmpDissi_z = tmpDissi_z+dtmp_z							!update integral
	
	if(zTick) then			!if the tip did a full cycle...
		DissiFilter_z(0)=tmpDissi_z
		tmpDissi_z=0.0d0
		!write(888,*)DissiFilter(0)
		tfq_z=1.0d0/Freq_z
		do i=1,2					!process the filters
			DissiFilter_z(i)=DissiFilter_z(i-1)*(tfq_z/(DFilterRC_z+tfq_z)) + (1.0d0-tfq_z/(DFilterRC_z+tfq_z))*DissiFiltero_z(i)
			DissiFiltero_z(i)=DissiFilter_z(i)    !back collector
		enddo
		DissipationI_z=DissiFilter_z(2)
	endif
	!----------------------------------------------------------------
	dtmp_y = 0.50d0*(Ftipsurf(2)+Ftipsurfo(2))*(x(2)-xo(2))*J2eV
	tmpDissi_y = tmpDissi_y+dtmp_y							!update integral
	
	if(yTick) then			!if the tip did a full cycle...
		DissiFilter_y(0)=tmpDissi_y
		tmpDissi_y=0.0d0
		!write(888,*)DissiFilter(0)
		tfq_y=1.0d0/Freq_y
		do i=1,2					!process the filters
			DissiFilter_y(i)=DissiFilter_y(i-1)*(tfq_y/(DFilterRC_y+tfq_y)) + (1.0d0-tfq_y/(DFilterRC_y+tfq_y))*DissiFiltero_y(i)
			DissiFiltero_y(i)=DissiFilter_y(i)    !back collector
		enddo
		DissipationI_y=DissiFilter_y(2)
	endif
	!zo=dtmp_z
	
end subroutine DissiCalc


!*** AGC Z Circuit Updater ***
!   Input variables:
!   (DBL)SignalIn   Signal coming from the sensor
!   (DBL)dt         timestep
!   (LOG)Tick       oscillation cycle indicator
!   (DBL)Freq       PLL measured frequency
!   Description:
!   first filters the signal to get the amplitude (At_z)
!   then, if the tip has completed a cycle, computes
!   the excitation needed to keep the preset amplitude (Aset_z)
subroutine AGC_UpdateZ(SignalIn,Freq)
  use consts
  implicit none
  double precision, intent(in) :: SignalIn,Freq
  double precision :: tfq
  integer :: i
  
  tfq=1.0d0/Freq
  !*** update the amplitude - Low pass Filter RC-mode - tick only mode
  AFilter_z(0)=SignalIn
	do i=1,AFilterOrder_z
		AFilter_z(i)=(tfq/(ARC_z+tfq))*AFilter_z(i-1) + (1.0d0-tfq/(ARC_z+tfq))*AFiltero_z(i)
		AFiltero_z(i)=AFilter_z(i)    !back collector
	enddo
	
	At_z=AFilter_z(AFilterOrder_z)*Aerr_z
	
  
  !*** update the gain
  if(AGCON_z) then
    FBi_z=FBi_z+KI_z*(Aset_z-At_Z)/Freq               !Integral part
    FBp_z=KP_z*(Aset_z-At_Z)                             !Proportional part
    FB_z=FBi_z+FBp_z+FBo_z
  endif
  
  
  
  !*** update old values
!   wvo_z=wv_z
!   Ato_z=At_z
  
  
end subroutine AGC_UpdateZ




!*** AGC Y Circuit Updater ***
!   Input variables:
!   (DBL)SignalIn   Signal coming from the sensor
!   (DBL)dt         timestep
!   (LOG)Tick       oscillation cycle indicator
!   (DBL)Freq       PLL measured frequency
!   Description:
!   first filters the signal to get the amplitude (At_y)
!   then compute the excitation needed to keep the preset
!   amplitude (Aset_y)

subroutine AGC_UpdateY(SignalIn,Freq)
  use consts
  implicit none
  double precision, intent(in) :: SignalIn,Freq
  double precision :: tfq
  integer :: i
  
  !*** update the amplitude - Low pass Filter RC-mode
  tfq=1.0d0/Freq
  AFilter_y(0)=SignalIn
	do i=1,AFilterOrder_y
		AFilter_y(i)=AFilter_y(i-1)*(tfq/(ARC_y+tfq)) + (1.0d0-tfq/(ARC_y+tfq))*AFiltero_y(i)
		AFiltero_y(i)=AFilter_y(i)    !back collector
	enddo
	
  At_y=AFilter_y(AFilterOrder_y)*Aerr_y     !set the frequency shift
  
  
  !*** update the gain
  if(AGCON_y) then
    FBi_y=FBi_y+KI_y*(Aset_y-At_y)/Freq                  !Integral part
    FBp_y=KP_y*(Aset_y-At_y)                             !Proportional part
    FB_y=FBi_y+FBp_y+FBo_y
  endif
  
  !*** update old values
!   wvo_y=wv_y
!   Ato_y=At_y
  
  
end subroutine AGC_UpdateY



!*** AGC Convergence checker ***
!   Input variables:
!   (DBL)dt         timestep
!   Description:
!   when the amplitudes At_z and At_y are the ones we set within some
!   tolerance RelaxTol_z and stays like that for a time
!   RelaxTimeSpan_z this function will return .TRUE.
logical function AmpConvergenceCheck(dt)
  use Consts
  use Cantilever
  implicit none
  double precision, intent(in) :: dt
  logical :: Conv_z, Conv_y
	
  AmpConvergenceCheck=.false.
	
	Conv_z = .true.
	Conv_y = .true.
	
	!if the feedbacks are activated check for convergence
	if(AGCON_z) Conv_z = (ABS(At_z-Aset_z) <= RelaxTol_z)
	if(AGCON_y) Conv_y = (ABS(At_y-Aset_y) <= RelaxTol_y)
	
  !*** check if the aplitude is within tolerance
  if( Conv_z .and. Conv_y ) then     !.and. (ABS(At_y-Aset_y) <= RelaxTol_y)
    RelaxTimer=RelaxTimer+dt          !increase the timer
  else
    RelaxTimer=0.0d0                  !reset the timer if out of tolerance
  endif

  !*** check if it was stable within tolerance for the required time span
  if(RelaxTimer>=RelaxTimeSpan_z) then
    AmpConvergenceCheck=.true.  
  endif

  return
  
end function AmpConvergenceCheck



end module AGC