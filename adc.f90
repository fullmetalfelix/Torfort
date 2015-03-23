module ADC

!*** Tip Height parameters
double precision :: zStart=50.0e-9      !Initial tip height
double precision :: zLand=10.5e-9       !Scanning tip height
double precision :: tLand=0.0d0         !Approaching time


!*** Dynamic Distance regulation
logical :: ADCON_z=.false.
double precision :: ADCSet_z=-30.0d0    !target frequency shift
double precision :: ADC_z=0.0d0         !Current tip height
double precision :: ADCo_z=50.0e-9      !Offset
double precision :: ADCp_z=0.0d0        !Proportional part
double precision :: ADCi_z=0.0d0        !Integral part

double precision :: AKP_z=0.1d0        !Proportional Gain
double precision :: AKI_z=1.0d0        !Integral Gain

!convergence parameters
double precision :: ARelaxTol_z=1.0d0        !tolerance for Frequency shift relaxation (Hz)
double precision :: ARelaxTimeSpan_z=0.01    !Time span during which the Amplitude has to stay within the tolerance (s)
double precision :: ARelaxTimer=0.0d0        !the relaxation timer

double precision :: zCrashLimit = dble(3.0e-10)

!backup variables
double precision :: bADCSet_z,bADC_z,bADCi_z

contains

!*** Initialization routine ***
subroutine ADCInit()
  implicit none

  
  ADC_z=zStart    !set the zPiezo at the starting point
  ADCo_z=ADC_z    !Set the offset at the same point for no real reason
  
  write(*,'(A)') 'ADC circuits initialized.'
  write(*,'(A)') 'ADC Feedback parameters:'
  write(*,'(A ES14.8)')'Offset: ',ADC_z
  write(*,'(A ES14.8)')'P Gain: ',AKP_z
  write(*,'(A ES14.8)')'I Gain: ',AKI_z
  write(*,*)
  
end subroutine ADCInit


subroutine ADC_Tune()
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
	write(*,*)'End tuning ADC.'
	exit
      case('enabled')
	write(*,'(A)',advance='no')'ADC circuit is now '
	ADCON_z=.not.ADCON_z
	if(ADCON_z) then
	  write(*,*)'ON'
	else
	  write(*,*)'OFF'
	endif
      case('pgain')
        write(*,'(A,ES16.8)')'The proportional gain is ',AKP_z
	write(*,'(A)')'new value:'
        read(*,*)AKP_z
        write(*,'(A,ES16.8)')'proportional gain is now ',AKP_z
      case('igain')
	write(*,'(A,ES16.8)')'The integral gain is ',AKI_z
	write(*,'(A)')'new value:'
	read(*,*)AKI_z
	write(*,'(A,ES16.8)')'integral gain is now ',AKI_z
      case DEFAULT
	write(*,'(A)')'Unknown command!'
      end select
    
  enddo
  
  
end subroutine ADC_Tune



!*** ADC circuit updater ***
!   Input variables:
!   (DBL)dFreq      frequency shift coming from the PLL
!   (DBL)dt         timestep constant
!   Description:
!   should be clear!
subroutine ADC_UpdateZ(dFreq,dt)
  use cantilever
  implicit none
  double precision, intent(in) :: dFreq,dt
  !double precision, intent(inout):: outposz
  
  !*** update the gain
  if(ADCON_z) then
    ADCi_z=ADCi_z+AKI_z*(ADCSet_z-dFreq)*dt   !Integral part
    ADCp_z=AKP_z*(ADCSet_z-dFreq)             !Proportional part
    ADC_z=ADCo_z+ADCi_z+ADCp_z
    HolderPosition(3)=ADC_z                   !set the HolderPosition Z to the ADC_z value
  endif
  
end subroutine ADC_UpdateZ


!*** ADC Convergence checker ***
!   Input variables:
!   (DBL)dFreq      frequency shift coming from the PLL
!   (DBL)dt         timestep
!   Description:
!   when the frequency shift dFreq is the one we set within some
!   tolerance ARelaxTol_z and stays like that for a time
!   ARelaxTimeSpan_z this function will return .TRUE.
logical function AConvergenceCheck(dFreq,dt)
  implicit none
  double precision, intent(in) :: dFreq,dt
  AConvergenceCheck=.false.
  
  !*** check if the aplitude is within tolerance
  if(ABS(dFreq-ADCset_z) <= ARelaxTol_z) then
    ARelaxTimer=ARelaxTimer+dt       !increase the timer
  else
    ARelaxTimer=0.0d0                  !reset the timer if out of tolerance
  endif
  
  !*** check if it was stable within tolerance for the required time span
  if(ARelaxTimer>=ARelaxTimeSpan_z)    AConvergenceCheck=.true.  
  
  return
  
end function AConvergenceCheck




end module ADC