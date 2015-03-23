!*** ****** PLL Torfort ****** ***!
Program Torfort

  Implicit None
  character*80 :: command
  integer :: wnum
  
  call MainINIT() 

  
  
contains


Subroutine MainINIT()
  use PLL
  use AGC
  use ADC
  use cantilever
  use input
  use timer
  use forcer
  use imager
  use output
  implicit none
	
	
  call LoadALL()
	
  !*** setting up the force field

  call ForceInit()
	
  call CLVInit()        !initialize the cantilever
  call AGCInit(tstep)   !initialize the Automatic Gain Control
  call PLLInit(tstep)   !initialize the PLL
  call ADCInit()        !initialize the Automatic Distance Control
  call ImagerInit()     !initialize the imager unit
	
  call InitOutput()     !initialize the output manager
	
	
	!main process
	
	if(WorkMode == 7) then
		call TheorImage()
		stop
	endif
	
  t=0.0d0
  !call WriteLateralMapXY(dble(6.9e-10))
!    call writeLateralMapYZ(0.0d0)
!    call writeForceLine(0.0d0,dble(5.0e-10))
!   call writeForceLine(0.0d0, dble(5.5e-10) )
!   stop
  
	!IgnoreDissip=.true.
 	
  	!call OscillDBG()
  	!stop
 	
	call InitFreeRelax()
  call FreeRelax()
  
  call InitRelax()
  call Relax()
	
  Phase=1     !use phase 1 for first approach
  !*** Approach the surface ***
  if(.true.) then
		if(tLand==0.0d0) then
			call Approach(.false.)
		else
			call Approach(.true.)
		endif
	else
		!insta approach!
		IntON_z=.true.
		IntON_y=.true.
		holderposition(3)=zLand
		call freeevo(0.10d0,.true.)
	endif

	!IgnoreDissip=.true.
	!phase=3
  call FreeEvo(0.30d0,.true.)	
  write(*,*)'relaxing all done!'
!   call ProbeSteadyDissip()
!   stop
!   call FreeEvo(1.00d0,.false.)
!   open(file='point.dat',unit=1005)
!   write(1005,*)holderposition(1:2),holderposition(3),dissipation_z,dissipationI_z,dissipation_y,dissipationI_y
!   close(1005)
!   stop

  !for single scans!
  !call ZSpectrLine(dble(11.28e-10), 10000, dble(0.05e-10), 100)
	!call ScanLine(dble(2.8e-9), 1000)
	!stop
	
	
	write(*,*)'RELAXATION DONE!'
  if(WorkMode == 1) call Commander()
  if(WorkMode == 2) call ScanLine()
  if(WorkMode == 3) call ScanImage()
  if(WorkMode == 4) call ZSpectrLine()
  if(WorkMode == 5) call SinglePoint()
  if(WorkMode == 6) call ForceCurve_Static()
  
  if(WorkMode == 8) call ProbeSteadyDissip()
	
end subroutine MainINIT


subroutine ProbeSteadyDissip()
	use cantilever
	use AGC
	use output
	use timer
	use PLL
	use Imager
	implicit none
	double precision :: dissip(3),z
	integer :: cnt=0
	
	open(file='steady.out',unit=101)
	
	cnt=0
	!scan z on Na
	Phase=9
	!HolderPosition(3)=dble(4.3e-10)+ASet_z
	!call FreeEvo(0.10d0,.true.)
	write(*,*)'Relaxed on Na.'
	
	do
		write(101,*)holderposition(3)-Aset_z,At_z,Dissipation_z,DissipationI_z, DeltaF_z	!write data
		write(*,*)  holderposition(3)-Aset_z,At_z,Dissipation_z,DissipationI_z, DeltaF_z	!write data
		HolderPosition(3)=HolderPosition(3)+ ZS_step						!move holder
		cnt = cnt + 1
		call FreeEvo(0.30d0,.false.)								!relax
		!if(HolderPosition(3)-ASet_z>=dble(5.8e-10)) exit
		if(cnt >= ZS_npts) exit
	enddo
	
	close(101)
	
end subroutine ProbeSteadyDissip

!*** ****** Phase 0: oscillation relaxation
!   Description:
!   This first relaxation just let some time pass with no external forces,
!   no excitation signal and no feedbacks. And also no friction!
!   The PLL is supposed to lock (plz ffs lock).
!   This is used to correct the gain loss from RC-filters
!
Subroutine InitFreeRelax()
  Use Consts
  Use Cantilever
  Use AGC
  use PLL
  use ADC
  implicit none

  IntON_z=.false.
  IntON_y=.false.
  ADCON_z=.false.
  AGCON_z=.false.
  AGCON_y=.false.
  ExciteON_z=.false.
  ExciteON_y=.false.
  ExciteSimple=.false.
  
  !initial position
  x(1)=0.0d0
  x(2)=Aset_y
  x(3)=Aset_z
  
  Gamma1=0.0d0
  Gamma2=0.0d0
  
  Xo=X
  Xoo=X
  
  v=0.0d0
  Call GetForce(A)
  
  HolderPosition(3)=ADC_z
  
  
end subroutine InitFreeRelax

subroutine FreeRelax()
  use consts
  use cantilever
  use timer
  use output
  use PLL
  use AGC
  implicit none
  integer :: cnt=0
  logical :: check=.false.
  
  t=0.0d0
  
  write(*,'(A)',advance='no')'Gauging...'
  do
    t=t+tstep
    cnt=cnt+1
    
    call Verlet(tstep)
    call UpdateCircuits()
    
    if(t>=0.050d0) exit
  enddo

  !now set the F detector offset
  RealFo_z=DeltaF_z
  RealFo_y=DeltaF_y
  
  !set the AGC Gain
  Aerr_z=Aset_z/At_z
  Aerr_y=Aset_y/At_y
  
  write(*,*)'done!'
  
 end subroutine FreeRelax


!*** ******* Phase 1: electronics relaxation
!   Description:
!   Switch on the circuits and relax them.
!   AGC will be on, and there is the excitation signal.
!
Subroutine InitRelax()
  Use Consts
  Use Cantilever
  use ADC
  Use AGC
  use PLL
  Implicit None

  Gamma1=W1/Q1
  Gamma2=W2/Q2

  IntON_z=.false.
  IntON_y=.false.
  
  ADCON_z=.false.
  ExciteSimple=.false.


	x(1) = 0.0d0
  !*** Vertical Initial Conditions
  if(Oscill_z) then   !switch on if it has to oscillate
  	ExciteON_z = .true.
  	AGCON_z    = .true.
  	FB_z       = FBo_z!/2.0d0
  	x(3)       = Aset_z
  else
  	ExciteON_z = .false.
  	AGCON_z    = .false.
  	FB_z       = 0.0d0
  	x(3)       = 0.0d0
  endif
  
  !*** Lateral Initial Conditions
  if(Oscill_y) then   !switch on if it has to oscillate
  	ExciteON_y = .true.
  	AGCON_y    = .true.
  	FB_y       = FBo_y
  	x(2)       = Aset_y
  else
  	ExciteON_y = .false.
  	AGCON_y    = .false.
  	FB_y       = 0.0d0
  	x(2)       = 0.0d0
  endif
  
  Xo=X
  Xoo=X
  v=0.0d0
  
  Call GetForce(A)
  
end subroutine InitRelax

subroutine Relax()
  use consts
  use timer
  use cantilever
  !use PLL
  use AGC
  !use ADC
  use output
  !use Imager
  implicit none
  double precision :: tEnd=2.0d0
  integer :: cnt=0
  logical :: Check=.false.

  write(*,*) "Relaxing the cantilever's oscillation..."
  Phase=0   !phase for output
	t=0.0d0
	
	open(file='p_relax.out',unit=888)
	
  do !T=tstep,Tend,tStep
    t=t+tstep
    cnt=cnt+1

    call Verlet(tstep)       !perform a Verlet step
    call UpdateCircuits()    !update electronics

    !*** write data every 'PointBuffer' time steps
    if((cnt==PointBuffer).and.DumpPhase(2)) then
      Call WriteAll_f(888)
      !write(*,*)t,At_z,DeltaF_z,FB_z
      cnt=0
    endif
    !write(111,*)t,x(:)
  
    if(t>=0.50d0) check = AmpConvergenceCheck(tstep)    !some time needs to pass before checking
    
    if(check .and. .true.) then    !if relaxed end this loop
      write(*,*) 'The oscillation amplitude was stable for the required time.'
      write(*,*) 'Relaxed in:',t, 'seconds'
      exit
    endif
  
  end do
	
	close(888)
	
	!set the excitation setpoints
	Aex0_z=FB_z
	Aex0_y=FB_y
	!write(*,*)'ex set',Aex0_z
	!stop
	
end subroutine Relax


!*** ****** Phase 2: Tip Approach ****** ***!
!   Input variables:
!   (LGC) DistMode          approach method
!   Description:
!   move the tip down until dF = dF_set (if DistMode=.false.) or
!   until it reaches zLand (when DistMode=.true.)
!
subroutine Approach(DistMode)
  use consts
  use timer
  use cantilever
  use PLL
  use AGC
  use ADC
  use output
  implicit none
  logical, intent(in) :: DistMode    !true for distance controlled landing, false for realistic ADC controlled landing
  double precision :: tEnd,tstart
  integer :: cnt=0
  logical :: Check
  logical :: approaching=.true., relaxing=.false.
  CHARACTER(LEN=30) :: FRMT= '(3F10.3)'
  
  !switch on the interaction (for the active oscillations)!
  IntON_z = Oscill_z
  IntON_y = Oscill_y

  tstart=t
  tEnd=t+3.0d0
  
  !*** info message ***
  write(*,'(A)',advance='no') 'Approaching the tip '
  if(DistMode)then
    write(*,*)'(Direct Distance Control)...'
  else
    write(*,*)'(Automatic Distance Control)...'
    ADCON_z=.true.        !switch on the ADC
  endif
  !*** ************ ***
  
  open(file='p_approach.out',unit=888)
  
  do 
    t=t+tstep
    cnt=cnt+1
    
    call Verlet(tstep)    !perform a Verlet step
    
    call UpdateCircuits() !operate the circuits
    
    !*** write data every 'PointBuffer' time steps
    if((cnt==PointBuffer).and.DumpPhase(3)) then
      Call WriteAll_f(888)
      write(*,*) DeltaF_z,At_z,ADC_z
      !write(*,*)(FB_z-Aex0_z)/Aex0_z
      cnt=0
    endif
  
    
    if(Approaching) then
      
      !*** select the right method for approaching...
      if(DistMode)then
      
        if(abs(ADC_z-zLand)>=1.0e-12) then      !if the tip is above the desired position...
          ADC_z=ADC_z+((zLand-zStart)/tland)*tstep     !move it down a bit
          HolderPosition(3)=ADC_z                      !update the holder position
        else
          Approaching=.false.                          !else, the approaching phase ends
          HolderPosition(3)=zLand
          ADC_z=zLand
          write(*,*)'The tip approached safely.'       !and the oscillation will be relaxed a bit again
          write(*,'(A)',advance='no') 'Now relaxing the amplitude again...'
        endif
        
      else    !*** if the ADC is used
        
        check = AConvergenceCheck(DeltaF_z,tstep)  !*** Check ADC convergence for approach
        if(check) then                             !if relaxed end this loop
          write(*,'(A F8.3 A)') 'The Tip approached safely in ',t-tstart, ' seconds'
          write(*,'(A21 F10.4 A)') 'The ZPiezo is now at ',ADC_z/nano,' nm'
          write(*,'(A)',advance='no') 'Now relaxing the amplitude again...'
          Approaching=.false.             !*** the approach is done
        endif
    
      endif
      
    else            !if the approach is finished
      
      !*** Relaxation criteria
      check=.true.
      if(.not.DistMode)then
        check=AConvergenceCheck(DeltaF_z,tstep)   !if we are in ADC mode relax the ZPiezo
      else
        check=(t>=tstart+tLand+0.050d0)           !if in DDC just let some time pass (0.05s) to relax the PLL
      endif
      check = check .and. AmpConvergenceCheck(tstep) !Also waits for the amplitude to relax
      !*** relaxation criteria are checked

      if(check) then
        write(*,*) 'Done!'
        write(*,'(A F8.4 A)') 'Amplitude is ',At_z/nano,' nm (vertical)'
        write(*,'(A F8.4 A)') 'Amplitude is ',At_y/ang,' å (lateral)'
        write(*,'(A F12.4 A)') 'Frequency Shift is ',DeltaF_z,' Hz (vertical)'
        
        if(DistMode)  then
          ADCSet_z=DeltaF_z       !if we were in DDC mode set the frequency shift setpoint to the actual one
          ADCo_z=ADC_z            !set the ADC offset to the actual value so that there wont be any jump in interactive mode
       endif
       
       ADCON_z=.false.                       !switch off the ADC
        
        exit
      endif
      
    endif

  end do
  
  !AGCON_z=.false.
	
	close(888)

end subroutine Approach



!*** ****** Phase 3b: Scanline (FOR INTERACTIVE MODE ONLY)
!   Input variables:
!   (DBL) LineLen         Length of the scanline
!   (INT) Points          Number of points to get on the line
!
!   Description:
!   this routine saves some scanlines for the interactive scanning mode.
!
subroutine ScanLine()
  use consts
  use timer
  use cantilever
  use PLL
  use AGC
  use ADC
  use output
  use Imager
  implicit none
  !integer, intent(in) :: Points
  !double precision, intent(in) :: LineLen
  integer :: Points
  double precision :: LineLen
  double precision :: pxDist,dist=0.0d0
  double precision, allocatable :: dataz(:),datay(:)
  integer :: cnt,pnt
  logical :: oldADCStatus
  
  oldADCStatus=ADCON_z
  
  LineLen = SL_Length
  Points = SL_npts
  
  allocate(dataz(0:Points),datay(0:Points))
  
  pxDist=LineLen/DBLE(Points)
  cnt=0
  pnt=0
  !t=0.0d0   !reset the time indicator

  call Freeevo(0.010d0,.false.) !let some time pass to adjust electronics
	phase=phase+1
	
  dataz(cnt)=DeltaF_z
  datay(cnt)=DeltaF_y
  call WriteScanLine()
  write(*,*)HolderPosition(:)
  SL_start(1:2) = HolderPosition(1:2) !save the starting position for recover

  do
    T=t+tstep     !increase the time indicator
        
    call Verlet(tstep)    !perform a Verlet step
    
    call UpdateCircuits() !operate the circuits
    
!     pnt=pnt+1
!     if((pnt==100)) then
!       !Call WriteAll()
!       !write(*,*)HolderPosition(2),At_z,DeltaF_z,Dissipation_z,FB_z
!       !write(*,*) DeltaF_z,At_z,ADC_z
!       !write(*,*)(FB_z-Aex0_z)/Aex0_z
!       pnt=0
!     endif

    !if the new image points was reached
    if(dist>=pxDist) then
      dist=0.0d0          !reset the travelled distance
      cnt=cnt+1
      
      call WriteScanLine()
      !write(*,*)HolderPosition(2),At_z,DeltaF_z,Dissipation_z,FB_z
      dataz(cnt)=DeltaF_z
      datay(cnt)=DeltaF_y
      
      if(cnt==Points) then     !*** if the line is done...
      
        ADCON_z=.false.
        call Lift(.true.,DBLE(5.0e-10))       !move the holder up of 5 angs
        call Reposition(.false.,.true.)       !move the holder back at the starting point
        call Lift(.false.,0.0d0)              !move the tip down to zLand
        call Freeevo(0.10d0,.false.)               !w8 some time
        ADCSet_z=DeltaF_z
        ADCON_z=oldADCStatus                  !reset the ADC as it was
        call Freeevo(0.10d0,.false.)               !w8 some time
	
        exit
      endif
      
    endif
    
    !move the tip
    HolderPosition(1:2) = HolderPosition(1:2) + ScanSpeed*tstep*SL_fdir(1:2)
    !HolderPosition(ScanAxis) = HolderPosition(ScanAxis) + ScanSpeed*tstep
    dist = dist + ScanSpeed*tstep
    
  end do

  !run some statistics on the scanline
  write(*,*)'Contrast Z:',Abs(MaxVal(dataz)-MinVal(dataz))
  write(*,*)'Contrast Y:',Abs(MaxVal(datay)-MinVal(datay))

end subroutine ScanLine


subroutine ZSpectrLine()
  use consts
  use timer
  use cantilever
  use PLL
  use AGC
  use ADC
  use output
  use Imager
  use Forcer
	implicit none
  double precision :: pxDist,dist=0.0d0
  integer :: cnt,pnt,znt
  logical :: oldADCStatus
  
  oldADCStatus=ADCON_z
  
  pxDist=SL_Length/DBLE(SL_npts)
  cnt=0
  pnt=0
  znt = 0
  !t=0.0d0   !reset the time indicator
  
  open(file="zSpec.out",unit=321)
	write(*,'(A)',advance='no')'Scanning: '
 SL_start(1:2) = HolderPosition(1:2) !save the starting position for recover
  do znt=0, ZS_npts
		!SCANLINE LOOP!*****************************************************
		call Freeevo(0.010d0,.false.) !let some time pass to adjust electronics
		!call WriteLateralMapXY(HolderPosition(3)-ASet_z,znt)
		write(321,*)HolderPosition(:),DeltaF_z,DeltaF_y
                
		cnt = 0
		do
			T=t+tstep     !increase the time indicator
			call Verlet(tstep)    !perform a Verlet step
			call UpdateCircuits() !operate the circuits
			!if the new image points was reached
			if(dist>=pxDist) then
				dist=0.0d0          !reset the travelled distance
				cnt=cnt+1
				write(321,*)HolderPosition(:),DeltaF_z,DeltaF_y
				if(cnt==SL_npts) then     !*** if the line is done...
					ADCON_z=.false.
					call Lift(.true.,DBLE(5.0e-10))       !move the holder up of 5 angs
					call Reposition(.false.,.true.)       !move the holder back at the starting point
					call Lift(.false.,0.0d0)              !move the tip down to zLand
					call Freeevo(0.10d0,.false.)               !w8 some time
					ADCSet_z=DeltaF_z
					ADCON_z=oldADCStatus                  !reset the ADC as it was
					call Freeevo(0.10d0,.false.)               !w8 some time
					exit
				endif
			endif
			!move the tip
                        HolderPosition(1:2) = HolderPosition(1:2) + ScanSpeed*tstep*SL_fdir(1:2)
			!HolderPosition(ScanAxis) = HolderPosition(ScanAxis) + ScanSpeed*tstep
			dist = dist + ScanSpeed*tstep
		end do
		!*****************************************************************
		call Lift(.true.,ZS_step)				!put the tip up a bit
		call Freeevo(0.10d0,.false.)	!wait a bit
		zLand = zLand + ZS_step
		dist = 0.0d0
		write(321,*)
		write(*,'(A)',advance='no')'.'
	enddo
	write(*,*)'done!'
	
end subroutine

!*** ****** Imaging Routine ****** ***!
!   Description:
!   this is the main routine to scan images.
!   It samples a scanline, and then drives the holder
!   to the next one.
!
subroutine ScanImage()
  use consts
  use Imager
  use cantilever
  use Timer
  implicit none
  integer :: i=0
  
  ScanQ(1) = 0 !reset the scanline counter

  do i=1, imgPixs      !loop on lines
     
     SL_start(1:2) = HolderPosition(1:2) !save the starting position for recover
     ScanQ(1) = ScanQ(1)+1  !increase the scanline counter
     
     write(*,*) 'Start line',HolderPosition(1),HolderPosition(2)
     Call ScanLoop(imgPix)                 !do the scanline sampling
     
     !reposition the tip
    
     call Reposition(.true.,.true.)
     
  enddo
  
  
  
end subroutine ScanImage


subroutine SinglePoint()
	use consts
	use Cantilever
	use Timer
	use PLL
	use AGC
	implicit none
	
	call FreeEvo(0.10d0,.false.)
	
	write(*,*) 'Single Point INFO (z,y):'
	write(*,*) 'Freqs', TotalF_z, TotalF_y
	write(*,*) 'DFs  ', DeltaF_z, DeltaF_y
	write(*,*) 'Amps ', At_z, At_y
	
	
end subroutine SinglePoint


!*** ****** Force curve calculator ****** ***!
!   Description:
!   The tip is moved towards zEnd with fixed size steps depending
!		on the given parameters, then the oscillation is relaxed and
!		data is stored to file.
!
subroutine ForceCurve_Static()
	use Cantilever
	use Output
	use Imager
	implicit none
	double precision :: FC_zEnd = dble(15.5e-9)
	integer :: FC_Points = 50
	double precision :: zStep
	integer :: i
	!logical :: up
	
	zStep = (FC_zEnd - HolderPosition(3)) / (FC_Points - 1)
	
	write(*,'(A)',advance='no')'Measuring df curve .'
	open(file='forcecurve.out',unit=555)
	
	call WriteForceCurve_Static(555)		!get the first point
	
	do i = 1, FC_Points - 1
		
		call DisplaceHolder(3, zStep)			!move the holder vertically
		call FreeEvo(0.10d0, .false.)			!wait for some time
		
		call WriteForceCurve_Static(555)	!write the data
		write(*,'(A)',advance='no')'.'
		
	enddo
	write(*,*)'done'
	close(555)
	
end subroutine ForceCurve_Static


!*** ****** Holder Repositioning ****** ***!
!   Input variables:
!   (LOG) down          .TRUE. if the tip has to slide down
!   (LOG) back          .TRUE. if the tip has to go backward
!   
subroutine Reposition(down,back)
  use cantilever
  use Timer
  use output
  use Imager
  implicit none
  logical, intent(in) :: down,back
  logical :: doneD=.false.,doneB=.false.
  double precision :: pxDist,xdist=0.0d0
  double precision :: pyDist,ydist=0.0d0
  double precision :: trlx=0.02d0   !time to wait after repositioning
  
  !write(*,'(A)',advance='no') 'Sliding to new line...'
  
  pxDist=imgSide/DBLE(imgPix)
  pyDist=imgSide/DBLE(imgPixs)
  
  doneD=.not.down
  doneB=.not.back
  xdist=0.0d0
  ydist=0.0d0
  trlx=0.02d0
  
  !loop to go back along the line
  if(back.and.(.not.doneB)) then
     
     do            !time evolution loop between scan points
        T=t+tstep
        call Verlet(tstep)    !perform a Verlet step
        call UpdateCircuits() !operate the circuits
        
        if(xdist<=pxDist*(imgPix-1)) then
           HolderPosition(1:2) = HolderPosition(1:2) - RepoSpeed*tstep*SL_fdir(1:2)
           !HolderPosition(ScanAxis) = HolderPosition(ScanAxis) - RepoSpeed*tstep
           xdist = xdist + RepoSpeed*tstep
        else
           doneB=.true.
           HolderPosition(1:2) = SL_start(1:2) !put the holder exactly at the beginning of the line
           exit
        endif
        
     end do
     
  end if
  
  if(down.and.(.not.doneD)) then
     
     do            !time evolution loop between scan points
        T=t+tstep
        call Verlet(tstep)    !perform a Verlet step
        call UpdateCircuits() !operate the circuits
        
        if(ydist <= pyDist) then !move the holder along the slowscan
           HolderPosition(1:2) = HolderPosition(1:2) + RepoSpeed*tstep*SL_sdir(1:2)
           !HolderPosition(ScanSlow) = HolderPosition(ScanSlow) + RepoSpeed*tstep
           ydist = ydist + RepoSpeed*tstep
        else
           doneD=.true.
           exit
        endif
       
    end do
    
    
 end if

 do            !time evolution loop between scan points
    T=t+tstep
    call Verlet(tstep)    !perform a Verlet step
    call UpdateCircuits() !operate the circuits
    
    !check if all itz done...
    if(doneD.and.doneB) then
       if(trlx>=0.0d0) then   !check if also some time has passed after repositioning
          trlx=trlx-tstep
       else
          exit
       endif
    endif
    
 enddo

  !write(*,*)'done.'


end subroutine Reposition


!*** ****** Tip lifter up/down ****** ***!
!	If up is .true. lift the tip up of height
! if up is .false. put the tip back to zLand!
subroutine Lift(up,height)
  use cantilever
  use timer
  use Imager
  use ADC
  implicit none
  logical, intent(in) :: up
  double precision, intent(in) :: height
  double precision :: dist
  
  dist=0.0d0
  ADCON_z=.false.
  
!   if(up) then
!     write(*,'(A)',advance='no')' Lifting tip...'
!   else
!     write(*,'(A)',advance='no')' ReApproaching tip...'
!   endif
  
  do            !time evolution loop between scan points
  
    T=t+tstep
    
    call Verlet(tstep)    !perform a Verlet step
  
    call UpdateCircuits() !operate the circuits
    
    
    if(up) then
      ADC_z=ADC_z+RepoSpeed*tstep
      HolderPosition(3)=ADC_z
      dist=dist+RepoSpeed*tstep   
      
      if(dist>=height) then     !check if it was lifted properly
! 				write(*,*) 'done!'
				exit
      endif
      
    else
      ADC_z=ADC_z-RepoSpeed*tstep
      HolderPosition(3)=ADC_z
      dist=dist+RepoSpeed*tstep
      
      ! in down lift mode, stop when the holder is at zLand
      if(ADC_z<=zLand) then
!         write(*,*)'done!'
        ADC_z=zLand
        HolderPosition(3)=ADC_z
        ADCo_z=ADC_z              !reset the ADC offset and integral part
        ADCi_z=0.0d0
        exit
      endif
      
    endif

    
  enddo
  
end subroutine Lift



!*** ****** Scan line sampler (from imager) ****** ***!
!   Input variables:
!   (INT) endPoint          the number of points to take on a scanline.
!
subroutine ScanLoop(endPoint)
  use Timer
  use imager
  use cantilever
  use output
  use PLL
  use AGC
  use NURBSDissiField
  implicit none
  integer, intent(in) :: endPoint
  double precision :: pxDist,dist=0.0d0
  double precision :: dataz(0:imgPix),datay(0:imgPix)
  integer :: cnt
  character*40 :: scanfile,linestr
  
  pxDist=imgSide/DBLE(imgPix)
  
  cnt=0
  
  write(linestr,*)ScanQ(1)
  scanfile='sl_'//trim(adjustl(linestr))//'.dat'
  open(file=scanfile, unit=888)
  
  call Freeevo(0.050d0,.false.)
  
  !sample the first point
  write(999,*) HolderPosition(1),HolderPosition(2),DeltaF_z,At_z,Dissipation_z,DissipationI_z, &
  DeltaF_y,At_y,Dissipation_y,DissipationI_y
  write(888,*) HolderPosition(1),HolderPosition(2),DeltaF_z,At_z,Dissipation_z,DissipationI_z, &
  DeltaF_y,At_y,Dissipation_y,DissipationI_y
  !write(*,*) HolderPosition(1),HolderPosition(2),DeltaF_z,At_z,Dissipation_z,DissipationI_z
  dataz(cnt)=DeltaF_z
  datay(cnt)=DeltaF_y
  
  do            !time evolution loop between scan points
    
    T=t+tstep
    
    call Verlet(tstep)    !perform a Verlet step
    call UpdateCircuits() !operate the circuits
    
    !if the new image points was reached
    if(dist>=pxDist) then
      dist=0.0d0          !reset the travelled distance
      cnt=cnt+1
      
      !write data for 3D plotting
      write(999,*) HolderPosition(1),HolderPosition(2),DeltaF_z,At_z,Dissipation_z,DissipationI_z, &
           DeltaF_y,At_y,Dissipation_y,DissipationI_y
      write(888,*) HolderPosition(1),HolderPosition(2),DeltaF_z,At_z,Dissipation_z,DissipationI_z, &
           DeltaF_y,At_y,Dissipation_y,DissipationI_y
      
      dataz(cnt)=DeltaF_z
      datay(cnt)=DeltaF_y
      
      if(cnt==endPoint) exit
      
    endif
    
    !move the tip
    HolderPosition(1) = HolderPosition(1) + ScanSpeed*tstep*SL_fdir(1)
    HolderPosition(2) = HolderPosition(2) + ScanSpeed*tstep*SL_fdir(2)
    !HolderPosition(ScanAxis) = HolderPosition(ScanAxis) + ScanSpeed*tstep
    dist = dist + ScanSpeed*tstep
    
  enddo
  
  write(999,*)      !write an empty line to separate the scanlines for 3D plotting
  close(888)
  

  !run some statistics on the scanline
  write(*,*)'Contrast Z:',Abs(MaxVal(dataz)-MinVal(dataz))
  write(*,*)'Contrast Y:',Abs(MaxVal(datay)-MinVal(datay))

end subroutine ScanLoop



subroutine DisplaceHolder(Axis,dx)
  use Cantilever
  use Imager
  use timer
  implicit none
  integer, intent(in) :: Axis
  double precision, intent(in) :: dx
  double precision :: dist
  double precision :: oldPos
  
  !t=0.0d0
  dist=0.0d0
  oldPos=HolderPosition(Axis)
  
  do
    T=t+tstep
    
    call Verlet(tstep)    !perform a Verlet step
  
    call UpdateCircuits() !operate the circuits
    
    !if the required distance has been travelled
    if(dist>=dx) then
      HolderPosition(Axis)=oldpos+dx  !place it in the right position (minimal correction)
      call FreeEvo(0.010d0,.false.)           !let 0.01 s pass just to be cool
      exit                            !get outta here
    endif
    
    !move the tip
    HolderPosition(Axis) = HolderPosition(Axis) + RepoSpeed*tstep
    dist = dist + RepoSpeed*tstep
    
  enddo
  

end subroutine DisplaceHolder



subroutine FreeEvo(timelen,writeDBG)
  use consts
  use cantilever
  use timer
  use Output
  use AGC
  use PLL
  implicit none
  double precision, intent(in) :: timelen
  logical, intent(in) :: writeDBG
  double precision :: tEnd
  integer :: cnt=0
  character*200 :: fname,tint
  
  tEnd=t+timelen
  cnt=0
  
  if(writeDBG) then
  	oFEVS = oFEVS + 1
  	write(tint,*) oFEVS
		fname = 'p_free_'//trim(adjustl(tint))//'.out'
		open(file=trim(adjustl(fname)), unit=888)
  endif
  
  do 
  
  	t=t+tstep
    cnt=cnt+1
    
    call Verlet(tstep)
    if(((Xo(3)>Xoo(3)).and.(Xo(3)>X(3))).and.(t>tEnd)) exit
    call UpdateCircuits()
  	!write(*,*)HolderPosition(2),At_z,DeltaF_z,HolderPosition(3)
  	
    if((cnt==PointBuffer).and.(writeDBG)) then
      Call WriteAll_f(888)
      !write(*,*)'AAASD'
      write(*,*)t,At_z,DeltaF_z,HolderPosition(3)
      cnt=0
    endif
    
    !stop only when the tip is at the top
    
  end do

	if(writeDBG) then
		close(888)
	endif

end subroutine FreeEvo




!*** Customized Verlet Velocity Algorithm ***
!   Input variables:
!   (DBL)dt         timestep
subroutine Verlet(dt)
  use consts
  use cantilever
  use timer
  implicit none
  double precision, intent(in) :: dt
  
  
  !*** Customized Verlet Velocity *********************************************************************
  X=X + V*dt + 0.50d0*A*dt**2                       !Calc X in the next time
  !V=V + 0.50d0*A*dt                                 !Calc V(T+Dt/2) USING A(T)
  V(3)= 0.50d0*A(3)*dt + V(3)*(1.0d0-0.50d0*dt*Gamma1)
  V(2)= 0.50d0*A(2)*dt + V(2)*(1.0d0-0.50d0*dt*Gamma2)

  !x(2)=0.0d0				!lock X,Y
  !x(1)=0.0d0!dble(2.820e-10)
  Call GetForce(A)                                  !Calc A(T+Dt)

  V(3)= 0.50d0*A(3)*dt + V(3)*(1.0d0-0.50d0*dt*Gamma1)
  V(2)= 0.50d0*A(2)*dt + V(2)*(1.0d0-0.50d0*dt*Gamma2)

  !V(3)=(V(3) + 0.50d0*A(3)*dt)/(1.0d0+0.50d0*dt*Gamma1) !Calc V(T+Dt)   USING A(T+Dt/2) (Z Component)
  !V(2)=(V(2) + 0.50d0*A(2)*dt)/(1.0d0+0.50d0*dt*Gamma2) !Calc V(T+Dt)   USING A(T+Dt/2) (Y Component)
  !V(1)=(V(1) + 0.50d0*A(1)*dt)                      !Calc V(T+Dt)   USING A(T+Dt/2) (X Component)
  !****************************************************************************************************
	
 	!X(3)=15.0e-9*Cos(dPi*t*F1)
	
end subroutine Verlet



!***************** ACCELERATION CALCULATION *****************!
!   Input variables
!   (DBL)Fout(3)    scrap vector to fill with Acceleration (no friction term included)
!
Subroutine GetForce(Fout)
  Use Consts
  use timer
  Use Cantilever
  Use AGC
  use PLL
  Implicit None
  Double Precision :: Fint(3)=0.0d0
  Double Precision,Intent(Out) :: Fout(3)

  Fout=0.0d0

  if(IntON_z.or.IntON_y) then
    Call Fext(Fint)
  endif

  !The Friction Terms Are Left Out for good, they will be added in the verlet algorithm!!!

  !vertcal part:    *******************************************
  Fout(3)=-W1**2*(X(3)) !cantilever spring term
  
  if(ExciteON_z) then   !if the external driver is on...
    if(ExciteSimple) then
      Fout(3) = Fout(3) + FB_z/mass*Sin((w1)*t+PiHalf)
    else
      Fout(3) = Fout(3) + FB_z/mass*SinPLL_z
    endif
  endif
  
  if(IntON_z) then      !if the tip-surface interaction is on...
    Fout(3) = Fout(3) + Fint(3)/mass
  endif
  
  !************************************************************
  
  !Lateral part (y): ******************************************
  Fout(2)=-W2**2*(X(2)) !cantilever torsional spring term
  
  if(ExciteON_y) then   !if the external driver is on...
    if(ExciteSimple) then
      Fout(2) = Fout(2) + FB_y*Sin(w2*t+PiHalf)/massTR!FB_y*SinPLL_y/massTR
    else
      Fout(2) = Fout(2) + FB_y*SinPLL_y/massTR
    endif
    
  endif
  
  if(IntON_y) then      !if the tip-surface interaction is on...
    Fout(2) = Fout(2) + Fint(2)/massTR
  endif

  !************************************************************  
  !Ftipsurf=Fout
  
End Subroutine GetForce
!***********************************************************************



!*** ****** Tip-Surface Interaction ****** ***!
!   Output variables:
!   (DBL)f(1:3)         scrap vector to fill with the tip-sample interaction
!   Description:
!   This calls Forcer routines to evaluate the force in the tip's position.
!
subroutine Fext(f)
  use consts
  use cantilever
  use timer
  use ADC
  use AGC
  use Forcer
  use NURBSDissiField
  implicit none
  double precision, intent(out) :: f(3)
  double precision :: pos(3)
  double precision :: vdW=0.0d0
  
  F=0.0d0
  
  pos(1)=HolderPosition(1)
  pos(2)=HolderPosition(2) + x(2)
  pos(3)=HolderPosition(3) + x(3)
  
  call CalcForce(pos,F)           !the forcer module knows what to do
  Ftipsurf=F
  call CalcVanDerWaals(pos,vdW)   !get the vdw interaction
  
  F(3)=F(3)+vdW                   !sum up short+long range interactions
  
  if(pos(3) < zCrashLimit)then
    write(*,*)'AAAAARGH! tip crash!'
    stop
  endif
  
  
end subroutine Fext



subroutine UpdateCircuits()
  use Cantilever
  use PLL
  use ADC
  use AGC
  use timer
  use forcer
  use output
  use NURBSDissiField
  implicit none
  double precision :: zSignal, ySignal,zPLLSignal,yPLLSignal,pos
  logical :: zTick=.false., yTick=.false.
  
  zSignal=(X(3)+Noise(3))
  ySignal=(X(2)+Noise(2))
  
  call FM_UpdateZ(tstep)
  call FM_UpdateY(tstep)

  zTick=(Xo(3)>Xoo(3)).and.(Xo(3)>X(3))
  if(zTick) call AGC_UpdateZ(Abs(zSignal),TotalF_z)
  
  yTick=(Xo(2)>Xoo(2)).and.(Xo(2)>X(2))
  if(yTick) call AGC_UpdateY(Abs(ySignal),TotalF_y)
  
  call ADC_UpdateZ(DeltaF_z,tstep)

  call Excite_Update()

!   if(zTick) then
!     !write(343,*)x(3),At_z,FB_z,At_y,FB_y
!   endif
	
  pos=HolderPosition(3)+X(3)
  
  !if entering the zone...
  if((.not.ATzone).and.((pos<ATup).and.(pos>ATdw)).and.(.not.IgnoreDissip)) then
    tstep=dtATzone						!change the time step
    call ForceNoiseInit(tstep)
    if(v(3)<0.0d0)PA_t=1.0d0	!if we are entering from above, set PA to 1
    ATzone=.true.							!activate the absolute terror switch
  elseif((ATzone).and.((pos>ATup)).and.(.not.IgnoreDissip)) then			!(pos<ATdw).or.
    !we are exiting the AT zone from the upper limit...
    tstep=dtNorm
    call ForceNoiseInit(tstep)
    ATzone=.false.
  endif
  
  !write(*,*)X(3)
!	call UpdateDissipation(X(3),Ftipsurf(3),TotalF_z,zTick,tstep)
	call DissiCalc(tstep,TotalF_z,TotalF_y,zTick,yTick)
	!call GetDissi(tstep)
	
! 	if(zTick) then
! 		if(phase>1) write(444,*)t,X(3),Noise(3)*J2eV
! 		Noise(3)=0.0d0
! 	endif
	
  !update the old positions
  Xoo(3)=Xo(3)
  Xo(3)=X(3)
  Xoo(2)=Xo(2)
  Xo(2)=X(2)
	ftipsurfo=ftipsurf
	ftipsurfoo=ftipsurfo
	
end subroutine UpdateCircuits


subroutine OscillDBG()
  use consts
  use cantilever
  use timer
  use output
  use PLL
  use AGC
  use ADC
  use NURBSDissiField
  use forcer
  implicit none
  integer :: cnt = 0,fn=300
  double precision :: tmpd = 0.0d0 
  
  t=0.0d0
  IntON_z=.true.
  IntON_y=.true.
  ADCON_z=.false.
  AGCON_z=.false.
  AGCON_y=.false.
  ExciteON_z=.false.
  ExciteON_y=.false.
  ExciteSimple=.false.
	!IgnoreDissip=.false.
  Gamma1=0.0d0
  Gamma2=0.0d0
  !initial position
  x(1)=0.0d0
  x(2)=0.0d0!dble(2.82e-10)!Aset_y
  x(3)=Aset_z
	Noise=0.0d0
  Ftipsurf=0.0d0
  !IgnoreDissip=.true.
  Xo=X
  Xoo=X
  
  v=0.0d0
  
  write(*,*)HolderPosition(:)
  HolderPosition(3)=zLand!Dble(15.440e-9)!Dble(15.50610e-9)
  tmpd = 0.0d0 
  Call GetForce(A)
  write(*,*)A
  phase=10
  
  DZold=Aset_z+HolderPosition(3)
  DYold=0.0d0
  DFold=A/mass
  do
    t=t+tstep

    call Verlet(tstep)
    X(2)=0.0d0
    X(1)=0.0d0
    
    if(X(3)+HolderPosition(3)<=8.0e-10) then
    	tmpd = tmpd+0.50d0*(Ftipsurf(3)+Ftipsurfo(3))*(x(3)-xo(3))
    	
    endif
    !           1         2               3     4        5   6   7   (8 9 10)
    write(151,*)t,X(3)+HolderPosition(3),PA_t,FmA(3),FmB(3),Wab,Wba,Ftipsurf(:)
    write(*,*)t,X(3)+HolderPosition(3),Ftipsurf(:)
    tmpd = tmpd+0.50d0*(Ftipsurf(3)+Ftipsurfo(3))*(x(3)-xo(3))
    
    if((Xo(3)>Xoo(3)).and.(Xo(3)>X(3))) then
    	cnt=cnt+1
    	if(cnt==1)	exit
    endif
    if((Xo(3)<Xoo(3)).and.(Xo(3)<X(3))) fn=fn+1
    
    call UpdateCircuits()
		!write(*,*) t,X(3)+HolderPosition(3),ATZone,ATup,Ftipsurf(3)
		
  enddo
!   write(*,*)'halting!',noise(3)*J2eV,YDissip*J2eV
  !call dbgNURBSz2()
  write(*,*)'Dissipated E',tmpd*J2eV
  
  
end subroutine OscillDBG


subroutine Commander()
  use timer
  use Imager
  use Cantilever
  use PLL
  use ADC
  use AGC
  use Output
  use Forcer
  implicit none
  character*80 :: command
  double precision :: dval
  integer :: ival
  do
    write(*,*)'what do u wanna do now?'
    read(5,'(A)') command
    command=trim(adjustl(command))
    
    if((command=='exit') .OR. (command=='gtfo'))then
      write(*,*)'u decided to go on with your life and leave me here! good bye then!'
      stop
    endif
  
    select case (command)
    case('help','wtf')
      call PrintHelp()
    case('status','sup')
      call PrintStatus()
      
    case('scanline','scan')           !get the scanline
      phase=phase+1
      write(*,'(A)')'ScanLine length:'
      read(*,*)SL_Length
      write(*,'(A)')'ScanLine points:'
      read(*,*)SL_npts
      call ScanLine()
    case('scanimage','image')
      phase=phase+1
      call ScanImage()
      
    case('movex')                     !move the holder along X
      write(*,*)'X displacement: '
      read(*,*)dval
      call DisplaceHolder(1,dval)
    
    case('movey')                     !move the holder along Y
      write(*,*)'Y displacement: '
      read(*,*)dval
      call DisplaceHolder(2,dval)
      
    case('zpiezoset')            !directly set the zpiezo
      phase=phase+1
      write(*,*)'Current zpiezo: ',ADC_z
      write(*,*)'type the desired zpiezo: '
      zStart=ADC_z
      read(*,*)zLand
      call Approach(.true.)
      
    case('adcset')         !set the adc frequency shift setpoint
      phase=phase+1
      write(*,*)'Current dF(z):',DeltaF_z
      write(*,*)'type the desired df(z): '
      read(*,*)ADCSet_z
      call Approach(.false.)
      zLand=HolderPosition(3)
      
    case('freeevo')
      write(*,*) 'Specify the time length: '
      read(*,*)dval
      call FreeEvo(dval,.true.)
    
    case('adctune')
      call ADC_Tune()
    case('agctuney')
    	call AGCy_Tune()
    
!     case('reload')
!       call LoadRelax()
!       Phase=1
!       if(tLand==0.0d0) then
!         call Approach(.false.)
!       else
!         call Approach(.true.)
!       endif
    case('temp')
      write(*,'(A,F8.3,A)')'Temperature is ', Temp,' K'
      write(*,'(A)')'new value:'
      read(*,*)Temp
      !call SetForceNoise()
      write(*,'(A,F8.3,A)')'Temperature is now ', Temp,' K'
    case('speed')
      write(*,'(A,ES16.8,A)')'The scan speed is ',ScanSpeed,' m/s'
      write(*,'(A)')'new value:'
      read(*,*)ScanSpeed
      write(*,'(A,ES16.8,A)')'The scan speed is now ',ScanSpeed,' m/s'
    
      
    end select
    
  enddo





end subroutine Commander

subroutine PrintHelp()
  implicit none
  

  write(*,*)'Avaiable commands: '
  write(*,*)'movex   - move the sample holder along the X axis'
  write(*,*)'movey   - move the sample holder along the Y axis'
  write(*,*)'adcset  - change the frequency shift and approach/relax the tip again'
  write(*,*)'scan    - run the holder and get a scanline of length 2La from the current position'
  write(*,*)'          and then go back to the starting point'
  write(*,*)'freeevo - let the cantilever go for the specified time'
  write(*,*)'reload  - reload the situation as it was after relaxing'
  write(*,*)'speed   - change the scanning speed'
  
end subroutine PrintHelp

subroutine PrintStatus()
  use consts
  use cantilever
  use ADC
  use AGC
  use PLL
  implicit none
  
  write(*,'(A, 3ES16.6)')'Holder position: ',HolderPosition
  write(*,'(A, ES16.6, ES16.6)')'Amplitude (z,y): ',At_z,At_y
  write(*,'(A, ES16.6, ES16.6)')'dFrequency (z,y): ',DeltaF_z,DeltaF_y
  
  write(*,'(A)',advance='no')'ADC circuit is '
  if(ADCON_z) then
    write(*,*)'ON'
  else
    write(*,*)'OFF'
  endif
  
  write(*,'(A)',advance='no')'AGC (z) circuit is '
  if(AGCON_z) then
    write(*,*)'ON'
  else
    write(*,*)'OFF'
  endif
  
  write(*,'(A)',advance='no')'AGC (y) circuit is '
  if(AGCON_y) then
    write(*,*)'ON'
  else
    write(*,*)'OFF'
  endif
  

end subroutine PrintStatus


end program Torfort
