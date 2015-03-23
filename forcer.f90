!*** ****** Forcer Module ****** ***!
!
!   Description:
!   It manages the tip-sample interaction forces. The starting point
!   is ReadGridFile() which reads the forcefield grid parameters and then builds
!   a 3D grid in the gridQ(:) arrays. There are also interpolation routines which are
!   good to get the forces between the specified points. The tip doesn't spend all
!   day jumping from grid point to grid point u know!?! Then you can read the forces
!   to the main matrix ForceField using the ReadForceField() routine. The forces are supposed to
!   die whenever you jump too far, i.e. above the highest grid point (zmax)!
!
module Forcer
implicit none
integer, parameter :: rndSize = 100
integer :: ForceType											!the type of force file: 1-> Map, 2->SinglePoint

character*40 :: GridFile                  !the grid parameters filename
integer :: pX,pY,pZ,pdZ                   !number of grid points along each axis + zdissimap
double precision :: xMin,xMax,xStep       !minimum and maximum X value
double precision, allocatable :: gridX(:) !the grid points coordinates along X
double precision :: yMin,yMax,yStep       !minimum and maximum Y value
double precision, allocatable :: gridY(:) !the grid points coordinates along Y
double precision :: zMin,zMax,zStep       !minimum and maximum Z value
double precision, allocatable :: gridZ(:) !the grid points coordinates along Z
double precision :: hx,hy,hz              !grid steps size
double precision :: dzMin,dzMax,dzStep    !grid range for dissipation map
double precision :: ATup,ATdw             !Absolute Terror zone heights

integer :: xsymmetry = 1                  !symmetry along x axis (1->pacman , 2->mirror)
integer :: ysymmetry = 1                  !symmetry along y axis (1->pacman , 2->mirror)
integer :: Geometry = 1                   !Surface geometry (1->Uniform , 2->Step)
double precision :: StepYPoint = 0.0d0    !the Y coord where the step is (last point on the topmost terrace)
double precision :: StepHeight = 0.0d0    !the step heigh
double precision :: StepBroad  = 0.0d0    !the step broadening (for an der Waals morphing)
double precision :: GridFactor            !the unit of the grid

character*40 :: ForceFile                               !the force field filename
character*40 :: DissiFile                               !the force field filename
logical :: IgnoreDissip = .false.
logical :: NoDissiEver  = .false.
double precision, allocatable :: ForceField(:,:,:,:)    !the force field values (xindex,yindex,zindex,component)
double precision :: Temp
double precision :: ZDissip = 0.0d0, DZold = 0.0d0, DFold(3)
double precision :: YDissip = 0.0d0,YDissipStore = 0.0d0, DYold = 0.0d0
integer :: YDissipN = 0

logical :: ForceNoise = .true.

double precision :: rndNums(2*rndSize),SDF_z,SDF_y
integer :: rndIndex = 1

contains


subroutine ForceInit()
  use NURBSForceField
  use timer
  implicit none
  
  
  write(*,*)
  call ReadGridFile()
  call GenerateRegularGrid()
  
  
  call ForceSetup()             !setup the forcefield for interpolation
  
  if(.not.NoDissiEver) call DissiSetup()
  
  call ForceNoiseInit(tstep)
  
  !call dbgdbg()

  !stop
end subroutine ForceInit

subroutine ForceNoiseInit(dt)
	use cantilever
	use consts
	implicit none
	double precision, intent(in) :: dt
	
	
	!Temp = 1000.0d0
	SDF_z=dsqrt((4.d0*kB*Temp*K_z/(Q1*w1))/dt) ! /dt?!
	SDF_y=dsqrt((4.d0*kB*Temp*K_y/(Q2*w2))/dt)
	!call RMARIN(123,456)
	!refill the random vector
	!call gauss(rndNums,rndSize)
	call FillRndVector()
	!write(*,*)rndSize
	!write(*,*)rndNums(:)
	!stop
end subroutine ForceNoiseInit

!*** setup the NURBSForce 
subroutine ForceSetup() 
  use NURBSForceField
  implicit none
  integer :: i,j,k,c
  double precision,allocatable :: tmparr(:)
  
  call NURBSfAllocate(pX+3,pY+3,pZ)
  
  allocate(tmparr(0:pX+3))
  if( pX >= 3) then
		tmparr(2:pX+1)=gridX(1:pX)
		tmparr(0) = -gridX(3)
		tmparr(1) = -gridX(2)
		tmparr(pX+2)=xMax
		tmparr(pX+3)=gridX(2)+xMax
	else
		do i = 0, pX+3
			tmparr(i) = i * dble(1.0e-10)
		enddo
	endif
	call NURBSfSetGrid(tmparr,pX+3,1)
	deallocate(tmparr)


  allocate(tmparr(0:pY+3))
  if( pY >= 3) then
		tmparr(2:pY+1)=gridY(1:pY)
		tmparr(0) = -GridY(3)
		tmparr(1) = -GridY(2)
		tmparr(pY+2)=yMax
		tmparr(pY+3)=gridY(2)+yMax
	else
		do i = 0, pY+3
			tmparr(i) = i * dble(1.0e-10)
		enddo
	endif
  call NURBSfSetGrid(tmparr,pY+3,2)
	deallocate(tmparr)
	
		
  allocate(tmparr(0:pZ))
  tmparr(0:pZ-1)=gridZ(1:pZ)
  tmparr(pZ)=zMax
  call NURBSfSetGrid(tmparr,pZ,3)
  deallocate(tmparr)
  
  ForceType = 1
  write(*,*)'Reading forcefield file...',ForceType

  if(ForceType == 1) call ReadForceFieldSimple(ForceFile)   !read the force field values
  if(ForceType == 2) call ReadForceFieldSingle(ForceFile)   !read the force field values
  
  !setup the NURBS interpolation
  call NURBSfSetup()
  !call dbgNURBSz()
  
end subroutine ForceSetup

!*** setup the Disspiation map
! This is called only if the DissiFile is given in the input.
subroutine DissiSetup()

  use NURBSDissiField
  implicit none
  integer :: i
  double precision, allocatable :: tmparr(:)
  double precision :: PA,ForceA(3),ForceB(3),x(3)
  
  
  call NURBSdAllocate(pX+3,pY+3,pdZ)      !allocate the dissipation data stuff
  
  
  allocate(tmparr(0:pX+3))
  if( pX >= 3) then
  tmparr(2:pX+1)=gridX(1:pX)
  tmparr(0) = -gridX(3)
  tmparr(1) = -gridX(2)
  tmparr(pX+2)=xMax
  tmparr(pX+3)=gridX(2)+xMax
  else
		do i = 0, pX+3
			tmparr(i) = i * dble(1.0e-10)
		enddo
	endif
	call NURBSdSetGrid(tmparr,pX+3,1)
  deallocate(tmparr)
  
  allocate(tmparr(0:pY+3))
  if( pY >= 3) then
		tmparr(2:pY+1)=gridY(1:pY)
		tmparr(0) = -GridY(3)
		tmparr(1) = -GridY(2)
		tmparr(pY+2)=yMax
		tmparr(pY+3)=gridY(2)+yMax
  else
		do i = 0, pY+3
			tmparr(i) = i * dble(1.0e-10)
		enddo
	endif
  call NURBSdSetGrid(tmparr,pY+3,2)
  deallocate(tmparr)
  

  allocate(tmparr(0:pdZ))
  do i=0,pdZ
    tmparr(i)=dzMin + dzStep*i !(dzMax-dzMin)*dble(i)/dble(pdZ-1)
  enddo
  call NURBSdSetGrid(tmparr,pdZ,3)
  deallocate(tmparr)
  write(*,*)'Reading dissipation map... ',DissiFile
  
  if(ForceType == 1) call ReadDissimapSimple(DissiFile)   !read the force field values
  if(ForceType == 2) call ReadDissimapSingle(DissiFile)   !read the force field values
	
  call NURBSdSetup(Temp)
!   call dbgAll()
!   stop
  
end subroutine DissiSetup

subroutine CalcForce(x,f)
  use NURBSForceField
  use NURBSDissiField
  use timer
  implicit none
  double precision, intent(in) :: x(3)
  double precision, intent(out) :: f(3)
  double precision :: fy,fz,fr(3),xc(3)
  double precision :: PA
  
  f=0.0d0
  if(x(3)>=GridZ(pZ)) then
  	call AddNoise(F)
  	return
  endif
  
  call CenterCursor(x,xc)
  
  !check if we have to use the non-conservative method
  if(ATzone) then
  	!write(*,*)'using dissi'
	 	call Evolve(tstep,xc,PA,fr)					!evolve the probabilities PA/PB
  	if(PA<0.0d0) then
  		call NURBSForce(xc,fr)   !use the conservative field
  	endif
  else
!   	write(*,*)'CONS'
		!ZDissip = 0.0d0
  	call NURBSForce(xc,fr)			!use the conservative field
  	!write(555,*)x(3),1
  endif
	
	f=0.0d0
	f(:)=fr(:)
	!f(2)=0.0d0
	f(1)=0.0d0
! 	call UpdateDissipationZ(xc(3),fr(3))
! 	call UpdateDissipationY(xc(2),fr(2))
	
	!add the force noise
	call AddNoise(F)
	
	!write(666,*)x(3),f(3)
	DZold=xc(3)
	DYold=xc(2)
	DFold=f
	
end subroutine CalcForce

subroutine CalcConsForce(x,f)
  use NURBSForceField
  use NURBSDissiField
  use timer
  implicit none
  double precision, intent(in) :: x(3)
  double precision, intent(out) :: f(3)
  double precision :: fy,fz,fr(3),xc(3)
  double precision :: PA
  
  
  if(x(3)>=GridZ(pZ)) return
  
  call CenterCursor(x,xc)
	call NURBSForce(xc,fr)			!use the conservative field
	f=0.0d0
	f(:)=fr(:)
	f(1)=0.0d0
! 	f(2)=f(2)+dble(rndNums(rndIndex))*SDF_y
!   f(3)=f(3)+dble(rndNums(rndIndex+rndSize))*SDF_z
! 	rndIndex = rndIndex + 1
! 	if(rndIndex == rndSize) then
! 		call FillRndVector()	!gauss(rndNums,rndSize)
! 		rndIndex = 1
! 	endif
	
end subroutine CalcConsForce

subroutine AddNoise(f)
  implicit none
  double precision, intent(inout) :: f(3)

	!call AddNoise(F)
	if(ForceNoise) then
		f(2)=f(2)+dble(rndNums(rndIndex))*SDF_y
		f(3)=f(3)+dble(rndNums(rndIndex+rndSize))*SDF_z
		rndIndex = rndIndex + 1
		if(rndIndex == rndSize) then
			call FillRndVector()	!gauss(rndNums,rndSize)
			rndIndex = 1
		endif
	endif

end subroutine AddNoise
subroutine FillRndVector()
  use ziggurat
  implicit none
  integer :: i

  do i=1,rndSize*2
    rndNums(i)=rnor()
  enddo

  return

end subroutine FillRndVector


subroutine UpdateDissipationZ(Z,ZForce)
	implicit none
	double precision, intent(in) :: Z, ZForce
	
	!ZDissip=ZDissip+ZForce*(Z-DZold)
	ZDissip= ZDissip + (ZForce+DFold(3))*(Z-DZold)/2.0d0
	
end subroutine UpdateDissipationZ

subroutine UpdateDissipationY(Y,YForce)
	use cantilever
	implicit none
	double precision, intent(in) :: Y, YForce
	logical :: yTick
	
	!update the integral
	YDissip = YDissip + (YForce+DFold(2))*(Y-DYold)/2.0d0
	
	!check if theres a ytick and ...
! 	yTick=(Xo(2)>Xoo(2)).and.(Xo(2)>X(2))
! 	if(ytick) then
! 		YDissipStore=YDissipStore+YDissip			!increment the total dissipation
! 		YDissip=0.0d0													!reset the dissipation on one cycle
! 		YDissipN=YDissipN+1										!increment the cycle counter
! 		
! 	endif
	
end subroutine UpdateDissipationY



subroutine CalcVanDerWaals(p,vdw)
  use cantilever
  implicit none
  double precision, intent(in) :: p(3)
  double precision, intent(out) :: vdw
  double precision :: vdWup=0.0d0,vdWdw=0.0d0
  double precision :: pd(3)
  
  
  if(Geometry==1) then        !if itz plain surface
    
    !*** van der Waals contribution (already in N) ***
    vdW=(TipHamaker*TipRadius**2)*(1.0d0-sing)*(TipRadius*sing-p(3)*sing-TipRadius-p(3))
    vdW=vdW/(6.0d0*(p(3)**2)*(TipRadius+p(3)-TipRadius*sing)**2)
    vdW=vdW-(TipHamaker*tang*(p(3)*sing+TipRadius*sing+TipRadius*cos2g))/(6.0d0*cosg*(TipRadius+p(3)-TipRadius*sing)**2)
    return
    
  elseif(Geometry==2) then    !if there is a step
    
    call CenterCursor(p,pd)   !center the cursor in the unit cell
    pd(3)=p(3)+StepHeight
    
    !vdW interaction on the top terrace
    vdWup=(TipHamaker*TipRadius**2)*(1.0d0-sing)*(TipRadius*sing-p(3)*sing-TipRadius-p(3))
    vdWup=vdWup/(6.0d0*(p(3)**2)*(TipRadius+p(3)-TipRadius*sing)**2)
    vdWup=vdWup-(TipHamaker*tang*(p(3)*sing+TipRadius*sing+TipRadius*cos2g))/(6.0d0*cosg*(TipRadius+p(3)-TipRadius*sing)**2)
    
    !vdW interaction on the lower terrace
    vdWdw=(TipHamaker*TipRadius**2)*(1.0d0-sing)*(TipRadius*sing-pd(3)*sing-TipRadius-p(3))
    vdWdw=vdWdw/(6.0d0*(pd(3)**2)*(TipRadius+pd(3)-TipRadius*sing)**2)
    vdWdw=vdWdw-(TipHamaker*tang*(pd(3)*sing+TipRadius*sing+TipRadius*cos2g))/(6.0d0*cosg*(TipRadius+pd(3)-TipRadius*sing)**2)
    
    vdw=vdWup/(1.0d0+exp((pd(2)-StepYPoint)/StepBroad))+vdWdw/(1.0d0+exp(-(pd(2)-StepYPoint)/StepBroad))
    
    return
    
  endif
  

end subroutine CalcVanDerWaals




!*** ****** Find the containing grid box lower indexes *** ******!
!   Input variables:
!   (DBL)x(3)         the tip position, already reshifted in the unit cell
!
!   Output variables:
!   (INT)idx(3)       the lower indexes of the grid volume containing x
!
!   Description:
!   Killroy was here!
!
subroutine FindIndex(x,idx)
  implicit none
  double precision, intent(in) :: x(3)
  integer, intent(out) :: idx(3)
  integer :: i
  
  idx=0
  
  !*** x axis
  if( (x(1)>=gridX(pX)) ) then
    idx(1)=pX
  else
    do i=1,pX
      if(gridX(i)>x(1)) then
        idx(1)=i-1
        exit
      endif
    enddo
  endif
  

  
  !*** ******
    
  !*** y axis
  if( x(2)>=gridY(pY) ) then
    idx(2)=pY
  else
    do i=1,pY
      if(gridY(i)>x(2)) then
        idx(2)=i-1
        exit
      endif
    enddo 
  endif

  !*** ******

  !*** z axis
  if(x(3)<=zmin) then
    idx(3)=1
    return
  endif
  
  if( x(3)>=gridZ(pZ) ) then
    idx(3)=pZ
  else
    do i=1,pZ
      !write(*,*)'comparing',x(3),'to',gridZ(i), (gridZ(i)>x(3))
      if(gridZ(i)>x(3)) then
        idx(3)=i-1
        exit
      endif
    enddo 
  endif
  
  

end subroutine FindIndex


subroutine CenterCursor(x,xc)
  implicit none
  double precision, intent(in) :: x(3)
  double precision, intent(out) :: xc(3)
  
  
  !center the cursor in the x direction
  xc(1)=x(1)
  if(xsymmetry==1)then
    !the pacman symmetry
    do
      if((xc(1)>xMax).or.(xc(1)<0.0d0)) then 
				xc(1)=xc(1)-(x(1)/Abs(x(1)))*xMax
      else
				exit
      endif
    enddo
    
!   else
!     !mirror geometry
!     
!     if(xc(1)<0) xc(1)=-xc(1)      !if the point is in the negative side, flip it!
!     
!     do
!       if((xc(1)>=2.0d0*GridX(pX))) then
!         xc(1)=xc(1)-(2.0d0*GridX(pX))      !this centers the point in the range [0,2L]
!       else
!         if( xc(1)>GridX(pX) ) then
!           xc(1)=2.0d0*GridX(pX)-xc(1)
! 				endif
!         exit
!       endif
!     enddo
  endif
  
  !center the cursor in the y direction
  xc(2)=x(2)
  if(ysymmetry==1)then    !use the right symmetry fcs
    !the pacman symmetry
    do
      if((xc(2)>yMax).or.(xc(2)<0.0d0)) then 
	xc(2)=xc(2)-(x(2)/abs(x(2)))*yMax
      else
	exit
      endif
    enddo
    
  else
    !the mirror symmetry
    
    if(xc(2)<0) xc(2)=-xc(2)      !if the point is in the negative side, flip it!
    
    do
      if((xc(2)>=2.0d0*GridY(pY))) then
        xc(2)=xc(2)-(2.0d0*GridY(pY))      !this centers the point in the range [0,2L]
      else
        if( xc(2)>GridY(pY) ) then
          xc(2)=2.0d0*GridY(pY)-xc(2)
	endif
        exit
      endif
    enddo
    
  endif
  
  xc(3)=x(3)
  
  return

end subroutine CenterCursor




!*** ****** reads the gridfile ****** ***!
subroutine ReadGridFile()
  implicit none
  integer :: idx,ifile=999
  character*256 :: line,tmpstr,strv
  
  GridFile=trim(ADJUSTL(GridFile))
  
  open(file=trim(ADJUSTL(GridFile)), unit=ifile)
  write(*,*)'reading grid file: ' // GridFile
  
  rewind(ifile)
  call find_str(ifile,'gridunits',idx,line)
  if(line/='000') then
    read(line,'(A9 A20)')tmpstr,strv
    strv=trim(adjustl(strv))
    if(strv(1:3) == 'ang') then
      GridFactor=dble(1.0e-10)
      write(*,*)'Grid parameters are in Angstroms.' 
    endif
      
    if(strv=='uni') GridFactor=1.0d0
  endif
  
  !*** read the geometric parameters **************
  rewind(ifile)  
  call find_str(ifile,'geometry',idx,line)
  if(line/='000') then
    read(line,'(A8 A)')tmpstr,strv
    if(trim(adjustl(strv))=='step') then
      write(*,'(A)',advance='no')'Surface with a step at point: '
      rewind(ifile)  
      call find_str(ifile,'stepYpoint',idx,line)
      if(line/='000') then
        read(line,'(A10 F14.8)')tmpstr,StepYPoint
        write(*,*)StepYPoint*GridFactor
        Geometry=2
      else
        write(*,'(A)')'FATAL! Unable to find the step Y gridpoint! Check your input!'
        stop
      endif
      write(*,'(A)',advance='no')'Step height: '
      rewind(ifile)  
      call find_str(ifile,'stepheight',idx,line)
      if(line/='000') then
        read(line,'(A10 F14.8)')tmpstr,StepHeight
        write(*,*)StepHeight*GridFactor
      else
        write(*,'(A)')' FATAL! Unable to find the step height! Check your input!'
        stop
      endif  
      write(*,'(A)',advance='no')'Step broadening: '
      rewind(ifile)  
      call find_str(ifile,'stepbroad',idx,line)
      if(line/='000') then
        read(line,'(A9 F14.8)')tmpstr,StepBroad
        write(*,*)StepBroad*GridFactor
      else
        write(*,'(A)')'FATAL! Unable to find the step width! Check your input!'
        stop
      endif
      
    else if(trim(adjustl(strv))=='unity') then
      write(*,*)'Plain surface geometry.'
      Geometry=1
    else
      write(*,*)'FATAL! Unrecognized surface geometry. Use step or unit!'
      stop
    endif
    
  endif
  !************************************************
  
  
  !*** read the grid size *************************
  rewind(ifile)  
  call find_str(ifile,'x.gridpoints',idx,line)
  if(line/='000') then
    read(line,'(A12 I8)')tmpstr,pX
  endif  
  rewind(ifile)  
  call find_str(ifile,'y.gridpoints',idx,line)
  if(line/='000') then
    read(line,'(A12 I8)')tmpstr,pY
  endif  
  rewind(ifile)  
  call find_str(ifile,'z.gridpoints',idx,line)
  if(line/='000') then
    read(line,'(A12 I8)')tmpstr,pZ
  endif
  call find_str(ifile,'d.zgridpoints',idx,line)
  if(line/='000') then
    read(line,'(A13, I8)')tmpstr,pdZ
  endif
  !************************************************
  
  !*** read the grid min vals *********************
  rewind(ifile)  
  call find_str(ifile,'x.gridMin',idx,line)
  if(line/='000') then
    read(line,'(A9 F14.8)')tmpstr,xMin
  endif  
  rewind(ifile)  
  call find_str(ifile,'y.gridMin',idx,line)
  if(line/='000') then
    read(line,'(A9 F14.8)')tmpstr,yMin
  endif  
  rewind(ifile)  
  call find_str(ifile,'z.gridMin',idx,line)
  if(line/='000') then
    read(line,'(A9 F14.8)')tmpstr,zMin
  endif
  rewind(ifile)
  call find_str(ifile,'dmap.zMin',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,dzMin
  endif
  !************************************************
  
  !*** read the grid Max vals *********************
  rewind(ifile)  
  call find_str(ifile,'x.gridMax',idx,line)
  if(line/='000') then
    read(line,'(A9 F14.8)')tmpstr,xMax
  endif  
  rewind(ifile)  
  call find_str(ifile,'y.gridMax',idx,line)
  if(line/='000') then
    read(line,'(A9 F14.8)')tmpstr,yMax
  endif  
  rewind(ifile)
  call find_str(ifile,'z.gridMax',idx,line)
  if(line/='000') then
    read(line,'(A9 F14.8)')tmpstr,zMax
  endif
  rewind(ifile)
  call find_str(ifile,'dmap.zMax',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,dzMax
  endif
  !************************************************

  rewind(ifile)
  call find_str(ifile,'dmap.step',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,dzStep
  endif
  rewind(ifile)
  call find_str(ifile,'x.fgstep',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,xStep
  endif
  rewind(ifile)
  call find_str(ifile,'y.fgstep',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,yStep
  endif
  rewind(ifile)
  call find_str(ifile,'z.fgstep',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,zStep
  endif

  !*** read the grid symmetry *********************
  rewind(ifile)  
  call find_str(ifile,'x.symmetry',idx,line)
  if(line/='000') then
    read(line,'(A10 A)')tmpstr,strv
    if(trim(adjustl(strv))=='mirror')then
      xsymmetry=2
    else
      xsymmetry=1
    endif
  endif  
  rewind(ifile)  
  call find_str(ifile,'y.symmetry',idx,line)
  if(line/='000') then
    read(line,'(A10 A)')tmpstr,strv
    if(trim(adjustl(strv))=='mirror')then
      ysymmetry=2
    else
      ysymmetry=1
    endif
  endif  
  !************************************************
  
  rewind(ifile)
  call find_str(ifile,'dmap.ATup',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,ATup
  endif
  rewind(ifile)
  call find_str(ifile,'dmap.ATdw',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,ATdw
  endif
  
  !*** Multiply everything by GridFactor **********
  xMin=xMin*GridFactor
  yMin=yMin*GridFactor
  zMin=zMin*GridFactor
  dzMin=dzMin*GridFactor
  
  xMax=xMax*GridFactor
  yMax=yMax*GridFactor
  zMax=zMax*GridFactor
  dzMax=dzMax*GridFactor
  
  ATup=ATup*GridFactor
  ATdw=ATdw*GridFactor
  
  dzStep = dzStep*GridFactor
  xStep = xStep*GridFactor
  yStep = yStep*GridFactor
  zStep = zStep*GridFactor
  
  StepYPoint=StepYPoint*GridFactor
  StepHeight=StepHeight*GridFactor
  StepBroad=StepBroad*GridFactor
  !************************************************
    
  close(ifile)
  
end subroutine ReadGridFile


!*** generate the grid points from the input parameters
subroutine GenerateRegularGrid()
  implicit none
  integer :: cnt
  double precision :: asd
  
  write(*,*)'Generating grid points...',pX,pY,pZ
  
  allocate(gridX(pX))
  allocate(gridY(pY))
  allocate(gridZ(pZ))
  
  hx=xStep!(xMax-xMin)/dble(pX)
  hy=yStep!(yMax-yMin)/dble(pY)
  hz=zStep!(zMax-zMin)/dble(pZ)
  
  do cnt=1,pX
    gridX(cnt)=xMin+(cnt-1)*xStep
  enddo
  do cnt=1,pY
    gridY(cnt)=yMin+(cnt-1)*yStep
  enddo
  do cnt=1,pZ
    gridZ(cnt)=zMin+(cnt-1)*zStep  !hz
  enddo 
  
  open(file="grid.out",unit=131)
  write(131,*)'Generated Grid'
  write(131,*)'xPoints:'
  write(131,*)gridX
  write(131,*)'yPoints:'
  write(131,*)gridY
  write(131,*)'zPoints:'
  write(131,*)gridZ 
  close(131)
  
  write(*,*)'Grid generated (grid.out)'
  if( (pX == 1) .and. (pY == 1) ) then
  	write(*,'(A)')'INFO: the grid has only one point in the XY plane!'
  endif
  
end subroutine GenerateRegularGrid








!*** Lateral Force Field Map ***!
!	Prints out the lateral force field at some height
!
subroutine writeLateralMapXY(height,indexer)
  use consts
  implicit none
  integer, intent(in) :: indexer
  double precision, intent(in) :: height
  double precision :: x,y,pos(3),f(3),f2(3),grd,pos2(3),pos3(3),f3(3)
  integer :: i,j
  double precision :: ystep = dble(1.0e-12)
  character*80 :: fidx,fname

  !write the 2D force map (X,Y) at fixed height
  write(fidx,*)indexer
  fname = "map_lateral_xy_"//trim(adjustl(fidx))//".dat"
  open(file=fname,unit=772)
  pos=0.0d0
  pos(3)=height
  do x=0.0d0, 2.0*xMax, xMax/16.0d0
    do y=0.0d0, 2.0*yMax, yMax/16.0d0
      pos(1)=x
      pos(2)=y
      pos2 = pos
      pos2(2) = pos(2) + ystep
      call CalcForce(pos,f)
      call CalcForce(pos2,f2)
      write(772,*)x,y,f(2),(f2(2)-f(2))/ystep
    enddo
    write(772,*)
  enddo
  close(772)
  
end subroutine writeLateralMapXY

!*** Lateral Force Field Map ***!
!	Prints out the lateral force field at some height
!
subroutine writeLateralMapYZ(xPosition)
  use consts
  implicit none
  double precision, intent(in) :: xPosition
  double precision :: z,y,pos(3),f(3),f2(3),pos2(3)
  double precision :: ystep = dble(1.0e-10), zstep = dble(0.05e-10)

  !write the 2D force map (X,Y) at fixed height
  open(file="map_lateral_yz.out",unit=772)
  pos=0.0d0
  pos(1)=xPosition
  do z=GridZ(30), zMin+dble(8.0e-10), zstep
  	pos(3)=z
    do y=0.0d0, 2.0*yMax, yMax/64.0d0
      pos(2)=y
      pos2 = pos
      pos2(2) = pos(2) + ystep
      call CalcConsForce(pos,f)
      call CalcConsForce(pos2,f2)
      write(772,*)y,z,f(2),(f2(2)-f(2))/ystep
    enddo
    write(772,*)
  enddo
  close(772)
  
  open(file="map_vertical_yz.out",unit=772)
  pos=0.0d0
  pos(1)=xPosition
  do z=GridZ(30), zMin+dble(8.0e-10), zstep
  	pos(3)=z
    do y=0.0d0, 2.0*yMax, yMax/64.0d0
      pos(2)=y
      pos2 = pos
      pos2(3) = pos(3) + ystep
      call CalcConsForce(pos,f)
      call CalcConsForce(pos2,f2)
      write(772,*)y,z,f(3),(f2(3)-f(3))/ystep
    enddo
    write(772,*)
  enddo
  close(772)
  
  
end subroutine writeLateralMapYZ

subroutine writeForceLine(x,z)
	use NURBSForceField
	implicit none
	double precision, intent(in) :: x,z
	double precision :: y, pos(3),f(3)
	integer :: i
	
	pos = 0.0d0
	pos(1) = x
	pos(3) = z
	
	open(file='forceline.out',unit=1234)
	do y=-yMax, 2.0d0*yMax, yMax/64.0d0
		pos(2) = y
		call CalcConsForce(pos,f)
		write(1234,*)y,f(:)
	enddo
	close(1234)
	open(file='forcedata.out',unit=111)
	do i=0,pY+3
		pos(2) = Gi(i)
		call CalcConsForce(pos,f)
		write(111,*)pos(2),Dijk(2,i,50,:)
		!write(111,*)pos(2)+ymax,f(:)
	enddo
	close(111)
	
	!write(*,*)Dijk(0,:,55,3)
	
end subroutine writeForceLine



!*** this is the same routine found in input.f90
subroutine find_str(fnum,str,idx,line)
  implicit none
  integer, intent(in) :: fnum         !file number
  character*8, intent(in) :: str      !string to find (itz always 8 characters!!!
  integer, intent(out) :: idx         !character index in 'line'
  character*256, intent(out):: line   !Line containing 'str'
  
  
  10  read(fnum, '(A)', end=20) line                !read a line
      !write(*,*)line
      idx = index (line, str)                       !look for the string
  		!write(*,*)line,'||',str
      if ((idx .ne. 0) .and. (line(1:1)/='#')) then !if the string is present...
	return                                      !exit the subroutine
      endif                                         !***
      goto 10                             !loop!
      
  20  line='000'                    !the end of the file was reached!
      return


end subroutine find_str





end module Forcer


