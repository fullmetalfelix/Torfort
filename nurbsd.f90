module NURBSDissiField
  
  !Private
  !*** this module contains all the stuff to get the non conservative force in (x,y,z)

!character*200, public :: filein

integer, public :: Order=1                                        !NURBS order
double precision, allocatable :: Gi(:),Gj(:),Gk(:)        !gridpoint space coords
double precision, allocatable :: Dijk(:,:,:,:)            !datapoints on the grid (i,j,k,c) c=ETA,VA,fA(1:3),ETB,VB,fB(1:3),ET
integer :: NPts(3),MPts(3),FDim=12

double precision, allocatable :: Pijk(:,:,:,:)            !NURBS control points
double precision, allocatable :: Pi(:),Pj(:),Pk(:)        !Control points space coords

double precision, allocatable :: usk(:),vsk(:),wsk(:)     !NURBS parameters vectors
double precision, allocatable :: uk(:),vk(:),wk(:)        !NURBS knots vectors
double precision, allocatable :: Nb(:)                    !NURBS basis functions
double precision, allocatable :: Nb3(:,:)

double precision,Public :: PA_t=1.0d0,PAI_t,PAIo_t=0.0d0,Wbao_t=0.0d0      !probability of state A at time t
double precision,Public :: sumw_t=0.0d0, sumwo_t=0.0d0, exp_h=1.0d0
double precision,Public :: Wab,Wba
double precision :: totaltime=0.0d0

double precision :: Beta      !beta=1/(kT)
logical, public :: inA
double precision, public :: FmA(3),FmB(3)
!double precision :: Dissipation_z=0.0d0,DissiFilter(0:2)=0.0d0,Fzo=0.0d0,zo=0.0d0,DissiFiltero(0:2)
!double precision :: DFilterRC=0.001590d0,tmpDissi=0.0d0

public :: ReadDissimap,NURBSdSetup, NURBSdSetGrid, NURBSdAllocate, NURBSdForce, dbgNURBSz,dbgNURBSz2,dbgAll
public :: Evolve!, UpdateDissipation,Dissipation_z,tmpDissi


contains


Subroutine ReadDissimapSimple(DissiFile)
  implicit none
  character*40, intent(in) :: DissiFile
  character*250 :: line
  integer :: i,j,k,info
  double precision :: tmparr(FDim),EA,EB,VA,VB,ET,fA(3),fB(3)
  
  Dijk=0.0d0
  
  open(file=trim(adjustl(DissiFile)),unit=123)
  do
    read(123,'(A)',iostat=info)line
    if(info /= 0) exit
    read(line,*)i,j,k,tmparr(:)
    !write(*,*)'read from file',i,j,k,tmparr(1:2)
    Dijk(i+1,j+1,k-1,:) = tmparr(:)
  enddo
  close(123)
	
	call RefineField()
  
  
end subroutine ReadDissimapSimple

Subroutine ReadDissimapSingle(DissiFile)
  implicit none
  character*40, intent(in) :: DissiFile
  character*450 :: line
  integer :: i,j,k,info
  double precision :: tmparr(FDim),EA,EB,VA,VB,ET,fA(3),fB(3)
  
  Dijk=0.0d0
  
  open(file=trim(adjustl(DissiFile)),unit=123)
  do
    read(123,'(A)',iostat=info)line
    if(info /= 0) exit
    read(line,*)k,tmparr(:)
    Dijk(2,2,k-1,:)=tmparr(:)
    write(*,*)'read from dfile',k,tmparr(:)
  enddo
  close(123)
	
	call RefineField()
  
  
end subroutine ReadDissimapSingle

subroutine RefineField()
	implicit none
	integer :: pX, pY, pdZ
	
	pX = Npts(1)-3
	pY = Npts(2)-3
	pdZ = Npts(3)
	
  !set the points replicas
  if( pX >= 3 ) then
		Dijk(0,2:pY+1,0:pdZ,:)=Dijk(pX,2:pY+1,0:pdZ,:)
		Dijk(1,2:pY+1,0:pdZ,:)=Dijk(pX+1,2:pY+1,0:pdZ,:)
		Dijk(pX+2,2:pY+1,0:pdZ,:)=Dijk(2,2:pY+1,0:pdZ,:)						!X periodicity (first point)
		Dijk(pX+3,2:pY+1,0:pdZ,:)=Dijk(3,2:pY+1,0:pdZ,:)					!X periodicity (second point)
		
		Dijk(0:pX+3,0,0:pdZ,:)   = Dijk(0:pX+3,pY+0,0:pdZ,:)
		Dijk(0:pX+3,1,0:pdZ,:)   = Dijk(0:pX+3,pY+1,0:pdZ,:)
		Dijk(0:pX+3,pY+2,0:pdZ,:)= Dijk(0:pX+3,2,0:pdZ,:)			!Y periodicity (first poit)
		Dijk(0:pX+3,pY+3,0:pdZ,:)= Dijk(0:pX+3,3,0:pdZ,:)			!Y periodicity (first poit)
	else
		!or just fill everything with the single given point
		Dijk(0,:,:,:) = Dijk(2,:,:,:)
		Dijk(1,:,:,:) = Dijk(2,:,:,:)
		Dijk(0, :, 0:pdZ, :) = Dijk(2, :, 0:pdZ, :)
		Dijk(1, :, 0:pdZ, :) = Dijk(2, :, 0:pdZ, :)
		Dijk(3, :, 0:pdZ, :) = Dijk(2, :, 0:pdZ, :)
		Dijk(4, :, 0:pdZ, :) = Dijk(2, :, 0:pdZ, :)
		
		Dijk(:, 0, 0:pdZ, :) = Dijk(:, 2, 0:pdZ, :)
		Dijk(:, 1, 0:pdZ, :) = Dijk(:, 2, 0:pdZ, :)
		Dijk(:, 3, 0:pdZ, :) = Dijk(:, 2, 0:pdZ, :)
		Dijk(:, 4, 0:pdZ, :) = Dijk(:, 2, 0:pdZ, :)
	endif
	
	
end subroutine RefineField


subroutine NURBSdAllocate(dx,dy,dz)
  implicit none
  integer, intent(in) :: dx,dy,dz
  
  NPts(1)=dx
  NPts(2)=dy
  NPts(3)=dz
  allocate(Dijk(0:dx,0:dy,0:dz,FDim))
  allocate(Pijk(0:dx,0:dy,0:dz,FDim))
  
  Dijk=0.0d0
  Pijk=0.0d0
  
  allocate(Gi(0:NPts(1)),Pi(0:NPts(1)))
  allocate(Gj(0:NPts(2)),Pj(0:NPts(2)))
  allocate(Gk(0:NPts(3)),Pk(0:NPts(3)))
  
  write(*,*)'Dissipative field allocated',dx,dy,dz
  
end subroutine NURBSdAllocate


subroutine NURBSdSetGrid(gridarray,sz,comp)
  implicit none
  double precision,intent(in) :: gridarray(0:sz)
  integer, intent(in) :: sz,comp
  
  if(comp==1) then
    Gi=gridarray
  endif
  if(comp==2) then
    Gj=gridarray
  endif
  if(comp==3) then
    Gk=gridarray
  endif

end subroutine NURBSdSetGrid


!*** ****** Full 4D NURBS interpolation ****** ***!
!   If this works i want double salary!
!
subroutine NURBSdSetup(temperature)
  use consts
  implicit none
  double precision, intent(in) :: temperature
  integer :: i,j,k,ii,jj,kk
  double precision :: x(3),f(3)
  double precision :: Rxjk(0:NPts(1),0:NPts(2),0:NPts(3),FDim),Rxyk(0:NPts(1),0:NPts(2),0:NPts(3),FDim)
  
  Beta=1.0d0/(temperature*kB)
  
  !*** create the parameters vectors
  allocate(usk(0:NPts(1)))     !there is one more point... itz the periodic image of the first one!
  call NURBSParams(usk,NPts(1))
  allocate(vsk(0:NPts(2)))     !there is one more point... itz the periodic image of the first one!
  call NURBSParams(vsk,NPts(2))
  allocate(wsk(0:NPts(3)))
  call NURBSParams(wsk,NPts(3))
  
  do i=1,3
    MPts(i)=NPts(i)+Order+1
  enddo
  
  !*** compute the knot vectors
  allocate(uk(0:MPts(1)))
  call NURBSKnots(uk,usk,NPts(1),MPts(1))
  allocate(vk(0:MPts(2)))
  call NURBSKnots(vk,vsk,NPts(2),MPts(2))
  allocate(wk(0:MPts(3)))
  call NURBSKnots(wk,wsk,NPts(3),MPts(3))
  
  allocate(Nb(0:Order))   !this will contain the non zero basis functions
  allocate(Nb3(3,0:Order))
  
  write(*,*)
  write(*,'(A)',advance='no')'Computing dissipation NURBS control points... '
  call NURBS_DoX(Rxjk)
  call NURBS_DoY(Rxjk,Rxyk)
  call NURBS_DoZ(Rxyk)
  write(*,'(A)')'done!'
	write(*,*)
	
  !call dbgAll()
  !stop
  !call dbgNURBSz2()
  
end subroutine NURBSdSetup


subroutine Evolve(dt,x,PA,Force)
  implicit none
  double precision, intent(in) :: dt,x(3)
  double precision, intent(out):: PA,Force(3)
  double precision :: exp_dh, t(3),inc,dec
  double precision :: Values(FDim)
  double precision :: preDer(2),EAgrad(3),EBgrad(3)
  integer :: i,j,k
  
  Force=0.0d0
  
  !write(*,*)'Evolve called!'
  !write(*,*)'z.position',x(3)
  !<Gk(0)
  !(x(3)<dble(4.6e-10)).or.
  if((x(3)<Gk(0)).or.(x(3)>=Gk(Npts(3)))) then
    PA = -1.0d0
    !write(*,*)'No dissi deal!', x(3)
    
    return
  endif
  
  !****evaluate ALL
  call NURBSdEvalX(x,Values,t)   !fetch the values for barriers and forces
  !write(*,*)'evolve2',values(:)
  
  !if the barrier is 0 there is only one min and it is B (tip close)
  if( Values(11) == 0.0d0 ) then
  	!write(*,*)'ZERO BARRIER FOUND!',x(3)
		PA_t = 0.0d0
		Force(1:3) = Values(8:10)  !!there is only B -> give the force of B!!!
		!stop
		return
	endif
  
  !check if we are in a dummy zone
!   if(isDummy(t)) then
! 		PA=-1.0d0
! 		!write(444,*)totaltime,x(3),0.0d0,0,0,0,0
! 		If(x(3)<5.0e-10)then
! 			PA_t=0.0d0
! 		else
! 			PA_t=1.0d0
! 		endif
! 		return
!   endif
  
!   write(*,*)'evolve3',t(:),beta
  
  !setup the transition rates
  Wab=Values(2)*Exp(-Beta*abs(Values(1)) )  !ABS IS HERE TEMPORARILY!
  Wba=Values(7)*Exp(-Beta*abs(Values(6)) )
  
  sumw_t = Wab + Wba   !the sum of the rates gives h(t)
	exp_dh = Exp(-(sumw_t+sumwo_t)*dt*0.50d0)
	exp_h  = exp_h * exp_dh
	
	PA_t = exp_h + PA_t*exp_dh + 0.50d0*dt*( Wba + Wbao_t*exp_dh )
  !Integrate the master equation
!   PAI_t = PAIo_t*exp_dh + 0.50d0*dt*(Wba + exp_dh*Wbao_t)
!   PA_t = exp_h+PAI_t
	
	!simple integrator...
	!PA_t = PAIo_t - 2.0d0*PA_t*sumw_t*dt + 2.0d0*Wba*dt


  if(PA_t>1.0d0) then   !1 is the maximum
  	PA_t=1.0d0
  	PAI_t=1.0d0-exp_h
  endif
  if(PA_t<0.0d0) then   !0 is the minimum
  	PA_t=0.0d0
  	PAI_t=-exp_h
  endif
  
  !totaltime=totaltime+dt
  !PAIo_t=PAI_t ORIGINAL!
  !PAIo_t = PA_t
  sumwo_t = (Wab+Wba)
  Wbao_t = Wba
  
  FmA = PA_t*values(3:5)
  FmB = (1.0d0-PA_t)*values(8:10)
  !FmB = values(3:5)
  !FmB = values(8:10)
  Force = FmA + FmB !-PA_t*EAgrad - (1.0d0-PA_t)*EBgrad
  PA=PA_t
  !write(*,*)'evolve6',Wab,Wba,PA_t,Force(:)
  !write(*,*)'Evolved!'
  !write(*,*)'RATES:',x(3),Wab,Wba, Wab-Wba,PA_t, -PA_t*sumw_t*dt + Wba*dt
  !write(444,*)totaltime,x(3),Values(2),Exp(-Beta*Values(1)),Values(7),Exp(-Beta*Values(6))
  
end subroutine Evolve


! subroutine UpdateDissipation(z,Fz,Freq,zTick,dt)
! 	use consts
! 	implicit none
! 	double precision, intent(in) :: z, Fz,Freq,dt
! 	logical, intent(in) :: zTick
! 	double precision :: tfq
! 	integer :: i
! 	
!   tmpDissi = tmpDissi+(Fz+Fzo)*(z-zo)/2.0d0
! 	
! 	if(zTick) then			!if the tip did a full cycle...
! 		DissiFilter(0)=tmpDissi*J2eV
! 		tmpDissi=0.0d0
! 		!write(888,*)DissiFilter(0)
! 		tfq=1.0d0/Freq
! 		do i=1,2					!process the filters
! 			DissiFilter(i)=DissiFilter(i-1)*(tfq/(DFilterRC+tfq)) + (1.0d0-tfq/(DFilterRC+tfq))*DissiFiltero(i)
! 			DissiFiltero(i)=DissiFilter(i)    !back collector
! 		enddo
! 		Dissipation_z=DissiFilter(2)
! 	endif
! 	
! 	zo=z
! 	Fzo=Fz
! 	
! end subroutine UpdateDissipation

subroutine NURBSdEvalX(x,vals,tout)
  implicit none
  double precision, intent(in) :: x(3)
  double precision, intent(out):: vals(FDim),tout(3)
  double precision :: t(3)
  integer :: idx(3),i,j,k,c
  
  
  t = 0.0d0
  if(Npts(1) > 4) t(1) = AccurateT(x(1),1)
  if(Npts(2) > 4) t(2) = AccurateT(x(2),2)
  
  t(1)=AccurateT(x(1),1)!(x(1)-Gi(0))/(Gi(Npts(1))-Gi(0))   !first estimate for tx
  t(2)=AccurateT(x(2),2)!(x(2)-Gj(0))/(Gj(Npts(2))-Gj(0))   !first estimate for ty
  t(3)=AccurateT(x(3),3)!(x(3)-Gk(0))/(Gk(Npts(3))-Gk(0))   !first estimate for tz
  !t(3)=AccurateTZ(x(3))         !refine the value for Z
  
  if((t(1)>1.0d0).or.(t(2)>1.0d0).or.(t(3)>1.0d0) ) then
		write(*,*)'Error in t!',t,x
		write(*,*)'Gk(f)',Gk(Npts(3))
		stop
	endif
  
  call NURBSdEval(t,vals)
  tout=t
  
end subroutine NURBSdEvalX

subroutine NURBSdEval(t,vals)
  implicit none
  double precision, intent(in) :: t(3)
  double precision, intent(out):: vals(FDim)
  double precision :: si(FDim),sj(FDim),sk(FDim)
  integer :: idx(3),i,j,k
  
  vals=0.0d0
  si=0.0d0
  sj=0.0d0
  sk=0.0d0
  idx=0
  
  if((t(1)>1.0d0).or.(t(2)>1.0d0).or.(t(3)>1.0d0) ) then
		write(*,*)'Error in t (eval)!'
		stop
	endif
	
  idx(1)=NURBSSpan(NPts(1),t(1),uk)       !find the span of tx
  idx(2)=NURBSSpan(NPts(2),t(2),vk)       !find the span of ty
  idx(3)=NURBSSpan(NPts(3),t(3),wk)       !find the span of tz
  
  call NURBSBasis(idx(1),t(1),1,uk)
  Nb3(1,0:Order)=Nb(0:Order)
  call NURBSBasis(idx(2),t(2),2,vk)
  Nb3(2,0:Order)=Nb(0:Order)
  call NURBSBasis(idx(3),t(3),3,wk)
  Nb3(3,0:Order)=Nb(0:Order)
  
  si=0.0d0
  sj=0.0d0
  sk=0.0d0
  do i=0,Order
    sj=0.0d0
    do j=0,Order
      sk=0.0d0
      do k=0,Order
        sk(:)=sk(:)+Nb3(3,k)*Pijk(idx(1)-Order+i,idx(2)-Order+j,idx(3)-Order+k,:)!zcoords(idx(3)-Order+k)!
      enddo
      sj(:)=sj(:)+Nb3(2,j)*sk(:)
    enddo
    vals(:)=vals(:)+Nb3(1,i)*sj(:)
  enddo
  
  
  
end subroutine NURBSdEval


subroutine NURBSdGrad(t0,E0,aorb,grad)
  implicit none
  double precision, intent(in) :: t0(3),E0
  double precision, intent(out):: grad(3)
  integer, intent(in) :: aorb
  double precision :: t(3),vals(FDim),tmpv(2),x0,x1
  
  !move to t(1)+h
!   write(*,*)'--- gradientize ---',E0,t0(:)
  t=t0
  x0 = GetX(t(1))
  !write(*,*)'G1',t(:)
  t(1)=t(1)+0.1*(Gi(NPts(1))-Gi(NPts(1)-1))
  x1 = GetX(t(1))
  !write(*,*)'G2',t(:)
  call NURBSdEval(t,vals)
  grad(1)= ((vals(5)-vals(aorb))-E0)/(x1-x0)
!   write(*,*)'grad:',vals(5)-vals(aorb),E0,vals(5)-vals(aorb)-E0,'on',(x1-x0)
  t=t0
  x0 = GetY(t(2))
  t(2)=t(2)+0.1*(Gj(NPts(2))-Gj(NPts(2)-1))
  x1 = GetY(t(2))
  call NURBSdEval(t,vals)
  grad(2)= ((vals(5)-vals(aorb))-E0)/(x1-x0)
!   write(*,*)'grad:',vals(5)-vals(aorb),E0,vals(5)-vals(aorb)-E0,'on',(x1-x0)
  t=t0
  x0 = GetZ(t(3))
  t(3)=t(3)+0.1*(Gk(NPts(3))-Gk(NPts(3)-1))
  x1 = GetZ(t(3))
  call NURBSdEval(t,vals)
  grad(3)= ((vals(5)-vals(aorb))-E0)/(x1-x0)
!   write(*,*)'grad:',vals(5)-vals(aorb),E0,vals(5)-vals(aorb)-E0,'on',(x1-x0)
!   write(*,*)'--- ---'
  
  
  
end subroutine NURBSdGrad

function ParamSpan(t)
  implicit none
  double precision, intent(in) :: t(3)
  integer :: i,j,k
  integer :: ParamSpan(3)
  
  do i=0,NPts(1)
    if(usk(i)>t(1)) exit
  enddo
  do j=0,NPts(2)
    if(vsk(j)>t(2)) exit
  enddo
  do k=0,NPts(3)
    if(wsk(k)>t(3)) exit
  enddo
  ParamSpan(1)=i-1
  ParamSpan(2)=j-1
  ParamSpan(3)=k-1
  
end function ParamSpan


logical function isDummy(t)
  implicit none
  double precision, intent(in) :: t(3)
  integer :: idx(3)
  double precision :: refval
  
  idx=ParamSpan(t)
  !write(*,*)idx(:)
  refval=Dijk(idx(1),idx(2),idx(3),1)
  
  isDummy=.false.
  !no dummy if the reference point has a non zero value
  if((refval/=0.0d0).or.(idx(3)==0)) then
    isDummy=.false.
    return
  else
  	isDummy=.true.
  	return
  endif
  
  if(refval==Dijk(idx(1),idx(2),idx(3)-2,1)) then
    isDummy=.true.
    return
  endif
  
end function isDummy

double precision function AccurateTZ(z)
  implicit none
  double precision, intent(in) :: z
  double precision :: t,zl,zm,zr
  double precision :: lef,mid,rig
  integer :: i,j,k,kk,side=0
  
  side=0
  t=0.0d0
  lef= 0.0d0
  rig= 1.0d0
  zl=Gk(0)-z
  zr=Gk(NPts(3))-z

  do
    
    mid = (zl*rig-zr*lef)/(zl-zr)
    t=mid
    
    zm=GetZ(t)-z
    
    if(zm*zr>0.0d0) then
      rig=mid
      zr=zm
      if(side==-1) zl=zl/2.0d0
      side=-1
    endif
    if(zl*zm>0.0d0) then
      lef=mid
      zl=zm
      if(side==1) zr=zr/2.0d0
      side=1
    endif
    if(Abs(zm)<1.0e-15) exit
    if(abs(rig-lef)<1.0e-8) exit
    
  enddo
  AccurateTZ=mid
  
end function AccurateTZ


double precision function AccurateT(z,comp)
  implicit none
  integer, intent(in) :: comp
  double precision, intent(in) :: z
  double precision :: t,zl,zm,zr
  double precision :: lef,mid,rig
  integer :: i,j,k,kk,side=0
  
  side=0
  t=0.0d0
  lef= 0.0d0
  rig= 1.0d0
  
  if(comp == 1) then
		zl=Gi(0)-z
		zr=Gi(NPts(1))-z
  endif
  if(comp == 2) then
		zl=Gj(0)-z
		zr=Gj(NPts(2))-z
  endif
  if(comp == 3) then
		zl=Gk(0)-z
		zr=Gk(NPts(3))-z
  endif


  do
    
    mid = (zl*rig-zr*lef)/(zl-zr)
    t=mid
    
    
    if(comp == 1)  zm=GetX(t)-z
    if(comp == 2)  zm=GetY(t)-z
    if(comp == 3)  zm=GetZ(t)-z
    
    if(zm*zr>0.0d0) then
      rig=mid
      zr=zm
      if(side==-1) zl=zl/2.0d0
      side=-1
    endif
    if(zl*zm>0.0d0) then
      lef=mid
      zl=zm
      if(side==1) zr=zr/2.0d0
      side=1
    endif
    if(Abs(zm)<1.0e-15) exit
    if(abs(rig-lef)<1.0e-8) exit
    
  enddo
  AccurateT=mid
  
end function AccurateT


double precision function GetX(t)

  implicit none
  double precision, intent(in) :: t
  double precision :: z
  integer :: idx,i
  
  z=0.0d0
  
  idx=NURBSSpan(NPts(1),t,uk)
  call NURBSBasis(idx,t,1,uk)
  
  do i=0,Order
    z=z+Nb(i)*Pi(idx-Order+i)
  enddo
  
  GetX=z
  
end function GetX
double precision function GetY(t)

  implicit none
  double precision, intent(in) :: t
  double precision :: z
  integer :: idx,i
  
  z=0.0d0
  
  idx=NURBSSpan(NPts(2),t,vk)
  call NURBSBasis(idx,t,2,vk)
  
  do i=0,Order
    z=z+Nb(i)*Pj(idx-Order+i)
  enddo
  
  GetY=z
  
end function GetY
double precision function GetZ(t)

  implicit none
  double precision, intent(in) :: t
  double precision :: z
  integer :: idx,i
  
  z=0.0d0
  
  idx=NURBSSpan(NPts(3),t,wk)
  call NURBSBasis(idx,t,3,wk)
  
  do i=0,Order
    z=z+Nb(i)*Pk(idx-Order+i)
  enddo
  
  GetZ=z
  
end function GetZ

!*** ****** Evaluate the basis functions ****** ***
!   Input variables:
!   (INT)idx           knot span inter
!   (INT)c             spatial component (1->x,2->y,3->z)
!   (DBL)ukv(:)        knot vector for this spatial direction
!
!   Description:
!   Computes the non zero basis functions in t along the c-space direction
!
subroutine NURBSBasis(idx,t,c,ukv)
  implicit none
  integer :: idx
  integer, intent(in) :: c
  double precision, intent(in) :: t
  double precision, intent(in) :: ukv(0:MPts(c))
  double precision :: left(0:Order),right(0:Order),saved,temp
  integer :: j,r
  
  left=0.0d0
  right=0.0d0
  Nb=0.0d0
  Nb(0)=1.0d0
  
  do j=1,Order
    left(j)=t-ukv(idx-j+1)
    right(j)= ukv(idx+j)-t
    saved=0.0d0
    do r=0,j-1
      temp=Nb(r)/(right(r+1)+left(j-r))
      Nb(r)=saved+right(r+1)*temp
      saved=left(j-r)*temp
    enddo
    Nb(j)=saved
  enddo
  
end subroutine NURBSBasis

!*** ****** Create a knot vector ****** ***!
!   Input variables:
!   (INT)n             dimension of usk - number of points
!   (INT)m             dimension of uk - m=n+p+1
!   (INT)p             degree of the NURBS
!   (DBL)uk(m)         empty knot vector (dimension m should be n+p+1)
!   (DBL)usk(n)        parameters vector of dimension n
!
!   Output variables:
!   DBL)uk(m)          the cool knot vector!
!
!   Description:
!   generates a knot vector by averaging the parameter vector.
!
subroutine NURBSKnots(ukv,uskv,n,m)
  implicit none
  integer, intent(in) :: n,m
  double precision :: ukv(0:m),uskv(0:n)
  integer :: i
  
  do i=0,Order
    ukv(i)=0.0d0
  enddo
  do i=m-Order,m
    ukv(i)=1.0d0
  enddo
  do i=1,n-Order
    ukv(i+Order)=sum(uskv(i:i+Order-1))/dble(Order)
  enddo

end subroutine NURBSKnots


!*** ****** Create Parameters Vector equally spaced!! ****** ***!
!   Input variables:
!   (INT)n             dimension of u - number of points
!   (DBL)u(n)          enpty parameters vector of dimension n
!
!   Output variables:
!   DBL)u(n)          the cool parameters vector!
!
!   Description:
!   generates a parameters vector using the equally spaced method.
!
subroutine NURBSParams(u,n)
  implicit none
  integer, intent(in) :: n
  double precision,intent(inout) :: u(0:n)
  integer :: i
  
  !*** create the parameters vectors - Equally spaced method
  u(0)=0.0d0
  do i=1,n-1
    u(i)=dble(i)/dble(n)
  end do
  u(n)=1.0d0
  
end subroutine NURBSParams


!*** ****** Find the knot span ****** ***!
!
integer function NURBSSpan(n,t,ukv)
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: t
  double precision, intent(in) :: ukv(0:n+Order+1)
  integer :: low, mid, hig
  

  if(t>=ukv(n+1)) then
    NURBSSpan=n
    return
  endif
  
  low=Order
  hig=n+1
  mid=(low+hig)/2
  do
    
    if((t>=ukv(mid)).and.(t<ukv(mid+1))) then
      NURBSSpan=mid
      exit
    endif
    
    if(t<ukv(mid)) then
      hig=mid
    else
      low=mid
    endif
    mid=(low+hig)/2
  enddo
  
  
end function NURBSSpan



subroutine NURBS_DoX(Rxjk)
  implicit none
  double precision, intent(out) :: Rxjk(0:NPts(1),0:NPts(2),0:NPts(3),FDim)
  double precision :: linsys(0:NPts(1),0:NPts(1)),mat(0:NPts(1),0:NPts(1))
  double precision :: b(0:NPts(1)),IPIV(0:NPts(1))
  integer :: INFO,i,j,k,c,ii
  
  
  IPIV=0.0d0
  b=0.0d0
  linsys=0.0d0
  mat=0.0d0
  do i=0,NPts(1)       !loop on the rows
    ii=NURBSSpan(NPts(1),usk(i),uk)    !find the span of usk(i)
    call NURBSBasis(ii,usk(i),1,uk)
    do j=0,Order
      linsys(i,ii-Order+j)=Nb(j)
      mat(i,ii-Order+j)=Nb(j)
    enddo
  enddo
  Pi(:)=Gi(:)
  call DGESV( NPts(1)+1, 1, mat, NPts(1)+1, IPIV, Pi, NPts(1)+1, INFO )
  
  do c=1,FDim    !loop on components
    do j=0,NPts(2)
      do k=0,NPts(3)
        
        b(0:Npts(1))=Dijk(0:NPts(1),j,k,c)
        IPIV=0.0d0
        INFO=0
        mat=linsys
        call DGESV( NPts(1)+1, 1, mat, NPts(1)+1, IPIV, b, NPts(1)+1, INFO )
        Rxjk(0:NPts(1),j,k,c)=b(0:NPts(1))
      enddo
      
    enddo
    
  enddo
  
  
end subroutine NURBS_DoX


subroutine NURBS_DoY(Rxjk,Rxyk)
  implicit none
  double precision, intent(in) :: Rxjk(0:Npts(1),0:Npts(2),0:Npts(3),FDim)
  double precision, intent(out):: Rxyk(0:Npts(1),0:Npts(2),0:Npts(3),FDim)
  double precision :: linsys(0:NPts(2),0:NPts(2)),mat(0:NPts(2),0:NPts(2))
  double precision :: b(0:NPts(2)),IPIV(0:NPts(2))
  integer :: INFO,i,j,k,c,ii
  
  IPIV=0.0d0
  b=0.0d0
  linsys=0.0d0
  mat=0.0d0
  do i=0,NPts(2)       !loop on the rows
    ii=NURBSSpan(NPts(2),vsk(i),vk)    !find the span of usk(i)
    call NURBSBasis(ii,vsk(i),2,vk)
    do j=0,Order
      linsys(i,ii-Order+j)=Nb(j)
      mat(i,ii-Order+j)=Nb(j)
    enddo
  enddo
  Pj(:)=Gj(:)
  call DGESV( NPts(2)+1, 1, mat, NPts(2)+1, IPIV, Pj, NPts(2)+1, INFO )
  
  do c=1,FDim    !loop on components
    do i=0,NPts(1)
      do k=0,NPts(3)
        
        b(0:Npts(2))=Rxjk(i,0:Npts(2),k,c)
        IPIV=0.0d0
        INFO=0
        mat=linsys
        call DGESV( NPts(2)+1, 1, mat, NPts(2)+1, IPIV, b, NPts(2)+1, INFO )
        Rxyk(i,0:NPts(2),k,c)=b(0:NPts(2))
      enddo
      
    enddo
    
  enddo
  
end subroutine NURBS_DoY


subroutine NURBS_DoZ(Rxyk)
  implicit none
  double precision, intent(in) :: Rxyk(0:Npts(1),0:Npts(2),0:Npts(3),FDim)
  double precision :: linsys(0:NPts(3),0:NPts(3)),mat(0:NPts(3),0:NPts(3))
  double precision :: b(0:NPts(3)),IPIV(0:NPts(3))
  integer :: INFO,i,j,k,c,ii
  
  IPIV=0.0d0
  b=0.0d0
  linsys=0.0d0
  mat=0.0d0
  do i=0,NPts(3)       !loop on the rows
    ii=NURBSSpan(NPts(3),wsk(i),wk)    !find the span of usk(i)
    call NURBSBasis(ii,wsk(i),3,wk)
    do j=0,Order
      linsys(i,ii-Order+j)=Nb(j)
      mat(i,ii-Order+j)=Nb(j)
    enddo
  enddo
  Pk(:)=Gk(:)

  call DGESV( NPts(3)+1, 1, mat, NPts(3)+1, IPIV, Pk, NPts(3)+1, INFO )

  do c=1,FDim    !loop on components
    do i=0,NPts(1)
      do j=0,NPts(2)
        b(0:Npts(3))=Rxyk(i,j,0:Npts(3),c)
        IPIV=0.0d0
        INFO=0
        mat=linsys
        call DGESV( NPts(3)+1, 1, mat, NPts(3)+1, IPIV, b, NPts(3)+1, INFO )
        Pijk(i,j,0:NPts(3),c)=b(0:NPts(3))
      enddo
    enddo
  enddo

end subroutine NURBS_DoZ


subroutine dbgNURBSz()
  implicit none
  double precision :: x(3), f(FDim),tt(3)
  double precision :: z=0.0d0,xc=0.0d0,yc=0.0d0
  integer :: i
  
  x=0.0d0
  x(3)=5.0e-10
  do yc=Gj(0),Gj(NPts(2)),5.640e-10/32.0d0
    do xc=Gi(0),Gi(NPts(1)),5.640e-10/32.0d0
      x(2) = yc
      x(1) = xc
      call NURBSdEvalX(x,f,tt)
      write(123,*)xc,yc,f(:)
    enddo
    write(123,*)
  enddo
  do i=0,NPts(3)
    write(124,*)Gk(i),Dijk(0,0,i,:)
  enddo
  
  
end subroutine dbgNURBSz

subroutine dbgNURBSz2()
  implicit none
  double precision :: x(3), f(FDim),tt(3)
  double precision :: z=0.0d0,xc=0.0d0,yc=0.0d0
  integer :: i
  
  x=0.0d0
  !x(3)=5.0e-10
  do yc=Gk(0),Gk(NPts(3)),0.01e-10
      x(3)=yc
      call NURBSdEvalX(x,f,tt)
      write(123,*)yc,f(:)
  enddo
!   do yc=Gk(NPts(3)),Gk(0),-0.01e-10
!       x(3)=yc
!       call NURBSdEvalX(x,f,tt)
!       write(124,*)yc,f(2)*Exp(-Beta*f(1)),f(4)*Exp(-Beta*f(3))
!   enddo
!   
  do i=0,NPts(3)
    write(124,*)Gk(i),Dijk(0,0,i,:)
  enddo
  
  
end subroutine dbgNURBSz2


subroutine dbgAll()
	implicit none
	double precision :: x(3), f(FDim),tt(3),xd(3)
	double precision :: xc,yc,zc
	character*100 :: fname
	integer :: i,j,k
	
	x=0.0d0
	!write(*,*)Gi(:)
	!write(*,*)Gj(:)
	i=2
	!do i=2,(NPts(1)-3)
		write(fname,'(A4 I1 A4)')'map_',i-2,'.dat'
		open(file=trim(adjustl(fname)),unit=321)
		do j=2, NPts(2)
		
			do k=0,(NPts(3))
				write(321,*)Gj(j),Gk(k),Dijk(i,j,k,:)
			enddo
			write(321,*)
			
		enddo
		close(321)
		write(fname,'(A5 I1 A4)')'imap_',i-2,'.dat'
		open(file=trim(adjustl(fname)),unit=322)
		x(1)=Gi(i)
		do yc=Gj(0),2.0d0*5.640e-10,2.0d0*5.640e-10/32
			x(2)=yc
			!write(*,*)'original pos',x(:)
			do zc=Gk(0),Gk(NPts(3)),0.01e-10
				
				x(3)=zc
				!write(*,*)'original pos',x(:)
				call dCenterCursor(x,xd)
				!write(*,*)'c pos',xd(:)
				call NURBSdEvalX(xd,f,tt)
				write(322,*)yc,zc,f(:)
			enddo
			write(*,*)'c pos',x(2),xd(2)
			write(322,*)
		enddo
		close(322)
		
	!enddo
! 	
end subroutine dbgAll

subroutine dCenterCursor(x,xc)
  implicit none
  double precision, intent(in) :: x(3)
  double precision, intent(out) :: xc(3)
  
  
  !center the cursor in the x direction
  xc(1)=x(1)
    !the pacman symmetry
    do
      if((xc(1)>dble(5.64e-10)).or.(xc(1)<0.0d0)) then 
				xc(1)=xc(1)-(x(1)/Abs(x(1)))*dble(5.64e-10)
      else
				exit
      endif
    enddo

  !center the cursor in the y direction
  xc(2)=x(2)
    !the pacman symmetry
    do
      if((xc(2)>dble(5.64e-10)).or.(xc(2)<0.0d0)) then 
				xc(2)=xc(2)-(x(2)/abs(x(2)))*dble(5.64e-10)
      else
				exit
      endif
    enddo
    
      
  xc(3)=x(3)
  
  return

end subroutine dCenterCursor



end module NURBSDissiField

