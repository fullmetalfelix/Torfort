module NURBSForceField
  
  !Private
  !*** this module contains all the stuff to get the conservative force in (x,y,z)

integer, public :: Order=1                                        !NURBS order
double precision, allocatable :: Gi(:),Gj(:),Gk(:)        !gridpoint space coords
double precision, allocatable, public :: Dijk(:,:,:,:)            !datapoints on the grid
integer :: NPts(3),MPts(3)

double precision, allocatable :: Pijk(:,:,:,:)            !NURBS control points
double precision, allocatable :: Pi(:),Pj(:),Pk(:)        !Control points space coords

double precision, allocatable :: usk(:),vsk(:),wsk(:)     !NURBS parameters vectors
double precision, allocatable :: uk(:),vk(:),wk(:)        !NURBS knots vectors
double precision, allocatable :: Nb(:)                    !NURBS basis functions
double precision, allocatable :: Nb3(:,:)

public :: NURBSfSetup,NURBSfSetGrid, NURBSfAllocate, NURBSForce, dbgNURBSz,dbgdbg


contains


subroutine NURBSfAllocate(dx,dy,dz)
  implicit none
  integer, intent(in) :: dx,dy,dz
  
  NPts(1)=dx
  NPts(2)=dy
  NPts(3)=dz
  allocate(Dijk(0:dx,0:dy,0:dz,3))
  allocate(Pijk(0:dx,0:dy,0:dz,3))
  
  Dijk=0.0d0
  
  allocate(Gi(0:NPts(1)),Pi(0:NPts(1)))
  allocate(Gj(0:NPts(2)),Pj(0:NPts(2)))
  allocate(Gk(0:NPts(3)),Pk(0:NPts(3)))
  
  
end subroutine NURBSfAllocate

subroutine NURBSfSetGrid(gridarray,sz,comp)
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

end subroutine NURBSfSetGrid



!*** ****** Full 4D NURBS interpolation ****** ***!
!   If this works i want double salary!
!
subroutine NURBSfSetup()
  implicit none
  integer :: i,j,k,ii,jj,kk
  double precision :: x(3),f(3),t
  double precision :: Rxjk(0:NPts(1),0:NPts(2),0:NPts(3),3),Rxyk(0:NPts(1),0:NPts(2),0:NPts(3),3)
  !write(*,*)'points:',px,py,pz
  !write(*,*)ForceField(1,1,:,3)
  
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
  write(*,'(A)',advance='no')'Computing forcefield NURBS control points... '
  call NURBS_DoX(Rxjk)
  call NURBS_DoY(Rxjk,Rxyk)
  call NURBS_DoZ(Rxyk)
  write(*,'(A)')'done!'
  write(*,*)
  
end subroutine NURBSfSetup


!*** ****** Evaluate the force in (x,y,z) ****** ***
!
subroutine NURBSForce(x,force)
  implicit none
  double precision, intent(in) :: x(3)
  double precision, intent(out):: force(3)
  double precision :: t(3),si(3),sj(3),sk(3)
  integer :: idx(3),i,j,k,c
  
  force=0.0d0
  si=0.0d0
  sj=0.0d0
  sk=0.0d0
  t=0.0d0
  idx=0
  !write(*,*)x
  
  !find the parameters that gives the right x
  !write(*,*)'starting accurate...'
  t(1)=AccurateT(x(1),1)!(x(1)-Gi(0))/(Gi(Npts(1))-Gi(0))   !first estimate for x
  t(2)=AccurateT(x(2),2)!(x(2)-Gj(0))/(Gj(Npts(2))-Gj(0))   !first estimate for y
  t(3)=AccurateT(x(3),3)!(x(3)-Gk(0))/(Gk(Npts(3))-Gk(0))   !first estimate for z
  !write(*,*)'done!'
!   write(*,*)'Accurate Z...'
  !t(3)=AccurateTZ(x(3))
  !t(1) = 0.0d0
  !t(2) = 0.0d0
   !write(*,*)'got Z!',t(:)
	if( (t(3)<0.0d0).or.(t(3)>1.0d0) ) then
		write(*,*)'ARRGH t out:',t(:),x(:)
  endif
   
  idx(1)=NURBSSpan(NPts(1),t(1),uk)       !find the span of tx
!   write(*,*)'spans 1',idx(1)
  idx(2)=NURBSSpan(NPts(2),t(2),vk)       !find the span of ty
!   write(*,*)'spans 2',idx(2)
  idx(3)=NURBSSpan(NPts(3),t(3),wk)       !find the span of tz
!   write(*,*)'spans 3',idx(3)
!   write(*,*)'spans ok'
  
  call NURBSBasis(idx(1),t(1),1,uk)
  Nb3(1,0:Order)=Nb(0:Order)
  call NURBSBasis(idx(2),t(2),2,vk)
  Nb3(2,0:Order)=Nb(0:Order)
  call NURBSBasis(idx(3),t(3),3,wk)
  Nb3(3,0:Order)=Nb(0:Order)
!   write(*,*)'basis ok'
  
  do c=1,3
    do i=0,Order
      sj(c)=0.0d0
      do j=0,Order
        sk(c)=0.0d0
        do k=0,Order
          sk(c)=sk(c)+Nb3(3,k)*Pijk(idx(1)-Order+i,idx(2)-Order+j,idx(3)-Order+k,c)!zcoords(idx(3)-Order+k)!
        enddo
        sj(c)=sj(c)+Nb3(2,j)*sk(c)
      enddo
      force(c)=force(c)+Nb3(1,i)*sj(c)
    enddo
  enddo
  !force(1) = 0.0d0
  !force(2) = 0.0d0
  !write(*,*)force
  
end subroutine NURBSForce


double precision function AccurateTZ(z)
  implicit none
  double precision, intent(in) :: z
  double precision :: t,zl,zm,zr
  double precision :: lef,mid,rig
  integer :: i,j,k,kk,side=0
  
  side=0
  t=0.0d0
  !write(*,*)z
  lef=0.0d0
  rig=1.0d0
  zl=Gk(0)-z
  zr=Gk(Npts(3))-z
  
  
  do
    
    mid = (zl*rig-zr*lef)/(zl-zr)
    t=mid
    !write(*,*)lef,rig,mid
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
  double precision, intent(out) :: Rxjk(0:NPts(1),0:NPts(2),0:NPts(3),3)
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
  
  do c=1,3    !loop on components
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
  double precision, intent(in) :: Rxjk(0:Npts(1),0:Npts(2),0:Npts(3),3)
  double precision, intent(out):: Rxyk(0:Npts(1),0:Npts(2),0:Npts(3),3)
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
  
  do c=1,3    !loop on components
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
  double precision, intent(in) :: Rxyk(0:Npts(1),0:Npts(2),0:Npts(3),3)
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
  
  do c=1,3    !loop on components
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


!*** ****** Read the ForceField File ****** ***!
!   Description:
!   this routine reads the forcefield input file when it is written in the
!   simple plain form. The first 3 columns are the matrix indexes of the gridpoints
!   and the last 3 are the force components. Very straightforward.
!
subroutine ReadForceFieldSimple(ForceFile)
	use consts
  implicit none
  character*40, intent(in) :: ForceFile
  integer :: ifile=999, ciph=1
  integer :: i,j,k,ioer
  character*200 :: line
  double precision :: fx,fy,fz
  
!   ForceFile=trim(adjustl(ForceFile))
!   allocate(ForceField(pX,pY,pZ+1,3))  !the is one more z-plane where all the forces are 0!!!
!   ForceField=0.0d0
  
  open(file=trim(adjustl(ForceFile)),unit=ifile)
  write(*,*)'  Reading simple: '//trim(adjustl(ForceFile))
  do
  	read(ifile, '(A)', iostat=ioer) line   !read a line
  	if(ioer /= 0) exit
		read(line,*)i,j,k,fx,fy,fz  !read the data in that line
		Dijk(i+1,j+1,k-1,1)=fx
		Dijk(i+1,j+1,k-1,2)=fy            !set the forcefield matrix
		Dijk(i+1,j+1,k-1,3)=fz*eVang2N
		!write(*,*)'ASDASD'
		!write(*,*)i,j,k,fz
	enddo
	
  close(ifile)      !the end of the file was reached!
  call RefineField()

end subroutine ReadForceFieldSimple

subroutine ReadForceFieldSingle(ForceFile)
	use consts
  implicit none
  character*40, intent(in) :: ForceFile
  integer :: ifile=999, ciph=1
  integer :: i,j,k,ioer
  character*200 :: line
  double precision :: fx,fy,fz

  open(file=trim(adjustl(ForceFile)),unit=ifile)
  write(*,*)'  Reading single: '//trim(adjustl(ForceFile))
  !write(*,*)Npts(:)
  do
  	read(ifile, '(A)', iostat=ioer) line   !read a line
  	if(ioer /= 0) exit
		read(line,*)k,fx,fy,fz  !read the data in that line
		!write(*,*)k,fz
		Dijk(2,2,k-1,1)=fx
		Dijk(2,2,k-1,2)=fy            !set the forcefield matrix
		Dijk(2,2,k-1,3)=fz*eVang2N
		
	enddo
	
  close(ifile)      !the end of the file was reached!
  call RefineField()

end subroutine ReadForceFieldSingle
!***************************************************************************








!setup the periodic boundary point replicas
subroutine RefineField()
	implicit none
	integer :: pX, pY, pZ
	
	pX = Npts(1)-3
	pY = Npts(2)-3
	pZ = Npts(3)
	
  !if the field is not on a single point, create periodic boundary
  if( pX >= 3 ) then
		Dijk(0,2:pY+1,0:pZ-1,1:3)=Dijk(pX,  2:pY+1,0:pZ-1,1:3)!ForceField(pX-1,1:pY,1:pZ,1:3)
		Dijk(1,2:pY+1,0:pZ-1,1:3)=Dijk(pX+1,2:pY+1,0:pZ-1,1:3)!ForceField(pX,1:pY,1:pZ,1:3)
		Dijk(pX+2,2:pY+1,0:pZ-1,1:3)=Dijk(2,2:pY+1,0:pZ-1,1:3)					!X periodicity (first point)
		Dijk(pX+3,2:pY+1,0:pZ-1,1:3)=Dijk(3,2:pY+1,0:pZ-1,1:3)					!X periodicity (second point)
		
		Dijk(0:pX+3,0,0:pZ-1,1:3)   = Dijk(0:pX+3,pY+0,0:pZ-1,1:3)
		Dijk(0:pX+3,1,0:pZ-1,1:3)   = Dijk(0:pX+3,pY+1,0:pZ-1,1:3)
		Dijk(0:pX+3,pY+2,0:pZ-1,1:3)= Dijk(0:pX+3,2,0:pZ-1,1:3)			!Y periodicity (first poit)
		Dijk(0:pX+3,pY+3,0:pZ-1,1:3)= Dijk(0:pX+3,3,0:pZ-1,1:3)			!Y periodicity (first poit)
	else
		!otherwise just create fill all the points with the same values
		Dijk(0, :, 0:pZ-1, 1:3) = Dijk(2, :, 0:pZ-1, 1:3)
		Dijk(1, :, 0:pZ-1, 1:3) = Dijk(2, :, 0:pZ-1, 1:3)
		Dijk(3, :, 0:pZ-1, 1:3) = Dijk(2, :, 0:pZ-1, 1:3)
		Dijk(4, :, 0:pZ-1, 1:3) = Dijk(2, :, 0:pZ-1, 1:3)
		
		Dijk(:, 0, 0:pZ-1, 1:3) = Dijk(:, 2, 0:pZ-1, 1:3)
		Dijk(:, 1, 0:pZ-1, 1:3) = Dijk(:, 2, 0:pZ-1, 1:3)
		Dijk(:, 3, 0:pZ-1, 1:3) = Dijk(:, 2, 0:pZ-1, 1:3)
		Dijk(:, 4, 0:pZ-1, 1:3) = Dijk(:, 2, 0:pZ-1, 1:3)
	endif
  
  Dijk(0:pX+3,0:pY+3,pZ,1:3)=0.0d0			!set a zero point at large distance
	
	
	
end subroutine RefineField



subroutine dbgNURBSz()
  implicit none
  double precision :: x(3), f(3)
  double precision :: z=0.0d0
  integer :: i
  
  x=0.0d0
  do z=dble(21.0e-10), dble(1.0e-10), dble(-0.1e-10)
    x(3) = z
    call NURBSForce(x,f)
    write(123,*)z,f(3)
  enddo
  do i=0,NPts(3)
    write(124,*)Gk(i),Dijk(0,0,i,3)
  enddo
  
  
end subroutine dbgNURBSz

subroutine dbgdbg()
  implicit none
  double precision :: x(3),f(3),t,z
  
  x=0.0d0
  x(3)=5.0e-10
  do z=0.0d0,5.60e-10,5.64e-10/32.0d0
    do t=0.0d0,5.60e-10,5.64e-10/32.0d0
      x(2)=t
      x(1)=z
      call NURBSForce(x,f)
      write(333,*)z,t,f
    enddo
      write(333,*)
  enddo
  x(3)=5.050e-10
  do z=0.0d0,5.60e-10,5.64e-10/32.0d0
    do t=0.0d0,5.60e-10,5.64e-10/32.0d0
      x(2)=t
      x(1)=z
      call NURBSForce(x,f)
      write(444,*)z,t,f
    enddo
      write(444,*)
  enddo

end subroutine dbgdbg


end module NURBSForceField

