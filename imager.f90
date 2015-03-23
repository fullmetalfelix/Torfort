!*** ****** IMAGING PARAMETERS ****** ***!
module Imager

integer :: imgPix                !Image Pixels per side
integer :: imgPixs               !Image pixels on the slow scan direction
double precision :: imgSide      !Real space Image side

double precision :: ScanSpeed     !Scan speed
double precision :: RepoSpeed     !Repositioning speed
integer :: ScanAxis=2             !Scan direction (1->X, 2->Y)
integer :: ScanSlow=1             !Slow scan direction (1->X, 2->Y)


character*8 :: ImageMode        !Imaging mode (FMAFM - TRFFM)
integer :: ImagePattern         !Pattern mode (1->lineline, 2->zigozago, 3->forwback)
character*200 :: ScriptFile			!Scanning Script filename
integer :: ScriptLines					!number of instructions in the script
!integer,allocatable :: ScriptCommands(:)			!command sequence


integer :: scanQ(2)=1           !pixel coordinates (x,y)
integer :: cursorStep=1         !pixel cursor step, used to do forward and backward scan (1 or -1)

!--- Parameters for Scanliner ---!
double precision :: SL_start(2)
double precision :: SL_fdir(2) =  (/ 0.0d0, 1.0d0 /)
double precision :: SL_sdir(2) =  (/ 1.0d0, 0.0d0 /)
double precision :: SL_Length = dble(2.82E-9)
integer :: SL_npts = 100

!--- Parameters for ZSpectrer ---!
double precision :: ZS_step = dble(0.1e-10)
integer :: ZS_npts = 30

!--- Parameters for SteadyDissi ---!
!double precision :: SD_step = dble(0.2e-10)

contains


subroutine ImagerInit()
  use consts
  implicit none
  character*80 :: chr1,chr2
  
  !normalize the slow scan directions
  SL_fdir = SL_fdir / sqrt( SL_fdir(1)**2 + SL_fdir(2)**2 )
  SL_sdir = SL_sdir / sqrt( SL_sdir(1)**2 + SL_sdir(2)**2 )

  write(*,*) 'fast scan direction', SL_fdir
  write(*,*) 'slow scan direction', SL_sdir

  write(chr1,'(F15.8)')imgSide/nano
  write(chr2,'(I8)')imgPix
  
  write(*,'(A)')'Imaging area: '//trim(adjustl(chr1))//' x '//trim(adjustl(chr1))//' nm^2'
  write(chr1,'(I8)')imgPixs
  write(*,'(A)')'scan resolution: '//trim(adjustl(chr1))//' x '//trim(adjustl(chr2))//' px^2.'
  write(*,'(A)')
  
end subroutine ImagerInit



subroutine ReadScript()
	implicit none
	integer :: r=909,ioer
	character*300 :: line
	
	open(file=ScriptFile, unit=r)
	ScriptLines=0
	do
		read(r,'(A)',IOSTAT=ioer)line
		if(ioer/=0) exit
		if(line(1:1) /= '#') ScriptLines=ScriptLines+1
		
	enddo
	
	
	
	close(r)
	
end subroutine ReadScript


!*** ****** Theoretical dF map ****** ***!
!
subroutine TheorImage()
  use consts
  use cantilever
  use AGC
  use ADC
  implicit none
  
  double precision :: scanX=0.0d0,scanY=0.0d0,scanZ
  double precision :: dfZ,dfYt,dfY
  integer :: i,j,k
  
  open(file='theorimage.out',unit=888)
  
  scanZ=zLand
  !loop ver all points
  do i=0,imgPixs
    scanX=i*imgSide/imgPixs
    do j=0, imgPix
      scanY=j*imgSide/imgPix
      
      call CalcDFz(scanX,scanY,scanZ,ASet_z,dfZ)
      dfZ=-dfZ*F1/(dPi*ASet_z*k_z)
      
      call CalcDFy(scanX,scanY,6.5d0*ang,ASet_y,dfY)
      dfY=-dfY*F2/(dPi*ASet_y*k_y)
      
      write(888,*)scanX,scanY,dfZ,dfY
      
    enddo
    write(888,*)
  enddo
  close(888)
  
  
end subroutine TheorImage

!*** ****** Theoretical Frequency Shift ****** ***!
!
subroutine CalcDFz(x,y,z,Amp,dfZ)
  use consts
  use forcer
  implicit none
  double precision, intent(in) :: x,y,z,Amp
  double precision, intent(out):: dfZ
  double precision :: teta,tetao,f(3),fo(3),pos(3),vdW
  double precision :: dx
  
  dx=dPi/10000.0d0
  pos(1)=x
  pos(2)=y
  pos(3)=z
  dfZ=0.0d0
  
  !set the old position value for trapezoid integration
  tetao=0.0d0
  call CalcForce(pos,fo)
  call CalcVanDerWaals(pos,vdW)
  fo(3)=fo(3)+vdW
  
  do teta=dx, dPi, dx
  
    pos(3)=z+Amp*Cos(teta)
    call CalcForce(pos,F)
    call CalcVanDerWaals(pos,vdW)
    f(3)=f(3)+vdW
    
    dfZ=dfZ+0.50d0*(fo(3)*Cos(tetao)+f(3)*Cos(teta))*dx
    
    tetao=teta
    fo=f
  enddo
  
  
end subroutine CalcDFz

subroutine CalcDFy(x,y,z,Amp,dfY)
  use consts
  use forcer
  implicit none
  double precision, intent(in) :: x,y,z,Amp
  double precision, intent(out):: dfY
  double precision :: teta,tetao,f(3),fo(3),pos(3)
  double precision :: dx
  
  dx=dPi/10000.0d0
  pos(1)=x
  pos(2)=y
  pos(3)=z
  dfY=0.0d0
  
  !set the old position value for trapezoid integration
  tetao=0.0d0
  call CalcForce(pos,fo)
  
  do teta=dx, dPi, dx
  
    pos(2)=y+Amp*Cos(teta)
    call CalcForce(pos,f)
    
    dfY=dfY+0.50d0*(fo(2)*Cos(tetao)+f(2)*Cos(teta))*dx
    
    tetao=teta
    fo=f
    
  enddo
  
  
end subroutine CalcDFy


end module Imager


