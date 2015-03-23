module Input

logical :: DBG=.false.
logical :: ForceSpectr
integer :: WorkMode = 1          !working operation (1->interactive , 2->single line, 3->image)

integer :: ifile=959
character*40 :: filein

character*80 :: Mode, vAFMfile

integer :: FSpectrAxis=0
integer :: FS_file=23
double precision :: FS_Ipos=0.0d0
double precision :: FS_Fpos=0.0d0
double precision :: FS_dpos=1.0d0
double precision :: FS_1pos=0.0d0
double precision :: FS_2pos=0.0d0
double precision :: FS_3pos=0.0d0


contains


subroutine LoadALL()
  implicit none
  
  call getarg(1,filein)
  
  write(*,*)'Loading parameters from ',trim(adjustl(filein)), ' ...'
  
  Open(file=filein, unit=ifile)
  
  call Load_SYS()
  
  call load_force()
  call Load_Imaging()
  call Load_ForceSpect()
  call Load_Scanliner()

  !if(WorkMode == 2) call Load_Scanliner()
  if(WorkMode == 4) then
  	!call Load_Scanliner()
  	call Load_ZSprecter()
  endif
  
  if(WorkMode == 8) call Load_ZSprecter()
  
  close(ifile)
  
  Open(file=trim(adjustl(vAFMfile)), unit=ifile)
  
  call Load_CLV()
  call Load_AGC()
  call Load_PLLz()
  call Load_PLLy()
  call Load_ADC()
  close(ifile)
  
  write(*,*)'Input file read! All parameters set!'

end subroutine LoadALL



subroutine Load_Force()
  use forcer
  implicit none
  integer :: idx
  character*256 :: line,tmpstr,strval
	
	
	
	!read the force field type
  rewind(ifile)
  call find_str(ifile,'forcetype',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,strval
    if(trim(adjustl(strval)) == 'single') then
    	write(*,*)'FORCES: the force field type is single point.'
    	ForceType = 2
    else
    	write(*,*)'FORCES: the force field type is volumetric.'
    	ForceType = 1
    endif
  endif
  
  !read the forcefield filename
  rewind(ifile)  
  call find_str(ifile,'forcefield',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,ForceFile
  endif
	!read the force grid filename
  rewind(ifile)
  call find_str(ifile,'forcegrids',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,GridFile
  endif
	!read the dissipation map if any
  rewind(ifile)
  NoDissiEver = .true.
  call find_str(ifile,'dissifield',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,DissiFile
    NoDissiEver = .false.
  else
  	write(*,*)'Dissipation Field was not given!'
  	NoDissiEver = .true.
  endif
  !read vAFM parameters filename
  rewind(ifile)  
  call find_str(ifile,'vmachine',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,vAFMfile
  endif
  
  rewind(ifile)  
  call find_str(ifile,'ignoredissi',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,IgnoreDissip
  endif
  
  rewind(ifile)  
  call find_str(ifile,'s.temperature',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,Temp
  endif
  
  rewind(ifile)  
  call find_str(ifile,'s.forcenoise',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,ForceNoise
  endif

	write(*,*)'FORCES: reading forcefield from  : '//trim(adjustl(ForceFile))
	write(*,*)'FORCES: reading grid shape from  : '//trim(adjustl(GridFile))
	if(.not.NoDissiEver) then
		write(*,*)'FORCES: reading dissipation from : '//trim(adjustl(DissiFile))
	endif
	write(*,*)'INFO: virtual machine parameters: '//trim(adjustl(vAFMfile))
	
	!stop
end subroutine Load_Force 


subroutine Load_SYS()
  use Timer
  use Imager
  use Output
  implicit none
  integer :: idx
  character*256 :: line,tmpstr,strval
  
!   rewind(ifile)
!   call find_str(ifile,'s.SYS.DB',idx,line)
!   if(line/='000') then
!     read(line,'(A8, L5)')tmpstr,DBG
!     write(*,*)'debug is',DBG
!   endif  
	DBG = .false.
	
  rewind(ifile)
  call find_str(ifile,'dump.ini',idx,line)
  if(line/='000') then
    read(line,'(A9, L5)')tmpstr,DumpPhase(1)
    write(*,*)'Write Initialization: ',DumpPhase(1)
  endif
  rewind(ifile)
  call find_str(ifile,'dump.rel',idx,line)
  if(line/='000') then
    read(line,'(A10, L5)')tmpstr,DumpPhase(2)
    write(*,*)'Write Relaxation: ',DumpPhase(2)
  endif
  rewind(ifile)
  call find_str(ifile,'dump.lnd',idx,line)
  if(line/='000') then
    read(line,'(A9, L5)')tmpstr,DumpPhase(3)
    write(*,*)'Write Approach: ',DumpPhase(3)
  endif
  rewind(ifile)
  call find_str(ifile,'dump.fvs',idx,line)
  if(line/='000') then
    read(line,'(A10, L5)')tmpstr,DumpPhase(4)
    write(*,*)'Write Free Oscillations: ',DumpPhase(4)
  endif
  
  !read the time steps (normal operation and KMC region)
  rewind(ifile)
  call find_str(ifile,'s.timestep',idx,line)
  if(line/='000') then
    read(line,'(A10, F15.8)')tmpstr,tstep
    dtNorm=tstep
  endif
  rewind(ifile)
  call find_str(ifile,'s.ATtimestep',idx,line)
  if(line/='000') then
    read(line,'(A12, F15.8)')tmpstr,dtATzone
  endif

  rewind(ifile)
  call find_str(ifile,'writepts',idx,line)
  if(line/='000') then
    read(line,'(A8, I8)')tmpstr,PointBuffer
  endif
  
  rewind(ifile)
  call find_str(ifile,'workmode',idx,line)
  if(line/='000') then
    read(line,'(A8, A)')tmpstr,strval
    strval=trim(adjustl(strval))
    write(*,*)
    write(*,'(A)', advance='no')'Working mode: '
    select case(strval)
    case('debug')
    	WorkMode=0
    	write(*,'(A)')'Debug'
    case('interactive')
      WorkMode=1
      write(*,'(A)')'Interactive operation'
    case('image')
      WorkMode=3
      write(*,'(A)')'Direct imaging'
    case('scanline')
      WorkMode=2
      write(*,'(A)')'Single scanline'
    case('zspectra')
    	WorkMode=4
    	write(*,'(A)')'Z Spectroscopy'
    case('point')
    	WorkMode=5
    	write(*,'(A)')'Single point'
    case('forcecurve')
    	WorkMode=6
    	write(*,'(A)')'Force curve'
  case('gissimage')
     WorkMode=7
     write(*,'(A)')'Theoretical image'
  case('steadycalc')
     WorkMode=8
     write(*,'(A)')'Steady state calculation'
  case default
      write(*,'(A)')'No valid operation mode selected, using interactive!'
      stop
    end select
    write(*,*)
  endif
  

end subroutine Load_SYS


subroutine Load_Scanliner()
	use imager
	implicit none
  integer :: idx
  character*256 :: line,tmpstr,strval

  rewind(ifile)
  call find_str(ifile,'sl.length',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,SL_Length
  endif
  rewind(ifile)
  call find_str(ifile,'sl.points',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,SL_npts
  endif
  rewind(ifile)
  call find_str(ifile,'sl.slowdir',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,SL_sdir(1),SL_sdir(2)
  endif
  rewind(ifile)
  call find_str(ifile,'sl.fastdir',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,SL_fdir(1),SL_fdir(2)
    write(*,*)'asddir',SL_fdir
  endif
 
end subroutine Load_Scanliner


subroutine Load_Debug()
	use imager
	implicit none
	
	
	
	
end subroutine Load_Debug



subroutine Load_ZSprecter()
	use imager
	implicit none
  integer :: idx
  character*256 :: line,tmpstr,strval
  
  rewind(ifile)
  call find_str(ifile,'zs.zstep',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,ZS_step
  endif
  rewind(ifile)
  call find_str(ifile,'zs.points',idx,line)
  if(line/='000') then
    read(line,*)tmpstr,ZS_npts
  endif
  
end subroutine Load_ZSprecter

subroutine Load_Imaging()
  use Imager
  implicit none
  character*256 :: line,tmpstr,strval
  integer :: idx

  rewind(ifile)
  call find_str(ifile,'image.ps',idx,line)
  if(line/='000') then
    read(line,'(A8 I8)')tmpstr,imgPixs
  endif
  rewind(ifile)
  call find_str(ifile,'image.pf',idx,line)
  if(line/='000') then
    read(line,'(A8 I8)')tmpstr,imgPix
  endif
  
  rewind(ifile)
  call find_str(ifile,'image.side',idx,line)
  if(line/='000') then
    read(line,'(A10 F15.8)')tmpstr,imgSide
  endif  
  rewind(ifile)
  call find_str(ifile,'image.mode',idx,line)
  if(line/='000') then
    read(line,'(A10 A10)')tmpstr,ImageMode
  endif
  
  rewind(ifile)
  call find_str(ifile,'image.pattern',idx,line)
  if(line/='000') then
    read(line,'(A13 A)')tmpstr,strval
    select case(trim(adjustl(strval)))
    case('zigozago')
      ImagePattern=2
    case('lineline')
      ImagePattern=1
    case('forwback')
      ImagePattern=3
    end select
  endif  
  
  rewind(ifile)
  call find_str(ifile,'image.speed',idx,line)
  if(line/='000') then
    read(line,'(A11 F15.8)')tmpstr,ScanSpeed
  endif
  rewind(ifile)
  call find_str(ifile,'image.rspeed',idx,line)
  if(line/='000') then
    read(line,'(A12 F15.8)')tmpstr,RepoSpeed
  endif
  
!!$  rewind(ifile)
!!$  call find_str(ifile,'image.fastdir',idx,line)
!!$  if(line/='000') then
!!$    read(line,'(A13 A)')tmpstr,strval
!!$    if(trim(adjustl(strval))=='y' ) ScanAxis=2
!!$    if(trim(adjustl(strval))=='x' ) ScanAxis=1
!!$    
!!$    !write(*,*)ScanAxis
!!$    !stop
!!$  endif
  
  
  !stop

end subroutine Load_Imaging


subroutine Load_ForceSpect()
  implicit none
  integer :: idx
  character*256 :: line,tmpstr
  
  !*** read force spectroscopy parameters
  rewind(ifile)
  call find_str(ifile,'FSP.AXIS',idx,line)
  if(line/='000') then
    read(line,'(A8 I5)')tmpstr,FSpectrAxis
  endif
  rewind(ifile)
  call find_str(ifile,'FSP.FILE',idx,line)
  if(line/='000') then
    read(line,'(A8 I5)')tmpstr,FS_file
  endif
  rewind(ifile)
  call find_str(ifile,'FSP.Ipos',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,FS_Ipos
  endif
  rewind(ifile)
  call find_str(ifile,'FSP.dpos',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,FS_dpos
  endif  
  rewind(ifile)
  call find_str(ifile,'FSP.Fpos',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,FS_Fpos
  endif
  rewind(ifile)
  call find_str(ifile,'FSP.1pos',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,FS_1pos
  endif
  rewind(ifile)
  call find_str(ifile,'FSP.2pos',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,FS_2pos
  endif
  rewind(ifile)
  call find_str(ifile,'FSP.3pos',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,FS_3pos
  endif

end subroutine Load_ForceSpect


subroutine Load_ADC()
  use ADC
  implicit none
  character*256 :: line,tmpstr
  integer :: idx

  rewind(ifile)
  call find_str(ifile,'z.ADC.Zi',idx,line)
  if(line/='000') then
    read(line,'(A8,F15.8)')tmpstr,zStart
  endif  
  rewind(ifile)
  call find_str(ifile,'z.ADC.Zf',idx,line)
  if(line/='000') then
    read(line,'(A8,F18.12)')tmpstr,zLand
  endif  
  rewind(ifile)
  call find_str(ifile,'z.ADC.LandTime',idx,line)
  if(line/='000') then
    read(line,'(A14,F15.8)')tmpstr,tLand
  endif  
  rewind(ifile)
  call find_str(ifile,'z.ADC.zLim',idx,line)
  if(line/='000') then
    read(line,'(A10,F15.8)')tmpstr,zCrashLimit
    write(*,*)'INFO: Crash limit set to:',zCrashLimit
  endif  

  rewind(ifile)
  call find_str(ifile,'z.ADC.dF',idx,line)
  if(line/='000') then
    read(line,'(A8,F15.8)')tmpstr,ADCSet_z
  endif
  rewind(ifile)
  call find_str(ifile,'z.ADC.KP',idx,line)
  if(line/='000') then
    read(line,'(A8,F15.8)')tmpstr,AKP_z
  endif  
  rewind(ifile)
  call find_str(ifile,'z.ADC.KI',idx,line)
  if(line/='000') then
    read(line,'(A8,F15.8)')tmpstr,AKI_z
  endif    

  !*** approach relaxation
  rewind(ifile)
  call find_str(ifile,'z.ADC.Re',idx,line)
  if(line/='000') then
    read(line,'(A8,F15.8)')tmpstr,ARelaxTol_z
  endif
  rewind(ifile)
  call find_str(ifile,'z.ADC.Rt',idx,line)
  if(line/='000') then
    read(line,'(A8,F15.8)')tmpstr,ARelaxTimeSpan_z
  endif 



end subroutine Load_ADC



subroutine Load_CLV()
  use cantilever
  implicit none
  character*256 :: line,tmpstr
  integer :: idx
  double precision :: dvl
  
  
  rewind(ifile)
  call find_str(ifile,'CLV.mass',idx,line)
  if(line/='000') then
    read(line,'(A8,F15.8)')tmpstr,mass
  endif 
  
  
  !*** van der Waals parameters
  rewind(ifile)
  call find_str(ifile,'TIP.hamk',idx,line)
  if(line/='000') then
    read(line,'(A8,F15.8)')tmpstr,TipHamaker
  endif 
  rewind(ifile)
  call find_str(ifile,'TIP.radi',idx,line)
  if(line/='000') then
    read(line,'(A8,F15.8)')tmpstr,TipRadius
  endif 
  rewind(ifile)
  call find_str(ifile,'TIP.angl',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,TipAngle
  endif
  
  rewind(ifile)
  call find_str(ifile,'z.CLV.F0',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,F1
  endif
  rewind(ifile)
  call find_str(ifile,'z.CLV.Kf',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,k_z
  endif
  rewind(ifile)
  call find_str(ifile,'z.CLV.Qf',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,Q1
  endif
  
  rewind(ifile)
  call find_str(ifile,'y.CLV.F0',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,F2
  endif
  rewind(ifile)
  call find_str(ifile,'y.CLV.Kf',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,k_y
  endif
  rewind(ifile)
  call find_str(ifile,'y.CLV.Qf',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,Q2
  endif
  
  rewind(ifile)
  call find_str(ifile,'holder.x',idx,line)
  if(line/='000') then
    read(line,'(A8, F15.8)')tmpstr,dvl
    HolderPosition(1) = dvl
  endif
  rewind(ifile)
  call find_str(ifile,'holder.y',idx,line)
  if(line/='000') then
    read(line,'(A8, F15.8)')tmpstr,dvl
    HolderPosition(2) = dvl
  endif
  !write(*,*)HolderPosition(:)
  !stop
  
  !read what oscillation to perform
  rewind(ifile)
  call find_str(ifile,'z.oscill',idx,line)
  if(line/='000') then
    read(line,'(A8, L5)')tmpstr, Oscill_z
  endif
  rewind(ifile)
  call find_str(ifile,'y.oscill',idx,line)
  if(line/='000') then
    read(line,'(A8, L5)')tmpstr, Oscill_y
  endif
  write(*,*) 'Vertical Oscillation is ', Oscill_z
  write(*,*) 'Lateral  Oscillation is ', Oscill_y
  
end subroutine Load_CLV

subroutine Load_AGC()
  use consts
  use AGC
  implicit none
  character*256 :: line,tmpstr
  integer :: idx
  
  rewind(ifile)
  call find_str(ifile,'z.AGC.setpoint',idx,line)
  if(line/='000') then
    read(line,'(A14, F15.8)')tmpstr,Aset_z
  endif
  rewind(ifile)
  call find_str(ifile,'z.AGC.offset',idx,line)
  if(line/='000') then
    read(line,'(A12,F15.8)')tmpstr,FBo_z
  endif
  rewind(ifile)
  call find_str(ifile,'z.AGC.KP',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,KP_z
  endif  
  rewind(ifile)
  call find_str(ifile,'z.AGC.KI',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,KI_z
  endif    
  
  rewind(ifile)
  call find_str(ifile,'z.AGC.nf',idx,line)
  if(line/='000') then
    read(line,'(A8, I8)')tmpstr,AFilterOrder_z
  endif

  rewind(ifile)
  call find_str(ifile,'z.AGC.fc',idx,line)
  if(line/='000') then
    read(line,'(A8, F15.8)')tmpstr,ARC_z
    ARC_z=1.0d0/(dPi*ARC_z)
  endif

  
  rewind(ifile)
  call find_str(ifile,'z.AGC.Re',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,RelaxTol_z
  endif
  rewind(ifile)
  call find_str(ifile,'y.AGC.Re',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,RelaxTol_y
  endif
  rewind(ifile)
  call find_str(ifile,'z.AGC.Rt',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,RelaxTimeSpan_z
  endif  
  
  
  rewind(ifile)
  call find_str(ifile,'y.AGC.setpoint',idx,line)
  if(line/='000') then
    read(line,'(A14, F15.8)')tmpstr,Aset_y
  endif
  rewind(ifile)
  call find_str(ifile,'y.AGC.offset',idx,line)
  if(line/='000') then
    read(line,'(A12, F15.8)')tmpstr,FBo_y
  endif
  rewind(ifile)
  call find_str(ifile,'y.AGC.KP',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,KP_y
  endif  
  rewind(ifile)
  call find_str(ifile,'y.AGC.KI',idx,line)
  if(line/='000') then
    read(line,'(A8 F15.8)')tmpstr,KI_y
  endif
  
  rewind(ifile)
  call find_str(ifile,'y.AGC.nf',idx,line)
  if(line/='000') then
    read(line,'(A8, I8)')tmpstr,AFilterOrder_y
  endif
  rewind(ifile)
  call find_str(ifile,'y.AGC.fc',idx,line)
  if(line/='000') then
    read(line,'(A8, F15.8)')tmpstr,ARC_y
    ARC_y=1.0d0/(dPi*ARC_y)
  endif
  
  
end subroutine Load_AGC

subroutine Load_PLLz()
  use PLL
  use consts
  implicit none
  character*256 :: line,tmpstr
  integer :: idx
  
  
  rewind(ifile)
  call find_str(ifile,'nfilterz',idx,line)
  if(line/='000') then
    read(line,'(A8, I8)')tmpstr,FilterOrder_z
  endif  
  
  rewind(ifile)
  call find_str(ifile,'z.PLL.fc',idx,line)
  if(line/='000') then
    read(line,'(A8, F15.8)')tmpstr,RC_z
  endif
  RC_z=1.0d0/(RC_z*dPi)


end subroutine Load_PLLz

subroutine Load_PLLy()
  use PLL
  use consts
  implicit none
  character*256 :: line,tmpstr
  integer :: idx
  
  
  rewind(ifile)
  call find_str(ifile,'nfiltery',idx,line)
  if(line/='000') then
    read(line,'(A8, I8)')tmpstr,FilterOrder_y
  endif  
  
  rewind(ifile)
  call find_str(ifile,'y.PLL.fc',idx,line)
  if(line/='000') then
    read(line,'(A8, F15.8)')tmpstr,RC_y
  endif
  RC_y=1.0d0/(RC_y*dPi)


end subroutine Load_PLLy


subroutine find_strl(fnum,str,ln,idx,line)
	implicit none
	integer, intent(in) :: fnum,ln
	integer, intent(out) :: idx
	character*20, intent(in) :: str
	character*256, intent(out) :: line
	character*256 :: newline
	character*256 :: newstr
	integer :: ioer
	
	!newline = trim(adjustl(str))
	!newstr = newline(1:ln)
	
	do
		read(fnum,'(A)',iostat=ioer) line
		if(ioer /= 0) exit
		newline = trim(adjustl(line))
		write(*,*)newline(1:ln)//'---'//str(1:ln)//'***'
		if(newline(1:ln) == trim(adjustl(str))) return
		if ((line(1:1)/='#')) return
		
	enddo
	line = '000'
	
end subroutine find_strl

subroutine find_str(fnum,str,idx,line)
  implicit none
  integer, intent(in) :: fnum         !file number
  character*8, intent(in) :: str      !string to find (itz always 8 characters!!!
  integer, intent(out) :: idx         !character index in 'line'
  character*256, intent(out):: line   !Line containing 'str'

  write(*,*)str

10  line = ""
    read(fnum, '(A)', end=20) line                !read a line
    idx = index(trim(adjustl(line)), trim(adjustl(str)))                       !look for the string
    
!    if((str(1:2)=='sl')) then
!       !idx = index('start.x asd','start.x')
!       write(*,*)trim(adjustl(line)),idx,trim(adjustl(str))
!       !stop
!    endif
    if ((idx .ne. 0) .and. (line(1:1)/='#')) then !if the string is present...
      !if(DBG)print *, line(idx:)                  !print the line if in debug mode
      return                                      !exit the subroutine
    endif                                         !***
    goto 10                             !loop!
    
20  line='000'                    !the end of the file was reached!
21  return

end subroutine find_str


end module Input
