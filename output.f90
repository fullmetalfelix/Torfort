!*** ****** Output Module ****** ***!
module OutPut

integer :: PointBuffer=10000
integer :: Phase = 0
integer :: oFEVS = 0
integer :: LogFile = 321
logical :: DumpPhase(5)=.false.


contains

subroutine InitOutput()



end subroutine InitOutput


subroutine DeInitOutput()
  
  implicit none
  
  
  close(LogFile)

end subroutine DeInitOutput

subroutine WriteAll()
  use Consts
  use Cantilever
  use Timer
  use PLL
  use ADC
  use AGC
  use NURBSDissiField
  implicit none
	
	
  write(1000+Phase,*) t,HolderPosition(3)+x(3),SinPLL_z*ASet_z,At_z,FB_z,FBp_z,FBi_z,DeltaF_z, &
    ADC_z,ADCp_z,ADCi_z, &
    x(2),SinPLL_y*ASet_y,At_y,FB_y,FBp_y,FBi_y,DeltaF_y,Dissipation_z,DissipationI_z
	
	!write(2000+Phase,*) t,phd_z,p1_z,p2_z,ut_z,pref_z,phd_y,p1_y,p2_y,ut_y,pref_y
	

end subroutine WriteAll


subroutine WriteAll_f(fileptr)
  use Consts
  use Cantilever
  use Timer
  use PLL
  use ADC
  use AGC
  use NURBSDissiField
  implicit none
  integer, intent(in) :: fileptr
	
  write(fileptr,*) t,HolderPosition(3)+x(3),SinPLL_z*ASet_z,At_z,FB_z,FBp_z,FBi_z,DeltaF_z, &
    ADC_z,ADCp_z,ADCi_z, &
    x(2),SinPLL_y*ASet_y,At_y,FB_y,FBp_y,FBi_y,DeltaF_y,Dissipation_z,DissipationI_z
	
	!write(2000+Phase,*) t,phd_z,p1_z,p2_z,ut_z,pref_z,phd_y,p1_y,p2_y,ut_y,pref_y
	

end subroutine WriteAll_f

subroutine WriteScanLine()
  use Consts
  use Cantilever
  use Timer
  use PLL
  use ADC
  use AGC
  use NURBSDissiField
  implicit none

  
!   write(1000+Phase,*) t,x(3),SinPLL_z*ASet_z,At_z,FB_z,FBp_z,FBi_z,DeltaF_z, &
!     ADC_z,ADCp_z,ADCi_z, &
!     x(2),SinPLL_y*ASet_y,At_y,FB_y,FBp_y,FBi_y,DeltaF_y,Dissipation_z, &
!     HolderPosition(1), HolderPosition(2)
  
	write(1000+Phase,*) t,HolderPosition(1), HolderPosition(2), &
		At_z,FB_z,DeltaF_z,Dissipation_z,DissipationI_z, &
		At_y,FB_y,DeltaF_y, Dissipation_y,DissipationI_y,ADC_z

end subroutine WriteScanLine

subroutine WriteForceCurve_Static(fileptr)
	use Cantilever
  use AGC
  use PLL
  implicit none
  integer, intent(in) :: fileptr
  
  write(fileptr,*)HolderPosition(3), At_z, DeltaF_z
  
  
end subroutine WriteForceCurve_Static


subroutine LogRelax()
  use Cantilever
  use ADC
  use AGC
  use PLL
  implicit none
  
  !write(LogFile,*) 
  !write(LogFile,*) 'Relaxing the oscillation...'




end subroutine LogRelax



end module OutPut