Module Cantilever
  Use Consts
  ! Cantilever Parameters
  double precision :: k_z=30.0d0
  double precision :: mass=0.0d0
  Double Precision :: W1=150000.0d0*dPi, W1e=150000.0d0*dPi, Q1=30000.0d0
  double precision :: F1=150000.0d0
  double precision :: Gamma1=0.0d0
  
  double precision :: k_y=1250.0d0
  double precision :: massTR=0.0d0
  double precision :: F2=1500000.0d0
  Double Precision :: W2=1500000.0d0*dPi, W2e=1500500.0d0*dPi, Q2=150000.0d0
  double precision :: Gamma2
  
  Double Precision :: A0z=0.0d0*Nano,P0z=0.0d0*Nano
  Double Precision :: A0y=Ang, P0y=0.0d0
  
  ! Positioning stuff
  double precision :: HolderPosition(3)=0.0d0 !cantilever holder position
  Double Precision :: X(3)=0.0d0              !tip position
  double precision :: Xo(3),Xoo(3)            !tip old positions
  double precision :: V(3)=0.0d0              !speed
  double precision :: A(3)=0.0d0              !acceleration
  double precision :: Noise(3)=0.0d0          !noise on position
  double precision :: Ftipsurf(3)=0.0d0,Ftipsurfo(3)=0.0d0,Ftipsurfoo(3)=0.0d0
  double precision :: fNoise(3)=0.0d0         !noise on force
  
  
  ! Tip Parameters
  double precision :: TipHamaker=0.0d0
  double precision :: TipRadius=0.0d0
  double precision :: TipAngle=0.0d0
  double precision :: sing,tang,cosg,cos2g
  
  Logical :: IntON_z=.true.    !switch for the tip-surface interaction (z)
  Logical :: IntON_y=.true.    !switch for the tip-surface interaction (y)  
 
	Logical :: Oscill_z = .true.
	Logical :: Oscill_y = .true.
 
  !backup variables
  double precision :: bHP(3),bX(3),bXo(3),bXoo(3),bV(3),bA(3)
 
contains

!*** cantilever initialization routine ***
subroutine CLVInit()
  use consts
  implicit none
  
  !compute the angular frequencies
  W1=dPi*F1
  W2=dPi*F2

  if(mass==0.0d0) then
    mass=k_z/(w1**2.0d0)      !calculate mass if not supplied
    write(*,'(A ES14.8 A)') "Cantilever's mass was estimated: ",mass, " Kg"
  endif
  
  massTR=k_y/(W2**2.0d0)
  write(*,'(A ES14.8 A)')"Effective torsional mass was estimated: ",massTR, " Kg"
  
  !van der waals stuff
  sing=Sin(TipAngle*g2r)
  tang=Tan(TipAngle*g2r)
  cosg=Cos(TipAngle*g2r)
  cos2g=Cos(TipAngle*g2r*2.0d0)
  
  write(*,*)
  
end subroutine CLVInit
!*****************************************

End Module Cantilever
