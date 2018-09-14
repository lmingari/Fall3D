  subroutine solvesuzuki
  !********************************************************
  !*
  !*    Suzuki source solution
  !*     x : LON
  !*     y : LAT
  !*     z : height a.s.l
  !*
  !********************************************************
  use KindType
  use Master
  implicit none
  !
  integer(ip) :: is,ic
  real(rp)    :: sum,deltaz,z
  real(rp), allocatable :: S(:)
  !
  !***  Horizontal source position
  !
  Xplum(1:ns) = X0
  Yplum(1:ns) = Y0
  !
  !***  Vertical position and total mass according to Suzuki distribution
  !
  allocate(S(ns))
  sum = 0.0
  deltaz = Hplume/ns   ! from the botttom to Hplum
  do is = 1,ns
     z = is*deltaz
     S(is) = ((1.0-z/Hplume)*exp(Asuzu*(z/Hplume-1.0)))**Lsuzu
     Zplum(is)= Z0 + dZ0 + z
     sum = sum + S(is)
  end do
  !
  !*** Normalization to MFR (SUMz=MFR)
  !
  S(1:ns) = M0*S(1:ns)/sum
  !
  !***  Mass distribution
  !
  do is = 1,ns
     do ic = 1,nc
        Mplum(ic,is) = fc(ic)*S(is)
     end do
  end do
  !
  deallocate(S)
  return
  end subroutine solvesuzuki
