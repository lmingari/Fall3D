!***************************************************************
!*
!*		Module for master operations
!*
!***************************************************************
MODULE Master
  use KindType
  IMPLICIT NONE
  SAVE
  !
  !***  Version of input file
  !
  real(rp),parameter :: version = 7.0
  !
  !***  Granulometry and magma-type density
  !
  integer(ip) ::  nc,ng
  real   (rp) ::  fimin,fimax,fimean(2),fidisp(2),rhomin,rhomax,sphemin,sphemax
  !
  !***  Type of distribution
  !
  character(len=25) :: TYPE_DIST
  integer(ip)       :: mdist = 0
  !
  !***  Memory storage
  !
  integer(ip) :: memo = 0
  !
  !***  Arrays
  !
  real(rp), allocatable :: fc  (:)      ! fc(nc)
  real(rp), allocatable :: rhop(:)      ! rhop(nc)
  real(rp), allocatable :: diam(:)      ! diam(nc)
  real(rp), allocatable :: sphe(:)      ! sphe(nc)
  !
CONTAINS
  !
  subroutine getmem
    !**********************************************************************
    !*
    !*  Allocates memory for arrays related to nc
    !*
    !**********************************************************************
    implicit none
    !
    allocate(fc(nc))
    fc = 0.
    memo = memo + rp*nc
    !
    allocate(rhop(nc))
    rhop = 0.
    memo = memo + rp*nc
    !
    allocate(diam(nc))
    diam = 0.
    memo = memo + rp*nc
    !
    allocate(sphe(nc))
    sphe = 0.
    memo = memo + rp*nc
    !
  end subroutine getmem
  !
END MODULE Master
