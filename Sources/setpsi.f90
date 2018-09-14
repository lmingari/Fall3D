  subroutine setpsi(psi,sphe,diam,modv,nc)
  !*********************************************************************
  !*
  !*      Calculates the particle shape factor psi depending on the
  !*      velocity model
  !*
  !*      modv = 1   ARASTOPOUR. psi = 1 (not used)
  !*      modv = 2   GANSER      psi = sphericity
  !*      modv = 3   WILSON      psi = (b+c)/2a    a>b>c semi-axes
  !*      modv = 4   DELLINO     psi = sphericity/circularity
  !*
  !**********************************************************************
  use KindType
  implicit none
  integer(ip) :: modv,nc
  real   (rp) :: psi(nc),sphe(nc),diam(nc)
  !
  integer(ip) :: ic
  real   (rp) :: pi,gama,circula
  !
  !***    Initializations
  !
  pi = 4.0_rp*atan(1.0_rp)
  !
  !***    Computes psi
  !
  if(modv == 1) then          ! Arastopour
     !
     psi(1:nc) = 1.0_rp
     !
  else if(modv == 2) then     ! Ganser
     !
     psi(1:nc) = sphe(1:nc)
     !
  else if(modv == 3) then     ! Wilson
     !
     do ic = 1,nc
        call get_gama(diam(ic),sphe(ic),gama)  ! Get a/c
        if(gama >= 1.0_rp) then                ! oblate
           psi(ic) = 0.5_rp*(1.0_rp+1.0_rp/gama)
        else                                ! prolate
           psi(ic) = gama
        end if
     end do
     !
  else if(modv == 4)  then    ! Dellino
     !
     do ic = 1,nc
        call get_gama(diam(ic),sphe(ic),gama)         ! Get a/c
        if(gama >= 1.0_rp) then                       ! oblate
           circula = 1.0_rp
           psi(ic) = sphe(ic)/circula
        else                                    ! prolate
           circula = sqrt(gama)                 ! Riley 1941
           psi(ic) = sphe(ic)/circula
        end if
     end do
  end if
  !
  return
end subroutine setpsi
!
!
!
subroutine get_gama(diam,sphe,gama)
  !**************************************************************************
  !*
  !*     Gets gama = a/c
  !*
  !*     NOTE: In all cases it is assumed that particles fall as
  !*     prolate ellipsoids
  !*
  !*            a = b < c  prolate   (gama < 1)
  !*
  !*     The inversion of the area of the ellipsoid is done numerically.
  !*     Area given by:
  !*
  !*     A = 2*pi*(a**2 + c**2*e/tan(e) )   e = acos(gama)
  !*     d = 2*c*gama**(2/3)               (prolate)
  !*
  !*     NOTE: particle diameter is multiplied by a factor. It does not affect
  !*           results (a/c) and is done to facilitate convergence and prevent
  !*           propagation of rounding errors (e.g. for micron size particles
  !*           diam of the order 1d-6 rised to 2 or 3)
  !*
  !***************************************************************************
  use KindType
  implicit none
  real(rp)    :: diam,sphe,gama
  !
  integer(ip) :: iiter,niter
  real   (rp) :: d,pi,gmin,gmax,Ao,toler,e
  real   (rp) :: Vp,Ap
  !
  !***   Initializations
  !
  d     = diam*1d3         ! see NOTE
  niter = 1000
  toler = 1d-8
  gmin  = 1d-3
  gmax  = 1.0_rp
  !
  !***   Volume and area
  !
  pi = 4.0_rp*atan(1.0_rp)
  Vp = 4.0_rp*pi*((0.5_rp*d)**3)/3.0_rp
  Ap = (pi**(1.0_rp/3.0_rp))*((6.0_rp*Vp)**(2.0_rp/3.0_rp))/sphe
  !
  !***   Iterates
  !
  do iiter = 1,niter
     gama = 0.5_rp*(gmin+gmax)
     e    = acos(gama)
     Ao   = 0.5_rp*pi*d*d*(gama**(-4.0_rp/3.0_rp))*(gama*gama + (e/tan(e)))
     if(Ao < Ap) then
        gmax = gama
     else
        gmin = gama
     end if
     if((iiter > 1).and.(abs(Ao-Ap) < toler)) goto 10
  end do
  call runend('Subroutine get_gama: convergence not achieved')
  !
  !***  convergence
  !
10 return
end subroutine get_gama
