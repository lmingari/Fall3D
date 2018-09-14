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
  integer(ip):: modv,nc
  real(rp)   :: psi(nc),sphe(nc),diam(nc)
  !
  integer(ip):: ic
  real(rp)   :: pi,gama,circula
  !
  !***  Initializations
  !
  pi = 4.0*atan(1.0)
  !
  !***  Computes psi
  !
  if(modv.eq.1) then          ! Arastopour
     !
     do ic=1,nc
        psi(ic) = 1.0
     end do
     !
  else if(modv.eq.2) then     ! Ganser
     !
     do ic=1,nc
        psi(ic) = sphe(ic)
     end do
     !
  else if(modv.eq.3) then     ! Wilson
     !
     do ic = 1,nc
        call get_gama(gama,diam(ic),sphe(ic))  ! Get a/c
        if(gama.ge.1.0) then                ! oblate
           psi(ic) = 0.5*(1.0+1.0/gama)
        else                                ! prolate
           psi(ic) = gama
        end if
     end do
     !
  else if(modv.eq.4)  then    ! Dellino
     !
     do ic = 1,nc
        ! arnau OJO
        !               call get_gama(gama,diam(ic),sphe(ic))         ! Get a/c
        !               if(gama.ge.1.0) then                       ! oblate
        circula = 1.69
        psi(ic) = sphe(ic)/circula
        !	       else                                    ! prolate
        !                 circula = sqrt(0.5*(gama+1.0/gama))
        !                 psi(ic) = sphe(ic)/circula
        !               end if
     end do
  end if
  !
  return
end subroutine setpsi
!
!
!
subroutine get_gama(gama,diam,sphe)
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
  real(rp) :: gama,diam,sphe
  !
  integer(ip):: iiter,niter
  real(rp)   :: d,pi,gmin,gmax,Ao,toler,e
  real(rp)   :: Vp,Ap
  !
  !***   Initializations
  !
  d     = diam*1d4         ! see NOTE
  niter = 1000
  toler = 1d-6
  gmin  = 1d-3
  gmax  = 1.0
  !
  !***   Volume and area
  !
  pi = 4.0*atan(1.0)
  Vp = 4.0*pi*((0.5*d)**3)/3.0
  Ap = (pi**(1.0/3.0))*((6.0*Vp)**(2.0/3.0))/sphe
  !
  !***   Iterates
  !
  do iiter = 1,niter
     gama = 0.5*(gmin+gmax)
     e    = acos(gama)
     Ao   = 0.5*pi*d*d*(gama**(-4.0/3.0))*(gama*gama + (e/tan(e)))
     if(Ao.lt.Ap) then
        gmax = gama
     else
        gmin = gama
     end if
     if((iiter.gt.1).and.(abs(Ao-Ap).lt.toler)) goto 10
  end do
  call wriwar('Subroutine get_gama: convergence not achieved')
  !
  !***  convergence
  !
10 return
end subroutine get_gama
