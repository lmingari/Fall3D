subroutine setvset
  !***********************************************************************
  !*
  !*    Set terminal fall velocity in the computational domain
  !*
  !*    dz       - Array of sizes of the cells in Z direction (Input)
  !*    rho      - Atmospheric density (Input Matrix)
  !*    vset     - Particle settling velocity (Output Matrix)
  !*
  !***********************************************************************
  use KindType
  use Numeric
  use Master
  implicit none
  !
  logical     :: conve
  integer(ip) :: ic,i,j,k
  real   (rp) :: z,visc,visc0,Tair0,T
  !
  !***  Standard air properies at sea level
  !
  visc0 = 1.827d-5  ! reference viscosity
  Tair0 = 291.15    ! reference temperature
  !
  do ic=1,np        ! only particles, gases have zero settling velocity
     do i=1,nx
        do j=1,ny
           !
           z = topg(i,j)-dz(0)    ! Initialize
           !
           !***           Begin loop on the layers
           !
           do k=1,nz
              z = z+dz(k-1)
              T = tempe(i,j,k)
              ! Sutherland's law
              visc = visc0*((Tair0+120.0_rp)/(T+120.0_rp))*((T/Tair0)**1.5_rp)
              call vsettl(diam(ic),rhop(ic),rho(i,j,k),visc,vset(i,j,k,ic), &
                   modv,psi(ic),sphe(ic),conve)
              if(.not.conve) call runend('No convergence in vsettl routine')
           enddo
        enddo
     enddo
  enddo
  !
  !***  Correct settling velocities of aggregates
  !
  if(classname(np).eq.'aggregate') then
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx
              vset(i,j,k,ic) = vset_factor*vset(i,j,k,ic) 
           end do
        end do
     end do
  end if
  !
  !***  Sign criteria
  !
  vset = -vset
  !
  return
end subroutine setvset
!
subroutine vsettl(diam,rhop,rhoa,visc,vset,model,psi,sphe,conve)
  !*****************************************************************************
  !*
  !*    Set fall velocity, as a function of particles diameter and
  !*    air density and viscosity.
  !*
  !*    diam     - Particle diameter in meters (Input)
  !*    rhop     - Particle density (Input)
  !*    rhoa     - Air density (Input vector)
  !*    visc     - Air viscosity (Input)
  !*    vset     - Particle settling velocity (Output)
  !*    model    - Settling velocity model (Input)
  !*
  !*****************************************************************************
  use KindType
  implicit none
  !
  integer(ip) :: model
  integer(ip) :: it
  logical :: conve
  !
  integer(ip) :: maxit
  real(rp)    :: diam,rhop,rhoa,visc,vset,psi,sphe
  real(rp)    :: gi,eps,rey,rk1,rk2,a,b,vold,cd100,dnd,gama   ! work variables
  real(rp), save :: cd = 1.0_rp
  !
  conve = .true.
  gi=9.81_rp              ! Gravity acceleration
  eps=1d-3                ! Tolerance
  maxit=100               ! Maximum number of iterations
  vold =1d5
  !
  !***  Begin iteration to compute vset
  !
  !
  !***  1. Model of Arastoopour et al. 1982
  !
  if(model == 1) then
     do it=1,maxit
        vset=sqrt(4.0_rp*gi*diam*rhop/(3.0_rp*cd*rhoa))
        rey=rhoa*vset*diam/visc
        if(rey <= 988.947_rp) then
           cd=24.0_rp/rey*(1.0_rp+0.15_rp*rey**0.687_rp)
        else
           cd=0.44_rp
        endif
        if(it > 1.and.abs(vset-vold) <= eps) goto 10
        vold=vset
     enddo
  endif
  !
  !***  2. Model of Ganser 1993
  !
  if(model == 2)  then
     call get_gama(diam,sphe,gama)  ! Get gama=a/c
     dnd=0.5_rp*(1.0_rp+gama)*gama**(-2.0/3.0)
     !
     do it=1,maxit
        vset=sqrt(4.0_rp*gi*diam*rhop/(3.0_rp*cd*rhoa))
        rey=rhoa*vset*diam/visc
        !		   rk1=1.0_rp/(1.0_rp/3.0_rp+2.0_rp/(3.0_rp*sqrt(psi)))

        rk1=1.0_rp/(dnd/3.0_rp+2.0_rp/(3.0_rp*sqrt(psi)))
        rk2=10.0_rp**(1.8148_rp*(-log10(psi))**0.5743_rp)
        cd =24.0_rp/(rey*rk1)*(1.0_rp+0.1118_rp*            &
             (rey*rk1*rk2)**0.6567_rp)+0.4305_rp*rk2/      &
             (1.0_rp+3305.0_rp/(rey*rk1*rk2))

        if(it > 1.and.abs(vset-vold) <= eps) goto 10
        vold=vset
     enddo
  endif
  !
  !***  3. Model of Wilson & Huang 1979
  !
  if(model == 3) then
     do it=1,maxit
        vset=sqrt(4.0_rp*gi*diam*rhop/(3.0_rp*cd*rhoa))
        rey=rhoa*vset*diam/visc

        if(rey <= 100.0_rp) then
           cd=24.0_rp/rey*psi**(-0.828_rp)+2.0_rp*sqrt(1.0_rp-psi)
        elseif(rey > 100.0_rp.and.rey < 1000.0_rp) then
           cd100=0.24_rp*psi**(-0.828_rp)+2.0_rp*sqrt(1.0_rp-psi)
           a=(1.0_rp-cd100)/900.0_rp
           b=1.0_rp-1000.0_rp*a
           cd=a*rey+b
        else
           cd=1.0_rp
        endif

        if(it > 1 .and. abs(vset-vold) <= eps) goto 10
        vold=vset
     enddo
  endif
  !
  !***  4. Model of dellino et al. 2005
  !
  if(model == 4) then
     vset = ((diam*diam*diam*gi*(rhop-rhoa)*rhoa*(psi**1.6_rp))/ &
          (visc*visc))**0.5206_rp
     vset = (1.2065_rp*visc*vset)/(diam*rhoa)
     return
  end if
  !
  !*** No convergence
  !
  conve = .false.
  return
  !
  !*** Convergence
  !
10 continue
  vset=sqrt(4.0_rp*gi*diam*rhop/(3.0_rp*cd*rhoa))
  !
  return
end subroutine vsettl
