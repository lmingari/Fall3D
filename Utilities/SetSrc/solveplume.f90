  subroutine solveplume(lexit)
  !************************************************************************
  !*
  !*    Solves the 1-D radially-averaged equations for a plume.
  !*
  !*      np        Number of plume points for output. Points are equally spaced
  !*                along the s direction (i.e. not in the vertical z due
  !*                to plume bent-over by wind). The points start at s=0
  !*                and finish at the neutral bouyancy level.
  !*      ns        Total number of points for output (Plume+Umbrella). Mass
  !*                in the umbrella region is distributed following a Gaussian.
  !*      npart      Number of particle classes.
  !*      nc         Number of classes (npart+ngas).
  !*                 Plume equations are solved for the particles only. Gas is
  !*                 distributed around the NBL.
  !*
  !*    OUTPUT
  !*
  !*      lexit          Returns .true. or .false. (collapsing column)
  !*      Xplum(ns)      X-coordinate (UTM)
  !*      Yplum(ns)      Y-coordinate (UTM)
  !*      Zplum(ns)      Z-coordinate (a.s.l.)
  !*      Hplum(ns)      tHeta angle
  !*      Uplum(ns)      Averaged velocity
  !*      Tplum(ns)      Temperature
  !*      Dplum(ns)      Bulk density
  !*      Rplum(ns)      Radius of the column
  !*      Mair (ns)      Mass of entrained air
  !*      Mplum(nc,ns)   Mass of class ic per unit of time in the plume
  !*                     OUTPUT: Mass that falls per unit time
  !*
  !*************************************************************************
  use KindType
  use InpOut
  use Plume
  use Master
  use Air
  use Coordinates
  implicit none
  !
  integer(ip), save :: ipass = 0
  !
  logical     :: lexit,go_on_s,go_on_MFR
  integer(ip) :: ic,ieq,is,istate,i,ismax
  real(rp)    :: s,so,ds,sb,rhorhoa,rhorhoamin,Mjet,n_MFR_min,n_MFR_max
  real(rp)    :: x,y,z,th,u,T,Ma,r,rho,rhoair
  real(rp)    :: Hb,Hc,dz,dMdz,Mumb
  !
  real(rp), save, allocatable :: fo(:),f(:),Mf0(:) ! work arrays for BC and solutions
  !
  !***  Allocates memory (depending on npart only)
  !
  if(ipass == 0) then
     ipass = 1
     allocate(fo (8+npart))
     allocate(f  (8+npart))
     allocate(Mf0(npart  ))
  end if
  !
  !***  Loop over MFR. It depends on the solving strategy (for MFR or height)
  !
  go_on_MFR = .true.
  MFR_iter  = 1
  n_MFR_min = n_MFR(1)   ! min MFR (e.g. 10**2 . From input file)
  n_MFR_max = n_MFR(2)   ! max MFR (e.g. 10**10. From input file)
  !
  do while(go_on_MFR)
     !
     select case(solve_plume_for)
     case('HEIGHT')
        M0 = M0                                 ! MFR given, a single iteration with MFR=M0 is sufficient
     case('MFR')
        M0 = 10.0_rp**(0.5_rp*(n_MFR_min+n_MFR_max))  ! M0 for this iteration
     end select
     !
     !***     Load boundary conditions for this MFR iteration
     !
     fo(1) = M0                                !  Eq(1). Q = pi*r*r*rho*u
     fo(2) = M0*u0                             !  Eq(2). P = pi*r*r*rho*u*u
     fo(3) = 0.5_rp*pi                         !  Eq(3). Theta
     fo(4) = M0*(Cbo*T0+g*(Z0+dZ0)+0.5*u0*u0)  !  Eq(4). J = pi*r*r*rho*u*E
     fo(5) = 0.0_rp                            !  Eq(5). Ma
     fo(6) = X0_UTM                            !  Eq(6). X in UTM, not in Lon
     fo(7) = Y0_UTM                            !  Eq(7). Y in UTM, not in Lat
     fo(8) = Z0 + dZ0                          !  Eq(8). Z  a.s.l.
     do ic = 1,npart
        fo(8+ic) = M0*(1.0_rp-w0)*fc(ic)       !  Eq(8+ic). Initial mass flux for ic
     end do
     !
     !***     Loop over s to find Sb, the Neutral Buoyancy Level point in s coordinate
     !
     go_on_s = .true.
     s       =  0.0_rp
     so      =  0.0_rp
     ds      = 10.0_rp
     do ieq  = 1,8+npart
        f(ieq) = fo(ieq)
     end do
     !
     rhorhoamin = 1d9
     rhomax = 0.0_rp
     do while(go_on_s)
        s  = s + ds
        call splum(so,f,s,8+npart,istate)
        !
        if(istate == 2) then
           call getplume(f,x,y,z,th,u,T,Ma,r,rho)
           rhoair = rhoa(z)
           rhorhoa = rho/rhoair

           if(rhorhoa <= rhorhoamin) rhorhoamin = rhorhoa

           sb = s
           if(rhorhoa <= (1.001_rp) .and. &
                rhorhoa >= (0.999_rp) .and. &
                rhorhoa > rhorhoamin) then  ! Skip the rho/rhoair=1 after jet
              go_on_s = .false.
              lexit = .true.
           end if
        else
           go_on_s = .false.
           lexit   = .false.   ! Wrong termination
        end if
     end do
     !
     select case(solve_plume_for)
     case('HEIGHT')
        if(lexit) then
           go_on_MFR = .false.   ! Done
        else
           return                ! Wrong termination
        end if
        !
     case('MFR')
        if(lexit) then
           !
           !***  Checks if the plume height is close to the required value.
           !***  Note that z is a.s.l.
           !
           if( abs((Z0+dZ0+(HPlume/c_umbrella)-8.0_rp*R0)-z) <= 10.0_rp) then
              go_on_MFR = .false.   ! Done
           else if( z > (Z0+dZ0+(HPlume/c_umbrella)-8.0_rp*R0) ) then
              n_MFR_max = 0.5_rp*(n_MFR_min+n_MFR_max)
           else
              n_MFR_min = 0.5_rp*(n_MFR_min+n_MFR_max)
           end if
           !
        else
           n_MFR_max = 0.5_rp*(n_MFR_min+n_MFR_max)
        end if
        !
        MFR_iter = MFR_iter + 1
        !
        if(MFR_iter.eq.25) call runend('MFR iterations = Itermax. Convergence not achieved')
        !
     end select

  end do              ! end do while(go_on_MFR)
  !
  !*** Integrates again up to the NBL and load the variables
  !
  s  = 0.0_rp
  so = 0.0_rp
  ds = sb/np            ! np points
  do ieq = 1,8+npart
     f(ieq) = fo(ieq)
  end do
  !
  rhomax = 0.0_rp
  do is=1,np
     s  = s + ds
     call splum(so,f,s,8+npart,istate)
     call getplume(f,x,y,z,th,u,T,Ma,r,rho)
     !
     Xplum(is) = x
     Yplum(is) = y
     Zplum(is) = z
     Lplum(is) = s
     Hplum(is) = th
     Uplum(is) = u
     Tplum(is) = T
     Mair (is) = Ma
     Rplum(is) = r
     Dplum(is) = rho
     do ic=1,npart
        Mplum(ic,is) = f(8+ic)
     end do

  end do
  !
  !*** If necessary, convert coordinates from UTM to LON-LAT
  !
  if(TRIM(coord_sys).eq.'LON-LAT') then
     do is=1,np
        x = Xplum(is)
        y = Yplum(is)
        call utm2ll(x,y,UTM_zone,LonPlume,LatPlume,WGS_84_DATUM,istate)
        if(istate.ne.0) then
           lexit = .false.
           return
        end if
        Xplum(is) = LonPlume
        Yplum(is) = LatPlume
     end do
  end if
  !
  !*** Distributes the remainder mass vertically along the umbrella
  !*** region using a linear distribution. The mass distribution is
  !*** stored from Mplum(ic,np) to Mplum(ic,ns), i.e. on ns-np+1 points.
  !*** The total column height Hc is assumed to be Hc=1.32*(Hb+8R0)=c_umbrella*(Hb+8R0)
  !
  Hb = Zplum(np)-Z0-dZ0         ! NBL   height (above vent)
  Hc = c_umbrella*(Hb+8.0_rp*R0)    ! Total height (above vent)
  dz = (Hc-Hb)/(ns-np)          ! ns-np+1 points
  do ic = 1,npart
     Mumb = Mplum(ic,np)        ! Remainder particle mass to distribute in umbrella
     dMdz = Mumb/(Hc-Hb)
     i = 0
     do is = np,ns
        z = i*dz
        i = i+1
        Mplum(ic,is) = Mumb - dMdz*z
     end do
  end do
  !
  do is = 1,ns-np
     Xplum(np+is) = Xplum(np)
     Yplum(np+is) = Yplum(np)
     Zplum(np+is) = Zplum(np)+(is)*dz
     Lplum(np+is) = Lplum(np+is-1) + sqrt( (Xplum(np+is)-Xplum(np+is-1))*(Xplum(np+is)-Xplum(np+is-1)) + &
                                           (Yplum(np+is)-Yplum(np+is-1))*(Yplum(np+is)-Yplum(np+is-1)) + &
                                           (Zplum(np+is)-Zplum(np+is-1))*(Zplum(np+is)-Zplum(np+is-1)) )
     Hplum(np+is) = 0.5_rp*pi
     Uplum(np+is) = Uplum(np)
     Tplum(np+is) = Tplum(np)
     Mair (np+is) = Mair (np)
     Rplum(np+is) = Rplum(np)
     Dplum(np+is) = Dplum(np)
  end do
  !
  !*    Modifies Mplum in order to store the mass that falls instead
  !*    of the mass in the plume
  !
  do ic = 1,npart
     Mf0(ic) = M0*(1.0_rp-w0)*fc(ic)
  end do
  !
  do ic=1,npart
     Mplum(ic,ns) =  Mplum(ic,ns-1)
  end do
  do is=ns-1,2,-1
     do ic=1,npart
        Mplum(ic,is)= Mplum(ic,is-1)-Mplum(ic,is)
     end do
  end do
  do ic=1,npart
     Mplum(ic,1)=Mf0(ic)-Mplum(ic,1)
  end do
  !
  !***  Put all the mass below zjet at z=zjet
  !
  ismax=1
  do ic=1,npart
     Mjet = 0.0_rp
     do is=1,ns
        if(ZPlum(is).lt.(Z0+dZ0+zjet)) then
           Mjet = Mjet + MPlum(ic,is)
           ismax = is
           MPlum(ic,is)= 0.0_rp
        end if
     end do
     MPlum(ic,ismax+1) = MPlum(ic,ismax+1) + Mjet
  end do
  !
  !***  If gas species exist in the granulometry file distribute (linearly) the
  !***  mass of each gas specie in the NBL
  !
  if(ngas > 0) then
     do ic = 1,ngas
        Mumb = M0*fc(npart+ic)          ! Total mass of gas specie
        i    = ns-np +1                 ! ns-np+1 points
        do is = np,ns                   ! umbrella region
           Mplum(npart+ic,is) = Mumb/i
        end do
     end do
  end if
  !
  !*** Writes the plume results
  !
  call wriplumeprop                        ! Writes plume properties for this interval
  call wriplumeheight                      ! Writes plume height for this interval
  call wriplumetem(time,min(time2,ieend1)) ! Writes plume temperature file
  !
  !*** If necessary, modifies the mass of each class to account for
  !*** aggregation accoding to Costa et al. (2010) & Folch et al. (2010)
  !
  if(type_aggr.eq.'COSTA') call costa
  call wriplumefmas                        ! Writes the mass distribution for this interval
  !
  return
end subroutine solveplume
!
!
!
subroutine getplume(f,x,y,z,th,u,T,Ma,r,rho)
  !********************************************************************
  !*
  !*    Extracts physical variables from values in the ODE's. at
  !*    a certain value of the arc parameter s.
  !*
  !*    OUTPUTS : x,y,z,th,u,T,Ma,r,rho
  !*
  !********************************************************************
  use KindType
  use Master
  use Plume
  use Air
  implicit none
  !
  real(rp) :: f(8+npart),x,y,z,th,u,T,Ma,r,rho
  !
  integer(ip) :: ic
  real(rp)    :: Mp,npp,rhopart
  !
  !***  Extracts values
  !
  x  = f(6)           ! Coordinates
  y  = f(7)
  z  = f(8)
  th = f(3)                                   ! Plume bent-over angle (in rad).
  u  = f(2)/f(1)                              ! Bulk velocity
  !        T  = ((f(4)/f(1))-g*z-0.5*u*u)/Cb        ! Bulk temperature
  Ma = f(5)                                   ! Mass flux of entrained air
  !
  Mp = 0.0_rp                 ! Total mass flux of particles
  do ic = 1,npart
     Mp = Mp + f(8+ic)
  end do
  !
  npp = Mp/(Mp+M0*w0+Ma)    ! Particle mass fraction (with air and volatiles). Q= Ma+Mp+M0*w0
  !
  rhopart= 0.0_rp             ! Averaged particle density
  do ic = 1,npart
     rhopart= rhopart + f(8+ic)/rhop(ic)
  end do
  rhopart = rhopart/Mp
  !
  Cb = (Ma*Ca+M0*w0*Cv+Mp*Cp)/(Mp+M0*w0+Ma)   ! Bulk heat capacity
  !
  T  = ((f(4)/f(1))-g*z-0.5_rp*u*u)/Cb        ! Bulk temperature
  rho = npp*rhopart + (1.0_rp-npp)*T/(rhoa(z)*Ta(z))  ! Averaged bulk density
  rho = 1.0_rp/rho
  rhomax=max(rho,rhomax)
  !
  r   = sqrt(f(1)/(pi*rho*u))   ! Plume radius
  !
  return
end subroutine getplume
!
!
!
subroutine splum(so,f,s,neq,istate)
  !************************************************************************
  !*
  !*    Returns the value of f(neq,s ) (at point s ) given the
  !*    initial conditions   f(neq,so) (at point so)
  !*
  !************************************************************************
  use KindType
  use odepack
  implicit none
  !
  integer  ::  neq,istate
  real(rp) ::  f(neq),so,s
  !
  integer  ::  nneq(1),itol,itask,iopt,lrw,liw,mf,iwork(20)
  real(rp) ::  rtol(1),atol(1),rwork(20+16*neq)
  external ::  dfds,jac
  !
  !***  Defines lsode parameters
  !
  nneq(1) = neq
  itol   = 1
  itask  = 1                 ! Normal computation (overshooting)
  istate = 1
  iopt   = 0
  lrw    = 20+16*neq
  liw    = 20
  mf     = 10                ! Jacobian is not supplied
  rtol(1) = 1e-6_rp          ! Relative error
  atol(1) = 1e-6_rp          ! Absolute error

  call lsode(dfds,nneq,f,so ,s,itol,rtol,atol,itask, &
       istate,iopt,rwork,lrw,iwork,liw,jac,mf)

  return
end subroutine splum

subroutine jac
  ! This is a dummy routine needed by lsode
end subroutine jac
!
!
subroutine dfds(neq,s,f,E)
  !*************************************************************************
  !*
  !*    Defines the system of equations for lsode loading E(i)
  !*
  !*    df
  !*    -- = E(s,f)    i = 1:neq
  !*    ds
  !*
  !*************************************************************************
  use KindType
  use Master
  use Plume
  use Air
  implicit none
  !
  integer     :: neq
  real(rp)    :: s,f(*),E(*)
  !
  logical     :: conve
  integer     :: ic
  real(rp)    :: x,y,z,th,u,T,Ma,r,rho                    ! Derivated variables
  real(rp)    :: Taair,rhoair,muair,Vaair,ue,vlimit,dvaadz
  real(rp)    :: zs,dlnAdz,A,Ri,a_shear,a_vortex,Ljet,Fb0,Fm0
  !
  !***  1. Current values of derived variables
  !
  call getplume(f,x,y,z,th,u,T,Ma,r,rho)
  !
  Taair  = Ta   (z)                   ! Air temperature
  rhoair = rhoa (z)                   ! Air density
  muair  = mua  (z)                   ! Air viscosity
  Vaair  = Va   (z)
  dVaadz = dVadz(z)
  !
  !***  Entrainment velocity. For strong plumes Carazzo et al. 2008
  !
  R0      = min(sqrt(M0/(pi*rhomax*U0)),350.0_rp)
  zs      = max(z/(2.0_rp*R0),0.0_rp)                 ! zstar = z/Dvent
  A       = 1.6_rp*(zs*zs+1657.0_rp)/(zs*zs+2411.0_rp)
  dlnAdz  = 1508.0_rp*zs/((zs*zs+2411.0_rp)*(zs*zs+1657.0_rp))
  Ri      = max(g*(rhoair-rho)*r/(rhoair*u*u+1e-6_rp),1e-6_rp)
  !
  Fb0 = pi*R0*R0*U0*g*(rhomax-1.0_rp)/1.0_rp           !rhoa(0)=1.0
  Fm0 = M0*U0
  Ljet= min(5.0_rp*Fm0**0.75_rp/sqrt(Fb0),1000.0_rp*R0)
  !
  a_shear = 0.0675_rp + (1.0_rp-(1.0_rp/A))*Ri + (0.5_rp*r*dlnAdz)/(2.0_rp*R0)
  a_vortex= 0.34_rp*( Vaair/(sqrt(abs(Ri))*U0) )**0.125_rp

  a_shear  = min(a_shear ,0.2_rp)
  a_vortex = min(a_vortex,0.8_rp)

  if(zs <= 5.0_rp) then
     a_shear = 0.008_rp  !0.0675
     a_vortex= 0.0_rp    ! 0.01
  elseif(zs > 5.0_rp .and. z <= Ljet) then
     a_shear = 0.0675_rp ! 0.06
     a_vortex= 0.0_rp
  endif
  !
  ue   = a_shear*abs(u-Vaair*cos(th)) + a_vortex*abs(Vaair*sin(th))
  !
  !***    Computes xi parameter (Bursik 2001 gives 0.23)
  !
  xi = 0.1_rp
  !
  !***  2. Defines the equations
  !
  do ic = 1,npart
     call vsettl(diam(ic),rhop(ic),rhoair,muair,vlimit,modv,psi(ic),conve)
     if(.not.conve) call runend('No convergence in vsettl routine')
     E(8+ic) = -xi*vlimit*f(8+ic)/(u*r)
  end do
  !
  E(1) =  2.0_rp*pi*r*rhoair*ue
  do ic = 1,npart
     E(1) = E(1) + E(8+ic)
  end do
  !
  E(2) = pi*r*r*(rhoair-rho)*g*sin(th) + Vaair*cos(th)*E(1)
  do ic = 1,npart
     E(2) = E(2) + u*E(8+ic)
  end do
  !
  E(3) = pi*r*r*(rhoair-rho)*g*cos(th) - Vaair*sin(th)*E(1)
  E(3) = E(3)/(pi*r*r*rho*u*u)

  E(4) = 2.0_rp*pi*r*rhoair*ue*(Ca*Taair + g*z + 0.5_rp*ue*ue)
  do ic = 1,npart
     !         E(4) = E(4) + (Cp*T+g*z+0.5*u*u)*E(8+ic)
     E(4) = E(4) + Cp*T*E(8+ic)
  end do
  !
  E(5) = 2.0_rp*pi*r*rhoair*ue
  !
  E(6) = cos(th)*cos(PSIa)
  !
  E(7) = cos(th)*sin(PSIa)
  !
  E(8) = sin(th)
  !
  return
end subroutine dfds
!
!
!
subroutine vsettl(diam,rhop,rhoa,visc,vset,model,psi,conve)
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
  !*    psi      - Model factor
  !*
  !*****************************************************************************
  use KindType
  implicit none
  !
  integer(ip) :: model, it
  logical :: conve
  !
  integer(ip) :: maxit
  real(rp) ::  diam,rhop,rhoa,visc,vset,psi
  real(rp) ::  gi,eps,rey,rk1,rk2,a,b,vold,cd100   ! work variables
  real(rp), save ::  cd = 1.0
  !
  conve = .true.
  gi=9.81_rp                ! Gravity acceleration
  eps=1e-3_rp               ! Tolerance
  maxit=100                 ! Maximum number of iterations
  vold =1e5_rp
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
        if(rey <= 988.947_rp) then  ! This is the actual transition point
           cd=24.0_rp/rey*(1.0_rp+0.15_rp*rey**0.687_rp)
        else
           cd=0.44_rp
        endif
        if(it.gt.1.and.abs(vset-vold).le.eps) goto 10
        vold=vset
     enddo
  endif
  !
  !***  2. Model of Ganser 1993
  !
  if(model == 2)  then
     do it=1,maxit
        vset=sqrt(4.0*gi*diam*rhop/(3.0*cd*rhoa))
        rey=rhoa*vset*diam/visc
        rk1=3.0_rp/(1.0_rp+2.0_rp/sqrt(psi))
        rk2=10.0**(1.8148*(-log10(psi))**0.5743)
        cd=24.0/(rey*rk1)*(1.0+0.1118*     &
             (rey*rk1*rk2)**0.6567)+0.4305*rk2/   &
             (1.0+3305.0/(rey*rk1*rk2))
        if(it.gt.1.and.abs(vset-vold).le.eps) goto 10
        vold=vset
     enddo
  endif
  !
  !**  3. Model of Wilson & Huang 1979
  !
  if(model == 3) then
     do it=1,maxit
        vset=sqrt(4.0*gi*diam*rhop/(3.0*cd*rhoa))
        rey=rhoa*vset*diam/visc
        if(rey.le.100.0) then
           cd=24.0/rey*psi**(-0.828)+2.0*sqrt(1.07-psi)
        elseif(rey.gt.100.0.and.rey.lt.1000.0) then
           cd100=0.24*psi**(-0.828)+2.0*sqrt(1.07-psi)
           a=(1.0-cd100)/900.0
           b=1.0-1000.0*a
           cd=a*rey+b
        else
           cd=1.0
        endif
        !
        if(it.gt.1.and.abs(vset-vold).le.eps) goto 10
        vold=vset
     enddo
  endif
  !
  !***  4. Model of Dellino et al.
  !
  if(model == 4) then
     vset=((diam*diam*diam*gi*(rhop-rhoa)*rhoa*(psi**1.6))/(visc*visc))**0.5206
     vset=(1.2065*visc*vset)/(diam*rhoa)
     return
  end if
  !
  conve = .false.
  return
  !
10 continue
  vset=sqrt(4.0*gi*diam*rhop/(3.0*cd*rhoa))
  !
  return
end subroutine vsettl
