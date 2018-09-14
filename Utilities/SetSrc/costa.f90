subroutine costa
  !********************************************************
  !*
  !*    Computes aggregation according to Costa et al. (2010)
  !*
  !********************************************************
  use KindType
  use InpOut
  use Master
  use Air
  implicit none
  !
  logical           :: conve
  integer(ip), save :: ipass = 0
  integer(ip)       :: is,ic,jc,iaggr
  real   (rp)       :: alfa_mean,sumfc,sumff,w,rhopmean,Stij,viscl,alfa,Stcr,qq
  real   (rp)       :: Ab,As,Ad,psi3,psi4,gamas,duds,ntot,dntot,vij
  real   (rp)       :: rhoair,muair,Taair,fi,expb,exps,expd,sumNi,ka
  character(len=6)  :: wphase
  !
  real   (rp)       :: time_plum,H_min,H_max,U_min,U_max
  !
  real   (rp), save, allocatable :: Ni  (:)
  real   (rp), save, allocatable :: conc(:)
  real   (rp), save, allocatable :: vset(:)
  real   (rp), save, allocatable :: M_in_plum(:,:)
  !
  !*** Initializations
  !
  H_min = 5.0d4
  H_max = 5.0d0
  U_min = 5.0d2
  U_max = 5.0d-2
  !
  !***  Allocates memory (first time only)
  !
  if(ipass == 0) then
     ipass = 1
     allocate(Ni  (npart))
     allocate(conc(npart))
     allocate(vset(npart))
     allocate(M_in_plum(npart,ns))
  end if
  !
  !*** Index of the first aggregating class
  !
  iaggr = 0
  do ic = npart-1,1,-1
     if(diam_aggr.gt.diam(ic)) then
        iaggr = ic
     end if
  end do
  if(iaggr.eq.0) iaggr = npart-1  ! ensure at least 1 aggregating class
  !
  !*** Previous computations
  !
  sumfc = 0.0_rp
  do ic = iaggr,npart-1
     sumfc = sumfc + fc(ic)
  end do
  !
  sumff = 0.0_rp
  do ic = iaggr,npart-1
     do jc = iaggr,npart-1
        sumff = sumff + fc(ic)*fc(jc)
     end do
  end do
  !
  rhopmean = 0.0_rp                ! mean density of primary particles
  do ic = iaggr,npart-1
     rhopmean = rhopmean + fc(ic)*rhop(ic)
  end do
  rhopmean = rhopmean/sumfc
  !
  !*** First, recover the mass in the plume (only the mass that falls
  !*** is stored in Mplum). This is needed to estimate particle mass
  !*** concentration.
  !
  M_in_plum = 0.0_rp
  !
  do ic = 1,npart
     M_in_plum(ic,ns) = Mplum(ic,ns)  ! last layer
  end do
  do is = ns-1,1,-1
     do ic = 1,npart
        M_in_plum(ic,is) = M_in_plum(ic,is+1) + Mplum(ic,is)
     end do
  end do
  !
  !*** Results file header
  !
  write(lures,1) (TRIM(classname(ic)),ic=1,npart-1)
1 format(/,  &
       'AGGREGATION SUMMARY (Costa model): fraction of aggregated mass (in %)',/,&
       '   z (km a.s.l.)   w_phase ',20(2x,a))
  !
  !*** Loop over all the plume slabs (ns-1) above the zjet...
  !
  do is = 2,ns
     !
     if(MPlum(1,is).gt.0.0_rp) then
        !
        !*** Mass concentration at each plume point. Top-hat profile and
        !*** cilindric slabs are assumed
        !
        do ic = 1,npart
           conc(ic) = M_in_plum(ic,is)/(pi*Rplum(is)*Rplum(is)*(Lplum(is)-Lplum(is-1)))
        end do
        !
        !*** Particle settling velocity at each plume point
        !
        rhoair = rhoa (Zplum(is))           ! Air density
        muair  = mua  (Zplum(is))           ! Air viscosity
        Taair  = Ta   (Zplum(is))           ! Air temperature
        do ic = 1,npart
           call vsettl(diam(ic),rhop(ic),rhoair,muair,vset(ic),modv,psi(ic),conve)
           if(.not.conve) vset(ic) = 0.0_rp
        end do
        !
        !*** Determine the water phase. In the Plume it is assumed that:
        !***   T < 255K       ice
        !***   373 > T > 255  water
        !***   T > 373        vapor
        !
        !*** Also computes time. Note that time does not appear in the steady plume model. However
        !*** time is needed by the aggregation model to compute particle decay. As a time scale
        !*** we compute the time that a particle spends in the ice-water window of the plume,
        !*** where aggregation can occurr.
        !
        H_min = min(Zplum(is-1),H_min)
        H_max = max(Zplum(is  ),H_max)
        U_max = max(Uplum(is-1),U_max)
        U_min = min(Uplum(is  ),U_min)
        !
        if(Tplum(is).lt.255.) then
           wphase = 'ice'
           !
        else if (Tplum(is).lt.373.) then
           wphase = 'water'
           !
        else
           wphase = 'vapour'
           !
        end if
        !
        !*** Collision frequency Brownian motion at each point
        !
        Ab = -4.0_rp*1.38d-23*Tplum(is)/(3.0_rp*muair)
        !
        !*** Collision frequency for laminar and turbulent fluid sheat at each point
        !*** Given the 1D radially-averaged approach, shear is estimated only from duds
        !*** and the dissipation rate from (dud2)**2
        !
        psi3  = 6.0_rp/pi
        duds  = (Uplum(is)-Uplum(is-1))/(Lplum(is)-Lplum(is-1))
        gamas = sqrt(1.3*2.0*duds*duds)  ! turbulent flow
        As    = -2.0_rp*gamas*psi3/3.0_rp
        !
        !*** Collision frequency for differential settling velocity at each point
        !
        psi4 =  psi3*((6.0_rp/pi)**(1.0_rp/3.0_rp))
        Ad   = -pi*(rhopmean-rhop(npart))*g*psi4/(48.0_rp*muair)
        !
        !*** Nummer of particles per unit volume that can potentially participate in
        !*** the aggreagtion process and solid volume fraction
        !
        !
        !*** Computes the mean sticking efficiency (alfa mean) at each point
        !
        SELECT CASE(wphase)
        case('ice')
           !
           alfa_mean = 0.09_rp
           time_plum = 2.0_rp*(H_max-H_min)/(U_max+U_min)
           !
        case('water')
           !
           time_plum = 2.0_rp*(H_max-H_min)/(U_max+U_min)
           !
           viscl = 2.414d-5*(10.0_rp**(247.7_rp/(Taair-140.0_rp)))  ! water viscosity
           Stcr  = 1.3_rp                                           ! Critical Stokes number
           qq    = 0.8_rp
           !
           alfa_mean = 0.0_rp
           do ic = iaggr,npart-1
              do jc = iaggr,npart-1
                 vij=abs(vset(ic)-vset(jc))-Ab/(2.0_rp*pi*diam(ic)*diam(jc))-As*diam(ic)*diam(jc)/(4.0_rp*pi)
                 w    = fc(ic)*fc(jc)/sumff
                 Stij = 8.0_rp*rhopmean*vij*diam(ic)*diam(jc)/ &
                      (9.0_rp*viscl*(diam(ic)+diam(jc)))
                 alfa = 1.0_rp+((Stij/Stcr)**qq)
                 alfa = 1.0_rp/alfa
                 alfa_mean = alfa_mean + w*alfa
              end do
           end do
           !
        case('vapour')
           !
           alfa_mean = 0.0_rp      ! No wet aggregation
           time_plum = 0.0_rp
           !
        END SELECT
        !
        ntot = 0.0_rp
        fi   = 0.0_rp
        do ic = iaggr,npart-1
           fi   = fi   + conc(ic)/rhop(ic)
           ntot = ntot + 6.0_rp*conc(ic)/(pi*rhop(ic)*diam(ic)*diam(ic)*diam(ic))
        end do
        !
        !*** Total particle decay per unit volume during this time interval
        !
        expb = 2.0_rp                ! Brownian     ntot exponent
        exps = 2.0_rp-(3.0_rp/Df)    ! Shear        ntot exponent
        expd = 2.0_rp-(4.0_rp/Df)    ! Differential ntot exponent
        !
        dntot = Ab*(ntot**expb) + As*(ntot**exps)*(fi**(3.0_rp/Df)) + Ad*(ntot**expd)*(fi**(4.0_rp/Df))
        dntot = alfa_mean*dntot*time_plum
        !
        sumNi = 0.0_rp
        Ni    = 0.0_rp
        ka    = 1.0_rp                ! Fractal prefactor
        do ic = iaggr,npart-1
           Ni(ic) = ka*( (diam(npart)/diam(ic))**Df )
           sumNi  = sumNi + Ni(ic)
        end do
        !
        !*** Number of particles of class j per unit volume that aggregate
        !*** (stored in Ni)
        !
        do ic = iaggr,npart-1
           Ni(ic) = dntot*Ni(ic)/sumNi
        end do
        !
        !*** Mass of aggregates (stored in Ni)
        !
        do ic = iaggr,npart-1
           Ni(ic) = pi*rhop(ic)*diam(ic)*diam(ic)*diam(ic)*Ni(ic)/6.0_rp   ! per unit volume and time
           Ni(ic) = Ni(ic)*pi*Rplum(is)*Rplum(is)*(Lplum(is)-Lplum(is-1))  ! per unit time
           Ni(ic) = -Ni(ic)
           Ni(ic) = min(Ni(ic),Mplum(ic,is)) ! limit the total amount of mass
        end do
        !
        !*** Writes aggregation summary to the results file
        !
        write(lures,10) zplum(is),TRIM(wphase),(100*(1.0-(Mplum(ic,is)-Ni(ic))/Mplum(ic,is)),ic=1,npart-1)
10      format(2x,f8.1,8x,a6,2x,100(3x,f7.2))
        !
        !*** Finally, modifies the distribution of mass that falls from the plume
        !
        do ic = iaggr,npart-1
           Mplum(ic   ,is) = Mplum(ic   ,is) - Ni(ic)
           Mplum(npart,is) = Mplum(npart,is) + Ni(ic)
           !
        end do
        !
     end if  ! z=zjet
     !
  end do
  !
  return
end subroutine costa
