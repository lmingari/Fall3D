subroutine agrsrc
  !*****************************************************************
  !*
  !*    This routine computes tarat, the aggregation source
  !*    term (only for particles).
  !*
  !*****************************************************************
  use KindType
  use Numeric
  use Master
  use Domain, only: collect_c,ic_loc, nc_loc, iz_loc, nz_loc
  implicit none
  !
  integer(ip) :: ix,iy,iz,ic
  real   (rp) :: fi,ntot,dntot,dnj
  real   (rp) :: expb,exps,expd,ka,sumNi
  real   (rp), allocatable :: Ni(:),d3(:)
  !
  !***  Initializations
  !
  ka   = 1.0_rp                ! Fractal prefactor
  expb = 2.0_rp                ! Brownian     ntot exponent
  exps = 2.0_rp-(3.0_rp/Df)    ! Shear        ntot exponent
  expd = 2.0_rp-(4.0_rp/Df)    ! Differential ntot exponent
  !
  allocate(Ni(nc))
  allocate(d3(nc))
  sumNi = 0.0_rp
  do ic = iaggr,np-1
     d3(ic) = pi*rhop(ic)*diam(ic)*diam(ic)*diam(ic)
     Ni(ic) = ka*( (diam(np)/diam(ic))**Df )
     sumNi  = sumNi + Ni(ic)
  end do
  !
  !***  Collect for classes (not for z-layer)
  !
  call collect_c( c(0:,0:,0:,:), nc )
  !
  !***  Loop over domain
  !
  do iz = iz_loc,iz_loc+nz_loc-1
     do iy = 1,ny
        do ix = 1,nx
           !
           ntot = 0.0_rp
           fi   = 0.0_rp
           do ic = iaggr,np-1
              fi   = fi   + c(ix,iy,iz,ic)/rhop(ic)
              ntot = ntot + 6.0_rp*c(ix,iy,iz,ic)/d3(ic)
           end do
           ntot  = ntot/Hm1(ix,iy)        ! concentration is scaled
           fi    = fi/Hm1(ix,iy)          ! concentration is scaled
           !
           dntot = Abagr(ix,iy,iz)*(ntot**expb)                   + &
                Asagr(ix,iy,iz)*(ntot**exps)*(fi**(3.0_rp/Df)) + &
                Adagr(ix,iy,iz)*(ntot**expd)*(fi**(4.0_rp/Df))
           dntot = Alagr(ix,iy,iz)*dntot
           !
           !***     Calculates source (particles that aggregate only)
           !
           tarat(ix,iy,iz,np) = 0.0_rp
           do ic = iaggr,np-1
              !
              !***  Scale source term (in kg/(m3*s)) and sign criteria
              !
              !
              !*** Number of particles that aggragate per unit volume and unit time
              !
              dnj = dntot*Ni(ic)/sumNi
              ! Negative sign (sink)
              ! tarat(ix,iy,iz,ic) = dnj*d3(ic)*Hm1(ix,iy)/6.0_rp

              ! arnau this has to be check. This is mass per unit volume and time (kg/m3*s)
              tarat(ix,iy,iz,ic) = dnj*d3(ic)/6.0_rp ! Negative sign (sink)
              !
              ! Class mass conservation
              if( (c(ix,iy,iz,ic)+dt*tarat(ix,iy,iz,ic)) < 0.0_rp) &
                   tarat(ix,iy,iz,ic) = -c(ix,iy,iz,ic)/dt
              ! Positive sign (source)
              tarat(ix,iy,iz,np) = tarat(ix,iy,iz,np) - tarat(ix,iy,iz,ic)
           end do
           !
        end do
     end do
  end do
  !
  deallocate(Ni)
  deallocate(d3)
  !
  return
end subroutine agrsrc
!
!
!
subroutine agrpar
  !*****************************************************************
  !*
  !*    This routine computes aggregation arrays related to
  !*    meteorological data and source term
  !*
  !*****************************************************************
  use KindType
  use Numeric
  use Master
  use InpOut
  use Parallel
  use Domain , only: collect_c,collect_cload,collect_cz, ic_loc, nc_loc      
  implicit none
  !
  integer(ip) :: ix,iy,iz,ic,jc    
  real   (rp) :: T,visc,visc0,Tair0,rhopmean,rhoa,wcloud,viscl
  real   (rp) :: w,sumfc,sumff,alfa,Stcr,Stij,qq,bb,B,kb,vup,vdw,gamas,psi3,psi4,z, & 
       g,dmax,dmin,dzeta,vol0,vij
  character(len=5) :: wphase 
  !
  !***  Collect for classes (not for z-layer)
  !          
  do ic = ic_loc, ic_loc + nc_loc - 1
     call collect_cz( c( :, :, :, ic ), nz )  !  return to a serial behaviour
  end do
  call collect_c( c(0:,0:,0:,:), nc )
  !
  !***  Initializations
  !
  g     = 9.81_rp
  kb    = 1.38d-23                 ! Boltzman constant
  visc0 = 1.827d-5                 ! reference viscosity
  Tair0 = 291.15                   ! reference temperature
  bb    = 0.44_rp
  B     = (72.0_rp*4.8d-4)/pi
  Stcr  = 1.3_rp                   ! Critical Stokes number
  qq    = 0.8_rp
  psi3  = 6.0_rp/pi
  psi4  = psi3*((6.0_rp/pi)**(1.0_rp/3.0_rp))
  !
  sumfc    = 0.0_rp
  do ic = iaggr,np-1           ! Only particles      
     sumfc = sumfc + fc(ic)
  end do
  !
  sumff    = 0.0_rp
  do ic = iaggr,np-1
     do jc = iaggr,np-1      
        sumff = sumff + fc(ic)*fc(jc)
     end do
  end do
  !
  rhopmean = 0.0_rp                ! mean density of primary particles
  do ic = iaggr,np-1      
     rhopmean = rhopmean + fc(ic)*rhop(ic)
  end do
  rhopmean = rhopmean/sumfc
  !
  !***  Initializations
  !
  Alagrmax_ice = 0.0_rp
  Alagrmax_wat = 0.0_rp
  Alagrmin_ice = 1d20
  Alagrmin_wat = 1d20
  !
  !***  Loop over domain
  !
  do iz = 1,nz
     do iy = 1,ny
        do ix = 1,nx
           !
           !***     Temperature (air or plume) and density
           !
           if(tplume(ix,iy,iz).gt.0.0_rp) then
              T = tplume(ix,iy,iz)
           else
              T = tempe(ix,iy,iz) 
           end if
           rhoa = rho  (ix,iy,iz)
           !
           !***     Air viscosity
           !       
           visc = visc0*((Tair0+120.0_rp)/(T+120.0_rp))*((T/Tair0)**1.5)   ! Sutherland's law
           !
           !***     Water viscosity
           !
           viscl = 2.414d-5*(10.0_rp**(247.7_rp/(T-140.0_rp)))
           !
           !***     Chek wether we are in the plume or in the cloud and determine
           !***     the water phase
           !
           if(tplume(ix,iy,iz).gt.0.0_rp) then
              !
              !           In the Plume it is assumed that 
              !           T < 255K       ice
              !           373 > T > 255  water
              !           T > 373        vapor
              !
              if(T.lt.255) then
                 wphase = 'ice'
              else if (T.lt.373) then
                 wphase = 'water'
              else
                 wphase = 'vapour'
              end if
              !
           else
              !
              !          In the cloud              
              !          If there is some gas component it is assumed that
              !          the first component np+1 is (magmatic) water. 
              !
              if(ng.gt.0) then 
                 wcloud = c(ix,iy,iz,np+1)/Hm1(ix,iy)    
              else
                 wcloud = 0.
              end if
              wcloud = wcloud/rhoa   ! mixing ratio in kg/kg
              !         
              call water_phase(qv(ix,iy,iz),P(ix,iy,iz),T,wcloud,wphase)
              !
           end if
           !
           !***     Computes A Brownian function at each point
           !
           Abagr(ix,iy,iz) = -4.0_rp*kb*T/(3.0_rp*visc)
           !
           !***     Computes A Shear at each point
           !        
           call get_gamas(ix,iy,iz,gamas) 
           Asagr(ix,iy,iz) = -2.0_rp*gamas*psi3/3.0_rp
           !
           !***     Computes A differential settling velocity at each point
           !
           Adagr(ix,iy,iz) = -pi*(rhopmean-rhoa)*g*psi4/(48.0_rp*visc)
           !
           !***     Computes the mean sticking efficiency (alfa mean) at each point
           !
           Alagr(ix,iy,iz) = 0.0_rp 
           !
           SELECT CASE (wphase)
           CASE('ice')
              !
              Alagr(ix,iy,iz) = 0.09_rp
              if(Alagr(ix,iy,iz).gt.Alagrmax_ice) Alagrmax_ice = Alagr(ix,iy,iz)
              if(Alagr(ix,iy,iz).lt.Alagrmin_ice) Alagrmin_ice = Alagr(ix,iy,iz)          
              !
           case('water')
              do ic = iaggr,np-1
                 do jc = iaggr,np-1
                    vij   =  abs(vset(ix,iy,iz,ic)-vset(ix,iy,iz,jc))- &
                         Abagr(ix,iy,iz)/(2.0_rp*pi*diam(ic)*diam(jc))- &
                         Asagr(ix,iy,iz)*diam(ic)*diam(jc)/(4.0_rp*pi)
                    w    =  fc(ic)*fc(jc)/sumff
                    Stij =  8.0_rp*rhopmean*vij*diam(ic)*diam(jc)/ &
                         (9.0_rp*viscl*(diam(ic)+diam(jc)))
                    alfa = 1.0_rp+((Stij/Stcr)**qq)
                    alfa = 1.0_rp/alfa                       
                    Alagr(ix,iy,iz) = Alagr(ix,iy,iz) + w*alfa
                 end do
              end do
              !
              if(Alagr(ix,iy,iz).gt.Alagrmax_wat) Alagrmax_wat = Alagr(ix,iy,iz)
              if(Alagr(ix,iy,iz).lt.Alagrmin_wat) Alagrmin_wat = Alagr(ix,iy,iz) 
              !
           case('vapor')
              continue      ! No wet aggregation
           END SELECT

           !                                                 
        end do
     end do
  end do
  !
  !***  Computes the maximum/minimum values of the aggregation parameters
  !
  Alagrmax = maxval(Alagr)
  Alagrmin = minval(Alagr)
  Alagrave = sum(Alagr)/(nx*ny*nz)
  !      
  Abagrmax = maxval(Abagr)
  Abagrmin = minval(Abagr)
  Abagrave = sum(Abagr)/(nx*ny*nz)
  !      
  Asagrmax = maxval(Asagr)
  Asagrmin = minval(Asagr)
  Asagrave = sum(Asagr)/(nx*ny*nz)
  !      
  Adagrmax = maxval(Adagr)
  Adagrmin = minval(Adagr)
  Adagrave = sum(Adagr)/(nx*ny*nz)
  !
  !***  (kernel order of magnitue only)
  !
  dmin = diam(np  )
  dmax = diam(iaggr)
  !
  Kbagrmax = Alagrmax*maxval(-Abagr)*(4.0_rp)
  kbagrmin = Alagrmin*minval(-Abagr)*(4.0_rp)
  !
  ksagrmax = Alagrmax*maxval(-Asagr)*((2.0_rp*dmax)**3.0)
  ksagrmin = Alagrmin*minval(-Asagr)*((2.0_rp*dmin)**3.0)
  !
  kdagrmax = Alagrmax*(maxval(-Adagr))*(3.0*dmax**4.0)
  kdagrmin = Alagrmin*(minval(-Adagr))*(3.0*dmin**4.0)
  !
  !***  Release memory
  !
  !      deallocate(mi)
  !
  !***  Corrects the setling velocity of the aggregate
  !
  do iz = 1,nz
     do iy = 1,ny
        do ix = 1,nx
           vset(ix,iy,iz,iaggr) = vset_factor*vset(ix,iy,iz,iaggr) 
        end do
     end do
  end do
  !
  !***  Writes to log file
  !
  if(mpime == root) then
     write(nlst,30) Alagrmax    ,Alagrmin    ,Alagrave, &
          Alagrmax_ice,Alagrmin_ice,          &
          Alagrmax_wat,Alagrmin_wat,          &
          Abagrmax    ,Abagrmin    ,Abagrave, &
          Asagrmax    ,Asagrmin    ,Asagrave, &
          Adagrmax    ,Adagrmin    ,Adagrave, &
          Kbagrmax    ,Kbagrmin, &
          Ksagrmax    ,Ksagrmin, &
          Kdagrmax    ,Kdagrmin                         
30   format(/, &
          'AGGREGATION PARAMETERS',/,  &    
          '  Alfa mean                   ',/, &
          '    Maximum                 : ',e12.4,/, &
          '    Minimum                 : ',e12.4,/, & 
          '    Average                 : ',e12.4,/, &
          '    Ice   maximum           : ',e12.4,/, &
          '    Ice   minimum           : ',e12.4,/, &                 
          '    Water maximum           : ',e12.4,/, &
          '    Water minimum           : ',e12.4,/, &
          '  A Brownian coefficient      ',/,       &
          '    Maximum                 : ',e12.4,/, &
          '    Minimum                 : ',e12.4,/, & 
          '    Average                 : ',e12.4,/, &
          '  A Shear coefficient         ',/,       &
          '    Maximum                 : ',e12.4,/, &
          '    Minimum                 : ',e12.4,/, & 
          '    Average                 : ',e12.4,/, &
          '  A Diff. settling coeff.     ',/,       &
          '    Maximum                 : ',e12.4,/, &
          '    Minimum                 : ',e12.4,/, & 
          '    Average                 : ',e12.4,/, &
          '  K Brownian kernel           ',/,       &
          '    Maximum                 : ',e12.4,/, &
          '    Minimum                 : ',e12.4,/, &              
          '  K shear kernel              ',/,       &
          '    Maximum                 : ',e12.4,/, &
          '    Minimum                 : ',e12.4,/, &  
          '  K differential kernel       ',/,       &
          '    Maximum                 : ',e12.4,/, &
          '    Minimum                 : ',e12.4)  
  end if
  !
  !
  return
end subroutine agrpar
!
!
!
subroutine water_phase(qv,P,T,wcloud,wphase)
  !*******************************************************************************
  !*
  !*    Computes the water phase
  !*
  !*******************************************************************************
  use KindType
  implicit none
  !
  real(rp) :: qv,P,T,wcloud
  real(rp) :: w,e,ew,ei
  character(len=5) :: wphase 
  !
  !***  Mixing ratio from specific humidity
  !
  w = qv/(1.0_rp-qv)
  !
  !***  Add the water cloud contribution
  !
  w = w + wcloud
  !
  !***  vapor pressure
  !
  e = P*w/(w+0.622_rp)
  !
  !***  saturation vapour pressure with respect to water
  !
  !      ew =  10.79574*(1-273.16/T)  & 
  !            -5.02800*Log10(T/273.16)  &
  !            +1.50475d-4*(1 - 10**(-8.2969*(T/273.16-1))) &
  !            +0.42873d-3*(10**(4.76955*(1-273.16/T)) - 1)  + 0.78614 
  !      ew = 100*(10.**ew)   ! hPa --> Pa
  !
  ew = 6.112_rp*exp(17.67_rp*(T-273.16_rp)/(T-273.16_rp+243.5_rp)) 
  ew = 100.0_rp*ew     ! hPa --> Pa
  !                   
  !***   saturation vapour pressure with respect to ice
  ! 
  ei = -9.09718_rp*(273.16_rp/T - 1.0_rp) - 3.56654_rp*Log10(273.16_rp/T) + 0.876793_rp*(1.0_rp-T/273.16_rp) + Log10(6.1071_rp)    
  ei = 100.0_rp*(10.0_rp**ei)   ! hPa --> Pa           
  !
  !***   stat 
  !
  if(e.gt.ew) then 
     ! 
     !***      Value of 255k based on Durant2009
     ! 
     if(T.gt.255.) then
        wphase = 'water'
     else
        wphase = 'ice'
     end if
  else
     if(e.gt.ei) then
        wphase = 'ice'
     else
        wphase = 'vapor'
     end if
  end if
  !
  return
end subroutine water_phase
!
!
!
subroutine get_gamas(i,j,k,gamas)
  !*************************************************************************
  !*
  !*    Computes Gama_s according to:
  !*
  !*    Laminar:    Gama_s = abs( dU/dz )
  !*    Turbulent:  Gama_s = sqrt( 1.3 SijSij)
  !*
  !*************************************************************************
  use KindType
  use Numeric
  use Master
  implicit none
  !
  integer(ip) :: i,j,k
  real   (rp) :: gamas
  real   (rp) :: dXX,dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz, &
       S11,S12,S13,S21,S22,S23,S31,S32,S33,epsilon
  !
  dXX = Hm1(i,j)*dX1
  !                                  
  if(i.eq.1) then
     dudx=(vx(2,j,k)-vx(1,j,k))/dXX
     dvdx=(vy(2,j,k)-vy(1,j,k))/dXX
     dwdx=(vz(2,j,k)-vz(1,j,k))/dXX
  else if(i.eq.nx) then
     dudx=(vx(nx,j,k)-vx(nx-1,j,k))/dXX
     dvdx=(vy(nx,j,k)-vy(nx-1,j,k))/dXX
     dwdx=(vz(nx,j,k)-vz(nx-1,j,k))/dXX
  else
     dudx=(vx(i+1,j,k)-vx(i-1,j,k))/(2.0_rp*dXX)
     dvdx=(vy(i+1,j,k)-vy(i-1,j,k))/(2.0_rp*dXX)
     dwdx=(vz(i+1,j,k)-vz(i-1,j,k))/(2.0_rp*dXX)
  end if
  !                                 
  if(j.eq.1) then
     dudy=(vx(i,2,k)-vx(i,1,k))/dX2
     dvdy=(vy(i,2,k)-vy(i,1,k))/dX2
     dwdy=(vz(i,2,k)-vz(i,1,k))/dX2
  else if(j.eq.ny) then
     dudy=(vx(i,ny,k)-vx(i,ny-1,k))/dX2
     dvdy=(vy(i,ny,k)-vy(i,ny-1,k))/dX2
     dwdy=(vz(i,ny,k)-vz(i,ny-1,k))/dX2
  else
     dudy=(vx(i,j+1,k)-vx(i,j-1,k))/(2.0_rp*dX2)
     dvdy=(vy(i,j+1,k)-vy(i,j-1,k))/(2.0_rp*dX2)
     dwdy=(vz(i,j+1,k)-vz(i,j-1,k))/(2.0_rp*dX2)
  end if
  !                                 
  if(k.eq.1) then
     dudz=(vx(i,j,2)-vx(i,j,1))/dz(1)
     dvdz=(vy(i,j,2)-vy(i,j,1))/dz(1)
     dwdz=(vz(i,j,2)-vz(i,j,1))/dz(1)
  else if(k.eq.nz) then
     dudz=(vx(i,j,nz)-vx(i,j,nz-1))/dz(nz-1)
     dvdz=(vy(i,j,nz)-vy(i,j,nz-1))/dz(nz-1)
     dwdz=(vz(i,j,nz)-vz(i,j,nz-1))/dz(nz-1)
  else
     dudz=(vx(i,j,k+1)-vx(i,j,k-1))/(dz(k-1)+dz(k))
     dvdz=(vy(i,j,k+1)-vy(i,j,k-1))/(dz(k-1)+dz(k))
     dwdz=(vz(i,j,k+1)-vz(i,j,k-1))/(dz(k-1)+dz(k))
  end if
  !      
  S11 = dudx
  S12 = 0.5*(dvdx+dudy)
  S13 = 0.5*(dwdx+dudz)
  S21 = S12
  S22 = dvdy
  S23 = 0.5*(dwdy+dvdz)
  S31 = S13
  S32 = S23
  S33 = dwdz
  !      
  epsilon = S11*S11+S12*S12+S13*S13+ &
       S21*S21+S22*S22+S23*S23+ &
       S31*S31+S32*S32+S33*S33
  !
  gamas = sqrt(1.3*2.0*epsilon)
  !
  !      This is for laminar flow
  !         if(iz.eq.1) then
  !           vup    = sqrt(vx(ix,iy,2)*vx(ix,iy,2)+vy(ix,iy,2)*vy(ix,iy,2))  
  !           vdw    = sqrt(vx(ix,iy,1)*vx(ix,iy,1)+vy(ix,iy,1)*vy(ix,iy,1)) 
  !           gamas = abs( (vup-vdw)/dz(1) )
  !         else if(iz.eq.nz) then
  !           vup    = sqrt(vx(ix,iy,nz  )*vx(ix,iy,nz  )+vy(ix,iy,nz  )*vy(ix,iy,nz  ))  
  !           vdw    = sqrt(vx(ix,iy,nz-1)*vx(ix,iy,nz-1)+vy(ix,iy,nz-1)*vy(ix,iy,nz-1)) 
  !           gamas = abs( (vup-vdw)/dz(nz-1) ) 
  !         else
  !           vup    = sqrt(vx(ix,iy,iz+1)*vx(ix,iy,iz+1)+vy(ix,iy,iz+1)*vy(ix,iy,iz+1))  
  !           vdw    = sqrt(vx(ix,iy,iz-1)*vx(ix,iy,iz-1)+vy(ix,iy,iz-1)*vy(ix,iy,iz-1)) 
  !           gamas = abs( (vup-vdw)/(dz(iz-1)+dz(iz)) )
  !         end if
  !
  return
end subroutine get_gamas
                
      
             
       
