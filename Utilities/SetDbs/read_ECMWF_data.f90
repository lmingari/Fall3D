subroutine read_ECMWF_data(fname,time_lag,time,timesec,nx,ny,nz,zlayer,lon,lat,topg,LDU,  &
     u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
  !********************************************************************
  !*
  !*   Reads ECMWF data in NetCDF format. This routine MUST be called
  !*   after read_ECMWF_grid
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use netcdf
  use InpOut
  use MathFun
  use ECMWF_nc
  use wind_rotation
  implicit none
  !
  character(len=*) :: fname
  integer(ip)      :: nx,ny,nz
  real(rp)         :: time_lag,time
  real(rp)         :: timesec                     ! time in sec after 0000UTC
  real(rp)         :: zlayer(nz),lat(nx,ny),lon(nx,ny),topg(nx,ny),LDU(nx,ny)
  real(rp)         :: u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz),T(nx,ny,nz),Tv(nx,ny,nz), &
       Tp(nx,ny,nz),qv(nx,ny,nz),ro(nx,ny,nz),p(nx,ny,nz)
  real(rp)         :: pblh(nx,ny),ust(nx,ny),L(nx,ny)
  !
  logical          :: found
  integer          :: itime1,itime2,ix,iy,iz
  real(rp)         :: time_s
  real(rp), allocatable :: usfc(:,:),vsfc(:,:),Tsfc(:,:),Psfc(:,:)
  !
  !*** Finds the time itime1 and itime2 and the interpolation factor time_s
  !
  if(nt_ECMWF.gt.1) then
     found  = .false.
     it_ECMWF = 0
     do while(.not.found)
        it_ECMWF = it_ECMWF + 1
        if( (it_ECMWF+1).gt.nt_ECMWF ) call runend('read_ECMWF : Time not found')
        if( (time_lag+timesec.ge.timesec_ECMWF(it_ECMWF  )).and. &
             (time_lag+timesec.le.timesec_ECMWF(it_ECMWF+1)) ) then
           found = .true.
           itime1 = it_ECMWF
           itime2 = it_ECMWF + 1
           time_s = (timesec_ECMWF(itime2)-time_lag-timesec)/(timesec_ECMWF(itime2)-timesec_ECMWF(itime1))
        end if
     end do
  else
     itime1 = 1
     itime2 = 1
     time_s =1.0_rp
  end if
  !
  write(lulog,10) time,time_ECMWF(itime1),time_ECMWF(itime2),time_s
10 format('Interpolating at ',f15.0,/, &
       '        ECMWF from ',f15.0,' to ',f15.0,' s = ',f6.3)
  !
  !*** Open netCDF file and get ncID
  !
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) call runend('read_ECMWF : Error in nf90_open')
  !
  !*** First of all reads topograpgy form ECMWF. This is necessary since interpolation in done
  !*** in terrain following
  !
  allocate(work_ECMWF1(nx_ECMWF,ny_ECMWF,1))
  allocate(work_ECMWF2(nx_ECMWF,ny_ECMWF,1))
  !
  !*** LDU (no given)
  !
  LDU = -99
  !
  !*** Z:sfc(nx_ECMWF,ny_ECMWF) Geopotential FI at surface -> topg(nx,ny) (time independent actually)
  !
  if( nf90_inq_varid(ncID,fis_ECMWF_name,fisID) /= 0) call runend('read_ECMWF : Error getting fisID')
  if( nf90_get_var(ncID,fisID,work_ECMWF1,start=(/1,1,itime1/)) /=0 ) call runend('read_ECMWF: Error reading fis')
  !
  work_ECMWF1 = work_ECMWF1 / 9.81_rp            ! topg = (fis)/g
  !
  !***  LAND(nx_ECMWF,ny_ECMWF) Land Mask (time independent actually)
  !***  0 water 1 land.
  !
  if( nf90_inq_varid(ncID,lmask_ECMWF_name,lmaskID) /= 0) call runend('read_ECMWF : Error getting lmaskID')
  if( nf90_get_var(ncID,lmaskID,work_ECMWF2,start=(/1,1,itime1/)) /=0 ) call runend('read_ECMWF: Error reading landmask')
  !
  work_ECMWF1 = work_ECMWF1*work_ECMWF2
  !
  call nc_interpola2d(nx_ECMWF,ny_ECMWF,nelem_ECMWF, &
       work_ECMWF1,lelpo_ECMWF,lnods_ECMWF,s_po_ECMWF,t_po_ECMWF, &
       nx,ny,npoin_DAT,topg)
  topg_ECMWF(:,:) = work_ECMWF1(:,:,1)           ! needed later
  deallocate(work_ECMWF1)
  deallocate(work_ECMWF2)
  !
  !*** Reads 3D variables and interpolates
  !
  allocate (work_ECMWF1(nx_ECMWF,ny_ECMWF,nz_ECMWF))
  allocate (work_ECMWF2(nx_ECMWF,ny_ECMWF,nz_ECMWF))
  !
  !*** z_ECMWF(nx_ECMWF,ny_ECMWF,nz_ECMWF) Heights. This is necessary because pressure levels change with time.
  !*** Heights are computed from the geopotential FI. The topography is substracted in order to interpolate
  !*** in terrain following
  !
  if( nf90_inq_varid(ncID,fi_ECMWF_name,fiID) /= 0) call runend('read_ECMWF : Error getting fiID')
  if( nf90_get_var(ncID,fiID,work_ECMWF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ECMWF: Error reading fi')
  work_ECMWF1 = work_ECMWF1 / 9.81_rp            ! topg = (fis)/g at time1
  if( nf90_get_var(ncID,fiID,work_ECMWF2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ECMWF: Error reading fi')
  work_ECMWF2 = work_ECMWF2 / 9.81_rp            ! topg = (fis)/g at time2
  !
  z_ECMWF = time_s*work_ECMWF1 + (1.0_rp-time_s)*work_ECMWF2     ! z at time
  do iz = 1,nz_ECMWF
     z_ECMWF(:,:,iz) =  z_ECMWF(:,:,iz) - topg_ECMWF(:,:)
  end do
  !
  if(MAXVAL(z_ECMWF(:,:,:)).lt.(zlayer(nz))) &
       call wriwar('Dbs max vertical layer is higher than ECMWF max layer. Values will be extrapolated')
  !
  !*** TMP(nx_ECMWF,ny_ECMWF,nz_ECMWF) --> T(nx,ny,nz)  Temperature
  !
  if( nf90_inq_varid(ncID,T_ECMWF_name,tID) /= 0) call runend('read_ECMWF : Error getting tID')
  if( nf90_get_var(ncID,tID,work_ECMWF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ECMWF: Error reading T')
  if( nf90_get_var(ncID,tID,work_ECMWF2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ECMWF: Error reading T')
  !
  work_ECMWF1 = time_s*work_ECMWF1 + (1.0_rp-time_s)*work_ECMWF2   ! T at time
  !
  call nc_interpola3d(nx_ECMWF,ny_ECMWF,nz_ECMWF,nelem_ECMWF, &
       work_ECMWF1,z_ECMWF,lelpo_ECMWF,lnods_ECMWF,s_po_ECMWF,t_po_ECMWF, &
       nx,ny,nz,npoin_DAT,T,zlayer)
  !
  !*** RH(nx_ECMWF,ny_ECMWF,nz_ECMWF) Relative humidity --> qv(nx,ny,nz)  Specific humidity is computed later
  !
  if( nf90_inq_varid(ncID,Q_ECMWF_name,qID) /= 0) call runend('read_ECMWF : Error getting tID')
  if( nf90_get_var(ncID,qID,work_ECMWF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ECMWF: Error reading QV')
  if( nf90_get_var(ncID,qID,work_ECMWF2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ECMWF: Error reading QV')
  !
  work_ECMWF1 = time_s*work_ECMWF1 + (1.0_rp-time_s)*work_ECMWF2   ! QV at time
  !
  call nc_interpola3d(nx_ECMWF,ny_ECMWF,nz_ECMWF,nelem_ECMWF, &
       work_ECMWF1,z_ECMWF,lelpo_ECMWF,lnods_ECMWF,s_po_ECMWF,t_po_ECMWF, &
       nx,ny,nz,npoin_DAT,qv,zlayer)
  !
  !*** U(nx_ECMWF,ny_ECMWF,nz_ECMWF) --> u(nx,ny,nz)  Velocity
  !
  if( nf90_inq_varid(ncID,u_ECMWF_name,uID) /= 0) call runend('read_ECMWF : Error getting uID')
  if( nf90_get_var(ncID,uID,work_ECMWF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ECMWF: Error reading u')
  if( nf90_get_var(ncID,uID,work_ECMWF2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ECMWF: Error reading u')
  !
  work_ECMWF1 = time_s*work_ECMWF1 + (1.0_rp-time_s)*work_ECMWF2   ! U at time
  !
  call nc_interpola3d(nx_ECMWF,ny_ECMWF,nz_ECMWF,nelem_ECMWF, &
       work_ECMWF1,z_ECMWF,lelpo_ECMWF,lnods_ECMWF,s_po_ECMWF,t_po_ECMWF, &
       nx,ny,nz,npoin_DAT,u,zlayer)
  !
  !*** V(nx_ECMWF,ny_ECMWF,nz_ECMWF) --> v(nx,ny,nz)  Velocity
  !
  if( nf90_inq_varid(ncID,v_ECMWF_name,vID) /= 0) call runend('read_ECMWF : Error getting vID')
  if( nf90_get_var(ncID,vID,work_ECMWF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ECMWF: Error reading v')
  if( nf90_get_var(ncID,vID,work_ECMWF2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ECMWF: Error reading v')
  !
  work_ECMWF1 = time_s*work_ECMWF1 + (1.0_rp-time_s)*work_ECMWF2   ! V at time
  !
  call nc_interpola3d(nx_ECMWF,ny_ECMWF,nz_ECMWF,nelem_ECMWF, &
       work_ECMWF1,z_ECMWF,lelpo_ECMWF,lnods_ECMWF,s_po_ECMWF,t_po_ECMWF, &
       nx,ny,nz,npoin_DAT,v,zlayer)
  !
  ! Rotate wind if request
  if(rotate_wind) then
     allocate(unew(nx,ny))
     allocate(vnew(nx,ny))
     do iz=1,nz
        do iy=1,ny
           do ix = 1,nx
              call rotate2d(u(ix,iy,iz),v(ix,iy,iz),unew(ix,iy),vnew(ix,iy),rotation_angle)
           end do
        end do
        u(:,:,iz) = unew      ! Copy
        v(:,:,iz) = vnew      ! Copy
     end do
     deallocate(unew)
     deallocate(vnew)
  end if
  !
  !*** W(nx_ECMWF,ny_ECMWF,nz_ECMWF) --> w(nx,ny,nz) Omega velocity
  !
  if( nf90_inq_varid(ncID,w_ECMWF_name,wID) /= 0) call runend('read_ECMWF : Error getting wID')
  if( nf90_get_var(ncID,wID,work_ECMWF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ECMWF: Error reading w')
  if( nf90_get_var(ncID,wID,work_ECMWF2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ECMWF: Error reading w')
  !
  work_ECMWF1 = time_s*work_ECMWF1 + (1.0_rp-time_s)*work_ECMWF2            ! Omega at time
  do iz = 1,nz_ECMWF
     work_ECMWF1(:,:,iz) = work_ECMWF1(:,:,iz)/p_ECMWF(iz)                ! Omega/P at time
  end do
  !
  call nc_interpola3d(nx_ECMWF,ny_ECMWF,nz_ECMWF,nelem_ECMWF, &
       work_ECMWF1,z_ECMWF,lelpo_ECMWF,lnods_ECMWF,s_po_ECMWF,t_po_ECMWF, &
       nx,ny,nz,npoin_DAT,w,zlayer)
  !
  !*** p(nx,ny,nz) pressure
  !
  do iz = 1,nz_ECMWF
     work_ECMWF1(:,:,iz) = p_ECMWF(iz)                                   ! P at time
  end do
  !
  call nc_interpola3d(nx_ECMWF,ny_ECMWF,nz_ECMWF,nelem_ECMWF, &
       work_ECMWF1,z_ECMWF,lelpo_ECMWF,lnods_ECMWF,s_po_ECMWF,t_po_ECMWF, &
       nx,ny,nz,npoin_DAT,p,zlayer)
  !
  !*** Release memory
  !
  deallocate(work_ECMWF1)
  deallocate(work_ECMWF2)
  !
  !*** Rest of 2D scalar variables
  !
  allocate(work_ECMWF1(nx_ECMWF,ny_ECMWF,1))
  allocate(work_ECMWF2(nx_ECMWF,ny_ECMWF,1))
  !
  !*** HPBL(nx_ECMWF,ny_ECMWF) --> PBLH(nx,ny) Planetary boundary layer height
  !
  if( nf90_inq_varid(ncID,pblh_ECMWF_name,pblhID) /= 0) call runend('read_ECMWF : Error getting pblhID')
  if( nf90_get_var(ncID,pblhID,work_ECMWF1,start=(/1,1,itime1/)) /=0 ) call runend('read_ECMWF: Error reading pblh')
  if( nf90_get_var(ncID,pblhID,work_ECMWF2,start=(/1,1,itime2/)) /=0 ) call runend('read_ECMWF: Error reading pblh')
  !
  work_ECMWF1 = time_s*work_ECMWF1 + (1.0_rp-time_s)*work_ECMWF2   ! pblh box
  !
  call nc_interpola2d(nx_ECMWF,ny_ECMWF,nelem_ECMWF, &
       work_ECMWF1,lelpo_ECMWF,lnods_ECMWF,s_po_ECMWF,t_po_ECMWF, &
       nx,ny,npoin_DAT,pblh)
  !
  deallocate(work_ECMWF1)
  deallocate(work_ECMWF2)
  !
  !*** Other derived 3D variables and transformations (already in the dbs grid!)
  !
  allocate(work_ECMWF1(nx,ny,nz))
  allocate(work_ECMWF2(nx,ny,nz))
  !
  !*** qv(nx,ny,nz) Specific humidity (ECMWF gives relative humidity, by now stored in qv)
  !*** Following Bolton (1980)
  !***      es = 6.112 * exp((17.67 * T)/(T + 243.5));  saturation vapor pressure in mb, T in C
  !***      e  = es * (RH/100.0);                       vapor pressure in mb
  !***      q  = (0.622 * e)/(p - (0.378 * e))          specific humidity in kg/kg, p in mb
  !
  work_ECMWF1(:,:,:) = 6.112*exp( (17.67*(T(:,:,:)-273.15))/(243.5+(T(:,:,:)-273.15)) )  ! Saturation Pressure in mb
  work_ECMWF1(:,:,:) = work_ECMWF1(:,:,:)*qv(:,:,:)/100.0_rp                               ! Vapor      Pressure in mb
  qv(:,:,:) = 0.622*work_ECMWF1(:,:,:)/((p(:,:,:)/101_rp)-0.378*work_ECMWF1(:,:,:))        ! Specific   humidity in kg/kg
  !
  !*** Tv(nx,ny,nz) virtual temperature
  !
  work_ECMWF1 = qv(:,:,:)/(1.0_rp-qv(:,:,:))      ! Mixing ratio (humidity ratio)
  !
  Tv(:,:,:) = T(:,:,:)*(1.0_rp + 1.6077_rp*work_ECMWF1(:,:,:))/(1.0_rp+work_ECMWF1(:,:,:))
  !
  !*** Tp(nx,ny,nz) potential temperature  Theta=T*(1bar/P))**(R/cp)
  !
  Tp(:,:,:) = T(:,:,:)*(1.01e5_rp/p(:,:,:))**(0.285_rp)
  !
  !*** ro(nx,ny,nz) density. Gas law using virtual temperature (account for water in air)
  !
  ro(:,:,:) = p(:,:,:)/(287.06*Tv(:,:,:))
  !
  !*** w(nx,ny,nz) vertical velocity  w = - omega/(rho*g) = - omega*R*Tv/(P*g)
  !
  w(:,:,:) =  - w(:,:,:)*287.06_rp*Tv(:,:,:)/9.81_rp
  !
  deallocate(work_ECMWF1)
  deallocate(work_ECMWF2)
  !
  !***  Computes other variables not given by ECMWF
  !
  !
  !*** ust(nx,ny) and L(nx,ny).
  !*** For this, use values at surface (u,v,T and P)
  !
  allocate(usfc(nx,ny))
  allocate(vsfc(nx,ny))
  allocate(Tsfc(nx,ny))
  allocate(Psfc(nx,ny))
  allocate(work_ECMWF1(nx_ECMWF,ny_ECMWF,1))
  allocate(work_ECMWF2(nx_ECMWF,ny_ECMWF,1))
  !
  if( nf90_inq_varid(ncID,TRIM(us_ECMWF_name),usID) /= 0) call runend('read_ECMWF : Error getting usID')
  if( nf90_get_var(ncID,usID,work_ECMWF1,start=(/1,1,itime1/)) /=0 ) call runend('read_ECMWF: Error reading us')
  if( nf90_get_var(ncID,usID,work_ECMWF2,start=(/1,1,itime2/)) /=0 ) call runend('read_ECMWF: Error reading us')
  work_ECMWF1 = time_s*work_ECMWF1 + (1.0_rp-time_s)*work_ECMWF2            ! u:sfc at time
  call nc_interpola2d(nx_ECMWF,ny_ECMWF,nelem_ECMWF, &
       work_ECMWF1,lelpo_ECMWF,lnods_ECMWF,s_po_ECMWF,t_po_ECMWF, &
       nx,ny,npoin_DAT,usfc)
  !
  if( nf90_inq_varid(ncID,vs_ECMWF_name,vsID) /= 0) call runend('read_ECMWF : Error getting vsID')
  if( nf90_get_var(ncID,vsID,work_ECMWF1,start=(/1,1,itime1/)) /=0 ) call runend('read_ECMWF: Error reading vs')
  if( nf90_get_var(ncID,vsID,work_ECMWF2,start=(/1,1,itime2/)) /=0 ) call runend('read_ECMWF: Error reading vs')
  work_ECMWF1 = time_s*work_ECMWF1 + (1.0_rp-time_s)*work_ECMWF2            ! v:sfc at time
  call nc_interpola2d(nx_ECMWF,ny_ECMWF,nelem_ECMWF, &
       work_ECMWF1,lelpo_ECMWF,lnods_ECMWF,s_po_ECMWF,t_po_ECMWF, &
       nx,ny,npoin_DAT,vsfc)
  !
  ! Rotate surface wind if request
  if(rotate_wind) then
     allocate(unew(nx,ny))
     allocate(vnew(nx,ny))
     do iy=1,ny
        do ix = 1,nx
           call rotate2d(usfc(ix,iy),vsfc(ix,iy),unew(ix,iy),vnew(ix,iy),rotation_angle)
        end do
     end do
     usfc = unew      ! Copy
     vsfc = vnew      ! Copy
     deallocate(unew)
     deallocate(vnew)
  end if
  !
  if( nf90_inq_varid(ncID,Ts_ECMWF_name,TsID) /= 0) call runend('read_ECMWF : Error getting TsID')
  if( nf90_get_var(ncID,TsID,work_ECMWF1,start=(/1,1,itime1/)) /=0 ) call runend('read_ECMWF: Error reading Ts')
  if( nf90_get_var(ncID,TsID,work_ECMWF2,start=(/1,1,itime2/)) /=0 ) call runend('read_ECMWF: Error reading Ts')
  work_ECMWF1 = time_s*work_ECMWF1 + (1.0_rp-time_s)*work_ECMWF2            ! T:sfc at time
  call nc_interpola2d(nx_ECMWF,ny_ECMWF,nelem_ECMWF, &
       work_ECMWF1,lelpo_ECMWF,lnods_ECMWF,s_po_ECMWF,t_po_ECMWF, &
       nx,ny,npoin_DAT,Tsfc)
  !
  do iy = 1,ny
     do ix = 1,nx
        call get_par_ABL(1.0_rp,zlayer(2),Tsfc(ix,iy),T(ix,iy,2),P(ix,iy,1),P(ix,iy,2), &
             usfc(ix,iy),vsfc(ix,iy),ust(ix,iy),L(ix,iy))
     end do
  end do
  !
  deallocate(usfc)
  deallocate(vsfc)
  deallocate(Tsfc)
  deallocate(Psfc)
  deallocate(work_ECMWF1)
  deallocate(work_ECMWF2)
  !
  !***  Closes
  !
  if( nf90_close(ncID) /= 0) call runend('read_ECMWF : Error in nf90_close')
  !
  return
end subroutine read_ECMWF_data
