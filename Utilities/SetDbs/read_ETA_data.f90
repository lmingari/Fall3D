subroutine read_ETA_data(fname,time_lag,time,timesec,nx,ny,nz,zlayer,lon,lat,topg,LDU,  &
     u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
  !********************************************************************
  !*
  !*   Reads ETA data in NetCDF format. This routine MUST be called
  !*   after read_ETA_grid
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use netcdf
  use InpOut
  use MathFun
  use ETA_nc
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
  !
  !*** Finds the time itime1 and itime2 and the interpolation factor time_s
  !
  if(nt_ETA.gt.1) then
     found  = .false.
     it_ETA = 0
     do while(.not.found)
        it_ETA = it_ETA + 1
        if( (it_ETA+1).gt.nt_ETA ) call runend('read_ETA : Time not found')
        if( (time_lag+timesec.ge.timesec_ETA(it_ETA  )).and. &
             (time_lag+timesec.le.timesec_ETA(it_ETA+1)) ) then
           found = .true.
           itime1 = it_ETA
           itime2 = it_ETA + 1
           time_s = (timesec_ETA(itime2)-time_lag-timesec)/(timesec_ETA(itime2)-timesec_ETA(itime1))
        end if
     end do
  else
     itime1 = 1
     itime2 = 1
     time_s =1.0_rp
  end if
  !
  write(lulog,10) time,time_ETA(itime1),time_ETA(itime2),time_s
10 format('Interpolating at ',f15.0,/, &
       '        ETA from ',f15.0,' to ',f15.0,' s = ',f6.3)
  !
  !*** Open netCDF file and get ncID
  !
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) call runend('read_ETA : Error in nf90_open')
  !
  !*** First of all reads topograpgy form ETA. This is necessary since interpolation in done
  !*** in terrain following
  !
  allocate(work_ETA1(nx_ETA,ny_ETA,1))
  allocate(work_ETA2(nx_ETA,ny_ETA,1))
  !
  !*** LDU (no given)
  !
  LDU = -99
  !
  !*** HGT:sfc(nx_ETA,ny_ETA) -> topg(nx,ny) (time independent actually)
  !*** NOTE: translation of origin (from 0 to lonmin) is necessary each time a
  !*** variable is read
  !
  if( nf90_inq_varid(ncID,fis_ETA_name,fisID) /= 0) call runend('read_ETA : Error getting fisID')
  if( nf90_get_var(ncID,fisID,work_ETA1,start=(/1,1,itime1/)) /=0 ) call runend('read_ETA: Error reading fis')
  !
  !***  LAND(nx_ETA,ny_ETA) Land Mask (time independent actually)
  !***  0 water 1 land.
  !
  if( nf90_inq_varid(ncID,lmask_ETA_name,lmaskID) /= 0) call runend('read_ETA : Error getting lmaskID')
  if( nf90_get_var(ncID,lmaskID,work_ETA2,start=(/1,1,itime1/)) /=0 ) call runend('read_ETA: Error reading landmask')
  !
  work_ETA1 = work_ETA1*work_ETA2
  !
  call nc_interpola2d(nx_ETA,ny_ETA,nelem_ETA, &
       work_ETA1,lelpo_ETA,lnods_ETA,s_po_ETA,t_po_ETA, &
       nx,ny,npoin_DAT,topg)
  topg_ETA(:,:) = work_ETA1(:,:,1)           ! needed later
  deallocate(work_ETA1)
  deallocate(work_ETA2)
  !
  !*** Reads 3D variables and interpolates
  !
  allocate (work_ETA1(nx_ETA,ny_ETA,nz_ETA))
  allocate (work_ETA2(nx_ETA,ny_ETA,nz_ETA))
  !
  !*** z_ETA(nx_ETA,ny_ETA,nz_ETA) Heights. This is necessary because pressure levels change with time.
  !*** ETA already gives the geopotential in m (i.e. height). The topography is substracted in order to interpolate
  !*** in terrain following
  !
  if( nf90_inq_varid(ncID,fi_ETA_name,fiID) /= 0) call runend('read_ETA : Error getting fiID')
  if( nf90_get_var(ncID,fiID,work_ETA1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ETA: Error reading fi')
  if( nf90_get_var(ncID,fiID,work_ETA2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ETA: Error reading fi')
  !
  z_ETA = time_s*work_ETA1 + (1.0_rp-time_s)*work_ETA2     ! z at time
  !
  do iz = 1,nz_ETA
     z_ETA(:,:,iz) =  z_ETA(:,:,iz) - topg_ETA(:,:)
  end do
  !
  if(MAXVAL(z_ETA(:,:,:)).lt.(zlayer(nz))) &
       call wriwar('Dbs max vertical layer is higher than ETA max layer. Values will be extrapolated')
  !
  !*** TMP(nx_ETA,ny_ETA,nz_ETA) --> T(nx,ny,nz)  Temperature
  !
  if( nf90_inq_varid(ncID,T_ETA_name,tID) /= 0) call runend('read_ETA : Error getting tID')
  if( nf90_get_var(ncID,tID,work_ETA1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ETA: Error reading T')
  if( nf90_get_var(ncID,tID,work_ETA2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ETA: Error reading T')
  !
  work_ETA1 = time_s*work_ETA1 + (1.0_rp-time_s)*work_ETA2   ! T at time
  !
  call nc_interpola3d(nx_ETA,ny_ETA,nz_ETA,nelem_ETA, &
       work_ETA1,z_ETA,lelpo_ETA,lnods_ETA,s_po_ETA,t_po_ETA, &
       nx,ny,nz,npoin_DAT,T,zlayer)
  !
  !*** RH(nx_ETA,ny_ETA,nz_ETA) Relative humidity --> qv(nx,ny,nz)  Specific humidity is computed later
  !
  if( nf90_inq_varid(ncID,Q_ETA_name,qID) /= 0) call runend('read_ETA : Error getting tID')
  if( nf90_get_var(ncID,qID,work_ETA1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ETA: Error reading QV')
  if( nf90_get_var(ncID,qID,work_ETA2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ETA: Error reading QV')
  !
  work_ETA1 = time_s*work_ETA1 + (1.0_rp-time_s)*work_ETA2   ! QV at time
  !
  call nc_interpola3d(nx_ETA,ny_ETA,nz_ETA,nelem_ETA, &
       work_ETA1,z_ETA,lelpo_ETA,lnods_ETA,s_po_ETA,t_po_ETA, &
       nx,ny,nz,npoin_DAT,qv,zlayer)
  !
  !*** U(nx_ETA,ny_ETA,nz_ETA) --> u(nx,ny,nz)  Velocity
  !
  if( nf90_inq_varid(ncID,u_ETA_name,uID) /= 0) call runend('read_ETA : Error getting uID')
  if( nf90_get_var(ncID,uID,work_ETA1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ETA: Error reading u')
  if( nf90_get_var(ncID,uID,work_ETA2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ETA: Error reading u')
  !
  work_ETA1 = time_s*work_ETA1 + (1.0_rp-time_s)*work_ETA2   ! U at time
  !
  call nc_interpola3d(nx_ETA,ny_ETA,nz_ETA,nelem_ETA, &
       work_ETA1,z_ETA,lelpo_ETA,lnods_ETA,s_po_ETA,t_po_ETA, &
       nx,ny,nz,npoin_DAT,u,zlayer)
  !
  !*** V(nx_ETA,ny_ETA,nz_ETA) --> v(nx,ny,nz)  Velocity
  !
  if( nf90_inq_varid(ncID,v_ETA_name,vID) /= 0) call runend('read_ETA : Error getting vID')
  if( nf90_get_var(ncID,vID,work_ETA1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ETA: Error reading v')
  if( nf90_get_var(ncID,vID,work_ETA2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ETA: Error reading v')
  !
  work_ETA1 = time_s*work_ETA1 + (1.0_rp-time_s)*work_ETA2   ! V at time
  !
  call nc_interpola3d(nx_ETA,ny_ETA,nz_ETA,nelem_ETA, &
       work_ETA1,z_ETA,lelpo_ETA,lnods_ETA,s_po_ETA,t_po_ETA, &
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
  !*** W(nx_ETA,ny_ETA,nz_ETA) --> w(nx,ny,nz) Omega velocity
  !
  if( nf90_inq_varid(ncID,w_ETA_name,wID) /= 0) call runend('read_ETA : Error getting wID')
  if( nf90_get_var(ncID,wID,work_ETA1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ETA: Error reading w')
  if( nf90_get_var(ncID,wID,work_ETA2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ETA: Error reading w')
  !
  work_ETA1 = time_s*work_ETA1 + (1.0_rp-time_s)*work_ETA2            ! Omega at time
  do iz = 1,nz_ETA
     work_ETA1(:,:,iz) = work_ETA1(:,:,iz)/p_ETA(iz)                ! Omega/P at time
  end do
  !
  call nc_interpola3d(nx_ETA,ny_ETA,nz_ETA,nelem_ETA, &
       work_ETA1,z_ETA,lelpo_ETA,lnods_ETA,s_po_ETA,t_po_ETA, &
       nx,ny,nz,npoin_DAT,w,zlayer)
  !
  !*** p(nx,ny,nz) pressure
  !
  do iz = 1,nz_ETA
     work_ETA1(:,:,iz) = p_ETA(iz)                                   ! P at time
  end do
  !
  call nc_interpola3d(nx_ETA,ny_ETA,nz_ETA,nelem_ETA, &
       work_ETA1,z_ETA,lelpo_ETA,lnods_ETA,s_po_ETA,t_po_ETA, &
       nx,ny,nz,npoin_DAT,p,zlayer)
  !
  !*** Release memory
  !
  deallocate(work_ETA1)
  deallocate(work_ETA2)
  !
  !*** Other derived 3D variables and transformations (already in the dbs grid!)
  !
  allocate(work_ETA1(nx,ny,nz))
  allocate(work_ETA2(nx,ny,nz))
  !
  !*** qv(nx,ny,nz) Specific humidity (ETA gives relative humidity, by now stored in qv)
  !*** Following Bolton (1980)
  !***      es = 6.112 * exp((17.67 * T)/(T + 243.5));  saturation vapor pressure in mb, T in C
  !***      e  = es * (RH/100.0);                       vapor pressure in mb
  !***      q  = (0.622 * e)/(p - (0.378 * e))          specific humidity in kg/kg, p in mb
  !
  work_ETA1(:,:,:) = 6.112*exp( (17.67*(T(:,:,:)-273.15))/(243.5+(T(:,:,:)-273.15)) )  ! Saturation Pressure in mb
  work_ETA1(:,:,:) = work_ETA1(:,:,:)*qv(:,:,:)/100.0_rp                               ! Vapor      Pressure in mb
  qv(:,:,:) = 0.622*work_ETA1(:,:,:)/((p(:,:,:)/101_rp)-0.378*work_ETA1(:,:,:))        ! Specific   humidity in kg/kg
  !
  !*** Tv(nx,ny,nz) virtual temperature
  !
  work_ETA1 = qv(:,:,:)/(1.0_rp-qv(:,:,:))      ! Mixing ratio (humidity ratio)
  !
  Tv(:,:,:) = T(:,:,:)*(1.0_rp + 1.6077_rp*work_ETA1(:,:,:))/(1.0_rp+work_ETA1(:,:,:))
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
  deallocate(work_ETA1)
  deallocate(work_ETA2)
  !
  !***  Closes
  !
  if( nf90_close(ncID) /= 0) call runend('read_ETA : Error in nf90_close')
  !
  !***  Computes other variables not given by ETA
  !
  !
  !*** PBL(nx,ny) boundary layer height. Not given by ETA. Value assumed
  !
  pblh(:,:) = 1500.0_rp
  !
  !*** ust(nx,ny) and L(nx,ny). Estimated
  !
  do iy = 1,ny
     do ix = 1,nx
        call get_par_ABL(1.0_rp,zlayer(2),T(ix,iy,1),T(ix,iy,2),P(ix,iy,1),P(ix,iy,2), &
             u(ix,iy,1),v(ix,iy,1),ust(ix,iy),L(ix,iy))
     end do
  end do
  !
  return
end subroutine read_ETA_data
