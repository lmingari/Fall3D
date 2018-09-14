subroutine read_NCEP2_data(fname,time_lag,time,timesec,nx,ny,nz,zlayer,lon,lat,topg,LDU,  &
     u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
  !********************************************************************
  !*
  !*   Reads NCEP2 data in NetCDF format. This routine MUST be called
  !*   after read_NCEP2_grid
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use netcdf
  use InpOut
  use MathFun
  use NCEP2_nc
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
  if(nt_NCEP2.gt.1) then
     found  = .false.
     it_NCEP2 = 0
     do while(.not.found)
        it_NCEP2 = it_NCEP2 + 1
        if( (it_NCEP2+1).gt.nt_NCEP2 ) call runend('read_NCEP2 : Time not found')
        if( (time_lag+timesec.ge.timesec_NCEP2(it_NCEP2  )).and. &
             (time_lag+timesec.le.timesec_NCEP2(it_NCEP2+1)) ) then
           found = .true.
           itime1 = it_NCEP2
           itime2 = it_NCEP2 + 1
           time_s = (timesec_NCEP2(itime2)-time_lag-timesec)/(timesec_NCEP2(itime2)-timesec_NCEP2(itime1))
        end if
     end do
  else
     itime1 = 1
     itime2 = 1
     time_s =1.0_rp
  end if
  !
  write(lulog,10) time,time_NCEP2(itime1),time_NCEP2(itime2),time_s
10 format('Interpolating at ',f15.0,/, &
       '        NCEP2 from ',f15.0,' to ',f15.0,' s = ',f6.3)
  !
  !*** Open netCDF file and get ncID
  !
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) call runend('read_NCEP2 : Error in nf90_open')
  !
  !*** First of all reads topograpgy form NCEP2. This is necessary since interpolation in done
  !*** in terrain following
  !
  allocate(work_NCEP21(nx_NCEP2,ny_NCEP2,1))
  allocate(work_NCEP22(nx_NCEP2,ny_NCEP2,1))
  !
  !*** LDU (no given)
  !
  LDU = -99
  !
  !*** HGT:surface(nx_NCEP2,ny_NCEP2) -> topg(nx,ny) (time independent actually)
  !*** NOTE: translation of origin (from 0 to lonmin) is necessary each time a
  !*** variable is read
  !
  if( nf90_inq_varid(ncID,fis_NCEP2_name,fisID) /= 0) call runend('read_NCEP2 : Error getting fisID')
  if( nf90_get_var(ncID,fisID,work_NCEP21,start=(/1,1,itime1/)) /=0 ) call runend('read_NCEP2: Error reading fis')
  call read_NCEP2_translate(work_NCEP21,nx_NCEP2,ny_NCEP2,1)
  !
  call nc_interpola2d(nx_NCEP2,ny_NCEP2,nelem_NCEP2, &
       work_NCEP21,lelpo_NCEP2,lnods_NCEP2,s_po_NCEP2,t_po_NCEP2, &
       nx,ny,npoin_DAT,topg)
  topg_NCEP2(:,:) = work_NCEP21(:,:,1)           ! needed later
  deallocate(work_NCEP21)
  deallocate(work_NCEP22)
  !
  !*** Reads 3D variables and interpolates
  !
  allocate (work_NCEP21(nx_NCEP2,ny_NCEP2,nz_NCEP2))
  allocate (work_NCEP22(nx_NCEP2,ny_NCEP2,nz_NCEP2))
  !
  !*** z_NCEP2(nx_NCEP2,ny_NCEP2,nz_NCEP2) Heights. This is necessary because pressure levels change with time.
  !*** NCEP2 already gives the geopotential in m (i.e. height). The topography is substracted in order to interpolate
  !*** in terrain following
  !
  if( nf90_inq_varid(ncID,fi_NCEP2_name,fiID) /= 0) call runend('read_NCEP2 : Error getting fiID')
  if( nf90_get_var(ncID,fiID,work_NCEP21,start=(/1,1,1,itime1/)) /=0 ) call runend('read_NCEP2: Error reading fi')
  if( nf90_get_var(ncID,fiID,work_NCEP22,start=(/1,1,1,itime2/)) /=0 ) call runend('read_NCEP2: Error reading fi')
  !
  call read_NCEP2_translate(work_NCEP21,nx_NCEP2,ny_NCEP2,nz_NCEP2)
  call read_NCEP2_translate(work_NCEP22,nx_NCEP2,ny_NCEP2,nz_NCEP2)
  z_NCEP2 = time_s*work_NCEP21 + (1.0_rp-time_s)*work_NCEP22     ! z at time
  !
  do iz = 1,nz_NCEP2
     z_NCEP2(:,:,iz) =  z_NCEP2(:,:,iz) - topg_NCEP2(:,:)
  end do
  !
  if(MAXVAL(z_NCEP2(:,:,:)).lt.(zlayer(nz))) &
       call wriwar('Dbs max vertical layer is higher than NCEP2 max layer. Values will be extrapolated')
  !
  !*** TMP(nx_NCEP2,ny_NCEP2,nz_NCEP2) --> T(nx,ny,nz)  Temperature
  !
  if( nf90_inq_varid(ncID,T_NCEP2_name,tID) /= 0) call runend('read_NCEP2 : Error getting tID')
  if( nf90_get_var(ncID,tID,work_NCEP21,start=(/1,1,1,itime1/)) /=0 ) call runend('read_NCEP2: Error reading T')
  if( nf90_get_var(ncID,tID,work_NCEP22,start=(/1,1,1,itime2/)) /=0 ) call runend('read_NCEP2: Error reading T')
  !
  call read_NCEP2_translate(work_NCEP21,nx_NCEP2,ny_NCEP2,nz_NCEP2)
  call read_NCEP2_translate(work_NCEP22,nx_NCEP2,ny_NCEP2,nz_NCEP2)
  work_NCEP21 = time_s*work_NCEP21 + (1.0_rp-time_s)*work_NCEP22   ! T at time
  !
  call nc_interpola3d(nx_NCEP2,ny_NCEP2,nz_NCEP2,nelem_NCEP2, &
       work_NCEP21,z_NCEP2,lelpo_NCEP2,lnods_NCEP2,s_po_NCEP2,t_po_NCEP2, &
       nx,ny,nz,npoin_DAT,T,zlayer)
  !
  !*** RH(nx_NCEP2,ny_NCEP2,nz_NCEP2) Relative humidity --> qv(nx,ny,nz)  Specific humidity is computed later
  !
  if( nf90_inq_varid(ncID,Q_NCEP2_name,qID) /= 0) call runend('read_NCEP2 : Error getting tID')
  if( nf90_get_var(ncID,qID,work_NCEP21,start=(/1,1,1,itime1/)) /=0 ) call runend('read_NCEP2: Error reading QV')
  if( nf90_get_var(ncID,qID,work_NCEP22,start=(/1,1,1,itime2/)) /=0 ) call runend('read_NCEP2: Error reading QV')
  !
  call read_NCEP2_translate(work_NCEP21,nx_NCEP2,ny_NCEP2,nz_NCEP2)
  call read_NCEP2_translate(work_NCEP22,nx_NCEP2,ny_NCEP2,nz_NCEP2)
  work_NCEP21 = time_s*work_NCEP21 + (1.0_rp-time_s)*work_NCEP22   ! QV at time
  !
  call nc_interpola3d(nx_NCEP2,ny_NCEP2,nz_NCEP2,nelem_NCEP2, &
       work_NCEP21,z_NCEP2,lelpo_NCEP2,lnods_NCEP2,s_po_NCEP2,t_po_NCEP2, &
       nx,ny,nz,npoin_DAT,qv,zlayer)
  !
  !*** U(nx_NCEP2,ny_NCEP2,nz_NCEP2) --> u(nx,ny,nz)  Velocity
  !
  if( nf90_inq_varid(ncID,u_NCEP2_name,uID) /= 0) call runend('read_NCEP2 : Error getting uID')
  if( nf90_get_var(ncID,uID,work_NCEP21,start=(/1,1,1,itime1/)) /=0 ) call runend('read_NCEP2: Error reading u')
  if( nf90_get_var(ncID,uID,work_NCEP22,start=(/1,1,1,itime2/)) /=0 ) call runend('read_NCEP2: Error reading u')
  !
  call read_NCEP2_translate(work_NCEP21,nx_NCEP2,ny_NCEP2,nz_NCEP2)
  call read_NCEP2_translate(work_NCEP22,nx_NCEP2,ny_NCEP2,nz_NCEP2)
  work_NCEP21 = time_s*work_NCEP21 + (1.0_rp-time_s)*work_NCEP22   ! U at time
  !
  call nc_interpola3d(nx_NCEP2,ny_NCEP2,nz_NCEP2,nelem_NCEP2, &
       work_NCEP21,z_NCEP2,lelpo_NCEP2,lnods_NCEP2,s_po_NCEP2,t_po_NCEP2, &
       nx,ny,nz,npoin_DAT,u,zlayer)
  !
  !*** V(nx_NCEP2,ny_NCEP2,nz_NCEP2) --> v(nx,ny,nz)  Velocity
  !
  if( nf90_inq_varid(ncID,v_NCEP2_name,vID) /= 0) call runend('read_NCEP2 : Error getting vID')
  if( nf90_get_var(ncID,vID,work_NCEP21,start=(/1,1,1,itime1/)) /=0 ) call runend('read_NCEP2: Error reading v')
  if( nf90_get_var(ncID,vID,work_NCEP22,start=(/1,1,1,itime2/)) /=0 ) call runend('read_NCEP2: Error reading v')
  !
  call read_NCEP2_translate(work_NCEP21,nx_NCEP2,ny_NCEP2,nz_NCEP2)
  call read_NCEP2_translate(work_NCEP22,nx_NCEP2,ny_NCEP2,nz_NCEP2)
  work_NCEP21 = time_s*work_NCEP21 + (1.0_rp-time_s)*work_NCEP22   ! V at time
  !
  call nc_interpola3d(nx_NCEP2,ny_NCEP2,nz_NCEP2,nelem_NCEP2, &
       work_NCEP21,z_NCEP2,lelpo_NCEP2,lnods_NCEP2,s_po_NCEP2,t_po_NCEP2, &
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
  !*** W(nx_NCEP2,ny_NCEP2,nz_NCEP2) --> w(nx,ny,nz) Omega velocity
  !
  if( nf90_inq_varid(ncID,w_NCEP2_name,wID) /= 0) call runend('read_NCEP2 : Error getting wID')
  if( nf90_get_var(ncID,wID,work_NCEP21,start=(/1,1,1,itime1/)) /=0 ) call runend('read_NCEP2: Error reading w')
  if( nf90_get_var(ncID,wID,work_NCEP22,start=(/1,1,1,itime2/)) /=0 ) call runend('read_NCEP2: Error reading w')
  !
  call read_NCEP2_translate(work_NCEP21,nx_NCEP2,ny_NCEP2,nz_NCEP2)
  call read_NCEP2_translate(work_NCEP22,nx_NCEP2,ny_NCEP2,nz_NCEP2)
  work_NCEP21 = time_s*work_NCEP21 + (1.0_rp-time_s)*work_NCEP22            ! Omega at time
  do iz = 1,nz_NCEP2
     work_NCEP21(:,:,iz) = work_NCEP21(:,:,iz)/p_NCEP2(iz)                ! Omega/P at time
  end do
  !
  call nc_interpola3d(nx_NCEP2,ny_NCEP2,nz_NCEP2,nelem_NCEP2, &
       work_NCEP21,z_NCEP2,lelpo_NCEP2,lnods_NCEP2,s_po_NCEP2,t_po_NCEP2, &
       nx,ny,nz,npoin_DAT,w,zlayer)
  !
  !*** p(nx,ny,nz) pressure
  !
  do iz = 1,nz_NCEP2
     work_NCEP21(:,:,iz) = p_NCEP2(iz)                                   ! P at time
  end do
  !
  call nc_interpola3d(nx_NCEP2,ny_NCEP2,nz_NCEP2,nelem_NCEP2, &
       work_NCEP21,z_NCEP2,lelpo_NCEP2,lnods_NCEP2,s_po_NCEP2,t_po_NCEP2, &
       nx,ny,nz,npoin_DAT,p,zlayer)
  !
  !*** Release memory
  !
  deallocate(work_NCEP21)
  deallocate(work_NCEP22)
  !
  !*** Other derived 3D variables and transformations (already in the dbs grid!)
  !
  allocate(work_NCEP21(nx,ny,nz))
  allocate(work_NCEP22(nx,ny,nz))
  !
  !*** qv(nx,ny,nz) Specific humidity (NCEP2 gives relative humidity, by now stored in qv)
  !*** Following Bolton (1980)
  !***      es = 6.112 * exp((17.67 * T)/(T + 243.5));  saturation vapor pressure in mb, T in C
  !***      e  = es * (RH/100.0);                       vapor pressure in mb
  !***      q  = (0.622 * e)/(p - (0.378 * e))          specific humidity in kg/kg, p in mb
  !
  work_NCEP21(:,:,:) = 6.112*exp( (17.67*(T(:,:,:)-273.15))/(243.5+(T(:,:,:)-273.15)) )  ! Saturation Pressure in mb
  work_NCEP21(:,:,:) = work_NCEP21(:,:,:)*qv(:,:,:)/100.0_rp                               ! Vapor      Pressure in mb
  qv(:,:,:) = 0.622*work_NCEP21(:,:,:)/((p(:,:,:)/101_rp)-0.378*work_NCEP21(:,:,:))        ! Specific   humidity in kg/kg
  !
  !*** Tv(nx,ny,nz) virtual temperature
  !
  work_NCEP21 = qv(:,:,:)/(1.0_rp-qv(:,:,:))      ! Mixing ratio (humidity ratio)
  !
  Tv(:,:,:) = T(:,:,:)*(1.0_rp + 1.6077_rp*work_NCEP21(:,:,:))/(1.0_rp+work_NCEP21(:,:,:))
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
  deallocate(work_NCEP21)
  deallocate(work_NCEP22)
  !
  !***  Closes
  !
  if( nf90_close(ncID) /= 0) call runend('read_NCEP2 : Error in nf90_close')
  !
  !***  Computes other variables not given by NCEP2
  !
  !
  !*** PBL(nx,ny) boundary layer height. Not given by ARPA. Value assumed
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
end subroutine read_NCEP2_data
!
!
!
subroutine read_NCEP2_translate(var,nx,ny,nz)
  !**********************************************************************************
  !*
  !*   Changes the origin
  !*
  !**********************************************************************************
  use KindType, ONLY :  ip,rp
  implicit none
  !
  integer(ip) :: nx,ny,nz
  real   (rp) :: var(nx,ny,nz)
  !
  integer(ip) :: ix,iy,iz,is
  real   (rp), allocatable :: work(:)
  !
  !*** Upside down
  !
  allocate(work(ny))
  is = ny/2 + mod(ny,2)   ! shift position
  !
  do iz = 1,nz
     do ix = 1,nx
        do iy = 1,ny
           work(iy) = var(ix,ny-iy+1,iz)
        end do
        var(ix,1:ny,iz) = work(1:ny)
     end do
  end do
  deallocate(work)
  !
  !*** Changes the origin
  !
  allocate(work(nx))
  is = nx/2 + mod(nx,2)   ! shift position
  !
  do iz = 1,nz
     do iy = 1,ny
        if(mod(nx,2).eq.1) work(is) = var(is,iy,iz)
        do ix = 1,nx/2
           work(ix   ) = var(is+ix,iy,iz)
           work(is+ix) = var(   ix,iy,iz)
        end do
        var(1:nx,iy,iz) = work(1:nx)
     end do
  end do
  deallocate(work)
  !
  return
end subroutine read_NCEP2_translate
