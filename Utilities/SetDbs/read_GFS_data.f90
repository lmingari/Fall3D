subroutine read_GFS_data(fname,time_lag,time,timesec,nx,ny,nz,zlayer,lon,lat,topg,LDU,  &
     u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
  !********************************************************************
  !*
  !*   Reads GFS data in NetCDF format. This routine MUST be called
  !*   after read_GFS_grid
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use netcdf
  use InpOut
  use MathFun
  use GFS_nc
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
  if(nt_GFS.gt.1) then
     found  = .false.
     it_GFS = 0
     do while(.not.found)
        it_GFS = it_GFS + 1
        if( (it_GFS+1).gt.nt_GFS ) call runend('read_GFS : Time not found')
        if( (time_lag+timesec.ge.timesec_GFS(it_GFS  )).and. &
             (time_lag+timesec.le.timesec_GFS(it_GFS+1)) ) then
           found = .true.
           itime1 = it_GFS
           itime2 = it_GFS + 1
           time_s = (timesec_GFS(itime2)-time_lag-timesec)/(timesec_GFS(itime2)-timesec_GFS(itime1))
        end if
     end do
  else
     itime1 = 1
     itime2 = 1
     time_s =1.0_rp
  end if
  !
  write(lulog,10) time,time_GFS(itime1),time_GFS(itime2),time_s
10 format('Interpolating at ',f15.0,/, &
       '        GFS from ',f15.0,' to ',f15.0,' s = ',f6.3)
  !
  !*** Open netCDF file and get ncID
  !
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) call runend('read_GFS : Error in nf90_open')
  !
  !*** First of all reads topograpgy form GFS. This is necessary since interpolation in done
  !*** in terrain following
  !
  allocate(work_GFS1(nx_GFS,ny_GFS,1))
  allocate(work_GFS2(nx_GFS,ny_GFS,1))
  !
  !*** LDU (no given)
  !
  LDU = -99
  !
  !*** HGT:surface(nx_GFS,ny_GFS) -> topg(nx,ny) (time independent actually)
  !*** NOTE: translation of origin (from 0 to lonmin) is necessary each time a
  !*** variable is read
  !
  if( nf90_inq_varid(ncID,fis_GFS_name,fisID) /= 0) call runend('read_GFS : Error getting fisID')
  if( nf90_get_var(ncID,fisID,work_GFS1,start=(/1,1,itime1/)) /=0 ) call runend('read_GFS: Error reading fis')
  call read_GFS_translate(work_GFS1,nx_GFS,ny_GFS,1)
  !
  !***  LAND(nx_GFS,ny_GFS) Land Mask (time independent actually)
  !***  0 water 1 land.
  !
  if( nf90_inq_varid(ncID,lmask_GFS_name,lmaskID) /= 0) call runend('read_GFS : Error getting lmaskID')
  if( nf90_get_var(ncID,lmaskID,work_GFS2,start=(/1,1,itime1/)) /=0 ) call runend('read_GFS: Error reading landmask')
  call read_GFS_translate(work_GFS2,nx_GFS,ny_GFS,1)
  !
  work_GFS1 = work_GFS1*work_GFS2
  !
  call nc_interpola2d(nx_GFS,ny_GFS,nelem_GFS, &
       work_GFS1,lelpo_GFS,lnods_GFS,s_po_GFS,t_po_GFS, &
       nx,ny,npoin_DAT,topg)
  topg_GFS(:,:) = work_GFS1(:,:,1)           ! needed later
  deallocate(work_GFS1)
  deallocate(work_GFS2)
  !
  !*** Reads 3D variables and interpolates
  !
  allocate (work_GFS1(nx_GFS,ny_GFS,nz_GFS))
  allocate (work_GFS2(nx_GFS,ny_GFS,nz_GFS))
  !
  !*** z_GFS(nx_GFS,ny_GFS,nz_GFS) Heights. This is necessary because pressure levels change with time.
  !*** GFS already gives the geopotential in m (i.e. height). The topography is substracted in order to interpolate
  !*** in terrain following
  !
  if( nf90_inq_varid(ncID,fi_GFS_name,fiID) /= 0) call runend('read_GFS : Error getting fiID')
  if( nf90_get_var(ncID,fiID,work_GFS1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_GFS: Error reading fi')
  if( nf90_get_var(ncID,fiID,work_GFS2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_GFS: Error reading fi')
  !
  call read_GFS_translate(work_GFS1,nx_GFS,ny_GFS,nz_GFS)
  call read_GFS_translate(work_GFS2,nx_GFS,ny_GFS,nz_GFS)
  z_GFS = time_s*work_GFS1 + (1.0_rp-time_s)*work_GFS2     ! z at time
  !
  do iz = 1,nz_GFS
     z_GFS(:,:,iz) =  z_GFS(:,:,iz) - topg_GFS(:,:)
  end do
  !
  if(MAXVAL(z_GFS(:,:,:)).lt.(zlayer(nz))) &
       call wriwar('Dbs max vertical layer is higher than GFS max layer. Values will be extrapolated')
  !
  !*** TMP(nx_GFS,ny_GFS,nz_GFS) --> T(nx,ny,nz)  Temperature
  !
  if( nf90_inq_varid(ncID,T_GFS_name,tID) /= 0) call runend('read_GFS : Error getting tID')
  if( nf90_get_var(ncID,tID,work_GFS1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_GFS: Error reading T')
  if( nf90_get_var(ncID,tID,work_GFS2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_GFS: Error reading T')
  !
  call read_GFS_translate(work_GFS1,nx_GFS,ny_GFS,nz_GFS)
  call read_GFS_translate(work_GFS2,nx_GFS,ny_GFS,nz_GFS)
  work_GFS1 = time_s*work_GFS1 + (1.0_rp-time_s)*work_GFS2   ! T at time
  !
  call nc_interpola3d(nx_GFS,ny_GFS,nz_GFS,nelem_GFS, &
       work_GFS1,z_GFS,lelpo_GFS,lnods_GFS,s_po_GFS,t_po_GFS, &
       nx,ny,nz,npoin_DAT,T,zlayer)
  !
  !*** RH(nx_GFS,ny_GFS,nz_GFS) Relative humidity --> qv(nx,ny,nz)  Specific humidity is computed later
  !
  if( nf90_inq_varid(ncID,Q_GFS_name,qID) /= 0) call runend('read_GFS : Error getting tID')
  if( nf90_get_var(ncID,qID,work_GFS1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_GFS: Error reading QV')
  if( nf90_get_var(ncID,qID,work_GFS2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_GFS: Error reading QV')
  !
  call read_GFS_translate(work_GFS1,nx_GFS,ny_GFS,nz_GFS)
  call read_GFS_translate(work_GFS2,nx_GFS,ny_GFS,nz_GFS)
  work_GFS1 = time_s*work_GFS1 + (1.0_rp-time_s)*work_GFS2   ! QV at time
  !
  call nc_interpola3d(nx_GFS,ny_GFS,nz_GFS,nelem_GFS, &
       work_GFS1,z_GFS,lelpo_GFS,lnods_GFS,s_po_GFS,t_po_GFS, &
       nx,ny,nz,npoin_DAT,qv,zlayer)
  !
  !*** U(nx_GFS,ny_GFS,nz_GFS) --> u(nx,ny,nz)  Velocity
  !
  if( nf90_inq_varid(ncID,u_GFS_name,uID) /= 0) call runend('read_GFS : Error getting uID')
  if( nf90_get_var(ncID,uID,work_GFS1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_GFS: Error reading u')
  if( nf90_get_var(ncID,uID,work_GFS2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_GFS: Error reading u')
  !
  call read_GFS_translate(work_GFS1,nx_GFS,ny_GFS,nz_GFS)
  call read_GFS_translate(work_GFS2,nx_GFS,ny_GFS,nz_GFS)
  work_GFS1 = time_s*work_GFS1 + (1.0_rp-time_s)*work_GFS2   ! U at time
  !
  call nc_interpola3d(nx_GFS,ny_GFS,nz_GFS,nelem_GFS, &
       work_GFS1,z_GFS,lelpo_GFS,lnods_GFS,s_po_GFS,t_po_GFS, &
       nx,ny,nz,npoin_DAT,u,zlayer)
  !
  !*** V(nx_GFS,ny_GFS,nz_GFS) --> v(nx,ny,nz)  Velocity
  !
  if( nf90_inq_varid(ncID,v_GFS_name,vID) /= 0) call runend('read_GFS : Error getting vID')
  if( nf90_get_var(ncID,vID,work_GFS1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_GFS: Error reading v')
  if( nf90_get_var(ncID,vID,work_GFS2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_GFS: Error reading v')
  !
  call read_GFS_translate(work_GFS1,nx_GFS,ny_GFS,nz_GFS)
  call read_GFS_translate(work_GFS2,nx_GFS,ny_GFS,nz_GFS)
  work_GFS1 = time_s*work_GFS1 + (1.0_rp-time_s)*work_GFS2   ! V at time
  !
  call nc_interpola3d(nx_GFS,ny_GFS,nz_GFS,nelem_GFS, &
       work_GFS1,z_GFS,lelpo_GFS,lnods_GFS,s_po_GFS,t_po_GFS, &
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
  !*** W(nx_GFS,ny_GFS,nz_GFS) --> w(nx,ny,nz) Omega velocity
  !
  if( nf90_inq_varid(ncID,w_GFS_name,wID) /= 0) call runend('read_GFS : Error getting wID')
  if( nf90_get_var(ncID,wID,work_GFS1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_GFS: Error reading w')
  if( nf90_get_var(ncID,wID,work_GFS2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_GFS: Error reading w')
  !
  call read_GFS_translate(work_GFS1,nx_GFS,ny_GFS,nz_GFS)
  call read_GFS_translate(work_GFS2,nx_GFS,ny_GFS,nz_GFS)
  work_GFS1 = time_s*work_GFS1 + (1.0_rp-time_s)*work_GFS2            ! Omega at time
  do iz = 1,nz_GFS
     work_GFS1(:,:,iz) = work_GFS1(:,:,iz)/p_GFS(iz)                ! Omega/P at time
  end do
  !
  call nc_interpola3d(nx_GFS,ny_GFS,nz_GFS,nelem_GFS, &
       work_GFS1,z_GFS,lelpo_GFS,lnods_GFS,s_po_GFS,t_po_GFS, &
       nx,ny,nz,npoin_DAT,w,zlayer)
  !
  !*** p(nx,ny,nz) pressure
  !
  do iz = 1,nz_GFS
     work_GFS1(:,:,iz) = p_GFS(iz)                                   ! P at time
  end do
  !
  call nc_interpola3d(nx_GFS,ny_GFS,nz_GFS,nelem_GFS, &
       work_GFS1,z_GFS,lelpo_GFS,lnods_GFS,s_po_GFS,t_po_GFS, &
       nx,ny,nz,npoin_DAT,p,zlayer)
  !
  !*** Release memory
  !
  deallocate(work_GFS1)
  deallocate(work_GFS2)
  !
  !*** Rest of 2D scalar variables
  !
  allocate(work_GFS1(nx_GFS,ny_GFS,1))
  allocate(work_GFS2(nx_GFS,ny_GFS,1))
  !
  !*** HPBL(nx_GFS,ny_GFS) --> PBLH(nx,ny) Planetary boundary layer height
  !
  if( nf90_inq_varid(ncID,pblh_GFS_name,pblhID) /= 0) call runend('read_GFS : Error getting pblhID')
  if( nf90_get_var(ncID,pblhID,work_GFS1,start=(/1,1,itime1/)) /=0 ) call runend('read_GFS: Error reading pblh')
  if( nf90_get_var(ncID,pblhID,work_GFS2,start=(/1,1,itime2/)) /=0 ) call runend('read_GFS: Error reading pblh')
  !
  call read_GFS_translate(work_GFS1,nx_GFS,ny_GFS,1)
  call read_GFS_translate(work_GFS2,nx_GFS,ny_GFS,1)
  work_GFS1 = time_s*work_GFS1 + (1.0_rp-time_s)*work_GFS2   ! pblh box
  !
  call nc_interpola2d(nx_GFS,ny_GFS,nelem_GFS, &
       work_GFS1,lelpo_GFS,lnods_GFS,s_po_GFS,t_po_GFS, &
       nx,ny,npoin_DAT,pblh)
  !
  deallocate(work_GFS1)
  deallocate(work_GFS2)
  !
  !*** Other derived 3D variables and transformations (already in the dbs grid!)
  !
  allocate(work_GFS1(nx,ny,nz))
  allocate(work_GFS2(nx,ny,nz))
  !
  !*** qv(nx,ny,nz) Specific humidity (GFS gives relative humidity, by now stored in qv)
  !*** Following Bolton (1980)
  !***      es = 6.112 * exp((17.67 * T)/(T + 243.5));  saturation vapor pressure in mb, T in C
  !***      e  = es * (RH/100.0);                       vapor pressure in mb
  !***      q  = (0.622 * e)/(p - (0.378 * e))          specific humidity in kg/kg, p in mb
  !
  work_GFS1(:,:,:) = 6.112*exp( (17.67*(T(:,:,:)-273.15))/(243.5+(T(:,:,:)-273.15)) )  ! Saturation Pressure in mb
  work_GFS1(:,:,:) = work_GFS1(:,:,:)*qv(:,:,:)/100.0_rp                               ! Vapor      Pressure in mb
  qv(:,:,:) = 0.622*work_GFS1(:,:,:)/((p(:,:,:)/101_rp)-0.378*work_GFS1(:,:,:))        ! Specific   humidity in kg/kg
  !
  !*** Tv(nx,ny,nz) virtual temperature
  !
  work_GFS1 = qv(:,:,:)/(1.0_rp-qv(:,:,:))      ! Mixing ratio (humidity ratio)
  !
  Tv(:,:,:) = T(:,:,:)*(1.0_rp + 1.6077_rp*work_GFS1(:,:,:))/(1.0_rp+work_GFS1(:,:,:))
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
  deallocate(work_GFS1)
  deallocate(work_GFS2)
  !
  !***  Closes
  !
  if( nf90_close(ncID) /= 0) call runend('read_GFS : Error in nf90_close')
  !
  !***  Computes other variables not given by GFS
  !
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
end subroutine read_GFS_data
!
!
!
subroutine read_GFS_translate(var,nx,ny,nz)
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
  !*** allocates memory
  !
  allocate(work(nx))
  !
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
  !
  deallocate(work)
  !
  return
end subroutine read_GFS_translate
