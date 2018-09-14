subroutine read_ARPA_data(fname,time_lag,time,timesec,nx,ny,nz,zlayer, &
     lon,lat,topg,LDU,u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
  !********************************************************************
  !*
  !*   Reads ARPA data in NetCDF format. This routine MUST be called
  !*   after read_ARPA_grid
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use netcdf
  use InpOut
  use MathFun
  use ARPA_nc
  use wind_rotation
  implicit none
  !
  character(len=*) :: fname
  integer(ip)      :: nx,ny,nz
  real(rp)         :: time_lag,time
  real(rp)         :: timesec                     ! time in sec after 0000UTC
  real(rp)         :: zlayer(nz),lat(nx,ny),lon(nx,ny),topg(nx,ny),LDU(nx,ny)
  real(rp)         ::  u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz),T(nx,ny,nz), &
       Tv(nx,ny,nz),Tp(nx,ny,nz),qv(nx,ny,nz),ro(nx,ny,nz),p(nx,ny,nz)
  real(rp)         :: pblh(nx,ny),ust(nx,ny),L(nx,ny)
  !
  logical          :: found
  integer          :: itime1,itime2,ix,iy,iz
  real(rp)         :: time_s
  integer :: status
  integer :: i
  !
  !*** Finds the time itime1 and itime2 and the interpolation factor time_s
  !
  if(nt_ARPA.gt.1) then
     found  = .false.
     it_ARPA = 0
     do while(.not.found)
        it_ARPA = it_ARPA + 1
        if( (it_ARPA+1).gt.nt_ARPA ) call runend('read_ARPA : Time not found')
        if( (time_lag+timesec.ge.timesec_ARPA(it_ARPA  )).and. &
             (time_lag+timesec.le.timesec_ARPA(it_ARPA+1)) ) then
           found = .true.
           itime1 = it_ARPA
           itime2 = it_ARPA + 1
           time_s = (timesec_ARPA(itime2)-time_lag-timesec)/(timesec_ARPA(itime2)-timesec_ARPA(itime1))
        end if
     end do
  else
     itime1 = 1
     itime2 = 1
     time_s =1.0_rp
  end if
  !
  write(lulog,10) time,time_ARPA(itime1),time_ARPA(itime2),time_s
10 format('Interpolating at ',f15.0,/, &
       '        ARPA from ',f15.0,' to ',f15.0,' s = ',f6.3)
  !
  !*** Open netCDF file and get ncID
  !
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) call runend('read_ARPA : Error in nf90_open')
  !
  !*** First of all reads topograpgy form ARPA.
  !    This is necessary since interpolation in done in terrain following
  !
  if(do_transpose) then
     allocate(work_ARPA1(ny_ARPA,nx_ARPA,1))
  else
     allocate(work_ARPA1(nx_ARPA,ny_ARPA,1))
  end if
  !
  !*** FI:sfc(nx_ARPA,ny_ARPA) Geopotential FI at surface -> topg(nx,ny) (time independent actually)
  !
  status = nf90_inq_varid(ncID,fis_ARPA_name,fisID)
  if(status == 0) then
     if( nf90_get_var(ncID,fisID,work_ARPA1,start=(/1,1,itime1/)) /=0 ) call runend('read_ARPA: Error reading fis')
     !
     work_ARPA1 = work_ARPA1 / 9.81_rp            ! topg = (fis)/g
     if(do_transpose) work_ARPA1(:,:,1) = transpose(work_ARPA1(:,:,1))
     call nc_interpola2d(nx_ARPA,ny_ARPA,nelem_ARPA, &
          work_ARPA1,lelpo_ARPA,lnods_ARPA,s_po_ARPA,t_po_ARPA, &
          nx,ny,npoin_DAT,topg)
     topg_ARPA(:,:) = work_ARPA1(:,:,1)           ! needed later
  else
     ! call runend('read_ARPA : Error getting fisID')
     call wriwar('Cannot find variable FI:sfc in netCDF file: setting TOPOG=0')
     topg_ARPA(:,:) = 0.0_rp  ! WARNING: SET=0
  end if
  deallocate(work_ARPA1)

  !
  !*** LDU (no given)
  !
  LDU = -99
  !
  !
  !*** Reads 3D variables and interpolates
  !
  !
  !*** z_ARPA(nx_ARPA,ny_ARPA,nz_ARPA) Heights.
  !    This is necessary because pressure levels change with time.
  !*** Heights are computed from the geopotential FI.
  !    The topography is subtracted in order to interpolate in terrain following
  !
  if(do_transpose) then
     allocate (work_ARPA1(ny_ARPA,nx_ARPA,nz_ARPA))
     allocate (work_ARPA2(ny_ARPA,nx_ARPA,nz_ARPA))
  else
     allocate (work_ARPA1(nx_ARPA,ny_ARPA,nz_ARPA))
     allocate (work_ARPA2(nx_ARPA,ny_ARPA,nz_ARPA))
  end if
  !
  if( nf90_inq_varid(ncID,fi_ARPA_name,fiID) /= 0) call runend('read_ARPA : Error getting fiID')
  if( nf90_get_var(ncID,fiID,work_ARPA1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ARPA: Error reading fi')
  work_ARPA1 =  work_ARPA1 / 9.81_rp      ! z = (fi)/g at time1
  if( nf90_get_var(ncID,fiID,work_ARPA2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ARPA: Error reading fi')
  work_ARPA2 = work_ARPA2 / 9.81_rp      ! z = (fi)/g at time2
  !
  if(do_transpose) then
     do i=1,nz_ARPA
        work_ARPA1(:,:,i) = transpose(work_ARPA1(:,:,i))
        work_ARPA2(:,:,i) = transpose(work_ARPA2(:,:,i))
     end do
  end if
  !
  z_ARPA = time_s*work_ARPA1 + (1.0_rp-time_s)*work_ARPA2     ! z at time
  do iz = 1,nz_ARPA
     z_ARPA(:,:,iz) =  z_ARPA(:,:,iz) - topg_ARPA(:,:)
  end do
  !
  if(MAXVAL(z_ARPA(:,:,:)).lt.(zlayer(nz))) &
       call wriwar('Dbs max vertical layer is higher than ARPA max layer. Values will be extrapolated')
  deallocate(work_ARPA1)
  deallocate(work_ARPA2)

  !
  !*** T(nx_ARPA,ny_ARPA,nz_ARPA) --> T(nx,ny,nz)  Temperature
  !
  if(do_transpose) then
     allocate (work_ARPA1(ny_ARPA,nx_ARPA,nz_ARPA))
     allocate (work_ARPA2(ny_ARPA,nx_ARPA,nz_ARPA))
  else
     allocate (work_ARPA1(nx_ARPA,ny_ARPA,nz_ARPA))
     allocate (work_ARPA2(nx_ARPA,ny_ARPA,nz_ARPA))
  end if
  if( nf90_inq_varid(ncID,T_ARPA_name,tID) /= 0) call runend('read_ARPA : Error getting tID')
  if( nf90_get_var(ncID,tID,work_ARPA1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ARPA: Error reading T')
  if( nf90_get_var(ncID,tID,work_ARPA2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ARPA: Error reading T')
  work_ARPA1 = time_s*work_ARPA1 + (1.0_rp-time_s)*work_ARPA2   ! T at time
  !
  if(do_transpose) then
     do i=1,nz_ARPA
        work_ARPA1(:,:,i) = transpose(work_ARPA1(:,:,i))
        work_ARPA2(:,:,i) = transpose(work_ARPA2(:,:,i))
     end do
  end if
  !
  call nc_interpola3d(nx_ARPA,ny_ARPA,nz_ARPA,nelem_ARPA, &
       work_ARPA1,z_ARPA,lelpo_ARPA,lnods_ARPA,s_po_ARPA,t_po_ARPA, &
       nx,ny,nz,npoin_DAT,T,zlayer)
  !
  deallocate(work_ARPA1)
  deallocate(work_ARPA2)

  !
  !*** Qv(nx_ARPA,ny_ARPA,nz_ARPA) --> qv(nx,ny,nz)  Specific humidity
  !
  if(do_transpose) then
     allocate (work_ARPA1(ny_ARPA,nx_ARPA,nz_ARPA))
     allocate (work_ARPA2(ny_ARPA,nx_ARPA,nz_ARPA))
  else
     allocate (work_ARPA1(nx_ARPA,ny_ARPA,nz_ARPA))
     allocate (work_ARPA2(nx_ARPA,ny_ARPA,nz_ARPA))
  end if
  !
  if( nf90_inq_varid(ncID,Q_ARPA_name,qID) /= 0) call runend('read_ARPA : Error getting tID')
  if( nf90_get_var(ncID,qID,work_ARPA1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ARPA: Error reading QV')
  if( nf90_get_var(ncID,qID,work_ARPA2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ARPA: Error reading QV')
  work_ARPA1 = time_s*work_ARPA1 + (1.0_rp-time_s)*work_ARPA2   ! QV at time
  !
  if(do_transpose) then
     do i=1,nz_ARPA
        work_ARPA1(:,:,i) = transpose(work_ARPA1(:,:,i))
        work_ARPA2(:,:,i) = transpose(work_ARPA2(:,:,i))
     end do
  end if
  !
  call nc_interpola3d(nx_ARPA,ny_ARPA,nz_ARPA,nelem_ARPA, &
       work_ARPA1,z_ARPA,lelpo_ARPA,lnods_ARPA,s_po_ARPA,t_po_ARPA, &
       nx,ny,nz,npoin_DAT,qv,zlayer)
  !
  deallocate(work_ARPA1)
  deallocate(work_ARPA2)
  !
  !*** U(nx_ARPA,ny_ARPA,nz_ARPA) --> u(nx,ny,nz)  Rotated velocity
  !
  if(do_transpose) then
     allocate (work_ARPA1(ny_ARPA,nx_ARPA,nz_ARPA))
     allocate (work_ARPA2(ny_ARPA,nx_ARPA,nz_ARPA))
  else
     allocate (work_ARPA1(nx_ARPA,ny_ARPA,nz_ARPA))
     allocate (work_ARPA2(nx_ARPA,ny_ARPA,nz_ARPA))
  end if
  !
  if( nf90_inq_varid(ncID,u_ARPA_name,uID) /= 0) call runend('read_ARPA : Error getting uID')
  if( nf90_get_var(ncID,uID,work_ARPA1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ARPA: Error reading u')
  if( nf90_get_var(ncID,uID,work_ARPA2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ARPA: Error reading u')
  work_ARPA1 = time_s*work_ARPA1 + (1.0_rp-time_s)*work_ARPA2   ! U at time
  !
  if(do_transpose) then
     do i=1,nz_ARPA
        work_ARPA1(:,:,i) = transpose(work_ARPA1(:,:,i))
        work_ARPA2(:,:,i) = transpose(work_ARPA2(:,:,i))
     end do
  end if
  !
  call nc_interpola3d(nx_ARPA,ny_ARPA,nz_ARPA,nelem_ARPA, &
       work_ARPA1,z_ARPA,lelpo_ARPA,lnods_ARPA,s_po_ARPA,t_po_ARPA, &
       nx,ny,nz,npoin_DAT,u,zlayer)
  !
  deallocate(work_ARPA1)
  deallocate(work_ARPA2)
  !
  !*** V(nx_ARPA,ny_ARPA,nz_ARPA) --> v(nx,ny,nz)  Rotated velocity
  !
  if(do_transpose) then
     allocate (work_ARPA1(ny_ARPA,nx_ARPA,nz_ARPA))
     allocate (work_ARPA2(ny_ARPA,nx_ARPA,nz_ARPA))
  else
     allocate (work_ARPA1(nx_ARPA,ny_ARPA,nz_ARPA))
     allocate (work_ARPA2(nx_ARPA,ny_ARPA,nz_ARPA))
  end if
  !
  if( nf90_inq_varid(ncID,v_ARPA_name,vID) /= 0) call runend('read_ARPA : Error getting vID')
  if( nf90_get_var(ncID,vID,work_ARPA1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ARPA: Error reading v')
  if( nf90_get_var(ncID,vID,work_ARPA2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ARPA: Error reading v')
  work_ARPA1 = time_s*work_ARPA1 + (1.0_rp-time_s)*work_ARPA2   ! V at time
  !
  if(do_transpose) then
     do i=1,nz_ARPA
        work_ARPA1(:,:,i) = transpose(work_ARPA1(:,:,i))
        work_ARPA2(:,:,i) = transpose(work_ARPA2(:,:,i))
     end do
  end if
  !
  call nc_interpola3d(nx_ARPA,ny_ARPA,nz_ARPA,nelem_ARPA, &
       work_ARPA1,z_ARPA,lelpo_ARPA,lnods_ARPA,s_po_ARPA,t_po_ARPA, &
       nx,ny,nz,npoin_DAT,v,zlayer)
  !
  deallocate(work_ARPA1)
  deallocate(work_ARPA2)
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
  !*** W(nx_ARPA,ny_ARPA,nz_ARPA) --> w(nx,ny,nz) Omega velocity
  !
  if(do_transpose) then
     allocate (work_ARPA1(ny_ARPA,nx_ARPA,nz_ARPA))
     allocate (work_ARPA2(ny_ARPA,nx_ARPA,nz_ARPA))
  else
     allocate (work_ARPA1(nx_ARPA,ny_ARPA,nz_ARPA))
     allocate (work_ARPA2(nx_ARPA,ny_ARPA,nz_ARPA))
  end if
  !
  if( nf90_inq_varid(ncID,w_ARPA_name,wID) /= 0) call runend('read_ARPA : Error getting wID')
  if( nf90_get_var(ncID,wID,work_ARPA1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_ARPA: Error reading w')
  if( nf90_get_var(ncID,wID,work_ARPA2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_ARPA: Error reading w')
  work_ARPA1 = time_s*work_ARPA1 + (1.0_rp-time_s)*work_ARPA2            ! Omega at time
  do iz = 1,nz_ARPA
     work_ARPA1(:,:,iz) = work_ARPA1(:,:,iz)/p_ARPA(iz)                ! Omega/P at time
  end do
  !
  if(do_transpose) then
     do i=1,nz_ARPA
        work_ARPA1(:,:,i) = transpose(work_ARPA1(:,:,i))
        work_ARPA2(:,:,i) = transpose(work_ARPA2(:,:,i))
     end do
  end if
  !
  call nc_interpola3d(nx_ARPA,ny_ARPA,nz_ARPA,nelem_ARPA, &
       work_ARPA1,z_ARPA,lelpo_ARPA,lnods_ARPA,s_po_ARPA,t_po_ARPA, &
       nx,ny,nz,npoin_DAT,w,zlayer)
  !
  deallocate(work_ARPA1)
  deallocate(work_ARPA2)
  !
  !*** p(nx,ny,nz) pressure
  !
  allocate (work_ARPA1(nx_ARPA,ny_ARPA,nz_ARPA))
  do iz = 1,nz_ARPA
     work_ARPA1(:,:,iz) = p_ARPA(iz)                  ! P at time
  end do
  !
  call nc_interpola3d(nx_ARPA,ny_ARPA,nz_ARPA,nelem_ARPA, &
       work_ARPA1,z_ARPA,lelpo_ARPA,lnods_ARPA,s_po_ARPA,t_po_ARPA, &
       nx,ny,nz,npoin_DAT,p,zlayer)
  deallocate(work_ARPA1)
  !
  !*** Release memory
  !
  !
  !*** Other derived 3D variables and transformations
  !
  allocate(work_ARPA1(nx,ny,nz))
  allocate(work_ARPA2(nx,ny,nz))
  !
  !*** Tv(nx,ny,nz) virtual temperature
  !
  work_ARPA1 = qv(:,:,:)/(1.0_rp-qv(:,:,:))      ! Mixing ratio (humidity ratio)
  !
  Tv(:,:,:) = T(:,:,:)*(1.0_rp + 1.6077_rp*work_ARPA1(:,:,:))/(1.0_rp+work_ARPA1(:,:,:))
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
  !*** Rotatate velocity components
  !
  work_ARPA1(:,:,:) = u(:,:,:)
  work_ARPA2(:,:,:) = v(:,:,:)
  !
  call rotate_ARPA_wind(work_ARPA1,work_ARPA2,nx,ny,nz,      &
       ARPA_cen_lon,ARPA_CEN_lat,lat,lon,u,v)
  !
  deallocate(work_ARPA1)
  deallocate(work_ARPA2)
  !
  !***  Closes
  !
  if( nf90_close(ncID) /= 0) call runend('read_ARPA : Error in nf90_close')
  !
  !***  Computes other variables not given by ARPA
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
end subroutine read_ARPA_data
!
!
!
subroutine rotate_ARPA_wind(u,v,d1,d2,d3,cen_lon,cen_lat,xlat,xlon,urot,vrot)
  !***************************************************************
  !*
  !*  Rotates ARPA wind
  !*
  !***************************************************************
  use KindType, ONLY :  ip,rp
  implicit none
  integer(ip), intent(in) ::  d1,d2,d3
  real(rp)  :: u(d1,d2,d3),v(d1,d2,d3)
  real(rp)  :: xlat(d1,d2),xlon(d1,d2)
  real(rp)  :: urot(d1,d2,d3),vrot(d1,d2,d3)
  real(rp)  :: cen_lon, cen_lat
  !
  integer(ip) :: i,j,k
  !
  real(rp), parameter :: pii = 3.14159265_rp
  real(rp), parameter :: radians_per_degree = pii/180_rp
  real(rp) :: stp0,ctp0,slon,clon,slat,clat,tph,cray,dray,dc
  !
  stp0 = sin( cen_lat * radians_per_degree)
  ctp0 = cos( cen_lat * radians_per_degree)
  !
  do j = 1,d2
     do i = 1,d1
        !
        slon = sin( (xlon(i,j)-cen_lon) * radians_per_degree)
        clon = cos( (xlon(i,j)-cen_lon) * radians_per_degree)
        slat = sin(  xlat(i,j)          * radians_per_degree)
        clat = cos(  xlat(i,j)          * radians_per_degree)
        tph  = asin(ctp0*slat - stp0*clat*clon)
        tph  = 1.0_rp/cos(tph)
        cray = (            stp0*slon     )*tph
        dray = (ctp0*clat + stp0*slat*clon)*tph
        dc   = 1.0_rp/(dray*dray+cray*cray)
        !
        do k =1,d3
           urot(i,j,k) = (dray*u(i,j,k)+cray*v(i,j,k))*dc
           vrot(i,j,k) = (dray*v(i,j,k)-cray*u(i,j,k))*dc
        end do
     end do
  end do
  !
  return
end subroutine rotate_ARPA_wind
