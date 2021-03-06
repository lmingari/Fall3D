subroutine read_WRF_data(fname,time_lag,time,timesec,nx,ny,nz,zlayer,lon,lat,topg,LDU, &
     SOIL,vegfra,u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,rmol,znt,spd10,L,smoi,prat)
  !********************************************************************
  !*
  !*   Reads WRF data in NetCDF format. This routine MUST be called
  !*   after read_WRF_grid
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use InpOut
  use MathFun
  use netcdf
  use WRF_nc
  use wind_rotation
  implicit none
  !
  character(len=*) :: fname
  integer(ip)      :: nx,ny,nz
  real(rp)         :: time_lag,time
  real(rp)         :: timesec                     ! time in sec after 0000UTC
  real(rp)         :: zlayer(nz),lat(nx,ny),lon(nx,ny),topg(nx,ny),LDU(nx,ny),SOIL(nx,ny),vegfra(nx,ny)
  real(rp)         :: u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz),T(nx,ny,nz),Tv(nx,ny,nz)
  real(rp)         :: Tp(nx,ny,nz),qv(nx,ny,nz),ro(nx,ny,nz),p(nx,ny,nz)
  real(rp)         :: pblh(nx,ny),ust(nx,ny),rmol(nx,ny),znt(nx,ny)
  real(rp)         :: spd10(nx,ny),L(nx,ny),smoi(nx,ny),prat(nx,ny)
  !
  logical          :: found
  integer          :: itime1,itime2,ix,iy,iz
  real(rp)         :: time_s
  !
  !*** Finds the time itime1 and itime2 and the interpolation factor time_s
  !
  if(nt_WRF.gt.1) then
     found  = .false.
     it_WRF = 0
     do while(.not.found)
        it_WRF = it_WRF + 1
        if( (it_WRF+1).gt.nt_WRF ) call runend('read_WRF : Time not found')
        if( (time_lag+timesec.ge.timesec_WRF(it_WRF  )).and. &
             (time_lag+timesec.le.timesec_WRF(it_WRF+1)) ) then
           found = .true.
           itime1 = it_WRF
           itime2 = it_WRF + 1
           time_s = (timesec_WRF(itime2)-time_lag-timesec)/(timesec_WRF(itime2)-timesec_WRF(itime1))
        end if
     end do
  else
     itime1 = 1
     itime2 = 1
     time_s = 1.0_rp
  end if
  !
  write(lulog,10) time,time_WRF(itime1),time_WRF(itime2),time_s
10 format('Interpolating at ',f15.0,/, &
       '        WRF from ',f15.0,' to ',f15.0,' s = ',f6.3)
  !
  !*** Open netCDF file and get ncID
  !
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) call runend('read_WRF : Error in nf90_open')
  !
  !*** First of all reads topograpgy form WRF. This is necessary since interpolation in done
  !*** in terrain following
  !
  allocate(work_WRF (nx_WRF,ny_WRF,1))
  allocate(work_WRF1(nx_WRF,ny_WRF,1))
  !
  !***  LDU(nx_WRF,ny_WRF) USGS Land Use Category (time independent actually)
  !
  if( nf90_inq_varid(ncID,ldu_WRF_name,lduID) /= 0) call runend('read_WRF : Error getting lduID')
  if( nf90_get_var(ncID,lduID,work_WRF,start=(/1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading LDU')
  !
  call nc_assign2d(nx_WRF,ny_WRF,nelem_WRF, &
       work_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,npoin_DAT,LDU)
  !
  !***  SOIL(nx_WRF,ny_WRF) Soil Category (time independent actually)
  !
  if( nf90_inq_varid(ncID,soil_WRF_name,soilID) /= 0) call runend('read_WRF : Error getting soilID')
  if( nf90_get_var(ncID,soilID,work_WRF,start=(/1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading SOIL')
  !
  call nc_assign2d(nx_WRF,ny_WRF,nelem_WRF, &
       work_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,npoin_DAT,SOIL)
  !
  !***  VEGFRA(nx_WRF,ny_WRF) Vegetation Fraction (time independent actually)
  !
  if( nf90_inq_varid(ncID,vegfra_WRF_name,vegfraID) /= 0) call runend('read_WRF : Error getting vegfraID')
  if( nf90_get_var(ncID,vegfraID,work_WRF,start=(/1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading VEGFRA')
  !
  call nc_interpola2d(nx_WRF,ny_WRF,nelem_WRF, &
       work_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,npoin_DAT,vegfra)
  !
  !*** HGT(nx_WRF,ny_WRF) --> topg(nx,ny) Topograpgy (time independent actually)
  !
  if( nf90_inq_varid(ncID,hgt_WRF_name,hgtID) /= 0) call runend('read_WRF : Error getting hgtID')
  if( nf90_get_var(ncID,hgtID,work_WRF,start=(/1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading hgt')
  !
  !***  LANDMASK(nx_WRF,ny_WRF) Land Mask (time independent actually)
  !***  0 water 1 land.
  !
  if( nf90_inq_varid(ncID,lmask_WRF_name,lmaskID) /= 0) call runend('read_WRF : Error getting lmaskID')
  if( nf90_get_var(ncID,lmaskID,work_WRF1,start=(/1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading landmask')
  !
  work_WRF = work_WRF*work_WRF1
  call nc_interpola2d(nx_WRF,ny_WRF,nelem_WRF, &
       work_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,npoin_DAT,topg)
  topg_WRF(:,:) = work_WRF(:,:,1)           ! needed later
  !
  !*** PRATEC(nx_WRF,ny_WRF) Precipitation rate mm s-1
  !
  if( nf90_inq_varid(ncID,prat_WRF_name,pratID) /= 0) then
     prat = 0.0
     !         call runend('read_WRF : Error getting pratID')
  else
     if( nf90_get_var(ncID,pratID,work_WRF,start=(/1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading pratec')
     !
     call nc_interpola2d(nx_WRF,ny_WRF,nelem_WRF, &
          work_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
          nx,ny,npoin_DAT,prat)
  end if
  !
  deallocate(work_WRF )
  deallocate(work_WRF1)
  !
  !*** Reads 3D variables and interpolates
  !
  !
  !*** z_WRF(nx_WRF,ny_WRF,nz_WRF) Heights. This is necessary because pressure levels change with time.
  !*** Heights are already in terrain following.
  !*** Reads geopotentials PHB, PH and computes z and zstag at time
  !
  allocate (work_WRF1(nx_WRF,ny_WRF,nzstag_WRF))
  allocate (work_WRF2(nx_WRF,ny_WRF,nzstag_WRF))
  !
  if( nf90_inq_varid(ncID,phb_WRF_name,phbID) /= 0) call runend('read_WRF : Error getting phbID')
  if( nf90_get_var(ncID,phbID,work_WRF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading phb')
  phb_WRF = work_WRF1  ! time independent
  !
  if( nf90_inq_varid(ncID,ph_WRF_name,phID) /= 0) call runend('read_WRF : Error getting phbID')
  if( nf90_get_var(ncID,phID,work_WRF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading ph')
  work_WRF1 = ( work_WRF1 + phb_WRF ) / 9.81_rp      ! zstag = (ph+phb)/g at time1
  if( nf90_get_var(ncID,phID,work_WRF2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading ph')
  work_WRF2 = ( work_WRF2 + phb_WRF ) / 9.81_rp      ! zstag = (ph+phb)/g at time2
  !
  zstag_WRF = time_s*work_WRF1 + (1.0_rp-time_s)*work_WRF2                          ! zstag at time
  z_WRF = 0.5_rp*( zstag_WRF(:,:,1:nzstag_WRF-1) + zstag_WRF(:,:,2:nzstag_WRF) )  ! z at time
  !
  do iz = 1,nz_WRF
     z_WRF(:,:,iz) =  z_WRF(:,:,iz) - topg_WRF(:,:)
  end do
  !
  if(MAXVAL(z_WRF(:,:,:)).lt.(zlayer(nz))) &
       call wriwar('Dbs max vertical layer is higher than WRF max layer. Values extrapolated')
  !
  deallocate(work_WRF1)
  deallocate(work_WRF2)
  !
  !*** Reads 3D variables and interpolates
  !
  !
  !*** U(nxstag_WRF,ny_WRF,nz_WRF) --> u(nx,ny,nz)
  !
  allocate(work_WRF (nx_WRF    ,ny_WRF,nz_WRF))
  allocate(work_WRF1(nxstag_WRF,ny_WRF,nz_WRF))
  allocate(work_WRF2(nxstag_WRF,ny_WRF,nz_WRF))
  !
  if( nf90_inq_varid(ncID,u_WRF_name,uID) /= 0) call runend('read_WRF : Error getting uID')
  if( nf90_get_var(ncID,uID,work_WRF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading u')
  if( nf90_get_var(ncID,uID,work_WRF2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading u')
  work_WRF1 = time_s*work_WRF1 + (1.0_rp-time_s)*work_WRF2   ! U arakawa at time
  !
  work_WRF = 0.5_rp*(work_WRF1(1:nxstag_WRF-1,1:ny_WRF,1:nz_WRF)+ work_WRF1(2:nxstag_WRF,1:ny_WRF,1:nz_WRF)) ! u box
  !
  call nc_interpola3d(nx_WRF,ny_WRF,nz_WRF,nelem_WRF, &
       work_WRF,z_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,nz,npoin_DAT,u,zlayer)
  !
  deallocate(work_WRF )
  deallocate(work_WRF1)
  deallocate(work_WRF2)
  !
  !*** V(nx_WRF,nystag_WRF,nz_WRF) --> v(nx,ny,nz)
  !
  allocate(work_WRF (nx_WRF,ny_WRF    ,nz_WRF))
  allocate(work_WRF1(nx_WRF,nystag_WRF,nz_WRF))
  allocate(work_WRF2(nx_WRF,nystag_WRF,nz_WRF))
  !
  if( nf90_inq_varid(ncID,v_WRF_name,vID) /= 0) call runend('read_WRF : Error getting vID')
  if( nf90_get_var(ncID,vID,work_WRF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading v')
  if( nf90_get_var(ncID,vID,work_WRF2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading v')
  work_WRF1 = time_s*work_WRF1 + (1.0_rp-time_s)*work_WRF2   ! V arakawa at time
  !
  work_WRF = 0.5_rp*(work_WRF1(1:nx_WRF,1:nystag_WRF-1,1:nz_WRF)+ work_WRF1(1:nx_WRF,2:nystag_WRF,1:nz_WRF)) ! v box
  !
  call nc_interpola3d(nx_WRF,ny_WRF,nz_WRF,nelem_WRF, &
       work_WRF,z_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,nz,npoin_DAT,v,zlayer)
  !
  deallocate(work_WRF )
  deallocate(work_WRF1)
  deallocate(work_WRF2)
  !
  !*** Rotate velocity components (WRF operation)
  !
  allocate(work_WRF1(nx,ny,nz))
  allocate(work_WRF2(nx,ny,nz))
  work_WRF1(:,:,:) = u(:,:,:)
  work_WRF2(:,:,:) = v(:,:,:)
  !
  call rotate_WRF_wind(work_WRF1,work_WRF2,nx,ny,nz,      &
       WRF_map_proj,WRF_cen_lon ,lat,lon, &
       WRF_truelat1,WRF_truelat2,u,v)
  deallocate(work_WRF1)
  deallocate(work_WRF2)
  !
  ! Rotate wind if request (additional rotation)
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
  !*** W(nx_WRF,ny_WRF,nzstag_WRF) --> w(nx,ny,nz)
  !
  allocate(work_WRF (nx_WRF,ny_WRF,nz_WRF    ))
  allocate(work_WRF1(nx_WRF,ny_WRF,nzstag_WRF))
  allocate(work_WRF2(nx_WRF,ny_WRF,nzstag_WRF))
  !
  if( nf90_inq_varid(ncID,w_WRF_name,wID) /= 0) call runend('read_WRF : Error getting wID')
  if( nf90_get_var(ncID,wID,work_WRF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading w')
  if( nf90_get_var(ncID,wID,work_WRF2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading w')
  work_WRF1 = time_s*work_WRF1 + (1.0_rp-time_s)*work_WRF2   ! W arakawa at time
  !
  work_WRF = 0.5_rp*(work_WRF1(1:nx_WRF,1:ny_WRF,1:nzstag_WRF-1)+ work_WRF1(1:nx_WRF,1:ny_WRF,2:nzstag_WRF)) ! w box
  !
  call nc_interpola3d(nx_WRF,ny_WRF,nz_WRF,nelem_WRF, &
       work_WRF,z_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,nz,npoin_DAT,w,zlayer)
  !
  deallocate(work_WRF )
  deallocate(work_WRF1)
  deallocate(work_WRF2)
  !
  !*** Rest of 3D scalar variables. These are already given at the box center
  !
  allocate(work_WRF (nx_WRF,ny_WRF,nz_WRF))
  allocate(work_WRF1(nx_WRF,ny_WRF,nz_WRF))
  allocate(work_WRF2(nx_WRF,ny_WRF,nz_WRF))
  !
  !*** p(nx_WRF,ny_WRF,nz_WRF) --> p(nx,ny,nz) Pressure (reference plus perturbed)
  !
  if( nf90_inq_varid(ncID,p_WRF_name,pID) /= 0) call runend('read_WRF : Error getting pID')
  if( nf90_get_var(ncID,pID,work_WRF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading p')
  if( nf90_get_var(ncID,pID,work_WRF2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading p')
  work_WRF = time_s*work_WRF1 + (1.0_rp-time_s)*work_WRF2   ! p preturbed box
  !
  if( nf90_inq_varid(ncID,pb_WRF_name,pbID) /= 0) call runend('read_WRF : Error getting pbID')
  if( nf90_get_var(ncID,pbID,work_WRF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading pb')
  if( nf90_get_var(ncID,pbID,work_WRF2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading pb')
  work_WRF = work_WRF + time_s*work_WRF1 + (1.0_rp-time_s)*work_WRF2   ! p+pb box
  !
  call nc_interpola3d(nx_WRF,ny_WRF,nz_WRF,nelem_WRF, &
       work_WRF,z_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,nz,npoin_DAT,p,zlayer)
  !
  !***  T(nx_WRF,ny_WRF,nz_WRF) --> Tp(nx,ny,nz) potential temperature and T(nx,ny,nz) temperature
  !***  It is also necessary to use pressure (stored in work_WRF) because T is given as theta temperature
  !***  (minus a reference value = 300)
  !***
  !***  T(:,:,:)=(Tp(:,:,:)+300)*((p(:,:,:)+pb(:,:,:))/pref_theta)**(r/cp)
  !
  if( nf90_inq_varid(ncID,t_WRF_name,tID) /= 0) call runend('read_WRF : Error getting tID')
  if( nf90_get_var(ncID,tID,work_WRF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading T')
  if( nf90_get_var(ncID,tID,work_WRF2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading T')
  work_WRF1 = time_s*work_WRF1 + (1.0_rp-time_s)*work_WRF2   ! T perturbated potential temperature box
  work_WRF1 = work_WRF1 + 300.0_rp                         ! T potential box
  !
  work_WRF =  work_WRF1*((work_WRF/1e5_rp)**(0.285_rp))    ! T box
  !
  call nc_interpola3d(nx_WRF,ny_WRF,nz_WRF,nelem_WRF, &
       work_WRF1,z_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,nz,npoin_DAT,Tp,zlayer)
  call nc_interpola3d(nx_WRF,ny_WRF,nz_WRF,nelem_WRF, &
       work_WRF,z_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,nz,npoin_DAT,T,zlayer)
  !
  !*** qv(nx_WRF,ny_WRF,nz_WRF) --> qv(nx,ny,nz) specific humidity. WRF gives the Water vapor mixing ratio (e),
  !*** so the conversion q = e/1+e is needed
  !
  if( nf90_inq_varid(ncID,qv_WRF_name,qvID) /= 0) call runend('read_WRF : Error getting qvID')
  if( nf90_get_var(ncID,qvID,work_WRF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading qv')
  if( nf90_get_var(ncID,qvID,work_WRF2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading qv')
  work_WRF = time_s*work_WRF1 + (1.0_rp-time_s)*work_WRF2   ! mixing ratio
  work_WRF = work_WRF/(1.0_rp + work_WRF)                 ! specific humidity
  !
  call nc_interpola3d(nx_WRF,ny_WRF,nz_WRF,nelem_WRF, &
       work_WRF,z_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,nz,npoin_DAT,qv,zlayer)
  !
  !*** Release memory
  !
  deallocate(work_WRF )
  deallocate(work_WRF1)
  deallocate(work_WRF2)
  !
  !*** 3d scalar variables depending on nsoil
  !
  allocate(work_WRF (nx_WRF,ny_WRF,1))
  allocate(work_WRF1(nx_WRF,ny_WRF,nsoil_WRF))
  allocate(work_WRF2(nx_WRF,ny_WRF,nsoil_WRF))
  !
  !***  SMOIS(nx_WRF,ny_WRF,1) soil moisture (in m3 m-3). Read only the first soil layer
  !
  if( nf90_inq_varid(ncID,smoi_WRF_name,smoiID) /= 0) call runend('read_WRF : Error getting smoiID')
  if( nf90_get_var(ncID,smoiID,work_WRF1,start=(/1,1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading smoi')
  if( nf90_get_var(ncID,smoiID,work_WRF2,start=(/1,1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading smoi')
  work_WRF(1:nx_WRF,1:ny_WRF,1) =         time_s *work_WRF1(1:nx_WRF,1:ny_WRF,1) + &
       (1.0_rp-time_s)*work_WRF2(1:nx_WRF,1:ny_WRF,1)   ! smoist box
  !
  call nc_interpola2d(nx_WRF,ny_WRF,nelem_WRF, &
       work_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,npoin_DAT,smoi)
  !
  !*** Release memory
  !
  deallocate(work_WRF )
  deallocate(work_WRF1)
  deallocate(work_WRF2)
  !
  !***  Rest of 2D scalar variables. These are already given at the box center
  !
  allocate(work_WRF (nx_WRF,ny_WRF,1))
  allocate(work_WRF1(nx_WRF,ny_WRF,1))
  allocate(work_WRF2(nx_WRF,ny_WRF,1))
  !
  !***  PBLH(nx_WRF,ny_WRF) Planetary boundary layer height
  !
  if( nf90_inq_varid(ncID,pblh_WRF_name,pblhID) /= 0) call runend('read_WRF : Error getting pblhID')
  if( nf90_get_var(ncID,pblhID,work_WRF1,start=(/1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading pblh')
  if( nf90_get_var(ncID,pblhID,work_WRF2,start=(/1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading pblh')
  work_WRF = time_s*work_WRF1 + (1.0_rp-time_s)*work_WRF2   ! pblh box
  !
  call nc_interpola2d(nx_WRF,ny_WRF,nelem_WRF, &
       work_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,npoin_DAT,pblh)
  !
  !***  UST(nx_WRF,ny_WRF) U* similarity theory (m/s)
  !
  if( nf90_inq_varid(ncID,ust_WRF_name,ustID) /= 0) call runend('read_WRF : Error getting ustID')
  if( nf90_get_var(ncID,ustID,work_WRF1,start=(/1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading ust')
  if( nf90_get_var(ncID,ustID,work_WRF2,start=(/1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading ust')
  work_WRF = time_s*work_WRF1 + (1.0_rp-time_s)*work_WRF2   ! ust box
  !
  call nc_interpola2d(nx_WRF,ny_WRF,nelem_WRF, &
       work_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,npoin_DAT,ust)
  !
  !*** RMOL(nx_WRF,ny_WRF) Inverse Monin-Obukhov length (1/m)
  !
  if( nf90_inq_varid(ncID,rmol_WRF_name,rmolID) /= 0) call runend('read_WRF :Error getting rmolID')
  if( nf90_get_var(ncID,rmolID,work_WRF1,start=(/1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading rmol')
  if( nf90_get_var(ncID,rmolID,work_WRF2,start=(/1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading rmol')
  work_WRF = time_s*work_WRF1 + (1.0_rp-time_s)*work_WRF2   ! rmol box
  !
  call nc_interpola2d(nx_WRF,ny_WRF,nelem_WRF, &
       work_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,npoin_DAT,rmol)
  !
  !*** ZNT(nx_WRF,ny_WRF) Roughness length (m)
  !
  if( nf90_inq_varid(ncID,znt_WRF_name,zntID) /= 0) call runend('read_WRF: Error getting zntID')
  if( nf90_get_var(ncID,zntID,work_WRF1,start=(/1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading znt')
  if( nf90_get_var(ncID,zntID,work_WRF2,start=(/1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading znt')
  work_WRF = time_s*work_WRF1 + (1.0_rp-time_s)*work_WRF2
  call nc_interpola2d(nx_WRF,ny_WRF,nelem_WRF, &
       work_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,npoin_DAT,znt)
  !
  !*** Wind speed at 10m
  !
  if( nf90_inq_varid(ncID,u10_WRF_name,u10ID) /= 0) call runend('read_WRF :Error getting u10ID')
  if( nf90_get_var(ncID,u10ID,work_WRF1,start=(/1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading u10')
  if( nf90_get_var(ncID,u10ID,work_WRF2,start=(/1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading u10')
  work_WRF = (time_s*work_WRF1 + (1.0_rp-time_s)*work_WRF2)**2   ! U10**2
  !
  if( nf90_inq_varid(ncID,v10_WRF_name,v10ID) /= 0) call runend('read_WRF :Error getting v10ID')
  if( nf90_get_var(ncID,v10ID,work_WRF1,start=(/1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading v10')
  if( nf90_get_var(ncID,v10ID,work_WRF2,start=(/1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading v10')
  work_WRF = work_WRF + (time_s*work_WRF1 + (1.0_rp-time_s)*work_WRF2)**2 !U10**2 + V10**2
  work_WRF = SQRT(work_WRF)
  !
  call nc_interpola2d(nx_WRF,ny_WRF,nelem_WRF, &
       work_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,npoin_DAT,spd10)
  !
  !***  HFX(nx_WRF,ny_WRF) Surface heat flux (W/m2)
  !***  Once interpolated is stored in L(nx,ny) to avoid allocating
  !***  unnecessary memory
  !
  if( nf90_inq_varid(ncID,hfx_WRF_name,hfxID) /= 0) call runend('read_WRF : Error getting hfxID')
  if( nf90_get_var(ncID,hfxID,work_WRF1,start=(/1,1,itime1/)) /=0 ) call runend('read_WRF: Error reading hfx')
  if( nf90_get_var(ncID,hfxID,work_WRF2,start=(/1,1,itime2/)) /=0 ) call runend('read_WRF: Error reading hfx')
  work_WRF = time_s*work_WRF1 + (1.0_rp-time_s)*work_WRF2   ! hfl box
  !
  call nc_interpola2d(nx_WRF,ny_WRF,nelem_WRF, &
       work_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,npoin_DAT,L)
  !
  deallocate(work_WRF )
  deallocate(work_WRF1)
  deallocate(work_WRF2)
  !
  !*** Other derived 3D variables and transformations (already in the dbs grid!)
  !
  allocate(work_WRF (nx,ny,nz))
  !
  !*** Tv(nx,ny,nz) virtual temperature
  !
  work_WRF = qv(:,:,:)/(1.0_rp-qv(:,:,:))      ! Mixing ratio (=humidity ratio)
  !
  Tv(:,:,:) = T(:,:,:)*(1.0_rp + 1.6077_rp*work_WRF(:,:,:))/(1.0_rp+work_WRF(:,:,:))
  !
  !*** ro(nx,ny,nz) density. Gas law using virtual temperature (account for water in air)
  !
  ro(:,:,:) = p(:,:,:)/(287.06*Tv(:,:,:))
  !
  deallocate(work_WRF )
  !
  !***  Closes
  !
  if( nf90_close(ncID) /= 0) call runend('read_WRF : Error in nf90_close')
  !
  !***  Computes other variables not given by WRF
  !
  !
  !***  L(nx,ny)
  !
  do iy = 1,ny
     do ix = 1,nx
        L(ix,iy) =  - (ust(ix,iy)*ust(ix,iy)*ust(ix,iy)*Tp(ix,iy,1)*ro(ix,iy,1)*1006.0)/ &
             (0.4*9.81*L(ix,iy))
        L(ix,iy) = min(L(ix,iy), 1e3_rp)
        L(ix,iy) = max(L(ix,iy),-1e3_rp)
     end do
  end do

  !
  !*** ust(nx,ny) and L(nx,ny). Estimated
  !
  !     do iy = 1,ny
  !        do ix = 1,nx
  !           call get_par_ABL(1.0_rp ,zlayer(2),T(ix,iy,1),T(ix,iy,2),P(ix,iy,1),P(ix,iy,2), &
  !                            u(ix,iy,1),v(ix,iy,1),ust(ix,iy),L(ix,iy))
  !        end do
  !     end do
  !
  return
end subroutine read_WRF_data
!
!
!
subroutine rotate_WRF_wind(u,v,d1,d2,d3,            &
     map_proj,cen_lon,xlat,xlon,     &
     truelat1,truelat2,urot,vrot)
  !***************************************************************
  !*
  !*  Rotates WRF wind
  !*
  !***************************************************************
  use KindType, ONLY :  ip,rp
  implicit none
  integer(ip), intent(in) ::  d1,d2,d3,map_proj
  real(rp)  :: u(d1,d2,d3),v(d1,d2,d3)
  real(rp)  :: xlat(d1,d2),xlon(d1,d2)
  real(rp)  :: urot(d1,d2,d3),vrot(d1,d2,d3)
  real(rp)  :: cen_lon, truelat1, truelat2
  !
  integer(ip) :: i,j,k
  real(rp)    :: cone
  real(rp), allocatable :: diff(:,:), alpha(:,:)
  !
  real(rp), parameter :: pii = 3.14159265_rp
  real(rp), parameter :: radians_per_degree = pii/180_rp
  !
  !*** allocate
  !
  allocate(diff(d1,d2))
  allocate(alpha(d1,d2))
  !
  cone =1.0_rp                          !  PS
  if( map_proj .eq. 1) then                          !  Lambert Conformal mapping
     IF (ABS(truelat1-truelat2) .GT. 0.1_rp) THEN
        cone=(LOG(COS(truelat1*radians_per_degree))-            &
             LOG(COS(truelat2*radians_per_degree))) /          &
             (LOG(TAN((90_rp-ABS(truelat1))*radians_per_degree*0.5_rp ))- &
             LOG(TAN((90_rp-ABS(truelat2))*radians_per_degree*0.5_rp )) )
     ELSE
        cone = SIN(ABS(truelat1)*radians_per_degree )
     ENDIF
  end if
  !
  diff(:,:) = xlon(:,:) - cen_lon
  !
  do i = 1, d1
     do j = 1, d2
        if(diff(i,j) .gt. 180_rp) then
           diff(i,j) = diff(i,j) - 360_rp
        end if
        if(diff(i,j) .lt. -180_rp) then
           diff(i,j) = diff(i,j) + 360_rp
        end if
     end do
  end do
  !
  do i = 1, d1
     do j = 1, d2
        if(xlat(i,j) .lt.0.0_rp) then
           alpha(i,j) = - diff(i,j) * cone * radians_per_degree
        else
           alpha(i,j) = diff(i,j) * cone * radians_per_degree
        end if
     end do
  end do
  !
  do k=1,d3
     urot(1:d1,1:d2,k) = v(1:d1,1:d2,k)*sin(alpha(1:d1,1:d2)) + u(1:d1,1:d2,k)*cos(alpha(1:d1,1:d2))
     vrot(1:d1,1:d2,k) = v(1:d1,1:d2,k)*cos(alpha(1:d1,1:d2)) - u(1:d1,1:d2,k)*sin(alpha(1:d1,1:d2))
  end do
  !
  deallocate(diff)
  deallocate(alpha)
  !
  return
end subroutine rotate_WRF_wind
