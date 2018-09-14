subroutine read_CAL_data(fname,time_lag,time,timesec,nx,ny,nz,zlayer,lon,lat,topg,LDU, &
     u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
  !********************************************************************
  !*
  !*   Interpolates CALMET data
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use InpOut
  use MathFun
  use CAL_nc
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
  real(rp)         :: time_s,pp,TT,rhoo
  !
  !*** Finds the time itime1 and itime2 and the interpolation factor time_s
  !
  if(nt_CAL.gt.1) then
     found  = .false.
     it_CAL = 0
     do while(.not.found)
        it_CAL = it_CAL + 1
        if( (it_CAL+1).gt.nt_CAL ) call runend('read_CAL : Time not found')
        if( (time_lag+timesec.le.timesec_CAL(1)) ) then
           found = .true.
           itime1 = 1
           itime2 = 1
           time_s = 1.0_rp
        end if
        if( (time_lag+timesec.gt.timesec_CAL(it_CAL  )).and. &
             (time_lag+timesec.le.timesec_CAL(it_CAL+1)) ) then
           found = .true.
           itime1 = it_CAL
           itime2 = it_CAL + 1
           time_s = (timesec_CAL(itime2)-time_lag-timesec)/(timesec_CAL(itime2)-timesec_CAL(itime1))
        end if
     end do
  else
     itime1 = 1
     itime2 = 1
     time_s = 1.0_rp
  end if
  !
  write(lulog,10) time,time_CAL(itime1),time_CAL(itime2),time_s
10 format('Interpolating at ',f15.0,/, &
       '        CALMET from ',f15.0,' to ',f15.0,' s = ',f6.3)
  !
  !*** Interpolates
  !
  !
  !***  LDU_CAL(nx_CAL,ny_CAL) --> ldu USGS Land Use Category (time independent actually)
  !
  call nc_assign2d(nx_CAL,ny_CAL,nelem_CAL, &
       1.0_rp*LDU_CAL,lelpo_CAL,lnods_CAL,s_po_CAL,t_po_CAL, &
       nx,ny,npoin_DAT,LDU)
  !
  !***  Convert from the CALMET 52-category system to the USGS 24-category
  !
  do iy = 1,ny
     do ix = 1,nx
        if(LDU(ix,iy).ge.10.and.LDU(ix,iy).le.17) then
           LDU(ix,iy) = 1
        else if(LDU(ix,iy).ge.20.and.LDU(ix,iy).le.24) then
           LDU(ix,iy) = 2
        else if(LDU(ix,iy).ge.-20.and.LDU(ix,iy).le.-24) then
           LDU(ix,iy) = 3
        else if(LDU(ix,iy).ge.30.and.LDU(ix,iy).le.33) then
           LDU(ix,iy) = 8
        else if(LDU(ix,iy).eq.40) then
           LDU(ix,iy) = 15
        else if(LDU(ix,iy).eq.41) then
           LDU(ix,iy) = 11
        else if(LDU(ix,iy).eq.42) then
           LDU(ix,iy) = 13
        else if(LDU(ix,iy).eq.43) then
           LDU(ix,iy) = 15
        else if(LDU(ix,iy).ge.50.and.LDU(ix,iy).le.55) then
           LDU(ix,iy) = 16
        else if(LDU(ix,iy).ge.60.and.LDU(ix,iy).le.61) then
           LDU(ix,iy) = 18
        else if(LDU(ix,iy).eq.62) then
           LDU(ix,iy) = 17
        else if(LDU(ix,iy).ge.70.and.LDU(ix,iy).le.77) then
           LDU(ix,iy) = 19
        else if(LDU(ix,iy).ge.80.and.LDU(ix,iy).le.81) then
           LDU(ix,iy) = 22
        else if(LDU(ix,iy).eq.82) then
           LDU(ix,iy) = 20
        else if(LDU(ix,iy).eq.83) then
           LDU(ix,iy) = 23
        else if(LDU(ix,iy).eq.84) then
           LDU(ix,iy) = 23
        else if(LDU(ix,iy).ge.90.and.LDU(ix,iy).le.92) then
           LDU(ix,iy) = 24
        else
           LDU(ix,iy) = -99
        end if
     end do
  end do
  !
  !*** topg_CAL(nx_CAL,ny_CAL) --> topg(nx,ny) Topograpgy (time independent actually)
  !
  call nc_interpola2d(nx_CAL,ny_CAL,nelem_CAL, &
       topg_CAL,lelpo_CAL,lnods_CAL,s_po_CAL,t_po_CAL, &
       nx,ny,npoin_DAT,topg)
  !
  !*** Time-dependent variables
  !
  !
  !*** U_CAL(nx_CAL,ny_CAL,nz_CAL) --> u(nx,ny,nz)
  !
  call nc_interpola3d(nx_CAL,ny_CAL,nz_CAL,nelem_CAL, &
       u_CAL(1,1,1,it_CAL),z_CAL,lelpo_CAL,lnods_CAL,s_po_CAL,t_po_CAL, &
       nx,ny,nz,npoin_DAT,u,zlayer)
  !
  !*** v_CAL(nx_CAL,ny_CAL,nz_CAL) --> v(nx,ny,nz)
  !
  call nc_interpola3d(nx_CAL,ny_CAL,nz_CAL,nelem_CAL, &
       v_CAL(1,1,1,it_CAL),z_CAL,lelpo_CAL,lnods_CAL,s_po_CAL,t_po_CAL, &
       nx,ny,nz,npoin_DAT,v,zlayer)
  !
  !*** w_CAL(nx_CAL,ny_CAL,nz_CAL) --> w(nx,ny,nz)
  !
  call nc_interpola3d(nx_CAL,ny_CAL,nz_CAL,nelem_CAL, &
       w_CAL(1,1,1,it_CAL),z_CAL,lelpo_CAL,lnods_CAL,s_po_CAL,t_po_CAL, &
       nx,ny,nz,npoin_DAT,w,zlayer)
  !
  !*** T_CAL(nx_CAL,ny_CAL,nz_CAL) --> T(nx,ny,nz)
  !
  call nc_interpola3d(nx_CAL,ny_CAL,nz_CAL,nelem_CAL, &
       T_CAL(1,1,1,it_CAL),z_CAL,lelpo_CAL,lnods_CAL,s_po_CAL,t_po_CAL, &
       nx,ny,nz,npoin_DAT,T,zlayer)
  !
  !*** UST_CAL(nx_CAL,ny_CAL) --> ust(nx,ny) ust*
  !
  call nc_interpola2d(nx_CAL,ny_CAL,nelem_CAL, &
       ust_CAL(1,1,it_CAL),lelpo_CAL,lnods_CAL,s_po_CAL,t_po_CAL, &
       nx,ny,npoin_DAT,ust)
  !
  !*** PBL_CAL(nx_CAL,ny_CAL) --> pbl(nx,ny) ust*
  !
  call nc_interpola2d(nx_CAL,ny_CAL,nelem_CAL, &
       pbl_CAL(1,1,it_CAL),lelpo_CAL,lnods_CAL,s_po_CAL,t_po_CAL, &
       nx,ny,npoin_DAT,pblh)
  !
  !*** L_CAL(nx_CAL,ny_CAL) --> L(nx,ny) ust*
  !
  call nc_interpola2d(nx_CAL,ny_CAL,nelem_CAL, &
       L_CAL(1,1,it_CAL),lelpo_CAL,lnods_CAL,s_po_CAL,t_po_CAL, &
       nx,ny,npoin_DAT,L)
  !
  !*** Rest variables not given by CALMET
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
  !*** qv(nx,ny)
  !
  qv = -999
  !
  !*** Standard density and pressure assumed
  !*** Potentail temperature computed
  !
  do iz = 1,nz
     do iy = 1,ny
        do ix = 1,nx
           call standard_atm(pp,TT,rhoo,topg(ix,iy)+zlayer(iz))
           p  (ix,iy,iz) = pp
           ro (ix,iy,iz) = rhoo
           Tp (ix,iy,iz) = T(ix,iy,iz)*((1.01d5/pp)**(0.285))
        end do
     end do
  end do
  !
  return
end subroutine read_CAL_data
