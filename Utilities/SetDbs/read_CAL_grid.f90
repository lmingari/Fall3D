subroutine read_CAL_grid (fname,npoin,xutm,yutm)
  !********************************************************************
  !*
  !*   Reads CALMET grid. It also creates all the
  !*   interpolation arrays used later by read_CAL_data
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use TimeFun
  use MathFun
  use CAL_nc
  implicit none
  !
  character(len=*) :: fname
  integer(ip)      :: npoin
  real   (rp)      :: xutm(npoin),yutm(npoin)
  !
  logical     :: go_on
  integer(ip) :: iyr,imo,idy,ihr,imi,ise
  integer(ip) :: ix,iy,it
  !
  !*** Initializations
  !
  npoin_DAT = npoin
  !
  !*** First of all, read all calmet data. All necessary data is stored in memory
  !*** (all time steps). Contrary to a netCDF file, CALMET file is binary sequential, so
  !*** it is not convenient computationally to scan all the file each time data
  !*** is accessed.
  !
  open(ludat,file=TRIM(fname),status='old',err=100,form='unformatted')
  !
  !*** Read CALMET headers
  !
  call readcal620
  !
  !*** Determines the CAL grid dimensions in 4D
  !
  ibyr_CAL  = icbyr
  ibmo_CAL  = icbmo
  ibdy_CAL  = icbdy
  ibhr_CAL  = icbhr
  ibsec_CAL = icbsec
  !
  nx_CAL = nxx
  ny_CAL = nyy
  nz_CAL = nzl + 1 ! Calmet gives values at the cell center. The value z=0 (ground) is added
  !
  !***  Averiguates the number of time intervals
  !***  WARNING: One hour time step is assumed
  !
  nt_CAL = 0
  dt_CAL = 3600. ! 1h assumed
  go_on = .true.
  do while(go_on)
     nt_CAL = nt_CAL + 1
     call addtime(ibyr_CAL,ibmo_CAL,ibdy_CAL,ibhr_CAL, &
          iyr,imo,idy,ihr,imi,ise,nt_CAL*dt_CAL)
     if( (iyr.eq.iceyr).and.(imo.eq.icemo).and.(idy.eq.icedy).and.(ihr.eq.icehr) ) go_on = .false.
     if( nt_CAL.gt.24*31 ) call runend('read_CAL_grid : unable to find dt_CAL')  ! 1 month max
  end do
  !
  !*** Determinates the nodal conectivities of the CALMET grid. This is necessary
  !*** in order to interpolate data onto a GENERIC set of points (which do not
  !*** necessary form a structurated mesh!).
  !
  nelem_CAL = (nx_CAL - 1)*(ny_CAL - 1)
  !
  !***  Allocates memory
  !
  call CAL_getmem
  !
  !*** Reads the rest of file
  !
  call readcal62
  close(ludat)
  !
  !*** Calculates timesec_CAL(:) in sec after 0000UTC
  !
  timesec_CAL(1) = 3600.*ibhr_CAL + ibsec_CAL + dt_CAL  ! First interval
  do it = 2,nt_CAL
     timesec_CAL(it) = timesec_CAL(1) + (it-1)*dt_CAL
  end do
  !
  !***  Calculates time_CAL(:) in format YYYYMMDDHHMMSS
  !
  do it = 1,nt_CAL
     call addtime(ibyr_CAL,ibmo_CAL,ibdy_CAL,ibhr_CAL, &
          iyr,imo,idy,ihr,imi,ise,it*dt_CAL)
     time_CAL(it) = ise*1d0 + &
          imi*1d2 + &
          ihr*1d4 + &
          idy*1d6 + &
          imo*1d8 + &
          iyr*1d10
  end do
  !
  !*** Builds grid coordinates. This is necessari to interpolate
  !
  xo_CAL = DBLE(xorigr4)
  yo_CAL = DBLE(yorigr4)
  dx_CAL = DBLE(dgrid4)
  dy_CAL = DBLE(dgrid4)
  do ix = 1,nx_CAL
     xutm_CAL(ix,1:ny_CAL) = xo_CAL + (ix-1)*dx_CAL
  end do
  do iy = 1,ny_CAL
     yutm_CAL(1:nx_CAL,iy) = yo_CAL + (iy-1)*dy_CAL
  end do
  !
  !*** Builds the nodal connectivities
  !
  call set_lnods(nx_CAL,ny_CAL,nelem_CAL,lnods_CAL)
  !
  !*** Builds elements to points connectivities and computes the interpolation
  !*** factors in 2D. It also verifies that the points lay within the WRF domain.
  !
  call set_lelpo(nelem_CAL,npoin_DAT,xutm_CAL,yutm_CAL,xutm,yutm,lnods_CAL, &
       lelpo_CAL,s_po_CAL,t_po_CAL)
  !
  return
  !
100 call runend('Error opening calmet62 file '//TRIM(fname))
  !
end subroutine read_CAL_grid
!
!
!
subroutine readcal620
  !*****************************************************************************
  !*
  !*    This subroutine reads (and skips) SCALAR data from the a CALMET 6.2 meteo
  !*    binary file.
  !*
  !***************************************************************************
  use CAL_nc
  implicit none
  !
  !***  Record #1 - File Declaration -- 24 words
  !
  read(ludat) dataset,dataver,datamod
  !
  !***  Record #2 - Number of comment lines -- 1 word
  !
  read(ludat) ncom
  if(ncom.gt.mxncom) call runend('Readcal620: too many comments')
  !
  !***  Record #3 to NCOM+2 (Comment record section) -- 33 words each
  !
  do icom=1,ncom
     read(ludat) comment(icom)
  end do
  !
  !***  Record #NCOM+3 - run control parameters -- 39 words
  !
  read(ludat) icbyr,icbmo,icbdy,icbhr,icbsec,                   &
       iceyr,icemo,icedy,icehr,icesec,                   &
       axtz,irlg,irtype,                                 &
       nxx,nyy,nzl,dgrid4,xorigr4,yorigr4,iwfcod,nssta,    &
       nusta,npsta,nowsta,nlu,iwat1,iwat2,lcalgrd,       &
       pmap,datum,daten,feast,fnorth,utmhem,iutmzn,      &
       rnlat0,relon0,xlat1,xlat2
  !
  !***  Builds ndatcb and ndatce in the format YYYYMMDDHHSSSS
  !
  write(ndatcb(1 :4 ),'(i4.4)') icbyr
  write(ndatcb(5 :6 ),'(i2.2)') icbmo
  write(ndatcb(7 :8 ),'(i2.2)') icbdy
  write(ndatcb(9 :10),'(i2.2)') icbhr
  write(ndatcb(11:14),'(i4.4)') icbsec
  !
  write(ndatce(1 :4 ),'(i4.4)') iceyr
  write(ndatce(5 :6 ),'(i2.2)') icemo
  write(ndatce(7 :8 ),'(i2.2)') icedy
  write(ndatce(9 :10),'(i2.2)') icehr
  write(ndatce(11:14),'(i4.4)') icesec
  !
  return
end subroutine readcal620
!
!
!
subroutine readcal62
  !******************************************************************************
  !*
  !*    This subroutine reads VECTOR data from the a CALMET 6.2 meteo binary
  !*    file. real*4 data is always converted to real*8 type
  !*    It must be called after readcal620
  !*
  !*    OUTPUTS
  !*
  !*    z_CAL   (nx_CAL,ny_CAL,nz_CAL)       = z-layer coordinates
  !*    LDU_CAL (nx_CAL,ny_CAL       )       = Land use (integer)
  !*    z0_CAL  (nx_CAL,ny_CAL       )       = Roughness length
  !*    topg_CAL(nx_CAL,ny_CAL       )       = topography
  !*    u_CAL   (nx_CAL,ny_CAL,nz_CAL,nt_CAL)= u-vel
  !*    v_CAL   (nx_CAL,ny_CAL,nz_CAL,nt_CAL)= v-vel
  !*    w_CAL   (nx_CAL,ny_CAL,nz_CAL,nt_CAL)= w-vel
  !*    T_CAL   (nx_CAL,ny_CAL,nz_CAL,nt_CAL)= air temperature
  !*    UST_CAL (nx_CAL,ny_CAL,       nt_CAL)= u*
  !*    PBL_CAL (nx_CAL,ny_CAL,       nt_CAL)= Mixing height
  !*    L_CAL   (nx_CAL,ny_CAL,       nt_CAL)= Monin-Obukhov length
  !*
  !*    All the units are in MKS
  !*
  !*****************************************************************
  use CAL_nc
  implicit none
  !
  integer(ip) :: iz,it
  !
  !***  Record #NCOM+4 - cell face heights (NZL + 1 = NZ words)
  !
  call readrrec(ludat,clab1,idum,idum,idum,idum,work4,nz_CAL)  ! ZLAYER(nz=nzl+1) in real*4
  call v4tov8  (work8,work4,nz_CAL)                            ! ZLAYER           in real*8
  do iz = 1,nz_CAL
     z_CAL(1:nx_CAL,1:ny_CAL,iz) = work8(iz)
  end do
  !
  !***  Redefine zlayers.First layer (z=0) is not shifted
  !
  do iz = nz_CAL,2,(-1)
     z_CAL(1:nx_CAL,1:ny_CAL,iz) = 0.5d0*(z_CAL(1:nx_CAL,1:ny_CAL,iz) + z_CAL(1:nx_CAL,1:ny_CAL,iz-1))
  end do
  !
  !***  Records #NCOM+5 & 6 - x, y coordinates of surface stations
  !*** (NSSTA words each record)
  !
  if(nssta.ge.1) then
     call readrdum(ludat,clabdu,nssta) ! XSSTA
     call readrdum(ludat,clabdu,nssta) ! YSSTA
  endif
  !
  !***  Records #NCOM+7 & 8 - x, y coordinates of upper air stations
  !***  (NUSTA words each record)
  !
  if(nusta.ge.1) then
     call readrdum(ludat,clabdu,nusta) ! XUSTA
     call readrdum(ludat,clabdu,nusta) ! YUSTA
  endif
  !
  !***  Records #NCOM+9 & 10 - x, y coordinates of precipitation stations
  !***  (NPSTA words each record)
  !
  if(npsta.ge.1) then
     call readrdum(ludat,clabdu,npsta) ! XPSTA
     call readrdum(ludat,clabdu,npsta) ! YPSTA
  endif
  !
  !***  Record #NCOM+11 - surface roughness lengths (NX * NY words)
  !
  call readrrec(ludat,clabdu,idum,idum,idum,idum,work4,nx_CAL*ny_CAL)
  call v4tov8  (Z0_CAL,work4,nx_CAL*ny_CAL)   !  Z0 in real*8
  !
  !***  Record #NCOM+12 - land use categories (NX * NY words)
  !***  Integer array
  !
  call readirec(ludat,clabdu,idum,idum,idum,idum,LDU_CAL,nx_CAL*ny_CAL)
  !
  !***  Record #NCOM+13 - topography (NX * NY words)
  !
  call readrrec(ludat,clab10,idum,idum,idum,idum,work4,nx_CAL*ny_CAL)   ! TOPG in real*4
  call v4tov8  (topg_CAL,work4,nx_CAL*ny_CAL)                           ! TOPG in real*8
  !
  !***  Record #NCOM+14 - leaf area index (NX * NY words)
  !
  call readrdum(ludat,clabdu,nx_CAL*ny_CAL)
  !
  !***  Record #NCOM+15 - nearest surface station no. to each
  !***  grid point (NX * NY words)
  !
  if(nssta.ge.1)then
     call readidum(ludat,clabdu,nx_CAL*ny_CAL) ! NEARS
  endif
  !
  !***  READ DATA RECORDS (normally ervery hour)
  !
  do it = 1,nt_CAL
     !
     !***    Wind components
     !
     do iz = 1,nz_CAL-1         ! NZL = NZ-1
        call readrrec(ludat,clabu,idatc0,isec0,idatc1,isec1,work4,nx_CAL*ny_CAL)    ! VX in real*4
        call v4tov8  (u_CAL(1,1,iz+1,it),work4,nx_CAL*ny_CAL)
        !
        call readrrec(ludat,clabv,idatc0,isec0,idatc1,isec1,work4,nx_CAL*ny_CAL)    ! VY in real*4
        call v4tov8  (v_CAL(1,1,iz+1,it),work4,nx_CAL*ny_CAL)
        !
        if(lcalgrd) then
           call readrrec(ludat,clabw,idatc0,isec0,idatc1,isec1,work4,nx_CAL*ny_CAL)   ! VZ in real*4
           call v4tov8(w_CAL(1,1,iz+1,it),work4,nx_CAL*ny_CAL)
        end if
     end do
     !
     !***    Read the 3-D temperature field (LCALGRD run only)
     !
     if(lcalgrd.and.irtype.eq.1) then
        do iz=1,nz_CAL-1
           call readrrec(ludat,clabt,idatc0,isec0,idatc1,isec1,work4,nx_CAL*ny_CAL)   ! TEMPE in real*4
           call v4tov8  (T_CAL(1,1,iz+1,it),work4,nx_CAL*ny_CAL)
        end do
     endif
     !
     !***   Read other 2-D meteorological fields:
     !***   PGT stability class, Friction velocity (m/s), Mixing height (m),
     !***   Monin-Obukhov length (m), Convective velocity scale (m/s),
     !***   Precip. rate (mm/hr)
     !
     if(irtype.eq.1) then
        call readidum(ludat,clabdu,nx_CAL*ny_CAL)               ! IPGT
        !
        call readrrec(ludat,clabus,idatc0,isec0,idatc1,isec1,work4,nx_CAL*ny_CAL)  ! USTAR in real*4
        call v4tov8(UST_CAL(1,1,it),work4,nx_CAL*ny_CAL)
        !
        call readrrec(ludat,clabzi,idatc0,isec0,idatc1,isec1,work4,nx_CAL*ny_CAL)  ! PBLH in real*4
        call v4tov8(PBL_CAL(1,1,it),work4,nx_CAL*ny_CAL)
        !
        call readrrec(ludat,clabrl,idatc0,isec0,idatc1,isec1,work4,nx_CAL*ny_CAL)  ! L in real*4
        call v4tov8(L_CAL(1,1,it),work4,nx_CAL*ny_CAL)
        !
        call readrdum(ludat,clabdu,nx_CAL*ny_CAL)               ! WSTAR
        !
        if(npsta.ne.0) then
           call readrdum(ludat,clabdu,nx_CAL*ny_CAL)            ! RMM
        endif
        !
        call readrrec(ludat,clabtk,idatc0,isec0,idatc1,isec1,work4,nx_CAL*ny_CAL)  ! TEMPK(nx,ny) in real*4
        call v4tov8(TSFC_CAL(1,1,it),work4,nx_CAL*ny_CAL)
        !
        call readrdum(ludat,clabdu,nx_CAL*ny_CAL) ! RHO
        call readrdum(ludat,clabdu,nx_CAL*ny_CAL) ! QSW
        call readidum(ludat,clabdu,nx_CAL*ny_CAL) ! IRH
        if(npsta.ne.0) then
           call readidum(ludat,clabdu,nx_CAL*ny_CAL) ! IPCODE
        endif
     endif
     !
  end do ! it = 1,nt_CAL
  !
  !***  Other computations
  !
  u_CAL(1:nx_CAL,1:ny_CAL,1,1:nt_CAL) = 0.  ! Set vx(z=0)=0
  v_CAL(1:nx_CAL,1:ny_CAL,1,1:nt_CAL) = 0.  ! Set vy(z=0)=0
  w_CAL(1:nx_CAL,1:ny_CAL,1,1:nt_CAL) = 0.  ! Set vz(z=0)=0
  T_CAL(1:nx_CAL,1:ny_CAL,1,1:nt_CAL) = TSFC_CAL(1:nx_CAL,1:ny_CAL,1:nt_CAL)  ! Set T T(z=0)
  !
  return
end subroutine readcal62
!
!
!
!**********************************************************************
!*
!*    Routines to read integer/real records
!*
!**********************************************************************
subroutine readrrec(ludat,clab,ndatb,isec0,ndate,isec1,a,n)
  implicit none
  integer(4)  :: ludat,ndatb,isec0,ndate,isec1,n,i
  character(len=8) :: clab
  real(4)     :: a(*)
  read(ludat) clab,ndatb,isec0,ndate,isec1,(a(i),i=1,n)
end subroutine readrrec
!**********************************************************************
subroutine readirec(ludat,clab,ndatb,isec0,ndate,isec1,a,n)
  implicit none
  integer(4)  :: ludat,ndatb,isec0,ndate,isec1,n,i
  character(len=8) :: clab
  integer(4)  :: a(*)
  read(ludat) clab,ndatb,isec0,ndate,isec1,(a(i),i=1,n)
end subroutine readirec
!**********************************************************************
subroutine readrdum(ludat,clab,n)
  implicit none
  integer(4)  :: ludat
  character(len=8) :: clab
  integer(4)  :: n,i,idummy
  real(4)     :: dummy
  read(ludat) clab,idummy,idummy,idummy,idummy,(dummy,i=1,n)
end subroutine readrdum
!**********************************************************************
subroutine readidum(ludat,clab,n)
  implicit none
  integer(4)  :: ludat
  character(len=8) :: clab
  integer(4)  :: n,i,idummy
  read(ludat) clab,idummy,idummy,idummy,idummy,(idummy,i=1,n)
end subroutine readidum
!
subroutine v4tov8(v8,v4,n)
  !****************************************************************
  !*
  !*    Converts a real*4 vector into a real*8
  !*
  !****************************************************************
  use kindType, ONLY: rp
  implicit none
  integer  :: n,i
  real(rp) :: v8(*)
  real(4)  :: v4(*)
  !
  if(rp.eq.8) then
     do i=1,n
        v8(i) = dble(v4(i))
     end do
  else
     do i=1,n
        v8(i) = v4(i)
     end do
  end if
  !
  return
end subroutine v4tov8
