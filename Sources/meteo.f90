  subroutine meteo
  !*****************************************************************
  !*
  !*    This routine reads meteorological data.
  !*
  !*****************************************************************
  use KindType
  use InpOut
  use Numeric
  use Master
  use Parallel
  implicit none
  !
  logical      :: vver_correction
  integer(ip)  :: iyr,imo,idy,ihr,imi,ise,ic,k
  !
  !***  Proceed ?
  !
  if(.not.meteotime) return
  !
  !*** Reads meteo data
  !
  call reamet
  !
  !***  Calculate the diffusion coefficients
  !
  call kappa3
  !
  !***  Writes meteo information to the log file
  !
  call addtime(ibyr,ibmo,ibdy,0,iyr,imo,idy,ihr,imi,ise,time)
  !
  if( mpime == root ) then
     write(nlst,10) idy ,month(imo ),iyr ,ihr ,imi ,ise,time
10   format(/,                                                                        &
          'METEO FILE READ               ',/,                                         &
          '  Current time              : ',i2,1x,a3,1x,i4,' at ',i2,':',i2,':',i2,/,  &
          '  Time (sec after 00UTC)    : ',f10.0)
     !
     write(nlst,20) maxval(vx)   ,minval(vx)             ,sum(vx)/(nx*ny*nz),                  &
          maxval(vy)             ,minval(vy)             ,sum(vy)/(nx*ny*nz),                  &
          maxval(vz)             ,minval(vz)             ,sum(vz)/(nx*ny*nz),                  &
          maxval(tempe)-273.     ,minval(tempe)-273.     ,sum(tempe)/(nx*ny*nz)-273.,          &
          maxval(rho)            ,minval(rho)            ,sum(rho)/(nx*ny*nz),                 &
          maxval(rkhor(:,:,2:nz)),minval(rkhor(:,:,2:nz)),sum(rkhor(:,:,2:nz))/(nx*ny*(nz-1)), &
          maxval(rkver(:,:,2:nz)),minval(rkver(:,:,2:nz)),sum(rkver(:,:,2:nz))/(nx*ny*(nz-1))
20   format( /,                                     &
          '  X-velocity values           ',/,       &
          '    Maximum                 : ',f12.4,/, &
          '    Minimum                 : ',f12.4,/, &
          '    Average                 : ',f12.4,/, &
          '  Y-velocity values           ',/,       &
          '    Maximum                 : ',f12.4,/, &
          '    Minimum                 : ',f12.4,/, &
          '    Average                 : ',f12.4,/, &
          '  Z-velocity values           ',/,       &
          '    Maximum                 : ',f12.4,/, &
          '    Minimum                 : ',f12.4,/, &
          '    Average                 : ',f12.4,/, &
          '  Temperature  C              ',/,       &
          '    Maximum                 : ',f12.4,/, &
          '    Minimum                 : ',f12.4,/, &
          '    Average                 : ',f12.4,/, &
          '  Air Density                 ',/,       &
          '    Maximum                 : ',f12.4,/, &
          '    Minimum                 : ',f12.4,/, &
          '    Average                 : ',f12.4,/, &
          '  Horizontal diffusion        ',/,       &
          '    Maximum                 : ',f12.4,/, &
          '    Minimum                 : ',f12.4,/, &
          '    Average                 : ',f12.4,/, &
          '  Vertical   diffusion        ',/,       &
          '    Maximum                 : ',f12.4,/, &
          '    Minimum                 : ',f12.4,/, &
          '    Average                 : ',f12.4)
  end if
  !
  !***  Calculate settling velocities of all classes. Gas classes are assumed
  !***  to have a zero settling velocity
  !
  call setvset
  !
  if( mpime == root ) then
     write(nlst,100) (ic,ic=1,nc)
     do k=1,nz
        write(nlst,110) k,(MAXVAL(vset(:,:,k,ic)),ic=1,nc)
        write(nlst,111)   (MINVAL(vset(:,:,k,ic)),ic=1,nc)
     end do
100  format(/, &
          'TERMINAL VELOCITIES (m/s)',/,  &
          '  Level/Class ',1x,50(8x,i2))
110  format(/,'   ',i2,'          max ',50(f9.4,1x))
111  format(  '   ',2x,'          min ',50(f9.4,1x))
  end if
  !
  !***  Compute here dry deposition velocities to account for effects
  !***  distinct from sedimentation velocity
  !
  v_drydep = .true.
  if(v_drydep) call setvdrydep      !    vset(1:nx,1:ny,1,1:nc) = vset(1:nx,1:ny,1,1:nc) + vdep(1:nx,1:ny,1:nc)
  if( mpime == root ) then
    write(nlst,*) "ZNT MIN MAX", minval(znt), maxval(znt)
  end if
  !
  !***  Computes scaled quantities (stored in the same arrays) except rkh1
  !***  because scaling breaks symmetry in horizontal diffusion
  !
  do k = 1,nz
     vx  (1:nx,1:ny,k) = vx   (1:nx,1:ny,k)/ Hm1(1:nx,1:ny)
     rho (1:nx,1:ny,k) = rho  (1:nx,1:ny,k)* Hm1(1:nx,1:ny)
     rkh1(1:nx,1:ny,k) = rkhor(1:nx,1:ny,k)/(Hm1(1:nx,1:ny)*Hm1(1:nx,1:ny))
  end do
  !
  !***  Correct vertical velocity (scaling coming from transformation to s-coordinates).
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     vver_correction = .true.  ! Modified                          
  case('UTM')
     vver_correction = .false.
  END SELECT
  if(vver_correction) then
     call corvver
  !  if( mpime == root ) then
  !     write(nlst,300) (ic,ic=1,nc)
  !     do k=1,nz
  !        write(nlst,310) k,(MAXVAL(vset(:,:,k,ic)),ic=1,nc)
  !        write(nlst,311)   (MINVAL(vset(:,:,k,ic)),ic=1,nc)
  !     end do
300     format(/, &
             'TERMINAL VELOCITIES IN SIGMA COORDINATES (m/s)',/,  &
             '  Level/Class ',1x,50(8x,i2))
310     format(/,'   ',i2,'          max ',50(f9.4,1x))
311     format(  '   ',2x,'          min ',50(f9.4,1x))
  !  end if
  end if
  !
  !*** For gravity currents, modify the wind field adding a null-divergence wind field
  !
  if(gravity_current) call addradialwind
  !
  !*** Computes the divergence of velocity. This is used to ensure mass-conservation
  !*** in non-null divergence velocity fields and must be called after addradialwind
  !
  call setdivu
  !
  !*** Calculates the critical time step (using already the scaled quantities)
  !
  call finddt
  !
  return
  end subroutine meteo
