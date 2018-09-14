  subroutine wridat
  !*****************************************************************
  !*
  !*    Writes input data to the lst file
  !*
  !*****************************************************************
  use KindType
  use InpOut
  use Numeric
  use Master
  implicit none
  !
  character(len=15) :: message
  integer(ip)       :: iyr,imo,idy,ihr,imi,ise,iz,ic,i
  !
  !***  Writes Grid data to the lst file
  !
  write(nlst,10)
10 format(/,                                                          &
       '----------------------------------------------------',/,   &
       '                                                    ',/,   &
       '               METEO AND GRID DATA                  ',/,   &
       '                                                    ',/,   &
       '----------------------------------------------------',/)
  !
  call addtime(idbsyr,idbsmo,idbsdy,0,iyr,imo,idy,ihr,imi,ise,1.0_rp*ibdbs)
  write(nlst,20) idy,month(imo),iyr,ihr,(imi*60)+ise
20 format('METEO DATA TIME RANGE',/,                      &
       '  Initial time       : ',2x,i2.2,1x,a3,1x,i4.4,' at ',i2.2,' h ',i4.4,' s')
  !
  call addtime(idbsyr,idbsmo,idbsdy,0,iyr,imo,idy,ihr,imi,ise,1.0_rp*iedbs)
  write(nlst,21) idy,month(imo),iyr,ihr,(imi*60)+ise
21 format('  Final   time       : ',2x,i2.2,1x,a3,1x,i4.4,' at ',i2.2,' h ',i4.4,' s')
  !
  write(nlst,22) DBLE(iedbs-ibdbs)/36d2,(iedbs-ibdbs)
22 format('  Meteo coverage     : ',f6.1,' h (',i9,' sec)')
  !
  write(nlst,23) INT(time_lag)/86400,INT(time_lag)
23 format('  Time lag           : ',i3,' days (',i9,' sec)',/)
  !
  SELECT CASE(coord_sys)
  case('UTM')
     write(nlst,30) TRIM(coord_sys), &
          xorigr,yorigr,xorigr+(nx-1)*dx,yorigr+(ny-1)*dy,nx,ny,dx,dy,  &
          tpgmin,tpgmax,ztop,nz
30   format('MESH',/, &
          '  System of coord.   :  ',a               ,/,      &
          '  Bottom-left corner : (',f9.1,2x,f9.1,')',/,      &
          '  Top-right   corner : (',f9.1,2x,f9.1,')',/,      &
          '  Number points x    : ',i9  ,/,                   &
          '  Number points y    : ',i9  ,/,                   &
          '  X  grid spacing    : ',f9.1,/,                   &
          '  Y  grid spacing    : ',f9.1,/,                   &
          '  Min. topography    : ',f9.1             ,/,      &
          '  Max. topography    : ',f9.1             ,/,      &
          '  Max. z domain      : ',f9.1             ,/,      &
          '  Number z-levels    : ',i9)
  case('LON-LAT')
     write(nlst,31) TRIM(coord_sys), &
          xorigr,yorigr,xorigr+(nx-1)*dx,yorigr+(ny-1)*dy,nx,ny,dx,dy,  &
          tpgmin,tpgmax,ztop,nz
31   format('MESH',/, &
          '  System of coord.   :  ',a               ,/,      &
          '  Bottom-left corner : (',f9.4,2x,f9.4,')',/,      &
          '  Top-right   corner : (',f9.4,2x,f9.4,')',/,      &
          '  Number points x    : ',i9  ,/,                   &
          '  Number points y    : ',i9  ,/,                   &
          '  Grid incr. (deg)   : ',f9.5,/,                   &
          '  Grid incr. (deg)   : ',f9.5,/,                   &
          '  Min. topography    : ',f9.1             ,/,      &
          '  Max. topography    : ',f9.1             ,/,      &
          '  Max. z domain      : ',f9.1             ,/,      &
          '  Number z-levels    : ',i9)
  END SELECT
  !
  do iz = 1,nz
     write(nlst,40) iz,zlayer(iz),dz(iz)
  end do
40 format('    z-layer ',i2,'       : ',f9.1,' (dz=',f6.1,')')
  !
  !***  Writes Input data to the lst file
  !
  write(nlst,100)
100 format(/,                                                           &
       '----------------------------------------------------',/,    &
       '                                                    ',/,    &
       '               FALL3D INPUT DATA                    ',/,    &
       '                                                    ',/,    &
       '----------------------------------------------------',/)
  !
  call addtime(ibyr,ibmo,ibdy,0,iyr,imo,idy,ihr,imi,ise,beg_time)
  write(nlst,250) idy,month(imo),iyr,ihr,(imi*60)+ise
250 format('FALL3D TIME RANGE',/,                      &
           '  Initial time    : ',2x,i2.2,1x,a3,1x,i4.4,' at ',i2.2,' h ',i4.4,' s')
  !
  call addtime(ibyr,ibmo,ibdy,0,iyr,imo,idy,ihr,imi,ise,end_time)
  write(nlst,251) idy,month(imo),iyr,ihr,(imi*60)+ise
251 format('  Final   time    : ',2x,i2.2,1x,a3,1x,i4.4,' at ',i2.2,' h ',i4.4,' s')

  !
  write(nlst,252) (end_time-beg_time)/36d2,(end_time-beg_time),  &
       safety_factor
252 format('  Time increment  : ',f6.1,' h (',f12.1,' s)'            ,/,    &
           '  Safety factor   : ',f6.2)
  !
  write(nlst,260)
260 format(/,                         &
       'NUMERICAL TREATMENT')
  !
  if(modkv == 0)  then
     write(nlst,261) rkv0
261  format('  Vertical turbulence      : Kv constant ( ',e12.4,' )')
  else if(modkv == 1) then
     write(nlst,262)
262  format('  Vertical turbulence      : Kv similarity theory   ')
  else if(modkv == 2) then
     write(nlst,263)
263  format('  Vertical turbulence      : Kv surface layer       ')
  end if
  !
  if(modkh == 0)  then
     write(nlst,2611) rkh0
2611 format('  Horizontal turbulence    : Kh constant ( ',e12.4,' )')
  else if(modkh == 1) then
     write(nlst,2622)
2622 format('  Horizontal turbulence    : Kh RAMS                ')
  else if(modkh == 2) then
     write(nlst,2623)
2623 format('  Horizontal turbulence    : Kh CMAQ                ')
  end if
  !
  if(lmeth == 0) message ='None           '
  if(lmeth == 1) message ='Minmod         '
  if(lmeth == 2) message ='Superbee       '
  if(lmeth == 3) message ='van Leer       '
  if(lmeth == 4) message ='Monotonized cen'
  if(lmeth == 5) message ='Beam Warming   '
  write(nlst,264) message
264 format('  Limiter method           : ',a15)
  !
  if(modv == 1) message = 'Arastoopour 82 '
  if(modv == 2) message = 'Ganser 93      '
  if(modv == 3) message = 'Wilson-Huang 79'
  if(modv == 4) message = 'Dellino 05     '
  write(nlst,266) message
266 format('  Settling velocity model  : ',a15)
  !
  write(nlst,320) np
  write(nlst,321) (TRIM(classname(ic)),ic=1,np)
  write(nlst,322) (diam(ic)*1d3,ic=1,np)
  write(nlst,323) (-log(1d3*diam(ic))/log(2.0_rp),ic=1,np)
  write(nlst,324) (rhop(ic)      ,ic=1,np)
  write(nlst,325) (sphe(ic)      ,ic=1,np)
  write(nlst,326) (psi (ic)      ,ic=1,np)
  write(nlst,327) (fc(ic)*1d2,ic=1,np)
  write(nlst,328) sum(fc(1:np))*1d2
  !
320 format(/,                                       &
       'GRANULOMETRIC DISTRIBUTION            ',/,     &
       '  NUMBER OF PARTICLE CLASSES        : ',i4)
321 format(                                         &
       '  CLASS NAME                        : ',50(a10,1x),/)
322 format(                                         &
       '  DIAMETER                  (mm   ) : ',50(f10.4,1x),/)
323 format(                                         &
       '  PHI                       (-    ) : ',50(f10.2,1x),/)
324 format(                                         &
       '  DENSITY                   (kg/m3) : ',50(f10.2,1x),/)
325 format(                                         &
       '  SPHERICITY                (-    ) : ',50(f10.2,1x),/)
326 format(                                         &
       '  MODEL FACTOR              (-    ) : ',50(f10.2,1x),/)
327 format(                                         &
       '  PERCENTAGE                (in % ) : ',50(f10.4,1x),/)
328 format(                                         &
       '  SUM                       (in % ) : ',1(f10.1,1x))
  !
  if( ng > 0 ) then
     write(nlst,329) ng
     write(nlst,3291) ( TRIM(classname(ic)) ,ic=np+1,nc)
     write(nlst,3292) ( fc(ic)*1d2  ,ic=np+1,nc)
     write(nlst,3293) sum(fc(np+1:nc))*1d2
  end if
329 format(/,                                       &
       '  NUMBER OF GAS SPECIES             : ',i4)
3291 format(                                         &
       '  CLASS NAME                        : ',50(a10,1x),/)
3292 format(                                         &
       '  PERCENTAGE (in % )                : ',50(f10.4,1x),/)
3293 format(                                         &
       '  SUM        (in % )                : ',1(f10.1,1x),/)
  !
  if(aggregation) then

     write(nlst,330) 'Costa', Df, vset_factor, ndtagr, &
          np-1, rhop(iaggr), &
          1d3*diam(MIN(iaggr,np-2)),1d3*diam(np-1),1d3*diam(np), &
          -log(1d3*diam(MIN(iaggr,np-2)))/log(2.0_rp), &
          -log(1d3*diam(          np-1 ))/log(2.0_rp),-log(1d3*diam(np))/log(2.0_rp)
330  format( &
          'AGGREGATION MODEL :  ',a,/, &
          '  FRACTAL EXPONENT : ',f7.4,/,&
          '  V.SETTLING FACTOR: ',f7.4,/,&
          '  N DT AGGREGATION : ',i7  ,/,&
          '  NUMBER OF AGGREGATING CLASSES : ',i4,/, &
          '  AVERAGED DENSITY : ',f7.2,/,&
          '  FROM DIAMETER (mm) ',f7.4,' TO ',f7.4, ' AGGREGATE TO ',f7.4,/, &
          '  FROM FI            ',f7.2,' TO ',f7.2, ' AGGREGATE TO ',f7.2,/)
  end if
  !
  write(nlst,340)
340 format( &
       'OUTPUT STRATEGY')
  !
  write(nlst,341) print_time/36d2,print_time
341 format('  Print time      : ',f6.1,' h (',f12.1,' s)',/, &
           '  Output formats  : netCDF')
  !
  if(track_points) then
     write(nlst,350) npts,(name_pts(i),xpts(i),ypts(i),use_pts(i),i=1,npts)
  else
     write(nlst,351)
  end if
350 format('  Track points    : yes'     ,/, &
           '  Num pts in file : ',i5     ,/, &
       1000(4x,a15,1x,'Coordinates: ',f13.4,1x,f13.4,'  In the domain (T/F): ',L,/) )
351 format('  Track points    : no')
  !
  !***  Writes memory requirement
  !
  write(nlst,400) memo/(1024.0_rp*1024.0_rp)
400 format(/,                                                           &
       '----------------------------------------------------',/,    &
       '                                                    ',/,    &
       '             MEMORY REQUIREMENTS                    ',/,    &
       '                                                    ',/,    &
       '----------------------------------------------------',/,    &
       '                                                    ',/,    &
       'DYNAMIC MEMORY ALLOCATED : ',f8.2,' MBy'             ,/)
  !
  !***  Writes time integration header
  !
  write(nlst,500)
500 format(/,                                                        &
       '----------------------------------------------------',/,   &
       '                                                    ',/,   &
       '               TIME INTEGRATION                     ',/,   &
       '                                                    ',/,   &
       '----------------------------------------------------',/)
  !
  return
end subroutine wridat
