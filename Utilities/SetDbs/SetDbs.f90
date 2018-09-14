program SetDbs
  !*****************************************************************************
  !*
  !*    AUTHOR  : A.Folch
  !*    VERSION : FALL3D-7.0
  !*    DATE    : DEC-2012
  !*    PURPOSE : This program generates the DBS files needed by Fall3d-7.0
  !*
  !*****************************************************************************
  use InpOut
  use Master
  implicit none
  !
  !***  Gets filenames from the call arguments
  !
  iarg = 1                          ! log   file
  call GETARG(iarg,lulogname)
  iarg = 2                          ! inp   file
  call GETARG(iarg,luinpname)
  iarg = 3                          ! meteo netCDF file, CALMET or profile
  call GETARG(iarg,ludatname)
  iarg = 4                          ! dbs output file (netCDF format)
  call GETARG(iarg,ludbsname)
  iarg = 5                          ! topography file
  call GETARG(iarg,lutopname)
  iarg = 6
  call GETARG(iarg,TYPE_DATA)
  !
  SELECT CASE (TYPE_DATA)
  CASE('debug')
     continue
  CASE('profile')
     continue
  CASE('calmet62')
     continue
  CASE('ncep1')
     continue
  CASE('wrf')
     continue
  CASE('arpa')
     continue
  CASE('gfs')
     continue
  CASE('ncep2')
     continue
  CASE('eta')
     continue
  CASE('ecmwf')
     continue
  CASE DEFAULT
     TYPE_DATA = 'wrf'  ! default value
  END SELECT
  !
  !***  Opens the log file
  !
  call openinp
  !
  !***  Gets database properties from input file
  !
  call readinp
  !
  !***  First time. Reads the grid and calculates the interpolation
  !***  factors
  !
  SELECT CASE (TYPE_DATA)
  CASE('debug')
     continue
  CASE('calmet62')
     call read_CAL_grid (ludatname,nx*ny,xutm,yutm)
     call checktime_CAL
  CASE('ncep1')
     call read_NCEP1_grid (ludatname,nx*ny,lonp,latp)
     call checktime_NCEP1
  CASE('wrf')
     call read_WRF_grid (ludatname,nx*ny,lonp,latp)
     call checktime_WRF
  CASE('arpa')
     call read_ARPA_grid (ludatname,nx*ny,lonp,latp)
     call checktime_ARPA
  CASE('gfs')
     call read_GFS_grid (ludatname,nx*ny,lonp,latp)
     call checktime_GFS
  CASE('ncep2')
     call read_NCEP2_grid (ludatname,nx*ny,lonp,latp)
     call checktime_NCEP2
  CASE('eta')
     call read_ETA_grid (ludatname,nx*ny,lonp,latp)
     call checktime_ETA
  CASE('ecmwf')
     call read_ECMWF_grid (ludatname,nx*ny,lonp,latp)
     call checktime_ECMWF
  END SELECT
  !
  call write_DBS_grid(ludbsname,nx,ny,nz,nt,lonp,latp,zlayer,time,timesec,xutm,yutm,zutm,coord_sys)
  !
  !***  Loop over time steps
  !
  do it = 1,nt
     !
     SELECT CASE (TYPE_DATA)
     CASE('debug')
        call read_DEBUG_data(ludatname,time_lag,time(it),timesec(it), &
             nx,ny,nz,zlayer,lonp,latp,topg,LDU,u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
     CASE('profile')
        call read_PROF_data(ludatname,time_lag,time(it),timesec(it), &
             nx,ny,nz,zlayer,lonp,latp,topg,LDU,u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
     CASE('calmet62')
        call read_CAL_data(ludatname,time_lag,time(it),timesec(it), &
             nx,ny,nz,zlayer,lonp,latp,topg,LDU,u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
     CASE('ncep1')
        call read_NCEP1_data(ludatname,time_lag,time(it),timesec(it), &
             nx,ny,nz,zlayer,lonp,latp,topg,LDU,u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
     CASE('wrf')
        call read_WRF_data(ludatname,time_lag,time(it),timesec(it), &
             nx,ny,nz,zlayer,lonp,latp,topg,LDU,SOIL,vegfra,u,v,w,T,Tv,Tp, &
             qv,ro,p,pblh,ust,rmol,znt,spd10,L,smoi,prat)
     CASE('arpa')
        call read_ARPA_data(ludatname,time_lag,time(it),timesec(it), &
             nx,ny,nz,zlayer,lonp,latp,topg,LDU,u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
     CASE('gfs')
        call read_GFS_data(ludatname,time_lag,time(it),timesec(it), &
             nx,ny,nz,zlayer,lonp,latp,topg,LDU,u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
     CASE('ncep2')
        call read_NCEP2_data(ludatname,time_lag,time(it),timesec(it), &
             nx,ny,nz,zlayer,lonp,latp,topg,LDU,u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
     CASE('eta')
        call read_ETA_data(ludatname,time_lag,time(it),timesec(it), &
             nx,ny,nz,zlayer,lonp,latp,topg,LDU,u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
     CASE('ecmwf')
        call read_ECMWF_data(ludatname,time_lag,time(it),timesec(it), &
             nx,ny,nz,zlayer,lonp,latp,topg,LDU,u,v,w,T,Tv,Tp,qv,ro,p,pblh,ust,L)
     END SELECT
     !
     call write_DBS_data(ludbsname,nx,ny,nz,it,topg,LDU,SOIL,vegfra,u,v,w,T,Tp,p, &
                         qv,ro,pblh,ust,rmol,znt,spd10,L,smoi,prat)
     !
  end do
  !
  call runend('OK')
  !
end program SetDbs
