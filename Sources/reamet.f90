   subroutine reamet
  !****************************************************************
  !*
  !*    Reads data from the meteo dbs file in netCDF format
  !*
  !****************************************************************
  use KindType
  use InpOut
  use Master
  use Numeric
  use Parallel
  implicit none
  !
  character(len=s_mess) :: message
  integer(ip)           :: istat,k,iout
  real   (rp)           :: z
  !
  !***  Gets 3D data: wind velocity field and temperature
  !
  do k = 1,nz
     z = zlayer(k)
     !
     if( mpime == root ) then
        call get_dbs_value_plane(fdbs,INT(time_lag+time),'U',z,nx,ny,vx(1,1,k),iout,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     !
     if( mpime == root ) then
        call get_dbs_value_plane(fdbs,INT(time_lag+time),'V',z,nx,ny,vy(1,1,k),iout,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     !
     if( mpime == root ) then
        call get_dbs_value_plane(fdbs,INT(time_lag+time),'W',z,nx,ny,vz(1,1,k),iout,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     !
     if( mpime == root ) then
        call get_dbs_value_plane(fdbs,INT(time_lag+time),'T',z,nx,ny,tempe(1,1,k),iout,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     !
     if( mpime == root ) then
        call get_dbs_value_plane(fdbs,INT(time_lag+time),'TP',z,nx,ny,vptemp(1,1,k),iout,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     !
     if( mpime == root ) then
        call get_dbs_value_plane(fdbs,INT(time_lag+time),'P',z,nx,ny,p(1,1,k),iout,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     !
     if( mpime == root ) then
        call get_dbs_value_plane(fdbs,INT(time_lag+time),'RHO',z,nx,ny,rho(1,1,k),iout,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     !
     if( mpime == root ) then
        call get_dbs_value_plane(fdbs,INT(time_lag+time),'QV',z,nx,ny,qv(1,1,k),iout,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     !
     call bcast( iout, 1, root )
     meteo_time = DBLE(iout)-time_lag
  end do
  !
  call bcast( vx    , SIZE(vx    ), root )
  call bcast( vy    , SIZE(vy    ), root )
  call bcast( vz    , SIZE(vz    ), root )
  call bcast( tempe , SIZE(tempe ), root )
  call bcast( vptemp, SIZE(vptemp), root )
  call bcast( p     , SIZE(p     ), root )
  call bcast( rho   , SIZE(rho   ), root )
  call bcast( qv    , SIZE(qv    ), root )
  !
  !***  Gets 2D data at surface
  !
  if( mpime == root ) then
     call get_dbs_value_plane(fdbs,INT(time_lag+time),'UST',0.0_rp,nx,ny,ustar(1,1),iout,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  !
  if( mpime == root ) then
     call get_dbs_value_plane(fdbs,INT(time_lag+time),'PBLH',0.0_rp,nx,ny,zi(1,1),iout,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) then
     call wriwar('reamet : no value for PBLH. Assuming PBLH = 1.5km')
     zi = 1500.0_rp
  end if
  !
  if( mpime == root ) then
     call get_dbs_value_plane(fdbs,INT(time_lag+time),'RMOL',0.0_rp,nx,ny,rl(1,1),iout,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  !
  if( mpime == root ) then
     call get_dbs_value_plane(fdbs,INT(time_lag+time),'ZNT',0.0_rp,nx,ny,znt(1,1),iout,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) then
     call wriwar('reamet : no value for ZNT. Assuming ZNT = 0.1m')
     znt = 0.1_rp
  end if
  !
  call bcast( ustar, SIZE(ustar), root )
  call bcast( zi   , SIZE(zi   ), root )
  call bcast( rl   , SIZE(rl   ), root )
  call bcast( znt  , SIZE(znt  ), root )
  !
  !*** Wet deposition
  !
  if(wetdeposition) then
     if( mpime == root ) then
        call get_dbs_value_plane(fdbs,INT(time_lag+time),'PRATE',0.0_rp,nx,ny,prate(1,1),iout,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) then
        call wriwar('reamet : no value for PRATE. Assuming PRATE = 0.0')
        prate = 0.0_rp
     end if
     call bcast( prate, SIZE(prate), root )
     prate(1:nx,1:ny) = prate(1:nx,1:ny)*3600.0_rp  ! mms-1 --> mmh-1
  end if
  !
  !*** Other computations
  !
  !vx(1:nx,1:ny,1) = 0.  ! Set vx(z=0)=0
  !vy(1:nx,1:ny,1) = 0.  ! Set vy(z=0)=0
  vz(1:nx,1:ny,1) = 0.  ! Set vz(z=0)=0
  !
  return
  end subroutine reamet
