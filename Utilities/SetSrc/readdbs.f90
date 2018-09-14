   subroutine readdbs
  !*********************************************************************
  !*
  !*    This routine reads meteo data from the dbs file at current time
  !*
  !*    Eruption mode:
  !*      It computes Tair(z) and Vair(z), vertical values above the vent.
  !*      It also computes Tair0 and PSIa (averaged wind direction)
  !*
  !*    Resuspension mode:
  !*      Reads UST and moisture from the dbs at the current time
  !*
  !*********************************************************************
  use KindType
  use InpOut
  use Master
  implicit none
  !
  character(len=s_mess) :: message
  integer(ip)           :: iz,istat
  real   (rp)           :: ux,uy,angle,z,Cao,dTdz
  !
  !***  Mode selection
  !
  SELECT CASE(type_source)
  case('POINT','SUZUKI','PLUME')
  !
  !***  Loop over model layers
  !
  z    = 0.0
  PSIa = 0.0
  do iz = 1,nz
     !
     !***     Gets air velocity (module) in terrain following --> Vair
     !
     call get_dbs_value_point(ludbsname,time,'U',x0,y0,zmodel(iz),ux,time2,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_dbs_value_point(ludbsname,time,'V',x0,y0,zmodel(iz),uy,time2,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     Vair(iz) = sqrt(ux*ux+uy*uy)
     !
     !***     Gets air temperature in terrain following --> Tair
     !
     call get_dbs_value_point(ludbsname,time,'T',x0,y0,zmodel(iz),Tair(iz),time2,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     !***     Gets air density in terrain following --> Rair
     !
     call get_dbs_value_point(ludbsname,time,'RHO',x0,y0,zmodel(iz),Rair(iz),time2,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     !***     Gets air direction in terrain following --> Aair
     !
     if(abs(ux).gt.1d-8) then
        angle = atan2(uy,ux)*180.0/pi
     else
        if(uy.gt.1d-8) then
           angle = 90.0
        else if(uy.lt.-1d-8) then
           angle = 270.0
        else
           angle = 0.0
        end if
     end if
     !
     Aair(iz) = angle
     if(Aair(iz).lt.0.0) Aair(iz) = 360.0 + Aair(iz)    ! Angle in deg. (0 to 360)
     !
     PSIa = PSIa + angle*(zmodel(iz)-z)
     z    = zmodel(iz)
     !
  end do
  !
  !***  Averaged wind direction
  !
  PSIa = PSIa/(zmodel(nz))                ! Angle in deg. (-180 to +180)
  if(PSIa.lt.0.0) PSIa = 360.0 + PSIa   ! Angle in deg. (0 to 360)
  PSIa = PSIa*pi/180.0                   ! Angle in rad.
  !
  !***  Computes Tair0 and Rair0
  !
  Tair0 = Tair(1)
  Rair0 = Rair(1)
  !
  !*** Computes the buoyancy frequency (squared)
  !
  Cao  = 998.    ! specific heat capacity at constant pressure of dry air (J kg^-1 K^-1)
  do iz = 1,nz
     !
     if(iz.eq.1) then
        dTdz = (Tair(2)-Tair(1))/(zmodel(2)-zmodel(1))
     else if (iz.eq.nz) then
        dTdz = (Tair(nz)-Tair(nz-1))/(zmodel(nz)-zmodel(nz-1))
     else
        dTdz = (Tair(iz+1)-Tair(iz-1))/(zmodel(iz+1)-zmodel(iz-1))
     end if
     !
     Nair(iz) = g*g*(1+Cao*dTdz/g)/(Cao*Tair0)
     !
  end do
  !
  !
  !
  case('RESUSPENSION')
  !
     if (recalculate_ust) then
       call get_dbs_value_plane(ludbsname,time,'SPD10',0.0d0,nx,ny,ust,time2,istat,message)
       ust = ust * 0.4 / LOG(10.0_rp/roughness_ust)
     else
       call get_dbs_value_plane(ludbsname,time,'UST' ,0.0d0,nx,ny,ust,time2,istat,message)
     end if
     call get_dbs_value_plane(ludbsname,time,'SMOI',0.0d0,nx,ny,smoi ,time2,istat,message)
     call get_dbs_value_plane(ludbsname,time,'SOIL',0.0d0,nx,ny,soil,time2,istat,message)
     if(istat.lt.0) call runend(message)
     call get_dbs_value_plane(ludbsname,time,'VEGFRA',0.0d0,nx,ny,vegfra,time2,istat,message)
     if(istat.lt.0) call runend(message)
     call get_dbs_value_plane(ludbsname,time,'RMOL',0.0d0,nx,ny,rmol,time2,istat,message)
     if(istat.lt.0) call runend(message)
     call get_dbs_value_plane(ludbsname,time,'RHO' ,0.0d0,nx,ny,densi,time2,istat,message)
  !
  END SELECT
  !
  !***  Bounds the value of time2
  !
  time2 = MIN(time2,ieend)
  !
  return
end subroutine readdbs
