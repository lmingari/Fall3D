module wind_rotation
  !
  ! This module performs wind rotations in 3D
  !
  ! Note: angles are expressed in degrees
  !
  ! Author: G.Macedonio
  !
  ! Version: 1.0
  ! Date: 22-MAY-2013
  !
  use kindtype
  implicit none
  !
  logical  :: rotate_wind=.false.               ! Default is false
  real(rp) :: rotation_lonpole,rotation_latpole,rotation_angle
  real(rp), allocatable :: unew(:,:),vnew(:,:)  ! Temporary storage
  !
  ! Degrees <-> radians
  real(4) :: to_rad_r4 = .0174533            ! PI/180
  real(8) :: to_rad_r8 = .017453292519943d0  ! PI/180
  real(4) :: to_deg_r4 = 57.29578            ! 180/PI
  real(8) :: to_deg_r8 = 57.29577951308d0    ! 180/PI

  interface llconv
     module procedure llconv_real4
     module procedure llconv_real8
  end interface

  interface rllconv
     module procedure rllconv_real4
     module procedure rllconv_real8
  end interface

  interface rotate2d
     module procedure rotate2d_real4
     module procedure rotate2d_real8
  end interface

  interface rotate3d
     module procedure rotate3d_real4
     module procedure rotate3d_real8
  end interface

contains

  subroutine llconv_real4(lonpole,latpole,lon,lat,newlon,newlat)
    !
    ! Trasforma le coordinate geografiche (lon,lat) di un punto, rispetto ad
    ! un sistema di riferimento geografico, in nuove coordinate (newlon,newlat)
    ! definite rispetto ad un riferimento che ha il polo nella posizione
    ! (lonpole,latpole).
    ! Lo zero della longitudine e' il meridiano che passa per il nuovo polo.
    !
    ! Le coordinate sono definite in gradi
    ! Le locazioni di memoria delle variabili di output (newlon, newlat)
    ! possono coincidere con quelle di input (lon, lat)
    !
    ! Note: Function asin returns a value in the range  [-PI, PI]
    ! Note: Function acos returns a value in the range  [  0, 2*PI]
    ! Note: Function atan2 returns a value in the range [-PI, PI]
    !
    implicit none
    !
    real(4), intent(in) :: lonpole,latpole
    real(4), intent(in) :: lon,lat
    real(4), intent(out):: newlon,newlat
    real(4) :: mylon,sinmylon,cosmylon,sinlat,coslat,sinlatpole,coslatpole
    !
    ! Prende come riferimento il meridiano passante per il nuovo polo
    !
    mylon = lon-lonpole
    !
    sinmylon = sin(to_rad_r4*mylon)
    cosmylon = cos(to_rad_r4*mylon)
    sinlat = sin(to_rad_r4*lat)
    coslat = cos(to_rad_r4*lat)
    sinlatpole = sin(to_rad_r4*latpole)
    coslatpole = cos(to_rad_r4*latpole)
    !
    ! Apply the Cosine theorem to find the new latitude
    newlat = to_deg_r4*asin(sinlat*sinlatpole + coslat*coslatpole*cosmylon)
    !
    ! Apply cotangent theorem to find the new longitude
    newlon = to_deg_r4*atan2(sinmylon*coslat,coslatpole*sinlat - &
         sinlatpole*coslat*cosmylon)  ! This is 180-newlon
    !
    ! Transform to newlon in the range [-180,180]
    if(newlon < 0.0) then
       newlon = -180.0 - newlon
    else
       newlon = 180.0 - newlon
    end if
    !
  end subroutine llconv_real4

  subroutine llconv_real8(lonpole,latpole,lon,lat,newlon,newlat)
    !
    ! Trasforma le coordinate geografiche (lon,lat) di un punto, rispetto ad
    ! un sistema di riferimento geografico, in nuove coordinate (newlon,newlat)
    ! definite rispetto ad un riferimento che ha il polo nella posizione
    ! (lonpole,latpole). Lo zero della latitudine e' il meridiano che
    ! passa per il nuovo polo.
    !
    ! Le coordinate sono definite in gradi
    !
    ! Le locazioni di memoria delle variabili di output (newlon, newlat)
    ! possono coincidere con quelle di input (lon, lat)
    !
    implicit none
    !
    real(8), intent(in) :: lonpole,latpole
    real(8), intent(in) :: lon,lat
    real(8), intent(out):: newlon,newlat
    real(8) :: mylon,sinmylon,cosmylon,sinlat,coslat,sinlatpole,coslatpole
    !
    ! Prende come riferimento il meridiano passante per il nuovo polo
    mylon = lon-lonpole
    !
    sinmylon = sin(to_rad_r8*mylon)
    cosmylon = cos(to_rad_r8*mylon)
    sinlat = sin(to_rad_r8*lat)
    coslat = cos(to_rad_r8*lat)
    sinlatpole = sin(to_rad_r8*latpole)
    coslatpole = cos(to_rad_r8*latpole)
    !
    ! Apply the Cosine theorem to find the new latitude
    newlat = to_deg_r8*asin(sinlat*sinlatpole + coslat*coslatpole*cosmylon)
    !
    ! Apply cotangent theorem to find the new longitude
    newlon = to_deg_r8*atan2(sinmylon*coslat,coslatpole*sinlat - &
         sinlatpole*coslat*cosmylon)  ! This is 180-newlon
    !
    ! Transform to newlon in the range [-180,180]
    if(newlon < 0d0) then
       newlon = -180d0 - newlon
    else
       newlon = 180d0 - newlon
    end if

  end subroutine llconv_real8

  subroutine rllconv_real4(lonpole,latpole,lon,lat,newlon,newlat)
    ! Trasformazione inversa di llconv
    !
    ! Le locazioni di memoria delle variabili di output (lon, lat)
    ! possono coincidere con quelle di input (newlon, newlat)
    !
    implicit none
    !
    real(4), intent(in)  :: lonpole,latpole
    real(4), intent(inout) :: lon,lat
    real(4), intent(in)  :: newlon,newlat
    real(4) :: csi,sincsi,coscsi,sinnewlat,cosnewlat,sinlatpole,coslatpole
    !
    csi = 180.0 - newlon
    sincsi = sin(to_rad_r4*csi)
    coscsi = cos(to_rad_r4*csi)
    sinnewlat = sin(to_rad_r4*newlat)
    cosnewlat = cos(to_rad_r4*newlat)
    sinlatpole = sin(to_rad_r4*latpole)
    coslatpole = cos(to_rad_r4*latpole)
    !
    ! Apply cosin theorem to find the latitude
    lat = to_deg_r4*asin(sinnewlat*sinlatpole + cosnewlat*coslatpole*coscsi)

    ! Apply cotangent theorem to find the longitude
    lon = to_deg_r4*atan2(sincsi*cosnewlat,sinnewlat*coslatpole - &
         cosnewlat*sinlatpole*coscsi)
    !
    ! Reset the origin
    lon = lon + lonpole
    !
    if(lon > 180.0) then
       lon = lon - 360.
    elseif(lon < -180.0) then
       lon = lon + 360.
    end if

  end subroutine rllconv_real4

  subroutine rllconv_real8(lonpole,latpole,lon,lat,newlon,newlat)
    ! Trasformazione inversa di llconv
    !
    ! Le locazioni di memoria delle variabili di output (lon, lat)
    ! possono coincidere con quelle di input (newlon, newlat)
    !
    implicit none
    !
    real(8), intent(in)  :: lonpole,latpole
    real(8), intent(inout) :: lon,lat
    real(8), intent(in)  :: newlon,newlat
    real(8) :: csi,sincsi,coscsi,sinnewlat,cosnewlat,sinlatpole,coslatpole
    !
    csi = 180d0 - newlon
    sincsi = sin(to_rad_r8*csi)
    coscsi = cos(to_rad_r8*csi)
    sinnewlat = sin(to_rad_r8*newlat)
    cosnewlat = cos(to_rad_r8*newlat)
    sinlatpole = sin(to_rad_r8*latpole)
    coslatpole = cos(to_rad_r8*latpole)
    !
    ! Applica il Teorema del Coseno
    lat = to_deg_r8*asin(sinnewlat*sinlatpole + cosnewlat*coslatpole*coscsi)

    ! Apply cotangent theorem to find the longitude
    lon = to_deg_r8*atan2(sincsi*cosnewlat,sinnewlat*coslatpole - &
         cosnewlat*sinlatpole*coscsi)
    !
    ! Ripristina origine
    lon = lon + lonpole
    !
    if(lon > 180d0) then
       lon = lon - 360d0
    elseif(lon < -180d0) then
       lon = lon + 360d0
    end if

  end subroutine rllconv_real8

  subroutine rotate2d_real4(xold,yold,xnew,ynew,angle)
    !
    ! Ruota un punto nel piano
    !
    ! Le locazioni di memoria delle variabili nuove (xnew,ynew) possono
    ! coincidere con quelle precedenti (xold,yold) 
    !
    ! L'angolo viene espresso in gradi
    !
    implicit none
    real(4), intent(in)  :: xold,yold
    real(4), intent(in)  :: angle
    real(4), intent(out) :: xnew,ynew
    real(4) :: xsave
    !
    xsave = xold            ! Copia coordinata xold
    xnew = xold*cos(to_rad_r4*angle) + yold*sin(to_rad_r4*angle)
    ynew = xsave*sin(to_rad_r4*angle) + yold*cos(to_rad_r4*angle)
    !
  end subroutine rotate2d_real4

  subroutine rotate2d_real8(xold,yold,xnew,ynew,angle)
    !
    ! Ruota un punto nel piano
    !
    !
    ! Le locazioni di memoria delle variabili nuove (xnew,ynew) possono
    ! coincidere con quelle precedenti (xold,yold) 
    !
    ! L'angolo viene espresso in gradi
    !
    implicit none
    real(8), intent(in)  :: xold,yold
    real(8), intent(in)  :: angle
    real(8), intent(out) :: xnew,ynew
    real(8) :: xsave
    !
    xsave = xold            ! Copia coordinata xold
    xnew = xold*cos(to_rad_r8*angle) + yold*sin(to_rad_r8*angle)
    ynew = xsave*sin(to_rad_r8*angle) + yold*cos(to_rad_r8*angle)
    !
  end subroutine rotate2d_real8

  subroutine rotate3d_real4(lonpole,latpole,angle,oldlon,oldlat,newlon,newlat)
    !
    ! Ruota un punto (xold,xnew) attorno a un asse passante per un polo
    !
    ! Le variabili nuove (xnew,ynew) possono eventualmente coincidere
    ! con quelle precedenti (xold,yold) 
    !
    ! L'angolo di rotazione viene espresso in gradi
    !
    implicit none
    real(4), intent(in)  :: lonpole,latpole
    real(4), intent(in)  :: oldlon,oldlat
    real(4), intent(in)  :: angle
    real(4), intent(out) :: newlon,newlat
    !
    ! Phase 1: Transform to coordinates respect new system
    call llconv(lonpole,latpole,oldlon,oldlat,newlon,newlat)
    !
    ! Phase 2: Rotate around the new pole
    newlon = newlon + angle
    !
    ! Fase 3: Transform back to original system
    call rllconv(lonpole,latpole,newlon,newlat,newlon,newlat)
    !
  end subroutine rotate3d_real4

  subroutine rotate3d_real8(lonpole,latpole,angle,oldlon,oldlat,newlon,newlat)
    !
    ! Ruota un punto (xold,xnew) attorno a un asse passante per un polo
    !
    ! Le variabili nuove (xnew,ynew) possono eventualmente coincidere
    ! con quelle precedenti (xold,yold) 
    !
    ! L'angolo di rotazione viene espresso in gradi
    !
    implicit none
    real(8), intent(in)  :: lonpole,latpole
    real(8), intent(in)  :: oldlon,oldlat
    real(8), intent(in)  :: angle
    real(8), intent(out) :: newlon,newlat
    !
    ! Phase 1: Transform to coordinates respect new system
    call llconv(lonpole,latpole,oldlon,oldlat,newlon,newlat)
    !
    ! Phase 2: Rotate around the new pole
    newlon = newlon + angle
    !
    ! Fase 3: Transform back to original system
    call rllconv(lonpole,latpole,newlon,newlat,newlon,newlat)
    !
  end subroutine rotate3d_real8

end module wind_rotation
