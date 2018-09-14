!***************************************************************
!*
!*    Module for coordinate operations
!*
!***************************************************************
MODULE Coordinates
  use KindType, ONLY :  ip,rp
  implicit none
  !
  integer(ip), parameter :: GRID_ZONE_LENGTH = 3
  integer(ip), parameter :: CLARKE_1866_DATUM = 0
  integer(ip), parameter :: GRS_80_DATUM = 1
  integer(ip), parameter :: WGS_84_DATUM = 2
  integer(ip), parameter :: NAD83_DATUM = 3
  integer(ip), parameter :: AIRY_1830_DATUM = 4
  !
  real(rp), parameter :: LOWER_EPS_LIMIT = 1d-14
  real(rp), parameter :: M_PI = 3.14159265358979323846     ! pi
  real(rp), parameter :: M_PI_2 = 1.57079632679489661923   ! pi/2
  !
CONTAINS
  !
  !
  !
  subroutine ll2utm (rlat,rlon, grid_zone, utmx, utmy, datum, istat)
    !***********************************************************************
    !*
    !*    Converts lat/long to UTM, using the specified datum
    !*
    !*    Input  : rlat,rlon,datum
    !*    Output : grid_zone,utmx,utmy,istat
    !*
    !***********************************************************************
    implicit none
    integer(ip) :: istat,datum
    real   (rp) :: rlat,rlon, utmx, utmy
    character(len=3) :: grid_zone
    !
    character(len=3) :: utmgz
    integer(ip) :: i
    real(rp) :: a,b,f,e,e2,e4,e6,phi,phi0,t,x,y, rho, cc,cc2,dlamb0
    real(rp) :: aa,aa2,aa3,aa4,aa5,aa6,lambda, k0, mm, mm0, nn, ep2,tt2,tt
    !
    if(datum.eq.CLARKE_1866_DATUM) then
       a = 6378206.4        ! semimajor axis of ellipsoid (meters)
       b = 6356583.8        ! semiminor axis of ellipsoid (meters)
    elseif(datum.eq.GRS_80_DATUM) then
       a = 6378137.0
       b = 6356752.3
    elseif(datum.eq.WGS_84_DATUM) then
       a = 6378137.0
       b = 6356752.31425 ! Used by GoogleEarth
    elseif(datum.eq.NAD83_DATUM) then
       a = 6378137.0
       b = 6356752.3141
    elseif(datum.eq.AIRY_1830_DATUM) then
       a = 6377563.4
       b = 6356256.9
    else
       istat = -1
       return
    endif
    !
    !***  Calculate flatness and eccentricity
    !
    f = 1.0 - (b/a)
    e2 = (2.0 - f) * f
    e = sqrt(e2)
    e4 = e2*e2
    e6 = e4*e2
    !
    !***  Convert latitude/longitude to radians
    !
    phi = rlat * M_PI / 180.0
    lambda = rlon * M_PI / 180.0
    !
    !***  Figure out the UTM zone, as well as dlamb0
    !
    call get_grid_zone (rlat, rlon, utmgz, dlamb0)
    phi0 = 0.0
    !
    !***  See if this will use UTM or UPS
    !
    if (rlat.gt.84.0) then
       !     use Universal Polar Stereographic Projection (north polar aspect)
       k0 = 0.994

       t = sqrt(((1.0-sin(phi))/(1.0+sin(phi))) *  &
            ((1.0+e*sin(phi))/(1.0-e*sin(phi)))**e)
       rho = 2.0*a*k0*t/sqrt((1.0+e)**(1.0+e)*(1.0-e)**(1.0-e))
       x = rho * sin(lambda - dlamb0)
       y = -rho * cos(lambda - dlamb0)
       !
       !***     Apply false easting/northing
       !
       x = x + 2e6
       y = y + 2e6

    elseif (rlat.lt.(-80.0)) then

       !     use Universal Polar Stereographic Projection (south polar aspect)

       phi = -phi
       lambda = -lambda
       dlamb0 = -dlamb0

       k0 = 0.994
       t = sqrt(((1.0-sin(phi))/(1.0+sin(phi))) *  &
            ((1.0+e*sin(phi))/(1.0-e*sin(phi)))**e)
       rho =2.0*a*k0*t/ sqrt((1.0+e)**(1.0+e)*(1.0-e)**(1.0-e))

       x = rho * sin(lambda - dlamb0)
       y = -rho * cos(lambda - dlamb0)
       x = -x
       y = -y
       !
       !***    Apply false easting/northing
       !
       x = x + 2e6
       y = y + 2e6

    else

       !     Use UTM
       !     set scale on central median (0.9996 for UTM)

       k0 = 0.9996
       mm = a * ((1.0-e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0) * phi -  &
            (3.0*e2/8.0+3.0*e4/32.0+45.0*e6/1024.0)*sin(2*phi)+    &
            (15.0*e4/256.0 + 45.0*e6/1024.0) * sin(4.0*phi) -   &
            (35.0*e6/3072.0) * sin(6.0*phi))

       mm0 = a * ((1.0-e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0)*phi0 -  &
            (3.0*e2/8.0+3.0*e4/32.0+45.0*e6/1024.0)*sin(2.0*phi0) +  &
            (15.0*e4/256.0 + 45.0*e6/1024.0) * sin(4.0*phi0) -        &
            (35.0*e6/3072.0) * sin(6.0*phi0))

       aa = (lambda - dlamb0) * cos(phi)
       aa2 = aa*aa
       aa3 = aa2*aa
       aa4 = aa2*aa2
       aa5 = aa4*aa
       aa6 = aa3*aa3

       ep2 = e2 / (1.0 - e2)
       nn = a / sqrt(1.0 - e2*sin(phi)**2)
       tt = tan(phi)**2
       tt2 = tt*tt
       cc = ep2 * cos(phi)**2
       cc2 = cc*cc

       x = k0*nn*(aa+(1.0-tt+cc)*aa3/6.0+(5.0-18.0*tt+tt2+ &
            72.0*cc-58.0*ep2)*aa5/120.0)
       y = k0*(mm-mm0+nn*tan(phi)*(aa2/2.0+(5.0-tt+9.0*cc+ &
            4.0*cc2)*aa4/24.0+(61.0-58.0*tt+tt2+600.0*cc-  &
            330.0*ep2)*aa6/720.0))
       !
       !***     Apply false easting and northing
       !
       x = x + 5e5
       if(y.lt.0.0) y = y + 1d7
    endif
    !
    !***  Set entries in UTM structure
    !
    do i=1,GRID_ZONE_LENGTH
       grid_zone(i:i) = utmgz(i:i)
    enddo
    utmx = x
    utmy = y
    !
    !***  done
    !
    istat = 0
    return
  end subroutine ll2utm
  !
  !
  !
  subroutine utm2ll(utmx,utmy,utmgz,rlon,rlat,datum,istat)
    !***********************************************************************
    !*
    !*    Converts UTM to lat/long, using the specified datum
    !*
    !*    Input  : utmgz,utmx,utmy,datum
    !*    Output : rlat,rlon,istat
    !*
    !************************************************************************
    implicit none
    integer(ip) :: istat,datum
    real(rp) :: utmx,utmy,rlon,rlat
    character(len=3) :: utmgz
    !
    integer(ip) :: lat_zone
    real(rp) :: a,b,f,e,e1,e2,e4,e6,e8,dlamb0,x,y,rho,t,chi,phit,phi,phi0
    real(rp) :: lambda, k0, mm, mm0, mu, nn1,e12,e13,e14,ep2,tt1
    real(rp) :: dd,dd2,dd3,dd4,dd5,dd6,cc1,rr1,phi1
    !
    if(datum.eq.CLARKE_1866_DATUM) then
       a = 6378206.4        ! semimajor axis of ellipsoid (meters)
       b = 6356583.8        ! semiminor axis of ellipsoid (meters)
    elseif(datum.eq.GRS_80_DATUM) then
       a = 6378137.0
       b = 6356752.3
    elseif(datum.eq.WGS_84_DATUM) then
       a = 6378137.0
       b = 6356752.31425
    else
       istat = -1
       return
    endif
    !
    !***  Calculate flatness and eccentricity
    !
    f = 1.0 - (b/a)
    e2 = (2.0 - f) * f
    e = sqrt(e2)
    e4 = e2*e2
    e6 = e4*e2
    e8 = e4*e4
    !
    !***  Given the UTM grid zone, generate a baseline dlamb0
    !
    call get_lambda0(utmgz, dlamb0,istat)
    if(istat.lt.0) then
       istat = -1
       return
    endif
    !
    lat_zone = ichar(utmgz(3:3))
    !
    !***  Take care of the polar regions first.
    !
    if(lat_zone.eq.ichar('Y').or.lat_zone.eq.ichar('Z')) then
       !     north polar aspect
       !     Subtract the false easting/northing
       x = utmx - 2e6
       y = utmy - 2e6
       !
       !***  Solve for inverse equations
       !
       k0 = 0.994
       rho = sqrt(x*x + y*y)
       t = rho*sqrt((1.0+e)**(1.0+e)*(1.0-e)**(1.0-e))/(2.0*a*k0)
       !
       !***     Solve for latitude and longitude
       !
       chi = M_PI_2 - 2.0 * atan(t)
       phit = chi + (e2/2.0+5.0*e4/24.0+e6/12.0+13.0*e8/360.0)* &
            sin(2.0*chi)+(7.0*e4/48.0+29.0*e6/240.0+811.0*e8/11520.0)* &
            sin(4.0*chi)+(7.0*e6/120.0+81.0*e8/1120.0)*sin(6.0*chi)+ &
            (4279.0*e8/161280.0)*sin(8.0*chi)
       !
10     continue
       phi = phit
       phit = M_PI_2 - 2.0*atan(t*((1.0-e*sin(phi))/(1.0+e*sin(phi)))**(e/2.0))

       if(abs(phi-phit).le.LOWER_EPS_LIMIT) goto 20
       goto 10
20     continue

       lambda = dlamb0 + atan2(x,-y)
    elseif(lat_zone.eq.ichar('A').or.lat_zone.eq.ichar('B')) then
       !       south polar aspect
       !       Subtract the false easting/northing

       x = -(utmx - 2e6)
       y = -(utmy - 2e6)
       !
       !***     Solve for inverse equations
       !
       k0 = 0.994
       rho = sqrt (x*x + y*y)
       t = rho*sqrt((1.0+e)**(1.0+e)*(1.0-e)**(1.0-e))/(2.0*a*k0)
       !
       !***    Solve for latitude and longitude
       !
       chi = M_PI_2 - 2.0 * atan(t)
       phit = chi + (e2/2.0 + 5.0*e4/24.0 + e6/12.0 + 13.0*e8/360.0)* &
            sin(2.0*chi) + (7.0*e4/48.0 + 29.0*e6/240.0 + &
            811.0*e8/11520.0) * sin(4.0*chi) + (7.0*e6/120.0 + &
            81.0*e8/1120.0) * sin(6.0*chi) + (4279.0*e8/161280.0)* &
            sin(8.0*chi)
       !
11     continue
       phi = phit
       phit = M_PI_2 - 2.0*atan(t*((1.0-e*sin(phi))/(1.0+e*sin(phi)))**(e/2.0))
       !
       if(abs(phi-phit) .le. LOWER_EPS_LIMIT) goto 21
       goto 11
21     continue
       phi = -phi
       lambda = -(-dlamb0 + atan2(x,-y))
    else
       !
       !     Now take care of the UTM locations
       !
       k0 = 0.99960
       !
       !***   Remove false eastings/northings
       !
       x = utmx - 5e5
       y = utmy
       !
       if(lat_zone.gt.ichar('B').and.lat_zone.lt.ichar('N')) &
            y = y - 1d7       ! southern hemisphere
       !
       !     Calculate the footpoint latitude
       !
       phi0 = 0.0
       e1 = (1.0 - sqrt(1.0-e2))/(1.0 + sqrt(1.0-e2))
       e12 = e1*e1
       e13 = e1*e12
       e14 = e12*e12
       !
       mm0 = a * ((1.0-e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0)*phi0- &
            (3.0*e2/8.0+3.0*e4/32.0+45.0*e6/1024.0)*sin(2.0*phi0)+ &
            (15.0*e4/256.0 + 45.0*e6/1024.0)*sin(4.0*phi0)- &
            (35.0*e6/3072.0) * sin(6.0*phi0))
       mm = mm0 + y/k0
       mu = mm/(a*(1.0-e2/4.0-3*e4/64.0-5*e6/256.0))
       !
       phi1 = mu + (3.0*e1/2.0 - 27.0*e13/32.0) * sin(2.0*mu)+ &
            (21.0*e12/16.0 - 55.0*e14/32.0) * sin(4.0*mu)+ &
            (151.0*e13/96.0) * sin(6.0*mu)+ &
            (1097.0*e14/512.0) * sin(8.0*mu)
       !
       !***    Now calculate lambda and phi
       !
       ep2 = e2/(1.0-e2)
       cc1 = ep2*cos(phi1)**2
       tt1 = tan(phi1)**2
       nn1 = a / sqrt(1.0-e2*sin(phi1)**2)
       rr1 = a * (1.0-e2)/(1.0-e2*sin(phi1)**2)**(1.5)
       dd = x / (nn1 * k0)
       !
       dd2 = dd*dd
       dd3 = dd*dd2
       dd4 = dd2*dd2
       dd5 = dd3*dd2
       dd6 = dd4*dd2
       !
       phi = phi1-(nn1*tan(phi1)/rr1)*(dd2/2.0-(5.0+3.0*tt1+  &
            10.0*cc1-4.0*cc1*cc1-9.0*ep2)*dd4/24.0+ &
            (61.0+90.0*tt1+298.0*cc1+45.0*tt1*tt1-252.0*ep2- &
            3.0*cc1*cc1)*dd6/720.0)
       lambda = dlamb0 + (dd - (1.0+2.0*tt1+cc1)*dd3/6.0 + &
            (5.0-2.0*cc1+28.0*tt1-3.0*cc1*cc1+8.0*ep2+24.0*tt1*tt1)* &
            dd5/120.0)/ cos(phi1)
    endif
    !
    !***  Convert phi/lambda to degrees
    !
    rlat = phi * 180.0 / M_PI
    rlon = lambda * 180.0 / M_PI
    !
    !***  All done
    !
    istat = 0
    return
  end subroutine utm2ll
  !
  !
  !
  logical function isdigit(a)
    implicit none
    character(len=1) :: a
    if(ichar(a).ge.ichar('0').and.ichar(a).le.ichar('9')) then
       isdigit = .true.
    else
       isdigit = .false.
    endif
    return
  end function isdigit
  !
  !
  !
  subroutine get_grid_zone (rlat,rlon,grid_zone,dlamb0)
    !******************************************************************
    !*
    !*    Solve for the grid zone, returns the central meridian
    !*
    !*****************************************************************
    implicit none
    real(rp) :: rlat,rlon,dlamb0
    character(len=3) :: grid_zone
    !
    integer(ip) :: long_zone,lat_zone
    !
    !*** First, let's take care of the polar regions
    !
    if(rlat.lt.(-80.0)) then
       if (rlon.lt.0.0) then
          grid_zone = '30A'
          dlamb0 = 0.0 * M_PI / 180.0
       else
          grid_zone = '31B'
          dlamb0 = 0.0 * M_PI / 180.0
       endif
       return
    else if(rlat.gt.84.0) then
       if(rlon.lt.0.0) then
          grid_zone = '30Y'
          dlamb0 = 0.0 * M_PI / 180.0
       else
          grid_zone = '31Z'
          dlamb0 = 0.0 * M_PI / 180.0
       endif
       return
    endif
    !
    !***  Now the special "X" grid
    !
    if(rlat.gt.72.0.and.rlon.gt.0.and.rlon.lt.42.0) then
       if(rlon.lt.9.0) then
          dlamb0 = 4.5
          grid_zone = '31X'
       else if(rlon.lt.21.0) then
          dlamb0 = 15.0 * M_PI / 180.0
          grid_zone = '33X'
       else if (rlon.lt.33.0) then
          dlamb0 = 27.0 * M_PI / 180.0
          grid_zone = '35X'
       else if (rlon.lt.42.0) then
          dlamb0 = 37.5 * M_PI / 180.0
          grid_zone = '37X'
       endif
       !
       return
    endif
    !
    !***  Handle the special "V" grid
    !
    if (rlat.gt.56.0.and.rlat.lt.64.0.and.rlon.gt.0.0 &
         .and.rlon.lt.12.0) then
       if (rlon.lt.3.0) then
          dlamb0 = 1.5 * M_PI / 180.0
          grid_zone = '31V'
       else if (rlon.lt.12.0) then
          dlamb0 = 7.5 * M_PI / 180.0
          grid_zone = '32V'
       endif
       return
    endif
    !
    !***  The remainder of the grids follow the standard rule
    !
    long_zone = int((rlon - (-180.0)) / 6.0) + 1
    dlamb0 = ((long_zone - 1)*6.0 + (-180.0)+3.0)*M_PI/180.0
    lat_zone = int((rlat-(-80.0))/8.0) + ichar('C')
    !
    if(lat_zone.gt.ichar('H')) lat_zone = lat_zone+1
    if(lat_zone.gt.ichar('N')) lat_zone = lat_zone+1
    !
    if(rlat.gt.80.0) lat_zone = ichar('X')
    !
    grid_zone(1:1) = char((long_zone/10) + ichar('0'))
    grid_zone(2:2) = char(mod(long_zone,10) + ichar('0'))
    grid_zone(3:3) = char(lat_zone)
    !
    return
  end subroutine get_grid_zone
  !
  !
  !
  subroutine get_lambda0 (grid_zone,dlamb0,istat)
    !********************************************************************
    !*
    !*    Given the grid zone, sets the central meridian, dlamb0
    !*
    !********************************************************************
    implicit none
    character(len=3) :: grid_zone
    integer(ip) :: istat
    real(rp)    :: dlamb0
    !
    integer(ip) :: long_zone,lat_zone
    !
    !***  Check the grid zone format
    !
    if (.not.isdigit(grid_zone(1:1)).or. &
         .not.isdigit(grid_zone(2:2))) then
       istat = -1
       return
    endif
    !
    long_zone = (ichar(grid_zone(1:1)) - ichar('0')) * 10 &
         + (ichar(grid_zone(2:2)) - ichar('0'))
    lat_zone  =  ichar(grid_zone(3:3))
    !
    !***  Take care of special cases
    !
    if(lat_zone.eq.ichar('A').or.lat_zone.eq.ichar('B').or. &
         lat_zone.eq.ichar('Y').or.lat_zone.eq.ichar('Z')) then
       dlamb0 = 0.0
       istat=0
       return
    elseif(lat_zone.eq.ichar('V')) then
       if(long_zone.eq.31) then
          dlamb0 = 1.5 * M_PI / 180.0
          istat=0
          return
       elseif(lat_zone.eq.32) then
          dlamb0 = 7.5 * M_PI / 180.0
          istat=0
          return
       endif
    elseif(lat_zone.eq.ichar('X')) then
       if(long_zone.eq.31) then
          dlamb0 = 4.5 * M_PI / 180.0
          istat=0
          return
       elseif(lat_zone.eq.33) then
          dlamb0 = 15.0 * M_PI / 180.0
          istat=0
          return
       elseif(lat_zone.eq.35) then
          dlamb0 = 27.0 * M_PI / 180.0
          istat = 0
          return
       elseif(lat_zone.eq.37) THEN
          dlamb0 = 37.5 * M_PI / 180.0
          istat = 0
          return
       elseif(lat_zone.eq.32.or.lat_zone.eq.34.or.lat_zone.eq.36) then
          istat = -1
          return
       endif
    endif
    !
    !***  Now handle standard cases
    !
    dlamb0 = ((long_zone-1)*6.0+(-180.0)+3.0) * M_PI / 180.0
    !
    !***  All done
    !
    istat = 0
    return
  end subroutine get_lambda0
  !
  !
  !
END MODULE Coordinates
