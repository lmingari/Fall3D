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
  real(rp), parameter :: LOWER_EPS_LIMIT = 1e-14_rp
  real(rp), parameter :: M_PI = 3.14159265358979323846_rp     ! pi
  real(rp), parameter :: M_PI_2 = 1.57079632679489661923_rp   ! pi/2
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
    !***********************************************************************
    implicit none
    integer(ip) :: istat,datum
    real(rp)    :: rlat,rlon, utmx, utmy
    character(len=3) :: grid_zone
    !
    character(len=3) :: utmgz
    integer(ip) :: i
    real(rp) :: a,b,f,e,e2,e4,e6,phi,phi0,t,x,y, rho, cc,cc2,dlamb0
    real(rp) :: aa,aa2,aa3,aa4,aa5,aa6,lambda, k0, mm, mm0, nn, ep2,tt2,tt
    !
    if(datum.eq.CLARKE_1866_DATUM) then
       a = 6378206.4_rp        ! semimajor axis of ellipsoid (meters)
       b = 6356583.8_rp        ! semiminor axis of ellipsoid (meters)
    elseif(datum.eq.GRS_80_DATUM) then
       a = 6378137._rp
       b = 6356752.3_rp
    elseif(datum.eq.WGS_84_DATUM) then
       a = 6378137._rp
       b = 6356752.31425_rp ! Used by GoogleEarth
    elseif(datum.eq.NAD83_DATUM) then
       a = 6378137._rp
       b = 6356752.3141_rp
    elseif(datum.eq.AIRY_1830_DATUM) then
       a = 6377563.4_rp
       b = 6356256.9_rp
    else
       istat = -1
       return
    endif
    !
    !***  Calculate flatness and eccentricity
    !
    f =1.0_rp - (b/a)
    e2 = (2.0_rp - f) * f
    e = sqrt(e2)
    e4 = e2*e2
    e6 = e4*e2
    !
    !***  Convert latitude/longitude to radians
    !
    phi = rlat * M_PI / 180.0_rp
    lambda = rlon * M_PI / 180.0_rp
    !
    !***  Figure out the UTM zone, as well as dlamb0
    !
    call get_grid_zone (rlat, rlon, utmgz, dlamb0)
    phi0 =0.0_rp
    !
    !***  See if this will use UTM or UPS
    !
    if (rlat.gt.84.0_rp) then
       !     use Universal Polar Stereographic Projection (north polar aspect)
       k0 = 0.994_rp

       t = sqrt(((1.0_rp-sin(phi))/(1.0_rp+sin(phi))) *  &
            ((1.0_rp+e*sin(phi))/(1.0_rp-e*sin(phi)))**e)
       rho = 2.0_rp*a*k0*t/sqrt((1.0_rp+e)**(1.0_rp+e)*(1.0_rp-e)**(1.0_rp-e))
       x = rho * sin(lambda - dlamb0)
       y = -rho * cos(lambda - dlamb0)
       !
       !***     Apply false easting/northing
       !
       x = x + 2d6
       y = y + 2d6

    elseif (rlat.lt.(-80.0_rp)) then

       !     use Universal Polar Stereographic Projection (south polar aspect)

       phi = -phi
       lambda = -lambda
       dlamb0 = -dlamb0

       k0 = 0.994_rp
       t = sqrt(((1.0_rp-sin(phi))/(1.0_rp+sin(phi))) *  &
            ((1.0_rp+e*sin(phi))/(1.0_rp-e*sin(phi)))**e)
       rho =2.0_rp*a*k0*t/ sqrt((1.0_rp+e)**(1.0_rp+e)*(1.0_rp-e)**(1.0_rp-e))

       x = rho * sin(lambda - dlamb0)
       y = -rho * cos(lambda - dlamb0)
       x = -x
       y = -y
       !
       !***    Apply false easting/northing
       !
       x = x + 2d6
       y = y + 2d6

    else

       !     Use UTM
       !     set scale on central median (0.9996 for UTM)

       k0 = 0.9996_rp
       mm = a * ((1.0_rp-e2/4.0_rp - 3.0_rp*e4/64.0_rp - 5.0_rp*e6/256.0_rp) * phi -  &
            (3.0_rp*e2/8.0_rp+3.0_rp*e4/32.0_rp+45.0_rp*e6/1024.0_rp)*sin(2*phi)+    &
            (15.0_rp*e4/256.0_rp + 45.0_rp*e6/1024.0_rp) * sin(4.0_rp*phi) -   &
            (35.0_rp*e6/3072.0_rp) * sin(6.0_rp*phi))

       mm0 = a * ((1.0_rp-e2/4._rp - 3._rp*e4/64._rp - 5._rp*e6/256._rp)*phi0 -  &
            (3._rp*e2/8._rp+3._rp*e4/32._rp+45._rp*e6/1024._rp)*sin(2._rp*phi0) +  &
            (15._rp*e4/256._rp + 45._rp*e6/1024._rp) * sin(4._rp*phi0) -        &
            (35._rp*e6/3072._rp) * sin(6._rp*phi0))

       aa = (lambda - dlamb0) * cos(phi)
       aa2 = aa*aa
       aa3 = aa2*aa
       aa4 = aa2*aa2
       aa5 = aa4*aa
       aa6 = aa3*aa3

       ep2 = e2 / (1.0_rp - e2)
       nn = a / sqrt(1.0_rp - e2*sin(phi)**2)
       tt = tan(phi)**2
       tt2 = tt*tt
       cc = ep2 * cos(phi)**2
       cc2 = cc*cc

       x = k0*nn*(aa+(1.0_rp-tt+cc)*aa3/6._rp+(5._rp-18._rp*tt+tt2+ &
            72._rp*cc-58._rp*ep2)*aa5/120._rp)
       y = k0*(mm-mm0+nn*tan(phi)*(aa2/2._rp+(5._rp-tt+9._rp*cc+ &
            4._rp*cc2)*aa4/24._rp+(61._rp-58._rp*tt+tt2+600._rp*cc-  &
            330._rp*ep2)*aa6/720._rp))
       !
       !***     Apply false easting and northing
       !
       x = x + 5d5
       if(y.lt.0._rp) y = y + 1d7
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
       a = 6378206.4_rp        ! semimajor axis of ellipsoid (meters)
       b = 6356583.8_rp        ! semiminor axis of ellipsoid (meters)
    elseif(datum.eq.GRS_80_DATUM) then
       a = 6378137._rp
       b = 6356752.3_rp
    elseif(datum.eq.WGS_84_DATUM) then
       a = 6378137.0_rp
       b = 6356752.31425_rp
    else
       istat = -1
       return
    endif
    !
    !***  Calculate flatness and eccentricity
    !
    f =1.0_rp - (b/a)
    e2 = (2._rp - f) * f
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
       x = utmx - 2d6
       y = utmy - 2d6
       !
       !***  Solve for inverse equations
       !
       k0 = 0.994_rp
       rho = sqrt(x*x + y*y)
       t = rho*sqrt((1.0_rp+e)**(1.0_rp+e)*(1.0_rp-e)**(1.0_rp-e))/(2._rp*a*k0)
       !
       !***     Solve for latitude and longitude
       !
       chi = M_PI_2 - 2._rp * atan(t)
       phit = chi + (e2/2._rp+5._rp*e4/24._rp+e6/12._rp+13._rp*e8/360._rp)* &
            sin(2._rp*chi)+(7._rp*e4/48._rp+29._rp*e6/240._rp+811._rp*e8/11520._rp)* &
            sin(4._rp*chi)+(7._rp*e6/120._rp+81._rp*e8/1120._rp)*sin(6._rp*chi)+ &
            (4279._rp*e8/161280._rp)*sin(8._rp*chi)
       !
10     continue
       phi = phit
       phit = M_PI_2 - 2._rp*atan(t*((1.0_rp-e*sin(phi))/(1.0_rp+e*sin(phi)))**(e/2._rp))

       if(abs(phi-phit).le.LOWER_EPS_LIMIT) goto 20
       goto 10
20     continue

       lambda = dlamb0 + atan2(x,-y)
    elseif(lat_zone.eq.ichar('A').or.lat_zone.eq.ichar('B')) then
       !       south polar aspect
       !       Subtract the false easting/northing

       x = -(utmx - 2d6)
       y = -(utmy - 2d6)
       !
       !***     Solve for inverse equations
       !
       k0 = 0.994_rp
       rho = sqrt (x*x + y*y)
       t = rho*sqrt((1.0_rp+e)**(1.0_rp+e)*(1.0_rp-e)**(1.0_rp-e))/(2._rp*a*k0)
       !
       !***    Solve for latitude and longitude
       !
       chi = M_PI_2 - 2._rp * atan(t)
       phit = chi + (e2/2._rp + 5._rp*e4/24._rp + e6/12._rp + 13._rp*e8/360._rp)* &
            sin(2._rp*chi) + (7._rp*e4/48._rp + 29._rp*e6/240._rp + &
            811._rp*e8/11520._rp) * sin(4._rp*chi) + (7._rp*e6/120._rp + &
            81._rp*e8/1120._rp) * sin(6._rp*chi) + (4279._rp*e8/161280._rp)* &
            sin(8._rp*chi)
       !
11     continue
       phi = phit
       phit = M_PI_2 - 2._rp*atan(t*((1.0_rp-e*sin(phi))/(1.0_rp+e*sin(phi)))**(e/2._rp))
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
       k0 = 0.9996_rp
       !
       !***   Remove false eastings/northings
       !
       x = utmx - 5d5
       y = utmy
       !
       if(lat_zone.gt.ichar('B').and.lat_zone.lt.ichar('N')) &
            y = y - 1d7       ! southern hemisphere
       !
       !     Calculate the footpoint latitude
       !
       phi0 = 0.0_rp
       e1 = (1.0_rp - sqrt(1.0_rp-e2))/(1.0_rp + sqrt(1.0_rp-e2))
       e12 = e1*e1
       e13 = e1*e12
       e14 = e12*e12
       !
       mm0 = a * ((1.0_rp-e2/4._rp - 3._rp*e4/64._rp - 5._rp*e6/256._rp)*phi0- &
            (3._rp*e2/8._rp+3._rp*e4/32._rp+45._rp*e6/1024._rp)*sin(2._rp*phi0)+ &
            (15._rp*e4/256._rp + 45._rp*e6/1024._rp)*sin(4._rp*phi0)- &
            (35._rp*e6/3072._rp) * sin(6._rp*phi0))
       mm = mm0 + y/k0
       mu = mm/(a*(1.0_rp-e2/4._rp-3*e4/64._rp-5*e6/256._rp))
       !
       phi1 = mu + (3._rp*e1/2._rp - 27._rp*e13/32._rp) * sin(2._rp*mu)+ &
            (21._rp*e12/16._rp - 55._rp*e14/32._rp) * sin(4._rp*mu)+ &
            (151._rp*e13/96._rp) * sin(6._rp*mu)+ &
            (1097._rp*e14/512._rp) * sin(8._rp*mu)
       !
       !***    Now calculate lambda and phi
       !
       ep2 = e2/(1.0_rp-e2)
       cc1 = ep2*cos(phi1)**2
       tt1 = tan(phi1)**2
       nn1 = a / sqrt(1.0_rp-e2*sin(phi1)**2)
       rr1 = a * (1.0_rp-e2)/(1.0_rp-e2*sin(phi1)**2)**(1.5_rp)
       dd = x / (nn1 * k0)
       !
       dd2 = dd*dd
       dd3 = dd*dd2
       dd4 = dd2*dd2
       dd5 = dd3*dd2
       dd6 = dd4*dd2
       !
       phi = phi1-(nn1*tan(phi1)/rr1)*(dd2/2._rp-(5._rp+3._rp*tt1+  &
            10._rp*cc1-4._rp*cc1*cc1-9._rp*ep2)*dd4/24._rp+ &
            (61._rp+90._rp*tt1+298._rp*cc1+45._rp*tt1*tt1-252._rp*ep2- &
            3._rp*cc1*cc1)*dd6/720._rp)
       lambda = dlamb0 + (dd - (1.0_rp+2._rp*tt1+cc1)*dd3/6._rp + &
            (5._rp-2._rp*cc1+28._rp*tt1-3._rp*cc1*cc1+8._rp*ep2+24._rp*tt1*tt1)* &
            dd5/120._rp)/ cos(phi1)
    endif
    !
    !***  Convert phi/lambda to degrees
    !
    rlat = phi * 180.0_rp / M_PI
    rlon = lambda * 180.0_rp / M_PI
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
    if(rlat.lt.(-80._rp)) then
       if (rlon.lt.0._rp) then
          grid_zone = '30A'
          dlamb0 =0.0_rp * M_PI / 180.0_rp
       else
          grid_zone = '31B'
          dlamb0 =0.0_rp * M_PI / 180.0_rp
       endif
       return
    else if(rlat.gt.84._rp) then
       if(rlon.lt.0._rp) then
          grid_zone = '30Y'
          dlamb0 =0.0_rp * M_PI / 180.0_rp
       else
          grid_zone = '31Z'
          dlamb0 =0.0_rp * M_PI / 180.0_rp
       endif
       return
    endif
    !
    !***  Now the special "X" grid
    !
    if(rlat.gt.72._rp.and.rlon.gt.0.and.rlon.lt.42._rp) then
       if(rlon.lt.9._rp) then
          dlamb0 = 4.5_rp
          grid_zone = '31X'
       else if(rlon.lt.21._rp) then
          dlamb0 = 15._rp * M_PI / 180.0_rp
          grid_zone = '33X'
       else if (rlon.lt.33._rp) then
          dlamb0 = 27._rp * M_PI / 180.0_rp
          grid_zone = '35X'
       else if (rlon.lt.42._rp) then
          dlamb0 = 37.5_rp * M_PI / 180.0_rp
          grid_zone = '37X'
       endif
       !
       return
    endif
    !
    !***  Handle the special "V" grid
    !
    if (rlat.gt.56._rp.and.rlat.lt.64._rp.and.rlon.gt.0._rp &
         .and.rlon.lt.12._rp) then
       if (rlon.lt.3._rp) then
          dlamb0 = 1.5_rp * M_PI / 180.0_rp
          grid_zone = '31V'
       else if (rlon.lt.12._rp) then
          dlamb0 = 7.5_rp * M_PI / 180.0_rp
          grid_zone = '32V'
       endif
       return
    endif
    !
    !***  The remainder of the grids follow the standard rule
    !
    long_zone = int((rlon - (-180.0_rp)) / 6._rp) + 1
    dlamb0 = ((long_zone - 1)*6._rp + (-180.0_rp)+3._rp)*M_PI/180.0_rp
    lat_zone = int((rlat-(-80._rp))/8._rp) + ichar('C')
    !
    if(lat_zone.gt.ichar('H')) lat_zone = lat_zone+1
    if(lat_zone.gt.ichar('N')) lat_zone = lat_zone+1
    !
    if(rlat.gt.80._rp) lat_zone = ichar('X')
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
       dlamb0 =0.0_rp
       istat=0
       return
    elseif(lat_zone.eq.ichar('V')) then
       if(long_zone.eq.31) then
          dlamb0 = 1.5_rp * M_PI / 180.0_rp
          istat=0
          return
       elseif(lat_zone.eq.32) then
          dlamb0 = 7.5_rp * M_PI / 180.0_rp
          istat=0
          return
       endif
    elseif(lat_zone.eq.ichar('X')) then
       if(long_zone.eq.31) then
          dlamb0 = 4.5_rp * M_PI / 180.0_rp
          istat=0
          return
       elseif(lat_zone.eq.33) then
          dlamb0 = 15._rp * M_PI / 180.0_rp
          istat=0
          return
       elseif(lat_zone.eq.35) then
          dlamb0 = 27._rp * M_PI / 180.0_rp
          istat = 0
          return
       elseif(lat_zone.eq.37) THEN
          dlamb0 = 37.5_rp * M_PI / 180.0_rp
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
    dlamb0 = ((long_zone-1)*6._rp+(-180.0_rp)+3._rp) * M_PI / 180.0_rp
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
