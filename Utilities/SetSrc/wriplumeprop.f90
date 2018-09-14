  subroutine wriplumeprop
  !**************************************************************
  !*
  !*    Outputs plume porperties
  !*
  !*************************************************************
  use KindType
  use InpOut
  use Master
  use Plume
  use Coordinates
  implicit none
  !
  character(len=s_work) :: message
  integer(ip) :: is,ic,iz
  integer(ip), save :: ipass = 0
  real   (rp) :: x,y,z,th
  !
  !***  1. Writes characteristics of the plume
  !
  if(ipass.eq.0) then
     ipass = 1
     if(modv.eq.1) message = 'Arastoopour 82 '
     if(modv.eq.2) message = 'Ganser 93      '
     if(modv.eq.3) message = 'Wilson-Huang 79'
     if(modv.eq.4) message = 'Dellino 05     '
     !
     write(lures,1)
1    format(/, &
          '----------------------------------------------------',/, &
          '           CHARACTERISTICS OF THE PLUME             ',/, &
          '----------------------------------------------------',/)
     write(lures,10) X0_LL,Y0_LL,X0_UTM,Y0_UTM,(z0+dZ0)/1d3,TRIM(message), &
          xi,nc
10   format(/, &
          'CHARACTERISTICS OF THE COLUMN         ',/    ,/, &
          '  LONGITUDE OF THE VENT             : ',f12.7,/, &
          '  LATITUDE  OF THE VENT             : ',f12.7,/, &
          '  X-COORDINATE OF THE VENT  (UTM  ) : ',f11.0,/, &
          '  Y-COORDINATE OF THE VENT  (UTM  ) : ',f11.0,/, &
          '  ALTITUDE OF THE VENT      (km   ) : ',f11.4,/, &
          '  VELOCITY    MODEL                 : ',a    ,/, &
          '  xi                                : ',f11.4,/, &
          'GRANULOMETRIC DISTRIBUTION            '      ,/, &
          '  NUMBER OF PARTICLE CLASSES        : ',i4)
     write(lures,11) (diam(ic)*1d3,ic=1,nc)
     write(lures,12) (1d2*fc(ic)  ,ic=1,nc)
     write(lures,13) (sphe(ic)    ,ic=1,nc)
     write(lures,14) (psi (ic)    ,ic=1,nc)
11   format( &
          '  DIAMETER                  (mm   ) : ',50(f7.4,1x))
12   format( &
          '  PERCENTAGE                (in % ) : ',50(f7.1,1x))
13   format( &
          '  SPHERICITY                (-    ) : ',50(f7.2,1x))
14   format( &
          '  MODEL FACTOR              (-    ) : ',50(f7.2,1x))
  end if
  !
  !***  Writes results along the plume
  !
  write(lures,20) time,time2
20 format(/, &
       '--------------------  ',/, &
       'SIMULATION INTERVAL : ',i7,'    to  ',i7,/, &
       '--------------------  ')
  !
  write(lures,21) M0, &
       Zplum(np)-Z0-dZ0,Zplum(np), &
       Zplum(ns)-Z0-dZ0,Zplum(ns), &
       PSIa*180./pi
21 format( &
       'MFR                           : ',e12.4,/,&
       'NBL    HEIGHT                 : ',f6.0,' (',f6.0,' a.s.l.)',/, &
       'COLUMN HEIGHT                 : ',f6.0,' (',f6.0,' a.s.l.)',/, &
       'AVERAGED WIND DIRECTION (DEG) : ',f6.1,/)
  !
  do iz = 1,nz
     write(lures,22) zmodel(iz),Z0+dZ0+zmodel(iz),Vair(iz),Aair(iz),Tair(iz)-273.,Rair(iz)
22   format( &
          'HEIGHT a.vent: ',f6.0,' (',f6.0,' a.s.l)   WIND SPEED: ',f6.1, &
          '   DIRECTION: ',f5.0,'   AIR TEMPERATURE(C):',f6.1, '   AIR DENSITY:',f6.3)
  end do
  !
  write(lures,100)
100 format(/, &
       'RESULTS :  x,y,z(km a.s.l),Th(deg),u(ms),r(m) T(C)',/)
  !
  do is=1,ns
     x = Xplum(is)
     y = Yplum(is)
     z = Zplum(is)
     th= Hplum(is)
     write(lures,110) x, &
          y, &
          z/1d3, &
          th*180.0/pi, &
          Uplum(is), &
          Rplum(is), &
          Tplum(is)-273.15
110  format(7(1x,e12.5))
  end do
  !
  return
  end subroutine wriplumeprop
