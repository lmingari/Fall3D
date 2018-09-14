subroutine wriplumeheight
  !**************************************************************
  !*
  !*    Outputs plume height
  !*
  !*************************************************************
  use KindType
  use Master
  use InpOut
  implicit none
  !
  integer(ip), save :: ipass = 0
  !
  !***  1. Writes characteristics of the plume
  !
  if(ipass.eq.0) then
     ipass = 1
     write(luhei,1)
1    format(/, &
          '----------------------------------------------------',/, &
          'Time 1  Time 2  NBL HEIGHT  COLUMN HEIGHT   MFR     ',/, &
          '  (h)     (h)   (km a.vent)  (km a.vent)   (kg/s)   ',/, &
          '----------------------------------------------------',/)
  end if
  !
  !***  Writes results along the plume
  !
  write(luhei,21) time/36e2, &
       min(time2,ieend1)/36e2, &
       (Zplum(np)-Z0-dZ0)/1d3, &
       (Zplum(ns)-Z0-dZ0)/1d3, &
       M0
21 format(f7.2,2x,f7.2,2x,f8.3,2x,f8.3,2x,e12.4)
  !
  return
end subroutine wriplumeheight
