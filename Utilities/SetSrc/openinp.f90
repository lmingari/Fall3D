  subroutine openinp
  !**************************************************************
  !*
  !*    Opens the control file
  !*
  !**************************************************************
  use InpOut
  use Master
  implicit none
  !
  character(len=10) :: rdate
  character(len=8 ) :: rtime
  !
  !***  Clock
  !
  call datem(rdate,rtime)
  !
  !***  Opens and writes the log file
  !
  open(lulog,file=TRIM(lulogname),status='unknown')
  write(lulog,1) version,rdate,rtime
   1  format('---------------------------------------------------------',/, &
             '                   PROGRAM SETSRC                        ',/, &
             '                   FALL3D Version ',f4.1,'               ',/, &
             '                                                         ',/, &
             '  Copyright (C) 2013 GNU General Public License version 3',/, &
             '  Arnau Folch, Antonio Costa, Giovanni Macedonio         ',/, &
             '  See licence for details                                ',/, &
             '---------------------------------------------------------',/, &
             '  Starting date         : ',a10,' at time: ',a8)
  !
  write(lulog,2) TRIM(luinpname),TRIM(lutgsname), &
                 TRIM(lugrnname),TRIM(lusrcname),TRIM(ludbsname),TRIM(ludepname)
2 format(/, &
       '  Input          file   : ',a,/, &
       '  TGSG           file   : ',a,/, &
       '  Granulometry   file   : ',a,/, &
       '  Source         file   : ',a,/, &
       '  Meteo Database file   : ',a,/, &
       '  Deposit        file   : ',a,' (For resuspension only)  ')
  !
  !***  Opens the source file
  !
  open(lusrc,file=TRIM(lusrcname),status='unknown',err=102)
  !
  return
102 call runend('Error opening file '//TRIM(lusrcname))
  return
  !
end subroutine openinp
!
subroutine datem(rdate,rtime)
  !***********************************************************************
  !*
  !*    Get system date and time from system clock
  !*
  !*    rtime*8    HH:MM:SS
  !*    rdate*10   MM-DD-YYYY
  !*
  !***********************************************************************
  implicit none
  !
  character(len=8)  :: rtime,sdate
  character(len=10) :: rdate,stime
  !
  call DATE_AND_TIME(sdate,stime)
  !
  rdate='  -  -    '
  rdate(1:2 )=sdate(5:6)
  rdate(4:5 )=sdate(7:8)
  rdate(7:10)=sdate(1:4)
  rtime='  :  :  '
  rtime(1:2 )=stime(1:2)
  rtime(4:5 )=stime(3:4)
  rtime(7:8 )=stime(5:6)
  !
  return
end subroutine datem
