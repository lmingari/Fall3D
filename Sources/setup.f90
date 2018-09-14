subroutine setup
  !**************************************************************************
  !*
  !*    Set up the run
  !*
  !**************************************************************************
  use Master
  use InpOut
  use Numeric
  use Parallel,  only: parallel_startup,parallel_setgroups,bcast,mpime,root,nproc
  implicit none
  !
  integer(ip)       :: info
  character(len=10) :: rdate,stime
  character(len=8 ) :: rtime,sdate
  !
  !***  Initializes parallel data
  !
  call parallel_startup
  !
  !***  Initializes data
  !
  call inidat
  !
  !***  Clock
  !
  call DATE_AND_TIME(sdate,stime)
  rdate='  -  -    '
  rdate(1:2)=sdate(5:6)
  rdate(4:5)=sdate(7:8)
  rdate(7:10)=sdate(1:4)
  rtime='  :  :  '
  rtime(1:2)=stime(1:2)
  rtime(4:5)=stime(3:4)
  rtime(7:8)=stime(5:6)
  !
  !***  Opens the lst file
  !
  if( mpime == root ) then
     open(nlst,file=TRIM(flst),iostat=info,status='unknown')
  end if
  call bcast( info, 1, root )
  if( info /= 0 ) goto 100
  !
  !***  Writes header
  !
  if( mpime == root ) then
     if( nproc == 1 ) then
        write(nlst,1) version,' SERIAL',nproc,rdate,rtime
     else
        write(nlst,1) version,' PARALLEL',nproc,rdate,rtime
     end if
1    format(&
          '---------------------------------------------------------',/,  &
          '                                                         ',/,  &
          '                          FALL3D                         ',/,  &
          '                                                         ',/, &
          '  Copyright (C) 2013 GNU General Public License version 3',/, &
          '  Arnau Folch, Antonio Costa, Giovanni Macedonio         ',/, &
          '  See licence for details                                ',/, &
          '---------------------------------------------------------',/,  &
          '  Code version          : ',f4.1,a                       ,/,  &
          '  Number of processors  : ',i4                           ,/,  &
          '  Starting date         : ',a10,' at time: ',a8)
     !
     write(nlst,2) TRIM(finp),TRIM(fsrc),TRIM(fgrn),TRIM(fdbs),TRIM(fpts), &
                   TRIM(flst),TRIM(fres_nc),TRIM(frst)
2    format(/, &
          '  INPUT files             '  ,/, &
          '  -----------             '  ,/, &
          '  FALL3D input    file  : ',a,/, &
          '  Source          file  : ',a,/, &
          '  Granulometry    file  : ',a,/, &
          '  Meteo           file  : ',a,/, &
          '  Tracking points file  : ',a,/, &
          '                          '  ,/, &
          '  OUTPUT files            '  ,/, &
          '  ------------            '  ,/, &
          '  Log             file  : ',a,/, &
          '  Results         file  : ',a,/, &
          '  Restart         file  : ',a)
  end if
  !
  !***  Reads data from the granulometry input file
  !
  call reagrn
  !
  !***  Reads data from the input file
  !
  call readat
  !
  !***  Reads grid data and allocates memory for arrays related to nx,ny,nz
  !
  call reagrd
  !
  !***  Reads output strategy
  !
  call reaout
  !
  !***  Checks input data
  !
  call chkdat
  !
  !***  Sets the parallel groups
  !
  call parallel_setgroups( number_of_cpu_groups )
  !
  !***  Computes topography bounds
  !
  tpgmax = MAXVAL(topg)
  tpgmin = MINVAL(topg)
  ztop   = tpgmax + zlayer(nz)
  !
  !***  If necessary, reads the restart file
  !
  if(restart) call rearst
  !
  !***  Writes inputs to the lst file
  !
  if( mpime == root ) call wridat
  !
  return
  !
  !*** List of errors
  !
100 call runend('Can not open file: '//TRIM(flst))
    end subroutine setup
