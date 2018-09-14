subroutine runend(message)
  !*************************************************************************
  !*
  !*    This routine stops the run and writes termination status
  !*
  !*************************************************************************
  use kindType
  use InpOut
  use Master
  use Parallel
  implicit none
  !
  integer(ip)      :: iwarn
  character(len=*) :: message
  !
  !***  Writes header
  !
  if( mpime == root ) write(nlst,10)
10 format('---------------------------------------------------',/,  &
       '                                                   ',/   &
       '             FALL3D : END OF SIMULATION            ',/,  &
       '                                                   ',/,  &
       '---------------------------------------------------',/)
  !
  !***  Writes EOF
  !
  call cpu_time(cpu_time_end)
  !
  !***  Writes warning list
  !
  if( mpime == root ) then
     write(nlst,21) nwarn
21   format('  WARNINGS : ',i2)
     do iwarn = 1,nwarn
        write(nlst,22) iwarn,TRIM(warning(iwarn))
     end do
22   format(3x,i2,' : ',a)
  end if
  !
  !*** Write message and stop the run.
  !
  if(TRIM(message) /= 'OK') then             ! Abnormal termination
     !
     if( mpime == root ) write(nlst,31) TRIM(message),(cpu_time_end-cpu_time_begin)
31   format(&
          '  STATUS   :  ABNORMAL TERMINATION',/,                &
          '  ERROR    : ',a                   ,/,                &
          '  CPU TIME : ',f9.1, 'sec')
  else                                     ! Normal termination
     !
     if( mpime == root ) write(nlst,32) (cpu_time_end-cpu_time_begin)
32   format(&
          '  STATUS   :  NORMAL TERMINATION',/,                  &
          '  CPU TIME : ',f9.1, 'sec')
  end if
  !
  !*** Closes LST file
  !
  if( mpime == root ) then
     close(nlst)
  end if
  !
  !*** Ends the run
  !
  call parallel_hangup
  !
  stop
end subroutine runend
!
!
!
subroutine wriwar(message)
  !*************************************************************************
  !*
  !*    This routine writes a warning message to the warnings list
  !*
  !*************************************************************************
  use InpOut
  implicit none
  character(len=*) :: message
  !
  nwarn = nwarn + 1
  nwarn = MIN(nwarn,maxwarn)
  warning(nwarn) = message(1:LEN_TRIM(message))
  !
  return
end subroutine wriwar
