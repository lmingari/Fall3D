  subroutine chkdat
  !*****************************************************************
  !*
  !*    Checks the input data
  !*
  !*****************************************************************
  use Parallel
  use Master, only: number_of_cpu_groups,nc
  implicit none
  !
  if( mod( nproc, number_of_cpu_groups ) /= 0 .or. number_of_cpu_groups <= 0 ) then
     call runend('number of CPUs and number of groups are incompatible')
  end if
  !
  if( number_of_cpu_groups > nc ) then
     call runend('number of groups greater than the number of particle classes')
  end if
  !
  return
  end subroutine chkdat
