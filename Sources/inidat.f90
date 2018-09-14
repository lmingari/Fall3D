subroutine inidat
  !****************************************************************
  !*
  !*    Initializes data
  !*
  !****************************************************************
  use KindType
  use Numeric
  use Master
  implicit none
  !
  integer(ip) :: ilen,jlen,is
  !
  !***  CPU time
  !
  call cpu_time(cpu_time_begin)
  !
  !***  Initialization of variables (default values)
  !
  end_time       = 0.0_rp
  beg_time       = 0.0_rp
  safety_factor  = 0.9_rp
  dt             = 1.0_rp
  !
  modkv = 0                    ! Turbulence model
  modkh = 0                    ! Turbulence model
  modv  = 2                    ! Settling velocity model (2=Ganser)
  lmeth = 1                    ! Limiter method
  !
  rkh0 = 0.0_rp           ! Horizontal diffusion
  rkv0 = 0.0_rp           ! Vertical   diffusion
  !
  erumass = 0.0_rp    ! Total erupted mass
  outmass = 0.0_rp    ! Output mass
  time    = 0.0_rp    ! Time
  !
  memo = 0             ! Required memory
  !
  !*** Gets the number of cpu groups
  !
  number_of_cpu_groups = 0
  jlen = LEN_TRIM(c_cpu_groups)
  do ilen = 1,jlen
     call s1toi1(c_cpu_groups(ilen:ilen),is)
     number_of_cpu_groups = number_of_cpu_groups + (10**(jlen-ilen))*is
  end do
  !
  return
end subroutine inidat
!
!
!
subroutine s1toi1(string1,i1)
  !********************************************************
  !*
  !*   Decodes a character*1 string
  !*
  !********************************************************
  use KindType
  implicit none
  integer(ip) :: i1
  character(len=1) :: string1
  !
  if(string1 == '0') then
     i1 = 0
  else if(string1 == '1') then
     i1 = 1
  else if(string1 == '2') then
     i1 = 2
  else if(string1 == '3') then
     i1 = 3
  else if(string1 == '4') then
     i1 = 4
  else if(string1 == '5') then
     i1 = 5
  else if(string1 == '6') then
     i1 = 6
  else if(string1 == '7') then
     i1 = 7
  else if(string1 == '8') then
     i1 = 8
  else if(string1 == '9') then
     i1 = 9
  end if
  !
  return
end subroutine s1toi1
