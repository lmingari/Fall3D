subroutine writim(wordin)
  !*****************************************************
  !*
  !*     Writes to the log file
  !*
  !*****************************************************
  use KindType
  use Master
  use InpOut
  implicit none
  !
  integer(ip) :: ilen
  character(len=*)      :: wordin
  character(len=s_word) :: word
  !
  ilen = MIN(LEN_TRIM(wordin),s_word)
  word(:) = ' '
  word(1:ilen) = wordin(1:ilen)
  !
  !***   Writes
  !
  write(lulog,'(a,f15.0)') 'Results for '//word//' at time ',time(it)
  call flush(lulog)
  !
  return
end subroutine writim
