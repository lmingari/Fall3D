!***************************************************************
!*
!*		Module for input/output
!*
!***************************************************************
MODULE InpOut
  use KindType
  IMPLICIT NONE
  SAVE
  !
  !***  File logical units
  !
  integer(ip), parameter  :: lulog = 10    ! log  file
  integer(ip), parameter  :: luinp = 11    ! inp  file
  integer(ip), parameter  :: lugrn = 12    ! grn file
  !
  !***  File names
  !
  character(len=s_file) :: lulogname,luinpname,lugrnname
  integer(ip)           :: iarg
  !
  !***  List of Warnings
  !
  integer(ip), parameter :: maxwarn = 100
  integer(ip)            :: nwarn = 0
  character(len=s_mess)  :: warning(maxwarn)
  !
END MODULE InpOut
