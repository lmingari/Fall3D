!***************************************************************
!*
!* Module for input/output
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
  integer(ip), parameter  :: ludat = 12    ! Data file
  integer(ip), parameter  :: ludbs = 13    ! DBS  file
  integer(ip), parameter  :: lutop = 14    ! top  file
  !
  !***  File names
  !
  character(len=s_file) :: lulogname,luinpname,ludatname,ludbsname,lutopname
  integer(ip)           :: iarg
  !
  !***  Type of data
  !
  character(len=25) :: TYPE_DATA
  !
  !***  List of Warnings
  !
  integer(ip), parameter :: maxwarn = 100
  integer(ip)            :: nwarn = 0
  character(len=s_mess)  :: warning(maxwarn)
  !
END MODULE InpOut
