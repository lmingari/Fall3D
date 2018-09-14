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
  integer(ip), parameter  :: lutgs = 12    ! tgs  file
  integer(ip), parameter  :: lugrn = 13    ! Data file
  integer(ip), parameter  :: lusrc = 14    ! DBS  file
  integer(ip), parameter  :: ludbs = 15    ! top  file
  integer(ip), parameter  :: lures = 16    ! res  file (plume)
  integer(ip), parameter  :: lumas = 17    ! mas  file (plume)
  integer(ip), parameter  :: luhei = 18    ! hei  file (plume)
  integer(ip), parameter  :: lutem = 19    ! tem  file (plume)
  !
  !***  File names
  !
  character(len=s_file) :: luinpname,lutgsname,ludbsname,ludepname          ! Input files
  character(len=s_file) :: lulogname,lugrnname,lusrcname                    ! Output files
  character(len=s_file) :: luresname,luheiname,lumasname,lutemname          ! Output files for plume case
  character(len=s_file) :: lusncname,lusgnname                              ! Output files for resuspension case
  integer(ip)           :: iarg
  !
  !***  List of Warnings
  !
  integer(ip), parameter :: maxwarn = 100
  integer(ip)            :: nwarn = 0
  character(len=s_mess)  :: warning(maxwarn)
  !
END MODULE InpOut
