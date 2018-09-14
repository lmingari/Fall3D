module odepack
  !
  ! Simple module for interfacing lsode with F90 programs (by G.Macedonio)
  !
  ! NOTE: For SetDbs the Real*4 version is commented
  !
  implicit none

  private :: lsode_real4
  private :: lsode_real8

  interface lsode

     module procedure lsode_real4
     module procedure lsode_real8

  end interface

contains

  subroutine lsode_real4 (f, neq, y, t, tout, itol, rtol, atol, itask, &
       istate, iopt, rwork, lrw, iwork, liw, jac, mf)
    implicit none
    external :: f, jac
    integer :: liw,lrw
    integer :: neq(*), itol, itask, istate, iopt, iwork(liw), mf
    real(4) :: y(*), t, tout, rtol(*), atol(*), rwork(lrw)
    !
    ! @@@@ STOP IF REAL*4 VERSION
    stop 'SLSODE NOT IMPLEMENTED IN ODEPACK.F90'
    !  UNCOMMENT IF NEEDED
    !
    !    call slsode (f, neq, y, t, tout, itol, rtol, atol, itask, &
    !         istate, iopt, rwork, lrw, iwork, liw, jac, mf)

  end subroutine lsode_real4

  subroutine lsode_real8 (f, neq, y, t, tout, itol, rtol, atol, itask, &
       istate, iopt, rwork, lrw, iwork, liw, jac, mf)
    implicit none
    external :: f, jac
    integer :: liw,lrw
    integer :: neq(*), itol, itask, istate, iopt, iwork(liw), mf
    real(8) :: y(*), t, tout, rtol(*), atol(*), rwork(lrw)
    call dlsode (f, neq, y, t, tout, itol, rtol, atol, itask, &
         istate, iopt, rwork, lrw, iwork, liw, jac, mf)
  end subroutine lsode_real8

end module odepack
