  subroutine writim
  !****************************************************************
  !*
  !*    Writes time evolution in the log file
  !*
  !****************************************************************
  use KindType
  use InpOut
  use Numeric
  use Master
  use Parallel
  implicit none
  !
  integer(ip) :: iyr,imo,idy,ihr,imi,ise
  !
  !***  Computes mass in the atmosphere and in the deposit
  !
  call cmass3d(c,cload,cdry,cwet,nx,ny,nz,nc,dX1,dX2,dz,massvol,massdep,massdry,masswet)
  !
  !***  Writes summary to the log file
  !
  call addtime(ibyr,ibmo,ibdy,0,iyr,imo,idy,ihr,imi,ise,time)
  !
  if( mpime == root )  &
       write(nlst,31) iiter,dtc,dt,int((time-beg_time)/60.0_rp),          &
       100.0_rp*((time-beg_time)/(end_time-beg_time)),     &
       idy ,month(imo ),iyr ,ihr ,imi ,ise,                &
       massvol,1d2*massvol/(massvol+massdep+outmass),   &
       outmass,1d2*outmass/(massvol+massdep+outmass),   &
!
       massdry,1d2*massdry/(massdep),   &
       masswet,1d2*masswet/(massdep),   &
!
       massdep,1d2*massdep/(massvol+massdep+outmass),   &
       (massvol+massdep+outmass),erumass
  !
31 format(/,                                                       &
       '-> Iteration                 : ',i8                   ,/,  &
       '   Critical time step        : ',e13.6                ,/,  &
       '   Time step                 : ',e13.6                ,/,  &
       '   Elapsed time              : ',i7,' min   (',f7.2,'%)',/,  &
       '   Current time              : ',i2,1x,a3,1x,i4,' at ',i2,':',i2.2,':',i2.2 ,/, &
       '   -------                     '                      ,/,  &
       '   (1) Approx. mass inside  domain : ',e13.6,' (',f7.2,'%)' ,/,  &
       '   (2) Approx. mass outside domain : ',e13.6,' (',f7.2,'%)' ,/,  &
       '   (3) Approx. mass at the deposit   '                      ,/,  &
       '       (3.1) Dry deposition        : ',e13.6,'             (',f7.2,'%)' ,/,  &
       '       (3.2) Wet deposition        : ',e13.6,'             (',f7.2,'%)' ,/,  &
       '       (3.3) Total                 : ',e13.6,' (',f7.2,'%)' ,/,  &
       '       ---------------------------  '                      ,/,  &
       '   Summ (1)+(2)+(3.3)              : ',e13.6                ,/,  &
       '   Erupted mass                    : ',e13.6)
  !
  !***  Forces writting
  !
  if( mpime == root ) call flush(nlst)
  !
  end subroutine writim
