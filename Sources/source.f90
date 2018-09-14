subroutine source
  !***********************************************************************
  !*
  !*    Reads sources from the source file. It gets also
  !*    the mass flow rate at the current interval.
  !*
  !***********************************************************************
  use KindType
  use InpOut
  use Master
  use Numeric
  use Parallel
  implicit none
  !
  integer(ip) :: isou,ic
  integer(ip) :: iyr,imo,idy,ihr,imi,ise
  !
  !***  Proceed
  !
  if(.not.sourcetime) return
  !
  !***  Reads source data : xs,ys,zs,isrc,jsrc,ksrc,tmrat and source_time
  !
  call reasrc
  !
  !***  Computes the mass flow rate
  !
  flowrate      = sum(tmrat(1:nsou,1:nc))
  flowrate_part = sum(tmrat(1:nsou,1:np))
  if( ng > 0 ) then
     flowrate_gas = sum(tmrat(1:nsou,np+1:nc))
  else
     flowrate_gas = 0.0_rp
  end if
  !
  !***  Writes source information to the lst file
  !
  call addtime(ibyr,ibmo,ibdy,0,iyr,imo,idy,ihr,imi,ise,time)
  if( mpime == root ) write(nlst,5) idy ,month(imo ),iyr ,ihr ,(imi*60)+ise
5 format(/,                                                       &
       'SOURCE FILE READ              '                       ,/,   &
       '  From time                 : ',2x,i2.2,1x,a3,1x,i4.4,' at ',i2.2,' h ',i4.4,' s')
  !
  call addtime(ibyr,ibmo,ibdy,0,iyr,imo,idy,ihr,imi,ise,source_time)
  if( mpime == root ) write(nlst,6) idy ,month(imo ),iyr ,ihr ,(imi*60)+ise,INT(source_time-time), &
       flowrate_part,flowrate_gas,flowrate,nsou
6 format(  &
       '  To   time                 : ',2x,i2.2,1x,a3,1x,i4.4,' at ',i2.2,' h ',i4.4,' s',/, &
       '  Time interval (in sec)    : ',i9                    ,/,   &
       '  Mass flow rate particles  : ',e14.6                 ,/,   &
       '  Mass flow rate gas        : ',e14.6                 ,/,   &
       '  Mass flow rate total      : ',e14.6                 ,/,   &
       '  Number of sources         : ',i5,/)
  !
  if( mpime == root ) then
     do isou = 1,nsou
        write(nlst,10) isrc(isou),jsrc(isou),ksrc(isou),                 &
             (zlayer(ksrc(isou))+topg(isrc(isou),jsrc(isou))), &
             sum(tmrat(isou,1:np))
10      format('  I:',i3,'  J:',i3,'  K:',i3,' (',f6.0,' m a.s.l.)  MFR:',e12.4)
     end do
     write(nlst,20)
20   format(/,                                                 &
          '      Class    Diameter(mm)    Density        MFR  ')
     !
     do ic = 1,nc
        write(nlst,21) TRIM(classname(ic)),diam(ic)*1d3,rhop(ic),sum(tmrat(1:nsou,ic))
21      format(2x,a9,4x,2x,f8.3,2x,4x,f7.0,4x,e11.4)
     end do
     !
  end if
  !
  !***  Scales the source term (point dependent scaling)
  !
  do isou = 1,nsou
     tmrat(isou,1:nc) = tmrat(isou,1:nc)*Hm1(isrc(isou),jsrc(isou))
  end do
  !
  !***  If necessary, reads plume temperature and water flow rate
  !***  Aggragation parameters are recalculated later
  !
  if(aggregation) call settemplume
  !
  return
end subroutine source
!
!
!
subroutine settemplume
  !*****************************************************************************
  !*     
  !*   Sets plume temperature and water rate at the affected nodes
  !*
  !*****************************************************************************      
  use KindType
  use InpOut
  use Master
  use Parallel
  implicit none
  !
  logical     :: go_on
  integer(ip) :: timesec,itime1,itime2,ntemp,itemp,ix,iy,iz,ivoid
  real   (rp) :: templume,waterplume,rvoid
  !
  !***  Initializations
  !
  timesec       = INT(time)
  tplume(:,:,:) = 0.0_rp
  wplume(:,:,:) = 0.0_rp
  !
  !***  Master reads plume temperature file for the current interval
  !
  if( mpime == root ) then
     !
     !***    Opens the file
     !
     open(ntem,FILE=TRIM(ftem),STATUS='unknown',ERR=100)
     !
     !***    Reads until the required time
     !
     go_on = .true.
     do while(go_on)
        read(ntem,*,ERR=100) itime1,itime2
        read(ntem,*,ERR=100) ntemp
        if((timesec.ge.itime1).and.(timesec.lt.itime2)) then
           go_on = .false.
           do itemp = 1,ntemp 
              read(ntem,*,ERR=100) ix,iy,iz,templume,waterplume
              !
              tplume(ix,iy,iz) = templume      ! t plume in K
              wplume(ix,iy,iz) = waterplume    ! w in Kg/s
           end do
        else
           do itemp = 1,ntemp
              read(ntem,*,ERR=100) ivoid,ivoid,ivoid,rvoid,rvoid
           end do
        end if
     end do
     !
100  continue      
  end if
  !
  !*** Broadcast
  !
  call bcast( tplume, SIZE(tplume), root )
  call bcast( wplume, SIZE(wplume), root )
  !
  return
end subroutine settemplume
