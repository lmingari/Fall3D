subroutine wriplumefmas
  !**************************************************************
  !*
  !*    Outputs mass fallen from the eruptive column
  !*
  !**************************************************************
  use KindType
  use Master
  use InpOut
  implicit none
  !
  integer(ip), save :: ipass = 0
  integer(ip) :: is,ic
  real   (rp) :: z,Mfi,Mtcon,Mtumb
  real   (rp), save, allocatable :: Mcon(:),Mumb(:)
  !
  !***  1. Writes mass falled from the plume
  !
  if(ipass.eq.0) then
     ipass = 1
     allocate(Mcon(nc))
     allocate(Mumb(nc))
     !
     write(lumas,1)
1    format(/, &
          '----------------------------------------------------',/, &
          '           MASS FALLEN FROM THE PLUME               ',/, &
          '----------------------------------------------------',/)
  end if
  !
  !***  Initializations
  !
  Mcon = 0.0
  Mumb = 0.0
  !
  !***  Writes results
  !
  write(lumas,10) time,time2
10 format(/, &
       '-------------------  ',/,  &
       'TIME INTERVAL (SEC): ',i7,' to ',i7,/, &
       '-------------------  ')
  !
  write(lumas,11)
11 format(/, &
       ' MASS FALLEN: z(km a.s.l.), Mtotal, Mfi for each Phi      ',/)
  !
  write(lumas,'(a)') 'BELOW NBL REGION'
  do is=1,np
     z = Zplum(is)
     Mfi= 0.0
     do ic=1,nc
        Mcon(ic) = Mcon(ic) + Mplum(ic,is)     ! Convective region fallen mass
        Mfi      = Mfi      + Mplum(ic,is)
     end do
     write(lumas,12) z/1d3,Mfi,(Mplum(ic,is),ic=1,nc)
12   format(1x,30(1x,e12.5))
  end do
  !
  write(lumas,'(a)') 'ABOVE NBL REGION'
  do is= np+1,ns
     z = Zplum(is)
     Mfi = 0.0
     do ic = 1,nc
        Mumb(ic) = Mumb(ic) + Mplum(ic,is)        ! Umbrella region fallen mass
        Mfi      = Mfi      + Mplum(ic,is)
     end do
     write(lumas,12) z/1d3,Mfi,(Mplum(ic,is),ic=1,nc)
  end do
  !
  !***  Write mass balance
  !
  Mtcon = 0.0
  Mtumb = 0.0
  do ic = 1,nc
     Mtcon = Mtcon + Mcon(ic)
     Mtumb = Mtumb + Mumb(ic)
  end do
  !
  write(lumas,20) Mtcon,1d2*(Mtcon/(Mtcon+Mtumb))
20 format(/, &
       ' TOTAL MASS FALLEN BELOW NBL ',1x,e12.5,1x,'(',f5.1,'%)')
  !
  write(lumas,21) Mtumb,1d2*(Mtumb/(Mtcon+Mtumb))
21 format( &
       ' TOTAL MASS FALLEN ABOVE NBL ',1x,e12.5,1x,'(',f5.1,'%)')
  !
  write(lumas,22) (Mtcon+Mtumb)
22 format( &
       ' TOTAL MASS FALLEN           ',1x,e12.5)
  !
  write(lumas,30) (-log(diam(ic)*1d3)/log(2.0),ic=1,nc)
30 format(/, &
       ' FI         ',50(1x,f12.2,9x))
  !
  write(lumas,31) (Mcon(ic),1d2*(Mcon(ic)/(Mcon(ic)+Mumb(ic))),ic=1,nc)
31 format( &
       ' BELOW NBL  ',50(1x,e12.5,1x,'(',f5.1,'%)'))
  !
  write(lumas,32) (Mumb(ic),1d2*(Mumb(ic)/(Mcon(ic)+Mumb(ic))),ic=1,nc)
32 format( &
       ' ABOVE NBL  ',50(1x,e12.5,1x,'(',f5.1,'%)'))
  !
  !*** 2. Writes granulometry file
  !
  !  if(type_aggr.eq.'COSTA') then
  write(lumas,40)
40 format(/,&
       ' GRANULOMETRY')
  do ic = 1,npart
     write(lumas,41) 1d3*diam(ic),rhop(ic),sphe(ic),(Mcon(ic)+Mumb(ic))/(Mtcon+Mtumb),TRIM(classname(ic))
41   format(f10.6,1x,f8.1,1x,f7.3,1x,e16.9,2x,a)
  end do
  !  end if
  !
  return
end subroutine wriplumefmas
