subroutine writps
  !**********************************************************
  !*
  !*    Write the tracking point files
  !*
  !*********************************************************
  use KindType
  use InpOut
  use Master
  use Numeric
  use Parallel
  implicit none
  !
  character(len=50) :: header1
  integer(ip), save :: ipass = 0
  integer(ip) :: iyr,imo,idy,ihr,imi,ise,i,ic,iz,info
  real   (rp) :: s,t,st,val(4),cload_pts,thick_pts,c_pts,pm5_pts, &
       pm10_pts,pm20_pts,cummc_pts,cp1,cp2
  !
  !***  Calculates time header1 in the form DDMMMM-HH:MM
  !
  call addtime(ibyr,ibmo,ibdy,0,iyr,imo,idy,ihr,imi,ise,time)
  header1 = '00aaa-00:00'
  write(header1(1 :2 ),'(i2.2)') idy
  header1(3 :5 ) = month(imo)
  write(header1(7 :8 ),'(i2.2)') ihr
  write(header1(10:11),'(i2.2)') imi
  !
  !***  Loop over points
  !
  do i = 1,npts
     if(use_pts(i)) then
        !
        !***       Opens file. First time (or not restart) write header
        !
        if(ipass == 0) then
           if(.not.restart) then
              if( mpime == root ) then
                 open(90,file=TRIM(name_file_pts(i))//'.res',iostat=info,status='unknown')
              end if
              call bcast( info, 1, root )
              if( info /= 0 ) goto 100
              !
              if( mpime == root ) write(90,10) TRIM(name_pts(i)),xpts(i),ypts(i)
10            format(&
                   'Tracking point file for : ',a              ,/, &
                   'Coordinates             : ',f13.4,1x,f13.4 ,/, &
                   '  Time    DDMMMM-HH:MM    load     thickness    conc.   conc.PM5   conc.PM10  conc.PM20   cumul',/, &
                   '                         ground   (compact=1)   ground   ground     ground    ground      conc. ',/, &
                   '  (min)      (--)        (kg/m2)     (cm)      (gr/m3)   (gr/m3)    (gr/m3)   (gr/m3)    (gr/m2)',/, &
                   '------------------------------------------------------------------------------------------------')
           else
              if( mpime == root ) then
                 open(90,file=TRIM(name_file_pts(i))//'.res',iostat=info,status='old',position='append')
              end if
              call bcast( info, 1, root )
              if( info /= 0 ) goto 100
           end if
        else
           if( mpime == root ) then
              open(90,file=TRIM(name_file_pts(i))//'.res',iostat=info,status='old',position='append')
           end if
           call bcast( info, 1, root )
           if( info /= 0 ) goto 100
        end if
        !
        !***       Interpolates results
        !
        s = spts(i)
        t = tpts(i)
        st = s*t
        val(1) = (1.0_rp-t-s+st)*0.25_rp                   !  4      3
        val(2) = (1.0_rp-t+s-st)*0.25_rp                   !
        val(3) = (1.0_rp+t+s+st)*0.25_rp                   !
        val(4) = (1.0_rp+t-s-st)*0.25_rp                   !  1      2
        !
        !***       cload (particles only)
        !
        cload_pts =  val(1)*(cload(ipts(i)  ,jpts(i)  ,0)/Hm1(ipts(i)  ,jpts(i)  )) + &
             val(2)*(cload(ipts(i)+1,jpts(i)  ,0)/Hm1(ipts(i)+1,jpts(i)  )) + &
             val(3)*(cload(ipts(i)+1,jpts(i)+1,0)/Hm1(ipts(i)+1,jpts(i)+1)) + &
             val(4)*(cload(ipts(i)  ,jpts(i)+1,0)/Hm1(ipts(i)  ,jpts(i)+1))
        !
        !***       thickness (particles only)
        !
        thick_pts = 0.0_rp
        do ic = 1,np ! particles
           thick_pts =  thick_pts + &
                val(1)*(cload(ipts(i)  ,jpts(i)  ,ic)/Hm1(ipts(i)  ,jpts(i)  )/rhop(ic)) + &
                val(2)*(cload(ipts(i)+1,jpts(i)  ,ic)/Hm1(ipts(i)+1,jpts(i)  )/rhop(ic)) + &
                val(3)*(cload(ipts(i)+1,jpts(i)+1,ic)/Hm1(ipts(i)+1,jpts(i)+1)/rhop(ic)) + &
                val(4)*(cload(ipts(i)  ,jpts(i)+1,ic)/Hm1(ipts(i)  ,jpts(i)+1)/rhop(ic))
        end do
        thick_pts = 1d2*thick_pts  ! cm
        !
        !***       particle concentration ground
        !
        c_pts = 0.0_rp
        do ic = 1,np
           c_pts =  c_pts + &
                val(1)*(c(ipts(i)  ,jpts(i)  ,1,ic)/Hm1(ipts(i)  ,jpts(i)  )) + &
                val(2)*(c(ipts(i)+1,jpts(i)  ,1,ic)/Hm1(ipts(i)+1,jpts(i)  )) + &
                val(3)*(c(ipts(i)+1,jpts(i)+1,1,ic)/Hm1(ipts(i)+1,jpts(i)+1)) + &
                val(4)*(c(ipts(i)  ,jpts(i)+1,1,ic)/Hm1(ipts(i)  ,jpts(i)+1))
        end do
        c_pts = c_pts*1d3    ! gr/m3
        !
        !***       PM5 ground
        !
        pm5_pts = 0.0_rp
        do ic = 1,np  ! particles
           if(diam(ic) <= 5d-6) then
              pm5_pts =  pm5_pts + &
                   val(1)*(c(ipts(i)  ,jpts(i)  ,1,ic)/Hm1(ipts(i)  ,jpts(i)  )) + &
                   val(2)*(c(ipts(i)+1,jpts(i)  ,1,ic)/Hm1(ipts(i)+1,jpts(i)  )) + &
                   val(3)*(c(ipts(i)+1,jpts(i)+1,1,ic)/Hm1(ipts(i)+1,jpts(i)+1)) + &
                   val(4)*(c(ipts(i)  ,jpts(i)+1,1,ic)/Hm1(ipts(i)  ,jpts(i)+1))
           end if
        end do
        pm5_pts = pm5_pts*1d3  ! gr/m3
        !
        !***       PM10 ground
        !
        pm10_pts = 0.0_rp
        do ic = 1,np  ! particles
           if(diam(ic) <= 1d-5) then
              pm10_pts =  pm10_pts + &
                   val(1)*(c(ipts(i)  ,jpts(i)  ,1,ic)/Hm1(ipts(i)  ,jpts(i)  )) + &
                   val(2)*(c(ipts(i)+1,jpts(i)  ,1,ic)/Hm1(ipts(i)+1,jpts(i)  )) + &
                   val(3)*(c(ipts(i)+1,jpts(i)+1,1,ic)/Hm1(ipts(i)+1,jpts(i)+1)) + &
                   val(4)*(c(ipts(i)  ,jpts(i)+1,1,ic)/Hm1(ipts(i)  ,jpts(i)+1))
           end if
        end do
        pm10_pts = pm10_pts*1d3  ! gr/m3
        !
        !***       PM20 ground
        !
        pm20_pts = 0.0_rp
        do ic = 1,np  ! particles
           if(diam(ic) <= 2e-5_rp) then
              pm20_pts =  pm20_pts + &
                   val(1)*(c(ipts(i)  ,jpts(i)  ,1,ic)/Hm1(ipts(i)  ,jpts(i)  )) + &
                   val(2)*(c(ipts(i)+1,jpts(i)  ,1,ic)/Hm1(ipts(i)+1,jpts(i)  )) + &
                   val(3)*(c(ipts(i)+1,jpts(i)+1,1,ic)/Hm1(ipts(i)+1,jpts(i)+1)) + &
                   val(4)*(c(ipts(i)  ,jpts(i)+1,1,ic)/Hm1(ipts(i)  ,jpts(i)+1))
           end if
        end do
        pm20_pts = pm20_pts*1d3  ! gr/m3
        !
        !***       particle cumulative concentration
        !
        cummc_pts = 0.0_rp
        do ic = 1,np
           do iz = 1,nz-1

              cp1 =  val(1)*(c(ipts(i)  ,jpts(i)  ,iz  ,ic)/Hm1(ipts(i)  ,jpts(i)  )) + &
                   val(2)*(c(ipts(i)+1,jpts(i)  ,iz  ,ic)/Hm1(ipts(i)+1,jpts(i)  )) + &
                   val(3)*(c(ipts(i)+1,jpts(i)+1,iz  ,ic)/Hm1(ipts(i)+1,jpts(i)+1)) + &
                   val(4)*(c(ipts(i)  ,jpts(i)+1,iz  ,ic)/Hm1(ipts(i)  ,jpts(i)+1))
              cp2 =  val(1)*(c(ipts(i)  ,jpts(i)  ,iz+1,ic)/Hm1(ipts(i)  ,jpts(i)  )) + &
                   val(2)*(c(ipts(i)+1,jpts(i)  ,iz+1,ic)/Hm1(ipts(i)+1,jpts(i)  )) + &
                   val(3)*(c(ipts(i)+1,jpts(i)+1,iz+1,ic)/Hm1(ipts(i)+1,jpts(i)+1)) + &
                   val(4)*(c(ipts(i)  ,jpts(i)+1,iz+1,ic)/Hm1(ipts(i)  ,jpts(i)+1))
              cummc_pts = cummc_pts + dz(iz)*0.5_rp*(cp1+cp2)
           end do
        end do
        cummc_pts = 1d3*cummc_pts    ! gr/m2
        !
        !***       Writes results
        !
        if( mpime == root ) write(90,20) int(time/60.0_rp),TRIM(header1),cload_pts,thick_pts, &
             c_pts, pm5_pts, pm10_pts, pm20_pts, cummc_pts
20      format(i7,1x,a,100(1x,e13.6))
        !
        !***       Close files
        !
        if( mpime == root ) close(90)
        !
     end if
  end do
  !
  ipass = 1
  return
  !
  !*** List of errors
  !
100 call runend('Can not open file: '//TRIM(name_file_pts(i))//'.res')
  !
end subroutine writps
!
!
!
subroutine writps_end
  !**********************************************************
  !*
  !*    Write the tracking point files
  !*
  !*********************************************************
  use KindType
  use InpOut
  use Master
  use Parallel
  implicit none
  !
  integer(ip) :: i,ic,info
  real   (rp) :: s,t,st,val(4),total
  real   (rp), allocatable :: cload_pts(:)
  !
  !***  Allocates
  !
  allocate(cload_pts(nc))  ! load for each class
  !
  !***  Loop over points
  !
  if( mpime == root ) then
     open(90,file=TRIM(fres)//'.tps',iostat=info,status='unknown')  ! all tps points
  end if
  call bcast( info, 1, root )
  if( info /= 0 ) goto 100
  !
  do i = 1,npts
     if(use_pts(i)) then
        !
        !***       Opens files
        !
        if( mpime == root ) then
           open(91,file=TRIM(name_file_pts(i))//'.final.res',iostat=info,status='unknown')   ! this tps
        end if
        call bcast( info, 1, root )
        if( info /= 0 ) goto 101
        !
        if( mpime == root ) write(91,10) TRIM(name_pts(i)),xpts(i),ypts(i)
10      format('Deposit file for : ',a          ,/, &
             'Coordinates : ',f13.4,1x,f13.4  ,/, &
             ' Diameter    fi       load     fraction ',/, &
             '  (mm)      (--)     (kg/m2)     (%)    ',/, &
             '-------------------------------------------')
        !
        !***       Interpolates results
        !
        s = spts(i)
        t = tpts(i)
        st = s*t
        val(1) = (1.0_rp-t-s+st)*0.25_rp                   !  4      3
        val(2) = (1.0_rp-t+s-st)*0.25_rp                   !
        val(3) = (1.0_rp+t+s+st)*0.25_rp                   !
        val(4) = (1.0_rp+t-s-st)*0.25_rp                   !  1      2
        !
        do ic = 1,np  ! particles
           cload_pts(ic) =  val(1)*(cload(ipts(i)  ,jpts(i)  ,ic)/Hm1(ipts(i)  ,jpts(i)  )) + &
                val(2)*(cload(ipts(i)+1,jpts(i)  ,ic)/Hm1(ipts(i)+1,jpts(i)  )) + &
                val(3)*(cload(ipts(i)+1,jpts(i)+1,ic)/Hm1(ipts(i)+1,jpts(i)+1)) + &
                val(4)*(cload(ipts(i)  ,jpts(i)+1,ic)/Hm1(ipts(i)  ,jpts(i)+1))
        end do
        total = sum(cload_pts)
        !
        !***       Writes results
        !
        if( mpime == root ) write(90,20) xpts(i),ypts(i),total
20      format(f13.4,1x,f13.4,1x,e13.6)

        do ic = 1,np  ! particles
           if( mpime == root ) &
                write(91,30) diam(ic)*1d3, -log(1d3*diam(ic))/log(2.0_rp), &
                cload_pts(ic)  , 1d2*cload_pts(ic)/total
30         format(f7.4,1x,f7.2,1x,e13.6,1x,f9.4)
        end do
        !
        !***       Close files
        !
        if( mpime == root ) close(91)
        !
     end if
  end do
  !
  if( mpime == root ) close(90)
  return
  !
  !*** List of errors
  !
100 call runend('Can not open file: '//TRIM(fres)//'.tps')
101 call runend('Can not open file: '//TRIM(name_file_pts(i))//'.final.res')
  !
end subroutine writps_end
