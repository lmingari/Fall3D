subroutine wriagr
  !**********************************************************
  !*
  !*    Write the aggregation summary
  !*
  !*********************************************************
  use KindType
  use InpOut
  use Master
  use Parallel
  use Domain, only: ic_loc,nc_loc
  implicit none
  integer(ip) :: ic
  real(rp)    :: mtot,mtota
  real(rp), allocatable :: mcl1(:),mcl2(:)
  !
  !***  Collect masses. Each processor has its own mass
  !
  allocate(mcl1(nc))
  allocate(mcl2(nc))
  mcl1(1:nc) = mclass(1:nc,1)              ! erupted mass       (already global)
  mcl2(1:nc) = mclass(1:nc,2)              ! my aggregated mass
  !
  call parallel_sum( mcl2, grp_inter )
  call parallel_sum( mcl2, grp_intra )
  mclass(1:nc,1) = ABS(mcl1(1:nc))
  mclass(1:nc,2) = ABS(mcl2(1:nc))
  !
  !***  Write file
  !
  if(mpime == root) then
     write(nlst,10)
10   format(/,                                                           &
          '----------------------------------------------------',/,    &
          '                                                    ',/,    &
          '             AGGREGATION SUMMARY                    ',/,    &
          '                                                    ',/,    &
          '----------------------------------------------------',/,    &
          '                                                    ',/,    &
          'Class    Fi      Erupted      Aggregated       relat. (%)   total(%)')
     !
     ! Total mass
     !
     mtot  = 0.0_rp
     do ic = 1,np
        mtot = mtot + mclass(ic,1)
     end do
     !
     ! Mass that potentially aggregates
     !
     mtota = 0.0_rp
     do ic = iaggr,np-1
        mtota = mtota + mclass(ic,1)
     end do
     !
     do ic = 1,np
        if(ic /= np) then
           write(nlst,20) ic, -log(1d3*diam(ic))/log(2.0_rp), &
                mclass(ic,1),mclass(ic,2),1d2*abs(mclass(ic,2))/ &
                mclass(ic,1),1d2*abs(mclass(ic,2))/mtot
        else
           write(nlst,21) ic, -log(1d3*diam(ic))/log(2.0_rp), &
                mclass(ic,1),mclass(ic,2),'-','-'
        end if
20      format(1x,i3,2x,f6.1,2x,e12.5,2x,e12.5,2x,f10.2,2x,f10.2)
21      format(1x,i3,2x,f6.1,2x,e12.5,2x,e12.5,2x,a10,2x,a10)
     end do
     !
     write(nlst,30) mtot,mclass(np,2),1d2*abs(mclass(np,2))/mtota,1d2*abs(mclass(np,2))/mtot
30   format('TOTAL ',8x,e12.5,2x,e12.5,2x,f10.2,2x,f10.2)
     !
  end if
  !
  return
end subroutine wriagr
