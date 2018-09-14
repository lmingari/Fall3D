  subroutine setsrc
  !************************************************************************
  !*
  !*    Sets the source term (already scaled)
  !*
  !************************************************************************
  use KindType
  use Numeric
  use Master
  use Domain, only: ic_loc, nc_loc, iz_loc, nz_loc
  use Parallel
  implicit none
  !
  integer(ip), save :: ipass = 0
  integer(ip) :: is,i,j,k,ic,isou
  real   (rp) :: dzeta,vol0
  !
  !***  Source term (as in the file)
  !
  do is=1,nsou
     i=isrc(is)
     j=jsrc(is)
     k=ksrc(is)
     if( k >= iz_loc .and. k <= ( iz_loc + nz_loc -1 ) ) then
        if(k==1) then
          dzeta = 0.5_rp*dz(k)+zlayer(k)
        else
          dzeta = 0.5_rp*(dz(k-1)+dz(k))
        endif
        vol0 = Hm1(i,j)*dX1*dX2*dzeta
        do ic = 1,nc
           c(i,j,k,ic) = c(i,j,k,ic) + dt*tmrat(is,ic)/vol0
        end do
     end if
  end do
  !
  !*** Account for aggregation. This is done for particles only.
  !
  if(aggregation) then
     !
     ipass = ipass + 1
     !
     !***  Compute tarat every ndtagr time steps only (avoid a collect_all call
     !***  every time step). This reduces the computational time
     !
     if (mod(ipass,ndtagr) == 0) then
        !
        call agrsrc
        !
        do ic = ic_loc, ic_loc + nc_loc - 1
           do k = iz_loc,iz_loc+nz_loc-1
              dzeta= 0.5_rp*(dz(k-1)+dz(k))
              do j = 1,ny
                 do i = 1,nx
                    c(i,j,k,ic) = c(i,j,k,ic) + ndtagr*dt*tarat(i,j,k,ic)
                    ! without /vol0 (rate is already per unit volume)
                 end do
              end do
           end do
        end do
        !
        !***  Accumulated mass of aggregates.
        !***  This is stored to visualize where the aggregates form
        !
        do k = 1, nz
           dzeta= 0.5_rp*(dz(k-1)+dz(k))
           if( ( k < iz_loc ) .or. ( k > (iz_loc + nz_loc - 1) ) ) then
              Maggr( :, :, k ) = 0.0_rp
           else
              do j = 1,ny
                 do i = 1,nx
                    Maggr(i,j,k) = Maggr(i,j,k) + ndtagr*dt*tarat(i,j,k,np)*Hm1(i,j)*dX1*dX2*dzeta
                 end do
              end do
           end if
        end do
        !
        !***  Class mass balance
        !***  mclass(ic,1) stores the source mass.
        !***  This is global because tmrat is global
        !
        do isou = 1,nsou
           do ic   = 1,nc
              mclass(ic,1) = mclass(ic,1) + ndtagr*dt*tmrat(isou,ic)/Hm1(isrc(isou),jsrc(isou)) ! source mass
           end do
        end do
        !
        !***  Class mass balance
        !***  mclass(ic,2) stores the mass of aggregates. This is local because tarat is local
        !
        do ic = ic_loc, ic_loc + nc_loc - 1
           do k = iz_loc,iz_loc+nz_loc-1
              dzeta= 0.5_rp*(dz(k-1)+dz(k))
              do j = 1,ny
                 do i = 1,nx
                    !
                    mclass(ic,2) = mclass(ic,2) + ndtagr*dt*tarat(i,j,k,ic)*Hm1(i,j)*dX1*dX2*dzeta ! Mass of aggragates.
                    !
                 end do
              end do
           end do
        end do
        !
     end if
  end if
  !
  return
  end subroutine setsrc
