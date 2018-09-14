  subroutine divcorr
  !************************************************************************
  !*
  !*    Correct the mass unbalance due to the nonnull divergence
  !*
  !************************************************************************
  use KindType
  use Numeric
  use Master
  use Domain,  only: ic_loc, nc_loc, nz_loc, iz_loc
  use Parallel
  implicit none
  !
  integer(ip) :: i,j,k,ic
  real   (rp) :: dive, corrdiv, eps
  real   (rp) :: compmass, eps1
  !
  call bcast( outmass, 1, root )
  call bcast( massvol, 1, root )
  call bcast( massdep, 1, root )
  !
  compmass = massvol+massdep+outmass
  eps  = 0.005_rp
  eps1= 10.0_rp*eps
  corrdiv = max(min((erumass-compmass)/(erumass),eps1),-eps1)
  !
  do k = iz_loc, iz_loc + nz_loc - 1
     do j=1,ny
        do i=1,nx
           dive = divu(i,j,k)
           if(abs(dive*dt).ge.eps) then
              do ic = ic_loc, ic_loc + nc_loc - 1
                 c(i,j,k,ic) = c(i,j,k,ic)*(1.0_rp+corrdiv)
              end do
            endif
        end do
     end do
  end do
  !
  return
  end subroutine divcorr
