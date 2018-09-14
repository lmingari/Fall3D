subroutine cmass3d(c,cload,cdry,cwet,nx,ny,nz,nc,dx,dy,dz,massvol,massdep,massdry,masswet)
  !****************************************************************
  !*
  !*   Computes the mass in the domain and in the deposit
  !*   Note that c and cload hold the scaled concentrations :
  !*       c*     = c     * Hm1
  !*       cload* = cload * Hm1
  !*   However, since dx = Hm1*dX1 = sin(gama)*Rearth*dgama the
  !*   terms cancel in the summ over the area. Consequently, for
  !*   spherical coordinates it is sufficient to call cmass3d with
  !*   dX1 and dX2 (without correcting for Hm1)
  !*
  !****************************************************************
  use KindType
  implicit none
  !
  integer(ip) :: nx,ny,nz,nc
  integer(ip) :: i,j,k,ic
  real   (rp) :: dx,dy,massvol,massdep,massdry,masswet
  real   (rp) :: c(0:nx+1,0:ny+1,0:nz+1,nc),dz(0:nz)
  real   (rp) :: cload(1:nx,1:ny,0:nc),cdry(1:nx,1:ny,0:nc),cwet(1:nx,1:ny,0:nc)
  real   (rp) :: fatx,faty,zmass
  !
  massvol  = 0.0_rp
  massdep  = 0.0_rp
  massdry  = 0.0_rp
  masswet  = 0.0_rp
  !
  do i=1,nx
     fatx=1.0_rp
     if(i == 1.or.i == nx) fatx=0.5_rp
     do j=1,ny
        faty=1.0_rp
        if(j == 1 .or. j == ny) faty=0.5_rp
        !
        !*** volume
        !
        do ic=1,nc
           zmass = 0.0_rp
           do k=1,nz-1
              zmass  = zmass + 0.5_rp*dz(k)*(c(i,j,k,ic)+c(i,j,k+1,ic))
           end do
           massvol = massvol + fatx*faty*zmass
        end do
        !
        !*** deposit (only particles)
        !
        massdep = massdep + fatx*faty*cload(i,j,0)
        massdry = massdry + fatx*faty*cdry (i,j,0)
        masswet = masswet + fatx*faty*cwet (i,j,0)
        !
     end do
  end do
  !
  massvol    = massvol*dx*dy      ! see note above
  massdep    = massdep*dx*dy
  massdry    = massdry*dx*dy
  masswet    = masswet*dx*dy
  !
  return
end subroutine cmass3d
