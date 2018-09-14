  subroutine cmass2d(c,vx,vy,vz,vset,nx,ny,nz,nc, &
     ic_loc,nc_loc,iz_loc,nz_loc,dx,dy,dz,dt,outmass)
  !****************************************************************
  !*
  !*    Computes the mass lost at the boundaries
  !*    Note that c and u hold the scaled quantities :
  !*       c* = c * Hm1
  !*       u* = u / Hm1
  !*   However, since dx = Hm1*dX1 = sin(gama)*Rearth*dgama the
  !*   terms cancel in the summ over the area. Consequently, for
  !*   spherical coordinates it is sufficient to call cmass2d with
  !*   dX1 and dX2 (without correcting for Hm1)
  !*
  !****************************************************************
  use KindType
  use parallel, only:  nproc,mpime,root,parallel_sum,group,bcast
  implicit none
  !
  integer(ip) :: nx,ny,nz,nc,ic_loc,nc_loc,iz_loc,nz_loc
  integer(ip) :: i,j,k,ic
  !
  real   (rp) :: dx,dy,dt,outmass
  real   (rp) :: c(0:nx+1,0:ny+1,0:nz+1,nc),dz(0:nz)
  real   (rp) :: vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz),vset(nx,ny,nz,nc)
  real   (rp) :: fatx,faty,area,velo,conc
  !
  integer(ip), save             :: ipass = 0
  real   (rp), save,allocatable :: mymass(:)   ! mymass(0:nproc-1)
  !
  !***  First time
  !
  if(ipass == 0) then
     ipass = 1
     allocate(mymass(0:nproc-1))
  end if
  call bcast( ipass, 1, root )  ! wait for the rest to allocate
  !
  mymass = 0.0_rp
  !
  !***  Left/Right (X). Sign criteria: influx < 0
  !
  do ic = ic_loc, ic_loc + nc_loc - 1
     do j=1,ny
        faty = 1.0_rp
        if(j == 1.or.j == ny) faty=0.5_rp
        !
        do k = iz_loc, iz_loc + nz_loc - 1
           area = faty*dy*0.5_rp*(dz(k)+dz(k-1))
           !
           velo    = vx(1,j,k   )
           conc    = c (0,j,k,ic)
           mymass(mpime) = mymass(mpime) - dt*area*velo*conc  ! Left. outflux > 0
           !
           velo    = vx(nx  ,j,k   )
           conc    = c (nx+1,j,k,ic)
           mymass(mpime) = mymass(mpime) + dt*area*velo*conc  ! Right. outflux > 0
        enddo
     enddo
  enddo
  !
  !***  Up/Down (Y). Sign criteria: influx < 0
  !
  do ic = ic_loc, ic_loc + nc_loc - 1
     do i=1,nx
        fatx = 1.0_rp
        if(i == 1.or.i == nx) fatx=0.5_rp
        !
        do k = iz_loc, iz_loc + nz_loc - 1
           area    = fatx*dx*0.5_rp*(dz(k)+dz(k-1))
           !
           velo    = vy(i,1,k   )
           conc    = c (i,0,k,ic)
           mymass(mpime) = mymass(mpime) - dt*area*velo*conc  ! Down. outflux > 0
           !
           velo    = vy(i,ny  ,k   )
           conc    = c (i,ny+1,k,ic)
           mymass(mpime) = mymass(mpime) + dt*area*velo*conc  ! Up. outflux > 0
        enddo
     enddo
  enddo
  !
  !***  Top(Z). Bottom is computed as diposit mass. Sign criteria: influx < 0
  !
  do ic = ic_loc, ic_loc + nc_loc - 1
     if( ( iz_loc + nz_loc - 1 ) == nz ) then
        do j=1,ny
           faty = 1.0_rp
           if(j == 1.or.j == ny) faty=0.5_rp
           do i=1,nx
              fatx=1.0_rp
              if(i == 1 .or. i == nx) fatx=0.5_rp
              area = fatx*faty*dx*dy
              !
              velo    = vz(i,j,nz     ) + vset(i,j,nz,ic)
              conc    = c (i,j,nz+1,ic)
              mymass(mpime) = mymass(mpime) + dt*area*velo*conc  ! Top. outflux > 0
           enddo
        enddo
     end if
  enddo
  !
  !***  Communicate
  !
  call parallel_sum(mymass,group)
  outmass = outmass + sum(mymass)
  !
  return
  end subroutine cmass2d
