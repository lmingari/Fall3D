  subroutine sizloop
  !*******************************************************************
  !*
  !*    Loop on particles classes
  !*
  !*    c        - Concentration (Input/Output matrix)
  !*    cload    - Ground loading (Output)
  !*    vx       - X-component of the wind direction (Input Matrix)
  !*    vy       - Y-component of the wind direction (Input Matrix)
  !*    vz       - Z-component of the wind direction (Input Matrix)
  !*    dz       - Dimensions of the cells in Z direction (Input Vector)
  !*    rkhor    - Horizontal diffusion coefficients (Input)
  !*    rkver    - Vertical diffusion coefficient (Input)
  !*    rho      - Atmospheric density (Input Matrix)
  !*    vset     - Particle settling velocity (Work Matrix)
  !*    work     - Auxiliary vector (Input)
  !*
  !**********************************************************************
  use KindType
  use Numeric
  use Master
  use Domain,   only: ic_loc, nc_loc, nz_loc, iz_loc, &
                      collect_cz, exchange_1layer, exchange_2layer
  use Parallel, only: mpime
  implicit none
  !
  integer(ip) :: ic
  !
  !***  Loop over classes
  !
  do ic = ic_loc, ic_loc + nc_loc - 1
     !
     !*** Set Boundary conditions
     !
     call setbcc(c(0,0,0,ic),nx,ny,nz,nz_loc,iz_loc,vx,vy,vz,vset(1,1,1,ic),vdep(1,1,ic),rkver,dz)
     !
     !*** Advection
     !
     call exchange_2layer( c( :, :, :, ic ), nz )  !  you must communicate boundary cells
     call advctzc(c(0,0,0,ic),nx,ny,nz,nz_loc,iz_loc,vz,vset(1,1,1,ic),dz,dt,lmeth)
     !
     if(mod(iiter,2) == 0) then
        call advctx(c(0,0,0,ic),nx,ny,nz,nz_loc,iz_loc,vx,dX1,dt,work,lmeth)
        call advcty(c(0,0,0,ic),nx,ny,nz,nz_loc,iz_loc,vy,dX2,dt,work,lmeth)
     else
        call advcty(c(0,0,0,ic),nx,ny,nz,nz_loc,iz_loc,vy,dX2,dt,work,lmeth)
        call advctx(c(0,0,0,ic),nx,ny,nz,nz_loc,iz_loc,vx,dX1,dt,work,lmeth)
     end if
     !
     !*** Diffusion (scaling breaks symmetry for rhhor!)
     !
     call exchange_1layer( c( :, :, :, ic ), nz )  !  you must communicate boundary cells
     call diffz(c(0,0,0,ic),nx,ny,nz,nz_loc,iz_loc,dz,dt,rkver,rho,work)
     !
     if(mod(iiter,2) == 0) then
        call diffx(c(0,0,0,ic),nx,ny,nz,nz_loc,iz_loc,dX1,dt,rkh1 ,rho,work)
        call diffy(c(0,0,0,ic),nx,ny,nz,nz_loc,iz_loc,dX2,dt,rkhor,rho,work)
     else
        call diffy(c(0,0,0,ic),nx,ny,nz,nz_loc,iz_loc,dX2,dt,rkhor,rho,work)
        call diffx(c(0,0,0,ic),nx,ny,nz,nz_loc,iz_loc,dX1,dt,rkh1 ,rho,work)
     end if
     !
     !*** Wet deposition (computed for particles only, not for gas species)
     !
     if( (wetdeposition).and.(ic.le.np) ) then
       call wetdep(c(0,0,0,ic),cwet,ic,nx,ny,nz,nz_loc,iz_loc,dt,diam,zlayer,dz,zi,prate)
     end if
     !
     !*** Accumulation from dry deposition and sedimentation
     !
     call accum(c(0,0,0,ic),cdry,ic,nx,ny,nz,nz_loc,iz_loc,vz,dt,vset(1,1,1,ic),vdep(1,1,ic))
     !
  end do
  !
  !***  Clip negative c values
  !
  c(0:nx+1,0:ny+1,0:nz+1,ic_loc:(ic_loc + nc_loc - 1)) = MAX(0.0_rp,  &
  c(0:nx+1,0:ny+1,0:nz+1,ic_loc:(ic_loc + nc_loc - 1)))
  !
  !***  Computes mass flowing out during this time interval
  !
  call cmass2d(c,vx,vy,vz,vset,nx,ny,nz,nc, &
       ic_loc,nc_loc,iz_loc,nz_loc,dX1,dX2,dz,dt,outmass)
  !
  return
  end subroutine sizloop
