  subroutine finddt
  !************************************************************************
  !*
  !*    Finds the time step as  dt = safety_factor*dtc , where dtc is the critical
  !*    time step
  !*
  !*    dz       - Dimensions of the cells in the Z direction (Input Vector)
  !*    vx,vy,vz - Wind field (Input Matrix)
  !*    rho      - Atmospheric density (Input Matrix)
  !*    rkhor    - Horizontal diffusion coefficients (Input Matrix)
  !*    rkver    - Vertical diffusion coefficients (Input Matrix)
  !*    rjac     - Jacobian (Input Matrix)
  !*    ztop     - Maximum height of the layers (Input)
  !*    z_map      - Apply Z-mapping: 0=NO; 1=YES (Input)
  !*
  !***********************************************************************
  use KindType
  use Numeric
  use Master
  implicit none
  !
  integer(ip) :: ic,k,i,j
  real   (rp) :: dtcn,ddz,dtcni,dvsetdz
  real   (rp) :: R8MIN = 1d-13
  !
  dtcn = 0.0_rp
  !
  do ic=1,nc
     do k=1,nz
        ddz = min(dz(k-1),dz(k))
        do j=1,ny
           do i=1,nx
              !
              if(k == 1) then
                 dvsetdz = (vset(i,j,2,ic)-vset(i,j,1,ic))/(dz(1))
              else if(k == nz) then
                 dvsetdz = (vset(i,j,nz,ic)-vset(i,j,nz-1,ic))/(dz(nz-1))
              else
                 dvsetdz = (vset(i,j,k+1,ic)-vset(i,j,k-1,ic))/(dz(k)+dz(k-1))
              end if
              !
              dtcni= (2.0_rp*rkh1 (i,j,k)/(dX1**2.0_rp))+(abs(vx(i,j,k))/dX1)   &
                   +(2.0_rp*rkhor(i,j,k)/(dX2**2.0_rp))+(abs(vy(i,j,k))/dX2)   &
                   +(2.0_rp*rkver(i,j,k)/(ddz**2.0_rp))     &
                   +(abs((vz(i,j,k)+vset(i,j,k,ic)))/ddz)    &
                   +(abs(dvsetdz))
              dtcn=max(dtcn,dtcni)
           enddo
        enddo
     enddo
  enddo
  !
  dtcn=max(dtcn,R8MIN)
  !
  !***  Critical time step
  !
  dtc = 1.0_rp/dtcn
  !
  !***  Time setp
  !
  dt = safety_factor*dtc
  !
  return
  end subroutine finddt
