subroutine read_PROF_grid(fname,nx,ny,lonp,latp,xutm,yutm,topg)
  !********************************************************************
  !*
  !*   Reads topography from a GRD file
  !*
  !********************************************************************
  use KindType, ONLY : ip,rp
  use InpOut,   ONLY : lutop
  use Master,   ONLY : coord_sys,xmin,xmax,ymin,ymax,lonmin,lonmax,latmin,latmax
  use PROF_nc
  implicit none
  !
  character(len=*) :: fname
  integer(ip) :: nx,ny,ix,iy
  real(rp) :: lonp(nx,ny),latp(nx,ny),xutm(nx,ny),yutm(nx,ny)
  real(rp) :: topg(nx,ny)
  !
  !*** Reads topography from a regional DEM file and
  !
  open(lutop,FILE=TRIM(fname),STATUS='old',ERR=101)
  read(lutop,*,err=102) cvoid
  read(lutop,*,err=102) nxdem,nydem
  read(lutop,*,err=102) xodem,xfdem
  read(lutop,*,err=102) yodem,yfdem
  read(lutop,*,err=102) rvoid,rvoid
  !
  !*** Checks that the computational domain lays within the regional DEM file
  !
  SELECT CASE(coord_sys)
  CASE('UTM')
     if(xodem.gt.xmin) call runend('read_PRO_grid: xmin of the domain is outside the DEM file')
     if(xfdem.lt.xmax) call runend('read_PRO_grid: xmax of the domain is outside the DEM file')
     if(yodem.gt.ymin) call runend('read_PRO_grid: ymin of the domain is outside the DEM file')
     if(yfdem.lt.ymax) call runend('read_PRO_grid: ymax of the domain is outside the DEM file')
  CASE('LON-LAT')
     if(xodem.gt.lonmin) call runend('read_PRO_grid: lonmin of the domain is outside the DEM file')
     if(xfdem.lt.lonmax) call runend('read_PRO_grid: lonmax of the domain is outside the DEM file')
     if(yodem.gt.latmin) call runend('read_PRO_grid: latmin of the domain is outside the DEM file')
     if(yfdem.lt.latmax) call runend('read_PRO_grid: latmax of the domain is outside the DEM file')
  END SELECT
  !
  !*** Reads
  !
  allocate(dem(nxdem,nydem))
  dxdem = (xfdem-xodem)/(nxdem-1)
  dydem = (yfdem-yodem)/(nydem-1)
  do iydem=1,nydem
     read(lutop,*,err=102) (dem(ixdem,iydem),ixdem=1,nxdem)
  end do
  close(lutop)
  !
  !*** Interpolates to get topography
  !
  SELECT CASE(coord_sys)
  CASE('LON-LAT')
     do iy = 1,ny
        do ix = 1,nx
           call get_topo_val(xodem,yodem,nxdem,nydem,dxdem,dydem,dem,lonp(ix,iy),latp(ix,iy),topg(ix,iy))
        end do
     end do
  CASE('UTM')
     do iy = 1,ny
        do ix = 1,nx
           call get_topo_val(xodem,yodem,nxdem,nydem,dxdem,dydem,dem,xutm(ix,iy),yutm(ix,iy),topg(ix,iy))
        end do
     end do
  END SELECT
  !
  !*** Releases memory
  !
  deallocate(dem)
  !
  return
  !
  !*** List of errors
  !
101 call runend('Error opening topography file'//TRIM(fname))
102 call runend('Error reading topography file'//TRIM(fname))
  !
  return
end subroutine read_PROF_grid
!
!
!
subroutine get_topo_val(xodem,yodem,nxdem,nydem,dxdem,dydem,dem,x,y,topg)
  !***********************************************************
  !*
  !*   Interpolates dem at the (x,y) point
  !*
  !***********************************************************
  use KindType, ONLY : ip,rp
  implicit none
  !
  integer(ip) :: nxdem,nydem
  real   (rp) :: xodem,yodem,dxdem,dydem,dem(nxdem,nydem)
  real   (rp) :: x,y,xdem,ydem,topg
  !
  logical         :: xfound,yfound
  integer(ip) :: iydem,ixdem,ixpt,iypt
  real   (ip) :: s,t,st,shape(4)
  !
  xfound = .false.
  yfound = .false.
  !
  do ixdem = 1,nxdem-1
     xdem = xodem + (ixdem-1)*dxdem
     if((x.ge.xdem).and.(x.le.(xdem+dxdem))) then
        ixpt   = ixdem
        xfound = .true.
     end if
  end do
  if(.not.xfound) call runend('get_topo_val: value not found in x')
  !
  do iydem = 1,nydem-1
     ydem = yodem + (iydem-1)*dydem
     if((y.ge.ydem).and.(y.le.(ydem+dydem))) then
        iypt   = iydem
        yfound = .true.
     end if
  end do
  if(.not.yfound) call runend('get_topo_val: value not found in y')
  !
  !*** Interpolates
  !
  s = (x-(xodem+(ixpt-1)*dxdem))/dxdem      ! parameter s in (0,1)
  s = 2d0*s-1d0                               ! parameter s in (-1,1)
  t = (y-(yodem+(iypt-1)*dydem))/dydem      ! parameter t in (0,1)
  t = 2d0*t-1d0                               ! parameter t in (-1,1)
  st = s*t
  !
  shape(1)=(1d0-t-s+st)*0.25d0                      ! 4         3
  shape(2)=(1d0-t+s-st)*0.25d0                      !
  shape(3)=(1d0+t+s+st)*0.25d0                      !
  shape(4)=(1d0+t-s-st)*0.25d0                      ! 1         2
  !
  topg = shape(1)*dem(ixpt  ,iypt  ) +  &
       shape(2)*dem(ixpt+1,iypt  ) +  &
       shape(3)*dem(ixpt+1,iypt+1) +  &
       shape(4)*dem(ixpt  ,iypt+1)
  !
  return
end subroutine get_topo_val
