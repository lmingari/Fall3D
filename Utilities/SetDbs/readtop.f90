subroutine readtop
  !***********************************************************************************
  !*
  !*    Reads topography from an external file
  !*
  !***********************************************************************************
  use KindType, ONLY :  ip,rp
  use Master
  use InpOut
  implicit none
  !
  character(len=4) :: cvoid
  integer  (ip)    :: nxdem,nydem,iy,ix,ixdem,iydem
  real     (rp)    :: xodem,yodem,xfdem,yfdem,dxdem,dydem,x,y,rvoid
  real     (rp), allocatable :: dem(:,:)
  !
  !***  Initialization
  !
  topg = 0_rp
  !
  open(lutop,FILE=TRIM(lutopname),STATUS='old',ERR=101)
  read(lutop,*,err=102) cvoid
  read(lutop,*,err=102) nxdem,nydem
  read(lutop,*,err=102) xodem,xfdem
  read(lutop,*,err=102) yodem,yfdem
  read(lutop,*,err=102) rvoid,rvoid
  !
  !***  Checks that the computational domain is within the regional DEM file
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     !
     if(xodem.gt.lonp(1 ,1 )) goto 103
     if(xfdem.lt.lonp(nx,ny)) goto 104
     if(yodem.gt.latp(1 ,1 )) goto 105
     if(yfdem.lt.latp(nx,ny)) goto 106
     !
  case('UTM')
     !
     if(xodem.gt.xutm(1 ,1 )) goto 203
     if(xfdem.lt.xutm(nx,ny)) goto 204
     if(yodem.gt.yutm(1 ,1 )) goto 205
     if(yfdem.lt.yutm(nx,ny)) goto 206
     !
  END SELECT
  !
  allocate(dem(nxdem,nydem))
  dxdem = (xfdem-xodem)/(nxdem-1)
  dydem = (yfdem-yodem)/(nydem-1)
  do iydem=1,nydem
     read(lutop,*,err=102) (dem(ixdem,iydem),ixdem=1,nxdem)
  end do
  close(lutop)
  !
  !***  Interpolates
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     !
     do iy = 1,ny
        do ix = 1,nx
           x = lonp(ix,iy)
           y = latp(ix,iy)
           topg(ix,iy) = get_val(xodem,yodem,nxdem,nydem,dxdem,dydem,dem,x,y)
        end do
     end do
     !
  case('UTM')
     !
     do iy = 1,ny
        do ix = 1,nx
           x = xutm(ix,iy)
           y = yutm(ix,iy)
           topg(ix,iy) = get_val(xodem,yodem,nxdem,nydem,dxdem,dydem,dem,x,y)
        end do
     end do
     !
  END SELECT
  !
  !***  Releases memory
  !
  deallocate(dem)
  return
  !
  !***   List of errors
  !
101 call wriwar('Unable to open file'//TRIM(lutopname)//'. Assuming no topography')
  return
102 call wriwar('Error reading file'//TRIM(lutopname)//'. Assuming no topography')
  return
103 call wriwar('lonmin of the domain is outside the DEM file. Assuming no topography')
  return
104 call wriwar('lonmax of the domain is outside the DEM file. Assuming no topography')
  return
105 call wriwar('latmin of the domain is outside the DEM file. Assuming no topography')
  return
106 call wriwar('latmax of the domain is outside the DEM file. Assuming no topography')
  return
203 call wriwar('xmin of the domain is outside the DEM file. Assuming no topography')
  return
204 call wriwar('xmax of the domain is outside the DEM file. Assuming no topography')
  return
205 call wriwar('ymin of the domain is outside the DEM file. Assuming no topography')
  return
206 call wriwar('ymax of the domain is outside the DEM file. Assuming no topography')
  return
  !
end subroutine readtop
