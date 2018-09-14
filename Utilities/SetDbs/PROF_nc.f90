!***************************************************************
!*
!*    Module for PROF operations
!*
!***************************************************************
     MODULE PROF_nc
     use KindType, ONLY :  ip,rp
     implicit none
     save
!
!*** DEM variables (GRD file)
!
     character(len=4) :: cvoid
     integer(ip) :: nxdem,nydem,ixdem,iydem
     real(rp) :: xodem,yodem,xfdem,yfdem,dxdem,dydem,x,y,rvoid
     real(rp), allocatable :: dem(:,:)
!
!*** PROF variables (profile file)
!
     integer(ip) :: ibyr_PROF,ibmo_PROF,ibdy_PROF,nt_PROF,it_PROF,nz_PROF,iz_PROF
     real(rp) :: xo_PROF,yo_PROF
!
     integer(ip), parameter :: nt_PROF_max = 100
     real(rp) :: time_PROF(nt_PROF_max)
     real(rp) :: timesec_PROF(nt_PROF_max)
!
     END MODULE PROF_nc
