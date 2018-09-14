   subroutine printres_nc
  !*********************************************************************
  !*
  !*    Writes results for postprocess in netCDF format
  !*
  !*********************************************************************
  use KindType
  use InpOut
  use Master
  use Numeric
  use netcdf
  use Parallel
  implicit none
  !
  character(len=18)     :: header
  integer(ip)           :: iyr,imo,idy,ihr,imi,ise,ic,iz,istat,iflevel,ieff,i,k
  integer(ip), save     :: ipass = 0
  real   (rp)           :: flm,outtime
  real   (rp)           :: reff(20),Qeff(20)
  !
  real(rp), allocatable :: work2d(:,:)
  real(rp), allocatable :: work3d(:,:,:)
  !
  !***  First time writes headers and topography
  !
  if(ipass == 0) call print_nc_grid
  ipass = ipass + 1
  !
  !***  Calculates time header in the form YYYY:MM:DD:HH:SSSS
  !***  It goes to the current time step of the netCDF variable time
  !
  outtime = beg_time + (ippfile-1)*print_time
  call addtime(ibyr,ibmo,ibdy,0,iyr,imo,idy,ihr,imi,ise,outtime)
  header = '0000:00:00:00:0000'
  write(header(1 :4 ),'(i4.4)') iyr
  write(header(6 :7 ),'(i2.2)') imo
  write(header(9 :10),'(i2.2)') idy
  write(header(12:13),'(i2.2)') ihr
  write(header(15:18),'(i4.4)') (imi*60)+ise
  !
  !***  Allocates memory
  !
  allocate(work2d(nx,ny))
  allocate(work3d(nx,ny,nz))
  !
  !***  Start to write results for the current time step
  !***  Open netCDF file and get ncID
  !
  if( mpime == root ) then
     istat = nf90_open(TRIM(fres_nc),NF90_WRITE, ncID)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('printres_nc : Error in nf90_open')
  !
  !***  First time writes topography
  !
  if(ipass == 1) then
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,topog_nc_name,topog_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : Error getting topog_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,topog_nc_ID, topg, start=(/1,1/), count=(/nx,ny/) )
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable topg')
  end if
  !
  !***  Time variable (current step)
  !
  outtime = outtime/3600.0_rp
  if( mpime == root ) then
     istat = nf90_inq_varid(ncID,time_nc_name,time_nc_ID)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('printres_nc : Error getting time_nc_ID')
  !
  if( mpime == root ) then
     istat = nf90_put_var(ncID,time_nc_ID,outtime,start=(/ipass/))
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable time')
  !
  if( mpime == root ) then
     istat = nf90_inq_varid(ncID,timeh_nc_name,timeh_nc_ID)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('printres_nc : Error getting timeh_nc_ID')
  !
  if( mpime == root ) then
     istat = nf90_put_var(ncID,timeh_nc_ID, header, start=(/ipass/), count=(/1/) )
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable timeh')
  !
  !***  Total particle deposit load (kg/m2)
  !
  if( out_load ) then
     work2d(1:nx,1:ny) = cload(1:nx,1:ny,0)/Hm1(1:nx,1:ny)
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,cload_nc_name,cload_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : Error getting cload_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,cload_nc_ID, work2d, start=(/1,1,ipass/), count=(/nx,ny,1/) )
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable cload')
  end if
  !
  !***  Particle class deposit load (in kg/m2)
  !
  if( out_load ) then
     if(print_class) then  ! load of particles
        do ic = 1,np
           work2d(1:nx,1:ny) = cload(1:nx,1:ny,ic)/Hm1(1:nx,1:ny)
           if( mpime == root ) then
              istat =  nf90_inq_varid(ncID,cload_ic_name(ic),cload_ic_ID(ic))
              if(istat /= 0) call wriwar(nf90_strerror(istat))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) call runend('printres_nc : Error getting cload_ic_ID')

           if( mpime == root ) then
              istat =  nf90_put_var(ncID,cload_ic_ID(ic), work2d, start=(/1,1,ipass/), count=(/nx,ny,1/) )
              if(istat /= 0) call wriwar(nf90_strerror(istat))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable cload')
        end do
     end if
  end if
  !
  !***  Total particle wet deposition deposit load (kg/m2)
  !
  if( out_load .and. wetdeposition ) then
     work2d(1:nx,1:ny) = cwet(1:nx,1:ny,0)/Hm1(1:nx,1:ny)
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,cwet_nc_name,cwet_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : Error getting cwet_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,cwet_nc_ID, work2d, start=(/1,1,ipass/), count=(/nx,ny,1/) )
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable cwet')
  end if
  !
  !***  Wet deposition particle class deposit load (in kg/m2)
  !
  if( out_load .and. wetdeposition ) then
     if(print_class) then  ! load of particles
        do ic = 1,np
           work2d(1:nx,1:ny) = cwet(1:nx,1:ny,ic)/Hm1(1:nx,1:ny)
           if( mpime == root ) then
              istat =  nf90_inq_varid(ncID,cwet_ic_name(ic),cwet_ic_ID(ic))
              if(istat /= 0) call wriwar(nf90_strerror(istat))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) call runend('printres_nc : Error getting cwet_ic_ID')

           if( mpime == root ) then
              istat =  nf90_put_var(ncID,cwet_ic_ID(ic), work2d, start=(/1,1,ipass/), count=(/nx,ny,1/) )
              if(istat /= 0) call wriwar(nf90_strerror(istat))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable cwet')
        end do
     end if
  end if
  !
  !***  Writes ground deposit thickness (in cm). Like for cload, thickness
  !***  is computed for particles only.
  !
  if( out_thick ) then
     work2d = 0.0_rp
     do ic=1,np    ! particles
        work2d(1:nx,1:ny) = work2d(1:nx,1:ny) + cload(1:nx,1:ny,ic)/rhop(ic)
     end do
     work2d(1:nx,1:ny) = 1d2*work2d(1:nx,1:ny)/Hm1(1:nx,1:ny)
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,thick_nc_name,thick_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : Error getting thick_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,thick_nc_ID, work2d, start=(/1,1,ipass/), count=(/nx,ny,1/) )
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable thick')
  end if
  !
  !***  Calculates and writes z-cumulative concentration for particles (in gr/m2)
  !
  if( out_cummc ) then
     work3d = 0.0_rp
     do ic = 1,np
        work3d = work3d + c(1:nx,1:ny,1:nz,ic)
     end do
     do iz = 1,nz
        work3d(1:nx,1:ny,iz) = work3d(1:nx,1:ny,iz)/Hm1(1:nx,1:ny)  ! concentration
     end do
     !
     work2d = 0.0_rp
     do iz = 1,nz-1
        work2d(1:nx,1:ny) =  work2d(1:nx,1:ny) + dz(iz)*0.5_rp*  &
             (work3d(1:nx,1:ny,iz) + work3d(1:nx,1:ny,iz+1))
     end do
     work2d(1:nx,1:ny) =  1d3*work2d(1:nx,1:ny)  ! gr/m2
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,cummc_nc_name,cummc_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 )  call runend('printres_nc : Error getting cummc_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,cummc_nc_ID, work2d, start=(/1,1,ipass/), count=(/nx,ny,1/))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable cummc')
  end if
  !
  !***  Calculates and writes z-cumulative concentration for each gas (in gr/m2)
  !
  if( out_cummc ) then
     do ic = 1,ng          ! gas
        work3d = c(1:nx,1:ny,1:nz,np+ic)
        do iz = 1,nz
           work3d(1:nx,1:ny,iz) = work3d(1:nx,1:ny,iz)/Hm1(1:nx,1:ny)  ! concentration
        end do
        !
        work2d = 0.0_rp
        do iz = 1,nz-1
           work2d(1:nx,1:ny) =  work2d(1:nx,1:ny) + dz(iz)*0.5_rp*  &
                (work3d(1:nx,1:ny,iz) + work3d(1:nx,1:ny,iz+1))
        end do
        work2d(1:nx,1:ny) =  1d3*work2d(1:nx,1:ny)  ! gr/m2
        !
        if( mpime == root ) then
           istat = nf90_inq_varid(ncID,cummc_ic_name(ic),cummc_ic_ID(ic))
           if(istat /= 0) call wriwar(nf90_strerror(istat))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 )  call runend('printres_nc : Error getting cummc_ic_ID')
        !
        if( mpime == root ) then
           istat = nf90_put_var(ncID,cummc_ic_ID(ic), work2d, start=(/1,1,ipass/), count=(/nx,ny,1/))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable cummc_ic')
     end do
  end if
  !
  !***  Calculates and writes PM05 z-cumulative concentration (in gr/m2)
  !***  Again, this refers to particulate matter only
  !
  if(  out_pm05c ) then
     work3d = 0.0_rp
     do ic = 1,np   ! particles
        if(diam(ic) <= 5d-6) then
           work3d = work3d + c(1:nx,1:ny,1:nz,ic)
        end if
     end do
     do iz = 1,nz
        work3d(1:nx,1:ny,iz) = work3d(1:nx,1:ny,iz)/Hm1(1:nx,1:ny)  ! concentration
     end do
     !
     work2d = 0.0_rp
     do iz = 1,nz-1
        work2d(1:nx,1:ny) =  work2d(1:nx,1:ny) + dz(iz)*0.5_rp*  &
             (work3d(1:nx,1:ny,iz) + work3d(1:nx,1:ny,iz+1))
     end do
     work2d(1:nx,1:ny) =  1d3*work2d(1:nx,1:ny)  ! gr/m2
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,pm05c_nc_name,pm05c_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 )  call runend('printres_nc : Error getting pm05c_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,pm05c_nc_ID, work2d, start=(/1,1,ipass/), count=(/nx,ny,1/))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable pm05c')
  end if
  !
  !***  Calculates and writes PM10 z-cumulative concentration (in gr/m2)
  !***  Again, this refers to particulate matter only
  !
  if( out_pm10c ) then
     work3d = 0.0_rp
     do ic = 1,np   ! particles
        if(diam(ic) <= 1d-5) then
           work3d = work3d + c(1:nx,1:ny,1:nz,ic)
        end if
     end do
     do iz = 1,nz
        work3d(1:nx,1:ny,iz) = work3d(1:nx,1:ny,iz)/Hm1(1:nx,1:ny)  ! concentration
     end do
     !
     work2d = 0.0_rp
     do iz = 1,nz-1
        work2d(1:nx,1:ny) =  work2d(1:nx,1:ny) + dz(iz)*0.5_rp*  &
             (work3d(1:nx,1:ny,iz) + work3d(1:nx,1:ny,iz+1))
     end do
     work2d(1:nx,1:ny) =  1d3*work2d(1:nx,1:ny)  ! gr/m2
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,pm10c_nc_name,pm10c_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 )  call runend('printres_nc : Error getting pm10c_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,pm10c_nc_ID, work2d, start=(/1,1,ipass/), count=(/nx,ny,1/))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable pm10c')
  end if
  !
  !***  Calculates and writes PM20 z-cumulative concentration (in gr/m2)
  !***  Again, this refers to particulate matter only
  !
  if( out_pm20c ) then
     work3d = 0.0_rp
     do ic = 1,np  ! particles
        if(diam(ic) <= 2d-5) then
           work3d = work3d + c(1:nx,1:ny,1:nz,ic)
        end if
     end do
     do iz = 1,nz
        work3d(1:nx,1:ny,iz) = work3d(1:nx,1:ny,iz)/Hm1(1:nx,1:ny)  ! concentration
     end do
     !
     work2d = 0.0_rp
     do iz = 1,nz-1
        work2d(1:nx,1:ny) =  work2d(1:nx,1:ny) + dz(iz)*0.5_rp*  &
             (work3d(1:nx,1:ny,iz) + work3d(1:nx,1:ny,iz+1))
     end do
     work2d(1:nx,1:ny) =  1d3*work2d(1:nx,1:ny)  ! gr/m2
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,pm20c_nc_name,pm20c_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 )  call runend('printres_nc : Error getting pm20c_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,pm20c_nc_ID, work2d, start=(/1,1,ipass/), count=(/nx,ny,1/))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable pm20c')
  end if
  !
  !***  Calculates and writes total particle concentration at ground (gr/m3)
  !
  if( out_concg ) then
     work2d = 0.0_rp
     do ic=1,np
        work2d(1:nx,1:ny) = work2d(1:nx,1:ny) + c(1:nx,1:ny,1,ic)
     end do
     work2d(1:nx,1:ny) = 1d3*work2d(1:nx,1:ny)/Hm1(1:nx,1:ny)            ! gr/m3
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,concg_nc_name,concg_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 )  call runend('printres_nc : Error getting concg_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,concg_nc_ID, work2d, start=(/1,1,ipass/), count=(/nx,ny,1/))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable concg')
  end if
  !
  !***  Calculates and writes PM05 concentration at the first layer (gr/m3)
  !***  Again, this refers to particulate matter only
  !
  if( out_pm05g ) then
     work2d = 0.0_rp
     do ic=1,np ! particles
        if(diam(ic) <= 5d-6) then                                        ! 05 Microns
           work2d(1:nx,1:ny) = work2d(1:nx,1:ny) + c(1:nx,1:ny,1,ic)
        end if
     end do
     work2d(1:nx,1:ny) = 1d3*work2d(1:nx,1:ny)/Hm1(1:nx,1:ny)            ! gr/m3
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,pm05g_nc_name,pm05g_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 )  call runend('printres_nc : Error getting pm05g_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,pm05g_nc_ID, work2d, start=(/1,1,ipass/), count=(/nx,ny,1/))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable pm05g')
  end if
  !
  !***  Calculates and writes PM10 concentration at the first layer (gr/m3)
  !***  Again, this refers to particulate matter only
  !
  if( out_pm10g ) then
     work2d = 0.0_rp
     do ic=1,np  ! particles
        if(diam(ic) <= 1d-5) then                                        ! 10 Microns
           work2d(1:nx,1:ny) = work2d(1:nx,1:ny) + c(1:nx,1:ny,1,ic)
        end if
     end do

     work2d(1:nx,1:ny) = 1d3*work2d(1:nx,1:ny)/Hm1(1:nx,1:ny)            ! gr/m3
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,pm10g_nc_name,pm10g_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 )  call runend('printres_nc : Error getting pm10g_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,pm10g_nc_ID, work2d, start=(/1,1,ipass/), count=(/nx,ny,1/))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable pm10g')
  end if
  !
  !***  Calculates and writes PM20 concentration at the first layer (gr/m3)
  !***  Again, this refers to particulate matter only
  !
  if( out_pm20g ) then
     work2d = 0.0_rp
     do ic=1,np  ! particles
        if(diam(ic) <= 2d-5) then                                        ! 20 Microns
           work2d(1:nx,1:ny) = work2d(1:nx,1:ny) + c(1:nx,1:ny,1,ic)
        end if
     end do
     work2d(1:nx,1:ny) = 1d3*work2d(1:nx,1:ny)/Hm1(1:nx,1:ny)            ! gr/m3
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,pm20g_nc_name,pm20g_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 )  call runend('printres_nc : Error getting pm20g_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,pm20g_nc_ID, work2d, start=(/1,1,ipass/), count=(/nx,ny,1/))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable pm20g')
  end if
  !
  !***  Calculates and writes AOT
  !***  Aerosol optical properties are discretizated in 20 radius intervals
  !***  ranging from 0.01 to 10 mic.
  !
  if( out_aot05 ) then
     work2d = 0.0_rp
     call get_AOT_prop(reff,Qeff)
     !
     do ic = 1,np  ! particles
        if(0.5_rp*diam(ic) >= reff(1 ).and. &
             0.5_rp*diam(ic) <= reff(20)) then                 ! This particle class is an aerosol
           !
           do i = 1,19                                ! 20 effective radius
              if(0.5_rp*diam(ic) >= reff(i  ).and. &
                   0.5_rp*diam(ic) <= reff(i+1)) then
                 ieff = i+1
              end if
           end do
           !
           work3d = c(1:nx,1:ny,1:nz,ic)
           do iz = 1,nz
              work3d(1:nx,1:ny,iz) = work3d(1:nx,1:ny,iz)/Hm1(1:nx,1:ny)  ! concentration
           end do
           !
           do iz = 1,nz-1
              work2d(1:nx,1:ny) = work2d(1:nx,1:ny) + dz(iz)*0.5_rp*  &
                   (work3d(1:nx,1:ny,iz) + work3d(1:nx,1:ny,iz+1))* &   ! Masss in kg/m2
                   0.75_rp*Qeff(ieff)/(rhop(ic)*reff(ieff))
           end do
        end if
     end do
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,aot05_nc_name,aot05_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 )  call runend('printres_nc : Error getting aot05_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,aot05_nc_ID, work2d, start=(/1,1,ipass/), count=(/nx,ny,1/))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable aot05')
  end if
  !
  !***  Calculates and writes concentration at selected flight levels (in gr/m3)
  !***  Concentration refers to particles only
  !
  if( out_fl ) then
     work3d = 0.0_rp
     do ic = 1,np ! particles
        work3d = work3d + c(1:nx,1:ny,1:nz,ic)
     end do
     do iz = 1,nz
        work3d(1:nx,1:ny,iz) = work3d(1:nx,1:ny,iz)/Hm1(1:nx,1:ny)  ! concentration
     end do
     !
     do iflevel = 1,nflevel
        flm =  flevel_value(iflevel)*30.48_rp  ! FL in m
        call setcut_FL(nx,ny,nz,zlayer,work2d,work3d,topg,p,flm)
        work2d(1:nx,1:ny) =  1d3*work2d(1:nx,1:ny)  ! gr/m3
        !
        if( mpime == root ) then
           istat = nf90_inq_varid(ncID,flevel_nc_name(iflevel),flevel_nc_ID(iflevel))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('printres_nc : Error getting iflelvel_nc_ID')
        !
        if( mpime == root ) then
           istat = nf90_put_var(ncID,flevel_nc_ID(iflevel),work2d,start=(/1,1,ipass/),count=(/nx,ny,1/))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable flevel')
     end do
  end if
  !
  !***  Calculates and writes 3D particle concentration (gr/m3)
  !
  if(print_3d_var) then
     work3d = 0.0_rp
     do ic = 1,np
        work3d = work3d + c(1:nx,1:ny,1:nz,ic)
     end do
     do iz = 1,nz
        work3d(1:nx,1:ny,iz) = work3d(1:nx,1:ny,iz)/Hm1(1:nx,1:ny)  ! concentration
     end do
     work3d(1:nx,1:ny,1:nz) = 1d3*work3d(1:nx,1:ny,1:nz)                 ! gr/m3
     if(.not.terrain_following) call ter2asl(work3d,topg,zlayer,nx,ny,nz)   ! a.s.l.
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,conce_nc_name,conce_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : Error getting conce_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,conce_nc_ID, work3d, start=(/1,1,1,ipass/), count=(/nx,ny,nz,1/) )
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable conce')
  end if
  !
  !***  Gravity current
  !
  if(print_3d_var.and.gravity_current) then
     do k = 1,nz
        work3d(1:nx,1:ny,k) = vx(1:nx,1:ny,k)*Hm1(1:nx,1:ny)
     end do
     if(.not.terrain_following) call ter2asl(work3d,topg,zlayer,nx,ny,nz)   ! a.s.l.
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,velu_nc_name,velu_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : Error getting divu_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,velu_nc_ID, work3d, start=(/1,1,1,ipass/), count=(/nx,ny,nz,1/) )
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable velu')
     !
     do k = 1,nz
        work3d(1:nx,1:ny,k) = vy(1:nx,1:ny,k)*Hm1(1:nx,1:ny)
     end do
     if(.not.terrain_following) call ter2asl(work3d,topg,zlayer,nx,ny,nz)   ! a.s.l.
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,velv_nc_name,velv_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : Error getting divu_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,velv_nc_ID, work3d, start=(/1,1,1,ipass/), count=(/nx,ny,nz,1/) )
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable velv')
     !
  end if
  !
  !***  Calculates and writes 3D particle class concentration (gr/m3)
  !
  if(print_3d_var.and.print_class) then
     do ic = 1,nc
        work3d = c(1:nx,1:ny,1:nz,ic)
        do iz = 1,nz
           work3d(1:nx,1:ny,iz) = work3d(1:nx,1:ny,iz)/Hm1(1:nx,1:ny)  ! concentration
        end do
        work3d(1:nx,1:ny,1:nz) = 1d3*work3d(1:nx,1:ny,1:nz)                 ! gr/m3
        if(.not.terrain_following) call ter2asl(work3d,topg,zlayer,nx,ny,nz)   ! a.s.l.
        !
        if( mpime == root ) then
           istat = nf90_inq_varid(ncID,conce_ic_name(ic),conce_ic_ID(ic))
           if(istat /= 0) call wriwar(nf90_strerror(istat))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 )  call runend('printres_nc : Error getting conce_nc_ID')

        if( mpime == root ) then
           istat = nf90_put_var(ncID,conce_ic_ID(ic), work3d, start=(/1,1,1,ipass/), count=(/nx,ny,nz,1/))
           if(istat /= 0) call wriwar(nf90_strerror(istat))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable conce')
     end do
  end if
  !
  !***  Writes 3D mass of aggregates (kg)
  !
  if(aggregation) then
     work3d = Maggr
     if(.not.terrain_following) call ter2asl(work3d,topg,zlayer,nx,ny,nz)   ! a.s.l.
     !
     if( mpime == root ) then
        istat = nf90_inq_varid(ncID,maggr_nc_name,maggr_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : Error getting maggr_nc_ID')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID,maggr_nc_ID, work3d, start=(/1,1,1,ipass/), count=(/nx,ny,nz,1/) )
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('printres_nc : error in nf90_put_var for variable maggr')
  end if
  !
  !*** Close the file
  !
  if( mpime == root ) then
     istat = nf90_close(ncID)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('printres_nc : Error in nf90_close')
  !
  !***  Releases memory
  !
  deallocate(work2d)
  deallocate(work3d)
  !
  return
end subroutine printres_nc
!
!
!
  subroutine print_nc_grid
  !**************************************************************************
  !*
  !*   Writes grid data (including coordinate variables) in netCDF format
  !*
  !**************************************************************************
  use KindType
  use InpOut
  use Master
  use Numeric
  use netcdf
  use Parallel
  implicit none
  !
  integer(ip) :: i,j,istat,iflevel,ic
  real(rp), allocatable :: workx(:),worky(:)
  ! character(len=2) :: ext      ! Valid for number of classes < 100
  !
  !*** Open netCDF file (define mode) and get ncID
  !
  if( mpime == root ) then
!     istat = nf90_create(TRIM(fres_nc),NF90_CLOBBER, ncID)
!
!    MN
     istat = nf90_create(TRIM(fres_nc),cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncID) 
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : Error in nf90_create')
  !
  !**  Define dimensions
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     if( mpime == root ) then
        istat = nf90_def_dim(ncID, lon_nc_name , nx, nx_nc_ID )
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_def_dim for lon')
     !
     if( mpime == root ) then
        istat = nf90_def_dim(ncID, lat_nc_name , ny, ny_nc_ID )
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_def_dim for lat')
     !
  case('UTM')
     if( mpime == root ) then
        istat = nf90_def_dim(ncID, x_nc_name, nx, nx_nc_ID )
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_def_dim for x')
     !
     if( mpime == root ) then
        istat = nf90_def_dim(ncID, y_nc_name, ny, ny_nc_ID )
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_def_dim for y')
     !
  END SELECT
  !
  if( mpime == root ) then
     istat = nf90_def_dim(ncID, alt_nc_name , nz, nz_nc_ID )
        if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_def_dim for alt')
  !
  if( mpime == root ) then
     istat =  nf90_def_dim(ncID, time_nc_name, NF90_UNLIMITED, nt_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_def_dim for time')
  !
  if( mpime == root ) then
     istat = nf90_def_dim(ncID, np_nc_name , np, np_nc_ID )
        if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_def_dim for np')
  !
  !*** Define coordinate variables
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     !
     if( mpime == root ) then
        istat = nf90_def_var(ncID, lon_nc_name ,NF90_FLOAT, (/nx_nc_ID/), lon_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_variable for variable '//TRIM(lon_nc_name))
     !
     if( mpime == root ) then
        istat = nf90_def_var(ncID, lat_nc_name ,NF90_FLOAT, (/ny_nc_ID/), lat_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_variable for variable '//TRIM(lat_nc_name))
     !
  case('UTM')
     !
     if( mpime == root ) then
        istat = nf90_def_var(ncID, x_nc_name ,NF90_FLOAT, (/nx_nc_ID/), x_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_variable for variable '//TRIM(x_nc_name))
     !
     if( mpime == root ) then
        istat =  nf90_def_var(ncID, y_nc_name ,NF90_FLOAT, (/ny_nc_ID/), y_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_variable for variable '//TRIM(y_nc_name))
     !
  END SELECT
  !
  if( mpime == root ) then
     istat = nf90_def_var(ncID, alt_nc_name ,NF90_FLOAT, (/nz_nc_ID/), alt_nc_ID)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) &
       call runend('print_nc_grid : error in nf90_def_variable for variable '//TRIM(alt_nc_name))
  !
  if( mpime == root ) then
     istat = nf90_def_var(ncID, time_nc_name,NF90_FLOAT, (/nt_nc_ID/), time_nc_ID)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) &
       call runend('print_nc_grid : error in nf90_def_variable for variable '//TRIM(time_nc_name))
  !
  !*** Define the rest of variables
  !
  !
  !*** Times
  !
  if( mpime == root ) then
     istat = nf90_def_var(ncID, timeh_nc_name, NF90_CHAR, (/nt_nc_ID/), timeh_nc_ID)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) &
       call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(timeh_nc_name))
  !
  !*** Particle size
  !
  if( mpime == root ) then
     istat = nf90_def_var(ncID, diam_nc_name, NF90_FLOAT, (/np_nc_ID/), diam_nc_ID)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) &
       call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(diam_nc_name))
  !
  !*** Particle density
  !
  if( mpime == root ) then
     istat = nf90_def_var(ncID, rhop_nc_name, NF90_FLOAT, (/np_nc_ID/), rhop_nc_ID)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) &
       call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(rhop_nc_name))
  !
  !*** Particle sphericity
  !
  if( mpime == root ) then
     istat = nf90_def_var(ncID, sphe_nc_name, NF90_FLOAT, (/np_nc_ID/), sphe_nc_ID)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) &
       call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(sphe_nc_name))
  !
  !*** Topo
  !
  if( mpime == root ) then
     istat = nf90_def_var(ncID, topog_nc_name, NF90_FLOAT, (/nx_nc_ID,ny_nc_ID/), topog_nc_ID)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) &
       call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(topog_nc_name))
  !
  !*** load
  !
  if( out_load ) then
     if( mpime == root ) then
        istat = nf90_def_var(ncID, cload_nc_name, NF90_FLOAT, (/nx_nc_ID,ny_nc_ID,nt_nc_ID/), cload_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(cload_nc_name))
  end if
  !
  !*** particle class load
  !
  if( out_load ) then
     if(print_class) then
        do ic = 1,np
           cload_ic_name(ic) = TRIM(cload_ic_name(ic))//TRIM(classname(ic))
           if( mpime == root ) then
              istat = nf90_def_var(ncID, cload_ic_name(ic), NF90_FLOAT, &
                   (/nx_nc_ID,ny_nc_ID,nt_nc_ID/), cload_ic_ID(ic))
              if(istat /= 0) call wriwar(nf90_strerror(istat))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) &
                call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(cload_ic_name(ic)))
        end do
     end if
  end if
  !
  !*** wet depo
  !
  if( out_load .and. wetdeposition ) then
     if( mpime == root ) then
        istat = nf90_def_var(ncID, cwet_nc_name, NF90_FLOAT, (/nx_nc_ID,ny_nc_ID,nt_nc_ID/), cwet_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(cwet_nc_name))
  end if
  !
  !*** particle class wet depo
  !
  if( out_load .and. wetdeposition ) then
     if(print_class) then
        do ic = 1,np
           cwet_ic_name(ic) = TRIM(cwet_ic_name(ic))//TRIM(classname(ic))
           if( mpime == root ) then
              istat = nf90_def_var(ncID, cwet_ic_name(ic), NF90_FLOAT, &
                   (/nx_nc_ID,ny_nc_ID,nt_nc_ID/), cwet_ic_ID(ic))
              if(istat /= 0) call wriwar(nf90_strerror(istat))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) &
                call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(cwet_ic_name(ic)))
        end do
     end if
  end if
  !
  !*** Thickness
  !
  if( out_thick ) then
     if( mpime == root ) then
        istat = nf90_def_var(ncID, thick_nc_name, NF90_FLOAT, (/nx_nc_ID,ny_nc_ID,nt_nc_ID/), thick_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(thick_nc_name))
  end if
  !
  !*** Concentration ground
  !
  if( out_concg ) then
     if( mpime == root ) then
        istat = nf90_def_var(ncID, concg_nc_name, NF90_FLOAT,(/nx_nc_ID,ny_nc_ID,nt_nc_ID/), concg_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(concg_nc_name))
  end if
  !
  !*** Concentration PM05 ground
  !
  if( out_pm05g ) then
     if( mpime == root ) then
        istat = nf90_def_var(ncID, pm05g_nc_name, NF90_FLOAT,(/nx_nc_ID,ny_nc_ID,nt_nc_ID/), pm05g_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(pm05g_nc_name))
  end if
  !
  !*** Concentration PM10 ground
  !
  if( out_pm10g ) then
     if( mpime == root ) then
        istat = nf90_def_var(ncID, pm10g_nc_name, NF90_FLOAT,(/nx_nc_ID,ny_nc_ID,nt_nc_ID/), pm10g_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(pm10g_nc_name))
  end if
  !
  !*** Concentration PM20 ground
  !
  if( out_pm20g ) then
     if( mpime == root ) then
        istat = nf90_def_var(ncID, pm20g_nc_name, NF90_FLOAT,(/nx_nc_ID,ny_nc_ID,nt_nc_ID/), pm20g_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(pm20g_nc_name))
  end if
  !
  !*** Column mass
  !
  if( out_cummc ) then
     if( mpime == root ) then
        istat = nf90_def_var(ncID, cummc_nc_name, NF90_FLOAT,(/nx_nc_ID,ny_nc_ID,nt_nc_ID/), cummc_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(cummc_nc_name))
  end if
  !
  !*** Column mass gas
  !
  if( out_cummc ) then
     if(ng > 0) then
        do ic = 1,ng
           cummc_ic_name(ic) =  TRIM(cummc_ic_name(ic))//TRIM(classname(np+ic))
           if( mpime == root ) then
              istat = nf90_def_var(ncID, cummc_ic_name(ic), NF90_FLOAT, &
                   (/nx_nc_ID,ny_nc_ID,nt_nc_ID/), cummc_ic_ID(ic))
              if(istat /= 0) call wriwar(nf90_strerror(istat))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) &
                call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(cummc_ic_name(ic)))
        end do
     end if
  end if
  !
  !*** Column mass PM05
  !
  if( out_pm05c ) then
     if( mpime == root ) then
        istat = nf90_def_var(ncID, pm05c_nc_name, NF90_FLOAT,(/nx_nc_ID,ny_nc_ID,nt_nc_ID/), pm05c_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(pm05c_nc_name))
  end if
  !
  !*** Column mass PM10
  !
  if( out_pm10c ) then
     if( mpime == root ) then
        istat = nf90_def_var(ncID, pm10c_nc_name, NF90_FLOAT,(/nx_nc_ID,ny_nc_ID,nt_nc_ID/), pm10c_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(pm10c_nc_name))
  end if
  !
  !*** Column mass PM20
  !
  if( out_pm20c ) then
     if( mpime == root ) then
        istat = nf90_def_var(ncID, pm20c_nc_name, NF90_FLOAT,(/nx_nc_ID,ny_nc_ID,nt_nc_ID/), pm20c_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(pm20c_nc_name))
  end if
  !
  !*** Con at FL
  !
  if( out_fl ) then
     do iflevel = 1,nflevel
        if( mpime == root ) then
           istat = nf90_def_var(ncID, flevel_nc_name(iflevel), NF90_FLOAT, &
                (/nx_nc_ID,ny_nc_ID,nt_nc_ID/), flevel_nc_ID(iflevel))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) &
             call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(flevel_nc_name(iflevel)))
     end do
  end if
  !
  !*** AOD
  !
  if( out_aot05 ) then
     if( mpime == root ) then
        istat = nf90_def_var(ncID, aot05_nc_name, NF90_FLOAT,(/nx_nc_ID,ny_nc_ID,nt_nc_ID/), aot05_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(aot05_nc_name))
  end if
  !
  !*** 3D variables
  !
  if(print_3d_var) then
     !
     !***    Concentration
     !
     if( mpime == root ) then
        istat = nf90_def_var(ncID,conce_nc_name,NF90_FLOAT,(/nx_nc_ID,ny_nc_ID,nz_nc_ID,nt_nc_ID/),conce_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
         call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(conce_nc_name))
     !
     !***    Gravity current
     !
     if(gravity_current) then
     if( mpime == root ) then
        istat = nf90_def_var(ncID,velu_nc_name,NF90_FLOAT,(/nx_nc_ID,ny_nc_ID,nz_nc_ID,nt_nc_ID/),velu_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
         call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(velu_nc_name))
     !
     if( mpime == root ) then
        istat = nf90_def_var(ncID,velv_nc_name,NF90_FLOAT,(/nx_nc_ID,ny_nc_ID,nz_nc_ID,nt_nc_ID/),velv_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
         call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(velv_nc_name))
     end if
     !
     !***    Concentration class
     !
     if(print_class) then
        do ic = 1,nc
           conce_ic_name(ic) = TRIM(conce_ic_name(ic))//TRIM(classname(ic))
           if( mpime == root ) then
              istat = nf90_def_var(ncID, conce_ic_name(ic), NF90_FLOAT, &
                   (/nx_nc_ID,ny_nc_ID,nz_nc_ID,nt_nc_ID/), conce_ic_ID(ic))
              if(istat /= 0) call wriwar(nf90_strerror(istat))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) &
                call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(conce_ic_name(ic)))
        end do
     end if
  end if           ! print_3d_var
  !
  if(aggregation) then
     if( mpime == root ) then
        istat = nf90_def_var(ncID,maggr_nc_name,NF90_FLOAT,(/nx_nc_ID,ny_nc_ID,nz_nc_ID,nt_nc_ID/),maggr_nc_ID)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid : error in nf90_def_var for variable '//TRIM(maggr_nc_name))
  end if
  !
  !*** Define attributes for coordinate variables
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     attr_desc  = 'longitude. East positive'
     attr_units = 'degrees_east'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, lon_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, lon_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     attr_desc  = 'latitude. North positive'
     attr_units = 'degrees_north'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, lat_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, lat_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
  case('UTM')
     attr_desc  = 'UTM. West-East distance'
     attr_units = 'm'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, x_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, x_nc_ID, 'description', attr_desc)
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     attr_desc  = 'UTM. South-North distance'
     attr_units = 'm'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, y_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, y_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
  END SELECT
  !
  if(terrain_following) then
     attr_desc  = 'heights (terrain following)'
  else
     attr_desc  = 'heights (a.s.l.)'
  end if
  attr_units = 'm'
  if( mpime == root ) then
     istat = nf90_put_att(ncID, alt_nc_ID, 'units', attr_units)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, alt_nc_ID, 'description', attr_desc)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  attr_desc  = 'time after 0000UTC'
  write(attr_units,  &
       '(''hours since '',i4.4,''-'',i2.2,''-'',i2.2,'' 00:00:0.0'')') &
       ibyr,ibmo,ibdy
  if( mpime == root ) then
     istat = nf90_put_att(ncID, time_nc_ID, 'units', attr_units)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, time_nc_ID, 'description', attr_desc)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  attr_desc  = 'Number of particles'
  attr_units = '-'
  if( mpime == root ) then
     istat = nf90_put_att(ncID, np_nc_ID, 'units', attr_units)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, np_nc_ID, 'description', attr_desc)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  !*** Define attributes for other variables
  !
  attr_desc  = 'Time slices'
  attr_units = 'YYYY:MM:DD:HH:SSSS'
  if( mpime == root ) then
     istat = nf90_put_att(ncID, timeh_nc_ID, 'units', attr_units)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, timeh_nc_ID, 'description', attr_desc)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  attr_desc  = 'Particle diameter'
  attr_units = 'mm'
  if( mpime == root ) then
     istat = nf90_put_att(ncID, diam_nc_ID, 'units', attr_units)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, diam_nc_ID, 'description', attr_desc)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  attr_desc  = 'Particle density'
  attr_units = 'kg/m3'
  if( mpime == root ) then
     istat = nf90_put_att(ncID, rhop_nc_ID, 'units', attr_units)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, rhop_nc_ID, 'description', attr_desc)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  attr_desc  = 'Particle sphericity'
  attr_units = '-'
  if( mpime == root ) then
     istat = nf90_put_att(ncID, sphe_nc_ID, 'units', attr_units)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, sphe_nc_ID, 'description', attr_desc)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  !*** topo
  !
  attr_desc  = 'Topography'
  attr_units = 'm'
  if( mpime == root ) then
     istat = nf90_put_att(ncID, topog_nc_ID, 'units', attr_units)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, topog_nc_ID, 'description', attr_desc)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  !*** load
  !
  if( out_load ) then
     attr_desc  = 'Particle ground deposit load'
     attr_units = 'kg/m2'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, cload_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, cload_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  end if
  !
  !*** class load
  !
  if( out_load ) then
     if(print_class) then
        do ic = 1,np
           attr_desc  = 'Deposit load for particle '//TRIM(classname(ic))
           attr_units = 'kg/m2'
           if( mpime == root ) then
              istat = nf90_put_att(ncID, cload_ic_ID(ic), 'units', attr_units)
              if(istat /= 0) call wriwar(nf90_strerror(istat))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
           !
           if( mpime == root ) then
              istat = nf90_put_att(ncID, cload_ic_ID(ic), 'description', attr_desc)
              if(istat /= 0) call wriwar(nf90_strerror(istat))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
        end do
     end if
  end if
  !
  !*** wet depo
  !
  if( out_load .and. wetdeposition ) then
     attr_desc  = 'Wet deposition particle ground deposit'
     attr_units = 'kg/m2'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, cwet_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, cwet_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  end if
  !
  !*** class load wet depo
  !
  if( out_load .and. wetdeposition ) then
     if(print_class) then
        do ic = 1,np
           attr_desc  = 'Wet deposition deposit load for particle '//TRIM(classname(ic))
           attr_units = 'kg/m2'
           if( mpime == root ) then
              istat = nf90_put_att(ncID, cwet_ic_ID(ic), 'units', attr_units)
              if(istat /= 0) call wriwar(nf90_strerror(istat))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
           !
           if( mpime == root ) then
              istat = nf90_put_att(ncID, cwet_ic_ID(ic), 'description', attr_desc)
              if(istat /= 0) call wriwar(nf90_strerror(istat))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
        end do
     end if
  end if
  !
  !*** thickness
  !
  if( out_thick ) then
     attr_desc  = 'Ground deposit thickness'
     attr_units = 'cm'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, thick_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, thick_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  end if
  !
  !*** Concentration at ground
  !
  if( out_concg ) then
     attr_desc  = 'Particle concentration at ground (1st layer)'
     attr_units = 'gr/m3'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, concg_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, concg_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  end if
  !
  !*** Concentration at ground PM05
  !
  if( out_pm05g ) then
     attr_desc  = 'PM05 concentration at ground (1st layer)'
     attr_units = 'gr/m3'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, pm05g_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, pm05g_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  end if
  !
  !*** Concentration at ground PM10
  !
  if( out_pm10g ) then
     attr_desc  = 'PM10 concentration at ground (1st layer)'
     attr_units = 'gr/m3'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, pm10g_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, pm10g_nc_ID, 'description', attr_desc)
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  end if
  !
  !*** Concentration at ground PM20
  !
  if( out_pm20g ) then
     attr_desc  = 'PM20 concentration at ground (1st layer)'
     attr_units = 'gr/m3'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, pm20g_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, pm20g_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  end if
  !
  !*** Col mass
  !
  if( out_cummc ) then
     attr_desc  = 'Particle column mass'
     attr_units = 'gr/m2'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, cummc_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, cummc_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  end if
  !
  !*** Gas column mass
  !
  if( out_cummc ) then
     do ic = 1,ng
        attr_desc  = 'Column mass for '//TRIM(classname(np+ic))
        attr_units = 'kg/m2'
        if( mpime == root ) then
           istat = nf90_put_att(ncID, cummc_ic_ID(ic), 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
        !
        if( mpime == root ) then
           istat = nf90_put_att(ncID, cummc_ic_ID(ic), 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     end do
  end if
  !
  !*** Column mass PM05
  !
  if( out_pm05c ) then
     attr_desc  = 'PM05 column mass'
     attr_units = 'gr/m2'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, pm05c_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, pm05c_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  end if
  !
  !*** Column mass PM10
  !
  if( out_pm10c ) then
     attr_desc  = 'PM10 column mass'
     attr_units = 'gr/m2'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, pm10c_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, pm10c_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  end if
  !
  !*** Column mass PM20
  !
  if( out_pm20c ) then
     attr_desc  = 'PM20 column mass'
     attr_units = 'gr/m2'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, pm20c_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, pm20c_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  end if
  !
  !*** C at FL
  !
  if( out_fl ) then
     do iflevel = 1,nflevel
        attr_desc  = 'Concentration at flight level '//TRIM(flevel_nc_name(iflevel))
        attr_units = 'gr/m3'
        if( mpime == root ) then
           istat = nf90_put_att(ncID, flevel_nc_ID(iflevel), 'units', attr_units)
           if(istat /= 0) call wriwar(nf90_strerror(istat))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
        !
        if( mpime == root ) then
           istat = nf90_put_att(ncID, flevel_nc_ID(iflevel), 'description', attr_desc)
           if(istat /= 0) call wriwar(nf90_strerror(istat))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     end do
  end if
  !
  !*** AOD
  !
  if(  out_aot05 ) then
     attr_desc  = 'Particle Aerosol Optical Depth at 0.5 micron'
     attr_units = '-'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, aot05_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, aot05_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  end if
  !
  !***  3D variables
  !
  if(print_3d_var) then
     !
     !***   concentration
     !
     attr_desc  = 'Total particle concentration'
     attr_units = 'gr/m3'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, conce_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, conce_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     !***   Gravity current
     !
     if(gravity_current) then
     attr_desc  = 'u-velocity'
     attr_units = 'm/s'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, velu_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, velu_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     attr_desc  = 'v-velocity'
     attr_units = 'm/s'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, velv_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, velv_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     end if ! gravity current
     !
     if(print_class) then
        !
        !***     particles
        !
        do ic = 1,nc
           attr_desc  = 'Concentration for class '//TRIM(classname(ic))
           attr_units = 'gr/m3'
           if( mpime == root ) then
              istat = nf90_put_att(ncID, conce_ic_ID(ic), 'units', attr_units)
              if(istat /= 0) call wriwar(nf90_strerror(istat))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
           !
           if( mpime == root ) then
              istat = nf90_put_att(ncID, conce_ic_ID(ic), 'description', attr_desc)
              if(istat /= 0) call wriwar(nf90_strerror(istat))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
        end do
     end if
  end if  ! print_3d_var
  !
  if(aggregation) then
     attr_desc  = 'Cumulative mass of aggregates'
     attr_units = 'kg'
     if( mpime == root ) then
        istat = nf90_put_att(ncID, maggr_nc_ID, 'units', attr_units)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, maggr_nc_ID, 'description', attr_desc)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  end if
  !
  !*** Put global attributes
  !
  attr_title = 'Fall3d 6.2 results'
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'TITLE', attr_title)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  attr_title = TRIM(coord_sys)
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'COORDINATES', attr_title)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'YEAR', ibyr )
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'MONTH', ibmo )
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'DAY', ibdy )
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'RUN_START', INT(beg_time))
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'RUN_END', INT(end_time))
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'LONMIN', lonmin)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'LATMIN', latmin)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'LONMAX', lonmax)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'LATMAX', latmax)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
  case('UTM')
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'XMIN', xmin)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'YMIN', ymin)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'XMAX', xmax)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'YMAX', ymax)
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
     !
  END SELECT
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'NUMBER_PARTICLES', np)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'NUMBER_GAS_SPECIES', ng)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : error in nf90_put_att')
  !
  !*** Leave the define mode
  !
  if( mpime == root ) then
     istat = nf90_enddef(ncID)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid: error in nf90_enddef')
  !
  !*** Write coordinate variables and time-independent variables
  !
  allocate(workx(nx))
  allocate(worky(ny))
  do i = 1,nx
     workx(i) = xorigr + (i-1)*dx
  end do
  do j = 1,ny
     worky(j) = yorigr + (j-1)*dy
  end do

  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     if( mpime == root ) then
        istat = nf90_put_var(ncID, lon_nc_ID, workx(1:nx))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid: error in nf90_put_var for varialbe lon')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID, lat_nc_ID, worky(1:ny))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid: error in nf90_put_var for varialbe lat')
     !
  case('UTM')
     if( mpime == root ) then
        istat = nf90_put_var(ncID, x_nc_ID, workx(1:nx))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid: error in nf90_put_var for varialbe x')
     !
     if( mpime == root ) then
        istat = nf90_put_var(ncID, y_nc_ID, worky(1:ny))
        if(istat /= 0) call wriwar(nf90_strerror(istat))
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) &
          call runend('print_nc_grid: error in nf90_put_var for varialbe y')
     !
  END SELECT
  !
  if( mpime == root ) then
     istat = nf90_put_var(ncID, alt_nc_ID, zlayer(1:nz))
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) &
       call runend('print_nc_grid: error in nf90_put_var for varialbe alt')
  !
  if( mpime == root ) then
     diam(1:np) = diam(1:np)*1d3  ! in mm
     istat = nf90_put_var(ncID, diam_nc_ID, diam(1:np))
     diam(1:np) = diam(1:np)/1d3  ! back to m
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) &
       call runend('print_nc_grid: error in nf90_put_var for varialbe diam')
  !
  if( mpime == root ) then
     istat = nf90_put_var(ncID, rhop_nc_ID, rhop(1:np))
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) &
       call runend('print_nc_grid: error in nf90_put_var for varialbe rhop')
  !
  if( mpime == root ) then
     istat = nf90_put_var(ncID, sphe_nc_ID, sphe(1:np))
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) &
       call runend('print_nc_grid: error in nf90_put_var for varialbe sphe')
  !
  if( mpime == root ) then
     istat = nf90_put_var(ncID, topog_nc_ID, topg(1:nx,1:ny))
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) &
       call runend('print_nc_grid: error in nf90_put_var for varialbe topg')
  !
  deallocate(workx)
  deallocate(worky)
  !
  !*** Close the file
  !
  if( mpime == root ) then
     istat = nf90_close(ncID)
     if(istat /= 0) call wriwar(nf90_strerror(istat))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('print_nc_grid : Error in nf90_close')
  !
  return
end subroutine print_nc_grid
!
!
!
subroutine get_AOT_prop(reff,Qeff)
  !************************************************************************
  !*
  !*    Sets Aerosol (<10mic) Optical properties
  !*
  !************************************************************************
  use KindType
  implicit none
  !
  real   (rp) :: reff(20),Qeff(20)
  !
  !***  Load optical properties
  !
  reff( 1)= 0.014_rp
  Qeff( 1)=  0.00213_rp
  reff( 2)= 0.020_rp
  Qeff( 2)=  0.00386_rp
  reff( 3)= 0.028_rp
  Qeff( 3)=  0.00864_rp
  reff( 4)= 0.040_rp
  Qeff( 4)=  0.02324_rp
  reff( 5)= 0.056_rp
  Qeff( 5)=  0.06674_rp
  reff( 6)= 0.079_rp
  Qeff( 6)=  0.18285_rp
  reff( 7)= 0.112_rp
  Qeff( 7)=  0.44613_rp
  reff( 8)= 0.158_rp
  Qeff( 8)=  0.93225_rp
  reff( 9)= 0.224_rp
  Qeff( 9)=  1.63000_rp
  reff(10)= 0.316_rp
  Qeff(10)=  2.35760_rp
  reff(11)= 0.447_rp
  Qeff(11)=  2.83130_rp
  reff(12)= 0.631_rp
  Qeff(12)=  2.90880_rp
  reff(13)= 0.891_rp
  Qeff(13)=  2.72340_rp
  reff(14)= 1.259_rp
  Qeff(14)=  2.50530_rp
  reff(15)= 1.778_rp
  Qeff(15)=  2.35970_rp
  reff(16)= 2.512_rp
  Qeff(16)=  2.27100_rp
  reff(17)= 3.548_rp
  Qeff(17)=  2.21020_rp
  reff(18)= 5.012_rp
  Qeff(18)=  2.16410_rp
  reff(19)= 7.079_rp
  Qeff(19)=  2.12830_rp
  reff(20)=10.000_rp
  Qeff(20)=  2.10050_rp
  !
  reff = reff*1e-6_rp    ! From mic to m
  !
  return
end subroutine get_AOT_prop
