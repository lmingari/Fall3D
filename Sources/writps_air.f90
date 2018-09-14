subroutine writps_air
  !**********************************************************
  !*
  !*    Write the tracking point files
  !*
  !*********************************************************
  use KindType
  use InpOut
  use Master
  use Numeric
  use Parallel
  use netcdf
  implicit none
  !
  character(len=50) :: header1
  integer(ip), save :: ipass = 0
  integer(ip) :: iyr,imo,idy,ihr,imi,ise,i,iz,ic,istat,ieff,jeff
  real   (rp) :: outtime,s,t,st,val(4),aod_pts,cc
  real   (rp) :: reff(20),Qeff(20)
  real   (rp), allocatable :: cz_pts(:),extz_pts(:),ncoz_pts(:)
  !
  !***  Calculates time header1 in the form DDMMMM-HH:MM
  !
  outtime = beg_time + (ippfile-1)*print_time
  outtime = outtime/3600.0_rp
  !
  call addtime(ibyr,ibmo,ibdy,0,iyr,imo,idy,ihr,imi,ise,time)
  header1 = '00aaa0000'
  write(header1(1 :2 ),'(i2.2)') idy
  header1(3 :5 ) = month(imo)
  write(header1(6 :7 ),'(i2.2)') ihr
  write(header1(8 :9 ),'(i2.2)') imi
  !
  !***  First time creates the netCDF files
  !
  if(ipass == 0) call writps_air0
  ipass = ipass + 1
  !
  !***  Allocates memory
  !
  allocate(cz_pts  (nz))
  allocate(extz_pts(nz))
  allocate(ncoz_pts(nz))
  !
  !***  Calculates aerosol optical properties.
  !***  These discretizated in 20 radius intervals ranging
  !***  from 0.01 to 10 mic (MODTRAN).
  !
  call get_AOT_prop(reff,Qeff)
  !
  !***  Loop over points
  !
  do i = 1,npts
     if(use_pts(i)) then
        !
        !***      Start to write results for the current time step
        !***      Open netCDF file and get ncID
        !
        if( mpime == root ) then
           istat = nf90_open(TRIM(name_file_pts(i))//'.res.nc',NF90_WRITE, ncID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air : Error in nf90_open')
        !
        !***      Writes time variable (current step)
        !
        if( mpime == root ) then
           istat = nf90_inq_varid(ncID,time_nc_name,time_nc_ID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air : Error getting time_nc_ID')
        !
        if( mpime == root ) then
           istat = nf90_put_var(ncID,time_nc_ID,outtime,start=(/ipass/))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air : error in nf90_put_var for variable time')
        !
        if( mpime == root ) then
           istat = nf90_inq_varid(ncID,timeh_nc_name,timeh_nc_ID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air : Error getting timeh_nc_ID')
        !
        if( mpime == root ) then
           istat = nf90_put_var(ncID,timeh_nc_ID, header1, start=(/ipass/), count=(/1/) )
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air : error in nf90_put_var for variable timeh')
        !
        !***      Interpolates results
        !
        s = spts(i)
        t = tpts(i)
        st = s*t
        val(1) = (1.0_rp-t-s+st)*0.25_rp                   !  4      3
        val(2) = (1.0_rp-t+s-st)*0.25_rp                   !
        val(3) = (1.0_rp+t+s+st)*0.25_rp                   !
        val(4) = (1.0_rp+t-s-st)*0.25_rp                   !  1      2
        !
        !***      Computes and writes concentration, extinction and AOD
        !
        cz_pts  (1:nz) = 0.0_rp
        extz_pts(1:nz) = 0.0_rp
        ncoz_pts(1:nz) = 0.0_rp
        aod_pts = 0.0_rp
        !
        do iz = 1,nz
           do ic = 1,np         ! loop over particles
              !               Conce
              cc =  val(1)*(c(ipts(i)  ,jpts(i)  ,iz  ,ic)/Hm1(ipts(i)  ,jpts(i)  )) + &
                   val(2)*(c(ipts(i)+1,jpts(i)  ,iz  ,ic)/Hm1(ipts(i)+1,jpts(i)  )) + &
                   val(3)*(c(ipts(i)+1,jpts(i)+1,iz  ,ic)/Hm1(ipts(i)+1,jpts(i)+1)) + &
                   val(4)*(c(ipts(i)  ,jpts(i)+1,iz  ,ic)/Hm1(ipts(i)  ,jpts(i)+1))
              cz_pts(iz) =  cz_pts(iz) + cc
              !               AOD and extinction
              if(0.5*diam(ic) >= reff(1 ).and. &
                   0.5*diam(ic) <= reff(20)) then               ! This particle class is an aerosol
                 !
                 do jeff = 1,19                            ! Select one of the 20 effective radius
                    if(0.5*diam(ic) >= reff(jeff  ).and. &
                         0.5*diam(ic) <= reff(jeff+1)) then
                       ieff = jeff+1
                    end if
                 end do
                 extz_pts(iz) = extz_pts(iz) + (       cc*0.75*Qeff(ieff))/(rhop(ic)*reff(ieff))
                 aod_pts      = aod_pts      + (dz(iz)*cc*0.75*Qeff(ieff))/(rhop(ic)*reff(ieff))
              end if
              !               number concentration
              if(0.5*diam(ic) <= 32e-6_rp) then    ! r < 32 mic
                 ncoz_pts(iz) = ncoz_pts(iz) + (6.0_rp*cc)/(rhop(ic)*pi*diam(ic)*diam(ic)*diam(ic))
              end if
              !
           end do
        end do
        cz_pts(1:nz)   =   cz_pts(1:nz)*1d3   ! gr/m3
        extz_pts(1:nz) = extz_pts(1:nz)*1d6   ! Mm-1
        ncoz_pts(1:nz) = ncoz_pts(1:nz)*1d-6  ! cm-3
        !
        !***      Writes
        !
        !         Concentration
        if( mpime == root ) then
           istat = nf90_inq_varid(ncID,concz_nc_name,concz_nc_ID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 )  call runend('writps_air : Error getting concg_nc_ID')
        !
        if( mpime == root ) then
           istat = nf90_put_var(ncID,concz_nc_ID, cz_pts, start=(/1,ipass/), count=(/nz,1/))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air : error in nf90_put_var for variable concz')
        !
        !         Extinction coefficient
        if( mpime == root ) then
           istat = nf90_inq_varid(ncID,EXT05_nc_name,EXT05_nc_ID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 )  call runend('writps_air : Error getting EXT05_nc_ID')
        !
        if( mpime == root ) then
           istat = nf90_put_var(ncID,EXT05_nc_ID, extz_pts, start=(/1,ipass/), count=(/nz,1/))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air : error in nf90_put_var for variable EXT05')
        !
        !         Number concentration
        if( mpime == root ) then
           istat = nf90_inq_varid(ncID,nconc_nc_name,nconc_nc_ID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 )  call runend('writps_air : Error getting nconc_nc_ID')
        !
        if( mpime == root ) then
           istat = nf90_put_var(ncID,nconc_nc_ID, ncoz_pts, start=(/1,ipass/), count=(/nz,1/))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air : error in nf90_put_var for variable nconc')
        !
        !         AOD
        if( mpime == root ) then
           istat = nf90_inq_varid(ncID,AOT05_nc_name,aot05_nc_ID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 )  call runend('writps_air : Error getting AOT05_nc_name')
        !
        if( mpime == root ) then
           istat = nf90_put_var(ncID,aot05_nc_ID, aod_pts, start=(/ipass/))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air : error in nf90_put_var for variable aot')
        !
        !***      Close the file
        !
        if( mpime == root ) then
           istat = nf90_close(ncID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air : Error in nf90_close')
        !
     end if
  end do
  !
  !***  Reslase memory
  !
  deallocate(cz_pts)
  deallocate(extz_pts)
  deallocate(ncoz_pts)
  !
  return
end subroutine writps_air
!
!
!
subroutine writps_air0
  !*******************************************************************************
  !*
  !*     Defines the netCDF tracking points file for ALL points
  !*
  !*******************************************************************************
  use KindType
  use InpOut
  use Master
  use Numeric
  use Parallel
  use netcdf
  implicit none
  !
  integer(ip) :: i,istat
  !
  if(restart) return   ! File already exists!
  !
  !***  Loop over points
  !
  do i = 1,npts
     if(use_pts(i)) then
        !
        !***     Open netCDF file (define mode) and get ncID
        !
        if( mpime == root ) then
           istat = nf90_create(TRIM(name_file_pts(i))//'.res.nc',NF90_CLOBBER, ncID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : Error in nf90_create')
        !
        !**      Define dimensions
        !
        if( mpime == root ) then
           istat = nf90_def_dim(ncID, alt_nc_name , nz, nz_nc_ID )
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_def_dim for alt')
        !
        if( mpime == root ) then
           istat =  nf90_def_dim(ncID, time_nc_name, NF90_UNLIMITED, nt_nc_ID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_def_dim for time')
        !
        !***     Define coordinate variables
        !
        if( mpime == root ) then
           istat = nf90_def_var(ncID, alt_nc_name ,NF90_FLOAT, (/nz_nc_ID/), alt_nc_ID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_def_variable for variable '//TRIM(alt_nc_name))
        !
        if( mpime == root ) then
           istat = nf90_def_var(ncID, time_nc_name,NF90_FLOAT, (/nt_nc_ID/), time_nc_ID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_def_variable for variable '//TRIM(time_nc_name))
        !
        !***     Define the rest of variables
        !
        !
        !***     AOD (550 nm)
        !
        if( mpime == root ) then
           istat = nf90_def_var(ncID, AOT05_nc_name, NF90_FLOAT, (/nt_nc_ID/), aot05_nc_ID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_def_var for variable '//TRIM(AOT05_nc_name))
        !
        !***     Extinction coefficient at (550nm)
        !
        if( mpime == root ) then
           istat = nf90_def_var(ncID, EXT05_nc_name, NF90_FLOAT, (/nz_nc_ID,nt_nc_ID/), ext05_nc_ID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_def_var for variable '//TRIM(EXT05_nc_name))
        !
        !***     Number concentration
        !
        if( mpime == root ) then
           istat = nf90_def_var(ncID, nconc_nc_name, NF90_FLOAT, (/nz_nc_ID,nt_nc_ID/), nconc_nc_ID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_def_var for variable '//TRIM(nconc_nc_name))
        !
        !***     Times
        !
        if( mpime == root ) then
           istat = nf90_def_var(ncID, timeh_nc_name, NF90_CHAR, (/nt_nc_ID/), timeh_nc_ID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_def_var for variable '//TRIM(timeh_nc_name))
        !
        !***     Concentration
        !
        if( mpime == root ) then
           istat = nf90_def_var(ncID, concz_nc_name, NF90_FLOAT,(/nz_nc_ID,nt_nc_ID/), concz_nc_ID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_def_var for variable '//TRIM(concz_nc_name))
        !
        !***     Define attributes for coordinate variables
        !
        attr_desc  = 'heights above vent'
        attr_units = 'm'
        if( mpime == root ) then
           istat = nf90_put_att(ncID, alt_nc_ID, 'units', attr_units)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        if( mpime == root ) then
           istat = nf90_put_att(ncID, alt_nc_ID, 'description', attr_desc)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        attr_desc  = 'time after 0000UTC'
        write(attr_units,  &
             '(''hours since '',i4.4,''-'',i2.2,''-'',i2.2,'' 00:00:0.0'')') &
             ibyr,ibmo,ibdy
        if( mpime == root ) then
           istat = nf90_put_att(ncID, time_nc_ID, 'units', attr_units)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        if( mpime == root ) then
           istat = nf90_put_att(ncID, time_nc_ID, 'description', attr_desc)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        !***     Define attributes for other variables
        !
        attr_desc  = 'Time instants'
        attr_units = 'DD MMMM HHMM'
        if( mpime == root ) then
           istat = nf90_put_att(ncID, timeh_nc_ID, 'units', attr_units)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        if( mpime == root ) then
           istat = nf90_put_att(ncID, timeh_nc_ID, 'description', attr_desc)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        !        AOD
        attr_desc  = 'AOT at 0.5 mic'
        attr_units = '(-)'
        if( mpime == root ) then
           istat = nf90_put_att(ncID, aot05_nc_ID, 'units', attr_units)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        if( mpime == root ) then
           istat = nf90_put_att(ncID, aot05_nc_ID, 'description', attr_desc)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        !        Extinction
        attr_desc  = 'Extinction coefficient at 0.5 mic'
        attr_units = 'Mm-1'
        if( mpime == root ) then
           istat = nf90_put_att(ncID, EXT05_nc_ID, 'units', attr_units)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        if( mpime == root ) then
           istat = nf90_put_att(ncID, EXT05_nc_ID, 'description', attr_desc)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        !        number concentration
        attr_desc  = 'Number concentration < 32 mic'
        attr_units = 'cm-3'
        if( mpime == root ) then
           istat = nf90_put_att(ncID, nconc_nc_ID, 'units', attr_units)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        if( mpime == root ) then
           istat = nf90_put_att(ncID, nconc_nc_ID, 'description', attr_desc)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        !        Concentration
        attr_desc  = 'Particle concentration'
        attr_units = 'gr/m3'
        if( mpime == root ) then
           istat = nf90_put_att(ncID, concz_nc_ID, 'units', attr_units)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')

        if( mpime == root ) then
           istat = nf90_put_att(ncID, concz_nc_ID, 'description', attr_desc)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        !***     Put global attributes
        !
        attr_title = 'Fall3d 6.2 results'
        if( mpime == root ) then
           istat = nf90_put_att(ncID, NF90_GLOBAL, 'TITLE', attr_title)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        attr_title = TRIM(name_pts(i))
        if( mpime == root ) then
           istat = nf90_put_att(ncID, NF90_GLOBAL, 'POINT_NAME', attr_title)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        attr_title = TRIM(coord_sys)
        if( mpime == root ) then
           istat = nf90_put_att(ncID, NF90_GLOBAL, 'COORDINATES', attr_title)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        SELECT CASE(coord_sys)
        case('LON-LAT')
           if( mpime == root ) then
              istat = nf90_put_att(ncID, NF90_GLOBAL, 'LONGITUDE', xpts(i))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
           !
           if( mpime == root ) then
              istat = nf90_put_att(ncID, NF90_GLOBAL, 'LATITUDE', ypts(i))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
           !
        case('UTM')
           if( mpime == root ) then
              istat = nf90_put_att(ncID, NF90_GLOBAL, 'X', xpts(i))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
           !
           if( mpime == root ) then
              istat = nf90_put_att(ncID, NF90_GLOBAL, 'Y', ypts(i))
           end if
           call bcast( istat, 1, root )
           if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
           !
        END SELECT
        !
        if( mpime == root ) then
           istat = nf90_put_att(ncID, NF90_GLOBAL, 'YEAR', ibyr )
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        if( mpime == root ) then
           istat = nf90_put_att(ncID, NF90_GLOBAL, 'MONTH', ibmo )
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        if( mpime == root ) then
           istat = nf90_put_att(ncID, NF90_GLOBAL, 'DAY', ibdy )
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        if( mpime == root ) then
           istat = nf90_put_att(ncID, NF90_GLOBAL, 'RUN_START', INT(beg_time))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        if( mpime == root ) then
           istat = nf90_put_att(ncID, NF90_GLOBAL, 'RUN_END', INT(end_time))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : error in nf90_put_att')
        !
        !***     Leave the define mode
        !
        if( mpime == root ) then
           istat = nf90_enddef(ncID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0: error in nf90_enddef')
        !
        !***     Write coordinate variables and time-independent variables
        !
        if( mpime == root ) then
           istat = nf90_put_var(ncID, alt_nc_ID, zlayer(1:nz))
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0: error in nf90_put_var for varialbe alt')
        !
        !***     Close the file
        !
        if( mpime == root ) then
           istat = nf90_close(ncID)
        end if
        call bcast( istat, 1, root )
        if( istat /= 0 ) call runend('writps_air0 : Error in nf90_close')
        !
     end if
  end do
  !
  return
end subroutine writps_air0
