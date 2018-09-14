program Fall3d
  !*****************************************************************************
  !*
  !*  FALL3D: 3-D Advection-Diffusion-Sedimentation model for tephra fallout
  !*          SERIAL AND PARALLEL code versions. The code is compiled in
  !*          parallel if the macro __MPI is defined and in serial otherwise.
  !*
  !*  NOTE: The total number of processors available is divided in groups.
  !*        Each group deals with one or several particle classes (the
  !*        number of groups must be equal or lower than the number of classes).
  !*        In turn, a second parallelization on geometry is performed by dividing
  !*        the domain in z-layers. Different layers are assigned to processors
  !*        belonging to a same group.
  !*
  !*    Authors:         Arnau Folch, Antonio Costa, Giovanni Macedonio
  !*    Current version: 7.0
  !*    Version date:    JUNE-2013
  !*
  !*    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN AND ARE
  !*    FURNISHED "AS IS." THE AUTHORS MAKE NO WARRANTY, EXPRESS OR IMPLIED,
  !*    AS TO THE USEFULNESS OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.
  !*    THEY ASSUME NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
  !*    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS
  !*
  !*****************************************************************************
  use KindType
  use InpOut
  use Master
  use Numeric, only: nz
  use Domain,  only: decomp
  implicit none
  !
  integer(ip) :: iarg
  integer(4)  :: iargc
  logical     :: go_on
  !
  !***  Gets filenames from the call arguments
  !
  iarg = 1                          ! inp file
  call GETARG(iarg,finp)
  iarg = 2                          ! src file
  call GETARG(iarg,fsrc)
  iarg = 3                          ! grn file
  call GETARG(iarg,fgrn)
  iarg = 4                          ! dbs file
  call GETARG(iarg,fdbs)
  iarg = 5                          ! log file
  call GETARG(iarg,flst)
  iarg = 6                          ! res file
  call GETARG(iarg,fres)
  iarg = 7                          ! pts file
  call GETARG(iarg,fpts)
  !
  c_cpu_groups(:) = '1'             ! Default=1
  if(iargc() >= 8) then             ! Optional argument (not needed for serial)
     iarg = 8
     call GETARG(iarg,c_cpu_groups)
  end if
  !
  frst    = TRIM(fres)//'.rst.nc'   ! restart file
  fres_nc = TRIM(fres)//'.nc'       ! res.nc file
  !
  !***  Initializes the run
  !
  call setup
  !
  !***  Define local dimension for parallel execution.
  !
  call decomp( nc, nz )
  !
  !***  Starts time integration
  !
  time        = beg_time
  meteo_time  = beg_time
  source_time = beg_time
  !
  go_on = .true.
  iiter = 0
  do while(go_on)
     !
     !*** If necessary reads source data
     !
     if(source_time <= time) then
        if(sourcetime) then
           call source
           if(aggregation) aggregationtime = .true.
        end if
        if(source_time >= end_time) sourcetime = .false.
     end if
     !
     !***  If necessary reads meteo data
     !
     if(meteo_time <= time) then
        if(meteotime) then
           call meteo
           if(aggregation) aggregationtime = .true.
        end if
        if(meteo_time >= end_time) meteotime = .false.
     end if
     !
     !*** If necessary computes aggregation parameters and corrects setling
     !*** velocity of the aggregates.
     !*** Aggregation parameters are recomputed each time either the 
     !*** meteo (air temperature) or the source (plume temperature and
     !*** plume water flow rate) vary. 
     !
     if(aggregationtime) then
        call agrpar
        aggregationtime = .false.
     end if
     !
     time    = time  + dt
     iiter   = iiter + 1
     erumass = erumass + flowrate*dt        ! Total mass (including gas species)
     !
     call setsrc
     call sizloop
     call endstep
     !
     if(time >= end_time) go_on = .false.
  end do
  !
  !***  Final deposit at specific locations
  !
  if(track_points) call writps_end
  !
  !***  Aggregation summary
  !
  if(aggregation) call wriagr
  !
  !***  Writes the rst file
  !
  call wrirst
  !
  !***  Ends the run
  !
  call runend('OK')
  !
end program Fall3d
