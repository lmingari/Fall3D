program Fall3d2GMT
  !*****************************************************************************
  !*
  !*    AUTHOR : A.Folch
  !*    VERSION: FALL3D-7.0
  !*    DATE   : MAY-2013
  !*    PURPOSE: This program generates the GMT Scripts to postprocess Fall3d
  !*             results
  !*
  !*****************************************************************************
  use KindType
  use InpOut
  use Master
  implicit none
  !
  integer(ip) :: ic
  !
  !***  Gets filenames from the call arguments
  !
  iarg = 1                          ! log   file
  call GETARG(iarg,lulogname)
  iarg = 2                          ! input file
  call GETARG(iarg,luinpname)
  iarg = 3                          ! nc    file
  call GETARG(iarg,luncname )
  iarg = 4                          ! problemname
  call GETARG(iarg,problemname)
  iarg = 5                          ! problem directory
  call GETARG(iarg,problemdir)
  iarg = 6                          ! output directory
  call GETARG(iarg,outdir)
  !
  !***  Opens the log file
  !
  call openinp
  !
  !***  Reads input data
  !
  call readinp
  !
  !*** Writes the script header
  !
  fname = TRIM(outdir)//'/Fall3d2GMT'
  open(luout,FILE=TRIM(fname),status='unknown')
  call system('chmod 750 '//TRIM(fname))
  write(luout,1)
1 format('#!/bin/sh                                           ',/, &
       '#                                                    ',/, &
       '# Script generated by the Fall3d2GMT utility program ',/, &
       '# Based on the GMT package                           ',/, &
       '#                                                    ')
  !
  !*** Time independent variables
  !
  !
  !    TOPOGRAPHY
  !
  if(pp_topog) then
     varname = 'TOPOGRAPHY'
     nval = nval_topog
     cval = cval_topog
     fact = fact_topog
     unit = unit_topog
     call write_script
  end if
  !
  !*** Deposit variables
  !
  !
  !    LOAD
  !
  if(pp_load0) then
     varname = 'LOAD'
     nval = nval_load0
     cval = cval_load0
     fact = fact_load0
     unit = unit_load0
     call write_script
  end if
  !
  !    WET
  !
  if(pp_loadw) then
     varname = 'WET'
     nval = nval_loadw
     cval = cval_loadw
     fact = fact_loadw
     unit = unit_loadw
     call write_script
  end if
  !
  !    CLASS-LOAD
  !
  if(pp_loadc) then
     nval = nval_loadc
     cval = cval_loadc
     fact = fact_loadc
     unit = unit_loadc
     do ic = 1,nclass
        varname = TRIM(classname(ic))
        call write_script
     end do
  end if
  !
  !    CLASS-WET
  !
  if(pp_wetc) then
     nval = nval_wetc
     cval = cval_wetc
     fact = fact_wetc
     unit = unit_wetc
     do ic = 1,nclass
        varname = TRIM(classnamew(ic))
        call write_script
     end do
  end if
  !
  !    THICKNESS
  !
  if(pp_thick) then
     varname = 'THICKNESS'
     nval = nval_thick
     cval = cval_thick
     fact = fact_thick
     unit = unit_thick
     call write_script
  end if
  !
  !*** Ground level variables
  !
  !
  !    C_GRND
  !
  if(pp_concg) then
     varname = 'C_GRND'
     nval = nval_concg
     cval = cval_concg
     fact = fact_concg
     unit = unit_concg
     call write_script
  end if
  !
  !    PMxx_GRND
  !
  if(pp_PMxxg) then
     nval = nval_PMxxg
     cval = cval_PMxxg
     fact = fact_PMxxg
     unit = unit_PMxxg
     varname = 'C_PM05_GRND'
     call write_script
     varname = 'C_PM10_GRND'
     call write_script
     varname = 'C_PM20_GRND'
     call write_script
  end if
  !
  !*** Airborne variables
  !
  !
  !    COL_MASS
  !
  if(pp_cumul) then
     varname = 'COL_MASS'
     nval = nval_cumul
     cval = cval_cumul
     fact = fact_cumul
     unit = unit_cumul
     call write_script
  end if
  !
  !    PMxx_COLUMN
  !
  if(pp_PMxxc) then
     nval = nval_PMxxc
     cval = cval_PMxxc
     fact = fact_PMxxc
     unit = unit_PMxxc
     varname = 'COL_MASSPM05'
     call write_script
     varname = 'COL_MASSPM10'
     call write_script
     varname = 'COL_MASSPM20'
     call write_script
  end if
  !
  !    FL
  !
  if(pp_fl) then
     nval = nval_fl
     cval = cval_fl
     fact = fact_fl
     unit = unit_fl
     varname = 'C_FL050'
     call write_script
     varname = 'C_FL100'
     call write_script
     varname = 'C_FL150'
     call write_script
     varname = 'C_FL200'
     call write_script
     varname = 'C_FL250'
     call write_script
     varname = 'C_FL300'
     call write_script
     varname = 'C_FL350'
     call write_script
     varname = 'C_FL400'
     call write_script
  end if
  !
  !    AOD
  !
  if(pp_aod05) then
     varname = 'AOD'
     nval = nval_aod05
     cval = cval_aod05
     fact = fact_aod05
     unit = unit_aod05
     call write_script
  end if
  !
  close(luout)
  !
  !*** Animation file (based on convert utility from ImageMagick)
  !
  fname = TRIM(outdir)//'/Fall3d2animation'
  open(luout,FILE=TRIM(fname),status='unknown')
  call system('chmod 750 '//TRIM(fname))
  write(luout,100) TRIM(outdir)
100 format('#!/bin/sh                                           ',/, &
       '#                                                    ',/, &
       '# Script generated by the Fall3d2GMT utility program ',/, &
       '# Based on the convert utility (ImageMagick)         ',/, &
       '#                                                    ',/, &
       'cd ',a)
  !
  if(pp_load0) write(luout,'(3a)') 'convert -delay 20 -loop 0 *.LOAD.*.jpg ',TRIM(problemname),'.LOAD.animated.gif'
  if(pp_loadw) write(luout,'(3a)') 'convert -delay 20 -loop 0 *.WET.*.jpg ',TRIM(problemname),'.WET.animated.gif'
  if(pp_loadc) then
     do ic = 1,nclass
        write(luout,'(7a)') 'convert -delay 20 -loop 0 *.',TRIM(classname(ic)),'.*.jpg ', &
             TRIM(problemname),'.',TRIM(classname(ic)),'.animated.gif'
     end do
  end if
  if(pp_wetc) then
     do ic = 1,nclass
        write(luout,'(7a)') 'convert -delay 20 -loop 0 *.',TRIM(classnamew(ic)),'.*.jpg ', &
             TRIM(problemname),'.',TRIM(classnamew(ic)),'.animated.gif'
     end do
  end if
  if(pp_thick) write(luout,'(3a)') 'convert -delay 20 -loop 0 *.THICKNESS.*.jpg ',TRIM(problemname),'.THICKNESS.animated.gif'
  !
  if(pp_concg) write(luout,'(3a)') 'convert -delay 20 -loop 0 *.C_GRND.*.jpg '   ,TRIM(problemname),'.C_GRND.animated.gif'
  if(pp_PMxxg) then
     write(luout,'(3a)') 'convert -delay 20 -loop 0 *.C_PM05_GRND.*.jpg '   ,TRIM(problemname),'.C_PM05_GRND.animated.gif'
     write(luout,'(3a)') 'convert -delay 20 -loop 0 *.C_PM10_GRND.*.jpg '   ,TRIM(problemname),'.C_PM10_GRND.animated.gif'
     write(luout,'(3a)') 'convert -delay 20 -loop 0 *.C_PM20_GRND.*.jpg '   ,TRIM(problemname),'.C_PM20_GRND.animated.gif'
  end if
  !
  if(pp_cumul) write(luout,'(3a)') 'convert -delay 20 -loop 0 *.COL_MASS.*.jpg ' ,TRIM(problemname),'.COL_MASS.animated.gif'
  if(pp_PMxxc) then
     write(luout,'(3a)') 'convert -delay 20 -loop 0 *.COL_MASSPM05.*.jpg '   ,TRIM(problemname),'.COL_MASSPM05.animated.gif'
     write(luout,'(3a)') 'convert -delay 20 -loop 0 *.COL_MASSPM10.*.jpg '   ,TRIM(problemname),'.COL_MASSPM10.animated.gif'
     write(luout,'(3a)') 'convert -delay 20 -loop 0 *.COL_MASSPM20.*.jpg '   ,TRIM(problemname),'.COL_MASSPM20.animated.gif'
  end if
  if(pp_fl) then
     write(luout,'(3a)') 'convert -delay 20 -loop 0 *.C_FL050.*.jpg '   ,TRIM(problemname),'.C_FL050.animated.gif'
     write(luout,'(3a)') 'convert -delay 20 -loop 0 *.C_FL100.*.jpg '   ,TRIM(problemname),'.C_FL100.animated.gif'
     write(luout,'(3a)') 'convert -delay 20 -loop 0 *.C_FL150.*.jpg '   ,TRIM(problemname),'.C_FL150.animated.gif'
     write(luout,'(3a)') 'convert -delay 20 -loop 0 *.C_FL200.*.jpg '   ,TRIM(problemname),'.C_FL200.animated.gif'
     write(luout,'(3a)') 'convert -delay 20 -loop 0 *.C_FL250.*.jpg '   ,TRIM(problemname),'.C_FL250.animated.gif'
     write(luout,'(3a)') 'convert -delay 20 -loop 0 *.C_FL300.*.jpg '   ,TRIM(problemname),'.C_FL300.animated.gif'
     write(luout,'(3a)') 'convert -delay 20 -loop 0 *.C_FL350.*.jpg '   ,TRIM(problemname),'.C_FL350.animated.gif'
     write(luout,'(3a)') 'convert -delay 20 -loop 0 *.C_FL400.*.jpg '   ,TRIM(problemname),'.C_FL400.animated.gif'
  end if
  if(pp_aod05) write(luout,'(3a)') 'convert -delay 20 -loop 0 *.AOD.*.jpg ' ,TRIM(problemname),'.AOD.animated.gif'
  !
  close(luout)
  !
  !*** Ends
  !
  call runend('OK')
  !
end program Fall3d2GMT