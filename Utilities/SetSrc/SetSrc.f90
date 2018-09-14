  program SetSrc
  !*****************************************************************************
  !*
  !*    AUTHOR  : A.Folch
  !*    VERSION : FALL3D-7.0
  !*    DATE    : DEC-2012
  !*    PURPOSE : This program generates the source files needed by Fall3d-7.0
  !*
  !*****************************************************************************
  use InpOut
  use Master
  implicit none
  !
  logical :: go_on,lexit
  !
  !***  Gets filenames from the call arguments
  !
  iarg = 1                          ! Model name
  call GETARG(iarg,model_name)
  call upcase(model_name)
  iarg = 2                          ! log   file
  call GETARG(iarg,lulogname)
  iarg = 3                          ! input file
  call GETARG(iarg,luinpname)
  iarg = 4                          ! tgs   file
  call GETARG(iarg,lutgsname)
  iarg = 5                          ! granulometry file
  call GETARG(iarg,lugrnname)
  iarg = 6                          ! source file
  call GETARG(iarg,lusrcname)
  iarg = 7                          ! Meteo database file
  call GETARG(iarg,ludbsname)
  iarg = 8                          ! Deposit results file. Only used for the resuspension case
  call GETARG(iarg,ludepname)
  !
  use_mesh = 'YES'
  !
  !***  Opens the log and source files
  !
  call openinp
  !
  !***  Reads input data depending on the mode (eruption or resuspension)
  !
  call readinp
  !
  !***  Loop over time instants (i.e. eruption phases for the eruption mode or
  !***  meteo data intervals for the resuspension mode).
  !
  do idt = 1,ndt
     !                             Time interval
     iebeg1 = idt_src(idt)
     if(idt.lt.ndt) then
        ieend1 = idt_src(idt+1)
     else
        ieend1 = ieend
     end if
     !
     !*** Source-dependent choice
     !
     SELECT CASE(type_source)
     case('POINT','SUZUKI')
         !
         !***  POINT and SUZUKI cases
         !
         SELECT CASE(MER_vs_h)
         case('NONE','ESTIMATE-MASTIN')
            !
            !*** No plume-wind coupling
            !
            M0     = M0_dt    (idt)
            HPlume = HPlume_dt(idt)
            !
            SELECT CASE(type_source)
            case('POINT')
               call solvepoint
            case('SUZUKI')
               Asuzu  = Asuzu_dt (idt)
               Lsuzu  = Lsuzu_dt (idt)
               call solvesuzuki
            END SELECT
            !
            SELECT CASE(use_mesh)
            case('YES')
               call wrisrc_mesh  (iebeg1,ieend1)
            case default
               call wrisrc_nomesh(iebeg1,ieend1)
            END SELECT
            !
            erumass = erumass + (ieend1-iebeg1)*M0                       ! erupted mass
            write(lulog,1) iebeg1,ieend1,M0,HPlume,Zplum(ns),erumass
            call flush(lulog)
            !
         case('ESTIMATE-DEGRUYTER','ESTIMATE-WOODHOUSE')
            !
            !*** Plume-wind coupling
            !
            go_on = .true.
            time  = iebeg1
            do while(go_on)
               !
               !*** Reads meteo data
               !
               call readdbs
               !
               !*** M0 depending on height and wind profile
               !
               HPlume = HPlume_dt(idt)
               call merwind
               !
               SELECT CASE(type_source)
               case('POINT')
                  call solvepoint
               case('SUZUKI')
                  Asuzu  = Asuzu_dt (idt)
                  Lsuzu  = Lsuzu_dt (idt)
                  call solvesuzuki
               END SELECT
               !
               SELECT CASE(use_mesh)
               case('YES')
                  call wrisrc_mesh  (time,min(time2,ieend1))
               case default
                  call wrisrc_nomesh(time,min(time2,ieend1))
               END SELECT
               !
               erumass = erumass + (min(time2,ieend1)-time)*M0    ! erupted mass
               write(lulog,1) time,min(time2,ieend1),M0,HPlume,Zplum(ns),erumass
               call flush(lulog)
               !
               !*** Time update
               !
               time = time2
               if(time.ge.ieend1) go_on = .false.
               !
            end do
            !
         END SELECT
         !
     case('PLUME')
         !
         !***  PLUME case
         !
         M0     = M0_dt    (idt)
         HPlume = HPlume_dt(idt)
         u0     = u0_dt    (idt)
         T0     = T0_dt    (idt)
         w0     = w0_dt    (idt)
         !
         go_on = .true.
         time  = iebeg1
         do while(go_on)
            !
            !*** Reads meteo data
            !
            call readdbs
            !
            call solveplume(lexit)
            if(.not.lexit) &
               call runend('Error when solving the plume equations') ! Checks output
            !
            SELECT CASE(use_mesh)
            case('YES')
               call wrisrc_mesh  (time,min(time2,ieend1))
            case default
               call wrisrc_nomesh(time,min(time2,ieend1))
            END SELECT
            !
            erumass = erumass + (min(time2,ieend1)-time)*M0    ! erupted mass
            write(lulog,1) time,min(time2,ieend1),M0,Zplum(np),Zplum(ns),erumass
            call flush(lulog)
            !
            !*** Time update
            !
            time = time2
            if(time.ge.ieend1) go_on = .false.
            !
         end do
         !
      case('RESUSPENSION')
         !
         time  = iebeg1
         !
         !*** Reads meteo data
         !
         call readdbs
         !
         call solvedust
         !
         call wrisrc_resu(iebeg1,ieend1)
         !
         erumass = erumass + (ieend1-iebeg1)*SUM(emis)                       ! resuspended mass
         write(lulog,2) iebeg1,ieend1,(ieend1-iebeg1)*SUM(emis),erumass
         call flush(lulog)
         !
    END SELECT
  end do          ! idt = 1,ndt
  !
  !***  Writes the rest of file (from the end of the eruption to the end of run)
  !
  SELECT CASE(type_source)
  case('POINT','SUZUKI')
     !
     SELECT CASE(MER_vs_h)
     case('NONE','ESTIMATE-MASTIN')
        if(irend.gt.ieend1) call wrisrc_void(ieend1,irend)
     case('ESTIMATE-DEGRUYTER','ESTIMATE-WOODHOUSE')
        if(irend.gt.time2)  call wrisrc_void(min(time2,ieend1),irend)
     END SELECT
     !
  case('PLUME')
     !
     if(irend.gt.time2) call wrisrc_void(min(time2,ieend1),irend)
     if(irend.gt.time2) call wriplumetem_void(min(time2,ieend1),irend)
     close(lures)
     close(lumas)
     close(luhei)
     close(lutem)
     !
  case('RESUSPENSION')
     !
     if(irend.gt.ieend1) call wrisrc_void(ieend1,irend)
     !
  END SELECT
  !
  !*** log file formats
  !
1 format(1x,i8,3x,i8,2x,e12.5,2x,f9.1,4x,f9.1,4x,e12.5)
2 format(1x,i8,3x,i8,2x,e12.5,2x,e12.5)
  !
  !***  End of program
  !
  close(lusrc)
  call runend('OK')
  !
  end program SetSrc
