   subroutine solvedust
  !********************************************************
  !*
  !*  Computes the source term for the resuspension case
  !*  depending on the emission scheme
  !*
  !********************************************************
  use Master
  use Deposit
  use KindType
  implicit none
  !
  integer(ip)        :: ix,iy,ic,ic_dep,ic2
  !Iteration indexes
  real(rp),parameter :: rho_water = 1d3
  real(rp)           :: rho_part,rho_bulk
  integer(ip)        :: isoil
  !rho_water: Density of water (1000 kg/m3)
  !rho_part: Mean Particle Density (kg/m3)
  !rho_bulk: Soil Bulk Density (kg/m3)
  !isoil: Soil Category
  real(rp)    :: wg,wv,f_w
  !wg: Gravimetric soil moisture
  !wv: Volumetric soil moisture
  !f_w: Factor correction (moisture)
  real(rp)    :: f_l, lambda
  !f_l: Factor correction (vegetation)
  !lambda: Frontal area index
  real(rp)    :: emis_step,sup0
  !emis_step: emitted mass on the current time step (kg m-2)
  !sup0: cell area (m2)
  real(rp)    :: d_cm,rhoa_cm,rhop_cm,Re,g_cm,K
  !Marticorena scheme variables
  real(rp)    :: gama,d_m,alfa,F_H(100),p(100)
  !Shao scheme variables
  !
  !***  Computes the threshold friction velocity of a dry bare surface
  !***  This is done for ALL particle classes IN THE DEPOSIT
  !***  (not only for the resuspended ones).
  !
  SELECT CASE(emission_scheme)
  case('WESTPHAL')
    !
    do iy = 1,ny
    do ix = 1,nx
    do ic = 1,nc_dep
      tust(ix,iy,ic) = tust_cte
    end do
    end do
    end do
    !
  case('MARTICORENA')
    !
    g_cm = 981.0_rp   ! cm/s2
    do iy = 1,ny
    do ix = 1,nx
    do ic = 1,nc_dep
      d_cm    = diam_dep(ic)*1d2   ! in cm
      rhoa_cm = densi(ix,iy)/1d3   ! g/cm3
      rhop_cm = rhop_dep(ic)/1d3
      Re      = 1331_rp*(d_cm**1.56_rp) + 0.38_rp
      K       = sqrt( (rhop_cm*g_cm*d_cm/rhoa_cm)*  &
                (1.0_rp+ (0.006_rp/(rhop_cm*g_cm*(d_cm**2.5_rp)))))
    if(Re.lt.0.03) then
      tust(ix,iy,ic) = (0.129_rp*K)/sqrt( (1.928_rp*(0.03_rp**0.092_rp))-1.0_rp )
    else if( re.le.10.0) then
      tust(ix,iy,ic) = (0.129_rp*K)/sqrt( (1.928_rp*(Re**0.092_rp))-1.0_rp )
    else
      tust(ix,iy,ic) = 0.129_rp*k*(1.0_rp-0.0858_rp*exp(-0.0617_rp*(Re-10.0_rp)))
    end if
      tust(ix,iy,ic) = tust(ix,iy,ic)/1d2  ! in m/s
    end do
    end do
    end do
    !
  case('SHAO')
    !
    !gama = 3d-4  ! kg s-2
    gama = 1.65d-4
    do iy = 1,ny
    do ix = 1,nx
    do ic = 1,nc_dep
      d_m  = diam_dep(ic)  ! in m
      tust(ix,iy,ic) = sqrt( 0.0123_rp*( (rhop_dep(ic)*g*d_m/densi(ix,iy)) + (gama/(densi(ix,iy)*d_m)) ))
    end do
    end do
    end do
    !
  END SELECT
  !
  !***  Applies soil moisture correction
  !
  if(moisture_correction) then
    !
     do iy = 1,ny
     do ix = 1,nx
        !
        !***  Compute mean particle density from deposit data
        !
        if(load(ix,iy).gt.mindep) then
          rho_part = 0.0d0
          do ic = 1,nc_dep
            rho_part = rho_part + rhop_dep(ic)*clas(ix,iy,ic)
          end do
          rho_part=rho_part/load(ix,iy)                          !Mean Particle Density
        else
          rho_part=2650.0_rp                                     !Default particle density
        end if
        !
        !***  Compute gravimetric soil moisture
        !
        isoil    = NINT(soil(ix,iy))                             !Soil Category
        rho_bulk = (1-maxsmc(isoil))*rho_part                    !Soil Bulk Density
        !
        if (isoil .le. 12) then                                  !Exclude water bodies, organic material, etc...
          wv = smoi(ix,iy)                                       !Volumetric soil moisture (First layer) 
          wg = 1d2*wv*rho_water/rho_bulk                         !Gravimetric soil moisture
        else
          wg = 0.0
        end if
        !
        !***  Compute factor correction
        !
        if(wg.le.w_threshold) then                               !Factor correction 
          f_w = 1.0_rp
        else
          f_w = sqrt(1.0_rp + 1.21_rp*((wg-w_threshold)**(0.68_rp)))
        end if
        !
        do ic = 1,nc_dep
          tust(ix,iy,ic) = f_w*tust(ix,iy,ic)
        end do
     end do
     end do
     !
  end if
  !
  !***  Applies vegetation correction
  !
  if(vegetation_correction) then
    !
     do iy = 1,ny
     do ix = 1,nx
       lambda = -0.35*LOG(1-0.01_rp*vegfra(ix,iy))
       if(lambda>1.0) lambda=1.0_rp
       !f_l = SQRT(1-0.5_rp*lambda)*SQRT(1+45.0_rp*lambda)
       f_l = SQRT(1-lambda)*SQRT(1+45.0_rp*lambda)
       !
       do ic = 1,nc_dep
         tust(ix,iy,ic) = f_l*tust(ix,iy,ic)
       end do
     end do
     end do
     !
  end if


  !
  !***  Computes the source term depending on the emission scheme (in kg m-2 s-1 )
  !
  emis(1:nx,1:ny,1:nc) = 0.0d0
  !
  SELECT CASE(emission_scheme)
  case('WESTPHAL')
    !
    do iy = 1,ny
    do ix = 1,nx
      isoil = NINT(soil(ix,iy))
      if (isoil .le. 12) then                 ! exclude water bodies (see, lakes, etc) where resuspension is impossible
      if (load(ix,iy).gt.mindep) then           ! It considers a minimum load mass
        do ic = 1,nc
          ic_dep = nc_dep-nc+ic
          if(clas(ix,iy,ic_dep).gt.emis_mas(ix,iy,ic)) then
          if(ust(ix,iy).ge.tust(ix,iy,ic_dep)) then               ! u* above threshold value
            emis(ix,iy,ic) = 1d-5*(ust(ix,iy)**4.0_rp)            !  kg m-2 s-1
          end if
          end if
        end do
      end if
      end if
    end do
    end do
    !
  case('MARTICORENA')
    !
    K = 5.4d-4  ! m-1
    do iy = 1,ny
    do ix = 1,nx
      isoil = NINT(soil(ix,iy))
      if (isoil .le. 12) then                 ! exclude water bodies (see, lakes, etc) where resuspension is impossible
      if (load(ix,iy).gt.mindep) then           ! It considers a minimum load mass
        do ic = 1,nc
          ic_dep = nc_dep-nc+ic
          if(clas(ix,iy,ic_dep).gt.emis_mas(ix,iy,ic)) then
          if(ust(ix,iy).ge.tust(ix,iy,ic_dep)) then                              ! u* above threshold value
            emis(ix,iy,ic) = K*densi(ix,iy)*ust(ix,iy)* &
              (ust(ix,iy)*ust(ix,iy)-tust(ix,iy,ic_dep)*tust(ix,iy,ic_dep))/g    !  kg m-2 s-1
          end if
          end if
        end do
      end if
      end if
    end do
    end do
    !
  case('SHAO')
    !
    do iy = 1,ny
    do ix = 1,nx
      isoil = NINT(soil(ix,iy))
      if (isoil .le. 12) then                 ! exclude water bodies (see, lakes, etc) where resuspension is impossible
      if (load(ix,iy).gt.mindep) then         ! It considers a minimum load mass
        !
        !*** Horizontal flux of saltating particles (all deposit classes)
        !*** Weight of each saltating class
        !
        F_H = 0.0_rp
        do ic = 1,nc_dep
          if(diam_dep(ic).ge.76d-6 .and. vset_dep(ic)>0.2*ust(ix,iy)) then
          if(ust(ix,iy).ge.tust(ix,iy,ic)) then
            F_H(ic) = densi(ix,iy)*ust(ix,iy)*ust(ix,iy)*ust(ix,iy)/g* &
              (1.0_rp-(tust(ix,iy,ic)*tust(ix,iy,ic)/(ust(ix,iy)*ust(ix,iy))))
          end if
          end if
          p(ic) = clas(ix,iy,ic)/load(ix,iy)
        end do
        !
        !*** Vertical flux (emission)
        !
        do ic = 1,nc
          ic_dep = nc_dep-nc+ic
          if(vset_dep(ic_dep)<0.2*ust(ix,iy)) then
          if(clas(ix,iy,ic_dep).gt.emis_mas(ix,iy,ic)) then
            do ic2 = 1,ic_dep
!              alfa = (0.6_rp*log10(diam_dep(ic2)*1d3)+1.6_rp)* &
!                exp(-140.0_rp*diam(ic)*1d3)                             ! d in mm
              if (diam_dep(ic2).ge.76d-6) then
                alfa = 1D-4 * (0.125_rp*LOG(diam_dep(ic2)*1D3)+0.328_rp) * &
                        EXP(-140.7_rp*diam(ic)*1D3 + 0.37)
!                       EXP(-100.0_rp*diam(ic)*1D3)
                alfa = alfa * 5.0/3.0 * rhop(ic) / densi(ix,iy) * 9.81_rp
                emis(ix,iy,ic) = emis(ix,iy,ic) + &
                  alfa*p(ic2)*p(ic_dep)*F_H(ic2)/ &
                  (tust(ix,iy,ic_dep)*tust(ix,iy,ic_dep))
              end if
            end do
          end if
          end if
        end do
        !
      end if
      end if
    end do
    end do
    !
  END SELECT
  !
  !*** Corrects the emission by a factor
  !
  do iy = 1,ny
  do ix = 1,nx
  do ic = 1,nc
     emis(ix,iy,ic) = emission_factor*emis(ix,iy,ic)  ! kg s-1 m-2
  end do
  end do
  end do
  !
  !*** Total emitted mass
  !
  do iy = 1,ny
  do ix = 1,nx
    do ic = 1,nc
      ic_dep = nc_dep-nc+ic
      emis_step = (ieend1-iebeg1)*emis(ix,iy,ic)                                   ! kg m-2
      if (clas(ix,iy,ic_dep) < emis_mas(ix,iy,ic) + emis_step) then
        emis_step          = clas(ix,iy,ic_dep) - emis_mas(ix,iy,ic)
        emis(ix,iy,ic)     = emis_step / (ieend1-iebeg1)                           ! kg s-1 m-2
        emis_mas(ix,iy,ic) = clas(ix,iy,ic_dep)                                    ! kg m-2
      else
        emis_mas(ix,iy,ic) = emis_mas(ix,iy,ic) + emis_step
      end if
    end do
  end do
  end do
  !
  !*** Scale de source term depending of the cell area (result in kg s-1)
  !
  do iy = 1,ny
  do ix = 1,nx
     sup0 = Hm1(ix,iy)*dX1*dX2
     do ic = 1,nc
        emis(ix,iy,ic) = emis(ix,iy,ic)*sup0     !  kg s-1
     end do
  end do
  end do
  !
  !*** Writes variables for this time step
  !
  call wriresu_res
  !
  return
  end subroutine solvedust
