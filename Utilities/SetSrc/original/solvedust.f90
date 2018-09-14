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
  integer(ip) :: ix,iy,ic,ic_dep
  real(rp)    :: sup0,gama,d_m,d_cm,rhoa_cm,rhop_cm,Re,g_cm,K
  real(rp)    :: rho_bulk,wg,f_w,alfa,F_H(100),p(100)
  !
  !***  Computes the threshold friction velocity of a dry bare surface
  !***  This is done for ALL particle classes IN THE DEPOSIT (not
  !***  only for the resuspended ones).
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
          tust(ix,iy,ic) = 0.12_rp*k*(1.0_rp-0.0858_rp*exp(-0.0617_rp*(Re-10.0_rp)))
        end if
        tust(ix,iy,ic) = tust(ix,iy,ic)/1d2  ! in m/s
     end do
     end do
     end do
     !
  case('SHAO')
     !
     gama = 3d-4  ! kg s-2
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
        rho_bulk = 0.0d0   ! bulk density
        do ic = 1,nc_dep
           rho_bulk = rho_bulk + rhop_dep(ic)*clas(ix,iy,ic)
        end do
        if( (smoi(ix,iy).le.0.9)) then                        ! exclude water bodies (see, lakes, etc) where resuspension is impossible
           if(load(ix,iy).ge.emis_mas(ix,iy)) then               ! emited mass can not be larger than the existing mass
               rho_bulk = rho_bulk/load(ix,iy)
           else
               rho_bulk = 1d3
           end if
           wg = 1d2*smoi(ix,iy)*1d3/rho_bulk
        else
           wg = w_threshold
        end if
        !
        if(wg.le.w_threshold) then
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
  !***  Computes the source term depending on the emission scheme (in kg m-2 s-1 )
  !
  emis(1:nx,1:ny,1:nc) = 0.0d0
  !
  SELECT CASE(emission_scheme)
  case('WESTPHAL')
     !
     do iy = 1,ny
     do ix = 1,nx
        if(load(ix,iy).ge.emis_mas(ix,iy)) then           ! emited mass can not be larger than the existing mass
        if(smoi(ix,iy).le.0.9  ) then                     ! exclude water bodies (see, lakes, etc) where resuspension is impossible
           do ic = 1,nc
              ic_dep = nc_dep-nc+ic

              if(ust(ix,iy).ge.tust(ix,iy,ic_dep)) then   ! u* above threshold value
                 !
                 emis(ix,iy,ic) = 1d-5*(ust(ix,iy)**4.0_rp)  !  kg m-2 s-1
                 !
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
        if(load(ix,iy).ge.emis_mas(ix,iy)) then               ! emited mass can not be larger than the existing mass
        if(smoi(ix,iy).le.0.9  ) then                         ! exclude water bodies (see, lakes, etc) where resuspension is impossible
           do ic = 1,nc
              ic_dep = nc_dep-nc+ic
              if(ust(ix,iy).ge.tust(ix,iy,ic_dep)) then       ! u* above threshold value
                 !
                 emis(ix,iy,ic) = K*densi(ix,iy)*ust(ix,iy)* &
                                  (ust(ix,iy)*ust(ix,iy)-tust(ix,iy,ic_dep)*tust(ix,iy,ic_dep))/g   !  kg m-2 s-1
                 !
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
        if(load(ix,iy).ge.emis_mas(ix,iy)) then           ! emited mass can not be larger than the existing mass
        if(smoi(ix,iy).le.0.9  ) then                     ! exclude water bodies (see, lakes, etc) where resuspension is impossible
           !
           !*** Horizontal flux of saltating particles (all deposit classes)
           !*** Weight of each saltating class
           !
           do ic = 1,nc_dep
              if(ust(ix,iy).ge.tust(ix,iy,ic)) then                        ! u* above threshold value
                 F_H(ic) = densi(ix,iy)*ust(ix,iy)*ust(ix,iy)*ust(ix,iy)/g* &
                           (1.0_rp-(tust(ix,iy,ic)*tust(ix,iy,ic)/(ust(ix,iy)*ust(ix,iy))))
              else
                 F_H(ic) = 0.0d0
              end if
              p(ic) = clas(ix,iy,ic)/load(ix,iy)
           end do
           !
           !*** Vertical flux (emission)
           !
           do ic = 1,nc
              do ic_dep = nc_dep-nc+ic,nc_dep
                 alfa = (0.6_rp*log10(diam_dep(ic_dep)*1d3)+1.6_rp)*exp(-140.0_rp*diam(ic)*1d3)  ! d in mm
                 emis(ix,iy,ic) = emis(ix,iy,ic) + alfa*p(ic_dep)*F_H(ic_dep)/(tust(ix,iy,nc_dep-nc+ic)*tust(ix,iy,nc_dep-nc+ic))
              end do
           end do
           !
        end if
        end if
     end do
     end do
     !
  END SELECT
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
  !*** Corrects the emission by a factor
  !
  do iy = 1,ny
  do ix = 1,nx
  do ic = 1,nc
     emis(ix,iy,ic) = emission_factor*emis(ix,iy,ic)  ! kg s-1
  end do
  end do
  end do
  !
  !*** Total emited mass
  !
  do iy = 1,ny
  do ix = 1,nx
     sup0 = Hm1(ix,iy)*dX1*dX2
     do ic = 1,nc
        emis_mas(ix,iy) = emis_mas(ix,iy) + (ieend1-iebeg1)*emis(ix,iy,ic)/sup0  ! kg m2
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
