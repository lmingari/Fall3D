subroutine setfrac
  !************************************************************
  !*
  !*    Gets the granulometry
  !*
  !*    outputs: fc(nc),rhop(nc),diam(nc),sphe(nc)
  !*
  !************************************************************
  use KindType
  use Master
  implicit none
  !
  integer(ip) :: ic,ig
  real   (rp) :: deltafi,fi,work
  !
  !***  Compute particle diameter in mm
  !
  deltafi = (fimax-fimin)/(nc-1)
  do ic = 1,nc
     fi = fimin + (ic-1)*deltafi
     diam(ic) = 2**(-fi)
  end do
  !
  !***  Density rhop(ic) and sphericity sphe(nc)
  !
  do ic = 1,nc
     fi = fimin + (ic-1)*deltafi
     rhop(ic) = (rhomax -rhomin )*(fi-fimin)/(fimax-fimin) + rhomin
     sphe(ic) = (sphemax-sphemin)*(fi-fimin)/(fimax-fimin) + sphemin
  end do
  !
  !***  Fraction fc(ic) (gaussian or bigaussian)
  !
  do ig = 1,ng
     do ic = 1,nc
        fi = fimin + (ic-1)*deltafi
        call M0_fi_gaussian(work,fi,fimin,fimax,deltafi,fimean(ig),fidisp(ig))
        if(ig>1) then
           fc(ic) = fc(ic) + 0.02*work
        else
           fc(ic) = fc(ic) + 0.98*work
        endif
     end do
  end do
!  fc = fc/ng
  !
  return
end subroutine setfrac
!
!
!
subroutine M0_fi_gaussian(fc,fi,fimin,fimax,deltafi,fimean,fidisp)
  !************************************************************************
  !*
  !*   Computes the mass for a given value of fi between fimin and fimax
  !*   (mass between fi-deltafi/2 and fi+deltafi/2) for a Gaussian distribution.
  !*   The ultimate values fimin and fimax are extended to +- inifinity
  !*   to ensure mass conservation.
  !*
  !*************************************************************************
  use KindType
  implicit none
  real(rp) :: fc,fimin,fimax,deltafi,fimean,fidisp
  !
  real(rp) :: fi,fi2
  real(rp) :: fer
  !
  !***  lower limit
  !
  if(fi.eq.fimin) then
     fi2 = fimin + 0.5_rp*deltafi
     fi2 = (fi2-fimean)/fidisp
     fc  = fer(fi2)
     return
  end if
  !
  !***  upper limit
  !
  if(fi.eq.fimax) then
     fi2 = fimax - 0.5_rp*deltafi
     fi2 = (fi2-fimean)/fidisp
     fc  = 1.0_rp - fer(fi2)
     return
  end if
  !
  !***  other values
  !
  fi2 = fi + 0.5_rp*deltafi
  fi2 = (fi2-fimean)/fidisp
  fc = fer(fi2)
  !
  fi2 = fi - 0.5_rp*deltafi
  fi2 = (fi2-fimean)/fidisp
  fc = fc - fer(fi2)
  !
  return
end subroutine M0_fi_gaussian
!
!
real(kind=8) function fer(x)
  !************************************************************************
  !*
  !*    Computes the area below the typifiyed normal distribution between
  !*    -infinity and x
  !*
  !************************************************************************
  use KindType
  implicit none
  logical  :: go_on
  real(rp) ::  x
  real(rp) ::  t1,t2,dt
  !
  !***  Computes integral between t=0 and t=abs(x)
  !
  fer = 0.0_rp
  dt  = abs(x)/ 10000.0_rp
  !
  t1 = 0.0_rp
  t2 = dt
  go_on = .true.
  do while(go_on)
     fer = fer + 0.5_rp*(exp(-t1*t1/2.0_rp)+exp(-t2*t2/2.0_rp))*dt
     t1 = t2
     t2 = t2 + dt
     if(t2.ge.abs(x)) go_on = .false.
  end do
  fer = fer/sqrt(8.0_rp*atan(1.0))
  !
  if(x.ge.0.0_rp) then
     fer = 0.5_rp + fer
  else
     fer = 0.5_rp - fer
  end if
  !
  return
end function fer
