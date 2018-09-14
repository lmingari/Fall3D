   subroutine getsrc
  !*****************************************************************
  !*
  !*    Gets the source points from the plume solution stored
  !*    in Mplum. The result is stored in the array src
  !*
  !*****************************************************************
  use KindType
  use Master
  implicit none
  !
  logical     :: go_on,found
  integer(ip) :: is,ic,ix,iy,i,iz(4)
  real(rp)    :: shape(4),zp(4),w(4)
  real(rp)    :: x,y,z,s,t
  !
  !***  Initializes the source array
  !
  src(:,:,:,:) = 0.0
  !
  !***  Loop over all the plume discrete points
  !
  do is = 1,ns
     !
     x = Xplum(is)
     found = .false.
     do i=1,nx
        if(x.ge.(xorigr+(i-1)*dx).and.x.le.(xorigr+(i  )*dx)) then
           found = .true.
           ix = i         ! index ix position (left point)
        end if
     end do
     if(.not.found) call runend('getsrc : source position not found in x')
     s = (x-(xorigr+(ix-1)*dx))/dx            ! parameter s in (0,1)
     s = 2.0*s-1.0                          ! parameter s in (-1,1)
     !
     y = Yplum(is)
     found = .false.
     do i=1,ny
        if(y.ge.(yorigr+(i-1)*dy).and.y.lt.(yorigr+(i  )*dy)) then
           found = .true.
           iy = i           ! index iy position (bottom point)
        end if
     end do
     if(.not.found) call runend('getsrc : source position not found in x')
     t = (y-(yorigr+(iy-1)*dy))/dy            ! parameter t in (0,1)
     t = 2.0*t-1.0                          ! parameter s in (-1,1)
     !
     call getshape(s,t,shape)                 ! Nodal values on a (x,y,z=0) plane
     !
     !***     Interpolation on a (x,y) plane. This is necessary because each of the
     !***     four points (ix,iy),(ix+1,iy),(ix+1,iy+1),(ix,iy+1) may have a different
     !***     vertical interpolation duo to topography (i.e. interpolation on a simple
     !***     8-nodes brick is not correct).
     !
     z = Zplum(is)                ! Absolute position, a.s.l. Vent correction already accounted for
     if(terrain_following) then
        zp(1) = topg(ix  ,iy  )   ! Point (ix  ,iy  ) elevation
        zp(2) = topg(ix+1,iy  )   ! Point (ix+1,iy  ) elevation
        zp(3) = topg(ix+1,iy+1)   ! Point (ix+1,iy+1) elevation
        zp(4) = topg(ix  ,iy+1)   ! Point (ix  ,iy+1) elevation
     else
        zp(:) = 0.0
     end if
     !
     do i = 1,4
        iz(i) = 0
        go_on = .true.
        do while(go_on)
           iz(i) = iz(i) + 1
           if((z-zp(i)-dZ0).lt.zmodel(1)) then
              go_on = .false.
              w(i)  = 0.0                                           ! parameter w in (0,1)
           else if((z-zp(i)-dZ0).ge.zmodel(iz(i)   ).and. &
                (z-zp(i)-dZ0).lt.zmodel(iz(i)+1)) then
              go_on = .false.
              w(i) = (z-zp(i)-dZ0-zmodel(iz(i)))/(zmodel(iz(i)+1)-zmodel(iz(i)))  ! parameter w in (0,1)
           else if(iz(i).eq.nz-1) then
              call runend('Source position not found in z')
           end if
        end do
     end do
     !
     do ic = 1,nc
        src(ic,ix,iy,iz(1)  ) = src(ic,ix,iy,iz(1)  )+ &
             (1.0-w(1))*shape(1)*Mplum(ic,is)
        src(ic,ix,iy,iz(1)+1) = src(ic,ix,iy,iz(1)+1)+ &
             (       w(1))*shape(1)*Mplum(ic,is)
        !
        src(ic,ix+1,iy,iz(2)  ) = src(ic,ix+1,iy,iz(2)  )+ &
             (1.0-w(2))*shape(2)*Mplum(ic,is)
        src(ic,ix+1,iy,iz(2)+1) = src(ic,ix+1,iy,iz(2)+1)+ &
             (     w(2))*shape(2)*Mplum(ic,is)
        !
        src(ic,ix+1,iy+1,iz(3)  ) = src(ic,ix+1,iy+1,iz(3)  )+ &
             (1.0-w(3))*shape(3)*Mplum(ic,is)
        src(ic,ix+1,iy+1,iz(3)+1) = src(ic,ix+1,iy+1,iz(3)+1)+ &
             (       w(3))*shape(3)*Mplum(ic,is)
        !
        src(ic,ix,iy+1,iz(4)  ) = src(ic,ix,iy+1,iz(4)  )+ &
             (1.0-w(4))*shape(4)*Mplum(ic,is)
        src(ic,ix,iy+1,iz(4)+1) = src(ic,ix,iy+1,iz(4)+1)+ &
             (       w(4))*shape(4)*Mplum(ic,is)
     end do
     !
  end do
  !
  return
end subroutine getsrc
!
!
!
subroutine getshape(s,t,shape)
  !***************************************************************
  !*
  !*    Evaluates the shape function for 4-nodes a brick at the
  !*    point (s,t).
  !*
  !***************************************************************
  use KindType
  implicit none
  !
  real(rp) :: s,t,shape(4)
  real(rp) :: sm,tm,sp,tp
  !
  sm = 0.5*(1.0-s)
  tm = 0.5*(1.0-t)
  sp = 0.5*(1.0+s)
  tp = 0.5*(1.0+t)
  !
  shape(1) = sm*tm
  shape(2) = sp*tm
  shape(3) = sp*tp
  shape(4) = sm*tp
  !
  return
end subroutine getshape
