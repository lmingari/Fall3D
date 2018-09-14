!***************************************************************
!*
!*    Module for math operations
!*
!***************************************************************
MODULE MathFun
  use kindType, ONLY :  ip,rp
  use wind_rotation
  implicit none
  !
CONTAINS
  !
  !
  !
  subroutine set_lnods(nx,ny,nelem,lnods)
    !**************************************************************
    !*
    !*   Builds the nodal connectivities
    !*
    !*   4----3
    !*   |    |
    !*   1----2
    !*
    !**************************************************************
    implicit none
    integer(ip) :: nx,ny,nelem
    integer(ip) :: lnods(4,nelem)
    !
    integer(ip) :: ix,iy,ielem
    !
    ielem = 0
    do iy = 1,ny-1
       do ix = 1,nx-1
          ielem = ielem + 1
          lnods(1,ielem) = (iy-1)*nx + ix
          lnods(2,ielem) = (iy-1)*nx + ix + 1
          lnods(3,ielem) = (iy  )*nx + ix + 1
          lnods(4,ielem) = (iy  )*nx + ix
       end do
    end do
    !
    if(ielem /= nelem) call runend('set_lnods: incorrect number of elements')
    return
  end subroutine set_lnods
  !
  !
  !
  subroutine set_lelpo(nelem,npoin,lon,lat,lonp,latp,lnods,lelpo,s_po,t_po)
    !**************************************************************
    !*
    !*   Builds:
    !*     lelpo(ipoin) = ielem to identify where a point lays
    !*     s_po(ipoin)  = s (interpolation factor)
    !*     t_po(ipoin)  = t (interpolation factor)
    !*
    !**************************************************************
    use InpOut
    use wind_rotation
    implicit none
    integer(ip) :: nelem,npoin
    integer(ip) :: lnods(4,nelem),lelpo(npoin)
    real(rp)    :: lat(*),lon(*),lonp(npoin),latp(npoin),s_po(npoin),t_po(npoin)
    !
    logical     :: found
    integer(ip) :: ipoin,ielem,i
    real   (rp) :: elcod(2,4),coglo(2),coloc(2)
    !
    !*** Loop over all points
    !
    do ipoin = 1,npoin
       coglo(1) = lonp(ipoin)
       coglo(2) = latp(ipoin)
       !
       ! Wind rotations (NOTE: here angle must be passed with sign changed)
       if(rotate_wind) then
          call rotate3d(rotation_lonpole,rotation_latpole,-rotation_angle, &
               coglo(1),coglo(2),coglo(1),coglo(2))
       end if
       !
       found = .false.
       ielem = 0
       do while(.not.found)
          !
          ielem = ielem + 1
          if(ielem == nelem+1) then
             write(lulog,10) lon(lnods(1,    1)),lat(lnods(1,    1)), &
                  lon(lnods(3,nelem)),lat(lnods(3,nelem)), &
                  ipoin, &
                  coglo(1),coglo(2)
10           format(' BL    Corner      : ',f9.1,1x,f9.1,/,&
                  ' TP    Corner      : ',f9.1,1x,f9.1,/,&
                  ' Poin  Number      : ',i8,/,&
                  ' Point Coordinates : ',f9.1,1x,f9.1)
             call runend('set_lelpo : point not found in meteo model domain')
          end if
          do i=1,4
             elcod(1,i) = lon(lnods(i,ielem))
             elcod(2,i) = lat(lnods(i,ielem))
          end do
          call elsest_elq1p1(4_ip,0.0_rp,1.0_rp,elcod,coglo,coloc,found)
          !
       end do
       !
       !***    Element found
       !
       lelpo(ipoin) = ielem
       s_po (ipoin) = MAX(-1.0_rp,MIN(1.0_rp,coloc(1)))
       t_po (ipoin) = MAX(-1.0_rp,MIN(1.0_rp,coloc(2)))
       !
    end do
    !
    return
  end subroutine set_lelpo
  !
  !
  !
  subroutine elsest_elq1p1(pnode,lmini,lmaxi,elcod,coglo,coloc,found)
    !**************************************************************
    !*
    !*    Check if point with global coordinates (x,y)=COGLO is inside
    !*    a triangle P1 or a quadrilateral Q1. The Q1 element is
    !*    divided into two P1 elements. Returns the local coordinates
    !*    (s,t)=COLOC
    !*
    !*    For P1 triangles we have:
    !*    x = (1-s-t)*x1 + s*x2 + t*x3
    !*    y = (1-s-t)*y1 + s*y2 + t*y3
    !*
    !*    This linear problem is solved for (s,t):
    !*         (x3-x1)(y-y1) -(y3-y1)(x-x1)
    !*    s =  ----------------------------
    !*         (x3-x1)(y2-y1)-(y3-y1)(x2-x1)
    !*
    !*         (x-x1)(y2-y1) -(y-y1)(x2-x1)
    !*    t =  ----------------------------
    !*         (x3-x1)(y2-y1)-(y3-y1)(x2-x1)
    !*
    !**************************************************************
    implicit none
    logical                  :: found
    integer(ip), intent(in)  :: pnode
    real(rp),    intent(in)  :: lmini,lmaxi,elcod(2,pnode),coglo(2)
    real(rp),    intent(out) :: coloc(2)
    !
    real(rp)                 :: deter,colo3,x2x1,y2y1,x3x1,y3y1,xx1,yy1
    !
    found = .false.
    !
    !*** P1 and Q1: Check if point is in first triangle: nodes 1-2-3
    !
    x2x1     = elcod(1,2)-elcod(1,1)
    y2y1     = elcod(2,2)-elcod(2,1)
    x3x1     = elcod(1,3)-elcod(1,1)
    y3y1     = elcod(2,3)-elcod(2,1)
    xx1      = coglo(1)  -elcod(1,1)
    yy1      = coglo(2)  -elcod(2,1)
    deter    = 1.0_rp/(x3x1*y2y1-y3y1*x2x1)
    coloc(1) = deter*(x3x1*yy1-y3y1*xx1)
    coloc(2) = deter*(y2y1*xx1-x2x1*yy1)
    !
    if(abs(coloc(1)).le.1d-10) coloc(1) = 0.0_rp
    if(abs(coloc(2)).le.1d-10) coloc(2) = 0.0_rp

    if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
       if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
          colo3 = 1.0_rp-coloc(1)-coloc(2)
          if(colo3>=lmini.and.colo3<=lmaxi) found = .true.
       end if
    end if

    if(pnode==4) then
       !
       ! Q1: Check if point is in second triangle: nodes 1-3-4
       !
       if(.not.found) then
          x2x1     = elcod(1,3)-elcod(1,1)
          y2y1     = elcod(2,3)-elcod(2,1)
          x3x1     = elcod(1,4)-elcod(1,1)
          y3y1     = elcod(2,4)-elcod(2,1)
          xx1      = coglo(1)  -elcod(1,1)
          yy1      = coglo(2)  -elcod(2,1)
          deter    = 1.0_rp/(x3x1*y2y1-y3y1*x2x1)
          coloc(1) = deter*(x3x1*yy1-y3y1*xx1)
          coloc(2) = deter*(y2y1*xx1-x2x1*yy1)
          !
          if(abs(coloc(1)).le.1d-10) coloc(1) = 0.0_rp
          if(abs(coloc(2)).le.1d-10) coloc(2) = 0.0_rp

          if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
             if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
                colo3 = 1.0_rp-coloc(1)-coloc(2)
                if(colo3>=lmini.and.colo3<=lmaxi) found = .true.
             end if
          end if
       end if
       if(found) then
          call elsest_newrap(coglo,coloc,2_ip,4_ip,elcod)
       end if
       !
    end if
    !
  end subroutine elsest_elq1p1
  !
  !
  !
  subroutine elsest_newrap(coglo,coloc,ndime,pnode,elcod)
    !**************************************************************
    !*
    !*    Calculate the inverse transformation (x,y,z)-->(s,t,r)
    !*
    !*    Iterate for f(s_i)=x_i: ds = J^{-1}.xpoin - J^{-1}.f(s_i)
    !*                               = J^{-1}dx
    !*                            ds = s_{i+1}-s_i  (deltas)
    !*                            dx = xpoin-f(s_i) (deltax)
    !*    where the s_i's are the coordinates in the local basis and
    !*    xpoin(idime)'s the real ones.
    !*
    !**************************************************************
    implicit none
    integer(ip), intent(in)  :: ndime,pnode
    real(rp),    intent(in)  :: coglo(ndime),elcod(ndime,pnode)
    real(rp),    intent(out) :: coloc(ndime)
    !
    integer(ip)              :: inode,iiter,jdime,maxit,idime,ntens
    real(rp)                 :: shapf(4),deriv(2,4)
    real(rp)                 :: xjacm(2,2),xjaci(2,2)
    real(rp)                 :: deltx(3),delts(3),xnorm,detja
    real(rp)                 :: rnode,diame,coocg(3)
    real(rp)                 :: shacg(4),dercg(2,4)
    !
    ! Initial condition
    !
    coloc(1)=0.0_rp
    coloc(2)=0.0_rp
    coloc(ndime)=0.0_rp
    if(ndime==1) then
       ntens=1
    else if(ndime==2) then
       ntens=3
    else
       ntens=6
    end if
    call shafun_q1(coloc(1),coloc(2),shapf,deriv)
    !
    ! Element diameter
    !
    do idime=1,ndime
       coocg(idime)=0.0_rp
    end do
    do inode=1,pnode
       do idime=1,ndime
          coocg(idime)=coocg(idime)+elcod(idime,inode)
       end do
    end do
    rnode=1.0_rp/real(pnode)
    do idime=1,ndime
       coocg(idime)=rnode*coocg(idime)
    end do
    call shafun_q1(coocg(1),coocg(2),shacg,dercg)
    call mbmabt(xjacm,elcod,dercg,ndime,ndime,pnode)
    call invmtx(xjacm,xjaci,diame,ndime)
    diame=1.0_rp/(diame**(2.0_rp/real(ndime)))
    !
    ! Initialize dx=coglo-f(s_1) with s_1=(0,0,0)
    !
    do idime=1,ndime
       deltx(idime)=0.0_rp
       do inode=1,pnode
          deltx(idime)=deltx(idime)&
               +shapf(inode)*elcod(idime,inode)
       end do
    end do
    xnorm=0.0_rp
    do idime=1,ndime
       deltx(idime) = coglo(idime)-deltx(idime)
       xnorm        = xnorm+deltx(idime)*deltx(idime)
    end do
    xnorm=xnorm*diame
    iiter=0
    maxit=10
    !
    ! Iterate for f(s_i)=x_i
    !
    do while( (xnorm>1e-8).and.(iiter<=maxit) )
       iiter=iiter+1
       !
       ! Compute J
       !
       do jdime=1,ndime
          do idime=1,ndime
             xjacm(idime,jdime)=0.0_rp
             do inode=1,pnode
                xjacm(idime,jdime)=xjacm(idime,jdime)&
                     +deriv(idime,inode)*elcod(jdime,inode)
             end do
          end do
       end do
       !
       ! Compute J^{-1}
       !
       call invmtx(xjacm,xjaci,detja,ndime)

       do idime=1,ndime
          delts(idime)=0.0_rp
          !
          ! Compute J^{-1}.dx
          !
          do jdime=1,ndime
             delts(idime)=delts(idime)+deltx(jdime)*xjaci(jdime,idime)
          end do
       end do
       do idime=1,ndime
          coloc(idime)=coloc(idime)+delts(idime)
       end do
       if ((coloc(1)>1d99).or.(coloc(2)>1d99).or.(coloc(ndime)>1d99)) then
          iiter=maxit+1
       else
          call shafun_q1(coloc(1),coloc(2),shapf,deriv)
       end if
       !
       ! Compute f_i
       !
       do idime=1,ndime
          deltx(idime)=0.0_rp
          do inode=1,pnode
             deltx(idime)=deltx(idime)&
                  +shapf(inode)*elcod(idime,inode)
          end do
       end do
       !
       ! Compute dx=coglo-f
       !         xnorm=sum ds^2
       !
       xnorm=0.0_rp
       do idime=1,ndime
          deltx(idime)=coglo(idime)-deltx(idime)
          xnorm=xnorm+delts(idime)*delts(idime)
       end do
       xnorm=xnorm*diame
    end do
    if(xnorm>1e-8) coloc(1)=2.0_rp

  end subroutine elsest_newrap
  !
  !
  !
  subroutine mbmabt(a,b,c,n1,n2,n3)
    !********************************************************************
    !*
    !*  This routine evaluates the matrix product A = B Ct, where
    !*  A -> Mat(n1,n2), B -> Mat(n1,n3), C -> Mat(n2,n3)
    !*
    !********************************************************************
    implicit none
    integer(ip) :: n1,n2,n3
    integer(ip) :: i,j,k
    real(rp)    :: a(n1,n2), b(n1,n3), c(n2,n3)
    !
    do i=1,n1
       do j=1,n2
          a(i,j)=0.0_rp
          do k=1,n3
             a(i,j)=a(i,j)+b(i,k)*c(j,k)
          end do
       end do
    end do
    !
    return
  end subroutine mbmabt
  !
  !
  !
  subroutine invmtx(a,b,deter,nsize)
    !********************************************************************
    !
    ! This routine inverts a square matrix A -> Mat(nsize,nsize). The
    ! inverse is stored in B. Its determinant is DETER
    !
    !********************************************************************
    implicit none
    integer(ip), intent(in)  :: nsize
    real(rp),    intent(in)  :: a(nsize,*)
    real(rp),    intent(out) :: deter,b(nsize,*)
    ! integer(ip)              :: jsize,isize
    real(rp)                 :: t1,t2,t3,t4,denom

    if(nsize==1) then
       !
       ! Inverse of a 1*1 matrix
       !
       deter=a(1,1)
       if(deter/=0.0_rp) return
       b(1,1) = 1.0_rp/a(1,1)

    else if(nsize==2) then
       !
       ! Inverse of a 2*2 matrix
       !
       deter=a(1,1)*a(2,2)-a(2,1)*a(1,2)
       if(deter==0.0_rp) return
       denom  = 1.0_rp/deter
       b(1,1) = a(2,2)*denom
       b(2,2) = a(1,1)*denom
       b(2,1) =-a(2,1)*denom
       b(1,2) =-a(1,2)*denom

    else if(nsize==3) then
       !
       ! Inverse of a 3*3 matrix
       !
       t1  = a(2,2)*a(3,3) - a(3,2)*a(2,3)
       t2  =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
       t3  = a(2,1)*a(3,2) - a(3,1)*a(2,2)
       deter = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
       if(deter==0.0_rp) return
       denom  = 1.0_rp/deter
       b(1,1) = t1*denom
       b(2,1) = t2*denom
       b(3,1) = t3*denom
       b(2,2) = ( a(1,1)*a(3,3) - a(3,1)*a(1,3))*denom
       b(3,2) = (-a(1,1)*a(3,2) + a(1,2)*a(3,1))*denom
       b(3,3) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2))*denom
       b(1,2) = (-a(1,2)*a(3,3) + a(3,2)*a(1,3))*denom
       b(1,3) = ( a(1,2)*a(2,3) - a(2,2)*a(1,3))*denom
       b(2,3) = (-a(1,1)*a(2,3) + a(2,1)*a(1,3))*denom

    else if(nsize==4) then
       !
       ! Inverse of a 4*4 matrix
       !
       t1=   a(2,2)*a(3,3)*a(4,4) + a(2,3)*a(3,4)*a(4,2)&
            +a(2,4)*a(3,2)*a(4,3) - a(2,3)*a(3,2)*a(4,4)&
            -a(2,2)*a(3,4)*a(4,3) - a(2,4)*a(3,3)*a(4,2)
       t2=  -a(2,1)*a(3,3)*a(4,4) - a(2,3)*a(3,4)*a(4,1)&
            -a(2,4)*a(3,1)*a(4,3) + a(2,4)*a(3,3)*a(4,1)&
            +a(2,3)*a(3,1)*a(4,4) + a(2,1)*a(3,4)*a(4,3)
       t3=   a(2,1)*a(3,2)*a(4,4) + a(2,2)*a(3,4)*a(4,1)&
            +a(2,4)*a(3,1)*a(4,2) - a(2,4)*a(3,2)*a(4,1)&
            -a(2,2)*a(3,1)*a(4,4) - a(2,1)*a(3,4)*a(4,2)
       t4=  -a(2,1)*a(3,2)*a(4,3) - a(2,2)*a(3,3)*a(4,1)&
            -a(2,3)*a(3,1)*a(4,2) + a(2,3)*a(3,2)*a(4,1)&
            +a(2,2)*a(3,1)*a(4,3) + a(2,1)*a(3,3)*a(4,2)
       deter= a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3 + a(1,4)*t4
       if(deter==0.0_rp) return
       denom = 1.0_rp/deter
       b(1,1) = t1*denom
       b(2,1) = t2*denom
       b(3,1) = t3*denom
       b(4,1) = t4*denom
       b(1,2) =(- a(1,2)*a(3,3)*a(4,4) - a(1,3)*a(3,4)*a(4,2)&
            &   - a(1,4)*a(3,2)*a(4,3) + a(1,3)*a(3,2)*a(4,4)&
            &   + a(1,2)*a(3,4)*a(4,3) + a(1,4)*a(3,3)*a(4,2))*denom
       b(2,2) =(  a(1,1)*a(3,3)*a(4,4) + a(1,3)*a(3,4)*a(4,1)&
            &   + a(1,4)*a(3,1)*a(4,3) - a(1,4)*a(3,3)*a(4,1)&
            &   - a(1,3)*a(3,1)*a(4,4) - a(1,1)*a(3,4)*a(4,3))*denom
       b(3,2) =(- a(1,1)*a(3,2)*a(4,4) - a(1,2)*a(3,4)*a(4,1)&
            &   - a(1,4)*a(3,1)*a(4,2) + a(1,4)*a(3,2)*a(4,1)&
            &   + a(1,2)*a(3,1)*a(4,4) + a(1,1)*a(3,4)*a(4,2))*denom
       b(4,2) =(  a(1,1)*a(3,2)*a(4,3) + a(1,2)*a(3,3)*a(4,1)&
            &   + a(1,3)*a(3,1)*a(4,2) - a(1,3)*a(3,2)*a(4,1)&
            &   - a(1,2)*a(3,1)*a(4,3) - a(1,1)*a(3,3)*a(4,2))*denom
       b(1,3) =(  a(1,2)*a(2,3)*a(4,4) + a(1,3)*a(2,4)*a(4,2)&
            &   + a(1,4)*a(2,2)*a(4,3) - a(1,3)*a(2,2)*a(4,4)&
            &   - a(1,2)*a(2,4)*a(4,3) - a(1,4)*a(2,3)*a(4,2))*denom
       b(2,3) =(- a(1,1)*a(2,3)*a(4,4) - a(1,3)*a(2,4)*a(4,1)&
            &   - a(1,4)*a(2,1)*a(4,3) + a(1,4)*a(2,3)*a(4,1)&
            &   + a(1,3)*a(2,1)*a(4,4) + a(1,1)*a(2,4)*a(4,3))*denom
       b(3,3) =(  a(1,1)*a(2,2)*a(4,4) + a(1,2)*a(2,4)*a(4,1)&
            &   + a(1,4)*a(2,1)*a(4,2) - a(1,4)*a(2,2)*a(4,1)&
            &   - a(1,2)*a(2,1)*a(4,4) - a(1,1)*a(2,4)*a(4,2))*denom
       b(4,3) =(- a(1,1)*a(2,2)*a(4,3) - a(1,2)*a(2,3)*a(4,1)&
            &   - a(1,3)*a(2,1)*a(4,2) + a(1,3)*a(2,2)*a(4,1)&
            &   + a(1,2)*a(2,1)*a(4,3) + a(1,1)*a(2,3)*a(4,2))*denom
       b(1,4) =(- a(1,2)*a(2,3)*a(3,4) - a(1,3)*a(2,4)*a(3,2)&
            &   - a(1,4)*a(2,2)*a(3,3) + a(1,4)*a(2,3)*a(3,2)&
            &   + a(1,3)*a(2,2)*a(3,4) + a(1,2)*a(2,4)*a(3,3))*denom
       b(2,4) =(  a(1,1)*a(2,3)*a(3,4) + a(1,3)*a(2,4)*a(3,1)&
            &   + a(1,4)*a(2,1)*a(3,3) - a(1,4)*a(2,3)*a(3,1)&
            &   - a(1,3)*a(2,1)*a(3,4) - a(1,1)*a(2,4)*a(3,3))*denom
       b(3,4) =(- a(1,1)*a(2,2)*a(3,4) - a(1,2)*a(2,4)*a(3,1)&
            &   - a(1,4)*a(2,1)*a(3,2) + a(1,4)*a(2,2)*a(3,1)&
            &   + a(1,2)*a(2,1)*a(3,4) + a(1,1)*a(2,4)*a(3,2))*denom
       b(4,4) =(  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1)&
            &   + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1)&
            &   - a(1,2)*a(2,1)*a(3,3) - a(1,1)*a(2,3)*a(3,2))*denom
    else
       !
       ! Inverse of a nsize*nsize matrix
       !
       !    do isize=1,nsize
       !        do jsize=1,nsize
       !           b(isize,jsize)=a(isize,jsize)
       !        enddo
       !     enddo
       !     call elsest_invert(b,nsize,nsize)
    end if

  end subroutine invmtx
  !
  !
  !
  subroutine shafun_q1(s,t,myshape,deriv)
    !************************************************************
    !*
    !*   Evaluates Q1 shape funcion and derivatives at (s,t)
    !*
    !************************************************************
    implicit none
    real(rp) :: s,t,myshape(4),deriv(2,4)
    real(rp) :: st
    !
    st=s*t
    myshape(1)=(1.-t-s+st)*0.25                           !  4         3
    myshape(2)=(1.-t+s-st)*0.25                           !
    myshape(3)=(1.+t+s+st)*0.25                           !
    myshape(4)=(1.+t-s-st)*0.25                           !
    deriv(1,1)=(-1.+t)*0.25                             !  1         2
    deriv(1,2)=(+1.-t)*0.25
    deriv(1,3)=(+1.+t)*0.25
    deriv(1,4)=(-1.-t)*0.25
    deriv(2,1)=(-1.+s)*0.25
    deriv(2,2)=(-1.-s)*0.25
    deriv(2,3)=(+1.+s)*0.25
    deriv(2,4)=(+1.-s)*0.25
    !
    return
  end subroutine shafun_q1
  !
  !
  !
  integer(ip) function stoi1(string1)
    !**************************************************************
    !*
    !*    Decodes a character*1 string
    !*
    !**************************************************************
    implicit none
    character(len=1) :: string1
    !
    if(string1.eq.'0') then
       stoi1 = 0
    else if(string1.eq.'1') then
       stoi1 = 1
    else if(string1.eq.'2') then
       stoi1 = 2
    else if(string1.eq.'3') then
       stoi1 = 3
    else if(string1.eq.'4') then
       stoi1 = 4
    else if(string1.eq.'5') then
       stoi1 = 5
    else if(string1.eq.'6') then
       stoi1 = 6
    else if(string1.eq.'7') then
       stoi1 = 7
    else if(string1.eq.'8') then
       stoi1 = 8
    else if(string1.eq.'9') then
       stoi1 = 9
    end if
    return
  end function stoi1
  !
  !
  !
  subroutine nc_interpola2d( &
       nx_WRF,ny_WRF,nelem_WRF, &
       var_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,npoin_DAT,var)
    !******************************************************************
    !*
    !*    Interploates
    !*
    !******************************************************************
    implicit none
    !
    integer(ip) :: nx_WRF,ny_WRF,nelem_WRF
    integer(ip) :: nx,ny,npoin_DAT
    integer(ip) :: lelpo_WRF(npoin_DAT),lnods_WRF(4,nelem_WRF)
    real(rp)    :: var_WRF(nx_WRF,ny_WRF)
    real(rp)    :: s_po_WRF(npoin_DAT),t_po_WRF(npoin_DAT)
    real(rp)    :: var(nx,ny)
    !
    integer(ip) :: ix,iy,ipoin_DAT,i
    integer(ip) :: ip_WRF(4),ix_WRF(4),iy_WRF(4)
    real(rp)    :: s,t,st
    real(rp)    :: myshape(4)
    !
    !*** Loop over points to interpolate
    !
    do iy = 1,ny
       do ix = 1,nx
          ipoin_DAT = (iy-1)*nx + ix
          !
          !***       Find the WRF indexes
          !
          do i=1,4
             ip_WRF(i) = lnods_WRF(i,lelpo_WRF(ipoin_DAT))
             !              iy_WRF(i) = ip_WRF(i)/nx_WRF + 1
             iy_WRF(i) = (ip_WRF(i)-1)/nx_WRF + 1           ! Modified
             ix_WRF(i) = ip_WRF(i)-(iy_WRF(i)-1)*nx_WRF
          end do
          !
          !***       shape funcions
          !
          s = s_po_WRF(ipoin_DAT)
          t = t_po_WRF(ipoin_DAT)
          st=s*t
          myshape(1)=(1.-t-s+st)*0.25_rp                           !  4     3
          myshape(2)=(1.-t+s-st)*0.25_rp                           !
          myshape(3)=(1.+t+s+st)*0.25_rp                           !
          myshape(4)=(1.+t-s-st)*0.25_rp                           !  1     2
          !
          !***      Interpolate
          !
          var(ix,iy) =0.0_rp
          do i=1,4
             var(ix,iy) = var(ix,iy) + myshape(i)*var_WRF(ix_WRF(i),iy_WRF(i))
          end do
          !
       end do
    end do
    !
    return
  end subroutine nc_interpola2d
  !
  !
  !
  subroutine nc_interpola3d( &
       nx_WRF,ny_WRF,nz_WRF,nelem_WRF, &
       var_WRF,zlayer_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,nz,npoin_DAT,var,zlayer)
    !******************************************************************
    !*
    !*    Interploates
    !*
    !******************************************************************
    implicit none
    !
    integer(ip) :: nx_WRF,ny_WRF,nz_WRF,nelem_WRF
    integer(ip) :: nx,ny,nz,npoin_DAT
    integer(ip) :: lelpo_WRF(npoin_DAT),lnods_WRF(4,nelem_WRF)
    real(rp)    :: var_WRF(nx_WRF,ny_WRF,nz_WRF)
    real(rp)    :: zlayer_WRF(nx_WRF,ny_WRF,nz_WRF)
    real(rp)    :: s_po_WRF(npoin_DAT),t_po_WRF(npoin_DAT)
    real(rp)    :: var(nx,ny,nz),zlayer(nz)
    !
    logical     :: found
    integer(ip) :: ix,iy,iz,ipoin_DAT,i,iz_WRF1,iz_WRF2
    integer(ip) :: ip_WRF(4),ix_WRF(4),iy_WRF(4)
    real(rp)    :: s,t,st,z,z_WRF1,z_WRF2,sz,var_WRF1,var_WRF2
    real(rp)    :: myshape(4)
    !
    !*** Loop over points to interpolate
    !
    do iy = 1,ny
       do ix = 1,nx
          ipoin_DAT = (iy-1)*nx + ix
          !
          !***       Find the WRF indexes
          !
          do i=1,4
             ip_WRF(i) = lnods_WRF(i,lelpo_WRF(ipoin_DAT))
             !              iy_WRF(i) = ip_WRF(i)/nx_WRF + 1
             iy_WRF(i) = (ip_WRF(i)-1)/nx_WRF + 1           ! Modified
             ix_WRF(i) = ip_WRF(i)-(iy_WRF(i)-1)*nx_WRF
          end do
          !
          !***       myshape funcions
          !
          s = s_po_WRF(ipoin_DAT)
          t = t_po_WRF(ipoin_DAT)
          st=s*t
          myshape(1)=(1.-t-s+st)*0.25_rp                           !  4     3
          myshape(2)=(1.-t+s-st)*0.25_rp                           !
          myshape(3)=(1.+t+s+st)*0.25_rp                           !
          myshape(4)=(1.+t-s-st)*0.25_rp                           !  1     2
          !
          do iz = 1,nz
             z = zlayer(iz)   ! Height to find
             !
             !***          First value
             !
             z_WRF1 =0.0_rp
             do i=1,4
                z_WRF1 = z_WRF1 + myshape(i)*zlayer_WRF(ix_WRF(i),iy_WRF(i),1)
             end do
             !
             if(z.le.z_WRF1) then    ! do not search
                iz_WRF1 = 1
                iz_WRF2 = 2
                sz =1.0_rp
             else                    ! search in z_WRF
                found = .false.
                iz_WRF1 = 0
                iz_WRF2 = 1
                do while(.not.found)
                   iz_WRF1 = iz_WRF1 + 1
                   iz_WRF2 = iz_WRF2 + 1
                   z_WRF2 =0.0_rp
                   do i=1,4
                      z_WRF2 = z_WRF2 + myshape(i)*zlayer_WRF(ix_WRF(i),iy_WRF(i),iz_WRF2)
                   end do
                   !
                   if( (z.ge.z_WRF1).and.(z.le.z_WRF2) ) then
                      found = .true.
                      sz    = (z_WRF2-z)/(z_WRF2-z_WRF1)
                   else
                      z_WRF1 = z_WRF2
                   end if
                   !                                                  No found (interpolate)
                   if( (iz_WRF2.eq.nz_WRF).and.(.not.found) ) then
                      found   = .true.
                      iz_WRF1 = nz_WRF
                      iz_WRF2 = nz_WRF
                      sz      =1.0_rp
                   end if
                   !
                end do
             end if
             !
             !***      Interpolate
             !
             var_WRF1 =0.0_rp   ! Value at z_WRF1
             var_WRF2 =0.0_rp   ! Value at z_WRF2
             do i=1,4
                var_WRF1 = var_WRF1 + myshape(i)*var_WRF(ix_WRF(i),iy_WRF(i),iz_WRF1)
                var_WRF2 = var_WRF2 + myshape(i)*var_WRF(ix_WRF(i),iy_WRF(i),iz_WRF2)
             end do
             !
             var(ix,iy,iz) = sz*var_WRF1 + (1.0_rp-sz)*var_WRF2
          end do
       end do
    end do
    !
    return
  end subroutine nc_interpola3d
  !
  !
  !
  subroutine nc_assign2d( &
       nx_WRF,ny_WRF,nelem_WRF, &
       var_WRF,lelpo_WRF,lnods_WRF,s_po_WRF,t_po_WRF, &
       nx,ny,npoin_DAT,var)
    !******************************************************************
    !*
    !*    Assigns for variables that cannot be interpolated
    !*    (e.g. for land use)
    !*
    !******************************************************************
    implicit none
    !
    integer(ip) :: nx_WRF,ny_WRF,nelem_WRF
    integer(ip) :: nx,ny,npoin_DAT
    integer(ip) :: lelpo_WRF(npoin_DAT),lnods_WRF(4,nelem_WRF)
    real(rp)    :: var_WRF(nx_WRF,ny_WRF)
    real(rp)    :: s_po_WRF(npoin_DAT),t_po_WRF(npoin_DAT)
    real(rp)    :: var(nx,ny)
    !
    integer(ip) :: ix,iy,ipoin_DAT,i,j
    integer(ip) :: ip_WRF(4),ix_WRF(4),iy_WRF(4)
    real(rp)    :: s,t,st,smax
    real(rp)    :: myshape(4)
    !
    !*** Loop over points to "interpolate"
    !
    do iy = 1,ny
       do ix = 1,nx
          ipoin_DAT = (iy-1)*nx + ix
          !
          !***       Find the WRF indexes
          !
          do i=1,4
             ip_WRF(i) = lnods_WRF(i,lelpo_WRF(ipoin_DAT))
             iy_WRF(i) = ip_WRF(i)/nx_WRF + 1
             ix_WRF(i) = ip_WRF(i)-(iy_WRF(i)-1)*nx_WRF
          end do
          !
          !***       myshape funcions
          !
          s = s_po_WRF(ipoin_DAT)
          t = t_po_WRF(ipoin_DAT)
          st=s*t
          myshape(1)=(1.-t-s+st)*0.25_rp                           !  4     3
          myshape(2)=(1.-t+s-st)*0.25_rp                           !
          myshape(3)=(1.+t+s+st)*0.25_rp                           !
          myshape(4)=(1.+t-s-st)*0.25_rp                           !  1     2
          !
          !***       Closest point (max value of myshape)
          !
          j = 1
          smax = 0.0
          do i=1,4
             if(myshape(i).gt.smax) then
                smax = myshape(i)
                j = i
             end if
          end do
          !
          !***      Assigns
          !
          var(ix,iy) = var_WRF(ix_WRF(j),iy_WRF(j))
          !
       end do
    end do
    !
    return
  end subroutine nc_assign2d
  !
  !
  !
END MODULE MathFun
