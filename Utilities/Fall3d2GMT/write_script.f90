subroutine write_script
  !************************************************************
  !*
  !*   Writes a script piece for a variable
  !*
  !************************************************************
  use KindType
  use InpOut
  use Master
  implicit none
  !
  character(len=3)      :: step
  character(len=1)      :: comma,bar
  character(len=s_mess) :: domain,ndeg
  integer(ip)           :: ival,it,ntt
  real(rp)              :: ndegx,ndegy
  character(len=s_file) :: ptsfile     ! Points file
  character(len=s_file) :: txtfile     ! Texts file
  !
  !*** Initializations
  !
  comma  ='"'
  bar    ='/'
  domain = ' '
  write(domain,'(f7.2,a,f7.2,a,f7.2,a,f7.2)') lonmin,bar,lonmax,bar, &
       latmin,bar,latmax
  call blank_string(domain)
  !
  ndegx = lonmax-lonmin
  ndegy = latmax-latmin
  if(ndegx.le.2) then
     ndeg = '-B01/'
  else if(ndegx.le.10) then
     ndeg = '-B02/'
  else
     ndeg = '-B05/'
  end if
  if(ndegy.le.2) then
     ndeg = TRIM(ndeg)//'01p:. '
  else if(ndegy.le.10) then
     ndeg = TRIM(ndeg)//'02p:. '
  else
     ndeg = TRIM(ndeg)//'05p:. '
  end if
  !
  !*** Writes the contours-VAR.txt file
  !
  fconto = TRIM(outdir)//'/contours-'//TRIM(varname)//'.txt'
  open(99,file=TRIM(fconto),status='unknown')
  do ival = 1,nval
     write(99,'(f12.6,1x,a)') cval(ival),'A'
  end do
  close(99)
  !
  !*** varname header
  !
  write(luout,1) TRIM(varname)
1 format(&
       '#',/, &
       '# ',a,/, &
       '#')
  !
  !*** Initializations
  !
  write(luout,'(''#'',/,''# Set Title font size'')')
  write(luout,'(''gmtset HEADER_FONT_SIZE 16'')')
  write(luout,'(''#'',/,''# Set axis annotation font size'')')
  write(luout,'(''gmtset ANNOT_FONT_SIZE 14'')')
  !
  !
  write(luout,'(''#'',/,''# Set output resolution for JPG files (dpi)'')')
  write(luout,'(''RESOLUTION=120'')')
  !
  write(luout,'(''#'',/,''# Set color used by psxy (optional points)'')')
  write(luout,'(''POINTS_COLOR=black'')')
  write(luout,'(''#'',/,''# Set points size for psxy (optional points)'')')
  write(luout,'(''POINTS_SIZE=5'')')
  !
  !*** Loop over time intstants (nt images for the animated gif)
  !
  if(TRIM(varname).eq.'TOPOGRAPHY') then
     ntt = 1
  else
     ntt = nt
  end if
  do it = 1,ntt
     !
     !***    Filename
     !
     step = '000'
     if(it.lt.10) then
        write(step(3:3),'(i1)') it
     else if(it.lt.100) then
        write(step(2:3),'(i2)') it
     else if(it.lt.1000) then
        write(step(1:3),'(i3)') it
     end if
     fname = TRIM(outdir)//'/'//TRIM(problemname)//'.'//TRIM(varname)//'.'//TRIM(step)//'.ps'
     !
     !***   Plot basemap
     !
     write(luout,10) TRIM(domain),TRIM(ndeg),comma,dates(it),TRIM(varname),TRIM(unit),comma,TRIM(fname)
10   format('#',/,'psbasemap -K -R',a,' -JX15cd/0cd ',a,a,' ',a, &
          ' ',a,' ',a,a,': -P > ',a)
     !
     !***   Plot coastlines
     !
     write(luout,20) TRIM(fname)
20   format('#',/,'pscoast -K -O -R -Na -Df -G220/220/220 -N1 -S0/192/255 -W -J -P >> ',a)

     !
     ! ***   Draw points (optional)
     !
     ptsfile=TRIM(problemdir)//'/points.txt'
     write(luout,21)
21   format('#'/'# Draw points (optional)')
     write(luout,22) trim(ptsfile)
22   format('if [ -r ',a,' ] ; then')
     write(luout,23) trim(ptsfile),trim(fname)
23   format('  psxy ',a,' -K -O -R -Sc${POINTS_SIZE}p -G$POINTS_COLOR -J -P >> ',a)
     write(luout,24)
24   format('fi')

     !
     ! ***   Write labels (optional)
     !
     txtfile=TRIM(problemdir)//'/texts.txt'
     write(luout,25)
25   format('#'/'# Write texts (optional)')
     write(luout,26) trim(txtfile)
26   format('if [ -r ',a,' ] ; then')
     write(luout,27) trim(txtfile),trim(fname)
27   format('  pstext ',a,' -K -O -Xa5p -R -J -P >> ',a) ! Shift labels by 5pt
     write(luout,28) trim(fname)
28   format('  echo | pstext -K -O -Xa-5p -R -J -P >> ',a) ! RESTORE SHIFT !!!!
     write(luout,29)
29   format('fi')

     !
     !***   Plot contours or shaded
     !
     if((it-1).lt.10) then
        write(luout,30) TRIM(luncname),TRIM(varname),it-1,TRIM(fact),TRIM(fconto),TRIM(fname)
     else
        write(luout,31) TRIM(luncname),TRIM(varname),it-1,TRIM(fact),TRIM(fconto),TRIM(fname)
     end if
30   format('#',/,'grdcontour -O ',a,'\?',a,'\[',i1,'\]',a, &
          ' -A+fHelvetica-Bold -R -C',a,' -J -P >> ',a)
31   format('#',/,'grdcontour -O ',a,'\?',a,'\[',i2,'\]',a, &
          ' -A+fHelvetica-Bold -R -C',a,' -J -P >> ',a)
     !
     !***   Convert to Jpeg (E120=dpi)
     !
     write(luout,90) TRIM(fname)
90   format('#',/,'ps2raster -A -Tj -E$RESOLUTION ',a)
     !
  end do
  !
  return
end subroutine write_script
!
!
!
subroutine blank_string(str1)
  use KindType
  !
  character(len=s_name) :: str1
  integer(ip)           :: i,ls1,ls2
  !
  ls1 = len_trim(str1)
  ls2 = 0
  do i = 1,ls1
     if(str1(i:i).ne.' ') then
        ls2 = ls2 + 1
        str1(ls2:ls2) = str1(i:i)
     endif
  enddo
  do i=ls2+1,ls1
     str1(i:i) =' '
  end do
  return
end subroutine blank_string
