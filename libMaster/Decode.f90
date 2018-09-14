!***************************************************************
!*
!*    Module for decode operations
!*
!***************************************************************
MODULE decode
  use kindtype
  implicit none
  !
CONTAINS
  !
  subroutine sdecode(card,words,param,nword,npar)
    !********************************************************************
    !*
    !*    This routine decodes a string card(s_long) into words and parameters
    !*
    !********************************************************************
    use kindtype
    implicit none
    integer(ip)  ::    nword,npar
    character(len=s_long) ::  card
    character(len=s_long) ::  words(nwormax)
    character(len=1)      ::  sstring(s_long)
    real(rp)              ::  param(nparmax)
    !
    integer(ip) ::  ipos,first,last,nstr,lflag,i,ii
    real(rp)    ::  digit
    !
    !***  Initializations
    !
    nword = 0
    npar  = 0
    ipos  = 0
    do while(1.ne.0)
       !                                    ! First position
       ipos = ipos + 1
       if(ipos.gt.s_long) return
10     if(card(ipos:ipos).eq.' '.or.card(ipos:ipos).eq.'=') then
          ipos = ipos + 1
          if(ipos.gt.s_long) return
          go to 10
       end if
       first = ipos
       !
       ipos = ipos + 1
       if(ipos.gt.s_long) return
20     if(card(ipos:ipos).ne.' '.and.card(ipos:ipos).ne.'=') then
          ipos = ipos + 1
          if(ipos.gt.s_long) return
          go to 20
       end if
       last = ipos-1
       !
       nstr = last-first+1

       ii = 0
       do i=first,last
          ii = ii + 1
          sstring(ii) = card(i:i)
       end do
       call decod1(sstring,nstr,lflag,digit)
       if(lflag.eq.0) then
          npar = npar + 1
          param(npar)= digit
       else if(lflag.eq.1) then
          nword = nword + 1
          words(nword)(:) = ' '
          words(nword)(1:nstr) = card(first:last)
       end if
       !
    end do
    return
  end subroutine sdecode
  !
  !
  !
  subroutine upcase(word)
    !***********************************************************************
    !*
    !*    This routine converts word to upper case
    !*
    !***********************************************************************
    implicit none
    character(len=*) :: word
    integer(ip)      :: iposi,ioctv,item1,item2,item3
    integer(ip)      :: ilen
    !
    item1 = int(o'141')
    item2 = int(o'172')
    item3 = int(o'40')
    !
    ilen=LEN_TRIM(word)
    !
    do iposi=1,ilen                                   ! process all positions
       ioctv=ichar(word(iposi:iposi))                 ! octal value
       if(item1.le.ioctv.and.item2.ge.ioctv) then     ! it is a lower case
          ioctv=ioctv-item3                           ! equivalent upper case
          word(iposi:iposi)=char(ioctv)               ! convert it to upcase
       end if
    end do ! iposi=1,ilen
    !
  end subroutine upcase
  !
  !
  !
  subroutine decod1(string,nstr,lflag,digit)
    !*******************************************************************
    !*
    !*    This subroutine decodes a single string(1:nstr)
    !*
    !*    If string(1:nstr) is a string returns lflag = 1
    !*    If string(1:nstr) is a number returns lflag = 0 and the digit
    !*
    !*******************************************************************
    implicit none
    integer(ip) :: nstr,lflag
    integer(ip) :: istr,decim
    real(rp)    :: digit
    character(len=1) :: string(s_long)
    !
    lflag = 0                                             ! Number by default
    istr  = 1
    !
    do while(istr.le.nstr)
       decim = ichar(string(istr))                         ! decimal value
       if(decim.lt.48.or.decim.gt.57) then                 ! It is not a num.
          if(decim.ne.43 .and.decim.ne.45 .and.  &         ! discard + -
               decim.ne.68 .and.decim.ne.69 .and.  &       ! discard D E
               decim.ne.100.and.decim.ne.101.and.  &       ! discard d e
               decim.ne.46) then                           ! discard .
             lflag = 1
             istr  = nstr
          end if
       end if
       istr = istr+1
    end do
    !
    if(lflag.eq.0) digit = stof(string,nstr)               ! It's a number
    !
    return
  end subroutine decod1
  !
  !
  !
  real(rp) function stof(string,nstr)
    !**************************************************************
    !*
    !*    This routine converts a real/integer number stored in a
    !*    string(1:nstr) into a real(rp)  digit format
    !*
    !**************************************************************
    implicit none
    integer(ip) ::   nstr
    character(len=1) :: string(*)
    !
    integer(ip) :: i,ipos,nsign,esign,nvalu
    integer(ip) :: expo,valu(s_long)
    logical     :: next
    !
    stof = 0.0_rp
    !
    !***  Sing decoding
    !
    ipos = 1
    if(ichar(string(ipos)).eq.43) then         !  + sign
       nsign = 1
       ipos  = ipos + 1
    else if(ichar(string(ipos)).eq.45) then    !  - sign
       nsign = -1
       ipos  = ipos + 1
    else                                       !  no sing (+)
       nsign = 1
       ipos  = ipos
    end if
    !
    !***  Base decoding
    !
    nvalu = 0
    next  = .true.
    do while(next)
       if((ichar(string(ipos)).eq.68 ).or. &         ! D
            (ichar(string(ipos)).eq.69 ).or. &       ! E
            (ichar(string(ipos)).eq.100).or. &       ! d
            (ichar(string(ipos)).eq.101).or. &       ! e
            (ichar(string(ipos)).eq.46 )) then       ! .
          next = .false.
       else
          nvalu = nvalu + 1
          valu(nvalu) = stof1(string(ipos))
          ipos = ipos + 1
          if(ipos.eq.(nstr+1)) then
             next = .false.
             ipos = ipos - 1
          end if
       end if
    end do
    do i = 1,nvalu
       stof = stof + valu(i)*1d1**(nvalu-i)
    end do
    !
    !***  Decimal decoding
    !
    if((ichar(string(ipos)).eq.46   ).and.  &
         ipos  .ne.nstr) then
       ipos = ipos + 1
       nvalu = 0
       next  = .true.
       do while(next)
          if((ichar(string(ipos)).eq.68 ).or. &         ! D
               (ichar(string(ipos)).eq.69 ).or. &       ! E
               (ichar(string(ipos)).eq.100).or. &       ! d
               (ichar(string(ipos)).eq.101)) then       ! e
             next = .false.
          else
             nvalu = nvalu + 1
             valu(nvalu) = stof1(string(ipos))
             ipos = ipos + 1
             if(ipos.eq.(nstr+1)) then
                next = .false.
                ipos = ipos - 1
             end if
          end if
       end do
       do i = 1,nvalu
          stof = stof + valu(i)*1d1**(-i)
       end do
    end if
    !
    !***  Exponent
    !
    if(((ichar(string(ipos)).eq.68 ).or. &         ! D
         (ichar(string(ipos)).eq.69 ).or. &        ! E
         (ichar(string(ipos)).eq.100).or. &        ! d
         (ichar(string(ipos)).eq.101)).and. &      ! e
         ipos  .ne.nstr) then
       ipos = ipos + 1
       if(ichar(string(ipos)).eq.43) then          !  + sign
          esign = 1
          ipos  = ipos + 1
       else if(ichar(string(ipos)).eq.45) then     !  - sign
          esign = -1
          ipos  = ipos + 1
       else                                        !  no sing (+)
          esign = 1
          ipos  = ipos
       end if
       !
       nvalu = 0
       next  = .true.
       do while(next)
          nvalu = nvalu + 1
          valu(nvalu) = stof1(string(ipos))
          ipos = ipos + 1
          if(ipos.eq.(nstr+1)) then
             next = .false.
             ipos = ipos - 1
          end if
       end do
       expo = 0
       do i = 1,nvalu
          expo = expo + valu(i)*10**(nvalu-i)
       end do
       !
       if(esign.eq.1) then
          stof = stof*(10.0_rp**expo)
       else if(esign.eq.-1) then
          stof = stof/(10.0_rp**expo)
       end if
       !
    end if
    !
    stof = nsign*stof
  end function stof
  !
  !
  !
  integer(ip) function stof1(string1)
    !**************************************************************
    !*
    !*    Decodes a character*1 string
    !*
    !**************************************************************
    implicit none
    character(len=1) :: string1
    !
    if(string1.eq.'0') then
       stof1 = 0
    else if(string1.eq.'1') then
       stof1 = 1
    else if(string1.eq.'2') then
       stof1 = 2
    else if(string1.eq.'3') then
       stof1 = 3
    else if(string1.eq.'4') then
       stof1 = 4
    else if(string1.eq.'5') then
       stof1 = 5
    else if(string1.eq.'6') then
       stof1 = 6
    else if(string1.eq.'7') then
       stof1 = 7
    else if(string1.eq.'8') then
       stof1 = 8
    else if(string1.eq.'9') then
       stof1 = 9
    end if
    return
  end function stof1
  !
  !
  !
END MODULE decode
