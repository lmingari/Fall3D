subroutine get_input_npar(fname,block,line,npar,istat,message)
  !**************************************************************************
  !*
  !*    Gets the number of parameters associated to a line
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the file
  !*    character*(*)   block    Block to search
  !*    character*(*)   line     Line block to search
  !*
  !*    OUTPUT:
  !*    integer         npar     Number of parameters
  !*    integer         istat    -1 ERROR  0 OK  1 WARNING
  !*    character       message  Output message with the error description
  !*
  !**************************************************************************
  use kindtype
  use Decode
  implicit none
  !
  character(len=*) :: message,fname,block,line
  integer(ip)      :: npar,istat
  !
  character(len=s_mess) :: mymessage
  character(len=s_long) :: card
  character(len=s_long) :: words(nwormax)
  logical               :: linefound,blockfound
  integer(ip)           :: ilen,nword
  real(rp)              :: param(nparmax),x0,xf,dx
  !
  !***  Initializations
  !
  ilen = LEN(message)
  message(:)  = ' '
  words(:)(:) = ' '
  istat = 0
  !
  !***  Opens the file
  !
  open(90,FILE=TRIM(fname),STATUS='old',ERR=101)
  !
  !***  Search the line
  !
  blockfound = .false.
  linefound  = .false.
  do while(.not.linefound)
     do while(.not.blockfound)
        read(90,'(a256)',END=102) card
        call sdecode(card,words,param,nword,npar)
        if(words(1)(1:LEN_TRIM(block)) == block(1:LEN_TRIM(block))) &
             blockfound=.true.
     end do
     read(90,'(a256)',END=103) card
     call sdecode(card,words,param,nword,npar)
     if(words(1)(1:LEN_TRIM(line)) == line(1:LEN_TRIM(line))) &
          linefound = .true.
  end do
  !
  !***  Line format  FROM x0 TO xf INCREMENT dx
  !***  Calculate npar
  !
  if(TRIM(words(2)) == 'FROM' .or. TRIM(words(2)) == 'from') then
     x0 = param(1)
     xf = param(2)
     if(x0.gt.xf) goto 104
     dx = min(xf-x0,param(3))
     npar = INT((xf-x0)/dx)+1
     if(npar.gt.nparmax) then
        npar = nparmax
        istat = 1
        mymessage = 'get_input_int4: warning. Too big number of parameters'
        message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
     end if
  end if
  !
  !***  Successful end
  !
  close(90)
  return
  !
  !***  List of errors
  !
101 istat = -1
  close(90)
  mymessage = 'get_input_npar : error opening the input file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
102 istat = -1
  close(90)
  mymessage = 'get_input_npar : block '//TRIM(block)// &
       ' not found in the input file'
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
103 istat = -1
  close(90)
  mymessage = 'get_input_npar : line '//TRIM(line)// &
       ' not found in the input file'
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
104 istat = -1
  close(90)
  mymessage = 'get_input_npar: error in line '//line(1:LEN_TRIM(line))
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
end subroutine get_input_npar
