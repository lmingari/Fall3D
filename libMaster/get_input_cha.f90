subroutine get_input_cha(fname,block,line,value,nval,istat,message)
  !**************************************************************************
  !*
  !*    Gets nval character inputs from the file fname
  !*
  !*    NOTE: words are converted to UPPER case
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the file
  !*    character*(*)   block    Block to search
  !*    character*(*)   line     Line block to search
  !*    integer         nval     Number of integers to read
  !*
  !*    OUTPUT:
  !*    integer          istat    -1 ERROR  0 OK  1 WARNING
  !*    character        value    Values of the nval integers read
  !*    character        message  Output message with the error description
  !*
  !*    CALLED BY: user
  !*
  !**************************************************************************
  use kindtype
  use decode
  implicit none
  integer(ip)           :: nval,istat
  character(len=*)      :: message,fname,block,line
  character(len=*)      :: value(nval)
  !
  character(len=s_mess) :: mymessage
  character(len=s_long) :: card
  character(len=s_long) :: words(nwormax)
  logical               :: linefound,blockfound
  integer(ip)           :: ilen,nword,npar,ival,j
  real(rp)              :: param(nparmax)
  !
  !***  Initializations
  !
  ilen = LEN(message)
  message(:)  = ' '
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
        if(words(1)(1:LEN_TRIM(block)).eq.block(1:LEN_TRIM(block))) &
             blockfound=.true.
     end do
     read(90,'(a256)',END=103) card
     call sdecode(card,words,param,nword,npar)
     if(words(1)(1:LEN_TRIM(line)).eq.line(1:LEN_TRIM(line))) &
          linefound = .true.
  end do
  !
  if((nword-1).lt.nval) goto 104
  !
  do ival = 1,nval
     value(ival)(1:LEN_TRIM(words(ival+1))) = &
          words(ival+1)(1:LEN_TRIM(words(ival+1)))
     !
     !***     Fill the rest with ' '
     !
     do j=LEN_TRIM(words(ival+1))+1,LEN(value(ival))
        value(ival)(j:j)=' '
     end do
  end do
  !
  !***  Convert to upper case
  !
  do ival = 1,nval
     call upcase(value(ival))
  end do
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
  mymessage = 'get_input_cha: error opening the input file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
102 istat = -1
  close(90)
  mymessage = 'get_input_cha: block '//TRIM(block)// &
       ' not found in the input file'
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
103 istat = -1
  close(90)
  mymessage = 'get_input_cha: line '//TRIM(line)// &
       ' not found in the input file'
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
104 istat = -1
  close(90)
  mymessage = 'get_input_cha: too few parameters in line '//TRIM(line)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
end subroutine get_input_cha
