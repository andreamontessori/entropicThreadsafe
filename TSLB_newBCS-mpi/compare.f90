
 program mio
 
  implicit none
  
  integer, parameter :: db=4 !kind(1.0)
  
    real(kind=db), parameter :: ZERO=real(0.d0,kind=db)
    real(kind=db), parameter :: HALF=real(0.5d0,kind=db)
    real(kind=db), parameter :: ONE=real(1.d0,kind=db)
    real(kind=db), parameter :: TWO=real(2.d0,kind=db)
    real(kind=db), parameter :: THREE=real(3.d0,kind=db)
    real(kind=db), parameter :: FOUR=real(4.d0,kind=db)
    real(kind=db), parameter :: FIVE=real(5.d0,kind=db)
    real(kind=db), parameter :: SIX=real(6.d0,kind=db)
    real(kind=db), parameter :: EIGHT=real(8.d0,kind=db)
    real(kind=db), parameter :: TEN=real(10.d0,kind=db)
    real(kind=db), parameter :: TWELVE=real(12.d0,kind=db)
    real(kind=db), parameter :: FOURTEEN=real(14.d0,kind=db)
    real(kind=db), parameter :: TWENTYFOUR=real(24.d0,kind=db)
  
  integer :: narg,inumchar
  
  integer, parameter :: mxln=120
  character(len=mxln) :: dataFile1,dataFile2,arg,directive
  
  integer :: ndim,nx,ny,nz,nard,i,j,k,l
  real(kind=db), allocatable, dimension(:,:,:,:) :: data1,data2
  real(kind=db) :: mymax
  
        narg = command_argument_count()
        if (narg /= 6) then
          write(6,*) 'error!'
          write(6,*) 'the command line should be'
          write(6,*) '[executable] [file1] [file2] [ndim] [nx] [ny] [nz]'
          write(6,*) 'file1  = name of first file'
          write(6,*) 'file2  = name of second file'
          write(6,*) 'ndim  = 1 for scalar and 3 for velocity fields'
          write(6,*) 'nx  = dimension along x'
          write(6,*) 'ny  = dimension along y'
          write(6,*) 'nz  = dimension along z'
          write(6,*) 'STOP!'
          stop
        endif
  
    
  
        do i = 1, narg
          call getarg(i, arg)
          if(i==1)then
            dataFile1=repeat(' ',mxln)
            call copystring(arg,dataFile1,mxln)
            write(6,*) 'file1  = ',trim(dataFile1)
          elseif(i==2)then
            dataFile2=repeat(' ',mxln)
            call copystring(arg,dataFile2,mxln)
            write(6,*) 'file2  = ',trim(dataFile2)
          elseif(i==3)then
            call copystring(arg,directive,mxln)
             ndim=intstr(directive,mxln,inumchar)
             write(6,*) 'ndim  = ',ndim
          elseif(i==4)then
            call copystring(arg,directive,mxln)
            nx=intstr(directive,mxln,inumchar)
            write(6,*) 'nx  = ',nx
          elseif(i==5)then
            call copystring(arg,directive,mxln)
            ny=intstr(directive,mxln,inumchar)
            write(6,*) 'ny  = ',ny
          elseif(i==6)then
            call copystring(arg,directive,mxln)
            nz=intstr(directive,mxln,inumchar)
            write(6,*) 'nz  = ',nz
          endif
        enddo
        
        allocate(data1(ndim,nx,ny,nz))
        allocate(data2(ndim,nx,ny,nz))
        
        open(unit=345,file=trim(dataFile1), &
         status='old',action='read',access='stream',form='unformatted')
        read(345)data1
        close(345)
        
        open(unit=346,file=trim(dataFile2), &
         status='old',action='read',access='stream',form='unformatted')
        read(346)data2
        close(346)
        
        mymax=real(0.d0,kind=db)
        
        do l=1,ndim
          do k=1,nz
            do j=1,ny
              do i=1,nx
                mymax=max(mymax,abs(data1(l,i,j,k)-data2(l,i,j,k)))
              enddo
            enddo
          enddo
        enddo
        
        write(6,'(g24.16)')mymax
        
   contains
 
   subroutine copystring(oldstring,newstring,lenstring)
 
!***********************************************************************
!     
!     LBsoft subroutine to copy one character string into another
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(in) :: oldstring
  character(len=*), intent(out) :: newstring
  integer, intent(in) :: lenstring
  
  integer :: i
  
  do i=1,lenstring
    newstring(i:i)=oldstring(i:i)
  enddo
  
  return
  
 end subroutine copystring
 
 function intstr(string,lenstring,laststring)
 
!***********************************************************************
!     
!     LBsoft subroutine for extracting integers from a character 
!     string
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(inout) :: string
  integer, intent(in) :: lenstring
  integer, intent(out) :: laststring
  
  integer :: intstr
  
  integer :: j,isn
  character*1, parameter, dimension(0:9) :: & 
   n=(/'0','1','2','3','4','5','6','7','8','9'/)
  logical :: flag,lcount,final
  character*1 :: ksn
  character*1, dimension(lenstring) :: word
  
  do j=1,lenstring
    word(j)=string(j:j)
  enddo
  
  isn=1
  laststring=0
  ksn='+'
  intstr=0
  flag=.false.
  final=.false.
  lcount=.false.
  
  
  do while(laststring<lenstring.and.(.not.final))
    
    laststring=laststring+1
    flag=.false.
    
    do j=0,9
      
      if(n(j)==word(laststring))then
        
        intstr=10*intstr+j
        lcount=.true.
        flag=.true.
        
      endif
    
    enddo
    
    if(lcount.and.(.not.flag))final=.true.
    if(flag .and. ksn=='-')isn=-1
    ksn=word(laststring)
    
  enddo

  intstr=isn*intstr

  do j=laststring,lenstring
    word(j-laststring+1)=word(j)
  enddo
  do j=lenstring-laststring+2,lenstring
    word(j)=' '
  enddo
  
  do j=1,lenstring
    string(j:j)=word(j)
  enddo
  
  return
  
  end function intstr
  
  function dblstr(string,lenstring,laststring)
  
!***********************************************************************
!     
!     LBsoft subroutine for extracting double precisions from a  
!     character string
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(inout) :: string
  integer, intent(in) :: lenstring
  integer, intent(out) :: laststring
  
  real(kind=db) :: dblstr
  
  logical :: flag,ldot,start,final
  integer :: iexp,idum,i,j,fail
  real(kind=db) :: sn,sten,sone

  character*1, parameter, dimension(0:9) :: & 
   n=(/'0','1','2','3','4','5','6','7','8','9'/)
  character*1, parameter :: dot='.'
  character*1, parameter :: d='d'
  character*1, parameter :: e='e'
  
  character*1 :: ksn
  character*1, dimension(lenstring) :: word
  character(len=lenstring) :: work
  
  do j=1,lenstring
    word(j)=string(j:j)
  enddo
  
  laststring=0
  sn= ONE
  ksn='+'
  sten= TEN
  sone= ONE
 
  dblstr = ZERO
  iexp=0
  idum=0
  start=.false.
  ldot=.false.
  final=.false.
  
  do while(laststring<lenstring .and. (.not.final))
    
    laststring=laststring+1
    flag=.false.
    
    do j=0,9
      
      if(n(j)==word(laststring))then
        
        dblstr=sten*dblstr+sone*real(j,kind=db)
        flag=.true.
        start=.true.
      endif
          
    enddo
    
    
    if(dot==word(laststring))then
          
      flag=.true.
      sten= ONE
      ldot=.true.
      start=.true.
          
    endif

    if(flag .and. ksn=='-') sn=- ONE
    if(ldot)sone= real(sone,kind=db)/ TEN
    ksn=word(laststring)
    if(ksn=="D")ksn="d"
    if(ksn=="E")ksn="e"
    
    if(start)then
      if(d==ksn .or. e==ksn)then
        do i=1,lenstring-laststring
          work(i:i)=word(i+laststring)
        enddo
        iexp=intstr(work,lenstring-laststring,idum)
        final=.true.
      endif
      if(.not.flag)final=.true.        
    endif
  enddo
  
  dblstr=sn*dblstr*( TEN ** iexp)
  laststring=laststring+idum
  
  do j=laststring,lenstring
    word(j-laststring+1)=word(j)
  enddo
  do j=lenstring-laststring+2,lenstring
    word(j)=' '
  enddo
  
  do j=1,lenstring
    string(j:j)=word(j)
  enddo
  
  return
  
 end function dblstr
  
 end program
