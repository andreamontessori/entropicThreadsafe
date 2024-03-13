
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

   real(kind=db), parameter :: infinity = huge(db)

   integer :: narg,inumchar

   integer, parameter :: mxln=120
   character(len=mxln) :: dataFile1,dataFile1b,dataFile2,arg,directive
   logical :: lelittle,lexist
   integer :: ncomp,nx,ny,nz,nard,i,j,k,l,ndataFile1,mystart,myend,iframe
   real(kind=db), allocatable, dimension(:,:,:) :: rho4
   real(kind=db), allocatable, dimension(:,:,:,:) :: vel4
   real(kind=db) :: mymax

   narg = command_argument_count()
   if (narg /= 7) then
      write(6,*) 'error!'
      write(6,*) 'the command line should be'
      write(6,*) '[executable] [filein] [ndim] [nx] [ny] [nz] [start] [end]'
      write(6,*) 'filein = name of first file'
      write(6,*) 'ndim   = 1 for scalar and 3 for velocity fields'
      write(6,*) 'nx     = dimension along x'
      write(6,*) 'ny     = dimension along y'
      write(6,*) 'nz     = dimension along z'
      write(6,*) 'start  = first frame to convert'
      write(6,*) 'end    = last frame to convert'
      write(6,*) 'STOP!'
      stop
   endif



   do i = 1, narg
      call getarg(i, arg)
      if(i==1)then
         dataFile1=repeat(' ',mxln)
         call copystring(arg,dataFile1,mxln)
         write(6,*) 'filein  = ',trim(dataFile1)
      elseif(i==2)then
         call copystring(arg,directive,mxln)
         ncomp=intstr(directive,mxln,inumchar)
         write(6,*) 'ndim  = ',ncomp
         if(ncomp.ne.1 .and. ncomp.ne.3)then
            write(6,*) 'ndim  wrong: ndim can be 1 or 3 !'
            stop
         endif
      elseif(i==3)then
         call copystring(arg,directive,mxln)
         nx=intstr(directive,mxln,inumchar)
         write(6,*) 'nx  = ',nx
      elseif(i==4)then
         call copystring(arg,directive,mxln)
         ny=intstr(directive,mxln,inumchar)
         write(6,*) 'ny  = ',ny
      elseif(i==5)then
         call copystring(arg,directive,mxln)
         nz=intstr(directive,mxln,inumchar)
         write(6,*) 'nz  = ',nz
      elseif(i==6)then
         call copystring(arg,directive,mxln)
         mystart=intstr(directive,mxln,inumchar)
         write(6,*) 'start = ',mystart
      elseif(i==7)then
         call copystring(arg,directive,mxln)
         myend=intstr(directive,mxln,inumchar)
         write(6,*) 'end   = ',myend
      endif
   enddo

   ndataFile1=firstspace(dataFile1,mxln)
   dataFile1b=dataFile1
   dataFile2=dataFile1
   dataFile2(ndataFile1-9:ndataFile1-5)='conv_'
   dataFile2(ndataFile1-4:ndataFile1+1)='######'
   dataFile2(ndataFile1+2:ndataFile1+5)='.vti'
   write(6,*) 'fileout  = ',trim(dataFile2)

   call test_little_endian(lelittle)

   if(ncomp==1)then

      allocate(rho4(nx,ny,nz))

      do iframe=mystart,myend

         dataFile1b(ndataFile1-9:ndataFile1-4)=write_fmtnumb(iframe)
         inquire(file=trim(dataFile1b),exist=lexist)
         if(.not. lexist)then
            write(6,*) 'file ',trim(dataFile1b),' not exist!'
            write(6,*) 'jump to next file!'
            cycle
         endif
         write(6,*) 'reading: ',trim(dataFile1b)
         open(unit=345,file=trim(dataFile1b), &
            status='old',action='read',access='stream',form='unformatted')
         read(345)rho4
         close(345)

         dataFile2(ndataFile1-4:ndataFile1+1)=write_fmtnumb(iframe)
         write(6,*) 'writing: ',trim(dataFile2)
         call write_rho_VTI(0,trim(dataFile2))

      enddo

   elseif(ncomp==3)then

      allocate(vel4(ncomp,nx,ny,nz))

      do iframe=mystart,myend

         dataFile1b(ndataFile1-9:ndataFile1-4)=write_fmtnumb(iframe)
         inquire(file=trim(dataFile1b),exist=lexist)
         if(.not. lexist)then
            write(6,*) 'file ',trim(dataFile1b),' not exist!'
            write(6,*) 'jump to next file!'
            cycle
         endif
         write(6,*) 'reading: ',trim(dataFile1b)
         open(unit=345,file=trim(dataFile1b), &
            status='old',action='read',access='stream',form='unformatted')
         read(345)vel4
         close(345)

         dataFile2(ndataFile1-4:ndataFile1+1)=write_fmtnumb(iframe)
         write(6,*) 'writing: ',trim(dataFile2)
         call write_vel_VTI(0,trim(dataFile2))

      enddo

   endif

   write(6,*) 'program correctly closed!'

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

   function firstspace(oldstring,lenstring)

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
      integer, intent(in) :: lenstring
      integer :: firstspace
      integer :: i

      firstspace=0
      do i=1,lenstring
         if(oldstring(i:i)==' ')exit
         firstspace=firstspace+1
      enddo

      return

   end function firstspace

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

   function space_fmtnumb(inum)

      !***********************************************************************
      !
      !     LBsoft function for returning the string of six characters
      !     with integer digits and leading spaces to the left
      !     originally written in JETSPIN by M. Lauricella et al.
      !
      !     licensed under Open Software License v. 3.0 (OSL-3.0)
      !     author: M. Lauricella
      !     last modification October 2019
      !
      !***********************************************************************

      implicit none

      integer,intent(in) :: inum
      character(len=6) :: space_fmtnumb
      integer :: numdigit,irest
      real(kind=8) :: tmp
      character(len=22) :: cnumberlabel

      numdigit=dimenumb(inum)
      irest=6-numdigit
      if(irest>0)then
         write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
         write(space_fmtnumb,fmt=cnumberlabel)repeat(' ',irest),inum
      else
         write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
         write(space_fmtnumb,fmt=cnumberlabel)inum
      endif

      return

   end function space_fmtnumb

   function space_fmtnumb12(inum)

      !***********************************************************************
      !
      !     LBsoft function for returning the string of six characters
      !     with integer digits and leading TWELVE spaces to the left
      !     originally written in JETSPIN by M. Lauricella et al.
      !
      !     licensed under Open Software License v. 3.0 (OSL-3.0)
      !     author: M. Lauricella
      !     last modification October 2019
      !
      !***********************************************************************

      implicit none

      integer,intent(in) :: inum
      character(len=12) :: space_fmtnumb12
      integer :: numdigit,irest
      real(kind=8) :: tmp
      character(len=22) :: cnumberlabel

      numdigit=dimenumb(inum)
      irest=12-numdigit
      if(irest>0)then
         write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
         write(space_fmtnumb12,fmt=cnumberlabel)repeat(' ',irest),inum
      else
         write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
         write(space_fmtnumb12,fmt=cnumberlabel)inum
      endif

      return

   end function space_fmtnumb12

   function dimenumb(inum)

      !***********************************************************************
      !
      !     LBsoft function for returning the number of digits
      !     of an integer number
      !     originally written in JETSPIN by M. Lauricella et al.
      !
      !     licensed under the 3-Clause BSD License (BSD-3-Clause)
      !     author: M. Lauricella
      !     last modification July 2018
      !
      !***********************************************************************

      implicit none

      integer,intent(in) :: inum
      integer :: dimenumb
      integer :: i
      real(kind=db) :: tmp

      i=1
      tmp=real(inum,kind=db)
      do
         if(tmp< 10.0_db )exit
         i=i+1
         tmp=tmp/ 10.0_db
      enddo

      dimenumb=i

      return

   end function dimenumb

   function write_fmtnumb(inum)

      !***********************************************************************
      !
      !     LBsoft function for returning the string of six characters
      !     with integer digits and leading zeros to the left
      !     originally written in JETSPIN by M. Lauricella et al.
      !
      !     licensed under the 3-Clause BSD License (BSD-3-Clause)
      !     author: M. Lauricella
      !     last modification July 2018
      !
      !***********************************************************************

      implicit none

      integer,intent(in) :: inum
      character(len=6) :: write_fmtnumb
      integer :: numdigit,irest
      !real*8 :: tmp
      character(len=22) :: cnumberlabel

      numdigit=dimenumb(inum)
      irest=6-numdigit
      if(irest>0)then
         write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
         write(write_fmtnumb,fmt=cnumberlabel)repeat('0',irest),inum
      else
         write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
         write(write_fmtnumb,fmt=cnumberlabel)inum
      endif

      return
   end function write_fmtnumb

   function write_fmtnumb2(inum)

      !***********************************************************************
      !
      !     LBsoft function for returning the string of six characters
      !     with integer digits and leading zeros to the left
      !     originally written in JETSPIN by M. Lauricella et al.
      !
      !     licensed under the 3-Clause BSD License (BSD-3-Clause)
      !     author: M. Lauricella
      !     last modification July 2018
      !
      !***********************************************************************

      implicit none

      integer,intent(in) :: inum
      character(len=2) :: write_fmtnumb2
      integer :: numdigit,irest
      !real*8 :: tmp
      character(len=22) :: cnumberlabel

      numdigit=dimenumb(inum)
      irest=2-numdigit
      if(irest>0)then
         write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
         write(write_fmtnumb2,fmt=cnumberlabel)repeat('0',irest),inum
      else
         write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
         write(write_fmtnumb2,fmt=cnumberlabel)inum
      endif

      return
   end function write_fmtnumb2

   subroutine write_rho_VTI(myframe,fnameFull)
      implicit none

      integer, intent(in) :: myframe

      logical, parameter :: textual=.false.
      character(len=120) :: extent
      integer i,j,k, iotest,iotest2
      integer(4)         :: length


      character(len=*) :: fnameFull
      character(len=*), parameter :: fname='dens'


      iotest = 55

      open(unit=iotest,file=trim(fnameFull),status='replace',action='write')

      extent =  trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(nx)) // ' ' &
         // trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(ny)) // ' ' &
         // trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(nz))

      if(lelittle)then
         write(iotest,*) '<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" >'
      else
         write(iotest,*) '<VTKFile type="ImageData" version="1.0" byte_order="BigEndian" >'
      endif
      write(iotest,*) ' <ImageData WholeExtent="' // trim(extent) // '" >'
      write(iotest,*) ' <Piece Extent="' // trim(extent) // '">'
      write(iotest,*) '   <PointData>'

      if (textual) then
         write(iotest,*) '    <DataArray type="Float32" Name="',trim(fname),'" format="ascii" >'

         do k=1,nz
            do j=1,ny
               do i=1,nx
                  write(iotest,fmt='("     ", F20.8)') rho4(i,j,k)
               enddo
            enddo
         enddo

         write(iotest,*) '    </DataArray>'
         write(iotest,*) '   </PointData>'
         write(iotest,*) ' </Piece>'
         write(iotest,*) ' </ImageData>'
      else
         write(iotest,*) '    <DataArray type="Float32" Name="',trim(fname),'" format="appended" offset="0" />'
         write(iotest,*) '   </PointData>'
         write(iotest,*) ' </Piece>'
         write(iotest,*) ' </ImageData>'
         write(iotest,*) ' <AppendedData encoding="raw">'
         close(iotest)

         open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
            access='stream',form='unformatted',action='write')
         write(iotest) '_'
         length = 4*nx*ny*nz
         write(iotest) length, rho4(1:nx,1:ny,1:nz)
         close(iotest)

         open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
            action='write')
         write(iotest,*) ''
         write(iotest,*) ' </AppendedData>'
      endif

      write(iotest,*) '</VTKFile >'
      close(iotest)



   end subroutine write_rho_VTI

   subroutine write_vel_VTI(myframe,fnameFull)
      implicit none

      integer, intent(in) :: myframe

      logical, parameter :: textual=.false.
      character(len=120) :: extent
      integer i,j,k, iotest
      integer(4)         :: length

      character(len=*) :: fnameFull


      iotest = 56

      open(unit=iotest,file=trim(fnameFull),status='replace',action='write')

      extent =  trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(nx)) // ' ' &
         // trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(ny)) // ' ' &
         // trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(nz))

      if(lelittle)then
         write(iotest,*) '<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" >'
      else
         write(iotest,*) '<VTKFile type="ImageData" version="1.0" byte_order="BigEndian" >'
      endif
      write(iotest,*) ' <ImageData WholeExtent="' // trim(extent) // '" >'
      write(iotest,*) ' <Piece Extent="' // trim(extent) // '">'
      write(iotest,*) '   <PointData>'

      if (textual) then
         write(iotest,*) '    <DataArray type="Float32" Name="vel" NumberOfComponents="3" format="ascii" >'

         do k=1,nz
            do j=1,ny
               do i=1,nx
                  write(iotest,fmt='("     ", F20.8)') vel4(1,i,j,k),vel4(2,i,j,k),vel4(3,i,j,k)
               enddo
            enddo
         enddo

         write(iotest,*) '    </DataArray>'
         write(iotest,*) '   </PointData>'
         write(iotest,*) ' </Piece>'
         write(iotest,*) ' </ImageData>'
      else
         write(iotest,*) '    <DataArray type="Float32" Name="vel" NumberOfComponents="3" format="appended" offset="0" />'
         write(iotest,*) '   </PointData>'
         write(iotest,*) ' </Piece>'
         write(iotest,*) ' </ImageData>'
         write(iotest,*) ' <AppendedData encoding="raw">'
         close(iotest)

         open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
            access='stream',form='unformatted',action='write')
         write(iotest) '_'
         length = 4*nx*ny*nz*3
         write(iotest) length, vel4(1:3, 1:nx,1:ny,1:nz)
         close(iotest)

         open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
            action='write')
         write(iotest,*) ''
         write(iotest,*) ' </AppendedData>'
      endif

      write(iotest,*) '</VTKFile >'
      close(iotest)



   end subroutine write_vel_VTI

   subroutine test_little_endian(ltest)

!***********************************************************************
!
!     LBsoft subroutine for checking if the computing architecture
!     is working in little-endian or big-endian
!
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!
!***********************************************************************

      implicit none
      integer, parameter :: ik1 = selected_int_kind(2)
      integer, parameter :: ik4 = selected_int_kind(9)

      logical, intent(out) :: ltest

      if(btest(transfer(int((/1,0,0,0/),ik1),1_ik4),0)) then
         !it is little endian
         ltest=.true.
      else
         !it is big endian
         ltest=.false.
      end if

      return

   end subroutine test_little_endian

end program
