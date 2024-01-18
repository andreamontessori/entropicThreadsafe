module vars
    !$if _OPENACC
    use openacc
    !$endif
    implicit none
    integer, parameter :: db=4 !kind(1.0)
    
    
    integer :: devNum
    integer(acc_device_kind) :: devType
    integer :: nprocs,myrank,lbecomm,localcomm
    integer:: xyplane, xzplane, yzplane, myxrank, yzcomm
    
    integer :: mydev, ndev 
    integer :: file_offset    
    integer :: proc_x,proc_y,proc_z
    integer :: mem_stop
    logical :: rreorder 
    integer, parameter::  mpid=2      ! mpi dimension
    logical :: periodic(mpid) 
    integer :: prgrid(mpid)
    integer:: coords(mpid)
    integer:: up(mpid),down(mpid),left(mpid)
    integer:: front(mpid),rear(mpid),right(mpid)
    integer, allocatable, dimension(:) :: yinidom,yfindom
    integer, allocatable, dimension(:) :: zinidom,zfindom
    
    !integer :: right_dest_x,left_source_x  
	!integer :: left_dest_x,right_source_x
    integer :: up_dest_y,down_source_y
	integer :: down_dest_y, up_source_y
	integer :: front_dest_z,rear_source_z
	integer :: rear_dest_z,front_source_z
	
	integer :: myoffset(3),lsizes(3),start_idx(3),gsizes(3),start_p(3)
	
	
	
    
    integer :: i,j,k,ll
    integer :: gi,gj,gk
    integer :: nx,ny,nz,step,stamp,stamp2D,nlinks,nsteps,ngpus
    integer :: lx,ly,lz
    integer :: istat,iframe,iframe2D
    
    logical :: lprint,lvtk,lasync,lpbc
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=db)  :: ts1,ts2,p0,p1,p2,p3,p1dcssq,p2dcssq,p3dcssq
    real(kind=db) :: visc_LB,uu,udotc,omega
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,fz,temp

    real(kind=db) :: fneq1,feq
    
    integer(kind=4), allocatable,dimension(:,:,:)   :: isfluid
    real(kind=db), allocatable, dimension(:,:,:) :: rho,rhoA,rhoB,psi,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
    real(kind=db), allocatable, dimension(:,:,:,:) :: f,g
    real(kind=db) :: mymemory,totmemory
	real(kind=db) :: uwall
	real(kind=4) :: rrx,rry,rrz

    !****************************print vars**************************************!
    
    integer, parameter :: mxln=120
    character(len=8), allocatable, dimension(:) :: namevarvtk
    character(len=500), allocatable, dimension(:) :: headervtk
    character(len=30), allocatable, dimension(:) :: footervtk
    integer, allocatable, dimension(:) :: ndimvtk
    integer, allocatable, dimension(:) :: vtkoffset
    integer, allocatable, dimension(:) :: ndatavtk
    integer, allocatable, dimension(:) :: nheadervtk
    integer :: nfilevtk
    integer, allocatable, dimension(:) :: varlistvtk
    character :: delimiter
    character(len=*), parameter :: filenamevtk='out'
    
    real(kind=4), allocatable, dimension(:,:,:) :: rhoprint
    real(kind=4), allocatable, dimension(:,:,:,:) :: velprint
    logical :: lelittle
    character(len=mxln) :: dir_out
    character(len=mxln) :: extentvtk
    character(len=mxln) :: sevt1,sevt2
    character(len=1), allocatable, dimension(:) :: head1,head2
    
    
    contains
    
    subroutine string_char(mychar,nstring,mystring)
    
      implicit none
      
      integer :: i
      character(1), allocatable, dimension(:) :: mychar
      integer, intent(in) :: nstring
      character(len=*), intent(in) :: mystring
      
      allocate(mychar(nstring))
      
      do i=1,nstring
        mychar(i)=mystring(i:i)
      enddo
      
    end subroutine string_char
  
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

endmodule
