module vars
#ifdef _OPENACC
    use openacc
#endif
    implicit none
    
#if PRC==4
    integer, parameter :: db=4 !kind(1.0)
#elif PRC==8
    integer, parameter :: db=8 !kind(1.0)
#else
//#error "ERROR in specifying PRC"
#endif    
    integer, parameter :: isf=1 !kind(1.0)

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

    integer :: i,j,k,ll,l
    integer :: gi,gj,gk
    integer :: nx,ny,nz,step,stamp,stamp2D,nsteps,ngpus
    integer :: lx,ly,lz
    integer :: istat,iframe,iframe2D
    
    logical :: lprint,lvtk,lasync,lpbc,lraw
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=db)  :: ts1,ts2,p0g,p1g
    real(kind=db) :: visc_LB,uu,udotc,omega,radius
    real(kind=db) :: tau,one_ov_nu,fx,fy,fz,temp,dumpYN

    real(kind=db) :: fneq1,feq,feqs,us,vs,ws
    
    integer, parameter :: nlinks=26  
    !lattice vectors
    integer, dimension(0:nlinks), parameter :: &
        !       0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26
        ex=   (/0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1/)
    integer, dimension(0:nlinks), parameter:: &
        ey=  (/ 0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1,-1, 1,-1, 1/)
    integer, dimension(0:nlinks), parameter:: &
        ez=  (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1/)
    integer, dimension(0:nlinks), parameter:: &
        opp =(/ 0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17,20,19,22,21,24,23,26,25/)
    
    real(kind=db), parameter :: p0=real(8.d0/27.d0,kind=db)
    real(kind=db), parameter :: p1=real(2.d0/27.d0,kind=db)
    real(kind=db), parameter :: p2=real(1.d0/54.d0,kind=db)
    real(kind=db), parameter :: p3=real(1.d0/216.d0,kind=db)
    real(kind=db), dimension(0:nlinks), parameter :: &
     p=(/p0,p1,p1,p1,p1,p1,p1,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p3,p3,p3,p3,p3,p3,p3,p3/)
    
    real(kind=db), parameter :: cssq=real(1.d0/3.d0,kind=db)
    real(kind=db), parameter :: p1dcssq=p1/cssq
    real(kind=db), parameter :: p2dcssq=p2/cssq
    real(kind=db), parameter :: p3dcssq=p3/cssq

    real(kind=db), dimension(0:nlinks), parameter :: dex=real(ex,kind=db)
    real(kind=db), dimension(0:nlinks), parameter :: dey=real(ey,kind=db)
    real(kind=db), dimension(0:nlinks), parameter :: dez=real(ez,kind=db)
    
    integer(kind=isf), allocatable,dimension(:,:,:)   :: isfluid
    real(kind=db), allocatable, dimension(:,:,:) :: rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
    real(kind=db), allocatable, dimension(:,:,:) :: phi,normx,normy,normz
    real(kind=db), allocatable, dimension(:,:,:,:) :: f,g
    real(kind=db) :: mymemory,totmemory
	real(kind=db) :: uwall
	real(kind=db) :: rrx,rry,rrz

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
    character(len=mxln) :: sevt1,sevt2,arg,directive
    character(len=1), allocatable, dimension(:) :: head1,head2
    
  contains
  
  pure function gauss_noseeded(i,j,k,l)
    !$acc routine gang
    !questo sopra serve per dire ad openacc di fare una copia sul device e gang significa che può essere chiamata da più threads/vector/worker/gangbang indipendentemente ad cazzum
    implicit none
    integer, intent(in) :: i,j,k,l
    integer :: kk
    real(kind=db) :: gauss_noseeded
    real(kind=db) :: dtemp1,dtemp2
    real(kind=db), parameter :: mylimit=real(1.d-30,kind=db)
    real(kind=db), parameter :: FIFTY = real(50.d0,kind=db)
    real(kind=db), parameter :: ONE = real(1.d0,kind=db)
    real(kind=db), parameter :: TWO = real(2.d0,kind=db)
    real(kind=db), parameter :: Pi = real(3.1415926535897932384626433832795028841971d0,kind=db)
    
    dtemp1=rand_noseeded(i,j,k,l)
    kk=nint(dtemp1*FIFTY)
    dtemp2=rand_noseeded(i,j,k,l+kk)
    
    dtemp1=dtemp1*(ONE-mylimit)+mylimit
    
  ! Box-Muller transformation
    gauss_noseeded=sqrt(- TWO *log(dtemp1))*cos(TWO*pi*dtemp2)
   end function gauss_noseeded
  
  pure function rand_noseeded(i,j,k,l)
    !$acc routine gang
    !questo sopra serve per dire ad openacc di fare una copia sul device e gang significa che può essere chiamata da più threads/vector/worker/gangbang indipendentemente ad cazzum
    implicit none
    integer, intent(in) :: i,j,k,l
    integer :: isub,jsub,ksub,lsub,msub
    integer ::ii,jj
    real(4) :: s,t,u33,u97,csub,uni
    real(kind=db) :: rand_noseeded
    
    real(4), parameter :: c =  362436.0/16777216.0
    real(4), parameter :: cd= 7654321.0/16777216.0
    real(4), parameter :: cm=16777213.0/16777216.0
  
    
  ! initial values of i,j,k must be in range 1 to 178 (not all 1)
  ! initial value of l must be in range 0 to 168.      
    isub=mod(i,178)
    isub=isub+1
    
    jsub=mod(j,178)
    jsub=jsub+1
    
    ksub=mod(k,178)
    ksub=ksub+1
    
    lsub=mod(l,169)
    
  ! initialization on fly  
    ii=97
    s=0.0
    t=0.5
    do jj=1,24
      msub=mod(mod(isub*jsub,179)*ksub,179)
      isub=jsub
      jsub=ksub
      ksub=msub
      lsub=mod(53*lsub+1,169)
      if(mod(lsub*msub,64).ge.32)s=s+t
      t=0.5*t
    enddo
    u97=s
    
    ii=33
    s=0.0
    t=0.5
    do jj=1,24
      msub=mod(mod(isub*jsub,179)*ksub,179)
      isub=jsub
      jsub=ksub
      ksub=msub
      lsub=mod(53*lsub+1,169)
      if(mod(lsub*msub,64).ge.32)s=s+t
      t=0.5*t
    enddo
    u33=s
    uni=u97-u33
    if (uni.lt.0.0) uni = uni + 1.0
    csub = c-cd
    if (csub.lt.0.0) csub = csub + cm
    uni = uni-csub
    if (uni.lt.0.0) uni = uni+1.0
    rand_noseeded = real(uni,kind=db)
   end function rand_noseeded

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

endmodule
