module vars

    implicit none
    integer, parameter :: db=4 !kind(1.0)
    integer :: i,j,k,ll
    integer :: nx,ny,nz,step,stamp,stamp2D,nsteps,ngpus
    integer :: istat,iframe,iframe2D
    
    logical :: lprint,lvtk,lasync,lpbc,lraw
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=db)  :: ts1,ts2,p0g,p1g
    real(kind=db) :: visc_LB,uu,udotc,omega
    real(kind=db) :: tau,one_ov_nu,fx,fy,fz,temp

    real(kind=db) :: fneq1,feq,feqs,us,vs,ws
    
    logical :: periodic(3)
    
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
    
    real(kind=db), parameter :: cssq=1.0_db/3.0_db
    real(kind=db), parameter :: p1dcssq=p1/cssq
    real(kind=db), parameter :: p2dcssq=p2/cssq
    real(kind=db), parameter :: p3dcssq=p3/cssq
    
    integer(kind=4), allocatable,dimension(:,:,:)   :: isfluid
    real(kind=db), allocatable, dimension(:,:,:) :: rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz
    real(kind=db), allocatable, dimension(:,:,:) :: phi,normx,normy,normz
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


endmodule
