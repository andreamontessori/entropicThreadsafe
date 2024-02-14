module vars

    implicit none
    integer, parameter :: db=4 !kind(1.0)
    integer :: i,j,k,ll
    integer :: nx,ny,nz,step,stamp,stamp2D,nlinks,nsteps,ngpus
    integer :: istat,iframe,iframe2D
    
    logical :: lprint,lvtk,lasync,lpbc
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=db)  :: ts1,ts2,p0,p1,p2,p3,p0g,p1g
    real(kind=db) :: visc_LB,uu,udotc,omega
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,fz,temp

    real(kind=db) :: fneq1,feq,feqs,us,vs,ws
    
    logical :: periodic(3)

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

endmodule