module vars
    !$if _OPENACC
    use openacc
    !$endif
    implicit none
    integer, parameter :: db=4 !kind(1.0)
    
    integer :: devNum
    integer(acc_device_kind) :: devType
    integer :: nprocs,myrank,lbecomm,localcomm
    integer:: rear_task, front_task
    integer:: left_task, right_task
    integer:: down_task, up_task
    integer:: xyplane, xzplane, yzplane, myxrank, yzcomm
    
    integer :: mydev, ndev 
    integer :: file_offset    
    integer :: proc_x,proc_y,proc_z
    integer :: mem_stop
    logical :: rreorder 
    integer, parameter::  mpid=2      ! mpi dimension
    logical :: periodic(mpid) 
    integer :: prgrid(mpid)
    integer:: mpicoords(mpid)
    integer:: up(2),down(2),left(2)
    integer:: front(2),rear(2),right(2)
    
    integer :: i,j,k,ll
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

endmodule
