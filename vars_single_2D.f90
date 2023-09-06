module vars

    implicit none
    integer, parameter :: db=4 !kind(1.0)
    integer :: i,j,ll,l,dumm
    integer :: nx,ny,step,stamp,nlinks,nsteps,ngpus
    integer :: istat,iframe
    integer, parameter :: nz=1
    
    logical :: lprint=.true.
    logical :: lvtk=.true.
    logical :: lasync=.false.
    logical :: lpbc=.true.
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    
    real(kind=4)  :: ts1,ts2 
    real(kind=db) :: visc_LB,omega,feq,one_pls_usq,one_pls_vsq
    real(kind=db) :: fneq1
    real(kind=db) :: qxx,qyy,qxy5_7,qxy6_8,pi2cssq1,pi2cssq2,pi2cssq0
    real(kind=db) :: tau,one_ov_nu,cssq,fx,fy,temp,dummy
    
    integer(kind=4), allocatable,  dimension(:,:)   :: isfluid
    
    real(kind=db), allocatable, dimension(:)     :: p
    real(kind=db), allocatable, dimension(:,:) :: rho,u,v,pxx,pyy,pxy
    real(kind=db), allocatable, dimension(:,:,:) :: f
    real(kind=db) :: mymemory,totmemory

endmodule