module vars

    implicit none
    integer, parameter :: db=4 !kind(1.0)
    integer :: i,j,ll,l,dumm
    integer :: nx,ny,step,stamp,nlinks,nsteps,stab_points,ngpus
    integer :: istat,iframe
    integer, parameter :: nz=1
    
    logical :: lprint=.true.
    logical :: lvtk=.true.
    logical :: lasync=.false.
    logical :: lpbc=.true.
    
    real(kind=db),parameter :: pi_greek=3.14159265359793234626433
    real(kind=db),parameter :: entropic_threshold=0.001
    
    real(kind=4)  :: ts1,ts2 
    real(kind=db) :: visc_LB,omega,feq,udotc,uu,fneq1
    real(kind=db) :: qxx,qyy,qxy5_7,qxy6_8,pi2cssq0,a2xxy,a2yyx
    real(kind=db) :: tau,one_ov_nu,cssq,cssq2,cssq3,fx,fy,temp,dummy
    real(kind=db) :: deltasmed
    
    integer(kind=4), allocatable,  dimension(:,:)   :: isfluid
    
    real(kind=db), allocatable, dimension(:)     :: p
    real(kind=db), allocatable, dimension(:,:) :: rho,u,v,pxx,pyy,pxy,alpha
    real(kind=db), allocatable, dimension(:,:,:) :: f
    real(kind=db) :: mymemory,totmemory

endmodule