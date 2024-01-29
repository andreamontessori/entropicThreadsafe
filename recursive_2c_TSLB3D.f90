program recursiveTSLB3D
    !$if _OPENACC
    use openacc
    !$endif
    use prints
    use vars
    use bcs3D
    
    implicit none
    logical :: dumpYN
    integer :: dumpstep
    real*4 :: grad_fix,grad_fiy,grad_fiz,mod_grad,curvature,sigma,sharp_c,ffx_s,ffy_s,ffz_s,indicator
    !$if _OPENACC
    integer :: devNum
    integer(acc_device_kind) :: devType
    devType = acc_get_device_type()
    devNum=acc_get_device_num(devType)
    !$endif

    nlinks=26 !pari!
    tau=0.8_db
    cssq=1.0_db/3.0_db
    visc_LB=cssq*(tau-0.5_db)
    one_ov_nu=1.0_db/visc_LB
    dumpYN=0
    !interface
    sigma=0.001
    sharp_c=0.2*3


#ifdef _OPENACC
    ngpus=acc_get_num_devices(acc_device_nvidia)
#else
    ngpus=0
#endif

    !*******************************user parameters and allocations**************************
        nx=128
        ny=128
        nz=128
        nsteps=50000
        stamp=2000000
        stamp2D=5000
        dumpstep=100000000
        fx=0.0_db*10.0**(-7)
        fy=0.0_db*10.0**(-5)
        fz=0.0_db*10.0**(-5)
		    uwall=0.05
        lprint=.false.
        lvtk=.false.
        lasync=.false.
        lpbc=.false.
        
        allocate(f(0:nx+1,0:ny+1,0:nz+1,0:nlinks),g(0:nx+1,0:ny+1,0:nz+1,0:6))
        allocate(rho(1:nx,1:ny,1:nz),u(1:nx,1:ny,1:nz),v(1:nx,1:ny,1:nz),w(1:nx,1:ny,1:nz))
		    allocate(normx(1:nx,1:ny,1:nz),normy(1:nx,1:ny,1:nz),normz(1:nx,1:ny,1:nz),phi(1:nx,1:ny,1:nz))
        allocate(pxx(1:nx,1:ny,1:nz),pxy(1:nx,1:ny,1:nz),pxz(1:nx,1:ny,1:nz),pyy(1:nx,1:ny,1:nz))
        allocate(pyz(1:nx,1:ny,1:nz),pzz(1:nx,1:ny,1:nz))
        allocate(isfluid(1:nx,1:ny,1:nz))
        
        if(lprint)then
          allocate(rhoprint(1:nx,1:ny,1:nz))
          allocate(velprint(1:3,1:nx,1:ny,1:nz))
          rhoprint(1:nx,1:ny,1:nz)=0.0
          velprint(1:3,1:nx,1:ny,1:nz)=0.0
        endif
        !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
        !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
        !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)

        p0=(8.0_db/27.0_db)
        p1=(2.0_db/27.0_db)
        p2=(1.0_db/54.0_db)
        p3=(1.0_db/216.0_db)
        p0g=1./4.
        p1g=1.0/8.0
        omega=1.0_db/tau
    !*****************************************geometry************************
        isfluid=1
        isfluid(1,:,:)=0 !left
        isfluid(nx,:,:)=0 !right
        isfluid(:,1,:)=0 !front 
        isfluid(:,ny,:)=0 !rear
        isfluid(:,:,1)=0 !bottom
        isfluid(:,:,nz)=0 !top
    !*************************************initial conditions ************************    
        u=0.0_db
        v=0.0_db
        w=0.0_db
        rho=1.0_db  !tot dens
		phi=0.0
		normx=0
		normy=0
		normz=0
        !do ll=0,nlinks
        if(dumpYN.eq.0)then
            do k=1,nz
                  do j=1,ny
                      do i=1,nx
                          if(isfluid(i,j,k).eq.1)then
                              !0
                              if((i-nx/2.)**2 + (j-ny/2.)**2 + (k-nz/2.)**2< 25**2)then
                                  phi(i,j,k)=1
                              endif
                              feq=(-4*rho(i,j,k)*(-2 + 3*u(i,j,k)**2 + 3*v(i,j,k)**2 + 3*w(i,j,k)**2))/27.
                              
                              f(i,j,k,0)=feq 
                              
                              !1
                              
                              feq=(rho(i,j,k)*(2 + 6*u(i,j,k) + 6*u(i,j,k)**2 - 3*v(i,j,k)**2 - 9*u(i,j,k)*v(i,j,k)**2 - 3*(1 + 3*u(i,j,k))*w(i,j,k)**2))/27.
                              
                              f(i,j,k,1)=feq 
                              
                              !2
                              feq=(rho(i,j,k)*(2 - 3*v(i,j,k)**2 - 3*w(i,j,k)**2 + 3*u(i,j,k)*(-2 + 2*u(i,j,k) + 3*v(i,j,k)**2 + 3*w(i,j,k)**2)))/27.
                              
                              f(i,j,k,2)=feq 
                              
                              !3
                              
                              feq=(rho(i,j,k)*(2 - 3*u(i,j,k)**2*(1 + 3*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(2 + 2*v(i,j,k) - 3*w(i,j,k)**2)))/27.
                              
                              f(i,j,k,3)=feq
                              
                              !4
                              feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(-2 + 2*v(i,j,k) + 3*w(i,j,k)**2)))/27.
                              
                              f(i,j,k,4)=feq
                              
                              !7
                              
                              feq=(2*rho(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) + u(i,j,k)*(3 + 9*v(i,j,k)*(1 + v(i,j,k)))) - 3*rho(i,j,k)*(1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                              
                              f(i,j,k,7)=feq 
                              
                              !8
                              feq=(2*rho(i,j,k)*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(-3 - 9*(-1 + v(i,j,k))*v(i,j,k))) + 3*rho(i,j,k)*(-1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                              
                              f(i,j,k,8)=feq 
                              
                              !10
                              
                              feq=(2*rho(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)))) + 3*rho(i,j,k)*(-1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                              
                              f(i,j,k,10)=feq
                              
                              !9
                              feq=(2*rho(i,j,k)*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(3 + 9*(-1 + v(i,j,k))*v(i,j,k))) - 3*rho(i,j,k)*(1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                              
                              f(i,j,k,9)=feq

                              !5
                              
                              feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k))))/27.
                              
                              f(i,j,k,5)=feq
                              
                              !6
                              feq=(rho(i,j,k)*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                              
                              f(i,j,k,6)=feq

                              !15
                              
                              feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                              
                              f(i,j,k,15)=feq
                              
                              !16
                              feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                              
                              f(i,j,k,16)=feq

                              !17
                              
                              feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                              
                              f(i,j,k,17)=feq
                              
                              !18
                              feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                              
                              f(i,j,k,18)=feq

                              !11
                              
                              feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                              
                              f(i,j,k,11)=feq
                              
                              !12
                              feq=(rho(i,j,k)*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                              
                              f(i,j,k,12)=feq

                              !13
                              
                              feq=(rho(i,j,k)*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                              
                              f(i,j,k,13)=feq
                              
                              !14
                              feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                              
                              f(i,j,k,14)=feq

                              !19
                              feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                              
                              f(i,j,k,19)=feq 
                              
                              !20
                              feq=(rho(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                              
                              f(i,j,k,20)=feq

                              !21
                              feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                              
                              f(i,j,k,21)=feq
                              
                              !22
                              feq=(rho(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                              
                              f(i,j,k,22)=feq

                              !23
                              feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                              
                              f(i,j,k,23)=feq

                              !24
                              feq=(rho(i,j,k)*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                              
                              f(i,j,k,24)=feq

                              !25
                              feq=(rho(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                              
                              f(i,j,k,25)=feq 

                              !26
                              feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                              
                              f(i,j,k,26)=feq 
							  !*******************************************************
                              !g0
                              feq=p0g*phi(i,j,k)
                              g(i,j,k,0)=feq
                              
                              !g1
                              
                              feq=p1g*phi(i,j,k)*(1 + 3*u(i,j,k))
                              g(i+1,j,k,1)=feq
                              
                              !g2
                              feq=p1g*phi(i,j,k)*(1 - 3*u(i,j,k))
                              g(i-1,j,k,2)=feq
                              
                              !g3
                              
                              feq=p1g*phi(i,j,k)*(1 + 3*v(i,j,k))
                              g(i,j+1,k,3)=feq
                              !g4
                              
                              feq=p1g*phi(i,j,k)*(1 - 3*v(i,j,k))
                              g(i,j-1,k,4)=feq
                              !g5
                              
                              feq=p1g*phi(i,j,k)*(1 + 3*w(i,j,k))
                              g(i,j,k+1,5)=feq
                              !g6
                              feq=p1g*phi(i,j,k)*(1 - 3*w(i,j,k))
                              g(i,j,k-1,6)=feq
                          endif
                      enddo
                  enddo
              enddo
          elseif(dumpYN.eq.1)then
              call read_distros_1c_3d
          endif
    
    !*************************************check data ************************ 
        write(6,*) '*******************LB data*****************'
        write(6,*) 'tau',tau
        write(6,*) 'omega',omega
        write(6,*) 'visc',visc_LB
        write(6,*) 'fx',fx
        write(6,*) 'fy',fy
        write(6,*) 'fz',fz
        write(6,*) 'cssq',cssq
        write(6,*) '*******************INPUT data*****************'
        write(6,*) 'nx',nx
        write(6,*) 'ny',ny
        write(6,*) 'ny',nz
        write(6,*) 'lpbc',lpbc
        write(6,*) 'lprint',lprint
        write(6,*) 'lvtk',lvtk
        write(6,*) 'lasync',lasync
        write(6,*) 'nsteps',nsteps
        write(6,*) 'stamp',stamp
        write(6,*) 'dumpYN', dumpYN
        write(6,*) 'max fx',huge(fx)
        write(6,*) 'max fx',huge(fy)
        write(6,*) 'max fx',huge(fz)
        write(6,*) '*******************************************' !,normx,normy,normz,arrx,arry,arrz)
    !$acc data copy(f,g,isfluid,p0,p1,p2,p3,p0g,p1g,&
             !$acc& pxx,pyy,pzz,pxy,pxz,pyz,rho,u,v,w,normx,normy,normz,phi) 
    !$if _OPENACC        
        call printDeviceProperties(ngpus,devNum,devType,6)
    !$endif
    iframe=0
    iframe2D=0
    write(6,'(a,i8,a,i8,3f16.4)')'start step : ',0,' frame ',iframe
    
    if(lprint)then  
        call init_output(nx,ny,nz,1,lvtk)
        call string_char(head1,nheadervtk(1),headervtk(1))
        call string_char(head2,nheadervtk(2),headervtk(2))
    endif
    
    if(lprint)then
      !$acc kernels present(rhoprint,velprint,rho,u,v,w) 
      !$acc loop independent collapse(3)  private(i,j,k)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rhoprint(i,j,k)=real(rho(i,j,k),kind=4)
            velprint(1,i,j,k)=real(u(i,j,k),kind=4)
            velprint(2,i,j,k)=real(v(i,j,k),kind=4)
            velprint(3,i,j,k)=real(w(i,j,k),kind=4)
          enddo
        enddo
      enddo
      !$acc end kernels 
      !$acc update host(rhoprint,velprint)
      if(lvtk)then
        call print_vtk_sync(iframe)
      else
        call print_raw_sync(iframe)
      endif
    endif
      
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !***********************************moments collision bbck + forcing************************ 
          !$acc kernels
          !$acc loop collapse(3)
          do k=1,nz
              do j=1,ny
                  do i=1,nx
                      if(isfluid(i,j,k).eq.1)then
                          
                          pxx(i,j,k)=0.0_db
                          pyy(i,j,k)=0.0_db
                          pzz(i,j,k)=0.0_db
                          pxy(i,j,k)=0.0_db
                          pxz(i,j,k)=0.0_db
                          pyz(i,j,k)=0.0_db
                          !
						  phi(i,j,k) = g(i,j,k,0)+g(i,j,k,1)+g(i,j,k,2)+g(i,j,k,3)+g(i,j,k,4)+g(i,j,k,5)+g(i,j,k,6)
                          rho(i,j,k) = f(i,j,k,0)+f(i,j,k,1)+f(i,j,k,2)+f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,5) &
                              +f(i,j,k,6)+f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11) &
                              +f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17) &
                              +f(i,j,k,18)+f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,24) &
                              +f(i,j,k,25) +f(i,j,k,26)
                              

                          u(i,j,k) = ((f(i,j,k,1)+f(i,j,k,7)+f(i,j,k,9)+f(i,j,k,15)+f(i,j,k,18)++f(i,j,k,19)+f(i,j,k,21)+f(i,j,k,24)+f(i,j,k,25)) &
                              -(f(i,j,k,2)+f(i,j,k,8)+f(i,j,k,10)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,20)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,26)) )/rho(i,j,k)
                          
                          v(i,j,k) = ((f(i,j,k,3)+f(i,j,k,7)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,13)+f(i,j,k,19)+f(i,j,k,22)+f(i,j,k,24)+f(i,j,k,26)) &
                              -(f(i,j,k,4)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,12)+f(i,j,k,14)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,23)+f(i,j,k,25)) )/rho(i,j,k)

                          w(i,j,k) = ((f(i,j,k,5)+f(i,j,k,11)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,17)+f(i,j,k,19)+f(i,j,k,21)+f(i,j,k,23)+f(i,j,k,26)) &
                              -(f(i,j,k,6)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,16)+f(i,j,k,18)+f(i,j,k,20)+f(i,j,k,22)+f(i,j,k,24)+f(i,j,k,25)) )/rho(i,j,k)                        
                          
                          !1-2
                          feq=(rho(i,j,k)*(2 + 6*u(i,j,k) + 6*u(i,j,k)**2 - 3*v(i,j,k)**2 - 9*u(i,j,k)*v(i,j,k)**2 - 3*(1 + 3*u(i,j,k))*w(i,j,k)**2))/27.
                          fneq1=f(i,j,k,1)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          feq=(rho(i,j,k)*(2 - 3*v(i,j,k)**2 - 3*w(i,j,k)**2 + 3*u(i,j,k)*(-2 + 2*u(i,j,k) + 3*v(i,j,k)**2 + 3*w(i,j,k)**2)))/27.
                          fneq1=f(i,j,k,2)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          !3-4
                          feq=(rho(i,j,k)*(2 - 3*u(i,j,k)**2*(1 + 3*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(2 + 2*v(i,j,k) - 3*w(i,j,k)**2)))/27.
                          fneq1=f(i,j,k,3)-feq
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(-2 + 2*v(i,j,k) + 3*w(i,j,k)**2)))/27.
                          fneq1=f(i,j,k,4)-feq
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          !5-6
                          feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k))))/27.
                          fneq1=f(i,j,k,5)-feq
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          feq=(rho(i,j,k)*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                          fneq1=f(i,j,k,6)-feq
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          !7-8
                          feq=(2*rho(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) + u(i,j,k)*(3 + 9*v(i,j,k)*(1 + v(i,j,k)))) - 3*rho(i,j,k)*(1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=f(i,j,k,7)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k) + fneq1
                          feq=(2*rho(i,j,k)*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(-3 - 9*(-1 + v(i,j,k))*v(i,j,k))) + 3*rho(i,j,k)*(-1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=f(i,j,k,8)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)+fneq1
                          !10-9
                          feq=(2*rho(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)))) + 3*rho(i,j,k)*(-1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=f(i,j,k,10)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)-fneq1
                          feq=(2*rho(i,j,k)*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(3 + 9*(-1 + v(i,j,k))*v(i,j,k))) - 3*rho(i,j,k)*(1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=f(i,j,k,9)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)-fneq1
                          !11-12
                          feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=f(i,j,k,11)-feq
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)+fneq1
                          feq=(rho(i,j,k)*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                          fneq1=f(i,j,k,12)-feq
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)+fneq1
                          !13-14
                          feq=(rho(i,j,k)*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                          fneq1=f(i,j,k,13) - feq
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)-fneq1
                          feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=f(i,j,k,14) - feq
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)-fneq1
                          !15-16
                          feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=f(i,j,k,15)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)+fneq1
                          feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                          fneq1=f(i,j,k,16)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)+fneq1
                          !17-18
                          feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=f(i,j,k,17)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)-fneq1
                          feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                          fneq1=f(i,j,k,18)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)-fneq1
                          !19-20
                          feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                          fneq1=f(i,j,k,19)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)+fneq1

                          feq=(rho(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                          fneq1=f(i,j,k,20)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)+fneq1

                          !21-22
                          feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                          fneq1=f(i,j,k,21)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)-fneq1
                          pxz(i,j,k)=pxz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)-fneq1

                          feq=(rho(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                          fneq1=f(i,j,k,22)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)-fneq1
                          pxz(i,j,k)=pxz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)-fneq1
                          !23-24
                          feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                          fneq1=f(i,j,k,23)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)-fneq1
                          pyz(i,j,k)=pyz(i,j,k)-fneq1

                          feq=(rho(i,j,k)*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                          fneq1=f(i,j,k,24)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)-fneq1
                          pyz(i,j,k)=pyz(i,j,k)-fneq1
                          !25-26
                          feq=(rho(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                          fneq1=f(i,j,k,25)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)-fneq1
                          pxz(i,j,k)=pxz(i,j,k)-fneq1
                          pyz(i,j,k)=pyz(i,j,k)+fneq1

                          feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                          fneq1=f(i,j,k,26)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)-fneq1
                          pxz(i,j,k)=pxz(i,j,k)-fneq1
                          pyz(i,j,k)=pyz(i,j,k)+fneq1

                          
                      endif
                  enddo
              enddo
          enddo
		  !$acc loop collapse(3) 
		  do k=1,nz
        do j=1,ny
          do i=1,nx
            if(isfluid(i,j,k).eq.1)then
              ! normals calculation
                  ! 0  1   2  3   4   5   6   7    8   9   10  11   12  13   14  15   16   17   18  19  20  21  22  23  24  25  26
              !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1,  1, -1,  1, -1, -1,  1,  1, -1/)
              !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0,  1  -1, -1,  1, -1,  1, -1,  1/)
              !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1,  1, -1,  1, -1,  1, -1, -1,  1/)
              
                !if(phi(i,j,k).gt.0.01 .and. phi(i,j,k).lt.0.96)then
                  grad_fix=3.0*(p1*(phi(i+1,j,k)-phi(i-1,j,k)) + &
                  p2*(phi(i+1,j+1,k)+phi(i+1,j-1,k)+phi(i+1,j,k+1)+phi(i+1,j,k-1)-(phi(i-1,j+1,k)+phi(i-1,j-1,k)+phi(i-1,j,k+1)+phi(i-1,j,k-1)) ) + &
                  p3*(phi(i+1,j+1,k+1)+phi(i+1,j-1,k-1)+phi(i+1,j-1,k+1)+phi(i+1,j+1,k-1) - (phi(i-1,j+1,k+1)+phi(i-1,j-1,k-1)+phi(i-1,j-1,k+1)+phi(i-1,j+1,k-1))) )
                  
                  grad_fiy=3.0*(p1*(phi(i,j+1,k)-phi(i,j-1,k)) + &
                  p2*(phi(i+1,j+1,k)+phi(i-1,j+1,k)+phi(i,j+1,k+1)+phi(i,j+1,k-1)-(phi(i+1,j-1,k)+phi(i-1,j-1,k)+phi(i,j-1,k+1)+phi(i,j-1,k-1)) ) + &
                  p3*(phi(i+1,j+1,k+1)+phi(i-1,j+1,k-1)+phi(i+1,j+1,k-1)+phi(i-1,j+1,k+1)-(phi(i+1,j-1,k+1)+phi(i-1,j-1,k-1)+phi(i+1,j-1,k-1)+phi(i-1,j-1,k+1))))
                  
                  grad_fiz=3.0*(p1*(phi(i,j,k+1)-phi(i,j,k-1)) + &
                  p2*(phi(i+1,j,k+1)+phi(i-1,j,k+1)+phi(i,j+1,k+1)+phi(i,j-1,k+1)-(phi(i+1,j,k-1)+phi(i-1,j,k-1)+phi(i,j+1,k-1)+phi(i,j-1,k-1)) ) + &
                  p3*(phi(i+1,j+1,k+1)+phi(i-1,j-1,k+1)+phi(i+1,j-1,k+1)+phi(i-1,j+1,k+1)-(phi(i+1,j+1,k-1)+phi(i-1,j-1,k-1)+phi(i+1,j-1,k-1)+phi(i-1,j+1,k-1))))

                  mod_grad= sqrt(grad_fix**2 + grad_fiy**2 + grad_fiz**2)

                  normx(i,j,k)=grad_fix/(mod_grad+10.0**(-9))
                  normy(i,j,k)=grad_fiy/(mod_grad+10.0**(-9))
                  normz(i,j,k)=grad_fiz/(mod_grad+10.0**(-9))
                !endif

            endif				
          enddo
        enddo
		  enddo
          !$acc end kernels
        !***********************************Print on files 3D************************
          if(mod(step,stamp).eq.0)write(6,'(a,i8)')'step : ',step
            if(lprint)then
              if(mod(step,stamp).eq.0)then
                iframe=iframe+1
                !$acc kernels present(rhoprint,velprint,rho,u,v,w) 
                !$acc loop independent collapse(3)  private(i,j,k)
                do k=1,nz
                  do j=1,ny
                    do i=1,nx
                        rhoprint(i,j,k)=real(rho(i,j,k),kind=4)
                        velprint(1,i,j,k)=real(u(i,j,k),kind=4)
                        velprint(2,i,j,k)=real(v(i,j,k),kind=4)
                        velprint(3,i,j,k)=real(w(i,j,k),kind=4)
                    enddo
                  enddo
                enddo
              !$acc end kernels 
              !$acc update host(rhoprint,velprint) 
              if(lvtk)then
                call print_vtk_sync(iframe)
              else
                call print_raw_sync(iframe)
              endif
            endif
          endif
        !***********************************Print on files 2D************************
          if(mod(step,stamp2D).eq.0)write(6,'(a,i8)')'step : ',step
            if(lprint)then
              if(mod(step,stamp2D).eq.0)then
                iframe2D=iframe2D+1
                !$acc update host(rho,phi,u,v,w) 
                call print_raw_slice_2c_sync(iframe2D)
            endif
          endif
        !***********************************dump f************************
          if(mod(step,dumpstep).eq.0) then
                write(6,'(a,i8)')'dump step at : ',step
                !$acc update host(f) 
                call dump_distros_1c_3d
          endif
        !***********************************collision + no slip + forcing: fused implementation*********
          !$acc kernels
          !$acc loop collapse(3) private(feq,feqs,curvature,ffx_s,ffy_s,ffz_s)
          do k=1,nz
              do j=1,ny
                  do i=1,nx
                      if(isfluid(i,j,k).eq.1)then
                          ! surface tension caluclation
                          curvature=0.0
                          indicator=(phi(i,j,k)*(1-phi(i,j,k))*(normx(i,j,k)**2+normy(i,j,k)**2+normz(i,j,k)**2))
                          if(indicator .gt. 0.2)then
                              curvature=-(3.0*(p1*(normx(i+1,j,k)-normx(i-1,j,k)) + &
                              p2*(normx(i+1,j+1,k)+normx(i+1,j-1,k)+normx(i+1,j,k+1)+normx(i+1,j,k-1)-(normx(i-1,j+1,k)+normx(i-1,j-1,k)+normx(i-1,j,k+1)+normx(i-1,j,k-1)) ) + &
                              p3*(normx(i+1,j+1,k+1)+normx(i+1,j-1,k-1)+normx(i+1,j-1,k+1)+normx(i+1,j+1,k-1) - (normx(i-1,j+1,k+1)+normx(i-1,j-1,k-1)+normx(i-1,j-1,k+1)+normx(i-1,j+1,k-1))) ) + &
                                            
											   3.0*(p1*(normy(i,j+1,k)-normy(i,j-1,k)) + &
                              p2*(normy(i+1,j+1,k)+normy(i-1,j+1,k)+normy(i,j+1,k+1)+normy(i,j+1,k-1)-(normy(i+1,j-1,k)+normy(i-1,j-1,k)+normy(i,j-1,k+1)+normy(i,j-1,k-1)) ) + &
                              p3*(normy(i+1,j+1,k+1)+normy(i-1,j+1,k-1)+normy(i+1,j+1,k-1)+normy(i-1,j+1,k+1)-(normy(i+1,j-1,k+1)+normy(i-1,j-1,k-1)+normy(i+1,j-1,k-1)+normy(i-1,j-1,k+1)))) + &
                                            
											   3.0*(p1*(normz(i,j,k+1)-normz(i,j,k-1)) + &
                              p2*(normz(i+1,j,k+1)+normz(i-1,j,k+1)+normz(i,j+1,k+1)+normz(i,j-1,k+1)-(normz(i+1,j,k-1)+normz(i-1,j,k-1)+normz(i,j+1,k-1)+normz(i,j-1,k-1)) ) + &
                              p3*(normz(i+1,j+1,k+1)+normz(i-1,j-1,k+1)+normz(i+1,j-1,k+1)+normz(i-1,j+1,k+1)-(normz(i+1,j+1,k-1)+normz(i-1,j-1,k-1)+normz(i+1,j-1,k-1)+normz(i-1,j+1,k-1)))))
             
                          endif
                          ffx_s=sigma*curvature*normx(i,j,k)
                          ffy_s=sigma*curvature*normy(i,j,k)
                          ffz_s=sigma*curvature*normz(i,j,k)
                          
						              us=u(i,j,k)+ffx_s
                          vs=v(i,j,k)+ffy_s
                          ws=w(i,j,k)+ffz_s
                          !0
                          !feq=(-4*rho(i,j,k)*(-2 + 3*u(i,j,k)**2 + 3*v(i,j,k)**2 + 3*w(i,j,k)**2))/27.
						              feqs=(-4*rho(i,j,k)*(-2 + 3*us**2 + 3*vs**2 + 3*ws**2))/27.
                          fneq1=(-3*(pxx(i,j,k) + pyy(i,j,k) + pzz(i,j,k)))/2.
                          f(i,j,k,0)=feqs + (1-omega)*fneq1*p0 
                          
                          !1
                          
                          !feq=(rho(i,j,k)*(2 + 6*u(i,j,k) + 6*u(i,j,k)**2 - 3*v(i,j,k)**2 - 9*u(i,j,k)*v(i,j,k)**2 - 3*(1 + 3*u(i,j,k))*w(i,j,k)**2))/27.
                          feqs=(rho(i,j,k)*(2 + 6*us + 6*us**2 - 3*vs**2 - 9*us*vs**2 - 3*(1 + 3*us)*ws**2))/27.
                          fneq1=(3*(2*pxx(i,j,k) - pzz(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - pyy(i,j,k)*(1 + 3*u(i,j,k)) - 6*pxy(i,j,k)*v(i,j,k) - 6*pxz(i,j,k)*w(i,j,k)))/2.
                          f(i+1,j,k,1)=feqs + (1-omega)*fneq1*p1 !+ (feqs-feq)
                          
                          !2
                          !feq=(rho(i,j,k)*(2 - 3*v(i,j,k)**2 - 3*w(i,j,k)**2 + 3*u(i,j,k)*(-2 + 2*u(i,j,k) + 3*v(i,j,k)**2 + 3*w(i,j,k)**2)))/27.
                          feqs=(rho(i,j,k)*(2 - 3*vs**2 - 3*ws**2 + 3*us*(-2 + 2*us + 3*vs**2 + 3*ws**2)))/27.
                          fneq1=(3*(2*pxx(i,j,k) - pzz(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + pyy(i,j,k)*(-1 + 3*u(i,j,k)) + 6*pxy(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*w(i,j,k)))/2.
                          f(i-1,j,k,2)=feqs + (1-omega)*fneq1*p1 !- (feqs-feq)
                          
                          !3
                          
                          !feq=(rho(i,j,k)*(2 - 3*u(i,j,k)**2*(1 + 3*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(2 + 2*v(i,j,k) - 3*w(i,j,k)**2)))/27.
                          feqs=(rho(i,j,k)*(2 - 3*us**2*(1 + 3*vs) - 3*ws**2 + 3*vs*(2 + 2*vs - 3*ws**2)))/27.
                          fneq1=(-3*(pxx(i,j,k) - 2*pyy(i,j,k) + pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) + 3*pxx(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i,j+1,k,3)=feqs + (1-omega)*fneq1*p1 !+ (feqs-feq)
                          
                          !4
                          !feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(-2 + 2*v(i,j,k) + 3*w(i,j,k)**2)))/27.
                          feqs=(rho(i,j,k)*(2 + us**2*(-3 + 9*vs) - 3*ws**2 + 3*vs*(-2 + 2*vs + 3*ws**2)))/27.
                          fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(-1 + 3*v(i,j,k)) + 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i,j-1,k,4)=feqs + (1-omega)*fneq1*p1 !- (feqs-feq)
                          
                          !7
                          
                          !feq=(2*rho(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) + u(i,j,k)*(3 + 9*v(i,j,k)*(1 + v(i,j,k)))) - 3*rho(i,j,k)*(1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                          feqs=(2*rho(i,j,k)*(1 + 3*vs*(1 + vs) + us**2*(3 + 9*vs) + us*(3 + 9*vs*(1 + vs))) - 3*rho(i,j,k)*(1 + 3*us + 3*vs)*ws**2)/108.
                          fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) + 6*pyy(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*(1 + 2*u(i,j,k) + 2*v(i,j,k)) + pxx(i,j,k)*(2 + 6*v(i,j,k)) - 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i+1,j+1,k,7)=feqs + (1-omega)*fneq1*p2 !+ (feqs-feq)
                          
                          !8
                          !feq=(2*rho(i,j,k)*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(-3 - 9*(-1 + v(i,j,k))*v(i,j,k))) + 3*rho(i,j,k)*(-1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                          feqs=(2*rho(i,j,k)*(1 + us**2*(3 - 9*vs) + 3*(-1 + vs)*vs + us*(-3 - 9*(-1 + vs)*vs)) + 3*rho(i,j,k)*(-1 + 3*us + 3*vs)*ws**2)/108.
                          fneq1=(-3*(-2*pyy(i,j,k) + pzz(i,j,k) + 6*pyy(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*(-1 + 2*u(i,j,k) + 2*v(i,j,k)) + pxx(i,j,k)*(-2 + 6*v(i,j,k)) - 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i-1,j-1,k,8)=feqs + (1-omega)*fneq1*p2 !- (feqs-feq)
                          
                          !10
                          
                          !feq=(2*rho(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)))) + 3*rho(i,j,k)*(-1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                          feqs=(2*rho(i,j,k)*(1 + 3*vs*(1 + vs) + us**2*(3 + 9*vs) - 3*us*(1 + 3*vs*(1 + vs))) + 3*rho(i,j,k)*(-1 + 3*us - 3*vs)*ws**2)/108.
                          fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) - 6*pyy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*(-1 + 2*u(i,j,k) - 2*v(i,j,k)) - 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(2 + 6*v(i,j,k)) + 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i-1,j+1,k,10)=feqs + (1-omega)*fneq1*p2 !+ (feqs-feq)
                          
                          !9
                          !feq=(2*rho(i,j,k)*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(3 + 9*(-1 + v(i,j,k))*v(i,j,k))) - 3*rho(i,j,k)*(1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                          feqs=(2*rho(i,j,k)*(1 + us**2*(3 - 9*vs) + 3*(-1 + vs)*vs + us*(3 + 9*(-1 + vs)*vs)) - 3*rho(i,j,k)*(1 + 3*us - 3*vs)*ws**2)/108.
                          fneq1=(-3*(-2*pyy(i,j,k) + pzz(i,j,k) - 6*pyy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*(1 + 2*u(i,j,k) - 2*v(i,j,k)) - 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(-2 + 6*v(i,j,k)) + 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i+1,j-1,k,9)=feqs + (1-omega)*fneq1*p2 !+ (feqs-feq)

                          !5
                          
                          !feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k))))/27.
                          feqs=(rho(i,j,k)*(2 + 6*ws*(1 + ws) - 3*us**2*(1 + 3*ws) - 3*vs**2*(1 + 3*ws)))/27.
                          fneq1=(-3*(pxx(i,j,k) + pyy(i,j,k) - 2*pzz(i,j,k) + 6*pxz(i,j,k)*u(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pxx(i,j,k)*w(i,j,k) + 3*pyy(i,j,k)*w(i,j,k)))/2.
                          f(i,j,k+1,5)=feqs + (1-omega)*fneq1*p1 !+ (feqs-feq)
                          
                          !6
                          !feq=(rho(i,j,k)*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                          feqs=(rho(i,j,k)*(2 + 6*(-1 + ws)*ws + us**2*(-3 + 9*ws) + vs**2*(-3 + 9*ws)))/27.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 6*pxz(i,j,k)*u(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(-1 + 3*w(i,j,k))))/2.
                          f(i,j,k-1,6)=feqs + (1-omega)*fneq1*p1 !- (feqs-feq)

                          !15
                          
                          !feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          feqs=(rho(i,j,k)*(2 + 6*ws*(1 + ws) + 6*us**2*(1 + 3*ws) - 3*vs**2*(1 + 3*ws) + 3*us*(2 - 3*vs**2 + 6*ws*(1 + ws))))/108.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) - 3*pyy(i,j,k)*u(i,j,k) + 6*pzz(i,j,k)*u(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) - 3*pyy(i,j,k)*w(i,j,k) + 6*pxz(i,j,k)*(1 + 2*u(i,j,k) + 2*w(i,j,k)) + pxx(i,j,k)*(2 + 6*w(i,j,k))))/2.
                          f(i+1,j,k+1,15)=feqs + (1-omega)*fneq1*p2 !+ (feqs-feq) 
                          
                          !16
                          !feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                          feqs=(rho(i,j,k)*(2 + us**2*(6 - 18*ws) + 6*(-1 + ws)*ws + vs**2*(-3 + 9*ws) + 3*us*(-2 + 3*vs**2 - 6*(-1 + ws)*ws)))/108.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(2 - 6*w(i,j,k)) + 3*pyy(i,j,k)*w(i,j,k) - 6*pxz(i,j,k)*(-1 + 2*u(i,j,k) + 2*w(i,j,k))))/2.
                          f(i-1,j,k-1,16)=feqs + (1-omega)*fneq1*p2 !- (feqs-feq)

                          !17
                          
                          !feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          feqs=(rho(i,j,k)*(2 + 6*ws*(1 + ws) + 6*us**2*(1 + 3*ws) - 3*vs**2*(1 + 3*ws) + 3*us*(-2 + 3*vs**2 - 6*ws*(1 + ws))))/108.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*(-1 + 2*u(i,j,k) - 2*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(2 + 6*w(i,j,k))))/2.
                          f(i-1,j,k+1,17)=feqs + (1-omega)*fneq1*p2 !+ (feqs-feq)
                          
                          !18
                          !feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                          feqs=(rho(i,j,k)*(2 + us**2*(6 - 18*ws) + 6*(-1 + ws)*ws + vs**2*(-3 + 9*ws) + 3*us*(2 - 3*vs**2 + 6*(-1 + ws)*ws)))/108.
                          fneq1=(-3*(pyy(i,j,k) - 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*(1 + 2*u(i,j,k) - 2*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(-2 + 6*w(i,j,k))))/2.
                          f(i+1,j,k-1,18)=feqs + (1-omega)*fneq1*p2 !+ (feqs-feq)

                          !11
                          
                          !feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          feqs=(rho(i,j,k)*(2 + 6*ws*(1 + ws) + 6*vs**2*(1 + 3*ws) - 3*us**2*(1 + 3*vs + 3*ws) + 2*vs*(3 + 9*ws*(1 + ws))))/108.
                          fneq1=(-3*pxx(i,j,k)*(1 + 3*v(i,j,k) + 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) - 3*pxy(i,j,k)*u(i,j,k) - 3*pxz(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*pyy(i,j,k)*w(i,j,k) + pyz(i,j,k)*(3 + 6*v(i,j,k) + 6*w(i,j,k)))
                          f(i,j+1,k+1,11)=feqs + (1-omega)*fneq1*p2 !+ (feqs-feq)
                          
                          !12
                          feq=(rho(i,j,k)*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                          feqs=(rho(i,j,k)*(2 + 2*vs**2*(3 - 9*ws) + 6*(-1 + ws)*ws + us**2*(-3 + 9*vs + 9*ws) + 2*vs*(-3 - 9*(-1 + ws)*ws)))/108.
                          fneq1=(3*pxx(i,j,k)*(-1 + 3*v(i,j,k) + 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) + 3*pxy(i,j,k)*u(i,j,k) + 3*pxz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + pyz(i,j,k)*(3 - 6*v(i,j,k) - 6*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k))
                          f(i,j-1,k-1,12)=feqs + (1-omega)*fneq1*p2 !- (feqs-feq)

                          !13
                          
                          !feq=(rho(i,j,k)*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                          feqs=(rho(i,j,k)*(2 - 6*ws + us**2*(-3 - 9*vs + 9*ws) + 6*(vs + vs**2*(1 - 3*ws) + 3*vs*(-1 + ws)*ws + ws**2)))/108.
                          fneq1=(-3*(pxx(i,j,k) - 2*pyy(i,j,k) + 6*pyz(i,j,k) - 2*pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) - 6*pxz(i,j,k)*u(i,j,k) + 3*pxx(i,j,k)*v(i,j,k) + 12*pyz(i,j,k)*v(i,j,k) - 6*pzz(i,j,k)*v(i,j,k) - 3*pxx(i,j,k)*w(i,j,k) + 6*pyy(i,j,k)*w(i,j,k) - 12*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i,j+1,k-1,13)=feqs + (1-omega)*fneq1*p2 !+ (feqs-feq)
                          
                          !14
                          !feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          feqs=(rho(i,j,k)*(2 + us**2*(-3 + 9*vs - 9*ws) + 6*ws*(1 + ws) + 6*vs**2*(1 + 3*ws) - 6*vs*(1 + 3*ws*(1 + ws))))/108.
                          fneq1=(3*pxx(i,j,k)*(-1 + 3*v(i,j,k) - 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) + 3*pxy(i,j,k)*u(i,j,k) - 3*pxz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + pyz(i,j,k)*(-3 + 6*v(i,j,k) - 6*w(i,j,k)) + 3*pyy(i,j,k)*w(i,j,k))
                          f(i,j-1,k+1,14)=feqs + (1-omega)*fneq1*p2 !+ (feqs-feq)

                          !19
                          !feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                          feqs=(rho(i,j,k)*(1 + 3*ws*(1 + ws) + vs**2*(3 + 9*ws) + us**2*(3 + 9*vs + 9*ws) + vs*(3 + 9*ws*(1 + ws)) + 3*us*(1 + 3*ws + 3*(vs + vs**2 + 3*vs*ws + ws**2))))/216.
                          fneq1=3*(pxx(i,j,k) + (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) + 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j+1,k+1,19)=feqs + (1-omega)*fneq1*p3 !+ (feqs-feq)
                          
                          !20
                          !feq=(rho(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                          feqs=(rho(i,j,k)*(1 + vs**2*(3 - 9*ws) + us**2*(3 - 9*vs - 9*ws) + 3*(-1 + ws)*ws + vs*(-3 - 9*(-1 + ws)*ws) - 3*us*(1 + 3*vs**2 + 3*(-1 + ws)*ws + vs*(-3 + 9*ws))))/216.
                          fneq1=-3*((pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(-1 + 3*u(i,j,k)) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + pxz(i,j,k)*(-3 + 6*u(i,j,k)) + 3*(2*pxy(i,j,k) + 3*pxz(i,j,k) + 2*pyz(i,j,k) + pzz(i,j,k))*v(i,j,k) + 3*(3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k) + pxx(i,j,k)*(-1 + 3*v(i,j,k) + 3*w(i,j,k)))
                          f(i-1,j-1,k-1,20)=feqs + (1-omega)*fneq1*p3 !- (feqs-feq)

                          !21
                          !feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                          feqs=(rho(i,j,k)*(1 + 3*ws*(1 + ws) + vs**2*(3 + 9*ws) + us**2*(3 - 9*vs + 9*ws) - 3*vs*(1 + 3*ws*(1 + ws)) + 3*us*(1 + 3*vs**2 + 3*ws*(1 + ws) - 3*vs*(1 + 3*ws))))/216.
						              fneq1=3*(pxx(i,j,k) - 3*pxy(i,j,k)*(1 + 2*u(i,j,k)) + (pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) - 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j-1,k+1,21)=feqs + (1-omega)*fneq1*p3 !+ (feqs-feq)
                          
                          !22
                          !feq=(rho(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                          feqs=(rho(i,j,k)*(1 + vs**2*(3 - 9*ws) + us**2*(3 + 9*vs - 9*ws) + 3*(-1 + ws)*ws + vs*(3 + 9*(-1 + ws)*ws) - 3*us*(1 - 3*ws + 3*(vs + vs**2 - 3*vs*ws + ws**2))))/216.
                          fneq1=3*(pxx(i,j,k) + 3*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k) - 3*(2*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) - 3*(pxx(i,j,k) - 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i-1,j+1,k-1,22)=feqs+ (1-omega)*fneq1*p3 !- (feqs-feq)

                          !23
                          !feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                          feqs=(rho(i,j,k)*(1 + 3*ws*(1 + ws) + vs**2*(3 + 9*ws) + us**2*(3 - 9*vs + 9*ws) - 3*vs*(1 + 3*ws*(1 + ws)) - 3*us*(1 + 3*vs**2 + 3*ws*(1 + ws) - 3*vs*(1 + 3*ws))))/216.
                          fneq1=3*(pxx(i,j,k) + 3*pxy(i,j,k) - 3*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k) - 3*(2*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k) - 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) + 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i-1,j-1,k+1,23)=feqs + (1-omega)*fneq1*p3 !+ (feqs-feq)

                          !24
                          !feq=(rho(i,j,k)*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                          feqs=(rho(i,j,k)*(1 + 3*vs - 3*ws+ 3*(us + us**2 + 3*us*vs + 3*us**2*vs + vs**2 + 3*us*vs**2 - 3*(us + us**2 + vs + 3*us*vs + vs**2)*ws + (1 + 3*us + 3*vs)*ws**2)))/216.
                          fneq1=3*(pxx(i,j,k) - 3*pxz(i,j,k)*(1 + 2*u(i,j,k)) + (pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) - 3*(pxx(i,j,k) + 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j+1,k-1,24)=feqs + (1-omega)*fneq1*p3 !+ (feqs-feq)

                          !25
                          !feq=(rho(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                          feqs=(rho(i,j,k)*(1 + vs**2*(3 - 9*ws) + us**2*(3 - 9*vs - 9*ws) + 3*(-1 + ws)*ws + vs*(-3 - 9*(-1 + ws)*ws) + 3*us*(1 + 3*vs**2 + 3*(-1 + ws)*ws + vs*(-3 + 9*ws))))/216.
                          fneq1=-3*(-pxx(i,j,k) - (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j-1,k-1,25)=feqs + (1-omega)*fneq1*p3 !+ (feqs-feq)

                          !26
                          !feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                          feqs=(rho(i,j,k)*(1 + 3*ws*(1 + ws) + vs**2*(3 + 9*ws) + us**2*(3 + 9*vs + 9*ws) + vs*(3 + 9*ws*(1 + ws)) - 3*us*(1 + 3*ws + 3*(vs + vs**2 + 3*vs*ws + ws**2))))/216.
                          fneq1=3*(pxx(i,j,k) - (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(-1 + 3*u(i,j,k)) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + pxz(i,j,k)*(-3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                          f(i-1,j+1,k+1,26)=feqs + (1-omega)*fneq1*p3 !+ (feqs-feq)
							
						              !g0
                          feq=p0g*phi(i,j,k)
                          g(i,j,k,0)=feq
                          !g1
                          feqs=p1g*sharp_c*phi(i,j,k)*(1-phi(i,j,k))*(normx(i,j,k))
                          feq=p1g*phi(i,j,k)*(1 + 3*u(i,j,k))
                          g(i+1,j,k,1)=feq+feqs
                          !g2
                          feq=p1g*phi(i,j,k)*(1 - 3*u(i,j,k))
                          g(i-1,j,k,2)=feq-feqs
                          !g3
                          feqs=p1g*sharp_c*phi(i,j,k)*(1-phi(i,j,k))*(normy(i,j,k))
                          feq=p1g*phi(i,j,k)*(1 + 3*v(i,j,k))
                          g(i,j+1,k,3)=feq+feqs
                          !g4
                          feq=p1g*phi(i,j,k)*(1 - 3*v(i,j,k))
                          g(i,j-1,k,4)=feq-feqs
                          !g5
                          feqs=p1g*sharp_c*phi(i,j,k)*(1-phi(i,j,k))*(normz(i,j,k))
                          feq=p1g*phi(i,j,k)*(1 + 3*w(i,j,k))
                          g(i,j,k+1,5)=feq+feqs
                          !g6
                          feq=p1g*phi(i,j,k)*(1 - 3*w(i,j,k))
                          g(i,j,k-1,6)=feq-feqs
                      endif
                  enddo
              enddo
          enddo
          !$acc end kernels
        !***********************************boundary conditions ********************************!
        call bcs_all_bback_2c
      
    enddo 
    !$acc end data
    call cpu_time(ts2)
    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
    write(6,*) 'glups: ',  real(nx)*real(ny)*real(nz)*real(nsteps)/1.0e9/(ts2-ts1)
    
    call get_memory_gpu(mymemory,totmemory)
    call print_memory_registration_gpu(6,'DEVICE memory occupied at the end', &
     'total DEVICE memory',mymemory,totmemory)
    
end program
