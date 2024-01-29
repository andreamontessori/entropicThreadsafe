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
    real(kind=db) :: one_ov_nu1, one_ov_nu2, st_coeff, sigma, beta, rtot, psi_x, psi_y,psi_z
    real(kind=db) :: mod_psi, mod_psi_sq, rprod, nu_avg, om1, om2, fpc, b0, b1, b2, b3
    !$if _OPENACC
    integer :: devNum
    integer(acc_device_kind) :: devType
    devType = acc_get_device_type()
    devNum=acc_get_device_num(devType)
    !$endif

    nlinks=26 !pari!
    tau=1.0_db
    cssq=1.0_db/3.0_db
    visc_LB=cssq*(tau-0.5_db)
    one_ov_nu1=1.0_db/visc_LB !fluid 1
    tau=1.0_db
    visc_LB=cssq*(tau-0.5_db)
    one_ov_nu2=1.0_db/visc_LB !fluid 2
    dumpYN=0


#ifdef _OPENACC
    ngpus=acc_get_num_devices(acc_device_nvidia)
#else
    ngpus=0
#endif

    !*******************************user parameters and allocations**************************
        nx=128
        ny=128
        nz=128
        nsteps=5000
        stamp=2000000000
        stamp2D=5000
        dumpstep=1000000000
        fx=0.0_db*10.0**(-7)
        fy=0.0_db*10.0**(-5)
        fz=0.0_db*10.0**(-5)
        lprint=.true.
        lvtk=.true.
        lasync=.false.
        lpbc=.true.
        
        allocate(f(0:nx+1,0:ny+1,0:nz+1,0:nlinks),g(0:nx+1,0:ny+1,0:nz+1,0:nlinks))
        allocate(rhoA(1:nx,1:ny,1:nz),rhoB(1:nx,1:ny,1:nz),u(1:nx,1:ny,1:nz),v(1:nx,1:ny,1:nz),w(1:nx,1:ny,1:nz))
        allocate(pxx(1:nx,1:ny,1:nz),pxy(1:nx,1:ny,1:nz),pxz(1:nx,1:ny,1:nz),pyy(1:nx,1:ny,1:nz))
        allocate(pyz(1:nx,1:ny,1:nz),pzz(1:nx,1:ny,1:nz))
        allocate(isfluid(1:nx,1:ny,1:nz),psi(1:nx,1:ny,1:nz))
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
        b0=-10.0/27.0
        b1=2.0/27.0
        b2=1.0/54.0
        b3=1.0/216.0
        p1dcssq=p1/cssq
        p2dcssq=p2/cssq
        p3dcssq=p3/cssq
        omega=1.0_db/tau

    !*********************************chromodynamics vars*****************************
        beta=0.95_db
        sigma=0.001_db
        st_coeff=(9.0_db/4.0_db)*sigma*omega
    !************************************geometry*************************************
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
        rhoA=0.0_db  
        rhoB=1.0_db
        psi=-1.0_db

        do k=1,nz
          do j=1,ny
            do i=1,nx
              if(sqrt((float(i)-nx/2.)**2 + (float(j)-ny/2.)**2 + (float(k)-nz/2.)**2)<20)then
                  rhoA(i,j,k)=1
                  rhoB(i,j,k)=0  
              endif
            enddo
          enddo
        enddo
        !do ll=0,nlinks
        if(dumpYN.eq.0)then
            do k=1,nz
                  do j=1,ny
                      do i=1,nx
                          if(isfluid(i,j,k).eq.1)then
                              !rhoA
                                !0
                                feq=(-4*rhoA(i,j,k)*(-2 + 3*u(i,j,k)**2 + 3*v(i,j,k)**2 + 3*w(i,j,k)**2))/27.
                                
                                f(i,j,k,0)=feq 
                                
                                !1
                                
                                feq=(rhoA(i,j,k)*(2 + 6*u(i,j,k) + 6*u(i,j,k)**2 - 3*v(i,j,k)**2 - 9*u(i,j,k)*v(i,j,k)**2 - 3*(1 + 3*u(i,j,k))*w(i,j,k)**2))/27.
                                
                                f(i,j,k,1)=feq 
                                
                                !2
                                feq=(rhoA(i,j,k)*(2 - 3*v(i,j,k)**2 - 3*w(i,j,k)**2 + 3*u(i,j,k)*(-2 + 2*u(i,j,k) + 3*v(i,j,k)**2 + 3*w(i,j,k)**2)))/27.
                                
                                f(i,j,k,2)=feq 
                                
                                !3
                                
                                feq=(rhoA(i,j,k)*(2 - 3*u(i,j,k)**2*(1 + 3*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(2 + 2*v(i,j,k) - 3*w(i,j,k)**2)))/27.
                                
                                f(i,j,k,3)=feq
                                
                                !4
                                feq=(rhoA(i,j,k)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(-2 + 2*v(i,j,k) + 3*w(i,j,k)**2)))/27.
                                
                                f(i,j,k,4)=feq
                                
                                !7
                                
                                feq=(2*rhoA(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) + u(i,j,k)*(3 + 9*v(i,j,k)*(1 + v(i,j,k)))) - 3*rhoA(i,j,k)*(1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                                
                                f(i,j,k,7)=feq 
                                
                                !8
                                feq=(2*rhoA(i,j,k)*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(-3 - 9*(-1 + v(i,j,k))*v(i,j,k))) + 3*rhoA(i,j,k)*(-1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                                
                                f(i,j,k,8)=feq 
                                
                                !10
                                
                                feq=(2*rhoA(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)))) + 3*rhoA(i,j,k)*(-1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                                
                                f(i,j,k,10)=feq
                                
                                !9
                                feq=(2*rhoA(i,j,k)*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(3 + 9*(-1 + v(i,j,k))*v(i,j,k))) - 3*rhoA(i,j,k)*(1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                                
                                f(i,j,k,9)=feq

                                !5
                                
                                feq=(rhoA(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k))))/27.
                                
                                f(i,j,k,5)=feq
                                
                                !6
                                feq=(rhoA(i,j,k)*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                                
                                f(i,j,k,6)=feq

                                !15
                                
                                feq=(rhoA(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                                
                                f(i,j,k,15)=feq
                                
                                !16
                                feq=(rhoA(i,j,k)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                                
                                f(i,j,k,16)=feq

                                !17
                                
                                feq=(rhoA(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                                
                                f(i,j,k,17)=feq
                                
                                !18
                                feq=(rhoA(i,j,k)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                                
                                f(i,j,k,18)=feq

                                !11
                                
                                feq=(rhoA(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                                
                                f(i,j,k,11)=feq
                                
                                !12
                                feq=(rhoA(i,j,k)*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                                
                                f(i,j,k,12)=feq

                                !13
                                
                                feq=(rhoA(i,j,k)*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                                
                                f(i,j,k,13)=feq
                                
                                !14
                                feq=(rhoA(i,j,k)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                                
                                f(i,j,k,14)=feq

                                !19
                                feq=(rhoA(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                                
                                f(i,j,k,19)=feq 
                                
                                !20
                                feq=(rhoA(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                                
                                f(i,j,k,20)=feq

                                !21
                                feq=(rhoA(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                                
                                f(i,j,k,21)=feq
                                
                                !22
                                feq=(rhoA(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                                
                                f(i,j,k,22)=feq

                                !23
                                feq=(rhoA(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                                
                                f(i,j,k,23)=feq

                                !24
                                feq=(rhoA(i,j,k)*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                                
                                f(i,j,k,24)=feq

                                !25
                                feq=(rhoA(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                                
                                f(i,j,k,25)=feq 

                                !26
                                feq=(rhoA(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                                
                                f(i,j,k,26)=feq 
                              !rhoB
                                !0
                              
                                feq=(-4*rhoB(i,j,k)*(-2 + 3*u(i,j,k)**2 + 3*v(i,j,k)**2 + 3*w(i,j,k)**2))/27.
                                
                                g(i,j,k,0)=feq 
                                
                                !1
                                
                                feq=(rhoB(i,j,k)*(2 + 6*u(i,j,k) + 6*u(i,j,k)**2 - 3*v(i,j,k)**2 - 9*u(i,j,k)*v(i,j,k)**2 - 3*(1 + 3*u(i,j,k))*w(i,j,k)**2))/27.
                                
                                g(i,j,k,1)=feq 
                                
                                !2
                                feq=(rhoB(i,j,k)*(2 - 3*v(i,j,k)**2 - 3*w(i,j,k)**2 + 3*u(i,j,k)*(-2 + 2*u(i,j,k) + 3*v(i,j,k)**2 + 3*w(i,j,k)**2)))/27.
                                
                                g(i,j,k,2)=feq 
                                
                                !3
                                
                                feq=(rhoB(i,j,k)*(2 - 3*u(i,j,k)**2*(1 + 3*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(2 + 2*v(i,j,k) - 3*w(i,j,k)**2)))/27.
                                
                                g(i,j,k,3)=feq
                                
                                !4
                                feq=(rhoB(i,j,k)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(-2 + 2*v(i,j,k) + 3*w(i,j,k)**2)))/27.
                                
                                g(i,j,k,4)=feq
                                
                                !7
                                
                                feq=(2*rhoB(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) + u(i,j,k)*(3 + 9*v(i,j,k)*(1 + v(i,j,k)))) - 3*rhoB(i,j,k)*(1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                                
                                g(i,j,k,7)=feq 
                                
                                !8
                                feq=(2*rhoB(i,j,k)*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(-3 - 9*(-1 + v(i,j,k))*v(i,j,k))) + 3*rhoB(i,j,k)*(-1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                                
                                g(i,j,k,8)=feq 
                                
                                !10
                                
                                feq=(2*rhoB(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)))) + 3*rhoB(i,j,k)*(-1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                                
                                g(i,j,k,10)=feq
                                
                                !9
                                feq=(2*rhoB(i,j,k)*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(3 + 9*(-1 + v(i,j,k))*v(i,j,k))) - 3*rhoB(i,j,k)*(1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                                
                                g(i,j,k,9)=feq

                                !5
                                
                                feq=(rhoB(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k))))/27.
                                
                                g(i,j,k,5)=feq
                                
                                !6
                                feq=(rhoB(i,j,k)*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                                
                                g(i,j,k,6)=feq

                                !15
                                
                                feq=(rhoB(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                                
                                g(i,j,k,15)=feq
                                
                                !16
                                feq=(rhoB(i,j,k)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                                
                                g(i,j,k,16)=feq

                                !17
                                
                                feq=(rhoB(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                                
                                g(i,j,k,17)=feq
                                
                                !18
                                feq=(rhoB(i,j,k)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                                
                                g(i,j,k,18)=feq

                                !11
                                
                                feq=(rhoB(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                                
                                g(i,j,k,11)=feq
                                
                                !12
                                feq=(rhoB(i,j,k)*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                                
                                g(i,j,k,12)=feq

                                !13
                                
                                feq=(rhoB(i,j,k)*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                                
                                g(i,j,k,13)=feq
                                
                                !14
                                feq=(rhoB(i,j,k)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                                
                                g(i,j,k,14)=feq

                                !19
                                feq=(rhoB(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                                
                                g(i,j,k,19)=feq 
                                
                                !20
                                feq=(rhoB(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                                
                                g(i,j,k,20)=feq

                                !21
                                feq=(rhoB(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                                
                                g(i,j,k,21)=feq
                                
                                !22
                                feq=(rhoB(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                                
                                g(i,j,k,22)=feq

                                !23
                                feq=(rhoB(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                                
                                g(i,j,k,23)=feq

                                !24
                                feq=(rhoB(i,j,k)*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                                
                                g(i,j,k,24)=feq

                                !25
                                feq=(rhoB(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                                
                                g(i,j,k,25)=feq 

                                !26
                                feq=(rhoB(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                                
                                g(i,j,k,26)=feq
                          endif
                      enddo
                  enddo
              enddo
          elseif(dumpYN.eq.1)then
              call read_distros_2c_3d
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
        write(6,*) '*******************************************'
    !$acc data copy(f,g,isfluid,p0,p1,p2,p3,&
             !$acc& pxx,pyy,pzz,pxy,pxz,pyz,rhoA,rhoB,psi,u,v,w,rhoprint,velprint)
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
      !$acc kernels present(rhoprint,velprint,rhoA,u,v,w) 
      !$acc loop independent collapse(3)  private(i,j,k)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rhoprint(i,j,k)=real(rhoA(i,j,k),kind=4)
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
          !$acc loop collapse(3) private(fneq1,feq,uu,udotc,rtot)
          do k=1,nz
              do j=1,ny
                  do i=1,nx
                      if(isfluid(i,j,k).eq.1)then
                          !
                          pxx(i,j,k)=0.0_db
                          pyy(i,j,k)=0.0_db
                          pzz(i,j,k)=0.0_db
                          pxy(i,j,k)=0.0_db
                          pxz(i,j,k)=0.0_db
                          pyz(i,j,k)=0.0_db
                          !
                          rhoA(i,j,k) = f(i,j,k,0)+f(i,j,k,1)+f(i,j,k,2)+f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,5) &
                              +f(i,j,k,6)+f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11) &
                              +f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17) &
                              +f(i,j,k,18)+f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,24) &
                              +f(i,j,k,25) +f(i,j,k,26)

                          rhoB(i,j,k) = g(i,j,k,0)+g(i,j,k,1)+g(i,j,k,2)+g(i,j,k,3)+g(i,j,k,4)+g(i,j,k,5) &
                              +g(i,j,k,6)+g(i,j,k,7)+g(i,j,k,8)+g(i,j,k,9)+g(i,j,k,10)+g(i,j,k,11) &
                              +g(i,j,k,12)+g(i,j,k,13)+g(i,j,k,14)+g(i,j,k,15)+g(i,j,k,16)+g(i,j,k,17) &
                              +g(i,j,k,18)+g(i,j,k,19)+g(i,j,k,20)+g(i,j,k,21)+g(i,j,k,22)+g(i,j,k,23)+g(i,j,k,24) &
                              +g(i,j,k,25) +g(i,j,k,26)
                              
                          rtot=rhoA(i,j,k)+rhoB(i,j,k)

                          u(i,j,k) = ((f(i,j,k,1)+f(i,j,k,7)+f(i,j,k,9)+f(i,j,k,15)+f(i,j,k,18)+f(i,j,k,19)+f(i,j,k,21)+f(i,j,k,24)+f(i,j,k,25) &
                                      + g(i,j,k,1)+g(i,j,k,7)+g(i,j,k,9)+g(i,j,k,15)+g(i,j,k,18)+g(i,j,k,19)+g(i,j,k,21)+g(i,j,k,24)+g(i,j,k,25)) &
                              -(f(i,j,k,2)+f(i,j,k,8)+f(i,j,k,10)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,20)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,26)  &
                              + g(i,j,k,2)+g(i,j,k,8)+g(i,j,k,10)+g(i,j,k,16)+g(i,j,k,17)+g(i,j,k,20)+g(i,j,k,22)+g(i,j,k,23)+g(i,j,k,26)))/rtot
                          
                          v(i,j,k) = ((f(i,j,k,3)+f(i,j,k,7)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,13)+f(i,j,k,19)+f(i,j,k,22)+f(i,j,k,24)+f(i,j,k,26) &
                                    + g(i,j,k,3)+g(i,j,k,7)+g(i,j,k,10)+g(i,j,k,11)+g(i,j,k,13)+g(i,j,k,19)+g(i,j,k,22)+g(i,j,k,24)+g(i,j,k,26)) &
                              -(f(i,j,k,4)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,12)+f(i,j,k,14)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,23)+f(i,j,k,25) &
                              + g(i,j,k,4)+g(i,j,k,8)+g(i,j,k,9)+g(i,j,k,12)+g(i,j,k,14)+g(i,j,k,20)+g(i,j,k,21)+g(i,j,k,23)+g(i,j,k,25)))/rtot

                          w(i,j,k) = ((f(i,j,k,5)+f(i,j,k,11)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,17)+f(i,j,k,19)+f(i,j,k,21)+f(i,j,k,23)+f(i,j,k,26) +&
                                      g(i,j,k,5)+g(i,j,k,11)+g(i,j,k,14)+g(i,j,k,15)+g(i,j,k,17)+g(i,j,k,19)+g(i,j,k,21)+g(i,j,k,23)+g(i,j,k,26) ) &
                              -(f(i,j,k,6)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,16)+f(i,j,k,18)+f(i,j,k,20)+f(i,j,k,22)+f(i,j,k,24)+f(i,j,k,25) + &
                                g(i,j,k,6)+g(i,j,k,12)+g(i,j,k,13)+g(i,j,k,16)+g(i,j,k,18)+g(i,j,k,20)+g(i,j,k,22)+g(i,j,k,24)+g(i,j,k,25)))/rtot                        
                        
                          
                          
                          psi(i,j,k)=(rhoA(i,j,k)-rhoB(i,j,k))/rtot
                          !1-2
                          feq=(rtot*(2 + 6*u(i,j,k) + 6*u(i,j,k)**2 - 3*v(i,j,k)**2 - 9*u(i,j,k)*v(i,j,k)**2 - 3*(1 + 3*u(i,j,k))*w(i,j,k)**2))/27.
                          fneq1=f(i,j,k,1)+g(i,j,k,1)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          feq=(rtot*(2 - 3*v(i,j,k)**2 - 3*w(i,j,k)**2 + 3*u(i,j,k)*(-2 + 2*u(i,j,k) + 3*v(i,j,k)**2 + 3*w(i,j,k)**2)))/27.
                          fneq1=f(i,j,k,2)+g(i,j,k,2)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          !3-4
                          feq=(rtot*(2 - 3*u(i,j,k)**2*(1 + 3*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(2 + 2*v(i,j,k) - 3*w(i,j,k)**2)))/27.
                          fneq1=f(i,j,k,3)+g(i,j,k,3)-feq
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          feq=(rtot*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(-2 + 2*v(i,j,k) + 3*w(i,j,k)**2)))/27.
                          fneq1=f(i,j,k,4)+g(i,j,k,4)-feq
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          !5-6
                          feq=(rtot*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k))))/27.
                          fneq1=f(i,j,k,5)+g(i,j,k,5)-feq
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          feq=(rtot*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                          fneq1=f(i,j,k,6)+g(i,j,k,6)-feq
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          !7-8
                          feq=(2*rtot*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) + u(i,j,k)*(3 + 9*v(i,j,k)*(1 + v(i,j,k)))) - 3*rtot*(1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=f(i,j,k,7)+g(i,j,k,7)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k) + fneq1
                          feq=(2*rtot*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(-3 - 9*(-1 + v(i,j,k))*v(i,j,k))) + 3*rtot*(-1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=f(i,j,k,8)+g(i,j,k,8)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)+fneq1
                          !10-9
                          feq=(2*rtot*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)))) + 3*rtot*(-1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=f(i,j,k,10)+g(i,j,k,10)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)-fneq1
                          feq=(2*rtot*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(3 + 9*(-1 + v(i,j,k))*v(i,j,k))) - 3*rtot*(1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=f(i,j,k,9)+g(i,j,k,9)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)-fneq1
                          !11-12
                          feq=(rtot*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=f(i,j,k,11)+g(i,j,k,11)-feq
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)+fneq1
                          feq=(rtot*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                          fneq1=f(i,j,k,12)+g(i,j,k,12)-feq
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)+fneq1
                          !13-14
                          feq=(rtot*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                          fneq1=f(i,j,k,13)+g(i,j,k,13) - feq
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)-fneq1
                          feq=(rtot*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=f(i,j,k,14)+g(i,j,k,14) - feq
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)-fneq1
                          !15-16
                          feq=(rtot*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=f(i,j,k,15)+g(i,j,k,15)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)+fneq1
                          feq=(rtot*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                          fneq1=f(i,j,k,16)+g(i,j,k,16)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)+fneq1
                          !17-18
                          feq=(rtot*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=f(i,j,k,17)+g(i,j,k,17)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)-fneq1
                          feq=(rtot*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                          fneq1=f(i,j,k,18)+g(i,j,k,18)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)-fneq1
                          !19-20
                          feq=(rtot*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                          fneq1=f(i,j,k,19)+g(i,j,k,19)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)+fneq1

                          feq=(rtot*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                          fneq1=f(i,j,k,20)+g(i,j,k,20)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)+fneq1

                          !21-22
                          feq=(rtot*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                          fneq1=f(i,j,k,21)+g(i,j,k,21)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)-fneq1
                          pxz(i,j,k)=pxz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)-fneq1

                          feq=(rtot*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                          fneq1=f(i,j,k,22)+g(i,j,k,22)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)-fneq1
                          pxz(i,j,k)=pxz(i,j,k)+fneq1
                          pyz(i,j,k)=pyz(i,j,k)-fneq1
                          !23-24
                          feq=(rtot*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                          fneq1=f(i,j,k,23)+g(i,j,k,23)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)-fneq1
                          pyz(i,j,k)=pyz(i,j,k)-fneq1

                          feq=(rtot*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                          fneq1=f(i,j,k,24)+g(i,j,k,24)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)+fneq1
                          pxz(i,j,k)=pxz(i,j,k)-fneq1
                          pyz(i,j,k)=pyz(i,j,k)-fneq1
                          !25-26
                          feq=(rtot*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                          fneq1=f(i,j,k,25)+g(i,j,k,25)-feq
                          pxx(i,j,k)=pxx(i,j,k)+fneq1
                          pyy(i,j,k)=pyy(i,j,k)+fneq1
                          pzz(i,j,k)=pzz(i,j,k)+fneq1
                          pxy(i,j,k)=pxy(i,j,k)-fneq1
                          pxz(i,j,k)=pxz(i,j,k)-fneq1
                          pyz(i,j,k)=pyz(i,j,k)+fneq1

                          feq=(rtot*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                          fneq1=f(i,j,k,26)+g(i,j,k,26)-feq
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
          !$acc end kernels
        !***********************************Print on files 3D************************
          if(mod(step,stamp).eq.0)write(6,'(a,i8)')'step : ',step
            if(lprint)then
              if(mod(step,stamp).eq.0)then
                iframe=iframe+1
                !$acc kernels present(rhoprint,velprint,rhoA,u,v,w) 
                !$acc loop independent collapse(3)  private(i,j,k)
                do k=1,nz
                  do j=1,ny
                    do i=1,nx
                        rhoprint(i,j,k)=real(rhoA(i,j,k),kind=4)
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
                !$acc update host(rhoA,u,v,w) 
                call print_raw_slice_2c_sync(iframe2D)
            endif
          endif
        !***********************************dump f************************
          if(mod(step,dumpstep).eq.0) then
                write(6,'(a,i8)')'dump step at : ',step
                !$acc update host(f,g) 
                call dump_distros_2c_3d
          endif
        !***********************************collision + no slip + forcing: fused implementation*********
          !$acc kernels
          !$acc loop collapse(3) private(fpc,om2,om1,nu_avg,st_coeff,rtot,mod_psi_sq,psi_z,rprod,mod_psi,psi_y,psi_x)
          do k=1,nz
              do j=1,ny
                  do i=1,nx
                      if(isfluid(i,j,k).eq.1)then
                        !*******************************chromodynamics***************************!
                          psi_x=(1.0_db/cssq)*(p1*(psi(i+1,j,k)-psi(i-1,j,k)) + &
                                              p2*(psi(i+1,j+1,k)+psi(i+1,j-1,k)+psi(i+1,j,k+1)+psi(i+1,j,k-1)) &
                                            -p2*(psi(i-1,j+1,k)+psi(i-1,j-1,k)+psi(i-1,j,k+1)+psi(i-1,j,k-1)) &
                                            + p3*(psi(i+1,j+1,k+1)+psi(i+1,j+1,k-1)+psi(i+1,j-1,k+1)+psi(i+1,j-1,k-1)) &
                                            -p3*(psi(i-1,j+1,k+1)+psi(i-1,j+1,k-1)+psi(i-1,j-1,k+1)+psi(i-1,j-1,k-1))) 


                          psi_y=(1.0_db/cssq)*(p1*(psi(i,j+1,k)-psi(i,j-1,k)) &
                                              +p2*(psi(i+1,j+1,k)+psi(i-1,j+1,k)+psi(i,j+1,k+1)+psi(i,j+1,k-1)) &
                                              -p2*(psi(i+1,j-1,k)+psi(i-1,j-1,k)+psi(i,j-1,k+1)+psi(i,j-1,k-1)) &
                                              +p3*(psi(i+1,j+1,k+1)+psi(i-1,j+1,k+1)+psi(i+1,j+1,k-1)+psi(i-1,j+1,k-1)) &
                                              -p3*(psi(i+1,j-1,k+1)+psi(i-1,j-1,k+1)+psi(i+1,j-1,k-1)+psi(i-1,j-1,k-1)))


                          psi_z=(1.0_db/cssq)*(p1*(psi(i,j,k+1)-psi(i,j,k-1)) + &
                                              p2*(psi(i+1,j,k+1)+psi(i-1,j,k+1)+ psi(i,j+1,k+1)+psi(i,j-1,k+1)) &
                                            -p2*(psi(i+1,j,k-1)+psi(i-1,j,k-1)+psi(i,j+1,k-1)+psi(i,j-1,k-1)) &
                                            +p3*(psi(i+1,j+1,k+1)+psi(i-1,j+1,k+1)+psi(i+1,j-1,k+1)+psi(i-1,j-1,k+1)) &
                                            -p3*(psi(i+1,j+1,k-1)+psi(i-1,j+1,k-1)+psi(i+1,j-1,k-1)+psi(i-1,j-1,k-1)))
                          
                          mod_psi=sqrt(psi_x**2+psi_y**2+psi_z**2)

                          mod_psi_sq=mod_psi**2 
                          
                          rtot=0.0_db
                          
                          rtot=rhoA(i,j,k)+rhoB(i,j,k)
                          
                          rprod=rhoA(i,j,k)*rhoB(i,j,k)
                          
                          nu_avg=1.0_db/(rhoA(i,j,k)*one_ov_nu1/rtot + rhoB(i,j,k)*one_ov_nu2/rtot)
                          
                          omega=2.0_db/(6.0_db*nu_avg + 1.0_db)
                          
                          st_coeff=(9.0_db/4.0_db)*sigma*omega
                          om1=0.0_db
                          om2=0.0_db
                          if (mod_psi>0.001)then
                            !0
                            om1=-st_coeff*mod_psi*b0 !(5.0_db*st_coeff*mod_psi)/27.0
                            feq=(-4*rtot*(-2 + 3*u(i,j,k)**2 + 3*v(i,j,k)**2 + 3*w(i,j,k)**2))/27.
                            fneq1=(-3*(pxx(i,j,k) + pyy(i,j,k) + pzz(i,j,k)))/2.
                            fpc=feq + (1-omega)*fneq1*p0 + om1
                            f(i,j,k,0)=fpc*rhoA(i,j,k)/rtot
                            g(i,j,k,0)=fpc*rhoB(i,j,k)/rtot
                            !1
                            om1=st_coeff*mod_psi*(p1*psi_x**2/mod_psi_sq - b1) !(st_coeff*(psi_x - mod_psi_sq))/(27.0*mod_psi)
                            om2=p1*(rtot)*(rprod*beta*(psi_x/mod_psi)/rtot**2) !(beta*rprod/(rtot*rtot))*(p1*rtot)*psi_x/mod_psi_sq
                            feq=(rtot*(2 + 6*u(i,j,k) + 6*u(i,j,k)**2 - 3*v(i,j,k)**2 - 9*u(i,j,k)*v(i,j,k)**2 - 3*(1 + 3*u(i,j,k))*w(i,j,k)**2))/27.
                            fneq1=(3*(2*pxx(i,j,k) - pzz(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - pyy(i,j,k)*(1 + 3*u(i,j,k)) - 6*pxy(i,j,k)*v(i,j,k) - 6*pxz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq + (1-omega)*fneq1*p1 + om1 + fx*p1dcssq
                            f(i+1,j,k,1)=fpc*rhoA(i,j,k)/rtot + om2
                            g(i+1,j,k,1)=fpc*rhoB(i,j,k)/rtot - om2
                            !2
                            !om1=st_coeff*mod_psi*(p1*psi_x**2/mod_psi_sq - b1)!-(st_coeff*(psi_x + mod_psi_sq))/(27.0*mod_psi)
                            !om2=-p1*(rtot)*(rprod*beta*(psi_x/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p1*rtot)*(-psi_x/mod_psi_sq)
                            feq=(rtot*(2 - 3*v(i,j,k)**2 - 3*w(i,j,k)**2 + 3*u(i,j,k)*(-2 + 2*u(i,j,k) + 3*v(i,j,k)**2 + 3*w(i,j,k)**2)))/27.
                            fneq1=(3*(2*pxx(i,j,k) - pzz(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + pyy(i,j,k)*(-1 + 3*u(i,j,k)) + 6*pxy(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq + (1-omega)*fneq1*p1 + om1 - fx*p1dcssq
                            f(i-1,j,k,2)=fpc*rhoA(i,j,k)/rtot - om2
                            g(i-1,j,k,2)=fpc*rhoB(i,j,k)/rtot + om2
                            !3
                            om2=p1*(rtot)*(rprod*beta*(psi_y/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p1*rtot)*(psi_y/mod_psi_sq)
                            om1=st_coeff*mod_psi*(p1*psi_y**2/mod_psi_sq - b1) !(st_coeff*(psi_y - mod_psi_sq))/(27.0*mod_psi)
                            feq=(rtot*(2 - 3*u(i,j,k)**2*(1 + 3*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(2 + 2*v(i,j,k) - 3*w(i,j,k)**2)))/27.
                            fneq1=(-3*(pxx(i,j,k) - 2*pyy(i,j,k) + pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) + 3*pxx(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq+ (1-omega)*fneq1*p1 + om1 + fy*p1dcssq
                            f(i,j+1,k,3)=fpc*rhoA(i,j,k)/rtot + om2
                            g(i,j+1,k,3)=fpc*rhoB(i,j,k)/rtot - om2
                            !4
                            !om2=-p1*(rtot)*(rprod*beta*(psi_y/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p1*rtot)*(-psi_y/mod_psi_sq)
                            !om1=st_coeff*mod_psi*(p1*psi_y**2/mod_psi_sq - b1)!-(st_coeff*(psi_y + mod_psi_sq))/(27.0*mod_psi)
                            feq=(rtot*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(-2 + 2*v(i,j,k) + 3*w(i,j,k)**2)))/27.
                            fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(-1 + 3*v(i,j,k)) + 6*pyz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq+ (1-omega)*fneq1*p1 + om1 - fy*p1dcssq
                            f(i,j-1,k,4)=fpc*rhoA(i,j,k)/rtot - om2
                            g(i,j-1,k,4)=fpc*rhoB(i,j,k)/rtot + om2
                            !7
                            om1=st_coeff*mod_psi*(p2*(psi_y+psi_x)**2/mod_psi_sq - b2)!(st_coeff*(psi_x + psi_y - mod_psi_sq))/(108.0*mod_psi)
                            om2=p2*(rtot)*(rprod*beta*((psi_x+psi_y)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p2*rtot)*((psi_x + psi_y)/(sqrt(2.0)*mod_psi_sq))
                            feq=(2*rtot*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) + u(i,j,k)*(3 + 9*v(i,j,k)*(1 + v(i,j,k)))) - 3*rtot*(1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                            fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) + 6*pyy(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*(1 + 2*u(i,j,k) + 2*v(i,j,k)) + pxx(i,j,k)*(2 + 6*v(i,j,k)) - 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq + (1-omega)*fneq1*p2 + om1 + (fx+fy)*p2dcssq 
                            f(i+1,j+1,k,7)=fpc*rhoA(i,j,k)/rtot + om2
                            g(i+1,j+1,k,7)=fpc*rhoB(i,j,k)/rtot - om2
                            !8
                            !om1=st_coeff*mod_psi*(p2*(psi_y+psi_x)**2/mod_psi_sq - b2) !-(st_coeff*(psi_x + psi_y + mod_psi_sq))/(108.0*mod_psi)
                            !om2=-p2*(rtot)*(rprod*beta*((psi_x+psi_y)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p2*rtot)*(-(psi_x + psi_y)/(sqrt(2.0)*mod_psi_sq))
                            feq=(2*rtot*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(-3 - 9*(-1 + v(i,j,k))*v(i,j,k))) + 3*rtot*(-1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                            fneq1=(-3*(-2*pyy(i,j,k) + pzz(i,j,k) + 6*pyy(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*(-1 + 2*u(i,j,k) + 2*v(i,j,k)) + pxx(i,j,k)*(-2 + 6*v(i,j,k)) - 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq + (1-omega)*fneq1*p2 + om1 - (fx+fy)*p2dcssq
                            f(i-1,j-1,k,8)=fpc*rhoA(i,j,k)/rtot - om2
                            g(i-1,j-1,k,8)=fpc*rhoB(i,j,k)/rtot + om2
                            !10
                            om1=st_coeff*mod_psi*(p2*(psi_y-psi_x)**2/mod_psi_sq - b2)!-((st_coeff*(psi_x - psi_y + mod_psi_sq))/(108.0*mod_psi))
                            om2=p2*(rtot)*(rprod*beta*((-psi_x+psi_y)/mod_psi)/rtot**2) !(beta*rprod/(rtot*rtot))*(p2*rtot)*(-(psi_x - psi_y)/(sqrt(2.0)*mod_psi_sq))
                            feq=(2*rtot*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)))) + 3*rtot*(-1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                            fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) - 6*pyy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*(-1 + 2*u(i,j,k) - 2*v(i,j,k)) - 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(2 + 6*v(i,j,k)) + 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq+ (1-omega)*fneq1*p2 + om1 + (fy-fx)*p2dcssq
                            f(i-1,j+1,k,10)=fpc*rhoA(i,j,k)/rtot + om2
                            g(i-1,j+1,k,10)=fpc*rhoB(i,j,k)/rtot - om2
                            !9
                            !om2=p2*(rtot)*(rprod*beta*((psi_x-psi_y)/mod_psi)/rtot**2) !(beta*rprod/(rtot*rtot))*(p2*rtot)*((psi_x - psi_y)/(sqrt(2.0)*mod_psi_sq))
                            !om1=st_coeff*mod_psi*(p2*(psi_y-psi_x)**2/mod_psi_sq - b2)!-((st_coeff*(-psi_x + psi_y + mod_psi_sq))/(108.0*mod_psi))
                            feq=(2*rtot*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(3 + 9*(-1 + v(i,j,k))*v(i,j,k))) - 3*rtot*(1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                            fneq1=(-3*(-2*pyy(i,j,k) + pzz(i,j,k) - 6*pyy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*(1 + 2*u(i,j,k) - 2*v(i,j,k)) - 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(-2 + 6*v(i,j,k)) + 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq+ (1-omega)*fneq1*p2 + om1 + (fx-fy)*p2dcssq
                            f(i+1,j-1,k,9)=fpc*rhoA(i,j,k)/rtot - om2
                            g(i+1,j-1,k,9)=fpc*rhoB(i,j,k)/rtot + om2

                            !5
                            
                            om1=st_coeff*mod_psi*(p1*(psi_z)**2/mod_psi_sq - b1)!(st_coeff*(psi_z - mod_psi_sq))/(27.0*mod_psi)
                            om2=p1*(rtot)*(rprod*beta*((psi_z)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p1*rtot)*(psi_z/mod_psi_sq)
                            feq=(rtot*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k))))/27.
                            fneq1=(-3*(pxx(i,j,k) + pyy(i,j,k) - 2*pzz(i,j,k) + 6*pxz(i,j,k)*u(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pxx(i,j,k)*w(i,j,k) + 3*pyy(i,j,k)*w(i,j,k)))/2.
                            fpc=feq+ (1-omega)*fneq1*p1 + om1 + fz*p1dcssq
                            f(i,j,k+1,5)=fpc*rhoA(i,j,k)/rtot + om2
                            g(i,j,k+1,5)=fpc*rhoB(i,j,k)/rtot - om2
                            
                            !6
                            !om1=st_coeff*mod_psi*(p1*(psi_z)**2/mod_psi_sq - b1)!-(st_coeff*(psi_z + mod_psi_sq))/(27.0*mod_psi)
                            !om2=p1*(rtot)*(rprod*beta*((-psi_z)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p1*rtot)*(-psi_z/mod_psi_sq)
                            feq=(rtot*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                            fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 6*pxz(i,j,k)*u(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(-1 + 3*w(i,j,k))))/2.
                            fpc=feq+ (1-omega)*fneq1*p1 + om1 - fz*p1dcssq
                            f(i,j,k-1,6)=fpc*rhoA(i,j,k)/rtot - om2
                            g(i,j,k-1,6)=fpc*rhoB(i,j,k)/rtot + om2

                            !15
                            om1=st_coeff*mod_psi*(p2*(psi_z+psi_x)**2/mod_psi_sq - b2)!(st_coeff*(psi_x + psi_z - mod_psi_sq))/(108.0*mod_psi)
                            om2=p2*(rtot)*(rprod*beta*((psi_z+psi_x)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p2*rtot)*((psi_x + psi_z)/(sqrt(2.0)*mod_psi_sq))
                            feq=(rtot*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) - 3*pyy(i,j,k)*u(i,j,k) + 6*pzz(i,j,k)*u(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) - 3*pyy(i,j,k)*w(i,j,k) + 6*pxz(i,j,k)*(1 + 2*u(i,j,k) + 2*w(i,j,k)) + pxx(i,j,k)*(2 + 6*w(i,j,k))))/2.
                            fpc=feq+ (1-omega)*fneq1*p2 + om1 + (fx+fz)*p2dcssq 
                            f(i+1,j,k+1,15)=fpc*rhoA(i,j,k)/rtot + om2
                            g(i+1,j,k+1,15)=fpc*rhoB(i,j,k)/rtot - om2
                            
                            !16
                            !om1=st_coeff*mod_psi*(p2*(psi_z+psi_x)**2/mod_psi_sq - b2)!-(st_coeff*(psi_x + psi_z + mod_psi_sq))/(108.0*mod_psi)
                            !om2=-p2*(rtot)*(rprod*beta*((psi_z+psi_x)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p2*rtot)*(-(psi_x + psi_z)/(sqrt(2.0)*mod_psi_sq))
                            feq=(rtot*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(2 - 6*w(i,j,k)) + 3*pyy(i,j,k)*w(i,j,k) - 6*pxz(i,j,k)*(-1 + 2*u(i,j,k) + 2*w(i,j,k))))/2.
                            fpc=feq+ (1-omega)*fneq1*p2 + om1 - (fx+fz)*p2dcssq
                            f(i-1,j,k-1,16)=fpc*rhoA(i,j,k)/rtot - om2
                            g(i-1,j,k-1,16)=fpc*rhoB(i,j,k)/rtot + om2

                            !17
                            om1=st_coeff*mod_psi*(p2*(psi_z-psi_x)**2/mod_psi_sq - b2)!-((st_coeff*(psi_x - psi_z + mod_psi_sq))/(108.0*mod_psi))
                            om2=p2*(rtot)*(rprod*beta*((psi_z-psi_x)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p2*rtot)*((-psi_x + psi_z)/(sqrt(2.0)*mod_psi_sq))
                            feq=(rtot*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*(-1 + 2*u(i,j,k) - 2*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(2 + 6*w(i,j,k))))/2.
                            fpc=feq+ (1-omega)*fneq1*p2 + om1 + (fz-fx)*p2dcssq
                            f(i-1,j,k+1,17)=fpc*rhoA(i,j,k)/rtot + om2
                            g(i-1,j,k+1,17)=fpc*rhoB(i,j,k)/rtot - om2
                            
                            !18
                            !om1=st_coeff*mod_psi*(p2*(psi_z-psi_x)**2/mod_psi_sq - b2)!-((st_coeff*(-psi_x + psi_z + mod_psi_sq))/(108.0*mod_psi))
                            !om2=p2*(rtot)*(rprod*beta*((-psi_z+psi_x)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p2*rtot)*((psi_x - psi_z)/(sqrt(2.0)*mod_psi_sq))
                            feq=(rtot*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(-3*(pyy(i,j,k) - 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*(1 + 2*u(i,j,k) - 2*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(-2 + 6*w(i,j,k))))/2.
                            fpc=feq+ (1-omega)*fneq1*p2 + om1 + (fx-fz)*p2dcssq
                            f(i+1,j,k-1,18)=fpc*rhoA(i,j,k)/rtot - om2
                            g(i+1,j,k-1,18)=fpc*rhoB(i,j,k)/rtot + om2

                            !11
                            om1=st_coeff*mod_psi*(p2*(psi_z+psi_y)**2/mod_psi_sq - b2)!!(st_coeff*(psi_y + psi_z - mod_psi_sq))/(108.0*mod_psi)
                            om2=p2*(rtot)*(rprod*beta*((psi_z+psi_y)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p2*rtot)*((psi_y + psi_z)/(sqrt(2.0)*mod_psi_sq))
                            feq=(rtot*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(-3*pxx(i,j,k)*(1 + 3*v(i,j,k) + 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) - 3*pxy(i,j,k)*u(i,j,k) - 3*pxz(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*pyy(i,j,k)*w(i,j,k) + pyz(i,j,k)*(3 + 6*v(i,j,k) + 6*w(i,j,k)))
                            fpc=feq+ (1-omega)*fneq1*p2 + om1 + (fy+fz)*p2dcssq
                            f(i,j+1,k+1,11)=fpc*rhoA(i,j,k)/rtot + om2
                            g(i,j+1,k+1,11)=fpc*rhoB(i,j,k)/rtot - om2
                            
                            !12
                            !om1=st_coeff*mod_psi*(p2*(psi_z+psi_y)**2/mod_psi_sq - b2)!-(st_coeff*(psi_y + psi_z + mod_psi_sq))/(108.0*mod_psi)
                            !om2=-p2*(rtot)*(rprod*beta*((psi_z+psi_y)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p2*rtot)*(-(psi_y + psi_z)/(sqrt(2.0)*mod_psi_sq))
                            feq=(rtot*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*pxx(i,j,k)*(-1 + 3*v(i,j,k) + 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) + 3*pxy(i,j,k)*u(i,j,k) + 3*pxz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + pyz(i,j,k)*(3 - 6*v(i,j,k) - 6*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k))
                            fpc=feq+ (1-omega)*fneq1*p2 + om1 - (fy+fz)*p2dcssq
                            f(i,j-1,k-1,12)=fpc*rhoA(i,j,k)/rtot - om2
                            g(i,j-1,k-1,12)=fpc*rhoB(i,j,k)/rtot + om2

                            !13
                            om1=st_coeff*mod_psi*(p2*(-psi_z+psi_y)**2/mod_psi_sq - b2)!-((st_coeff*(-psi_y + psi_z + mod_psi_sq))/(108.0*mod_psi))
                            om2=p2*(rtot)*(rprod*beta*((-psi_z+psi_y)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p2*rtot)*((psi_y - psi_z)/(sqrt(2.0)*mod_psi_sq))
                            feq=(rtot*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                            fneq1=(-3*(pxx(i,j,k) - 2*pyy(i,j,k) + 6*pyz(i,j,k) - 2*pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) - 6*pxz(i,j,k)*u(i,j,k) + 3*pxx(i,j,k)*v(i,j,k) + 12*pyz(i,j,k)*v(i,j,k) - 6*pzz(i,j,k)*v(i,j,k) - 3*pxx(i,j,k)*w(i,j,k) + 6*pyy(i,j,k)*w(i,j,k) - 12*pyz(i,j,k)*w(i,j,k)))/2.
                            fpc= feq+ (1-omega)*fneq1*p2 + om1 + (fy-fz)*p2dcssq
                            f(i,j+1,k-1,13)=fpc*rhoA(i,j,k)/rtot + om2
                            g(i,j+1,k-1,13)=fpc*rhoB(i,j,k)/rtot - om2
                            
                            !14
                            !om1=st_coeff*mod_psi*(p2*(-psi_z+psi_y)**2/mod_psi_sq - b2)!-((st_coeff*(psi_y - psi_z + mod_psi_sq))/(108.0*mod_psi))
                            !om2=p2*(rtot)*(rprod*beta*((psi_z-psi_y)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p2*rtot)*(-(psi_y - psi_z)/(sqrt(2.0)*mod_psi_sq))
                            feq=(rtot*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*pxx(i,j,k)*(-1 + 3*v(i,j,k) - 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) + 3*pxy(i,j,k)*u(i,j,k) - 3*pxz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + pyz(i,j,k)*(-3 + 6*v(i,j,k) - 6*w(i,j,k)) + 3*pyy(i,j,k)*w(i,j,k))
                            fpc=feq+ (1-omega)*fneq1*p2 + om1 + (fz-fy)*p2dcssq
                            f(i,j-1,k+1,14)=fpc*rhoA(i,j,k)/rtot - om2
                            g(i,j-1,k+1,14)=fpc*rhoB(i,j,k)/rtot + om2

                            !19
                            om1=st_coeff*mod_psi*(p3*(psi_z+psi_y+psi_x)**2/mod_psi_sq - b3)!(st_coeff*(psi_x + psi_y + psi_z - mod_psi_sq))/(432.0*mod_psi)
                            om2=p3*(rtot)*(rprod*beta*((psi_z+psi_y+psi_x)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p3*rtot)*((psi_x + psi_y + psi_z)/(sqrt(3.0)*mod_psi_sq))
                            feq=(rtot*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i,j,k) + (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) + 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                            fpc=feq + (1-omega)*fneq1*p3 + om1 + (fz+fy+fx)*p3dcssq
                            f(i+1,j+1,k+1,19)=fpc*rhoA(i,j,k)/rtot + om2
                            g(i+1,j+1,k+1,19)=fpc*rhoB(i,j,k)/rtot - om2
                            
                            !20
                            !om1=st_coeff*mod_psi*(p3*(psi_z+psi_y+psi_x)**2/mod_psi_sq - b3)!-(st_coeff*(psi_x + psi_y + psi_z + mod_psi_sq))/(432.0*mod_psi)
                            !om2=-p3*(rtot)*(rprod*beta*((psi_z+psi_y+psi_x)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p3*rtot)*(-(psi_x + psi_y + psi_z)/(sqrt(3.0)*mod_psi_sq))
                            feq=(rtot*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*((pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(-1 + 3*u(i,j,k)) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + pxz(i,j,k)*(-3 + 6*u(i,j,k)) + 3*(2*pxy(i,j,k) + 3*pxz(i,j,k) + 2*pyz(i,j,k) + pzz(i,j,k))*v(i,j,k) + 3*(3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k) + pxx(i,j,k)*(-1 + 3*v(i,j,k) + 3*w(i,j,k)))
                            fpc=feq+ (1-omega)*fneq1*p3 + om1 - (fz+fy+fx)*p3dcssq
                            f(i-1,j-1,k-1,20)=fpc*rhoA(i,j,k)/rtot - om2
                            g(i-1,j-1,k-1,20)=fpc*rhoB(i,j,k)/rtot + om2

                            !21
                            om1=st_coeff*mod_psi*(p3*(psi_z-psi_y+psi_x)**2/mod_psi_sq - b3)!(st_coeff*(psi_x - psi_y + psi_z - mod_psi_sq))/(432.0*mod_psi)
                            om2=p3*(rtot)*(rprod*beta*((psi_z-psi_y+psi_x)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p3*rtot)*((psi_x - psi_y + psi_z)/(sqrt(3.0)*mod_psi_sq))
                            feq=(rtot*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i,j,k) - 3*pxy(i,j,k)*(1 + 2*u(i,j,k)) + (pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) - 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                            fpc=feq+ (1-omega)*fneq1*p3 + om1 + (fx-fy+fz)*p3dcssq
                            f(i+1,j-1,k+1,21)=fpc*rhoA(i,j,k)/rtot + om2
                            g(i+1,j-1,k+1,21)=fpc*rhoB(i,j,k)/rtot - om2
                            
                            !22
                            !om1=st_coeff*mod_psi*(p3*(psi_z-psi_y+psi_x)**2/mod_psi_sq - b3)!-(st_coeff*(psi_x - psi_y + psi_z + mod_psi_sq))/(432.0*mod_psi)
                            !om2=p3*(rtot)*(rprod*beta*((-psi_z+psi_y-psi_x)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p3*rtot)*(-(psi_x - psi_y + psi_z)/(sqrt(3.0)*mod_psi_sq))
                            feq=(rtot*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i,j,k) + 3*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k) - 3*(2*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) - 3*(pxx(i,j,k) - 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                            fpc=feq+ (1-omega)*fneq1*p3 + om1 - (fx-fy+fz)*p3dcssq
                            f(i-1,j+1,k-1,22)=fpc*rhoA(i,j,k)/rtot - om2
                            g(i-1,j+1,k-1,22)=fpc*rhoB(i,j,k)/rtot + om2

                            !23
                            om1=st_coeff*mod_psi*(p3*(psi_z-psi_y-psi_x)**2/mod_psi_sq - b3)!-((st_coeff*(psi_x + psi_y - psi_z + mod_psi_sq))/(432.0*mod_psi))
                            om2=p3*(rtot)*(rprod*beta*((psi_z-psi_y-psi_x)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p3*rtot)*(-(psi_x + psi_y - psi_z)/(sqrt(3.0)*mod_psi_sq))
                            feq=(rtot*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i,j,k) + 3*pxy(i,j,k) - 3*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k) - 3*(2*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k) - 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) + 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                            fpc=feq+ (1-omega)*fneq1*p3 + om1 + (fz-fy-fx)*p3dcssq
                            f(i-1,j-1,k+1,23)=fpc*rhoA(i,j,k)/rtot + om2
                            g(i-1,j-1,k+1,23)=fpc*rhoB(i,j,k)/rtot - om2

                            !24
                            !om1=st_coeff*mod_psi*(p3*(psi_z-psi_y-psi_x)**2/mod_psi_sq - b3)!(st_coeff*(psi_x + psi_y - psi_z - mod_psi_sq))/(432.0*mod_psi)
                            !om2=p3*(rtot)*(rprod*beta*((-psi_z+psi_y+psi_x)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p3*rtot)*((psi_x + psi_y - psi_z)/(sqrt(3.0)*mod_psi_sq))
                            feq=(rtot*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                            fneq1=3*(pxx(i,j,k) - 3*pxz(i,j,k)*(1 + 2*u(i,j,k)) + (pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) - 3*(pxx(i,j,k) + 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                            fpc=feq+ (1-omega)*fneq1*p3 + om1 + (-fz+fy+fx)*p3dcssq
                            f(i+1,j+1,k-1,24)=fpc*rhoA(i,j,k)/rtot - om2
                            g(i+1,j+1,k-1,24)=fpc*rhoB(i,j,k)/rtot + om2

                            !25
                            om1=st_coeff*mod_psi*(p3*(-psi_z-psi_y+psi_x)**2/mod_psi_sq - b3)!-((st_coeff*(-psi_x + psi_y + psi_z + mod_psi_sq))/(432.0*mod_psi))
                            om2=p3*(rtot)*(rprod*beta*((-psi_z-psi_y+psi_x)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p3*rtot)*((psi_x - psi_y - psi_z)/(sqrt(3.0)*mod_psi_sq))
                            feq=(rtot*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*(-pxx(i,j,k) - (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                            fpc=feq + (1-omega)*fneq1*p3 + om1 + (-fz-fy+fx)*p3dcssq
                            f(i+1,j-1,k-1,25)=fpc*rhoA(i,j,k)/rtot + om2
                            g(i+1,j-1,k-1,25)=fpc*rhoB(i,j,k)/rtot - om2

                            !26
                            !om1=st_coeff*mod_psi*(p3*(-psi_z - psi_y + psi_x)**2/mod_psi_sq - b3)!((st_coeff*(-psi_x + psi_y + psi_z - mod_psi_sq))/(432.0*mod_psi))
                            !om2=p3*(rtot)*(rprod*beta*((psi_z+psi_y-psi_x)/mod_psi)/rtot**2)!(beta*rprod/(rtot*rtot))*(p3*rtot)*(-(psi_x - psi_y - psi_z)/(sqrt(3.0)*mod_psi_sq))
                            feq=(rtot*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i,j,k) - (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(-1 + 3*u(i,j,k)) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + pxz(i,j,k)*(-3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                            fpc=feq + (1-omega)*fneq1*p3 + om1 + (fz+fy-fx)*p3dcssq
                            f(i-1,j+1,k+1,26)=fpc*rhoA(i,j,k)/rtot - om2
                            g(i-1,j+1,k+1,26)=fpc*rhoB(i,j,k)/rtot + om2
                          else
                            !0
                            feq=(-4*rtot*(-2 + 3*u(i,j,k)**2 + 3*v(i,j,k)**2 + 3*w(i,j,k)**2))/27.
                            fneq1=(-3*(pxx(i,j,k) + pyy(i,j,k) + pzz(i,j,k)))/2.
                            fpc=feq + (1-omega)*fneq1*p0 
                            f(i,j,k,0)=fpc*rhoA(i,j,k)/rtot
                            g(i,j,k,0)=fpc*rhoB(i,j,k)/rtot
                            !1
                            feq=(rtot*(2 + 6*u(i,j,k) + 6*u(i,j,k)**2 - 3*v(i,j,k)**2 - 9*u(i,j,k)*v(i,j,k)**2 - 3*(1 + 3*u(i,j,k))*w(i,j,k)**2))/27.
                            fneq1=(3*(2*pxx(i,j,k) - pzz(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - pyy(i,j,k)*(1 + 3*u(i,j,k)) - 6*pxy(i,j,k)*v(i,j,k) - 6*pxz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq + (1-omega)*fneq1*p1 + fx*p1dcssq
                            f(i+1,j,k,1)=fpc*rhoA(i,j,k)/rtot 
                            g(i+1,j,k,1)=fpc*rhoB(i,j,k)/rtot 
                            !2
                            feq=(rtot*(2 - 3*v(i,j,k)**2 - 3*w(i,j,k)**2 + 3*u(i,j,k)*(-2 + 2*u(i,j,k) + 3*v(i,j,k)**2 + 3*w(i,j,k)**2)))/27.
                            fneq1=(3*(2*pxx(i,j,k) - pzz(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + pyy(i,j,k)*(-1 + 3*u(i,j,k)) + 6*pxy(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq + (1-omega)*fneq1*p1  - fx*p1dcssq
                            f(i-1,j,k,2)=fpc*rhoA(i,j,k)/rtot
                            g(i-1,j,k,2)=fpc*rhoB(i,j,k)/rtot
                            !3
                            feq=(rtot*(2 - 3*u(i,j,k)**2*(1 + 3*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(2 + 2*v(i,j,k) - 3*w(i,j,k)**2)))/27.
                            fneq1=(-3*(pxx(i,j,k) - 2*pyy(i,j,k) + pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) + 3*pxx(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq + (1-omega)*fneq1*p1 + fy*p1dcssq
                            f(i,j+1,k,3)=fpc*rhoA(i,j,k)/rtot
                            g(i,j+1,k,3)=fpc*rhoB(i,j,k)/rtot
                            !4
                            feq=(rtot*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(-2 + 2*v(i,j,k) + 3*w(i,j,k)**2)))/27.
                            fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(-1 + 3*v(i,j,k)) + 6*pyz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq+ (1-omega)*fneq1*p1 - fy*p1dcssq
                            f(i,j-1,k,4)=fpc*rhoA(i,j,k)/rtot 
                            g(i,j-1,k,4)=fpc*rhoB(i,j,k)/rtot 
                            !7
                            feq=(2*rtot*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) + u(i,j,k)*(3 + 9*v(i,j,k)*(1 + v(i,j,k)))) - 3*rtot*(1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                            fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) + 6*pyy(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*(1 + 2*u(i,j,k) + 2*v(i,j,k)) + pxx(i,j,k)*(2 + 6*v(i,j,k)) - 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq + (1-omega)*fneq1*p2 + (fx+fy)*p2dcssq 
                            f(i+1,j+1,k,7)=fpc*rhoA(i,j,k)/rtot 
                            g(i+1,j+1,k,7)=fpc*rhoB(i,j,k)/rtot 
                            !8
                            feq=(2*rtot*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(-3 - 9*(-1 + v(i,j,k))*v(i,j,k))) + 3*rtot*(-1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                            fneq1=(-3*(-2*pyy(i,j,k) + pzz(i,j,k) + 6*pyy(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*(-1 + 2*u(i,j,k) + 2*v(i,j,k)) + pxx(i,j,k)*(-2 + 6*v(i,j,k)) - 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq + (1-omega)*fneq1*p2 - (fx+fy)*p2dcssq
                            f(i-1,j-1,k,8)=fpc*rhoA(i,j,k)/rtot 
                            g(i-1,j-1,k,8)=fpc*rhoB(i,j,k)/rtot 
                            !10
                            feq=(2*rtot*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)))) + 3*rtot*(-1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                            fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) - 6*pyy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*(-1 + 2*u(i,j,k) - 2*v(i,j,k)) - 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(2 + 6*v(i,j,k)) + 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq+ (1-omega)*fneq1*p2 + (fy-fx)*p2dcssq
                            f(i-1,j+1,k,10)=fpc*rhoA(i,j,k)/rtot 
                            g(i-1,j+1,k,10)=fpc*rhoB(i,j,k)/rtot 
                            !9
                            feq=(2*rtot*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(3 + 9*(-1 + v(i,j,k))*v(i,j,k))) - 3*rtot*(1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                            fneq1=(-3*(-2*pyy(i,j,k) + pzz(i,j,k) - 6*pyy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*(1 + 2*u(i,j,k) - 2*v(i,j,k)) - 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(-2 + 6*v(i,j,k)) + 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                            fpc=feq+ (1-omega)*fneq1*p2 + (fx-fy)*p2dcssq
                            f(i+1,j-1,k,9)=fpc*rhoA(i,j,k)/rtot 
                            g(i+1,j-1,k,9)=fpc*rhoB(i,j,k)/rtot 

                            !5
                            feq=(rtot*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k))))/27.
                            fneq1=(-3*(pxx(i,j,k) + pyy(i,j,k) - 2*pzz(i,j,k) + 6*pxz(i,j,k)*u(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pxx(i,j,k)*w(i,j,k) + 3*pyy(i,j,k)*w(i,j,k)))/2.
                            fpc=feq+ (1-omega)*fneq1*p1 + fz*p1dcssq
                            f(i,j,k+1,5)=fpc*rhoA(i,j,k)/rtot 
                            g(i,j,k+1,5)=fpc*rhoB(i,j,k)/rtot 
                            
                            !6
                            feq=(rtot*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                            fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 6*pxz(i,j,k)*u(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(-1 + 3*w(i,j,k))))/2.
                            fpc=feq+ (1-omega)*fneq1*p1 - fz*p1dcssq
                            f(i,j,k-1,6)=fpc*rhoA(i,j,k)/rtot 
                            g(i,j,k-1,6)=fpc*rhoB(i,j,k)/rtot 

                            !15
                            feq=(rtot*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) - 3*pyy(i,j,k)*u(i,j,k) + 6*pzz(i,j,k)*u(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) - 3*pyy(i,j,k)*w(i,j,k) + 6*pxz(i,j,k)*(1 + 2*u(i,j,k) + 2*w(i,j,k)) + pxx(i,j,k)*(2 + 6*w(i,j,k))))/2.
                            fpc=feq+ (1-omega)*fneq1*p2 + (fx+fz)*p2dcssq 
                            f(i+1,j,k+1,15)=fpc*rhoA(i,j,k)/rtot 
                            g(i+1,j,k+1,15)=fpc*rhoB(i,j,k)/rtot 
                            
                            !16
                            feq=(rtot*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(2 - 6*w(i,j,k)) + 3*pyy(i,j,k)*w(i,j,k) - 6*pxz(i,j,k)*(-1 + 2*u(i,j,k) + 2*w(i,j,k))))/2.
                            fpc=feq+ (1-omega)*fneq1*p2 - (fx+fz)*p2dcssq
                            f(i-1,j,k-1,16)=fpc*rhoA(i,j,k)/rtot 
                            g(i-1,j,k-1,16)=fpc*rhoB(i,j,k)/rtot 

                            !17
                            feq=(rtot*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*(-1 + 2*u(i,j,k) - 2*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(2 + 6*w(i,j,k))))/2.
                            fpc=feq+ (1-omega)*fneq1*p2 + (fz-fx)*p2dcssq
                            f(i-1,j,k+1,17)=fpc*rhoA(i,j,k)/rtot 
                            g(i-1,j,k+1,17)=fpc*rhoB(i,j,k)/rtot
                            
                            !18
                            feq=(rtot*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(-3*(pyy(i,j,k) - 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*(1 + 2*u(i,j,k) - 2*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(-2 + 6*w(i,j,k))))/2.
                            fpc=feq+ (1-omega)*fneq1*p2 + (fx-fz)*p2dcssq
                            f(i+1,j,k-1,18)=fpc*rhoA(i,j,k)/rtot
                            g(i+1,j,k-1,18)=fpc*rhoB(i,j,k)/rtot 

                            !11
                            feq=(rtot*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(-3*pxx(i,j,k)*(1 + 3*v(i,j,k) + 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) - 3*pxy(i,j,k)*u(i,j,k) - 3*pxz(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*pyy(i,j,k)*w(i,j,k) + pyz(i,j,k)*(3 + 6*v(i,j,k) + 6*w(i,j,k)))
                            fpc=feq+ (1-omega)*fneq1*p2 + (fy+fz)*p2dcssq
                            f(i,j+1,k+1,11)=fpc*rhoA(i,j,k)/rtot 
                            g(i,j+1,k+1,11)=fpc*rhoB(i,j,k)/rtot 
                            
                            !12
                            feq=(rtot*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*pxx(i,j,k)*(-1 + 3*v(i,j,k) + 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) + 3*pxy(i,j,k)*u(i,j,k) + 3*pxz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + pyz(i,j,k)*(3 - 6*v(i,j,k) - 6*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k))
                            fpc=feq+ (1-omega)*fneq1*p2 - (fy+fz)*p2dcssq
                            f(i,j-1,k-1,12)=fpc*rhoA(i,j,k)/rtot
                            g(i,j-1,k-1,12)=fpc*rhoB(i,j,k)/rtot 

                            !13
                            feq=(rtot*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                            fneq1=(-3*(pxx(i,j,k) - 2*pyy(i,j,k) + 6*pyz(i,j,k) - 2*pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) - 6*pxz(i,j,k)*u(i,j,k) + 3*pxx(i,j,k)*v(i,j,k) + 12*pyz(i,j,k)*v(i,j,k) - 6*pzz(i,j,k)*v(i,j,k) - 3*pxx(i,j,k)*w(i,j,k) + 6*pyy(i,j,k)*w(i,j,k) - 12*pyz(i,j,k)*w(i,j,k)))/2.
                            fpc= feq+ (1-omega)*fneq1*p2 + (fy-fz)*p2dcssq
                            f(i,j+1,k-1,13)=fpc*rhoA(i,j,k)/rtot 
                            g(i,j+1,k-1,13)=fpc*rhoB(i,j,k)/rtot 
                            
                            !14
                            feq=(rtot*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*pxx(i,j,k)*(-1 + 3*v(i,j,k) - 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) + 3*pxy(i,j,k)*u(i,j,k) - 3*pxz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + pyz(i,j,k)*(-3 + 6*v(i,j,k) - 6*w(i,j,k)) + 3*pyy(i,j,k)*w(i,j,k))
                            fpc=feq+ (1-omega)*fneq1*p2 + (fz-fy)*p2dcssq
                            f(i,j-1,k+1,14)=fpc*rhoA(i,j,k)/rtot 
                            g(i,j-1,k+1,14)=fpc*rhoB(i,j,k)/rtot 

                            !19
                            feq=(rtot*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i,j,k) + (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) + 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                            fpc=feq + (1-omega)*fneq1*p3  + (fz+fy+fx)*p3dcssq
                            f(i+1,j+1,k+1,19)=fpc*rhoA(i,j,k)/rtot 
                            g(i+1,j+1,k+1,19)=fpc*rhoB(i,j,k)/rtot 
                            
                            !20
                            feq=(rtot*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*((pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(-1 + 3*u(i,j,k)) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + pxz(i,j,k)*(-3 + 6*u(i,j,k)) + 3*(2*pxy(i,j,k) + 3*pxz(i,j,k) + 2*pyz(i,j,k) + pzz(i,j,k))*v(i,j,k) + 3*(3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k) + pxx(i,j,k)*(-1 + 3*v(i,j,k) + 3*w(i,j,k)))
                            fpc=feq+ (1-omega)*fneq1*p3  - (fz+fy+fx)*p3dcssq
                            f(i-1,j-1,k-1,20)=fpc*rhoA(i,j,k)/rtot 
                            g(i-1,j-1,k-1,20)=fpc*rhoB(i,j,k)/rtot 

                            !21
                            feq=(rtot*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i,j,k) - 3*pxy(i,j,k)*(1 + 2*u(i,j,k)) + (pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) - 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                            fpc=feq+ (1-omega)*fneq1*p3 + (fx-fy+fz)*p3dcssq
                            f(i+1,j-1,k+1,21)=fpc*rhoA(i,j,k)/rtot 
                            g(i+1,j-1,k+1,21)=fpc*rhoB(i,j,k)/rtot 
                            
                            !22
                            feq=(rtot*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i,j,k) + 3*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k) - 3*(2*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) - 3*(pxx(i,j,k) - 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                            fpc=feq+ (1-omega)*fneq1*p3 - (fx-fy+fz)*p3dcssq
                            f(i-1,j+1,k-1,22)=fpc*rhoA(i,j,k)/rtot 
                            g(i-1,j+1,k-1,22)=fpc*rhoB(i,j,k)/rtot 

                            !23
                            feq=(rtot*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i,j,k) + 3*pxy(i,j,k) - 3*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k) - 3*(2*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k) - 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) + 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                            fpc=feq+ (1-omega)*fneq1*p3 + (fz-fy-fx)*p3dcssq
                            f(i-1,j-1,k+1,23)=fpc*rhoA(i,j,k)/rtot 
                            g(i-1,j-1,k+1,23)=fpc*rhoB(i,j,k)/rtot 

                            !24
                            feq=(rtot*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                            fneq1=3*(pxx(i,j,k) - 3*pxz(i,j,k)*(1 + 2*u(i,j,k)) + (pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) - 3*(pxx(i,j,k) + 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                            fpc=feq+ (1-omega)*fneq1*p3 + (-fz+fy+fx)*p3dcssq
                            f(i+1,j+1,k-1,24)=fpc*rhoA(i,j,k)/rtot 
                            g(i+1,j+1,k-1,24)=fpc*rhoB(i,j,k)/rtot

                            !25
                            feq=(rtot*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*(-pxx(i,j,k) - (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                            fpc=feq + (1-omega)*fneq1*p3 + (-fz-fy+fx)*p3dcssq
                            f(i+1,j-1,k-1,25)=fpc*rhoA(i,j,k)/rtot
                            g(i+1,j-1,k-1,25)=fpc*rhoB(i,j,k)/rtot 

                            !26
                            feq=(rtot*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i,j,k) - (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(-1 + 3*u(i,j,k)) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + pxz(i,j,k)*(-3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                            fpc=feq + (1-omega)*fneq1*p3 + (fz+fy-fx)*p3dcssq
                            f(i-1,j+1,k+1,26)=fpc*rhoA(i,j,k)/rtot 
                            g(i-1,j+1,k+1,26)=fpc*rhoB(i,j,k)/rtot 
                          endif
                      endif
                  enddo
              enddo
          enddo
          !$acc end kernels
        !***********************************boundary conditions********************************!
        call  bcs_allbback_multicomp          
    enddo 
    !$acc end data
    call cpu_time(ts2)
    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
    write(6,*) 'glups: ',  real(nx)*real(ny)*real(nz)*real(nsteps)/1.0e9/(ts2-ts1)
    
    call get_memory_gpu(mymemory,totmemory)
    call print_memory_registration_gpu(6,'DEVICE memory occupied at the end', &
     'total DEVICE memory',mymemory,totmemory)
    
end program
