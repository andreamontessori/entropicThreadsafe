program recursiveTSLB3D
    !$if _OPENACC
    use openacc
    !$endif
    use prints
    use vars
    
    implicit none
    
    !$if _OPENACC
    integer :: devNum
    integer(acc_device_kind) :: devType
    devType = acc_get_device_type()
    devNum=acc_get_device_num(devType)
    !$endif

    nlinks=26 !pari!
    tau=0.6_db
    cssq=1.0_db/3.0_db
    visc_LB=cssq*(tau-0.5_db)
    one_ov_nu=1.0_db/visc_LB
#ifdef _OPENACC
    ngpus=acc_get_num_devices(acc_device_nvidia)
#else
    ngpus=0
#endif

    !*******************************user parameters and allocations**************************m
        nx=64
        ny=64
        nz=64
        nsteps=20000
        stamp=20000
        fx=0.0_db*10.0**(-5)
        fy=1.0_db*10.0**(-5)
        fz=0.0_db*10.0**(-5)
        lprint=.true.
        lvtk=.true.
        lasync=.false.
        lpbc=.true.
        lpbcfull=.false.
        
        allocate(f(0:nx+1,0:ny+1,0:nz+1,0:nlinks))
        allocate(rho(1:nx,1:ny,1:nz),u(1:nx,1:ny,1:nz),v(1:nx,1:ny,1:nz),w(1:nx,1:ny,1:nz))
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
        p1dcssq=p1/cssq
        p2dcssq=p2/cssq
        p3dcssq=p3/cssq
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
        !do ll=0,nlinks
        f(1:nx,1:ny,1:nz,0)=rho(1:nx,1:ny,1:nz)*p0
        do ll=1,6
          f(1:nx,1:ny,1:nz,ll)=rho(1:nx,1:ny,1:nz)*p1
        enddo
        do ll=7,18
          f(1:nx,1:ny,1:nz,ll)=rho(1:nx,1:ny,1:nz)*p2
        enddo
        do ll=19,26
          f(1:nx,1:ny,1:nz,ll)=rho(1:nx,1:ny,1:nz)*p3
        enddo
    !enddo
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
        write(6,*) 'max fx',huge(fx)
        write(6,*) 'max fx',huge(fy)
        write(6,*) 'max fx',huge(fz)
        write(6,*) '*******************************************'
    !$acc data copy(f,isfluid,p0,p1,p2,p3,&
             !$acc& pxx,pyy,pzz,pxy,pxz,pyz,rho,u,v,w,rhoprint,velprint)
    !$if _OPENACC        
        call printDeviceProperties(ngpus,devNum,devType,6)
    !$endif
    iframe=0
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
          !$acc loop collapse(3) private(fneq1,feq,temp,uu,udotc)
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
                          rho(i,j,k) = f(i,j,k,0)+f(i,j,k,1)+f(i,j,k,2)+f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,5) &
                              +f(i,j,k,6)+f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11) &
                              +f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17) &
                              +f(i,j,k,18)+f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,24) &
                              +f(i,j,k,25) +f(i,j,k,26)
                              

                          u(i,j,k) = ((f(i,j,k,1)+f(i,j,k,7)+f(i,j,k,9)+f(i,j,k,15)+f(i,j,k,18)++f(i,j,k,19)+f(i,j,k,21)+f(i,j,k,24)+f(i,j,k,25)) &
                              -(f(i,j,k,2)+f(i,j,k,8)+f(i,j,k,10)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,20)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,26)))/rho(i,j,k)
                          
                          v(i,j,k) = ((f(i,j,k,3)+f(i,j,k,7)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,13)+f(i,j,k,19)+f(i,j,k,22)+f(i,j,k,24)+f(i,j,k,26)) &
                              -(f(i,j,k,4)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,12)+f(i,j,k,14)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,23)+f(i,j,k,25)))/rho(i,j,k)

                          w(i,j,k) = ((f(i,j,k,5)+f(i,j,k,11)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,17)+f(i,j,k,19)+f(i,j,k,21)+f(i,j,k,23)+f(i,j,k,26)) &
                              -(f(i,j,k,6)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,16)+f(i,j,k,18)+f(i,j,k,20)+f(i,j,k,22)+f(i,j,k,24)+f(i,j,k,25)))/rho(i,j,k)                        
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
          !$acc end kernels
        !***********************************Print on files************************
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
        
        !***********************************collision + no slip + forcing: fused implementation*********
         !$acc kernels
          !$acc loop collapse(3) private(feq,uu,temp,udotc)
          do k=1,nz
              do j=1,ny
                  do i=1,nx
                      if(isfluid(i,j,k).eq.1)then
                          !0
                          feq=(-4*rho(i,j,k)*(-2 + 3*u(i,j,k)**2 + 3*v(i,j,k)**2 + 3*w(i,j,k)**2))/27.
                          fneq1=(-3*(pxx(i,j,k) + pyy(i,j,k) + pzz(i,j,k)))/2.
                          f(i,j,k,0)=feq + (1-omega)*fneq1*p0
                          
                          !1
                          
                          feq=(rho(i,j,k)*(2 + 6*u(i,j,k) + 6*u(i,j,k)**2 - 3*v(i,j,k)**2 - 9*u(i,j,k)*v(i,j,k)**2 - 3*(1 + 3*u(i,j,k))*w(i,j,k)**2))/27.
                          fneq1=(3*(2*pxx(i,j,k) - pzz(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - pyy(i,j,k)*(1 + 3*u(i,j,k)) - 6*pxy(i,j,k)*v(i,j,k) - 6*pxz(i,j,k)*w(i,j,k)))/2.
                          f(i+1,j,k,1)=feq + (1-omega)*fneq1*p1 + fx*p1dcssq
                          
                          !2
                          feq=(rho(i,j,k)*(2 - 3*v(i,j,k)**2 - 3*w(i,j,k)**2 + 3*u(i,j,k)*(-2 + 2*u(i,j,k) + 3*v(i,j,k)**2 + 3*w(i,j,k)**2)))/27.
                          fneq1=(3*(2*pxx(i,j,k) - pzz(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + pyy(i,j,k)*(-1 + 3*u(i,j,k)) + 6*pxy(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*w(i,j,k)))/2.
                          f(i-1,j,k,2)=feq + (1-omega)*fneq1*p1 - fx*p1dcssq
                          
                          !3
                          
                          feq=(rho(i,j,k)*(2 - 3*u(i,j,k)**2*(1 + 3*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(2 + 2*v(i,j,k) - 3*w(i,j,k)**2)))/27.
                          fneq1=(-3*(pxx(i,j,k) - 2*pyy(i,j,k) + pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) + 3*pxx(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i,j+1,k,3)=feq+ (1-omega)*fneq1*p1 + fy*p1dcssq
                          
                          !4
                          feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(-2 + 2*v(i,j,k) + 3*w(i,j,k)**2)))/27.
                          fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(-1 + 3*v(i,j,k)) + 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i,j-1,k,4)=feq+ (1-omega)*fneq1*p1 - fy*p1dcssq
                          
                          !7
                          
                          feq=(2*rho(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) + u(i,j,k)*(3 + 9*v(i,j,k)*(1 + v(i,j,k)))) - 3*rho(i,j,k)*(1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) + 6*pyy(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*(1 + 2*u(i,j,k) + 2*v(i,j,k)) + pxx(i,j,k)*(2 + 6*v(i,j,k)) - 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i+1,j+1,k,7)=feq + (1-omega)*fneq1*p2 + (fx+fy)*p2dcssq 
                          
                          !8
                          feq=(2*rho(i,j,k)*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(-3 - 9*(-1 + v(i,j,k))*v(i,j,k))) + 3*rho(i,j,k)*(-1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=(-3*(-2*pyy(i,j,k) + pzz(i,j,k) + 6*pyy(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*(-1 + 2*u(i,j,k) + 2*v(i,j,k)) + pxx(i,j,k)*(-2 + 6*v(i,j,k)) - 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i-1,j-1,k,8)=feq + (1-omega)*fneq1*p2 - (fx+fy)*p2dcssq
                          
                          !10
                          
                          feq=(2*rho(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)))) + 3*rho(i,j,k)*(-1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) - 6*pyy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*(-1 + 2*u(i,j,k) - 2*v(i,j,k)) - 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(2 + 6*v(i,j,k)) + 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i-1,j+1,k,10)=feq+ (1-omega)*fneq1*p2 +(fy-fx)*p2dcssq
                          
                          !9
                          feq=(2*rho(i,j,k)*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(3 + 9*(-1 + v(i,j,k))*v(i,j,k))) - 3*rho(i,j,k)*(1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=(-3*(-2*pyy(i,j,k) + pzz(i,j,k) - 6*pyy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*(1 + 2*u(i,j,k) - 2*v(i,j,k)) - 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(-2 + 6*v(i,j,k)) + 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i+1,j-1,k,9)=feq+ (1-omega)*fneq1*p2 + (fx-fy)*p2dcssq

                          !5
                          
                          feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k))))/27.
                          fneq1=(-3*(pxx(i,j,k) + pyy(i,j,k) - 2*pzz(i,j,k) + 6*pxz(i,j,k)*u(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pxx(i,j,k)*w(i,j,k) + 3*pyy(i,j,k)*w(i,j,k)))/2.
                          f(i,j,k+1,5)=feq+ (1-omega)*fneq1*p1 + fz*p1dcssq
                          
                          !6
                          feq=(rho(i,j,k)*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 6*pxz(i,j,k)*u(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(-1 + 3*w(i,j,k))))/2.
                          f(i,j,k-1,6)=feq+ (1-omega)*fneq1*p1 - fz*p1dcssq

                          !15
                          
                          feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) - 3*pyy(i,j,k)*u(i,j,k) + 6*pzz(i,j,k)*u(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) - 3*pyy(i,j,k)*w(i,j,k) + 6*pxz(i,j,k)*(1 + 2*u(i,j,k) + 2*w(i,j,k)) + pxx(i,j,k)*(2 + 6*w(i,j,k))))/2.
                          f(i+1,j,k+1,15)=feq+ (1-omega)*fneq1*p2 + (fx+fz)*p2dcssq 
                          
                          !16
                          feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(2 - 6*w(i,j,k)) + 3*pyy(i,j,k)*w(i,j,k) - 6*pxz(i,j,k)*(-1 + 2*u(i,j,k) + 2*w(i,j,k))))/2.
                          f(i-1,j,k-1,16)=feq+ (1-omega)*fneq1*p2 - (fx+fz)*p2dcssq

                          !17
                          
                          feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*(-1 + 2*u(i,j,k) - 2*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(2 + 6*w(i,j,k))))/2.
                          f(i-1,j,k+1,17)=feq+ (1-omega)*fneq1*p2 + (fz-fx)*p2dcssq
                          
                          !18
                          feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                          fneq1=(-3*(pyy(i,j,k) - 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*(1 + 2*u(i,j,k) - 2*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(-2 + 6*w(i,j,k))))/2.
                          f(i+1,j,k-1,18)=feq+ (1-omega)*fneq1*p2 + (fx-fz)*p2dcssq

                          !11
                          
                          feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=(-3*pxx(i,j,k)*(1 + 3*v(i,j,k) + 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) - 3*pxy(i,j,k)*u(i,j,k) - 3*pxz(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*pyy(i,j,k)*w(i,j,k) + pyz(i,j,k)*(3 + 6*v(i,j,k) + 6*w(i,j,k)))
                          f(i,j+1,k+1,11)=feq+ (1-omega)*fneq1*p2 + (fy+fz)*p2dcssq
                          
                          !12
                          feq=(rho(i,j,k)*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                          fneq1=(3*pxx(i,j,k)*(-1 + 3*v(i,j,k) + 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) + 3*pxy(i,j,k)*u(i,j,k) + 3*pxz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + pyz(i,j,k)*(3 - 6*v(i,j,k) - 6*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k))
                          f(i,j-1,k-1,12)=feq+ (1-omega)*fneq1*p2 - (fy+fz)*p2dcssq

                          !13
                          
                          feq=(rho(i,j,k)*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                          fneq1=(-3*(pxx(i,j,k) - 2*pyy(i,j,k) + 6*pyz(i,j,k) - 2*pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) - 6*pxz(i,j,k)*u(i,j,k) + 3*pxx(i,j,k)*v(i,j,k) + 12*pyz(i,j,k)*v(i,j,k) - 6*pzz(i,j,k)*v(i,j,k) - 3*pxx(i,j,k)*w(i,j,k) + 6*pyy(i,j,k)*w(i,j,k) - 12*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i,j+1,k-1,13)=feq+ (1-omega)*fneq1*p2 + (fy-fz)*p2dcssq
                          
                          !14
                          feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=(3*pxx(i,j,k)*(-1 + 3*v(i,j,k) - 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) + 3*pxy(i,j,k)*u(i,j,k) - 3*pxz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + pyz(i,j,k)*(-3 + 6*v(i,j,k) - 6*w(i,j,k)) + 3*pyy(i,j,k)*w(i,j,k))
                          f(i,j-1,k+1,14)=feq+ (1-omega)*fneq1*p2 + (fz-fy)*p2dcssq

                          !19
                          feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                          fneq1=3*(pxx(i,j,k) + (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) + 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j+1,k+1,19)=feq + (1-omega)*fneq1*p3 + (fz+fy+fx)*p3dcssq
                          
                          !20
                          feq=(rho(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                          fneq1=-3*((pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(-1 + 3*u(i,j,k)) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + pxz(i,j,k)*(-3 + 6*u(i,j,k)) + 3*(2*pxy(i,j,k) + 3*pxz(i,j,k) + 2*pyz(i,j,k) + pzz(i,j,k))*v(i,j,k) + 3*(3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k) + pxx(i,j,k)*(-1 + 3*v(i,j,k) + 3*w(i,j,k)))
                          f(i-1,j-1,k-1,20)=feq+ (1-omega)*fneq1*p3 - (fz+fy+fx)*p3dcssq

                          !21
                          feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                          fneq1=3*(pxx(i,j,k) - 3*pxy(i,j,k)*(1 + 2*u(i,j,k)) + (pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) - 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j-1,k+1,21)=feq+ (1-omega)*fneq1*p3 + (fx-fy+fz)*p3dcssq
                          
                          !22
                          feq=(rho(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                          fneq1=3*(pxx(i,j,k) + 3*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k) - 3*(2*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) - 3*(pxx(i,j,k) - 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i-1,j+1,k-1,22)=feq+ (1-omega)*fneq1*p3 - (fx-fy+fz)*p3dcssq

                          !23
                          feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                          fneq1=3*(pxx(i,j,k) + 3*pxy(i,j,k) - 3*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k) - 3*(2*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k) - 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) + 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i-1,j-1,k+1,23)=feq+ (1-omega)*fneq1*p3 + (fz-fy-fx)*p3dcssq

                          !24
                          feq=(rho(i,j,k)*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                          fneq1=3*(pxx(i,j,k) - 3*pxz(i,j,k)*(1 + 2*u(i,j,k)) + (pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) - 3*(pxx(i,j,k) + 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j+1,k-1,24)=feq+ (1-omega)*fneq1*p3 + (-fz+fy+fx)*p3dcssq

                          !25
                          feq=(rho(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                          fneq1=-3*(-pxx(i,j,k) - (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j-1,k-1,25)=feq + (1-omega)*fneq1*p3 + (-fz-fy+fx)*p3dcssq

                          !26
                          feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                          fneq1=3*(pxx(i,j,k) - (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(-1 + 3*u(i,j,k)) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + pxz(i,j,k)*(-3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                          f(i-1,j+1,k+1,26)=feq + (1-omega)*fneq1*p3 + (fz+fy-fx)*p3dcssq
                      endif
                  enddo
              enddo
          enddo
        !***********************************boundary conditions no slip everywhere********************************!
            !$acc loop independent 
            do k=1,nz
                !$acc loop independent 
                do j=1,ny
                    !$acc loop independent 
                    do i=1,nx
                        if(isfluid(i,j,k).eq.0)then
                                ! 0  1   2  3   4   5   6   7    8   9   10  11   12  13   14  15   16   17   18  19  20  21  22  23  24  25  26
                            !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1,  1, -1,  1, -1, -1,  1,  1, -1/)
                            !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0,  1  -1, -1,  1, -1,  1, -1,  1/)
                            !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1,  1, -1,  1, -1,  1, -1, -1,  1/)
                            f(i+1,j-1,k-1,25)=f(i,j,k,26) !gpc 
                            f(i-1,j+1,k+1,26)=f(i,j,k,25) !hpc

                            f(i-1,j-1,k+1,23)=f(i,j,k,24) !gpc 
                            f(i+1,j+1,k-1,24)=f(i,j,k,23) !hpc

                            f(i+1,j-1,k+1,21)=f(i,j,k,22) !gpc 
                            f(i-1,j+1,k-1,22)=f(i,j,k,21) !hpc

                            f(i+1,j+1,k+1,19)=f(i,j,k,20) !gpc 
                            f(i-1,j-1,k-1,20)=f(i,j,k,19) !hpc

                            f(i+1,j,k-1,18)=f(i,j,k,17) !gpc 
                            f(i-1,j,k+1,17)=f(i,j,k,18) !hpc

                            f(i-1,j,k-1,16)=f(i,j,k,15) !gpc 
                            f(i+1,j,k+1,15)=f(i,j,k,16) !hpc

                            f(i,j-1,k+1,14)=f(i,j,k,13)!gpc 
                            f(i,j+1,k-1,13)=f(i,j,k,14)!hpc
                            
                            f(i,j-1,k-1,12)=f(i,j,k,11)!gpc 
                            f(i,j+1,k+1,11)=f(i,j,k,12)!hpc

                            f(i-1,j+1,k,10)=f(i,j,k,9)!gpc 
                            f(i+1,j-1,k,9)=f(i,j,k,10)!hpc

                            f(i-1,j-1,k,8)=f(i,j,k,7)!gpc 
                            f(i+1,j+1,k,7)=f(i,j,k,8)!hpc

                            f(i,j,k-1,6)=f(i,j,k,5)!gpc 
                            f(i,j,k+1,5)=f(i,j,k,6)!hpc 

                            f(i,j-1,k,4)=f(i,j,k,3)!gpc 
                            f(i,j+1,k,3)=f(i,j,k,4)!hpc 

                            f(i-1,j,k,2)=f(i,j,k,1)!gpc 
                            f(i+1,j,k,1)=f(i,j,k,2)!hpc 
                        endif
                    enddo
                enddo
            enddo
        
         !$acc end kernels
        
        !***********************************call other bcs:PERIODIC ************************
          if(lpbc)then      
            !periodic along x 
            !$acc kernels 
                ! 0  1   2  3   4   5   6   7    8   9   10  11   12  13   14  15   16   17   18  19  20  21  22  23  24  25  26

            !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1,  1, -1,  1, -1, -1,  1,  1, -1/)
            !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0,  1  -1, -1,  1, -1,  1, -1,  1/)
            !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1,  1, -1,  1, -1,  1, -1, -1,  1/)
            
            f(2,2:ny-1,2:nz-1,1)=f(nx-1,2:ny-1,2:nz-1,1)
            f(2,2:ny-1,2:nz-1,7)=f(nx-1,2:ny-1,2:nz-1,7)
            f(2,2:ny-1,2:nz-1,9)=f(nx-1,2:ny-1,2:nz-1,9)
            f(2,2:ny-1,2:nz-1,15)=f(nx-1,2:ny-1,2:nz-1,15)
            f(2,2:ny-1,2:nz-1,18)=f(nx-1,2:ny-1,2:nz-1,18)
            f(2,2:ny-1,2:nz-1,19)=f(nx-1,2:ny-1,2:nz-1,19)
            f(2,2:ny-1,2:nz-1,21)=f(nx-1,2:ny-1,2:nz-1,21)
            f(2,2:ny-1,2:nz-1,24)=f(nx-1,2:ny-1,2:nz-1,24)
            f(2,2:ny-1,2:nz-1,25)=f(nx-1,2:ny-1,2:nz-1,25)
            
            f(nx-1,2:ny-1,2:nz-1,2)=f(2,2:ny-1,2:nz-1,2)
            f(nx-1,2:ny-1,2:nz-1,8)=f(2,2:ny-1,2:nz-1,8)
            f(nx-1,2:ny-1,2:nz-1,10)=f(2,2:ny-1,2:nz-1,10)
            f(nx-1,2:ny-1,2:nz-1,16)=f(2,2:ny-1,2:nz-1,16)
            f(nx-1,2:ny-1,2:nz-1,17)=f(2,2:ny-1,2:nz-1,17)
            f(nx-1,2:ny-1,2:nz-1,20)=f(2,2:ny-1,2:nz-1,20)
            f(nx-1,2:ny-1,2:nz-1,22)=f(2,2:ny-1,2:nz-1,22)
            f(nx-1,2:ny-1,2:nz-1,23)=f(2,2:ny-1,2:nz-1,23)
            f(nx-1,2:ny-1,2:nz-1,26)=f(2,2:ny-1,2:nz-1,26)
            !$acc end kernels
            
            !periodic along y
            !$acc kernels 
                ! 0  1   2  3   4   5   6   7    8   9   10  11   12  13   14  15   16   17   18  19  20  21  22  23  24  25  26
                
            !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1,  1, -1,  1, -1, -1,  1,  1, -1/)
            !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0,  1  -1, -1,  1, -1,  1, -1,  1/)
            !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1,  1, -1,  1, -1,  1, -1, -1,  1/)
            
            f(2:nx-1,2,2:nz-1,3)=f(2:nx-1,ny-1,2:nz-1,3)
            f(2:nx-1,2,2:nz-1,7)=f(2:nx-1,ny-1,2:nz-1,7)
            f(2:nx-1,2,2:nz-1,10)=f(2:nx-1,ny-1,2:nz-1,10)
            f(2:nx-1,2,2:nz-1,11)=f(2:nx-1,ny-1,2:nz-1,11)
            f(2:nx-1,2,2:nz-1,13)=f(2:nx-1,ny-1,2:nz-1,13)
            f(2:nx-1,2,2:nz-1,19)=f(2:nx-1,ny-1,2:nz-1,19)
            f(2:nx-1,2,2:nz-1,22)=f(2:nx-1,ny-1,2:nz-1,22)
            f(2:nx-1,2,2:nz-1,24)=f(2:nx-1,ny-1,2:nz-1,24)
            f(2:nx-1,2,2:nz-1,26)=f(2:nx-1,ny-1,2:nz-1,26)
      
            f(2:nx-1,ny-1,2:nz-1,4)=f(2:nx-1,2,2:nz-1,4)
            f(2:nx-1,ny-1,2:nz-1,8)=f(2:nx-1,2,2:nz-1,8)
            f(2:nx-1,ny-1,2:nz-1,9)=f(2:nx-1,2,2:nz-1,9)
            f(2:nx-1,ny-1,2:nz-1,12)=f(2:nx-1,2,2:nz-1,12)
            f(2:nx-1,ny-1,2:nz-1,14)=f(2:nx-1,2,2:nz-1,14)
            f(2:nx-1,ny-1,2:nz-1,20)=f(2:nx-1,2,2:nz-1,20)
            f(2:nx-1,ny-1,2:nz-1,21)=f(2:nx-1,2,2:nz-1,21)
            f(2:nx-1,ny-1,2:nz-1,23)=f(2:nx-1,2,2:nz-1,23)
            f(2:nx-1,ny-1,2:nz-1,25)=f(2:nx-1,2,2:nz-1,25)

			      !$acc end kernels
			
          endif      
          
          if(lpbcfull)then
          
            !pbc side along x
            !$acc kernels

            
            f(1,1:ny,1:nz,1)=f(nx+1,1:ny,1:nz,1)
            f(1,1:ny,1:nz,7)=f(nx+1,1:ny,1:nz,7)
            f(1,1:ny,1:nz,9)=f(nx+1,1:ny,1:nz,9)
            f(1,1:ny,1:nz,15)=f(nx+1,1:ny,1:nz,15)
            f(1,1:ny,1:nz,18)=f(nx+1,1:ny,1:nz,18)
            f(1,1:ny,1:nz,19)=f(nx+1,1:ny,1:nz,19)
            f(1,1:ny,1:nz,21)=f(nx+1,1:ny,1:nz,21)
            f(1,1:ny,1:nz,24)=f(nx+1,1:ny,1:nz,24)
            f(1,1:ny,1:nz,25)=f(nx+1,1:ny,1:nz,25)
            
            f(nx,1:ny,1:nz,2)=f(0,1:ny,1:nz,2)
            f(nx,1:ny,1:nz,8)=f(0,1:ny,1:nz,8)
            f(nx,1:ny,1:nz,10)=f(0,1:ny,1:nz,10)
            f(nx,1:ny,1:nz,16)=f(0,1:ny,1:nz,16)
            f(nx,1:ny,1:nz,17)=f(0,1:ny,1:nz,17)
            f(nx,1:ny,1:nz,20)=f(0,1:ny,1:nz,20)
            f(nx,1:ny,1:nz,22)=f(0,1:ny,1:nz,22)
            f(nx,1:ny,1:nz,23)=f(0,1:ny,1:nz,23)
            f(nx,1:ny,1:nz,26)=f(0,1:ny,1:nz,26)
            !$acc end kernels   
            
            !pbc side along y
            !$acc kernels
            
            f(1:nx,1,1:nz,3)=f(1:nx,ny+1,1:nz,3)
            f(1:nx,1,1:nz,7)=f(1:nx,ny+1,1:nz,7)
            f(1:nx,1,1:nz,10)=f(1:nx,ny+1,1:nz,10)
            f(1:nx,1,1:nz,11)=f(1:nx,ny+1,1:nz,11)
            f(1:nx,1,1:nz,13)=f(1:nx,ny+1,1:nz,13)
            f(1:nx,1,1:nz,19)=f(1:nx,ny+1,1:nz,19)
            f(1:nx,1,1:nz,22)=f(1:nx,ny+1,1:nz,22)
            f(1:nx,1,1:nz,24)=f(1:nx,ny+1,1:nz,24)
            f(1:nx,1,1:nz,26)=f(1:nx,ny+1,1:nz,26)
            
            f(1:nx,ny,1:nz,4)=f(1:nx,0,1:nz,4)
            f(1:nx,ny,1:nz,8)=f(1:nx,0,1:nz,8)
            f(1:nx,ny,1:nz,9)=f(1:nx,0,1:nz,9)
            f(1:nx,ny,1:nz,12)=f(1:nx,0,1:nz,12)
            f(1:nx,ny,1:nz,14)=f(1:nx,0,1:nz,14)
            f(1:nx,ny,1:nz,20)=f(1:nx,0,1:nz,20)
            f(1:nx,ny,1:nz,21)=f(1:nx,0,1:nz,21)
            f(1:nx,ny,1:nz,23)=f(1:nx,0,1:nz,23)
            f(1:nx,ny,1:nz,25)=f(1:nx,0,1:nz,25)
            
            !$acc end kernels  
            
            !pbc side along z
            !$acc kernels
            
                ! 0  1   2  3   4   5   6   7    8   9   10  11   12  13   14  15   16   17   18  19  20  21  22  23  24  25  26

            !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1,  1, -1,  1, -1, -1,  1,  1, -1/)
            !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0,  1  -1, -1,  1, -1,  1, -1,  1/)
            !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1,  1, -1,  1, -1,  1, -1, -1,  1/)
            
            f(1:nx,1:ny,1,5)=f(1:nx,1:ny,nz+1,5)
            f(1:nx,1:ny,1,11)=f(1:nx,1:ny,nz+1,11)
            f(1:nx,1:ny,1,14)=f(1:nx,1:ny,nz+1,14)
            f(1:nx,1:ny,1,15)=f(1:nx,1:ny,nz+1,15)
            f(1:nx,1:ny,1,17)=f(1:nx,1:ny,nz+1,17)
            f(1:nx,1:ny,1,19)=f(1:nx,1:ny,nz+1,19)
            f(1:nx,1:ny,1,21)=f(1:nx,1:ny,nz+1,21)
            f(1:nx,1:ny,1,23)=f(1:nx,1:ny,nz+1,23)
            f(1:nx,1:ny,1,26)=f(1:nx,1:ny,nz+1,26)
            
            f(1:nx,1:ny,nz,6)=f(1:nx,1:ny,0,6)
            f(1:nx,1:ny,nz,12)=f(1:nx,1:ny,0,12)
            f(1:nx,1:ny,nz,13)=f(1:nx,1:ny,0,13)
            f(1:nx,1:ny,nz,16)=f(1:nx,1:ny,0,16)
            f(1:nx,1:ny,nz,18)=f(1:nx,1:ny,0,18)
            f(1:nx,1:ny,nz,20)=f(1:nx,1:ny,0,20)
            f(1:nx,1:ny,nz,22)=f(1:nx,1:ny,0,22)
            f(1:nx,1:ny,nz,24)=f(1:nx,1:ny,0,24)
            f(1:nx,1:ny,nz,25)=f(1:nx,1:ny,0,25)
            
            !$acc end kernels  
      
            
            !pbc edge along x
            !$acc kernels
            
                ! 0  1   2  3   4   5   6   7    8   9   10  11   12  13   14  15   16   17   18  19  20  21  22  23  24  25  26

            !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1,  1, -1,  1, -1, -1,  1,  1, -1/)
            !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0,  1  -1, -1,  1, -1,  1, -1,  1/)
            !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1,  1, -1,  1, -1,  1, -1, -1,  1/)
            
            !1:nx,1,1,1 <= 1:nx,ny+1,nz+1   *  +1  +1
            f(1:nx,1,1,11)=f(1:nx,ny+1,nz+1,11)
            f(1:nx,1,1,19)=f(1:nx,ny+1,nz+1,19)
            f(1:nx,1,1,26)=f(1:nx,ny+1,nz+1,26)
            
            !1:nx,ny,nz <= 1:nx,0,0         *  -1  -1
            f(1:nx,ny,nz,12)=f(1:nx,0,0,12)
            f(1:nx,ny,nz,20)=f(1:nx,0,0,20)
            f(1:nx,ny,nz,25)=f(1:nx,0,0,25)
            
            !1:nx,1,nz <= 1:nx,ny+1,0       *  +1  -1
            f(1:nx,1,nz,13)=f(1:nx,ny+1,0,13)
            f(1:nx,1,nz,22)=f(1:nx,ny+1,0,22)
            f(1:nx,1,nz,24)=f(1:nx,ny+1,0,24)
            
            !1:nx,ny,1 <= 1:nx,0,nz+1       *  -1  +1
            f(1:nx,ny,1,14)=f(1:nx,0,nz+1,14)
            f(1:nx,ny,1,21)=f(1:nx,0,nz+1,21)
            f(1:nx,ny,1,23)=f(1:nx,0,nz+1,23)
            
            !$acc end kernels   
            
            !pbc edge along y
            !$acc kernels
            
                ! 0  1   2  3   4   5   6   7    8   9   10  11   12  13   14  15   16   17   18  19  20  21  22  23  24  25  26

            !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1,  1, -1,  1, -1, -1,  1,  1, -1/)
            !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0,  1  -1, -1,  1, -1,  1, -1,  1/)
            !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1,  1, -1,  1, -1,  1, -1, -1,  1/)
            
            !1,1:ny,1 <= nx+1,1:ny,nz+1   +1  *  +1
            f(1,1:ny,1,15)=f(nx+1,1:ny,nz+1,15)
            f(1,1:ny,1,19)=f(nx+1,1:ny,nz+1,19)
            f(1,1:ny,1,21)=f(nx+1,1:ny,nz+1,21)
            
            !nx,1:ny,nz <= 0,1:ny,0   -1  *  -1
            f(nx,1:ny,nz,16)=f(0,1:ny,0,16)
            f(nx,1:ny,nz,20)=f(0,1:ny,0,20)
            f(nx,1:ny,nz,22)=f(0,1:ny,0,22)
            
            !1,1:ny,nz <= nx+1,1:ny,0   +1  *  -1
            f(1,1:ny,nz,18)=f(nx+1,1:ny,0,18)
            f(1,1:ny,nz,24)=f(nx+1,1:ny,0,24)
            f(1,1:ny,nz,25)=f(nx+1,1:ny,0,25)
            
            !nx,1:ny,1 <= 0,1:ny,nz+1   -1  *  +1
            f(nx,1:ny,1,17)=f(0,1:ny,nz+1,17)
            f(nx,1:ny,1,23)=f(0,1:ny,nz+1,23)
            f(nx,1:ny,1,26)=f(0,1:ny,nz+1,26)
            
            !$acc end kernels  
            
            !pbc edge along z
            !$acc kernels 
            
                ! 0  1   2  3   4   5   6   7    8   9   10  11   12  13   14  15   16   17   18  19  20  21  22  23  24  25  26

            !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1,  1, -1,  1, -1, -1,  1,  1, -1/)
            !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0,  1  -1, -1,  1, -1,  1, -1,  1/)
            !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1,  1, -1,  1, -1,  1, -1, -1,  1/)
            
            !1,1,1:nz <= nx+1,ny+1,1:nz   +1  +1  *
            f(1,1,1:nz,7)=f(nx+1,ny+1,1:nz,7)
            f(1,1,1:nz,19)=f(nx+1,ny+1,1:nz,19)
            f(1,1,1:nz,24)=f(nx+1,ny+1,1:nz,24)
            
            !nx,ny,1:nz <= 0,0,1:nz   -1  -1  *
            f(nx,ny,1:nz,8)=f(0,0,1:nz,8)
            f(nx,ny,1:nz,20)=f(0,0,1:nz,20)
            f(nx,ny,1:nz,23)=f(0,0,1:nz,23)
            
            !1,ny,1:nz <= nx+1,0,1:nz   +1  -1  *
            f(1,ny,1:nz,9)=f(nx+1,0,1:nz,9)
            f(1,ny,1:nz,21)=f(nx+1,0,1:nz,21)
            f(1,ny,1:nz,25)=f(nx+1,0,1:nz,25)
            
            !nx,1,1:nz <= 0,ny+1,1:nz   -1  +1  *
            f(nx,1,1:nz,10)=f(0,ny+1,1:nz,10)
            f(nx,1,1:nz,22)=f(0,ny+1,1:nz,22)
            f(nx,1,1:nz,26)=f(0,ny+1,1:nz,26)
            
            !$acc end kernels  
            
            !pbc corner
            !$acc kernels 
            
                ! 0  1   2  3   4   5   6   7    8   9   10  11   12  13   14  15   16   17   18  19  20  21  22  23  24  25  26

            !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1,  1, -1,  1, -1, -1,  1,  1, -1/)
            !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0,  1  -1, -1,  1, -1,  1, -1,  1/)
            !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1,  1, -1,  1, -1,  1, -1, -1,  1/)
            
            !1,1,1, <= nx+1,ny+1,nz+1   +1  +1  +1
            f(1,1,1,19)=f(nx+1,ny+1,nz+1,19)
            
            !nx,ny,nz <= 0,0,0   -1  -1  -1
            f(nx,ny,nz,20)=f(0,0,0,20)
            
            !1,ny,1 <= nx+1,0,nz+1  +1  -1  +1
            f(1,ny,1,21)=f(nx+1,0,nz+1,21)
            
            !nx,1,nz <= 0,ny+1,0    -1  +1  -1
            f(nx,1,nz,22)=f(0,ny+1,0,22)
            
            !nx,ny,1 <= 0,0,nz+1    -1  -1  +1
            f(nx,ny,1,23)=f(0,0,nz+1,23)
            
            !1,1,nz <= nx+1,ny+1,0  +1  +1  -1
            f(1,1,nz,24)=f(nx+1,ny+1,0,24)
            
            !1,ny,nz <= nx+1,0,0  +1  -1  -1
            f(1,ny,nz,25)=f(nx+1,0,0,25)
            
            !nx,1,1 <= 0,ny+1,nz+1  -1  +1  +1
            f(nx,1,1,26)=f( 0,ny+1,nz+1,26)
            
            !$acc end kernels  
          
          endif   
            
    enddo 
    !$acc end data
    call cpu_time(ts2)
    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
    write(6,*) 'glups: ',  real(nx)*real(ny)*real(nz)*real(nsteps)/1.0e9/(ts2-ts1)
    
    call get_memory_gpu(mymemory,totmemory)
    call print_memory_registration_gpu(6,'DEVICE memory occupied at the end', &
     'total DEVICE memory',mymemory,totmemory)
    
end program
