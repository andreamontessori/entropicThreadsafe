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

    nlinks=18 !pari!
    tau=1.0_db
    cssq=1.0_db/3.0_db
    cssq2=1.0_db/cssq**2
    cssq3=cssq**3
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
        nsteps=1000
        stamp=1000
        fx=1.0_db*10.0**(-5)
        fy=0.0_db*10.0**(-5)
        fz=0.0_db*10.0**(-5)
        lprint=.true.
        lvtk=.true.
        lasync=.false.
        lpbc=.true.
        
        allocate(f(0:nx+1,0:ny+1,0:nz+1,0:nlinks))
        allocate(rho(1:nx,1:ny,1:nz),u(1:nx,1:ny,1:nz),v(1:nx,1:ny,1:nz),w(1:nx,1:ny,1:nz))
        allocate(pxx(1:nx,1:ny,1:nz),pxy(1:nx,1:ny,1:nz),pxz(1:nx,1:ny,1:nz),pyy(1:nx,1:ny,1:nz))
        allocate(pyz(1:nx,1:ny,1:nz),pzz(1:nx,1:ny,1:nz))
        allocate(isfluid(1:nx,1:ny,1:nz)) !,omega_2d(1:nx,1:ny)) 
        if(lprint)then
          allocate(rhoprint(1:nx,1:ny,1:nz))
          allocate(velprint(1:3,1:nx,1:ny,1:nz))
          rhoprint(1:nx,1:ny,1:nz)=0.0
          velprint(1:3,1:nx,1:ny,1:nz)=0.0
        endif
        !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
        !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
        !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)

        p0=(1.0_db/3.0_db)
        p1=(1.0_db/18.0_db)
        p2=(1.0_db/36.0_db)
        p1dcssq=p1/cssq
        p2dcssq=p2/cssq
        omega=1.0_db/tau
    !*****************************************geometry************************
        isfluid=1
        isfluid(1,:,:)=0 !left
        isfluid(nx,:,:)=0 !right
        isfluid(:,1,:)=0 !front 
        isfluid(:,ny,:)=0 !rear
        isfluid(:,:,1)=0 !bottom
        isfluid(:,:,nz)=0 !top
    !****************************************hermite projection vars**********
        pi2cssq0=p0/(2.0_db*cssq**2)
        pi2cssq1=p1/(2.0_db*cssq**2)
        pi2cssq2=p2/(2.0_db*cssq**2)

        qxx=1.0_db-cssq
        qyy=1.0_db-cssq
        qzz=1.0_db-cssq
        qxy_7_8=1.0_db
        qxy_9_10=-1.0_db
        qxz_15_16=1.0_db
        qxz_17_18=-1.0_db
        qyz_11_12=1.0_db
        qyz_13_14=-1.0_db

    !*************************************initial conditions ************************    
        u=0.0_db
        v=0.0_db
        w=0.0_db
        rho=1.0_db  !tot dens
        !do ll=0,nlinks
        f(1:nx,1:ny,1:nz,0)=rho(1:nx,1:ny,1:nz)*p0
        f(1:nx,1:ny,1:nz,1)=rho(1:nx,1:ny,1:nz)*p1
        f(1:nx,1:ny,1:nz,2)=rho(1:nx,1:ny,1:nz)*p1
        f(1:nx,1:ny,1:nz,3)=rho(1:nx,1:ny,1:nz)*p1
        f(1:nx,1:ny,1:nz,4)=rho(1:nx,1:ny,1:nz)*p1
        f(1:nx,1:ny,1:nz,5)=rho(1:nx,1:ny,1:nz)*p1
        f(1:nx,1:ny,1:nz,6)=rho(1:nx,1:ny,1:nz)*p1
        f(1:nx,1:ny,1:nz,7)=rho(1:nx,1:ny,1:nz)*p2
        f(1:nx,1:ny,1:nz,8)=rho(1:nx,1:ny,1:nz)*p2
        f(1:nx,1:ny,1:nz,9)=rho(1:nx,1:ny,1:nz)*p2
        f(1:nx,1:ny,1:nz,10)=rho(1:nx,1:ny,1:nz)*p2
        f(1:nx,1:ny,1:nz,11)=rho(1:nx,1:ny,1:nz)*p2
        f(1:nx,1:ny,1:nz,12)=rho(1:nx,1:ny,1:nz)*p2
        f(1:nx,1:ny,1:nz,13)=rho(1:nx,1:ny,1:nz)*p2
        f(1:nx,1:ny,1:nz,14)=rho(1:nx,1:ny,1:nz)*p2
        f(1:nx,1:ny,1:nz,15)=rho(1:nx,1:ny,1:nz)*p2
        f(1:nx,1:ny,1:nz,16)=rho(1:nx,1:ny,1:nz)*p2
        f(1:nx,1:ny,1:nz,17)=rho(1:nx,1:ny,1:nz)*p2
        f(1:nx,1:ny,1:nz,18)=rho(1:nx,1:ny,1:nz)*p2
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
    !$acc data copy(f,isfluid,p0,p1,p2,&
             !$acc& pxx,pyy,pzz,pxy,pxz,pyz,rho,u,v,w,rhoprint,velprint) async(1)
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
    
    !$acc wait(1)
    if(lprint)then
      !$acc kernels present(rhoprint,velprint,rho,u,v,w) async(1)
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
      !$acc wait(1)
      if(lasync)then
        !$acc update host(rhoprint,velprint) async(2)
        continue
      else
        !$acc update host(rhoprint,velprint) async(2)
        !$acc wait(2)
        if(lvtk)then
          call print_vtk_sync(iframe)
        else
          call print_raw_sync(iframe)
        endif
      endif
    endif
      
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !***********************************moments collision bbck + forcing************************ 
        !$acc kernels async(1)
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
                            +f(i,j,k,18)
                            

                        u(i,j,k) = ((f(i,j,k,1)+f(i,j,k,7)+f(i,j,k,9)+f(i,j,k,15)+f(i,j,k,18)) &
                             -(f(i,j,k,2)+f(i,j,k,8)+f(i,j,k,10)+f(i,j,k,16)+f(i,j,k,17)))/rho(i,j,k)
                        
                        v(i,j,k) = ((f(i,j,k,3)+f(i,j,k,7)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,13)) &
                            -(f(i,j,k,4)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,12)+f(i,j,k,14)))/rho(i,j,k)

                        w(i,j,k) = ((f(i,j,k,5)+f(i,j,k,11)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,17)) &
                            -(f(i,j,k,6)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,16)+f(i,j,k,18)))/rho(i,j,k)
                        
                        uu=(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))*cssq
                        !1-2
                        udotc=u(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p1*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=f(i,j,k,1)-feq
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        feq=p1*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        fneq1=f(i,j,k,2)-feq
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        !3-4
                        udotc=v(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p1*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=f(i,j,k,3)-feq
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        feq=p1*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        fneq1=f(i,j,k,4)-feq
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        !5-6
                        udotc=w(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p1*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=f(i,j,k,5)-feq
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        feq=p1*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        fneq1=f(i,j,k,6)-feq
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        !7-8
                        udotc=u(i,j,k)+v(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=f(i,j,k,7)-feq
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pxy(i,j,k)=pxy(i,j,k) + fneq1
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        fneq1=f(i,j,k,8)-feq
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pxy(i,j,k)=pxy(i,j,k)+fneq1
                        !10-9
                        udotc=-u(i,j,k)+v(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=f(i,j,k,10)-feq
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pxy(i,j,k)=pxy(i,j,k)-fneq1
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        fneq1=f(i,j,k,9)-feq
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pxy(i,j,k)=pxy(i,j,k)-fneq1
                        !11-12
                        udotc=v(i,j,k)+w(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=f(i,j,k,11)-feq
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pyz(i,j,k)=pyz(i,j,k)+fneq1
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        fneq1=f(i,j,k,12)-feq
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pyz(i,j,k)=pyz(i,j,k)+fneq1
                        !13-14
                        udotc=v(i,j,k)-w(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=f(i,j,k,13) - feq
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pyz(i,j,k)=pyz(i,j,k)-fneq1
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        fneq1=f(i,j,k,14) - feq
                        pyy(i,j,k)=pyy(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pyz(i,j,k)=pyz(i,j,k)-fneq1
                        !15-16
                        udotc=u(i,j,k)+w(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=f(i,j,k,15)-feq
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pxz(i,j,k)=pxz(i,j,k)+fneq1
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        fneq1=f(i,j,k,16)-feq
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pxz(i,j,k)=pxz(i,j,k)+fneq1
                        !17-18
                        udotc=-u(i,j,k)+w(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=f(i,j,k,17)-feq
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pxz(i,j,k)=pxz(i,j,k)-fneq1
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        fneq1=f(i,j,k,18)-feq
                        pxx(i,j,k)=pxx(i,j,k)+fneq1
                        pzz(i,j,k)=pzz(i,j,k)+fneq1
                        pxz(i,j,k)=pxz(i,j,k)-fneq1
                    endif
                enddo
            enddo
        enddo
        !$acc end kernels
        !***********************************PRINT************************
        if(mod(step,stamp).eq.0)write(6,'(a,i8)')'step : ',step
        if(lprint)then
          if(mod(step,stamp).eq.0)then
            iframe=iframe+1
            !$acc wait(1)
            !$acc kernels present(rhoprint,velprint,rho,u,v,w) async(1)
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
           !$acc wait(1)
           if(lasync)then
              call close_print_async
              !$acc update host(rhoprint,velprint) async(2)
           else
              !$acc update host(rhoprint,velprint) async(2)
              !$acc wait(2)
              if(lvtk)then
                call print_vtk_sync(iframe)
              else
                call print_raw_sync(iframe)
              endif
            endif
          endif
          if(mod(step-stamp/4,stamp).eq.0 .and. lasync)then
            !write(6,*)'ciao 2',step,iframe
            !$acc wait(2)  
            if(lvtk)then
              call print_vtk_async(iframe)
            else
              call print_raw_async(iframe)
            endif
          endif
        endif
        
        !***********************************collision + no slip + forcing: fused implementation*********
        !$acc kernels async(1)
        !$acc loop collapse(3) private(feq,uu,temp,udotc)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(isfluid(i,j,k).eq.1)then
                        uu=(u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k))*cssq
                        a2xxy=2.0_db*u(i,j,k)*pxy(i,j,k)+pxx(i,j,k)*v(i,j,k)
                        a2xxz=2.0_db*u(i,j,k)*pxz(i,j,k)+pxx(i,j,k)*w(i,j,k)
                        a2xyy=2.0_db*v(i,j,k)*pxy(i,j,k)+pyy(i,j,k)*u(i,j,k)
                        a2xzz=2.0_db*w(i,j,k)*pxz(i,j,k)+pzz(i,j,k)*u(i,j,k)
                        a2yzz=2.0_db*w(i,j,k)*pyz(i,j,k)+v(i,j,k)*pzz(i,j,k)
                        a2yyz=2.0_db*v(i,j,k)*pyz(i,j,k)+w(i,j,k)*pyy(i,j,k)
                        a2xyz=u(i,j,k)pyz(i,j,k) + v(i,j,k)*pxz(i,j,k) + w(i,j,k)*pxy(i,j,k)
                        !0
                        feq=p0*rho(i,j,k)*(1.0_db-0.5_db*uu*cssq2)
                        f(i,j,k,0)=feq + (1.0_db-omega)*pi2cssq0*(-cssq*(pyy(i,j,k)+pxx(i,j,k)+pzz(i,j,k)))
                        
                        !1
                        udotc=u(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p1*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=(1.0_db-omega)*pi2cssq1*(qxx*pxx(i,j,k)-cssq*(pyy(i,j,k)+pzz(i,j,k)))
                        f(i+1,j,k,1)=feq + fneq1 + fx*p1dcssq
                        
                        !2
                        feq=p1*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        f(i-1,j,k,2)=feq + fneq1 - fx*p1dcssq
                        
                        !3
                        udotc=v(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p1*rho(i,j,k)*(1.0_db  + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=(1.0_db-omega)*pi2cssq1*(qyy*pyy(i,j,k)-cssq*(pxx(i,j,k)+pzz(i,j,k)))
                        f(i,j+1,k,3)=feq+fneq1 + fy*p1dcssq
                        
                        !4
                        feq=p1*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        f(i,j-1,k,4)=feq+fneq1 - fy*p1dcssq
                        
                        !7
                        udotc=u(i,j,k)+v(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_7_8*pxy(i,j,k))
                        f(i+1,j+1,k,7)=feq + fneq1 + (fx+fy)*p2dcssq 
                        
                        !8
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        f(i-1,j-1,k,8)=feq + fneq1 - (fx+fy)*p2dcssq
                        
                        !10
                        udotc=-u(i,j,k)+v(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qyy*pyy(i,j,k)-cssq*pzz(i,j,k)+2.0_db*qxy_9_10*pxy(i,j,k))
                        f(i-1,j+1,k,10)=feq+fneq1 +(fy-fx)*p2dcssq
                        
                        !9
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        f(i+1,j-1,k,9)=feq+fneq1 + (fx-fy)*p2dcssq

                        !5
                        udotc=w(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p1*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=(1.0_db-omega)*pi2cssq1*(qzz*pzz(i,j,k)-cssq*(pxx(i,j,k)+pyy(i,j,k)))
                        f(i,j,k+1,5)=feq+fneq1 + fz*p1dcssq
                        
                        !6
                        feq=p1*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        f(i,j,k-1,6)=feq+fneq1 - fz*p1dcssq

                        !15
                        udotc=u(i,j,k)+w(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_15_16*pxz(i,j,k))
                        f(i+1,j,k+1,15)=feq+fneq1 + (fx+fz)*p2dcssq 
                        
                        !16
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        f(i-1,j,k-1,16)=feq+fneq1 - (fx+fz)*p2dcssq

                        !17
                        udotc=-u(i,j,k)+w(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=(1.0_db-omega)*pi2cssq2*(qxx*pxx(i,j,k)+qzz*pzz(i,j,k)-cssq*pyy(i,j,k)+2.0_db*qxz_17_18*pxz(i,j,k))
                        f(i-1,j,k+1,17)=feq+fneq1 +(fz-fx)*p2dcssq
                        
                        !18
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 - temp)
                        f(i+1,j,k-1,18)=feq+fneq1 + (fx-fz)*p2dcssq

                        !11
                        udotc=v(i,j,k)+w(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_11_12*pyz(i,j,k))
                        f(i,j+1,k+1,11)=feq+fneq1+(fy+fz)*p2dcssq
                        
                        !12
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2- temp)
                        f(i,j-1,k-1,12)=feq+fneq1 - (fy+fz)*p2dcssq

                        !13
                        udotc=v(i,j,k)-w(i,j,k)
                        temp= udotc*3.0_db + udotc*(udotc*udotc-3.0_db*uu)/(6.0_db*cssq3)
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 + temp)
                        fneq1=(1.0_db-omega)*pi2cssq2*(qyy*pyy(i,j,k)+qzz*pzz(i,j,k)-cssq*pxx(i,j,k)+2.0_db*qyz_13_14*pyz(i,j,k))
                        f(i,j+1,k-1,13)=feq+fneq1 + (fy-fz)*p2dcssq
                        
                        !14
                        feq=p2*rho(i,j,k)*(1.0_db + 0.5_db*(udotc*udotc-uu)*cssq2 -temp)
                        f(i,j-1,k+1,14)=feq+fneq1 + (fz-fy)*p2dcssq
                    endif
                enddo
            enddo
        enddo
        !********************************boundary conditions no slip everywhere********************************!
            !$acc loop independent 
            do k=1,nz
                !$acc loop independent 
                do j=1,ny
                    !$acc loop independent 
                    do i=1,nx
                        if(isfluid(i,j,k).eq.0)then
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
        
        !******************************************call other bcs:PERIODIC ************************
          if(lpbc)then      
            !periodic along x 
            !$acc kernels async(1)
            !$acc loop independent 
            do k=2,nz-1
                !$acc loop independent 
                do j=2,ny-1
                  if(j>2 .and. j<ny-1 .and. k>2 .and. k<nz-1)then
                    f(2,j,k,1)=f(nx,j,k,1)
                    f(2,j,k,7)=f(nx,j,k,7)
                    f(2,j,k,9)=f(nx,j,k,9)
                    f(2,j,k,15)=f(nx,j,k,15)
                    f(2,j,k,18)=f(nx,j,k,18)
                    f(nx-1,j,k,2)=f(1,j,k,2)
                    f(nx-1,j,k,8)=f(1,j,k,8)
                    f(nx-1,j,k,10)=f(1,j,k,10)
                    f(nx-1,j,k,16)=f(1,j,k,16)
                    f(nx-1,j,k,17)=f(1,j,k,17)
                  else
                    if(j==2)then
                      if(k==2)then
                        f(2,j,k,1)=f(nx,j,k,1)
                        f(2,j,k,9)=f(nx,j,k,9)
                        f(2,j,k,18)=f(nx,j,k,18)
                        f(nx-1,j,k,2)=f(1,j,k,2)
                        f(nx-1,j,k,8)=f(1,j,k,8)
                        f(nx-1,j,k,16)=f(1,j,k,16)
                      elseif(k==nz-1)then
                        f(2,j,k,1)=f(nx,j,k,1)
                        f(2,j,k,9)=f(nx,j,k,9)
                        f(2,j,k,15)=f(nx,j,k,15)
                        f(nx-1,j,k,2)=f(1,j,k,2)
                        f(nx-1,j,k,8)=f(1,j,k,8)
                        f(nx-1,j,k,17)=f(1,j,k,17)
                      else
                        f(2,j,k,1)=f(nx,j,k,1)
                        f(2,j,k,9)=f(nx,j,k,9)
                        f(2,j,k,15)=f(nx,j,k,15)
                        f(2,j,k,18)=f(nx,j,k,18)
                        f(nx-1,j,k,2)=f(1,j,k,2)
                        f(nx-1,j,k,8)=f(1,j,k,8)
                        f(nx-1,j,k,16)=f(1,j,k,16)
                        f(nx-1,j,k,17)=f(1,j,k,17)
                      endif
                    elseif(j==ny-1)then
                      if(k==2)then
                        f(2,j,k,1)=f(nx,j,k,1)
                        f(2,j,k,7)=f(nx,j,k,7)
                        f(2,j,k,18)=f(nx,j,k,18)
                        f(nx-1,j,k,2)=f(1,j,k,2)
                        f(nx-1,j,k,10)=f(1,j,k,10)
                        f(nx-1,j,k,16)=f(1,j,k,16)
                      elseif(k==nz-1)then
                        f(2,j,k,1)=f(nx,j,k,1)
                        f(2,j,k,7)=f(nx,j,k,7)
                        f(2,j,k,15)=f(nx,j,k,15)
                        f(nx-1,j,k,2)=f(1,j,k,2)
                        f(nx-1,j,k,10)=f(1,j,k,10)
                        f(nx-1,j,k,17)=f(1,j,k,17)
                      else
                        f(2,j,k,1)=f(nx,j,k,1)
                        f(2,j,k,7)=f(nx,j,k,7)
                        f(2,j,k,15)=f(nx,j,k,15)
                        f(2,j,k,18)=f(nx,j,k,18)
                        f(nx-1,j,k,2)=f(1,j,k,2)
                        f(nx-1,j,k,10)=f(1,j,k,10)
                        f(nx-1,j,k,16)=f(1,j,k,16)
                        f(nx-1,j,k,17)=f(1,j,k,17)
                      endif
                    endif
                  endif 
                enddo
			      enddo
            !$acc end kernels
            
            !periodic along y
            !$acc kernels async(1)
            !$acc loop independent 
            do k=2,nz-1
              !$acc loop independent 
              do i=2,nx-1
                if(i>2 .and. i<nx-1 .and. k>2 .and. k<nz-1)then
                  f(i,2,k,3)=f(i,ny,k,3)
                  f(i,2,k,7)=f(i,ny,k,7)
                  f(i,2,k,10)=f(i,ny,k,10)
                  f(i,2,k,11)=f(i,ny,k,11)
                  f(i,2,k,13)=f(i,ny,k,13)
                  f(i,ny-1,k,4)=f(i,1,k,4)
                  f(i,ny-1,k,8)=f(i,1,k,8)
                  f(i,ny-1,k,9)=f(i,1,k,9)
                  f(i,ny-1,k,12)=f(i,1,k,12)
                  f(i,ny-1,k,14)=f(i,1,k,14)
                else
                  if(i==2)then
                    if(k==2)then
                      f(i,2,k,3)=f(i,ny,k,3)
                      f(i,2,k,10)=f(i,ny,k,10)
                      f(i,2,k,13)=f(i,ny,k,13)
                      f(i,ny-1,k,4)=f(i,1,k,4)
                      f(i,ny-1,k,8)=f(i,1,k,8)
                      f(i,ny-1,k,12)=f(i,1,k,12)
                    elseif(k==nz-1)then
                      f(i,2,k,3)=f(i,ny,k,3)
                      f(i,2,k,10)=f(i,ny,k,10)
                      f(i,2,k,11)=f(i,ny,k,11)
                      f(i,ny-1,k,4)=f(i,1,k,4)
                      f(i,ny-1,k,8)=f(i,1,k,8)
                      f(i,ny-1,k,14)=f(i,1,k,14)
                    else
                      f(i,2,k,3)=f(i,ny,k,3)
                      f(i,2,k,10)=f(i,ny,k,10)
                      f(i,2,k,11)=f(i,ny,k,11)
                      f(i,2,k,13)=f(i,ny,k,13)
                      f(i,ny-1,k,4)=f(i,1,k,4)
                      f(i,ny-1,k,8)=f(i,1,k,8)
                      f(i,ny-1,k,12)=f(i,1,k,12)
                      f(i,ny-1,k,14)=f(i,1,k,14)
                    endif
                  elseif(i==nx-1)then
                    if(k==2)then
                      f(i,2,k,3)=f(i,ny,k,3)
                      f(i,2,k,7)=f(i,ny,k,7)
                      f(i,2,k,13)=f(i,ny,k,13)
                      f(i,ny-1,k,4)=f(i,1,k,4)
                      f(i,ny-1,k,9)=f(i,1,k,9)
                      f(i,ny-1,k,12)=f(i,1,k,12)
                    elseif(k==nz-1)then 
                      f(i,2,k,3)=f(i,ny,k,3)
                      f(i,2,k,7)=f(i,ny,k,7)
                      f(i,2,k,11)=f(i,ny,k,11)
                      f(i,ny-1,k,4)=f(i,1,k,4)
                      f(i,ny-1,k,9)=f(i,1,k,9)
                      f(i,ny-1,k,14)=f(i,1,k,14)
                    else
                      f(i,2,k,3)=f(i,ny,k,3)
                      f(i,2,k,7)=f(i,ny,k,7)
                      f(i,2,k,11)=f(i,ny,k,11)
                      f(i,2,k,13)=f(i,ny,k,13)
                      f(i,ny-1,k,4)=f(i,1,k,4)
                      f(i,ny-1,k,9)=f(i,1,k,9)
                      f(i,ny-1,k,12)=f(i,1,k,12)
                      f(i,ny-1,k,14)=f(i,1,k,14)
                    endif
                  endif
                endif
			        enddo
			      enddo
			      !$acc end kernels
			
          endif
            ! !$acc kernels async(1)
            ! !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
            ! !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
            ! !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
            ! ! f(2,:,:,1)=f(nx,:,:,1)
            ! ! f(2,:,:,7)=f(nx,:,:,7)
            ! ! f(2,:,:,9)=f(nx,:,:,9)
            ! ! f(2,:,:,15)=f(nx,:,:,15)
            ! ! f(2,:,:,18)=f(nx,:,:,18)
            ! ! !x=nx 
            ! ! f(nx-1,:,:,2)=f(1,:,:,2)
            ! ! f(nx-1,:,:,8)=f(1,:,:,8)
            ! ! f(nx-1,:,:,10)=f(1,:,:,10)
            ! ! f(nx-1,:,:,16)=f(1,:,:,16)
            ! ! f(nx-1,:,:,17)=f(1,:,:,17)

            ! ! !y=1
            ! ! f3(:,2,:)=f3(:,ny,:)
            ! ! f7(:,2,:)=f7(:,ny,:)
            ! ! f10(:,2,:)=f10(:,ny,:)
            ! ! f11(:,2,:)=f11(:,ny,:)
            ! ! f13(:,2,:)=f13(:,ny,:)
        
            ! ! !y=ny
            ! ! f4(:,ny-1,:)=f4(:,1,:)
            ! ! f8(:,ny-1,:)=f8(:,1,:)
            ! ! f9(:,ny-1,:)=f9(:,1,:)
            ! ! f12(:,ny-1,:)=f12(:,1,:)
            ! ! f14(:,ny-1,:)=f14(:,1,:)
            ! !$acc end kernels
       
        
    enddo 
    !$acc wait
    if(lasync)then
      if(lvtk)then
        call print_vtk_sync(iframe)
      else
        call print_raw_sync(iframe)
      endif
    endif
    !$acc update host(rho,u,v,w)
    !$acc end data
    call cpu_time(ts2)
    write(6,*) 'u=',u(nx/2,ny/2,nz/2),'v=',v(nx/2,ny/2,nz/2),'w=',w(nx/2,ny/2,nz/2),'rho=',rho(nx/2,ny/2,nz/2)
    write(6,*) 'u=',u(nx/2,ny/2,1),'v=',v(nx/2,ny/2,1),'w=',w(nx/2,ny/2,1),'rho=',rho(nx/2,ny/2,1)
    write(6,*) 'u=',u(2,ny/2,nz/2),'v=',v(2,ny/2,nz/2),'w=',w(2,ny/2,nz/2),'rho=',rho(2,ny/2,nz/2)
    write(6,*) 'u=',u(2,ny/2,nz/2),'v=',v(2,ny/2,nz/2),'w=',w(2,ny/2,nz/2),'rho=',rho(2,ny/2,nz/2)
    write(6,*) 'u=',u(12,ny/2,nz/2),'v=',v(12,ny/2,nz/2),'w=',w(12,ny/2,nz/2),'rho=',rho(12,ny/2,nz/2)
    write(6,*) 'u=',u(24,ny/2,nz/2),'v=',v(24,ny/2,nz/2),'w=',w(24,ny/2,nz/2),'rho=',rho(24,ny/2,nz/2)
    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
    write(6,*) 'glups: ',  real(nx)*real(ny)*real(nz)*real(nsteps)/1.0e9/(ts2-ts1)
    
    call get_memory_gpu(mymemory,totmemory)
    call print_memory_registration_gpu(6,'DEVICE memory occupied at the end', &
     'total DEVICE memory',mymemory,totmemory)

  contains
  !$if _OPENACC  
    subroutine printDeviceProperties(ngpus,dev_Num,dev_Type,iu)
    
    
    use openacc
    
    integer :: ngpus,dev_Num
    integer(acc_device_kind) :: dev_Type
  
    integer,intent(in) :: iu 
    integer :: tot_mem,shared_mem
    character(len=255) :: myname,myvendor,mydriver
    
    call acc_get_property_string(dev_num,dev_type,acc_property_name,myname)
    tot_mem = acc_get_property(dev_num,dev_type,acc_property_memory)
    call acc_get_property_string(dev_num,dev_type,acc_property_vendor,myvendor)
    call acc_get_property_string(dev_num,dev_type,acc_property_driver,mydriver)
    
    write(iu,907)"                                                                               "
    write(iu,907)"*****************************GPU FEATURE MONITOR*******************************"
    write(iu,907)"                                                                               "
    
    write (iu,900) "Device Number: "      ,ngpus
    write (iu,901) "Device Name: "        ,trim(myname)
    write (iu,903) "Total Global Memory: ",real(tot_mem)/1e9," Gbytes"
    write (iu,901) "Vendor: "        ,trim(myvendor)
    write (iu,901) "Driver: "        ,trim(mydriver)
    
    write(iu,907)"                                                                               "
    write(iu,907)"*******************************************************************************"
    write(iu,907)"                                                                               "
    
    900 format (a,i0)
    901 format (a,a)
    902 format (a,i0,a)
    903 format (a,f16.8,a)
    904 format (a,2(i0,1x,'x',1x),i0)
    905 format (a,i0,'.',i0)
    906 format (a,l0)
    907 format (a)
    
    return
    
    end subroutine printDeviceProperties
  !$endif  
  subroutine print_raw_sync(iframe)
  
   implicit none
   
   integer, intent(in) :: iframe
  
   sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   open(unit=345,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted')
   write(345)rhoprint
   close(345)
   open(unit=346,file=trim(sevt2), &
    status='replace',action='write',access='stream',form='unformatted')
   write(346)velprint
   close(346)
   
  end subroutine print_raw_sync
  
  subroutine print_vtk_sync(iframe)
   implicit none
   
   integer, intent(in) :: iframe
   
   sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
    '_'//trim(write_fmtnumb(iframe)) // '.vti'
   sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
    '_'//trim(write_fmtnumb(iframe)) // '.vti'
   open(unit=345,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted')
   write(345)head1,ndatavtk(1),rhoprint,footervtk(1)
   close(345)
   open(unit=346,file=trim(sevt2), &
    status='replace',action='write',access='stream',form='unformatted')
   write(346)head2,ndatavtk(2),velprint,footervtk(2)
   close(346)
   
  end subroutine print_vtk_sync
  
  subroutine print_raw_async(iframe)
  
   implicit none
   
   integer, intent(in) :: iframe
  
   sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   open(unit=345,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted',&
    asynchronous='yes')
   write(345,asynchronous='yes')rhoprint
   
   open(unit=346,file=trim(sevt2), &
    status='replace',action='write',access='stream',form='unformatted',&
    asynchronous='yes')
   write(346,asynchronous='yes')velprint
   
   
  end subroutine print_raw_async
  
  subroutine print_vtk_async(iframe)
   implicit none
   
   integer, intent(in) :: iframe
   
   sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
    '_'//trim(write_fmtnumb(iframe)) // '.vti'
   sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
    '_'//trim(write_fmtnumb(iframe)) // '.vti'
    
   open(unit=345,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted',&
    asynchronous='yes')
   write(345,asynchronous='yes')head1,ndatavtk(1),rhoprint
   
   
   open(unit=780,file=trim(sevt2), &
    status='replace',action='write',access='stream',form='unformatted',&
    asynchronous='yes')
   write(780,asynchronous='yes')head2,ndatavtk(2),velprint
   
  end subroutine print_vtk_async
  
  subroutine close_print_async
  
   implicit none
   
   wait(345)
   if(lvtk)write(345)footervtk(1)
   close(345)
   
   
   wait(780)
   if(lvtk)write(780)footervtk(2)
   close(780) 
   
  end subroutine close_print_async
    
end program
