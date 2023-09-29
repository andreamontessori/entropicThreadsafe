program recursiveTSLB  
    !$if _OPENACC
    use openacc
    !$endif
    use prints
    use vars

    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!se possibile modulo con subs a chiusura
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!definizione equilibri entropici
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!newton-raphson per algoritmo di ricerca omega locale
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!benchmark: canale piatto turbolento, isotropic turbulence,
    !$if _OPENACC
    integer :: devNum
    integer(acc_device_kind) :: devType
    devType = acc_get_device_type()
    devNum=acc_get_device_num(devType)
    !$endif
       
    nlinks=8 !pari!
    tau=1._db
    cssq=1.0_db/3.0_db
    visc_LB=cssq*(tau-0.5_db)
    one_ov_nu=1.0_db/visc_LB
#ifdef _OPENACC
    ngpus=acc_get_num_devices(acc_device_nvidia)
#else
    ngpus=0
#endif

    !*******************************user parameters**************************
    nx=64
    ny=64
    nsteps=2000000
    stamp=2000000
    fx=1.0_db*10.0**(-6)
    fy=0.0_db
    allocate(p(0:nlinks))
    allocate(f(0:nx+1,0:ny+1,0:8))
    allocate(rho(1:nx,1:ny),u(1:nx,1:ny),v(1:nx,1:ny),pxx(1:nx,1:ny),pyy(1:nx,1:ny),pxy(1:nx,1:ny))
    allocate(isfluid(1:nx,1:ny)) !,omega_2d(1:nx,1:ny)) 
    if(lprint)then
      allocate(rhoprint(1:nx,1:ny,1:nz))
      allocate(velprint(1:3,1:nx,1:ny,1:nz))
      rhoprint(1:nx,1:ny,1:nz)=0.0
      velprint(1:3,1:nx,1:ny,1:nz)=0.0
    endif
    !ex=(/0,1,0,-1,0,1,-1,-1,1/)
    !ey=(/0,0,1,0,-1,1,1,-1,-1/)

    p=(/4.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/9.0_db,1.0_db/36.0_db, &
    1.0_db/36.0_db,1.0_db/36.0_db,1.0_db/36.0_db/)
    omega=1.0_db/tau
    
    !*****************************************geometry************************
    isfluid=1
    isfluid(1,:)=0 !EAST
    isfluid(nx,:)=0 !WEST
    isfluid(:,1)=0 !SOUTH 
    isfluid(:,ny)=0 !NORTH
    !*************************************initial conditions ************************    
    u=0.0_db
    v=0.0_db
    rho=1.0_db     !rho!
    !do ll=0,nlinks
    f(1:nx,1:ny,0)=p(0)*rho(:,:)!0.0_db
    f(1:nx,1:ny,1)=p(1)*rho(:,:)
    f(1:nx,1:ny,2)=p(2)*rho(:,:)
    f(1:nx,1:ny,3)=p(3)*rho(:,:)
    f(1:nx,1:ny,4)=p(4)*rho(:,:)
    f(1:nx,1:ny,5)=p(5)*rho(:,:)
    f(1:nx,1:ny,6)=p(6)*rho(:,:)
    f(1:nx,1:ny,7)=p(7)*rho(:,:)
    f(1:nx,1:ny,8)=p(8)*rho(:,:)
    !enddo
    !*************************************check data ************************ 
    write(6,*) '*******************LB data*****************'
    write(6,*) 'tau',tau
    write(6,*) 'omega',omega
    write(6,*) 'visc',visc_LB
    write(6,*) 'fx',fx
    write(6,*) 'cssq',cssq
    write(6,*) '*******************INPUT data*****************'
    write(6,*) 'nx',nx
    write(6,*) 'ny',ny
    write(6,*) 'lpbc',lpbc
    write(6,*) 'lprint',lprint
    write(6,*) 'lvtk',lvtk
    write(6,*) 'lasync',lasync
    write(6,*) 'nsteps',nsteps
    write(6,*) 'stamp',stamp
    write(6,*) 'max fx',huge(fx)
    write(6,*) 'available gpus',ngpus
    write(6,*) 'stabilizer pts',stab_points
    write(6,*) '*******************************************'


    !$acc data copy(p,rho,u,v,pxx,pxy,pyy,f,isfluid,rhoprint,velprint) async(1)
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
      !$acc kernels present(rhoprint,velprint,rho,u,v) async(1)
      !$acc loop independent collapse(2)  private(i,j)
      do j=1,ny
        do i=1,nx
          rhoprint(i,j,1)=real(rho(i,j),kind=4)
          velprint(1,i,j,1)=real(u(i,j),kind=4)
          velprint(2,i,j,1)=real(v(i,j),kind=4)
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
        !***********************************moment + neq pressor*********
        !$acc kernels async(1)
        !$acc loop collapse(2) private(fneq1,feq) 
        do j=1,ny
            do i=1,nx
                if(isfluid(i,j).eq.1)then
                    pxx(i,j)=0.0_db
                    pyy(i,j)=0.0_db
                    pxy(i,j)=0.0_db
                    rho(i,j) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
                    u(i,j) = (f(i,j,1) +f(i,j,5) +f(i,j,8)-(f(i,j,3) +f(i,j,6) +f(i,j,7)))/rho(i,j)
                    v(i,j) = (f(i,j,5) +f(i,j,2) +f(i,j,6)-(f(i,j,7) +f(i,j,4) +f(i,j,8)))/rho(i,j)
                    ! non equilibrium pressor components
                    !1-3
                    feq=p(1)*rho(i,j)*0.5_db*(2.0_db+6.0_db*u(i,j)**2 - 3.0_db*v(i,j)**2+u(i,j)*(6.0_db-9.0_db*v(i,j)**2))
                    fneq1=f(i,j,1)-feq
                    pxx(i,j)=pxx(i,j)+fneq1
                    feq=p(3)*rho(i,j)*0.5_db*(2.0_db+6.0_db*u(i,j)**2 - 3.0_db*v(i,j)**2+u(i,j)*(-6.0_db+9.0_db*v(i,j)**2))
                    fneq1=f(i,j,3)-feq
                    pxx(i,j)=pxx(i,j)+fneq1
                    !2-4
                    feq=p(2)*rho(i,j)*0.5_db*(2.0_db + 6.0_db*v(i,j) + 6.0_db*v(i,j)**2-3.0_db*u(i,j)**2*(1.0_db+3.0_db*v(i,j)))
                    fneq1=f(i,j,2)-feq
                    pyy(i,j)=pyy(i,j)+fneq1
                    feq=p(4)*rho(i,j)*0.5_db*(2.0_db - 6.0_db*v(i,j) + 6.0_db*v(i,j)**2 + u(i,j)**2*(-3.0_db+9.0_db*v(i,j)))
                    fneq1=f(i,j,4)-feq
                    pyy(i,j)=pyy(i,j)+fneq1
                    !5-7
                    feq=p(5)*rho(i,j)*(1.0_db + 3.0_db*v(i,j) + 3.0_db*v(i,j)**2 + u(i,j)**2*(3.0_db+9.0_db*v(i,j)) + u(i,j)*(3.0_db+9.0_db*v(i,j)+9.0_db*v(i,j)**2))
                    fneq1=f(i,j,5)- feq
                    pxx(i,j)=pxx(i,j)+fneq1
                    pyy(i,j)=pyy(i,j)+fneq1
                    pxy(i,j)=pxy(i,j)+fneq1
                    feq=-p(7)*rho(i,j)*(-1.0_db+3.0_db*v(i,j) - 3.0_db*v(i,j)**2 + u(i,j)**2*(-3.0_db+9.0_db*v(i,j)) + u(i,j)*(3.0_db-9.0_db*v(i,j)+9.0_db*v(i,j)**2))
                    fneq1=f(i,j,7)- feq
                    pxx(i,j)=pxx(i,j)+fneq1
                    pyy(i,j)=pyy(i,j)+fneq1
                    pxy(i,j)=pxy(i,j)+fneq1
                    !6-8
                    feq=p(6)*rho(i,j)*(1.0_db + 3.0_db*v(i,j) + 3.0_db*v(i,j)**2 + u(i,j)**2*(3.0_db+9.0_db*v(i,j)) - 3.0_db*u(i,j)*(1.0_db+3.0_db*v(i,j)+3.0_db*v(i,j)**2))
                    fneq1=f(i,j,6)-feq
                    pxx(i,j)=pxx(i,j)+fneq1
                    pyy(i,j)=pyy(i,j)+fneq1
                    pxy(i,j)=pxy(i,j)-fneq1
                    feq=p(8)*rho(i,j)*(1.0_db + u(i,j)**2*(3.0_db-9.0_db*v(i,j)) - 3.0_db*v(i,j)+3.0_db*v(i,j)**2 + u(i,j)*(3.0_db-9.0_db*v(i,j)+9.0_db*v(i,j)**2))
                    fneq1=f(i,j,8)-feq
                    pxx(i,j)=pxx(i,j)+fneq1
                    pyy(i,j)=pyy(i,j)+fneq1
                    pxy(i,j)=pxy(i,j)-fneq1
                endif
            enddo 
        enddo
        !$acc end kernels
      !***********************************PRINT************************
        if(mod(step,stamp).eq.0) write(6,'(a,i8)')'step : ',step
        if(lprint)then
            if(mod(step,stamp).eq.0)then
                iframe=iframe+1
                !$acc wait(1)
                !$acc kernels present(rhoprint,velprint,rho,u,v) async(1)
                !$acc loop independent collapse(2)  private(i,j)
                do j=1,ny
                    do i=1,nx
                      rhoprint(i,j,1)=real(rho(i,j),kind=4)
                      velprint(1,i,j,1)=real(u(i,j),kind=4)
                      velprint(2,i,j,1)=real(v(i,j),kind=4)
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
      !***********collision + no slip + forcing: fused implementation*********
        !$acc kernels async(1)
            !$acc loop collapse(2)
            do j=1,ny
                do i=1,nx 
                    if(isfluid(i,j).eq.1)then   
                        !oneminusuu= -uu !1.0_db - uu
                        !0
                        feq=p(0)*rho(i,j)*(1.0_db -3.0_db*0.5_db*(u(i,j)**2+v(i,j)**2))
                        fneq1=-3.0_db*0.5_db*(pxx(i,j)+pyy(i,j))
                        f(i,j,0)=feq + (1.0_db-omega)*p(0)*fneq1
                        !1
                        fneq1=3.0_db*pxx(i,j)-1.5_db*(pyy(i,j) + 3.0_db*pyy(i,j)*u(i,j)+6.0_db*pxy(i,j)*v(i,j))
                        feq=p(1)*rho(i,j)*0.5_db*(2.0_db+6.0_db*u(i,j)**2 - 3.0_db*v(i,j)**2+u(i,j)*(6.0_db-9.0_db*v(i,j)**2))
                        f(i+1,j,1)= feq + (1.0_db-omega)*p(1)*fneq1 + fx*p(1)/cssq 
                        !3
                        feq=p(3)*rho(i,j)*0.5_db*(2.0_db+6.0_db*u(i,j)**2 - 3.0_db*v(i,j)**2+u(i,j)*(-6.0_db+9.0_db*v(i,j)**2))
                        fneq1=3.0_db*pxx(i,j) + 1.5_db*pyy(i,j)*(-1.0_db+3.0_db*u(i,j)) + 9.0_db*pxy(i,j)*v(i,j)
                        f(i-1,j,3)= feq + (1.0_db-omega)*p(3)*fneq1  - fx*p(3)/cssq 
                        !2
                        fneq1=-1.5_db*(pxx(i,j) - 2.0_db*pyy(i,j) + 6.0_db*pxy(i,j)*u(i,j) +3.0_db*pxx(i,j)*v(i,j))
                        feq=p(2)*rho(i,j)*0.5_db*(2.0_db + 6.0_db*v(i,j) + 6.0_db*v(i,j)**2-3.0_db*u(i,j)**2*(1.0_db+3.0_db*v(i,j)))
                        f(i,j+1,2)= feq + (1.0_db-omega)*p(2)*fneq1  + fy*p(2)/cssq 
                        !4
                        feq=p(4)*rho(i,j)*0.5_db*(2.0_db - 6.0_db*v(i,j) + 6.0_db*v(i,j)**2 + u(i,j)**2*(-3.0_db+9.0_db*v(i,j)))
                        fneq1=1.5_db*(-pxx(i,j) + 2.0_db*pyy(i,j) + 6.0_db*pxy(i,j)*u(i,j) + 3.0_db*pxx(i,j)*v(i,j))
                        f(i,j-1,4)= feq + (1.0_db-omega)*p(4)*fneq1  - fy*p(4)/cssq 
                        !5
                        fneq1=3.0_db*(pxx(i,j) + pyy(i,j) + 3.0_db*pyy(i,j)*u(i,j) + 3.0_db*pxx(i,j)*v(i,j) + pxy(i,j)*(3.0_db+6.0_db*u(i,j)+6.0_db*v(i,j)))
                        feq=p(5)*rho(i,j)*(1.0_db + 3.0_db*v(i,j) + 3.0_db*v(i,j)**2 + u(i,j)**2*(3.0_db+9.0_db*v(i,j)) + u(i,j)*(3.0_db+9.0_db*v(i,j)+9.0_db*v(i,j)**2))
                        f(i+1,j+1,5)= feq + (1.0_db-omega)*p(5)*fneq1 + fx*p(5)/cssq + fy*p(5)/cssq
                        !7
                        fneq1=-3.0_db*(pyy(i,j)*(-1.0_db + 3.0_db*u(i,j)) + pxx(i,j)*(-1.0_db+3.0_db*v(i,j)) + pxy(i,j)*(-3.0_db+6.0_db*u(i,j)+6.0_db*v(i,j)))
                        feq=-p(7)*rho(i,j)*(-1.0_db+3.0_db*v(i,j) - 3.0_db*v(i,j)**2 + u(i,j)**2*(-3.0_db+9.0_db*v(i,j)) + u(i,j)*(3.0_db-9.0_db*v(i,j)+9.0_db*v(i,j)**2))
                        f(i-1,j-1,7)=feq + (1.0_db-omega)*p(7)*fneq1 - fx*p(7)/cssq - fy*p(7)/cssq 
                        !6
                        fneq1=3.0_db*(pxx(i,j)+pyy(i,j)-3.0_db*pyy(i,j)*u(i,j)+pxy(i,j)*(-3.0_db+6.0_db*u(i,j)-6.0_db*v(i,j))+3.0_db*pxx(i,j)*v(i,j))  
                        feq=p(6)*rho(i,j)*(1.0_db + 3.0_db*v(i,j) + 3.0_db*v(i,j)**2 + u(i,j)**2*(3.0_db+9.0_db*v(i,j)) - 3.0_db*u(i,j)*(1.0_db+3.0_db*v(i,j)+3.0_db*v(i,j)**2))
                        f(i-1,j+1,6)= feq + (1.0_db-omega)*p(6)*fneq1 - fx*p(6)/cssq + fy*p(6)/cssq 
                        !8
                        fneq1=3.0_db*(pxx(i,j) + pyy(i,j) + 3.0_db*pyy(i,j)*u(i,j) - 3.0_db*pxx(i,j)*v(i,j) + pxy(i,j)*(-3.0_db-6.0_db*u(i,j)+6.0_db*v(i,j)))
                        feq=p(8)*rho(i,j)*(1.0_db + u(i,j)**2*(3.0_db-9.0_db*v(i,j)) - 3.0_db*v(i,j)+3.0_db*v(i,j)**2 + u(i,j)*(3.0_db-9.0_db*v(i,j)+9.0_db*v(i,j)**2))
                        f(i+1,j-1,8)=feq + (1.0_db-omega)*p(8)*fneq1 + fx*p(8)/cssq - fy*p(8)/cssq 
                    endif
                enddo
            enddo
        !********************************************bcs no slip*****************************************!
            !$acc loop independent 
            do j=1,ny
                !$acc loop independent 
                do i=1,nx
                    if(isfluid(i,j).eq.0)then

                        f(i+1,j-1,8)=f(i,j,6)!gpc 
                        f(i-1,j-1,7)=f(i,j,5)!hpc

                        f(i-1,j+1,6)=f(i,j,8)!gpc 
                        f(i+1,j+1,5)=f(i,j,7)!hpc 


                        f(i,j-1,4)=f(i,j,2)!gpc 
                        f(i-1,j,3)=f(i,j,1)!hpc 

                        f(i,j+1,2)=f(i,j,4)!gpc 
                        f(i+1,j,1)=f(i,j,3)!hpc 
                    endif
                enddo
            enddo
        !$acc end kernels
        
      !******************************************call periodic bcs: always after fused************************
        
          if(lpbc)then    
              !periodic along x  
              !$acc kernels async(1)
              !$acc loop independent 
              do j=2,ny-1
                  if(j>2 .and. j<ny-1)then
                      f(2,j,1)=f(nx,j,1)
                      f(2,j,5)=f(nx,j,5)
                      f(2,j,8)=f(nx,j,8)	
                      f(nx-1,j,3)=f(1,j,3)
                      f(nx-1,j,6)=f(1,j,6)
                      f(nx-1,j,7)=f(1,j,7)
                  else
                      if(j==2)then
                          f(2,j,1)=f(nx,j,1)
                          f(2,j,8)=f(nx,j,8)
                          f(nx-1,j,3)=f(1,j,3)
                          f(nx-1,j,7)=f(1,j,7)
                      endif
                      if(j==ny-1)then
                          f(2,j,1)=f(nx,j,1)
                          f(2,j,5)=f(nx,j,5)
                          f(nx-1,j,3)=f(1,j,3)
                          f(nx-1,j,6)=f(1,j,6)
                      endif
                  endif
              enddo
              !$acc end kernels
          endif
    enddo 
    
    !$acc wait
    call cpu_time(ts2)
    !$acc update host(rho,u,v)
    !$acc end data


    !************************************************test points**********************************************!
    write(6,*) 'u=',u(nx/2,ny/2) ,'v=',v(nx/2,ny/2),'rho',rho(nx/2,ny/2) !'rho=',rho(nx/2,1+(ny-1)/2),nx/2,1+(ny-1)/2
    write(6,*) 'u=',u(2,ny/2) ,'v=',v(2,ny/2),'rho',rho(2,ny/2)
    write(6,*) 'u=',u(1,ny/2) ,'v=',v(1,ny/2),'rho',rho(1,ny/2)
    
    write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
    write(6,*) 'glups: ',  real(nx)*real(ny)*real(nsteps)*real(1.d-9,kind=4)/(ts2-ts1)

    open(101, file = 'v.out', status = 'replace')
    do j=1,ny
        do i=1,nx
            write(101,*) v(i,j) 
        enddo
    enddo
    close(101) 
    
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
