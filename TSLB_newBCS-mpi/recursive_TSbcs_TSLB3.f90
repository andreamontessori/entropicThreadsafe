program recursiveTSLB3D
#ifdef _OPENACC
    use openacc
#endif
    use prints
    use vars
    use bcs3D
    use profiling_m,   only : timer_init,itime_start, &
                        startPreprocessingTime,print_timing_partial, &
                        reset_timing_partial,printSimulationTime, &
                        print_timing_final,itime_counter,idiagnostic, &
                        ldiagnostic,start_timing2,end_timing2, &
                        set_value_ldiagnostic,set_value_idiagnostic, &
                        startSimulationTime,print_memory_registration, &
                        get_memory, &
#ifdef CUDA
                        get_totram,get_memory_cuda,print_memory_registration_cuda
#else
                        get_totram
#endif 
    
    implicit none
    integer :: dumpYN
    integer :: dumpstep,subchords(3),narg,inumchar
    logical :: mydiagnostic
    integer :: tdiagnostic
    real(kind=db) :: smemory,sram
    
    integer :: nplanes
    integer, allocatable, dimension(:) :: ndir,npoint

#ifdef _OPENACC
    integer :: devNum
    integer(acc_device_kind) :: devType
    devType = acc_get_device_type()
    devNum=acc_get_device_num(devType)
#endif

   
    tau=0.5004_db
    
    visc_LB=cssq*(tau-0.5_db)
    one_ov_nu=1.0_db/visc_LB
    dumpYN=0


#ifdef _OPENACC
    ngpus=acc_get_num_devices(acc_device_nvidia)
#else
    ngpus=0
#endif

    !*******************************user parameters and allocations**************************
        lx=256
        ly=256
        lz=1024
        nsteps=5000
        stamp=1000
        stamp2D=1000
        dumpstep=100000000
        fx=0.0_db*10.0**(-7)
        fy=0.0_db*10.0**(-5)
        fz=0.0_db*10.0**(-5)
		uwall=0.05
        lprint=.true.
        lvtk=.false.
        lraw=.true.
        lasync=.false.
        lpbc=.true.
        
        nplanes=3       !numero di piani da stampare
        allocate(ndir(nplanes),npoint(nplanes))
        ndir(1)=1       !perpendicolare all'asse x
        ndir(2)=2       !perpendicolare all'asse y
        ndir(3)=3       !perpendicolare all'asse z
        npoint(1)=lx/2  !nodo lungo l'asse x
        npoint(2)=ly/2  !nodo lungo l'asse y
        npoint(3)=lz/2  !nodo lungo l'asse z
        
        
        pbc_x=1  !(0=false 1=true)
        pbc_y=1
        pbc_z=0
        
        !!DECIDI DECOMPOSIZIONE MPI
        proc_x=1
        proc_y=1
        proc_z=1
#ifdef MPI        
        !leggi decomposizione da riga di comando
        narg = command_argument_count()
        if (narg /= 3) then
          write(6,*) 'error!'
          write(6,*) 'the command line should be'
          write(6,*) '[executable] [proc_x] [proc_y] [proc_z]'
          write(6,*) 'proc_x = decomposition along x'
          write(6,*) 'proc_y = decomposition along y'
          write(6,*) 'proc_z = decomposition along z'
          write(6,*) 'STOP!'
          stop
        endif
  
        do i = 1, narg
          call getarg(i, arg)
          if(i==1)then
            call copystring(arg,directive,mxln)
             proc_x=intstr(directive,mxln,inumchar)
             write(6,*) 'proc_x  = ',proc_x
          elseif(i==2)then
            call copystring(arg,directive,mxln)
             proc_y=intstr(directive,mxln,inumchar)
             write(6,*) 'proc_y  = ',proc_y
          elseif(i==3)then
            call copystring(arg,directive,mxln)
             proc_z=intstr(directive,mxln,inumchar)
             write(6,*) 'proc_z  = ',proc_z
          endif
        enddo
#endif
        
        !!!!!!! START MPI!!!!!!!!!!!!!!!!!!!
        call start_mpi
        
        ! start diagnostic if requested
        mydiagnostic=.true.
        tdiagnostic=1
        call set_value_ldiagnostic(mydiagnostic)
        call set_value_idiagnostic(tdiagnostic)
        if(ldiagnostic)then
         call timer_init()
         call startPreprocessingTime()  
       endif
        
        !!!!!!!SETUP MPI!!!!!!!!!!!!!!!!!!!
        call setup_mpi() 
        
        allocate(f(0:nx+1,0:ny+1,0:nz+1,0:nlinks))
        allocate(rho(1:nx,1:ny,1:nz),u(1:nx,1:ny,1:nz),v(1:nx,1:ny,1:nz),w(1:nx,1:ny,1:nz))
        allocate(pxx(1:nx,1:ny,1:nz),pxy(1:nx,1:ny,1:nz),pxz(1:nx,1:ny,1:nz),pyy(1:nx,1:ny,1:nz))
        allocate(pyz(1:nx,1:ny,1:nz),pzz(1:nx,1:ny,1:nz))
        allocate(isfluid(0:nx+1,0:ny+1,0:nz+1))
        if(lprint)then
          allocate(rhoprint(1:nx,1:ny,1:nz))
          allocate(velprint(1:3,1:nx,1:ny,1:nz))
          rhoprint(1:nx,1:ny,1:nz)=0.0
          velprint(1:3,1:nx,1:ny,1:nz)=0.0
        endif
        !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
        !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
        !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)

        omega=1.0_db/tau
    !*****************************************geometry************************
        isfluid=1
        do k=1,nz
            gk=nz*coords(3)+k
            do j=1,ny
              gj=ny*coords(2)+j
              do i=1,nx
                gi=nx*coords(1)+i
!                if(gi==1)isfluid(i,j,k)=0   !left
!                if(gi==lx)isfluid(i,j,k)=0  !right
!                if(gj==1)isfluid(i,j,k)=0   !front 
!                if(gj==ly)isfluid(i,j,k)=0  !rear
                if(gk==1)  isfluid(i,j,k)=0   !bottom
                if(gk==lz) isfluid(i,j,k)=0  !top
              enddo
            enddo
          enddo
          
        !setup domain decomposition omong MPI process 
        !******************************(solo all'inizio)************************  
        call exchange_isf_sendrecv
        call exchange_isf_intpbc
        call exchange_isf_wait
        
    !*************************************initial conditions ************************    
        u=0.0_db
        v=0.0_db
        w=0.0_db
        !devo trovare quali processi hanno in carico i nodi lungo il piano gk=1
        gk=1
        !subchords(1)=(gi-1)/nx
        !subchords(2)=(gj-1)/ny
        subchords(3)=(gk-1)/nz !!subchord coordinate della decomposizione che ha al suo interno gk=1
        !sono buoni tutti i processi che hanno subchords(3)==coords(3)
        if(subchords(3)==coords(3))then
          do j=1,ny
            gj=ny*coords(2)+j
            do i=1,nx
              gi=nx*coords(1)+i
              if((float(gi)-lx/2.0)**2 + (float(gj)-ly/2.0)**2<=10**2)then
                !call random_number(rrx)
                !call random_number(rry)
                !sto sul piano gk=1
                gk=1
                !myoffset(3) è il mio offset lungo z del mio sottodominio MPI e mi ridà il valore di k nel sottodominio
                k=gk-myoffset(3)
                !è uno pseudo generatore che da un numero randomico partendo da 4 integer come seed
                !devi fare in modo che ogni lattice point ad ogni time step abbia seed diversi
                !quindi uso come seed la posizione i j k e il timestep come quarto seed lo metto ad cazzum
                !il fatto che tutti i seed siano diversi è perchè può essere chiamata da più threads contemporaneamente
                !invece se tu hai un unico seed lo devi mettere in save per tutti i threads e poi dipende da chi chiama prima (dipende dall'ordine di chiamata)
                rrx=rand_noseeded(gi,gj,gk,524)
                rry=rand_noseeded(gi,gj,gk,1732)
                w(i,j,k)=uwall + 0.02*sqrt(-2.0*log(rry))*cos(2*3.1415926535897932384626433832795028841971*rrx)
              endif
            enddo
          enddo
        endif
        rho=1.0_db  !tot dens
        
        !rho(nx/2,ny/2,nz-10)=1.2
        !do ll=0,nlinks
        if(dumpYN.eq.0)then
            do k=1,nz
                  gk=nz*coords(3)+k
                  do j=1,ny
                      gj=ny*coords(2)+j
                      do i=1,nx
                          gi=nx*coords(1)+i
                          if(isfluid(i,j,k).eq.1)then
                              !0
                              
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
                          endif
                      enddo
                  enddo
              enddo
          elseif(dumpYN.eq.1)then
              !call read_distros_1c_3d
              call dostop('read_distros_1c_3d NOT implemented')
          endif
    
    !*************************************check data ************************ 
    if(myrank==0)then
#ifdef MPI
        write(6,*) 'MPI VERSION COMPILED'
#else
        write(6,*) 'SERIAL VERSION COMPILED'
#endif
        write(6,*) '*******************LB data*****************'
        write(6,*) 'tau',tau
        write(6,*) 'omega',omega
        write(6,*) 'visc',visc_LB
        write(6,*) 'fx',fx
        write(6,*) 'fy',fy
        write(6,*) 'fz',fz
        write(6,*) 'cssq',cssq
        write(6,*) '*******************INPUT data*****************'
        write(6,*) 'lx',lx
        write(6,*) 'ly',ly
        write(6,*) 'ly',lz
        write(6,*) 'lpbc',lpbc
        write(6,*) 'lprint',lprint
        write(6,*) 'lvtk',lvtk
        write(6,*) 'lasync',lasync
        write(6,*) 'nsteps',nsteps
        write(6,*) 'stamp',stamp
        write(6,*) 'dumpYN', dumpYN
        write(6,*) '*******************************************'
        ! info gpu
        call get_memory_gpu(mymemory,totmemory)
        call print_memory_registration_gpu(6,'DEVICE memory occupied at the start', &
        'total DEVICE memory',mymemory,totmemory)
        call flush(6)
     endif
     
    !$acc data copy(step,lx,ly,lz,nx,ny,nz,coords,myoffset,f,isfluid, &
       !$acc& pxx,pyy,pzz,pxy,pxz,pyz,rho,u,v,w,rhoprint,velprint, &
       !$acc& intpbc_dir,num_links_pops,links_pops,datampi,f_datampi,uwall, &
       !$acc& dest_extr,source_extr,f_dest_extr,f_source_extr) & 
       !$acc& create(dest_buffmpi,source_buffmpi,f_dest_buffmpi, &
       !$acc& f_source_buffmpi)
       
  
   
    
#ifdef _OPENACC      
        call printDeviceProperties(ngpus,devNum,devType,6)
#endif
    iframe=0
    iframe2D=0
    if(myrank==0)then
      write(6,'(a,i8,a,i8,3f16.4)')'start step : ',0,' frame ',iframe
      call flush(6)
    endif
    
    if(lprint)then  
        call init_output(1,lvtk,lraw)
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
        call driver_print_vtk_sync(iframe)
      endif
      if(lraw)then
        call driver_print_raw_sync(iframe)
        !scrivi i piani 2D
        call driver_print_raw_sync2D(iframe2D,nplanes,ndir,npoint)
      endif
    endif
    
    ! start diagnostic if requested
        if(ldiagnostic)then
           !call print_timing_partial(1,1,itime_start,6)
           !call reset_timing_partial()
          call startSimulationTime()
          call get_memory(smemory) 
          call get_totram(sram)
          call print_memory_registration(6,&
          'Occupied memory after setup MPI','Total memory',smemory,sram)
        endif
      
    !*************************************time loop************************  
    call cpu_time(ts1)
    do step=1,nsteps 
        !***********************************moments collision bbck + forcing************************ 
          if(ldiagnostic)call start_timing2("LB","moments")
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
          if(ldiagnostic)call end_timing2("LB","moments")
          ! thread-safe boundary condition setup
          if(ldiagnostic)call start_timing2("LB","bcs_TSLB")
          call bcs_TSLB_only_z_turbojet
          if(ldiagnostic)call end_timing2("LB","bcs_TSLB")
          !call bcs_TSLB_turbojet
          !
        !***********************************Print on files 3D************************
          if(mod(step,stamp).eq.0 .or. mod(step,stamp2D).eq.0)then
            if(myrank==0)write(6,'(a,i8)')'stamp step : ',step
          endif
          if(lprint)then
            if(mod(step,stamp).eq.0 .or. mod(step,stamp2D).eq.0)then
                if(ldiagnostic)call start_timing2("IO","print")
                
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
            endif
            
            if(mod(step,stamp).eq.0)then
              iframe=iframe+1
              if(lvtk)then
                call driver_print_vtk_sync(iframe)
              endif
              if(lraw)then
                call driver_print_raw_sync(iframe)
              endif
            endif
            !***********************************Print on files 2D************************
            if(mod(step,stamp2D).eq.0 .and. lraw)then
                iframe2D=iframe2D+1
                call driver_print_raw_sync2D(iframe2D,nplanes,ndir,npoint)
            endif
            if(mod(step,stamp).eq.0 .or. mod(step,stamp2D).eq.0)then
              if(ldiagnostic)call end_timing2("IO","print")
            endif
          endif
        !***********************************dump f************************
          if(mod(step,dumpstep).eq.0) then
                write(6,'(a,i8)')'dump step at : ',step
                if(ldiagnostic)call start_timing2("IO","dump_distros")
                !$acc update host(f) 
                call dump_distros_1c_3d
                if(ldiagnostic)call end_timing2("IO","dump_distros")
          endif
        !***********************************collision + no slip + forcing: fused implementation*********
          if(ldiagnostic)call start_timing2("LB","fused")
          !$acc kernels
          !$acc loop collapse(3) private(feq,uu,temp,udotc)
          do k=1,nz
              do j=1,ny
                  do i=1,nx
                      !if(isfluid(i,j,k).eq.1)then
                          !0
                          feq=(-4*rho(i,j,k)*(-2 + 3*u(i,j,k)**2 + 3*v(i,j,k)**2 + 3*w(i,j,k)**2))/27.
                          fneq1=(-3*(pxx(i,j,k) + pyy(i,j,k) + pzz(i,j,k)))/2.
                          f(i,j,k,0)=feq + (1-omega)*fneq1*p0
                          
                          !1
                          
                          feq=(rho(i,j,k)*(2 + 6*u(i,j,k) + 6*u(i,j,k)**2 - 3*v(i,j,k)**2 - 9*u(i,j,k)*v(i,j,k)**2 - 3*(1 + 3*u(i,j,k))*w(i,j,k)**2))/27.
                          fneq1=(3*(2*pxx(i,j,k) - pzz(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - pyy(i,j,k)*(1 + 3*u(i,j,k)) - 6*pxy(i,j,k)*v(i,j,k) - 6*pxz(i,j,k)*w(i,j,k)))/2.
                          f(i+1,j,k,1)=feq + (1-omega)*fneq1*p1 
                          
                          !2
                          feq=(rho(i,j,k)*(2 - 3*v(i,j,k)**2 - 3*w(i,j,k)**2 + 3*u(i,j,k)*(-2 + 2*u(i,j,k) + 3*v(i,j,k)**2 + 3*w(i,j,k)**2)))/27.
                          fneq1=(3*(2*pxx(i,j,k) - pzz(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + pyy(i,j,k)*(-1 + 3*u(i,j,k)) + 6*pxy(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*w(i,j,k)))/2.
                          f(i-1,j,k,2)=feq + (1-omega)*fneq1*p1 
                          
                          !3
                          
                          feq=(rho(i,j,k)*(2 - 3*u(i,j,k)**2*(1 + 3*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(2 + 2*v(i,j,k) - 3*w(i,j,k)**2)))/27.
                          fneq1=(-3*(pxx(i,j,k) - 2*pyy(i,j,k) + pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) + 3*pxx(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i,j+1,k,3)=feq+ (1-omega)*fneq1*p1 
                          
                          !4
                          feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k)) - 3*w(i,j,k)**2 + 3*v(i,j,k)*(-2 + 2*v(i,j,k) + 3*w(i,j,k)**2)))/27.
                          fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(-1 + 3*v(i,j,k)) + 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i,j-1,k,4)=feq+ (1-omega)*fneq1*p1 
                          
                          !7
                          
                          feq=(2*rho(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) + u(i,j,k)*(3 + 9*v(i,j,k)*(1 + v(i,j,k)))) - 3*rho(i,j,k)*(1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) + 6*pyy(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*(1 + 2*u(i,j,k) + 2*v(i,j,k)) + pxx(i,j,k)*(2 + 6*v(i,j,k)) - 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i+1,j+1,k,7)=feq + (1-omega)*fneq1*p2 
                          
                          !8
                          feq=(2*rho(i,j,k)*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(-3 - 9*(-1 + v(i,j,k))*v(i,j,k))) + 3*rho(i,j,k)*(-1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=(-3*(-2*pyy(i,j,k) + pzz(i,j,k) + 6*pyy(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*(-1 + 2*u(i,j,k) + 2*v(i,j,k)) + pxx(i,j,k)*(-2 + 6*v(i,j,k)) - 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i-1,j-1,k,8)=feq + (1-omega)*fneq1*p2 
                          
                          !10
                          
                          feq=(2*rho(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)*(1 + v(i,j,k)))) + 3*rho(i,j,k)*(-1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) - 6*pyy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*(-1 + 2*u(i,j,k) - 2*v(i,j,k)) - 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(2 + 6*v(i,j,k)) + 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i-1,j+1,k,10)=feq+ (1-omega)*fneq1*p2 
                          
                          !9
                          feq=(2*rho(i,j,k)*(1 + u(i,j,k)**2*(3 - 9*v(i,j,k)) + 3*(-1 + v(i,j,k))*v(i,j,k) + u(i,j,k)*(3 + 9*(-1 + v(i,j,k))*v(i,j,k))) - 3*rho(i,j,k)*(1 + 3*u(i,j,k) - 3*v(i,j,k))*w(i,j,k)**2)/108.
                          fneq1=(-3*(-2*pyy(i,j,k) + pzz(i,j,k) - 6*pyy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*(1 + 2*u(i,j,k) - 2*v(i,j,k)) - 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(-2 + 6*v(i,j,k)) + 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i+1,j-1,k,9)=feq+ (1-omega)*fneq1*p2 

                          !5
                          
                          feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k))))/27.
                          fneq1=(-3*(pxx(i,j,k) + pyy(i,j,k) - 2*pzz(i,j,k) + 6*pxz(i,j,k)*u(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pxx(i,j,k)*w(i,j,k) + 3*pyy(i,j,k)*w(i,j,k)))/2.
                          f(i,j,k+1,5)=feq+ (1-omega)*fneq1*p1 
                          
                          !6
                          feq=(rho(i,j,k)*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 6*pxz(i,j,k)*u(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(-1 + 3*w(i,j,k))))/2.
                          f(i,j,k-1,6)=feq+ (1-omega)*fneq1*p1 

                          !15
                          
                          feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) - 3*pyy(i,j,k)*u(i,j,k) + 6*pzz(i,j,k)*u(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) - 3*pyy(i,j,k)*w(i,j,k) + 6*pxz(i,j,k)*(1 + 2*u(i,j,k) + 2*w(i,j,k)) + pxx(i,j,k)*(2 + 6*w(i,j,k))))/2.
                          f(i+1,j,k+1,15)=feq+ (1-omega)*fneq1*p2 
                          
                          !16
                          feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(2 - 6*w(i,j,k)) + 3*pyy(i,j,k)*w(i,j,k) - 6*pxz(i,j,k)*(-1 + 2*u(i,j,k) + 2*w(i,j,k))))/2.
                          f(i-1,j,k-1,16)=feq+ (1-omega)*fneq1*p2 

                          !17
                          
                          feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*(-1 + 2*u(i,j,k) - 2*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(2 + 6*w(i,j,k))))/2.
                          f(i-1,j,k+1,17)=feq+ (1-omega)*fneq1*p2 
                          
                          !18
                          feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                          fneq1=(-3*(pyy(i,j,k) - 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*(1 + 2*u(i,j,k) - 2*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(-2 + 6*w(i,j,k))))/2.
                          f(i+1,j,k-1,18)=feq+ (1-omega)*fneq1*p2 

                          !11
                          
                          feq=(rho(i,j,k)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=(-3*pxx(i,j,k)*(1 + 3*v(i,j,k) + 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) - 3*pxy(i,j,k)*u(i,j,k) - 3*pxz(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*pyy(i,j,k)*w(i,j,k) + pyz(i,j,k)*(3 + 6*v(i,j,k) + 6*w(i,j,k)))
                          f(i,j+1,k+1,11)=feq+ (1-omega)*fneq1*p2 
                          
                          !12
                          feq=(rho(i,j,k)*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                          fneq1=(3*pxx(i,j,k)*(-1 + 3*v(i,j,k) + 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) + 3*pxy(i,j,k)*u(i,j,k) + 3*pxz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + pyz(i,j,k)*(3 - 6*v(i,j,k) - 6*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k))
                          f(i,j-1,k-1,12)=feq+ (1-omega)*fneq1*p2 

                          !13
                          
                          feq=(rho(i,j,k)*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                          fneq1=(-3*(pxx(i,j,k) - 2*pyy(i,j,k) + 6*pyz(i,j,k) - 2*pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) - 6*pxz(i,j,k)*u(i,j,k) + 3*pxx(i,j,k)*v(i,j,k) + 12*pyz(i,j,k)*v(i,j,k) - 6*pzz(i,j,k)*v(i,j,k) - 3*pxx(i,j,k)*w(i,j,k) + 6*pyy(i,j,k)*w(i,j,k) - 12*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i,j+1,k-1,13)=feq+ (1-omega)*fneq1*p2 
                          
                          !14
                          feq=(rho(i,j,k)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                          fneq1=(3*pxx(i,j,k)*(-1 + 3*v(i,j,k) - 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) + 3*pxy(i,j,k)*u(i,j,k) - 3*pxz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + pyz(i,j,k)*(-3 + 6*v(i,j,k) - 6*w(i,j,k)) + 3*pyy(i,j,k)*w(i,j,k))
                          f(i,j-1,k+1,14)=feq+ (1-omega)*fneq1*p2 

                          !19
                          feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                          fneq1=3*(pxx(i,j,k) + (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) + 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j+1,k+1,19)=feq + (1-omega)*fneq1*p3 
                          
                          !20
                          feq=(rho(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                          fneq1=-3*((pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(-1 + 3*u(i,j,k)) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + pxz(i,j,k)*(-3 + 6*u(i,j,k)) + 3*(2*pxy(i,j,k) + 3*pxz(i,j,k) + 2*pyz(i,j,k) + pzz(i,j,k))*v(i,j,k) + 3*(3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k) + pxx(i,j,k)*(-1 + 3*v(i,j,k) + 3*w(i,j,k)))
                          f(i-1,j-1,k-1,20)=feq+ (1-omega)*fneq1*p3 

                          !21
                          feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                          fneq1=3*(pxx(i,j,k) - 3*pxy(i,j,k)*(1 + 2*u(i,j,k)) + (pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) - 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j-1,k+1,21)=feq+ (1-omega)*fneq1*p3 
                          
                          !22
                          feq=(rho(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                          fneq1=3*(pxx(i,j,k) + 3*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k) - 3*(2*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) - 3*(pxx(i,j,k) - 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i-1,j+1,k-1,22)=feq+ (1-omega)*fneq1*p3 

                          !23
                          feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                          fneq1=3*(pxx(i,j,k) + 3*pxy(i,j,k) - 3*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k) - 3*(2*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k) - 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) + 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i-1,j-1,k+1,23)=feq+ (1-omega)*fneq1*p3 

                          !24
                          feq=(rho(i,j,k)*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                          fneq1=3*(pxx(i,j,k) - 3*pxz(i,j,k)*(1 + 2*u(i,j,k)) + (pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) - 3*(pxx(i,j,k) + 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j+1,k-1,24)=feq+ (1-omega)*fneq1*p3 

                          !25
                          feq=(rho(i,j,k)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                          fneq1=-3*(-pxx(i,j,k) - (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j-1,k-1,25)=feq + (1-omega)*fneq1*p3 

                          !26
                          feq=(rho(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                          fneq1=3*(pxx(i,j,k) - (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(-1 + 3*u(i,j,k)) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + pxz(i,j,k)*(-3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                          f(i-1,j+1,k+1,26)=feq + (1-omega)*fneq1*p3 
                      !endif
                  enddo
              enddo
          enddo
          !$acc end kernels
          if(ldiagnostic)call end_timing2("LB","fused")
        !***********************************pbcs boundary conditions ********************************!
        !call pbcs
        if(ldiagnostic)call start_timing2("LB","pbcs")
        call exchange_pops_sendrecv
        call exchange_pops_intpbc
        call exchange_pops_wait
        if(ldiagnostic)call end_timing2("LB","pbcs") 
    enddo 
    !$acc end data
    call cpu_time(ts2)
    
    if(ldiagnostic)then
      call printSimulationTime()
      call print_timing_final(idiagnostic,itime_counter, &
       itime_start,1,1,6)
      call get_memory(smemory) 
      call get_totram(sram)
      call print_memory_registration(6,&
      'Occupied memory after setup MPI','Total memory',smemory,sram)
    endif
    
    if(myrank==0)then
      write(6,*) 'time elapsed: ', ts2-ts1, ' s of your life time' 
      write(6,*) 'glups: ',  real(lx)*real(ly)*real(lz)*real(nsteps)/1.0e9/(ts2-ts1)
    endif
    
    call get_memory(smemory) 
    call get_totram(sram)
    call print_memory_registration(6,&
      'Occupied memory on HOST at the end','Total HOST memory',smemory,sram)
          
    call get_memory_gpu(mymemory,totmemory)
    call print_memory_registration_gpu(6,'DEVICE memory occupied at the end', &
     'total DEVICE memory',mymemory,totmemory)
    
end program
