program recursiveTSLB3D
   use mpi_template
#ifdef _OPENACC
   use openacc
#endif
   use prints
   use vars
   use bcs3D
   use lb_kernels
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

   integer :: dumpstep,narg,inumchar
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

   !relaxation time
   tau=0.5004_db
   omega=1.0_db/tau
   ! viscosity
   visc_LB=cssq*(tau-0.5_db)
   one_ov_nu=1.0_db/visc_LB
   !dump in
   dumpYN=0


#ifdef _OPENACC
   ngpus=acc_get_num_devices(acc_device_nvidia)
#else
   ngpus=0
#endif

   !*******************************user parameters and allocations**************************
   lx=500
   ly=500
   lz=2048
   nsteps=1000000
   stamp=10000
   stamp2D=2000
   dumpstep=500000
   fx=0.0_db*10.0**(-7)
   fy=0.0_db*10.0**(-5)
   fz=0.0_db*10.0**(-7)
   uwall=0.05
   radius=50.0
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

   !!DECIDI DECOMPOSIZIONE MPI: don't change, if mpi proc_j is set at command line
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

   !!!!!!! START MPI!!!!!!!!!
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
   !allocate(phi(0:nx+1,0:ny+1,0:nz+1))
   !allocate(g(0:nx+1,0:ny+1,0:nz+1,0:nlinks_advc))
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
   !*****************************geometry************************************************
   isfluid=1
   do k=1,nz
      gk=nz*coords(3)+k
      do j=1,ny
         gj=ny*coords(2)+j
         do i=1,nx
            gi=nx*coords(1)+i
                !if(gi==1)isfluid(i,j,k)=0   !left
                !if(gi==lx)isfluid(i,j,k)=0  !right
!                if(gj==1)isfluid(i,j,k)=0   !front
!                if(gj==ly)isfluid(i,j,k)=0  !rear
            if(gk==1)  isfluid(i,j,k)=0   !bottom
            if(gk==lz) isfluid(i,j,k)=0  !top
         enddo
      enddo
   enddo

   !setup domain decomposition among MPI process
   !******************************(only once at the beginning)************************
   call exchange_isf_sendrecv
   call exchange_isf_intpbc
   call exchange_isf_wait

   !*************************************initial conditions ************************
   call initial_conditions
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
   !*************************************copy data on device*****************
   !$acc data copy(step,lx,ly,lz,nx,ny,nz,coords,myoffset,f,isfluid, &
   !$acc& pxx,pyy,pzz,pxy,pxz,pyz,rho,u,v,w,rhoprint,velprint, &
   !$acc& intpbc_dir,num_links_pops,links_pops,datampi,f_datampi,uwall, &
   !$acc& send_extr,recv_extr,f_send_extr,f_recv_extr) &
   !$acc& create(send_buffmpi,recv_buffmpi,f_send_buffmpi, &
   !$acc& f_recv_buffmpi)
	! quali sono i buff effettivamente da tenere?



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
         !2d planes print
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
      call moments_TSLB
      if(ldiagnostic)call end_timing2("LB","moments")
      !
      !if(ldiagnostic)call start_timing2("LB","pbcs_phi")
      !call exchange_float_sendrecv
      !call exchange_float_intpbc
      !call exchange_float_wait
      !if(ldiagnostic)call end_timing2("LB","pbcs_phi")
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
      call fused_TSLB_1c
      if(ldiagnostic)call end_timing2("LB","fused")
      !***********************************pbcs boundary conditions ********************************!
      !call pbcs
      if(ldiagnostic)call start_timing2("LB","pbcs")
	  call exchange_pops_sendrecv
	  call exchange_pops_intpbc
	  call exchange_pops_wait
      if(ldiagnostic)call end_timing2("LB","pbcs")
      !if(ldiagnostic)call start_timing2("LB","pbcs_advc")
      !call exchange_pops_advc_sendrecv
      !call exchange_pops_advc_intpbc
      !call exchange_pops_advc_wait
      !if(ldiagnostic)call end_timing2("LB","pbcs_advc")
      ! thread-safe boundary condition setup
      if(ldiagnostic)call start_timing2("LB","bcs_TSLB")
      call bcs_turbulent_jet_meso
      if(ldiagnostic)call end_timing2("LB","bcs_TSLB")
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
