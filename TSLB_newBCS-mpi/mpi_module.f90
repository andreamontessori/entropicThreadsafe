
 module mpi_template
  
  use vars, only: db,isf,i,j,k,nx,ny,nz,lx,ly,lz,rho,f,isfluid,l,ll,opp,&
   ex,ey,ez,nlinks,filenamevtk,namevarvtk,sevt1,sevt2,dir_out,write_fmtnumb2, &
   write_fmtnumb,headervtk,nheadervtk,vtkoffset,ndatavtk,footervtk, &
   rhoprint,velprint,space_fmtnumb
#ifdef _OPENACC
    use openacc
#endif

#ifdef MPI
   use mpi
#endif

  implicit none
  
  integer :: MYMPIREAL
  integer :: MYMPIINTS
  
  integer, save :: nprocs,myrank,lbecomm,localcomm
  
  integer :: mydev, ndev 
  integer :: file_offset    
  integer :: proc_x,proc_y,proc_z
  integer :: pbc_x,pbc_y,pbc_z 
  integer :: mem_stop
  logical :: rreorder 
  integer, parameter::  mpid=3     ! mpi dimension
  logical :: periodic(mpid) 
  integer :: prgrid(mpid)
  integer:: coords(mpid)
  integer:: up(mpid),down(mpid),left(mpid)
  integer:: front(mpid),rear(mpid),right(mpid)
  integer, allocatable, dimension(:) :: xinidom,xfindom
  integer, allocatable, dimension(:) :: yinidom,yfindom
  integer, allocatable, dimension(:) :: zinidom,zfindom
    
  integer :: right_dest_x,left_source_x  
  integer :: left_dest_x,right_source_x
  integer :: up_dest_y,down_source_y
  integer :: down_dest_y, up_source_y
  integer :: front_dest_z,rear_source_z
  integer :: rear_dest_z,front_source_z
	
  integer, dimension(mpid) :: myoffset,lsizes,start_idx,gsizes,end_idx
  
  integer, parameter :: nlinksmpi=26
                                     
  integer, dimension(0:nlinksmpi), parameter :: &
      ! 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
   exmpi=(/0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1/)
  integer, dimension(0:nlinksmpi), parameter :: &
   eympi=(/0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1,-1, 1,-1, 1/)
  integer, dimension(0:nlinksmpi), parameter :: &
   ezmpi=(/0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1/)
  integer, dimension(0:nlinksmpi), parameter ::&
  oppmpi=(/0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17,20,19,22,21,24,23,26,25/)
  
  
  integer, save :: nlinks_faces,nlinks_edges,nlinks_corners,nlinks_max
  integer, save :: nfaces,nedges,ncorners
  logical, save :: lintbb=.false.
  integer, dimension(nlinksmpi) :: dest_dir,source_dir
  logical, dimension(nlinksmpi) :: ldestpop_dir,lsourcepop_dir
  logical, dimension(nlinksmpi) :: ldest_dir,lsource_dir
  logical, dimension(nlinksmpi) :: lintpbc_dir,lintpbcpop_dir
  logical, dimension(nlinksmpi) :: lintbb_dir
  logical, dimension(3,nlinksmpi) :: intpbc_dir
  integer, dimension(3,nlinksmpi) :: dest_dir_coord,source_dir_coord
  
  integer, allocatable, save, dimension(:,:) :: links_faces,links_edges, &
   links_corners,links_pops
  integer, allocatable, save, dimension(:,:) :: dest_extr,source_extr
  integer, allocatable, save, dimension(:,:) :: f_dest_extr,f_source_extr
  integer, allocatable, save, dimension(:,:) :: b_dest_extr,b_source_extr
  
  integer, save, dimension(nlinksmpi) :: num_extr,f_num_extr,b_num_extr,i_num_extr
  integer, save :: numtot_extr,f_numtot_extr,b_numtot_extr,i_numtot_extr
  integer, allocatable, save, dimension(:) :: num_links_pops
  
  integer, dimension(13) :: datampi,f_datampi,b_datampi,i_datampi
  integer, parameter :: num_f_datampi=1
  integer, parameter :: num_b_datampi=1
  
  real(kind=db), allocatable, save, dimension(:) :: dest_buffmpi,source_buffmpi
  integer, dimension(nlinksmpi) :: nbuffmpi_dest,nbuffmpi_source
  
  real(kind=db), allocatable, save, dimension(:) :: f_dest_buffmpi,f_source_buffmpi
  integer, dimension(nlinksmpi), save :: f_nbuffmpi_dest,f_nbuffmpi_source
  
  real(kind=db), allocatable, save, dimension(:) :: b_dest_buffmpi,b_source_buffmpi
  integer, dimension(nlinksmpi), save :: b_nbuffmpi_dest,b_nbuffmpi_source
  
  integer(kind=isf), allocatable, save, dimension(:) :: i_dest_buffmpi,i_source_buffmpi
  integer, dimension(nlinksmpi), save :: i_nbuffmpi_dest,i_nbuffmpi_source
  
  integer, dimension(nlinksmpi), save :: reqs_dest,reqs_source,mpitag
  
  integer, dimension(nlinksmpi), save :: f_reqs_dest,f_reqs_source,f_mpitag
  integer, dimension(nlinksmpi), save :: b_reqs_dest,b_reqs_source,b_mpitag
  integer, dimension(nlinksmpi), save :: i_reqs_dest,i_reqs_source,i_mpitag
  
  integer, parameter :: nbuff=2  
  logical, parameter :: lbuff=.true.
  
  
  
 contains
      
      subroutine start_mpi
      
      implicit none
      
      integer:: ierr
      
#ifdef MPI
!
      call mpi_init(ierr)
      call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)
      call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
      
#else
      nprocs=1
      myrank=0
      proc_x=1
      proc_y=1
      proc_z=1
#endif      

      
      end subroutine start_mpi
!
      subroutine setup_mpi


      implicit none
       
     
!
      integer:: uni,lopp,idrank,oi,oj,ok
      integer:: ierr, ijlen                    ! mpi variables
      character*15 hname
      character*17 file_name1
      character*17 file_name2
      character*17 file_name3
      character*15 file_name5
      integer,dimension(mpid) :: temp_coord
      logical :: lcheck=.false.
!
      real(db):: knorm
!     
#ifdef MPI
      if(db==4)then
        MYMPIREAL = MPI_REAL
      elseif(db==8)then
        MYMPIREAL = MPI_DOUBLE_PRECISION
      else
        write(6,*)'ERROR db not defined'
        stop
      endif
      if(isf==4)then
        MYMPIINTS=MPI_INTEGER
      elseif(isf==1)then
        MYMPIINTS=MPI_INTEGER1
      else
        write(6,*)'ERROR isf not defined'
        stop
      endif
#endif
      knorm = 1.0/1024.0
      
      
      nx = lx/proc_x
      ny = ly/proc_y
      nz = lz/proc_z
      
!
! some check
      lcheck=.false.
      if((nx*proc_x).NE.lx) then
         write(6,*) "ERROR: global and local size along x not valid!!" & 
     &                      , lx, nx, proc_x
         lcheck=.true.
      endif
!
      if((ny*proc_y).NE.ly) then
         write(6,*) "ERROR: global and local size along y not valid!!" &
     &                      , ly, ny, proc_y
         lcheck=.true.
      endif
!
      if((nz*proc_z).NE.lz) then
         write(6,*) "ERROR: global and local size along z not valid!!" &
     &                      , lz, nz, proc_z
         lcheck=.true.
      endif
      
      if(lcheck)then
#ifdef MPI      
         call MPI_finalize(ierr)
#endif
         stop
      endif
!
#ifdef _OPENACC
      ndev= acc_get_num_devices(acc_device_nvidia)
      if(ndev == 0) then
         if(myrank==0)write(6,*) 'WARNINIG: No GPUs found:', ndev
         call dostop
      endif
      mydev = mod(myrank, ndev)
      if(myrank==0)write(6,'(i4,a,i4,a)') ndev,' GPUs available for ', &
       nprocs,' MPI processes'
      call acc_set_device_num(mydev,acc_device_nvidia)
      write(6,'(a,i4,a,i4)') 'The MPI process ',myrank,' is using the GPU ',mydev
#endif
!

      
      
      rreorder=.false.
        
      periodic(1) = (pbc_x==1)
      periodic(2) = (pbc_y==1)
      periodic(3) = (pbc_z==1)
        
      prgrid(1) = proc_x!proc_x
      prgrid(2) = proc_y!proc_y
      prgrid(3) = proc_z
      
      if(myrank==0)write(6,'(a,3i4)')'MPI processes= ',proc_x,proc_y,proc_z   
      if(myrank==0)write(6,'(a,3i4)')'pbc applied= ',pbc_x,pbc_y,pbc_z
      call flush(6)
      
      
#ifdef MPI
      call mpi_barrier(MPI_COMM_WORLD,ierr)
! 
! set the gpu to the task id
      call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                                 MPI_INFO_NULL, localcomm, ierr)
      call MPI_Comm_rank(localcomm, mydev, ierr)
      call MPI_get_processor_name(hname,ijlen,ierr)

!
! check
      if ((proc_x*proc_y*proc_z).ne.nprocs) then
        if (myrank.eq.0) then
          write(6,'(a,3i4)') 'ERROR: decomposed for x y z procs ', & 
                             proc_x,proc_y,proc_z
          write(6,'(a,i4,a)') 'ERROR: decomposed for', & 
                             proc_x*proc_y*proc_z, 'procs'
          write(6,'(a,i4,a)') 'ERROR: launched on', nprocs, 'processes'
        end if
        call dostop('ERROR mpi job lunched with wrong number of processes')
      end if
!
!
! building virtual topology
      call MPI_cart_create(mpi_comm_world, mpid, prgrid, &
                              periodic,rreorder,lbecomm,ierr)

      call MPI_comm_rank(lbecomm, myrank, ierr)
      call MPI_cart_coords(lbecomm, myrank, mpid, &
                            coords, ierr)
                          
!
      call mpi_barrier(MPI_COMM_WORLD,ierr)
!    

      !pause
      !********************************************************!
	  !*******************neighbour processes******************!
	  !********************************************************!
      !call MPI_ERRHANDLER_SET(lbecomm, MPI_ERRORS_RETURN, ierr)
      do l=1,nlinksmpi
        mpitag(l) = 400 + l
        f_mpitag(l) = 500 + l
        b_mpitag(l) = 600 + l
        i_mpitag(l) = 700 + l
        temp_coord(1) = coords(1) + exmpi(l) 
	    temp_coord(2) = coords(2) + eympi(l) 
	    temp_coord(3) = coords(3) + ezmpi(l) 
        !call MPI_Cart_rank(lbecomm, temp_coord, dest_dir(l),ierr)
        
        oi=temp_coord(1)
        oj=temp_coord(2)
        ok=temp_coord(3)
        !oi=mod(oi+nx-1,nx)+1
        if(periodic(1))then 
          oi=mod(oi+proc_x,proc_x)
        else
          oi=min(max(oi,0),proc_x-1)
        endif
        if(periodic(2))then
          oj=mod(oj+proc_y,proc_y)
        else
          oj=min(max(oj,0),proc_y-1)
        endif
        if(periodic(3))then
          ok=mod(ok+proc_z,proc_z)
        else
          ok=min(max(ok,0),proc_z-1)
        endif

        dest_dir_coord(1:3,l)=[oi,oj,ok]
        call MPI_Cart_rank(lbecomm, [oi,oj,ok], dest_dir(l),ierr)

        
        lopp=opp(l)
        temp_coord(1) = coords(1) + exmpi(lopp) 
	    temp_coord(2) = coords(2) + eympi(lopp) 
	    temp_coord(3) = coords(3) + ezmpi(lopp) 
        !call MPI_Cart_rank(lbecomm, temp_coord, source_dir(l),ierr)
        
        oi=temp_coord(1)
        oj=temp_coord(2)
        ok=temp_coord(3)
!        if(periodic(1))oi=mod(oi+proc_x,proc_x)
!        if(periodic(2))oj=mod(oj+proc_y,proc_y)
!        if(periodic(3))ok=mod(ok+proc_z,proc_z)
        if(periodic(1))then 
          oi=mod(oi+proc_x,proc_x)
        else
          oi=min(max(oi,0),proc_x-1)
        endif
        if(periodic(2))then
          oj=mod(oj+proc_y,proc_y)
        else
          oj=min(max(oj,0),proc_y-1)
        endif
        if(periodic(3))then
          ok=mod(ok+proc_z,proc_z)
        else
          ok=min(max(ok,0),proc_z-1)
        endif

        source_dir_coord(1:3,l)=[oi,oj,ok]
        call MPI_Cart_rank(lbecomm, [oi,oj,ok],source_dir(l),ierr)

      enddo
		
      file_offset = 0    !to check



#else

      coords=0
      
      do l=1,nlinksmpi
          !vado sempre a me stesso perche sono il solo processo 
          dest_dir_coord(1:3,l)=coords(1:3)
          dest_dir(l)=myrank
          
          source_dir_coord(1:3,l)=coords(1:3)
          source_dir(l)=myrank
       enddo

#endif
      !gestisci se fare o no send e receive (deve andare su nodi diversi e non sfondare il range coords se non periodico)
      do l=1,nlinksmpi
          
          ldest_dir(l)=(myrank .ne. dest_dir(l))
          if(ldest_dir(l))then
            temp_coord(1) = coords(1) + exmpi(l) 
	        temp_coord(2) = coords(2) + eympi(l) 
	        temp_coord(3) = coords(3) + ezmpi(l) 
        
            oi=temp_coord(1)
            oj=temp_coord(2)
            ok=temp_coord(3)
 
            if(periodic(1))then 
              oi=mod(oi+proc_x,proc_x)
            endif
            if(periodic(2))then
              oj=mod(oj+proc_y,proc_y)
            endif
            if(periodic(3))then
              ok=mod(ok+proc_z,proc_z)
            endif
            !se sfondo allora non devo inviare
            if(oi<0 .or. oj<0 .or. ok<0 .or. &
             oi>=proc_x .or. oj>=proc_y .or. ok>=proc_z)then
               ldest_dir(l)=.false.
            endif
          endif
            
          lsource_dir(l)=(myrank .ne. source_dir(l))
          if(lsource_dir(l))then
            lopp=opp(l)
            temp_coord(1) = coords(1) + exmpi(lopp) 
	        temp_coord(2) = coords(2) + eympi(lopp) 
	        temp_coord(3) = coords(3) + ezmpi(lopp)  
  
            oi=temp_coord(1)
            oj=temp_coord(2)
            ok=temp_coord(3)
            if(periodic(1))then 
              oi=mod(oi+proc_x,proc_x)
            endif
            if(periodic(2))then
              oj=mod(oj+proc_y,proc_y)
            endif
            if(periodic(3))then
              ok=mod(ok+proc_z,proc_z)
            endif
            !se sfondo allora non devo ricevere
            if(oi<0 .or. oj<0 .or. ok<0 .or. &
             oi>=proc_x .or. oj>=proc_y .or. ok>=proc_z)then
               lsource_dir(l)=.false.
            endif
          endif
            
     enddo
   
     
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      
      !gestisci se fare o no le pbc interne (deve andare sullo stesso nodo in dest e source e non sfondare il range coords se non periodico)
      do l=1,nlinksmpi
        lintpbc_dir(l)=((myrank==dest_dir(l)) .and. (myrank==source_dir(l)))
        intpbc_dir(1,l)=.false.
        intpbc_dir(2,l)=.false.
        intpbc_dir(3,l)=.false.
        if(lintpbc_dir(l))then
            temp_coord(1) = coords(1) + exmpi(l) 
	        temp_coord(2) = coords(2) + eympi(l) 
	        temp_coord(3) = coords(3) + ezmpi(l) 
        
            oi=temp_coord(1)
            oj=temp_coord(2)
            ok=temp_coord(3)
 
            if(periodic(1))then 
              oi=mod(oi+proc_x,proc_x)
            endif
            if(periodic(2))then
              oj=mod(oj+proc_y,proc_y)
            endif
            if(periodic(3))then
              ok=mod(ok+proc_z,proc_z)
            endif
            !non devo sfondare se non periodico nella direzione specifica l
            if(oi<0 .or. oj<0 .or. ok<0 .or. &
             oi>=proc_x .or. oj>=proc_y .or. ok>=proc_z)then
               lintpbc_dir(l)=.false.
            endif
        endif
        !se sono periodico nel mio processo allora non devo sfondare neanche quando ricevo
        if(lintpbc_dir(l))then
            lopp=opp(l)
            temp_coord(1) = coords(1) + exmpi(lopp) 
	        temp_coord(2) = coords(2) + eympi(lopp) 
	        temp_coord(3) = coords(3) + ezmpi(lopp)  
  
            oi=temp_coord(1)
            oj=temp_coord(2)
            ok=temp_coord(3)
            if(periodic(1))then 
              oi=mod(oi+proc_x,proc_x)
            endif
            if(periodic(2))then
              oj=mod(oj+proc_y,proc_y)
            endif
            if(periodic(3))then
              ok=mod(ok+proc_z,proc_z)
            endif
            !non devo sfondare se non periodico nella direzione specifica l
            if(oi<0 .or. oj<0 .or. ok<0 .or. &
             oi>=proc_x .or. oj>=proc_y .or. ok>=proc_z)then
               lintpbc_dir(l)=.false.
            endif
        endif
        !sto in ballo devo fare periodico interno se true
        !allora storo se sono periodico per direzione l e condizioni lungo i tre assi
        if(lintpbc_dir(l))then
          !se lungo il vettore l mi muovo lungo un asse e sono periodico allora intpbc_dir è true
          intpbc_dir(1,l)=(abs(exmpi(l))==1 .and.periodic(1))
          intpbc_dir(2,l)=(abs(eympi(l))==1 .and.periodic(2))
          intpbc_dir(3,l)=(abs(ezmpi(l))==1 .and.periodic(3))
        endif
      enddo
      
      
      !setto variabili utili in particolare per MPI-IO
      !offset in coordinate globali di ogni processo MPI
      myoffset(1) = coords(1)*nx
      myoffset(2) = coords(2)*ny
      myoffset(3) = coords(3)*nz
      !dimensione locale grid di ogni processo MPI
      lsizes(1)=nx
      lsizes(2)=ny
      lsizes(3)=nz
      !dimensione globale grid di ogni processo MPI
      gsizes(1)=lx
      gsizes(2)=ly
      gsizes(3)=lz
      !inizio dimensione locale grid in coordinate globali di ogni processo MPI
      start_idx(1)=myoffset(1)+1
      start_idx(2)=myoffset(2)+1
      start_idx(3)=myoffset(3)+1
      !fine dimensione locale grid in coordinate globali di ogni processo MPI
      end_idx(1) = myoffset(1)+lsizes(1)
      end_idx(2) = myoffset(2)+lsizes(2)
      end_idx(3) = myoffset(3)+lsizes(3)
      
      !questo mi serve per la funzione GET_COORD_POINT 
      !funzione di debug che dalle coordinate generali x y z mi da le coordinate della decomposizione MPI
      !così so su quale nodo vengono lavorate le coordinate x y z
      allocate(xinidom(0:proc_x-1))
      allocate(xfindom(0:proc_x-1))
      xinidom(:)=0
      xfindom(:)=-1
      xinidom(0)=1
      xfindom(0)=nx+xinidom(0)-1
      do i=1,proc_x-1
        xinidom(i)=xfindom(i-1)+1
        xfindom(i)=xinidom(i)+nx-1
      enddo
      
      allocate(yinidom(0:proc_y-1))
      allocate(yfindom(0:proc_y-1))
      yinidom(:)=0
      yfindom(:)=-1
      yinidom(0)=1
      yfindom(0)=ny+yinidom(0)-1
      do i=1,proc_y-1
        yinidom(i)=yfindom(i-1)+1
        yfindom(i)=yinidom(i)+ny-1
      enddo
      
      allocate(zinidom(0:proc_z-1))
      allocate(zfindom(0:proc_z-1))
      zinidom(:)=0
      zfindom(:)=-1
      zinidom(0)=1
      zfindom(0)=nz+zinidom(0)-1
      do i=1,proc_z-1
        zinidom(i)=zfindom(i-1)+1
        zfindom(i)=zinidom(i)+nz-1
      enddo
     
   
      call flush(6) 
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      do idrank=0,nprocs
          if(idrank==myrank)then  !
            write(6,'(a,i4,a,9i8)')'DEC myrank ',myrank,' coords ',coords,&
             start_idx(1),end_idx(1),&
             start_idx(2),end_idx(2),&
             start_idx(3),end_idx(3)
            call flush(6)
          endif
#ifdef MPI
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      enddo
      call flush(6)
#ifdef VERBOSE 

#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

      do idrank=0,nprocs-1
        do l=1,nlinksmpi
          if(idrank==myrank)then
            write(6,'(a,i4,a,4i3,a,i4,a,i4,2l2,a,3l2)')'Myrank ',myrank,' l ',l,exmpi(l),eympi(l),ezmpi(l),&
             ' source ',source_dir(l),' dest ',dest_dir(l),lsource_dir(l),ldest_dir(l),' intpbc ',intpbc_dir(1:3,l)
            call flush(6)
          endif
#ifdef MPI
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
        enddo
      enddo
      call flush(6)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

#endif
        !trovo per ogni direzione l quali popolazioni devono essere inviate e le storo in links_faces
        !occhio che devo vedere quali popolazioni mandare ma il lattice può essere
        !tipo d3q15 o d3q19 e quindi non devo mandare nulla
        !exmpi eympi ezmpi sono le 26 direzioni di MPI 
        !ex ey ez sono le direzioni del lattice d3q15 o d3q19 o d3q27 
        allocate(num_links_pops(1:nlinksmpi))
        !faces
        nfaces=6
        nlinks_faces=9
        allocate(links_faces(1:nlinks_faces,1:6))
        do l=1,6
          nlinks_faces=0
          do ll=1,nlinks
            if((abs(exmpi(l))==1 .and. exmpi(l)==ex(ll)) .or. &
             (abs(eympi(l))==1 .and. eympi(l)==ey(ll)) .or. &
             (abs(ezmpi(l))==1 .and. ezmpi(l)==ez(ll)))then
               nlinks_faces=nlinks_faces+1
               links_faces(nlinks_faces,l)=ll
#ifdef VERBOSE  
              if(myrank==0)write(6,'(a,i3,a,3i3,a,a2,a,3i3)')'dir l ',l,' disp ',&
               exmpi(l),eympi(l),ezmpi(l),&
               ' f',write_fmtnumb2(ll),' dir ',ex(ll),ey(ll),ez(ll)
              call flush(6)
#endif
            endif
          enddo
          num_links_pops(l)=nlinks_faces
        enddo
        !edges
        nedges=12
        nlinks_edges=3
        allocate(links_edges(1:nlinks_edges,7:18))
        do l=7,18
          nlinks_edges=0
          do ll=1,nlinks
            if((abs(exmpi(l))==1 .and. abs(eympi(l))==1 .and. exmpi(l)==ex(ll) .and. eympi(l)==ey(ll)) .or. &
             (abs(exmpi(l))==1 .and. abs(ezmpi(l))==1 .and. exmpi(l)==ex(ll) .and. ezmpi(l)==ez(ll)) .or. &
             (abs(eympi(l))==1 .and. abs(ezmpi(l))==1 .and. eympi(l)==ey(ll) .and. ezmpi(l)==ez(ll)))then
               nlinks_edges=nlinks_edges+1
               links_edges(nlinks_edges,l)=ll
#ifdef VERBOSE  
              if(myrank==0)write(6,'(a,i3,a,3i3,a,a2,a,3i3)')'dir l ',l,' disp ',&
               exmpi(l),eympi(l),ezmpi(l),&
               ' f',write_fmtnumb2(ll),' dir ',ex(ll),ey(ll),ez(ll)
              call flush(6)
#endif
            endif
          enddo
          num_links_pops(l)=nlinks_edges
        enddo
        !corner
        ncorners=8
        nlinks_corners=1
        allocate(links_corners(1:nlinks_corners,19:nlinksmpi))
        do l=19,nlinksmpi
          nlinks_corners=0
          do ll=1,nlinks
            if(exmpi(l)==ex(ll) .and. eympi(l)==ey(ll) .and. ezmpi(l)==ez(ll))then
               nlinks_corners=nlinks_corners+1
               links_corners(nlinks_corners,l)=ll
#ifdef VERBOSE  
              if(myrank==0)write(6,'(a,i3,a,3i3,a,a2,a,3i3)')'dir l ',l,' disp ',&
               exmpi(l),eympi(l),ezmpi(l),&
               ' f',write_fmtnumb2(ll),' dir ',ex(ll),ey(ll),ez(ll)
              call flush(6)
#endif
            endif
          enddo
          num_links_pops(l)=nlinks_corners
        enddo
      
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      if(myrank==0)then
        write(6,'(a,3i3)')'dir nlinks ',nlinks_faces,nlinks_edges,nlinks_corners
      endif

      !riempio la lista links_pops con le popolazioni da inviare per ogni direzione l
      !il primo indice di lista links_pops è preso dal massimo delle pops da mandare tra faces edges e corners (ovviamente è sempre faces)
      !nota che num_links_pops è il numero di poplazioni da inviare per direzione l
      nlinks_max=max(nlinks_faces,nlinks_edges,nlinks_corners)
      allocate(links_pops(1:nlinks_max,1:nlinksmpi))
      do l=1,6
        do ll=1,num_links_pops(l)
          links_pops(ll,l)=links_faces(ll,l)
        enddo
      enddo
      do l=7,18
        do ll=1,num_links_pops(l)
          links_pops(ll,l)=links_edges(ll,l)
        enddo
      enddo
      do l=19,nlinksmpi
        do ll=1,num_links_pops(l)
          links_pops(ll,l)=links_corners(ll,l)
        enddo
      enddo
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      !solo se ci sono popolazioni da mandare 
      !allora metti ldestpop_dir e lsourcepop_dir true
      !con d3q19 o d3q15 spesso non devo mandare nulla 
      do l=1,nlinksmpi
        ldestpop_dir(l)=(ldest_dir(l) .and. num_links_pops(l)>0)
        lsourcepop_dir(l)=(lsource_dir(l) .and. num_links_pops(l)>0)
      enddo
      !solo se ci sono popolazioni da fare pbc interno lo fai 
      do l=1,nlinksmpi
        lintpbcpop_dir(l)=(lintpbc_dir(l) .and. num_links_pops(l)>0)
      enddo
      
      !mi storo gli estremi i j k che devono essere inviati e ricevuti lungo ogni direzione l
      allocate(dest_extr(6,nlinksmpi))
      allocate(source_extr(6,nlinksmpi))
      allocate(f_dest_extr(6,nlinksmpi))
      allocate(f_source_extr(6,nlinksmpi))
      allocate(b_dest_extr(6,nlinksmpi))
      allocate(b_source_extr(6,nlinksmpi))
      
      !faces
      do l=1,6
        if(exmpi(l)==1)then
          dest_extr(1,l)=nx+1
          dest_extr(2,l)=nx+1
          dest_extr(3,l)=1
          dest_extr(4,l)=ny
          dest_extr(5,l)=1
          dest_extr(6,l)=nz
          
          source_extr(1,l)=1
          source_extr(2,l)=1
          source_extr(3,l)=1
          source_extr(4,l)=ny
          source_extr(5,l)=1
          source_extr(6,l)=nz
          
          f_dest_extr(1,l)=nx
          f_dest_extr(2,l)=nx
          f_dest_extr(3,l)=1
          f_dest_extr(4,l)=ny
          f_dest_extr(5,l)=1
          f_dest_extr(6,l)=nz
          
          f_source_extr(1,l)=0
          f_source_extr(2,l)=0
          f_source_extr(3,l)=1
          f_source_extr(4,l)=ny
          f_source_extr(5,l)=1
          f_source_extr(6,l)=nz
          
          b_dest_extr(1,l)=nx-nbuff+1
          b_dest_extr(2,l)=nx
          b_dest_extr(3,l)=1
          b_dest_extr(4,l)=ny
          b_dest_extr(5,l)=1
          b_dest_extr(6,l)=nz
          
          b_source_extr(1,l)=1-nbuff
          b_source_extr(2,l)=0
          b_source_extr(3,l)=1
          b_source_extr(4,l)=ny
          b_source_extr(5,l)=1
          b_source_extr(6,l)=nz
        endif
        if(exmpi(l)==-1)then
          dest_extr(1,l)=0
          dest_extr(2,l)=0
          dest_extr(3,l)=1
          dest_extr(4,l)=ny
          dest_extr(5,l)=1
          dest_extr(6,l)=nz
          
          source_extr(1,l)=nx
          source_extr(2,l)=nx
          source_extr(3,l)=1
          source_extr(4,l)=ny
          source_extr(5,l)=1
          source_extr(6,l)=nz
          
          f_dest_extr(1,l)=1
          f_dest_extr(2,l)=1
          f_dest_extr(3,l)=1
          f_dest_extr(4,l)=ny
          f_dest_extr(5,l)=1
          f_dest_extr(6,l)=nz
          
          f_source_extr(1,l)=nx+1
          f_source_extr(2,l)=nx+1
          f_source_extr(3,l)=1
          f_source_extr(4,l)=ny
          f_source_extr(5,l)=1
          f_source_extr(6,l)=nz
          
          b_dest_extr(1,l)=1
          b_dest_extr(2,l)=nbuff
          b_dest_extr(3,l)=1
          b_dest_extr(4,l)=ny
          b_dest_extr(5,l)=1
          b_dest_extr(6,l)=nz
          
          b_source_extr(1,l)=nx+1
          b_source_extr(2,l)=nx+nbuff
          b_source_extr(3,l)=1
          b_source_extr(4,l)=ny
          b_source_extr(5,l)=1
          b_source_extr(6,l)=nz
        endif
        if(eympi(l)==1)then
          dest_extr(1,l)=1
          dest_extr(2,l)=nx
          dest_extr(3,l)=ny+1
          dest_extr(4,l)=ny+1
          dest_extr(5,l)=1
          dest_extr(6,l)=nz
          
          source_extr(1,l)=1
          source_extr(2,l)=nx
          source_extr(3,l)=1
          source_extr(4,l)=1
          source_extr(5,l)=1
          source_extr(6,l)=nz
          
          f_dest_extr(1,l)=1
          f_dest_extr(2,l)=nx
          f_dest_extr(3,l)=ny
          f_dest_extr(4,l)=ny
          f_dest_extr(5,l)=1
          f_dest_extr(6,l)=nz
          
          f_source_extr(1,l)=1
          f_source_extr(2,l)=nx
          f_source_extr(3,l)=0
          f_source_extr(4,l)=0
          f_source_extr(5,l)=1
          f_source_extr(6,l)=nz
          
          b_dest_extr(1,l)=1
          b_dest_extr(2,l)=nx
          b_dest_extr(3,l)=ny-nbuff+1
          b_dest_extr(4,l)=ny
          b_dest_extr(5,l)=1
          b_dest_extr(6,l)=nz
          
          b_source_extr(1,l)=1
          b_source_extr(2,l)=nx
          b_source_extr(3,l)=1-nbuff
          b_source_extr(4,l)=0
          b_source_extr(5,l)=1
          b_source_extr(6,l)=nz
        endif
        if(eympi(l)==-1)then
          dest_extr(1,l)=1
          dest_extr(2,l)=nx
          dest_extr(3,l)=0
          dest_extr(4,l)=0
          dest_extr(5,l)=1
          dest_extr(6,l)=nz
          
          source_extr(1,l)=1
          source_extr(2,l)=nx
          source_extr(3,l)=ny
          source_extr(4,l)=ny
          source_extr(5,l)=1
          source_extr(6,l)=nz
          
          f_dest_extr(1,l)=1
          f_dest_extr(2,l)=nx
          f_dest_extr(3,l)=1
          f_dest_extr(4,l)=1
          f_dest_extr(5,l)=1
          f_dest_extr(6,l)=nz
          
          f_source_extr(1,l)=1
          f_source_extr(2,l)=nx
          f_source_extr(3,l)=ny+1
          f_source_extr(4,l)=ny+1
          f_source_extr(5,l)=1
          f_source_extr(6,l)=nz
          
          b_dest_extr(1,l)=1
          b_dest_extr(2,l)=nx
          b_dest_extr(3,l)=1
          b_dest_extr(4,l)=nbuff
          b_dest_extr(5,l)=1
          b_dest_extr(6,l)=nz
          
          b_source_extr(1,l)=1
          b_source_extr(2,l)=nx
          b_source_extr(3,l)=ny+1
          b_source_extr(4,l)=ny+nbuff
          b_source_extr(5,l)=1
          b_source_extr(6,l)=nz
        endif
        if(ezmpi(l)==1)then
          dest_extr(1,l)=1
          dest_extr(2,l)=nx
          dest_extr(3,l)=1
          dest_extr(4,l)=ny
          dest_extr(5,l)=nz+1
          dest_extr(6,l)=nz+1
          
          source_extr(1,l)=1
          source_extr(2,l)=nx
          source_extr(3,l)=1
          source_extr(4,l)=ny
          source_extr(5,l)=1
          source_extr(6,l)=1
          
          f_dest_extr(1,l)=1
          f_dest_extr(2,l)=nx
          f_dest_extr(3,l)=1
          f_dest_extr(4,l)=ny
          f_dest_extr(5,l)=nz
          f_dest_extr(6,l)=nz
          
          f_source_extr(1,l)=1
          f_source_extr(2,l)=nx
          f_source_extr(3,l)=1
          f_source_extr(4,l)=ny
          f_source_extr(5,l)=0
          f_source_extr(6,l)=0
          
          b_dest_extr(1,l)=1
          b_dest_extr(2,l)=nx
          b_dest_extr(3,l)=1
          b_dest_extr(4,l)=ny
          b_dest_extr(5,l)=nz-nbuff+1
          b_dest_extr(6,l)=nz
          
          b_source_extr(1,l)=1
          b_source_extr(2,l)=nx
          b_source_extr(3,l)=1
          b_source_extr(4,l)=ny
          b_source_extr(5,l)=1-nbuff
          b_source_extr(6,l)=0
        endif
        if(ezmpi(l)==-1)then
          dest_extr(1,l)=1
          dest_extr(2,l)=nx
          dest_extr(3,l)=1
          dest_extr(4,l)=ny
          dest_extr(5,l)=0
          dest_extr(6,l)=0
          
          source_extr(1,l)=1
          source_extr(2,l)=nx
          source_extr(3,l)=1
          source_extr(4,l)=ny
          source_extr(5,l)=nz
          source_extr(6,l)=nz
          
          f_dest_extr(1,l)=1
          f_dest_extr(2,l)=nx
          f_dest_extr(3,l)=1
          f_dest_extr(4,l)=ny
          f_dest_extr(5,l)=1
          f_dest_extr(6,l)=1
          
          f_source_extr(1,l)=1
          f_source_extr(2,l)=nx
          f_source_extr(3,l)=1
          f_source_extr(4,l)=ny
          f_source_extr(5,l)=nz+1
          f_source_extr(6,l)=nz+1
          
          b_dest_extr(1,l)=1
          b_dest_extr(2,l)=nx
          b_dest_extr(3,l)=1
          b_dest_extr(4,l)=ny
          b_dest_extr(5,l)=1
          b_dest_extr(6,l)=nbuff
          
          b_source_extr(1,l)=1
          b_source_extr(2,l)=nx
          b_source_extr(3,l)=1
          b_source_extr(4,l)=ny
          b_source_extr(5,l)=nz+1
          b_source_extr(6,l)=nz+nbuff
        endif
      enddo
      !edges
      do l=7,18
        !!!!   x   y
        if(exmpi(l)==1 .and. eympi(l)==1)then
          dest_extr(1,l)=nx+1
          dest_extr(2,l)=nx+1
          dest_extr(3,l)=ny+1
          dest_extr(4,l)=ny+1
          dest_extr(5,l)=1
          dest_extr(6,l)=nz
          
          source_extr(1,l)=1
          source_extr(2,l)=1
          source_extr(3,l)=1
          source_extr(4,l)=1
          source_extr(5,l)=1
          source_extr(6,l)=nz
          
          f_dest_extr(1,l)=nx
          f_dest_extr(2,l)=nx
          f_dest_extr(3,l)=ny
          f_dest_extr(4,l)=ny
          f_dest_extr(5,l)=1
          f_dest_extr(6,l)=nz
          
          f_source_extr(1,l)=0
          f_source_extr(2,l)=0
          f_source_extr(3,l)=0
          f_source_extr(4,l)=0
          f_source_extr(5,l)=1
          f_source_extr(6,l)=nz
          
          b_dest_extr(1,l)=nx-nbuff+1
          b_dest_extr(2,l)=nx
          b_dest_extr(3,l)=ny-nbuff+1
          b_dest_extr(4,l)=ny
          b_dest_extr(5,l)=1
          b_dest_extr(6,l)=nz
          
          b_source_extr(1,l)=1-nbuff
          b_source_extr(2,l)=0
          b_source_extr(3,l)=1-nbuff
          b_source_extr(4,l)=0
          b_source_extr(5,l)=1
          b_source_extr(6,l)=nz
        endif
        if(exmpi(l)==-1 .and. eympi(l)==-1)then
          dest_extr(1,l)=0
          dest_extr(2,l)=0
          dest_extr(3,l)=0
          dest_extr(4,l)=0
          dest_extr(5,l)=1
          dest_extr(6,l)=nz
          
          source_extr(1,l)=nx
          source_extr(2,l)=nx
          source_extr(3,l)=ny
          source_extr(4,l)=ny
          source_extr(5,l)=1
          source_extr(6,l)=nz
          
          f_dest_extr(1,l)=1
          f_dest_extr(2,l)=1
          f_dest_extr(3,l)=1
          f_dest_extr(4,l)=1
          f_dest_extr(5,l)=1
          f_dest_extr(6,l)=nz
          
          f_source_extr(1,l)=nx+1
          f_source_extr(2,l)=nx+1
          f_source_extr(3,l)=ny+1
          f_source_extr(4,l)=ny+1
          f_source_extr(5,l)=1
          f_source_extr(6,l)=nz
          
          b_dest_extr(1,l)=1
          b_dest_extr(2,l)=nbuff
          b_dest_extr(3,l)=1
          b_dest_extr(4,l)=nbuff
          b_dest_extr(5,l)=1
          b_dest_extr(6,l)=nz
          
          b_source_extr(1,l)=nx+1
          b_source_extr(2,l)=nx+nbuff
          b_source_extr(3,l)=ny+1
          b_source_extr(4,l)=ny+nbuff
          b_source_extr(5,l)=1
          b_source_extr(6,l)=nz
        endif
        if(exmpi(l)==1 .and. eympi(l)==-1)then
          dest_extr(1,l)=nx+1
          dest_extr(2,l)=nx+1
          dest_extr(3,l)=0
          dest_extr(4,l)=0
          dest_extr(5,l)=1
          dest_extr(6,l)=nz
          
          source_extr(1,l)=1
          source_extr(2,l)=1
          source_extr(3,l)=ny
          source_extr(4,l)=ny
          source_extr(5,l)=1
          source_extr(6,l)=nz
          
          f_dest_extr(1,l)=nx
          f_dest_extr(2,l)=nx
          f_dest_extr(3,l)=1
          f_dest_extr(4,l)=1
          f_dest_extr(5,l)=1
          f_dest_extr(6,l)=nz
          
          f_source_extr(1,l)=0
          f_source_extr(2,l)=0
          f_source_extr(3,l)=ny+1
          f_source_extr(4,l)=ny+1
          f_source_extr(5,l)=1
          f_source_extr(6,l)=nz
          
          b_dest_extr(1,l)=nx-nbuff+1
          b_dest_extr(2,l)=nx
          b_dest_extr(3,l)=1
          b_dest_extr(4,l)=nbuff
          b_dest_extr(5,l)=1
          b_dest_extr(6,l)=nz
          
          b_source_extr(1,l)=1-nbuff
          b_source_extr(2,l)=0
          b_source_extr(3,l)=ny+1
          b_source_extr(4,l)=ny+nbuff
          b_source_extr(5,l)=1
          b_source_extr(6,l)=nz
        endif
        if(exmpi(l)==-1 .and. eympi(l)==1)then
          dest_extr(1,l)=0
          dest_extr(2,l)=0
          dest_extr(3,l)=ny+1
          dest_extr(4,l)=ny+1
          dest_extr(5,l)=1
          dest_extr(6,l)=nz
          
          source_extr(1,l)=nx
          source_extr(2,l)=nx
          source_extr(3,l)=1
          source_extr(4,l)=1
          source_extr(5,l)=1
          source_extr(6,l)=nz
          
          f_dest_extr(1,l)=1
          f_dest_extr(2,l)=1
          f_dest_extr(3,l)=ny
          f_dest_extr(4,l)=ny
          f_dest_extr(5,l)=1
          f_dest_extr(6,l)=nz
          
          f_source_extr(1,l)=nx+1
          f_source_extr(2,l)=nx+1
          f_source_extr(3,l)=0
          f_source_extr(4,l)=0
          f_source_extr(5,l)=1
          f_source_extr(6,l)=nz
          
          b_dest_extr(1,l)=1
          b_dest_extr(2,l)=nbuff
          b_dest_extr(3,l)=ny-nbuff+1
          b_dest_extr(4,l)=ny
          b_dest_extr(5,l)=1
          b_dest_extr(6,l)=nz
          
          b_source_extr(1,l)=nx+1
          b_source_extr(2,l)=nx+nbuff
          b_source_extr(3,l)=1-nbuff
          b_source_extr(4,l)=0
          b_source_extr(5,l)=1
          b_source_extr(6,l)=nz
        endif
        !!!!   x   z
        if(exmpi(l)==1 .and. ezmpi(l)==1)then
          dest_extr(1,l)=nx+1
          dest_extr(2,l)=nx+1
          dest_extr(3,l)=1
          dest_extr(4,l)=ny
          dest_extr(5,l)=nz+1
          dest_extr(6,l)=nz+1
          
          source_extr(1,l)=1
          source_extr(2,l)=1
          source_extr(3,l)=1
          source_extr(4,l)=ny
          source_extr(5,l)=1
          source_extr(6,l)=1
          
          f_dest_extr(1,l)=nx
          f_dest_extr(2,l)=nx
          f_dest_extr(3,l)=1
          f_dest_extr(4,l)=ny
          f_dest_extr(5,l)=nz
          f_dest_extr(6,l)=nz
          
          f_source_extr(1,l)=0
          f_source_extr(2,l)=0
          f_source_extr(3,l)=1
          f_source_extr(4,l)=ny
          f_source_extr(5,l)=0
          f_source_extr(6,l)=0
          
          b_dest_extr(1,l)=nx-nbuff+1
          b_dest_extr(2,l)=nx
          b_dest_extr(3,l)=1
          b_dest_extr(4,l)=ny
          b_dest_extr(5,l)=nz-nbuff+1
          b_dest_extr(6,l)=nz
          
          b_source_extr(1,l)=1-nbuff
          b_source_extr(2,l)=0
          b_source_extr(3,l)=1
          b_source_extr(4,l)=ny
          b_source_extr(5,l)=1-nbuff
          b_source_extr(6,l)=0
        endif
        if(exmpi(l)==-1 .and. ezmpi(l)==-1)then
          dest_extr(1,l)=0
          dest_extr(2,l)=0
          dest_extr(3,l)=1
          dest_extr(4,l)=ny
          dest_extr(5,l)=0
          dest_extr(6,l)=0
          
          source_extr(1,l)=nx
          source_extr(2,l)=nx
          source_extr(3,l)=1
          source_extr(4,l)=ny
          source_extr(5,l)=nz
          source_extr(6,l)=nz
          
          f_dest_extr(1,l)=1
          f_dest_extr(2,l)=1
          f_dest_extr(3,l)=1
          f_dest_extr(4,l)=ny
          f_dest_extr(5,l)=1
          f_dest_extr(6,l)=1
          
          f_source_extr(1,l)=nx+1
          f_source_extr(2,l)=nx+1
          f_source_extr(3,l)=1
          f_source_extr(4,l)=ny
          f_source_extr(5,l)=nz+1
          f_source_extr(6,l)=nz+1
          
          b_dest_extr(1,l)=1
          b_dest_extr(2,l)=nbuff
          b_dest_extr(3,l)=1
          b_dest_extr(4,l)=ny
          b_dest_extr(5,l)=1
          b_dest_extr(6,l)=nbuff
          
          b_source_extr(1,l)=nx+1
          b_source_extr(2,l)=nx+nbuff
          b_source_extr(3,l)=1
          b_source_extr(4,l)=ny
          b_source_extr(5,l)=nz+1
          b_source_extr(6,l)=nz+nbuff
        endif
        if(exmpi(l)==1 .and. ezmpi(l)==-1)then
          dest_extr(1,l)=nx+1
          dest_extr(2,l)=nx+1
          dest_extr(3,l)=1
          dest_extr(4,l)=ny
          dest_extr(5,l)=0
          dest_extr(6,l)=0
          
          source_extr(1,l)=1
          source_extr(2,l)=1
          source_extr(3,l)=1
          source_extr(4,l)=ny
          source_extr(5,l)=nz
          source_extr(6,l)=nz
          
          f_dest_extr(1,l)=nx
          f_dest_extr(2,l)=nx
          f_dest_extr(3,l)=1
          f_dest_extr(4,l)=ny
          f_dest_extr(5,l)=1
          f_dest_extr(6,l)=1
          
          f_source_extr(1,l)=0
          f_source_extr(2,l)=0
          f_source_extr(3,l)=1
          f_source_extr(4,l)=ny
          f_source_extr(5,l)=nz+1
          f_source_extr(6,l)=nz+1
          
          b_dest_extr(1,l)=nx-nbuff+1
          b_dest_extr(2,l)=nx
          b_dest_extr(3,l)=1
          b_dest_extr(4,l)=ny
          b_dest_extr(5,l)=1
          b_dest_extr(6,l)=nbuff
          
          b_source_extr(1,l)=1-nbuff
          b_source_extr(2,l)=0
          b_source_extr(3,l)=1
          b_source_extr(4,l)=ny
          b_source_extr(5,l)=nz+1
          b_source_extr(6,l)=nz+nbuff
        endif
        if(exmpi(l)==-1 .and. ezmpi(l)==1)then
          dest_extr(1,l)=0
          dest_extr(2,l)=0
          dest_extr(3,l)=1
          dest_extr(4,l)=ny
          dest_extr(5,l)=nz+1
          dest_extr(6,l)=nz+1
          
          source_extr(1,l)=nx
          source_extr(2,l)=nx
          source_extr(3,l)=1
          source_extr(4,l)=ny
          source_extr(5,l)=1
          source_extr(6,l)=1
          
          f_dest_extr(1,l)=1
          f_dest_extr(2,l)=1
          f_dest_extr(3,l)=1
          f_dest_extr(4,l)=ny
          f_dest_extr(5,l)=nz
          f_dest_extr(6,l)=nz
          
          f_source_extr(1,l)=nx+1
          f_source_extr(2,l)=nx+1
          f_source_extr(3,l)=1
          f_source_extr(4,l)=ny
          f_source_extr(5,l)=0
          f_source_extr(6,l)=0
          
          b_dest_extr(1,l)=1
          b_dest_extr(2,l)=nbuff
          b_dest_extr(3,l)=1
          b_dest_extr(4,l)=ny
          b_dest_extr(5,l)=nz-nbuff+1
          b_dest_extr(6,l)=nz
          
          b_source_extr(1,l)=nx+1
          b_source_extr(2,l)=nx+nbuff
          b_source_extr(3,l)=1
          b_source_extr(4,l)=ny
          b_source_extr(5,l)=1-nbuff
          b_source_extr(6,l)=0
        endif
        !!!!   y   z
        if(eympi(l)==1 .and. ezmpi(l)==1)then
          dest_extr(1,l)=1
          dest_extr(2,l)=nx
          dest_extr(3,l)=ny+1
          dest_extr(4,l)=ny+1
          dest_extr(5,l)=nz+1
          dest_extr(6,l)=nz+1
          
          source_extr(1,l)=1
          source_extr(2,l)=nx
          source_extr(3,l)=1
          source_extr(4,l)=1
          source_extr(5,l)=1
          source_extr(6,l)=1
          
          f_dest_extr(1,l)=1
          f_dest_extr(2,l)=nx
          f_dest_extr(3,l)=ny
          f_dest_extr(4,l)=ny
          f_dest_extr(5,l)=nz
          f_dest_extr(6,l)=nz
          
          f_source_extr(1,l)=1
          f_source_extr(2,l)=nx
          f_source_extr(3,l)=0
          f_source_extr(4,l)=0
          f_source_extr(5,l)=0
          f_source_extr(6,l)=0
          
          b_dest_extr(1,l)=1
          b_dest_extr(2,l)=nx
          b_dest_extr(3,l)=ny-nbuff+1
          b_dest_extr(4,l)=ny
          b_dest_extr(5,l)=nz-nbuff+1
          b_dest_extr(6,l)=nz
          
          b_source_extr(1,l)=1
          b_source_extr(2,l)=nx
          b_source_extr(3,l)=1-nbuff
          b_source_extr(4,l)=0
          b_source_extr(5,l)=1-nbuff
          b_source_extr(6,l)=0
        endif
        if(eympi(l)==-1 .and. ezmpi(l)==-1)then
          dest_extr(1,l)=1
          dest_extr(2,l)=nx
          dest_extr(3,l)=0
          dest_extr(4,l)=0
          dest_extr(5,l)=0
          dest_extr(6,l)=0
          
          source_extr(1,l)=1
          source_extr(2,l)=nx
          source_extr(3,l)=ny
          source_extr(4,l)=ny
          source_extr(5,l)=nz
          source_extr(6,l)=nz
          
          f_dest_extr(1,l)=1
          f_dest_extr(2,l)=nx
          f_dest_extr(3,l)=1
          f_dest_extr(4,l)=1
          f_dest_extr(5,l)=1
          f_dest_extr(6,l)=1
          
          f_source_extr(1,l)=1
          f_source_extr(2,l)=nx
          f_source_extr(3,l)=ny+1
          f_source_extr(4,l)=ny+1
          f_source_extr(5,l)=nz+1
          f_source_extr(6,l)=nz+1
          
          b_dest_extr(1,l)=1
          b_dest_extr(2,l)=nx
          b_dest_extr(3,l)=1
          b_dest_extr(4,l)=nbuff
          b_dest_extr(5,l)=1
          b_dest_extr(6,l)=nbuff
          
          b_source_extr(1,l)=1
          b_source_extr(2,l)=nx
          b_source_extr(3,l)=ny+1
          b_source_extr(4,l)=ny+nbuff
          b_source_extr(5,l)=nz+1
          b_source_extr(6,l)=nz+nbuff
        endif
        if(eympi(l)==1 .and. ezmpi(l)==-1)then
          dest_extr(1,l)=1
          dest_extr(2,l)=nx
          dest_extr(3,l)=ny+1
          dest_extr(4,l)=ny+1
          dest_extr(5,l)=0
          dest_extr(6,l)=0
          
          source_extr(1,l)=1
          source_extr(2,l)=nx
          source_extr(3,l)=1
          source_extr(4,l)=1
          source_extr(5,l)=nz
          source_extr(6,l)=nz
          
          f_dest_extr(1,l)=1
          f_dest_extr(2,l)=nx
          f_dest_extr(3,l)=ny
          f_dest_extr(4,l)=ny
          f_dest_extr(5,l)=1
          f_dest_extr(6,l)=1
          
          f_source_extr(1,l)=1
          f_source_extr(2,l)=nx
          f_source_extr(3,l)=0
          f_source_extr(4,l)=0
          f_source_extr(5,l)=nz+1
          f_source_extr(6,l)=nz+1
          
          b_dest_extr(1,l)=1
          b_dest_extr(2,l)=nx
          b_dest_extr(3,l)=ny-nbuff+1
          b_dest_extr(4,l)=ny
          b_dest_extr(5,l)=1
          b_dest_extr(6,l)=nbuff
          
          b_source_extr(1,l)=1
          b_source_extr(2,l)=nx
          b_source_extr(3,l)=1-nbuff
          b_source_extr(4,l)=0
          b_source_extr(5,l)=nz+1
          b_source_extr(6,l)=nz+nbuff
        endif
        if(eympi(l)==-1 .and. ezmpi(l)==1)then
          dest_extr(1,l)=1
          dest_extr(2,l)=nx
          dest_extr(3,l)=0
          dest_extr(4,l)=0
          dest_extr(5,l)=nz+1
          dest_extr(6,l)=nz+1
          
          source_extr(1,l)=1
          source_extr(2,l)=nx
          source_extr(3,l)=ny
          source_extr(4,l)=ny
          source_extr(5,l)=1
          source_extr(6,l)=1
          
          f_dest_extr(1,l)=1
          f_dest_extr(2,l)=nx
          f_dest_extr(3,l)=1
          f_dest_extr(4,l)=1
          f_dest_extr(5,l)=nz
          f_dest_extr(6,l)=nz
          
          f_source_extr(1,l)=1
          f_source_extr(2,l)=nx
          f_source_extr(3,l)=ny+1
          f_source_extr(4,l)=ny+1
          f_source_extr(5,l)=0
          f_source_extr(6,l)=0
          
          b_dest_extr(1,l)=1
          b_dest_extr(2,l)=nx
          b_dest_extr(3,l)=1
          b_dest_extr(4,l)=nbuff
          b_dest_extr(5,l)=nz-nbuff+1
          b_dest_extr(6,l)=nz
          
          b_source_extr(1,l)=1
          b_source_extr(2,l)=nx
          b_source_extr(3,l)=ny+1
          b_source_extr(4,l)=ny+nbuff
          b_source_extr(5,l)=1-nbuff
          b_source_extr(6,l)=0
        endif
        
      enddo
      !corner
      do l=19,nlinksmpi
        if(exmpi(l)==1)then
          dest_extr(1,l)=nx+1
          dest_extr(2,l)=nx+1
          
          source_extr(1,l)=1
          source_extr(2,l)=1
          
          f_dest_extr(1,l)=nx
          f_dest_extr(2,l)=nx
          
          f_source_extr(1,l)=0
          f_source_extr(2,l)=0
          
          b_dest_extr(1,l)=nx-nbuff+1
          b_dest_extr(2,l)=nx
          
          b_source_extr(1,l)=1-nbuff
          b_source_extr(2,l)=0
        else
          dest_extr(1,l)=0
          dest_extr(2,l)=0
          
          source_extr(1,l)=nx
          source_extr(2,l)=nx
          
          f_dest_extr(1,l)=1
          f_dest_extr(2,l)=1
          
          f_source_extr(1,l)=nx+1
          f_source_extr(2,l)=nx+1
          
          b_dest_extr(1,l)=1
          b_dest_extr(2,l)=nbuff
          
          b_source_extr(1,l)=nx+1
          b_source_extr(2,l)=nx+nbuff
        endif
        if(eympi(l)==1)then
          dest_extr(3,l)=ny+1
          dest_extr(4,l)=ny+1
          
          source_extr(3,l)=1
          source_extr(4,l)=1
          
          f_dest_extr(3,l)=ny
          f_dest_extr(4,l)=ny
          
          f_source_extr(3,l)=0
          f_source_extr(4,l)=0
          
          b_dest_extr(3,l)=ny-nbuff+1
          b_dest_extr(4,l)=ny
          
          b_source_extr(3,l)=1-nbuff
          b_source_extr(4,l)=0
        else
          dest_extr(3,l)=0
          dest_extr(4,l)=0
          
          source_extr(3,l)=ny
          source_extr(4,l)=ny
          
          f_dest_extr(3,l)=1
          f_dest_extr(4,l)=1
          
          f_source_extr(3,l)=ny+1
          f_source_extr(4,l)=ny+1
          
          b_dest_extr(3,l)=1
          b_dest_extr(4,l)=nbuff
          
          b_source_extr(3,l)=ny+1
          b_source_extr(4,l)=ny+nbuff
        endif
        if(ezmpi(l)==1)then
          dest_extr(5,l)=nz+1
          dest_extr(6,l)=nz+1
          
          source_extr(5,l)=1
          source_extr(6,l)=1
          
          f_dest_extr(5,l)=nz
          f_dest_extr(6,l)=nz
          
          f_source_extr(5,l)=0
          f_source_extr(6,l)=0
          
          b_dest_extr(5,l)=nz-nbuff+1
          b_dest_extr(6,l)=nz
          
          b_source_extr(5,l)=1-nbuff
          b_source_extr(6,l)=0
        else
          dest_extr(5,l)=0
          dest_extr(6,l)=0
          
          source_extr(5,l)=nz
          source_extr(6,l)=nz
          
          f_dest_extr(5,l)=1
          f_dest_extr(6,l)=1
          
          f_source_extr(5,l)=nz+1
          f_source_extr(6,l)=nz+1
          
          b_dest_extr(5,l)=1
          b_dest_extr(6,l)=nbuff
          
          b_source_extr(5,l)=nz+1
          b_source_extr(6,l)=nz+nbuff
        endif
      enddo
      !calcolo le quantita complessive da movimentare per ogni direzione l
       do l=1,nlinksmpi
         num_extr(l)=(source_extr(2,l)-source_extr(1,l)+1)* &
          (source_extr(4,l)-source_extr(3,l)+1)* &
          (source_extr(6,l)-source_extr(5,l)+1)*num_links_pops(l)
         i_num_extr(l)=(f_source_extr(2,l)-f_source_extr(1,l)+1)* &
          (f_source_extr(4,l)-f_source_extr(3,l)+1)* &
          (f_source_extr(6,l)-f_source_extr(5,l)+1)
         f_num_extr(l)=(f_source_extr(2,l)-f_source_extr(1,l)+1)* &
          (f_source_extr(4,l)-f_source_extr(3,l)+1)* &
          (f_source_extr(6,l)-f_source_extr(5,l)+1)*num_f_datampi
         b_num_extr(l)=(b_source_extr(2,l)-b_source_extr(1,l)+1)* &
          (b_source_extr(4,l)-b_source_extr(3,l)+1)* &
          (b_source_extr(6,l)-b_source_extr(5,l)+1)*num_b_datampi
       enddo
      
        
      
#ifdef VERBOSE 
       !stampo per debug 
       if(myrank==0)write(6,'(a)')'#######################   dest_extr    source_extr #######################'    
       do l=1,nlinksmpi
         if(myrank==0)write(6,'(a,i3,a,3i3,a,6i4,a,6i4,a,i4)')'dir l ',l,' disp ',&
             exmpi(l),eympi(l),ezmpi(l),' extremes d',dest_extr(1:6,l),&
             ' s ',source_extr(1:6,l),' num ',num_extr(l)
             call flush(6)
#ifdef MPI
             call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
       enddo
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      if(myrank==0)write(6,'(a)')'####################### f_dest_extr  f_source_extr #######################'
      do l=1,nlinksmpi
         if(myrank==0)write(6,'(a,i3,a,3i3,a,6i4,a,6i4,a,i4)')'dir l ',l,' disp ',&
             exmpi(l),eympi(l),ezmpi(l),' extremes d',f_dest_extr(1:6,l),&
             ' s ',f_source_extr(1:6,l),' num ',f_num_extr(l)
             call flush(6)
#ifdef MPI
             call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
       enddo
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
    if(lbuff)then
      if(myrank==0)write(6,'(a)')'####################### b_dest_extr  b_source_extr #######################'
      do l=1,nlinksmpi
         if(myrank==0)write(6,'(a,i3,a,3i3,a,6i4,a,6i4,a,i4)')'dir l ',l,' disp ',&
             exmpi(l),eympi(l),ezmpi(l),' extremes d',b_dest_extr(1:6,l),&
             ' s ',b_source_extr(1:6,l),' num ',b_num_extr(l)
             call flush(6)
#ifdef MPI
             call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
       enddo
     endif
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      
#endif

      !creo i tipi MPI contigui che mi servono per i send e receive
      !lo faccio su 13 direzioni perchè per l ed lopp le quantità da muovere sono uguali
#ifdef MPI      
      do l=1,nlinksmpi,2
        ll=(l+1)/2
        call MPI_type_contiguous(num_extr(l), MYMPIREAL, datampi(ll), ierr)
        call MPI_type_commit(datampi(ll),ierr)
#ifdef VERBOSE
        if(myrank.eq.0) then
          write(6,'(a,2i4,a,f16.8)') 'CREATE BUFFER: datampi',ll*2-1,ll*2,' (KB)-->',&
           real(num_extr(l),kind=db) *4 / 1024
          call flush(6)
        endif
#endif
      enddo
     
      !dir
      do l=1,nlinksmpi,2
        ll=(l+1)/2
        call MPI_type_contiguous(f_num_extr(l), MYMPIREAL, f_datampi(ll), ierr)
        call MPI_type_commit(f_datampi(ll),ierr)
#ifdef VERBOSE
        if(myrank.eq.0) then
          write(6,'(a,2i4,a,f16.8)') 'CREATE BUFFER: f_datampi',ll*2-1,ll*2,' (KB)-->', &
          real(f_num_extr(l),kind=db) *4 / 1024
          call flush(6)
        endif
#endif
      enddo
      
      !dir
    if(lbuff)then
      do l=1,nlinksmpi,2
        ll=(l+1)/2
        call MPI_type_contiguous(b_num_extr(l), MYMPIREAL, b_datampi(ll), ierr)
        call MPI_type_commit(b_datampi(ll),ierr)
#ifdef VERBOSE
        if(myrank.eq.0) then
          write(6,'(a,2i4,a,f16.8)') 'CREATE BUFFER: b_datampi',ll*2-1,ll*2,' (KB)-->', &
          real(b_num_extr(l),kind=db) *4 / 1024
          call flush(6)
        endif
#endif
      enddo          
    endif
    
      !dir
      do l=1,nlinksmpi,2
        ll=(l+1)/2
        call MPI_type_contiguous(i_num_extr(l), MYMPIINTS, i_datampi(ll), ierr)
        call MPI_type_commit(i_datampi(ll),ierr)
#ifdef VERBOSE
        if(myrank.eq.0) then
          write(6,'(a,2i4,a,f16.8)') 'CREATE BUFFER: i_datampi',ll*2-1,ll*2,' (KB)-->', &
          real(i_num_extr(l),kind=db) *1 / 1024
          call flush(6)
        endif
#endif
      enddo

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
    
    !alloca i buffer per mandare e ricevere
    numtot_extr=sum(num_extr)
    ll=0
    do l=1,nlinksmpi
      nbuffmpi_dest(l)=ll+1
      if(ldestpop_dir(l))ll=ll+num_extr(l)
    enddo
    allocate(dest_buffmpi(ll))
    dest_buffmpi=real(0.d0,kind=db)
    
    ll=0
    do l=1,nlinksmpi
      nbuffmpi_source(l)=ll+1
      if(lsourcepop_dir(l))ll=ll+num_extr(l)
    enddo
    allocate(source_buffmpi(ll))  
    source_buffmpi=real(0.d0,kind=db)
    
    f_numtot_extr=sum(f_num_extr)
    ll=0
    do l=1,nlinksmpi
      f_nbuffmpi_dest(l)=ll+1
      if(ldest_dir(l))ll=ll+f_num_extr(l)
    enddo
    allocate(f_dest_buffmpi(ll))
    f_dest_buffmpi=real(0.d0,kind=db)
    
    ll=0
    do l=1,nlinksmpi
      f_nbuffmpi_source(l)=ll+1
      if(lsource_dir(l))ll=ll+f_num_extr(l)
    enddo
    allocate(f_source_buffmpi(ll)) 
    f_source_buffmpi=real(0.d0,kind=db)
    
    if(lbuff)then
      b_numtot_extr=sum(b_num_extr)
      ll=0
      do l=1,nlinksmpi
        b_nbuffmpi_dest(l)=ll+1
        if(ldest_dir(l))ll=ll+b_num_extr(l)
      enddo
      allocate(b_dest_buffmpi(ll))
      b_dest_buffmpi=real(0.d0,kind=db)
      
      ll=0
      do l=1,nlinksmpi
        b_nbuffmpi_source(l)=ll+1
        if(lsource_dir(l))ll=ll+b_num_extr(l)
      enddo
      allocate(b_source_buffmpi(ll))
      b_source_buffmpi=real(0.d0,kind=db)
    endif
    
    i_numtot_extr=sum(i_num_extr)
    ll=0
    do l=1,nlinksmpi
      i_nbuffmpi_dest(l)=ll+1
      if(ldest_dir(l))ll=ll+i_num_extr(l)
    enddo
    allocate(i_dest_buffmpi(ll))
    i_dest_buffmpi=int(0,kind=isf)
    
    ll=0
    do l=1,nlinksmpi
      i_nbuffmpi_source(l)=ll+1
      if(lsource_dir(l))ll=ll+i_num_extr(l)
    enddo
    allocate(i_source_buffmpi(ll)) 
    i_source_buffmpi=int(0,kind=isf)
    
#ifdef MPI
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif    
    
    
!
#ifdef MEM_CHECK
      if(myrank == 0) then
         mem_stop = get_mem();
         write(6,*) "MEM_CHECK: after sub. input mem =", mem_stop
      endif
#endif
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. input"
        endif
#endif
      

      
      end subroutine setup_mpi
      
      subroutine setup_intbb
      
       implicit none
       
       integer(kind=isf), allocatable, dimension(:,:,:) :: myisfluid
       integer :: lmio
       
       !controllo se ho almeno un isfluid==0 nelle cornici
       !se si lintbb_dir è true lungo la direzione lmio
       do lmio=1,nlinksmpi
         allocate(myisfluid(dest_extr(1,lmio):dest_extr(2,lmio),&
          dest_extr(3,lmio):dest_extr(4,lmio),&
          dest_extr(5,lmio):dest_extr(6,lmio)))
         
         myisfluid(dest_extr(1,lmio):dest_extr(2,lmio),&
          dest_extr(3,lmio):dest_extr(4,lmio),&
          dest_extr(5,lmio):dest_extr(6,lmio)) = &
         isfluid(dest_extr(1,lmio):dest_extr(2,lmio),&
          dest_extr(3,lmio):dest_extr(4,lmio),&
          dest_extr(5,lmio):dest_extr(6,lmio))
         !lungo il settore di estremi lungo lmio (1-6 facce, 7-18 edges, 19-26 corner)
         !solo se ho almeno un isfluid==0 e ci sono popolazioni entranti fai bounce back 
         lintbb_dir(lmio)=((any(myisfluid==0)) .and. num_links_pops(lmio)>0)
         deallocate(myisfluid)
       enddo
       
      end subroutine setup_intbb
      
      subroutine setup_bb
      
       implicit none
       
       integer(kind=isf), allocatable, dimension(:,:,:) :: myisfluid
       
       allocate(myisfluid(1:nx,1:ny,1:nz))
       myisfluid(1:nx,1:ny,1:nz)  = isfluid(1:nx,1:ny,1:nz)
       lintbb=(any(myisfluid==0))
       deallocate(myisfluid)
       
       
      end subroutine setup_bb
      
      subroutine perform_pops_bb
      
       implicit none
       
       integer :: lmio,lopp,ll,l,i,j,k,ii,jj,kk
       
      !fai bounce-back nel bulk se necessario
      if(lintbb)then
        !$acc kernels present(isfluid,f) 
        !$acc loop independent collapse(3)  private(i,j,k,ii,jj,kk,lopp,l)
        do k=1,nz
          do j=1,ny
            do i=1,nx
              if(isfluid(i,j,k).ne.0)cycle
              do l=1,nlinks
                lopp=opp(l)
                ii=i+ex(l)
                jj=j+ey(l)
                kk=k+ez(l)
                f(ii,jj,kk,l)=f(i,j,k,lopp)
              enddo
            enddo
          enddo
        enddo
       !$acc end kernels
     endif
     
       !fai bounce back sulla cornice se necessario
       do lmio=1,nlinksmpi
         !lintbb_dir è true se esiste almeno in isfluid zero nell'intervallo dest_extr
         !lintbb_dir è true se nella cornice lungo la direzione lmio c'è almeno un isfluid==0
         if(.not. lintbb_dir(lmio))cycle
         do ll=1,num_links_pops(lmio)
           !trovo la popolazioni da gestire dalla lista per la direzione lmio
           !sono le stesse che uso per il send ma in direzione opposta perchè bounce-back
           l=links_pops(ll,lmio)
           lopp=opp(l)
         !$acc kernels present(dest_extr,num_links_pops,links_pops,f,isfluid) 
         !$acc loop independent collapse(3)  private(i,j,k,ii,jj,kk,ll,lopp)
           do k=dest_extr(5,lmio),dest_extr(6,lmio)
             do j=dest_extr(3,lmio),dest_extr(4,lmio)
               do i=dest_extr(1,lmio),dest_extr(2,lmio)
                 if(isfluid(i,j,k).ne. 0)cycle
                 ii=i+ex(lopp)
                 jj=j+ey(lopp)
                 kk=k+ez(lopp)
                 f(ii,jj,kk,lopp)=f(i,j,k,l) 
               enddo
             enddo
           enddo
         !$acc end kernels
         enddo
       enddo
       
      end subroutine perform_pops_bb
       
      subroutine exchange_pops_intpbc
       
       implicit none
       
       integer :: l,ll,lmio,oi,oj,ok
       !faccio le pbc per le popolazioni se interne allo stesso processo MPI, lintpbcpop_dir(lmio)=true
       !occhio se non ci sono popolazioni da mandare lintpbcpop_dir è falso
       do lmio=1,nlinksmpi
         if(.not. lintpbcpop_dir(lmio))cycle
         !scorro sul numero di popolazioni da prendere per la direzione lmio
         !$acc kernels present(dest_extr,intpbc_dir,num_links_pops,links_pops,f) async
         !$acc loop independent collapse(4) private(i,j,k,oi,oj,ok,ll,l)
         do ll=1,num_links_pops(lmio)
           do k=dest_extr(5,lmio),dest_extr(6,lmio)
             do j=dest_extr(3,lmio),dest_extr(4,lmio)
               do i=dest_extr(1,lmio),dest_extr(2,lmio)
                 !trovo la popolazioni da gestire dalla lista per la direzione lmio
                 l=links_pops(ll,lmio)
                 oi=i
                 oj=j
                 ok=k
                 if(intpbc_dir(1,lmio))oi=mod(oi+nx-1,nx)+1
                 if(intpbc_dir(2,lmio))oj=mod(oj+ny-1,ny)+1
                 if(intpbc_dir(3,lmio))ok=mod(ok+nz-1,nz)+1
                 f(oi,oj,ok,l)=f(i,j,k,l)      
               enddo
             enddo
           enddo
         
         enddo
         !$acc end kernels
       enddo
       
       !$acc wait
       
      end subroutine exchange_pops_intpbc
      
      subroutine exchange_pops_sendrecv
       
       implicit none
       
       integer :: l,ll,myoffset,tag,ierr
       
       do l=1,nlinksmpi
         !se devo mandare lungo l allora impacchetto
         if(ldestpop_dir(l))call packaging_buffmpi(l)
       enddo
       !source_buffmpi=dest_buffmpi
#ifdef MPI       
       do l=1,nlinksmpi
         ll=(l+1)/2
         if(ldestpop_dir(l))then
           myoffset=nbuffmpi_dest(l)
           !$acc host_data use_device(dest_buffmpi)
           call mpi_isend(dest_buffmpi(myoffset),1,datampi(ll),dest_dir(l), &
            mpitag(l),lbecomm,reqs_dest(l),ierr)
           !$acc end host_data
         endif
         if(lsourcepop_dir(l))then
           myoffset=nbuffmpi_source(l)
           !$acc host_data use_device(source_buffmpi)
           call mpi_irecv(source_buffmpi(myoffset),1,datampi(ll),source_dir(l), &
            mpitag(l),lbecomm,reqs_source(l),ierr)
           !$acc end host_data
         endif
       enddo
#endif
       
      
      end subroutine exchange_pops_sendrecv
      
      subroutine exchange_pops_wait
       
       implicit none
       
       integer :: l,ll,myoffset,tag,ierr
       integer, dimension(nlinksmpi) :: ierr_dest,ierr_source
#ifdef MPI 
       integer, dimension(MPI_STATUS_SIZE) :: status_dest,status_source
#endif       
      
       ierr_dest=0
       ierr_source=0
#ifdef MPI 
       do l=1,nlinksmpi
         ll=(l+1)/2
         if(ldestpop_dir(l))then
           call mpi_wait(reqs_dest(l),status_dest,ierr_dest(l))
         endif
         if(lsourcepop_dir(l))then
           call mpi_wait(reqs_source(l),status_source,ierr_source(l))
         endif
       enddo
#endif      
       if(any(ierr_dest.ne.0))call doerror(6,'ERROR in mpi_wait dest') 
       if(any(ierr_source.ne.0))call doerror(6,'ERROR in mpi_wait source')
       
       do l=1,nlinksmpi
         if(lsourcepop_dir(l))call depackaging_buffmpi(l)
       enddo
       
      end subroutine exchange_pops_wait
      
      subroutine packaging_buffmpi(lmio)
      
       implicit none
       
       integer, intent(in) :: lmio
       integer :: myoffset
       
       integer :: i,j,k,l,ll,m1,m2,m3
       
       integer :: idx
       
       myoffset=nbuffmpi_dest(lmio)
       m1=dest_extr(2,lmio)-dest_extr(1,lmio)+1
       m2=dest_extr(4,lmio)-dest_extr(3,lmio)+1
       m3=dest_extr(6,lmio)-dest_extr(5,lmio)+1
       !$acc kernels present(dest_buffmpi,num_links_pops,links_pops,f,dest_extr) 
       !$acc loop independent collapse(4)  private(i,j,k,idx,l,ll)
       !scorro sul numero di popolazioni da prendere per la direzione lmio
       do ll=1,num_links_pops(lmio)
         !trovo la popolazioni da gestire dalla lista per la direzione lmio
         do k=dest_extr(5,lmio),dest_extr(6,lmio)
           do j=dest_extr(3,lmio),dest_extr(4,lmio)
             do i=dest_extr(1,lmio),dest_extr(2,lmio)
               l=links_pops(ll,lmio)
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione
               idx=myoffset+(i-dest_extr(1,lmio))+(j-dest_extr(3,lmio))*m1+(&
                k-dest_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               dest_buffmpi(idx)=f(i,j,k,l)
             enddo
           enddo
         enddo
       enddo
       !$acc end kernels
       
      end subroutine packaging_buffmpi
      
      subroutine depackaging_buffmpi(lmio)
      
       implicit none
       
       integer, intent(in) :: lmio
       integer :: myoffset
       
       integer :: i,j,k,l,ll,m1,m2,m3
       
       integer :: idx
       
       myoffset=nbuffmpi_source(lmio)
       m1=source_extr(2,lmio)-source_extr(1,lmio)+1
       m2=source_extr(4,lmio)-source_extr(3,lmio)+1
       m3=source_extr(6,lmio)-source_extr(5,lmio)+1
       !$acc kernels present(source_buffmpi,num_links_pops,links_pops,f,source_extr) 
       !$acc loop independent collapse(4)  private(i,j,k,idx,l,ll)
       !scorro sul numero di popolazioni da prendere per la direzione lmio
       do ll=1,num_links_pops(lmio)
         do k=source_extr(5,lmio),source_extr(6,lmio)
           do j=source_extr(3,lmio),source_extr(4,lmio)
             do i=source_extr(1,lmio),source_extr(2,lmio)
               !trovo la popolazioni da gestire dalla lista per la direzione lmio
               l=links_pops(ll,lmio)
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione
               idx=myoffset+(i-source_extr(1,lmio))+(j-source_extr(3,lmio))*m1+(&
                k-source_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               f(i,j,k,l)=source_buffmpi(idx)
             enddo
           enddo
         enddo
       enddo
       !$acc end kernels
       
      end subroutine depackaging_buffmpi
      
      subroutine exchange_float_intpbc
       
       implicit none
       
       integer :: lmio,oi,oj,ok
       
       do lmio=1,nlinksmpi
         if(.not. lintpbc_dir(lmio))cycle
         !$acc kernels present(intpbc_dir,rho) 
         !$acc loop independent collapse(3)  private(i,j,k,oi,oj,ok)
         do k=f_source_extr(5,lmio),f_source_extr(6,lmio)
           do j=f_source_extr(3,lmio),f_source_extr(4,lmio)
             do i=f_source_extr(1,lmio),f_source_extr(2,lmio)
               oi=i
               oj=j
               ok=k
               if(intpbc_dir(1,lmio))oi=mod(oi+nx-1,nx)+1
               if(intpbc_dir(2,lmio))oj=mod(oj+ny-1,ny)+1
               if(intpbc_dir(3,lmio))ok=mod(ok+nz-1,nz)+1
               rho(i,j,k)=rho(oi,oj,ok)
             enddo
           enddo
         enddo
         !$acc end kernels
       enddo
       
      end subroutine exchange_float_intpbc
      
      subroutine exchange_float_sendrecv
       
       implicit none
       
       integer :: l,ll,myoffset,tag,ierr
       
       do l=1,nlinksmpi
         if(ldest_dir(l))call packaging_float_buffmpi(l)
       enddo
       !f_source_buffmpi=f_dest_buffmpi
#ifdef MPI       
       do l=1,nlinksmpi
         ll=(l+1)/2
         if(ldest_dir(l))then
           myoffset=f_nbuffmpi_dest(l)
           call mpi_isend(f_dest_buffmpi(myoffset),1,f_datampi(ll),dest_dir(l), &
            f_mpitag(l),lbecomm,f_reqs_dest(l),ierr)
         endif
         if(lsource_dir(l))then
           myoffset=f_nbuffmpi_source(l)
           call mpi_irecv(f_source_buffmpi(myoffset),1,f_datampi(ll),source_dir(l), &
            f_mpitag(l),lbecomm,f_reqs_source(l),ierr)
         endif
       enddo
#endif
       
      
      end subroutine exchange_float_sendrecv
      
      subroutine exchange_float_wait
       
       implicit none
       
       integer :: l,ll,myoffset,tag,ierr
       integer, dimension(nlinksmpi) :: ierr_dest,ierr_source
#ifdef MPI 
       integer, dimension(MPI_STATUS_SIZE) :: status_dest,status_source
#endif       
      
       ierr_dest=0
       ierr_source=0
#ifdef MPI 
       do l=1,nlinksmpi
         ll=(l+1)/2
         if(ldest_dir(l))then
           call mpi_wait(f_reqs_dest(l),status_dest,ierr_dest(l))
         endif
         if(lsource_dir(l))then
           call mpi_wait(f_reqs_source(l),status_source,ierr_source(l))
         endif
       enddo
#endif      
       if(any(ierr_dest.ne.0))call doerror(6,'ERROR in mpi_wait dest') 
       if(any(ierr_source.ne.0))call doerror(6,'ERROR in mpi_wait source')
       
       do l=1,nlinksmpi
         if(lsource_dir(l))call depackaging_float_buffmpi(l)
       enddo
       
      end subroutine exchange_float_wait
      
      subroutine packaging_float_buffmpi(lmio)
      
       implicit none
       
       integer, intent(in) :: lmio
       integer :: myoffset
       
       integer :: i,j,k,l,ll,m1,m2,m3
       
       integer :: idx
       
       myoffset=f_nbuffmpi_dest(lmio)
       m1=f_dest_extr(2,lmio)-f_dest_extr(1,lmio)+1
       m2=f_dest_extr(4,lmio)-f_dest_extr(3,lmio)+1
       m3=f_dest_extr(6,lmio)-f_dest_extr(5,lmio)+1
       !scorro sul numero di campi da prendere (per scalare = 1)
       ll=1
         !$acc kernels present(f_dest_buffmpi,rho,f_dest_extr) 
         !$acc loop independent collapse(3)  private(i,j,k,idx)
         do k=f_dest_extr(5,lmio),f_dest_extr(6,lmio)
           do j=f_dest_extr(3,lmio),f_dest_extr(4,lmio)
             do i=f_dest_extr(1,lmio),f_dest_extr(2,lmio)
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione
               idx=myoffset+(i-f_dest_extr(1,lmio))+(j-f_dest_extr(3,lmio))*m1+(&
                k-f_dest_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               
               f_dest_buffmpi(idx)=rho(i,j,k)
             enddo
           enddo
         enddo
         !$acc end kernels
       
       
      end subroutine packaging_float_buffmpi
      
      subroutine depackaging_float_buffmpi(lmio)
      
       implicit none
       
       integer, intent(in) :: lmio
       integer :: myoffset
       
       integer :: i,j,k,l,ll,m1,m2,m3
       
       integer :: idx
       
       myoffset=f_nbuffmpi_source(lmio)
       m1=f_source_extr(2,lmio)-f_source_extr(1,lmio)+1
       m2=f_source_extr(4,lmio)-f_source_extr(3,lmio)+1
       m3=f_source_extr(6,lmio)-f_source_extr(5,lmio)+1
       !scorro sul numero di campi da prendere (per scalare = 1)
       ll=1
         !$acc kernels present(f_source_buffmpi,rho,f_source_extr) 
         !$acc loop independent collapse(3)  private(i,j,k,idx)
         do k=f_source_extr(5,lmio),f_source_extr(6,lmio)
           do j=f_source_extr(3,lmio),f_source_extr(4,lmio)
             do i=f_source_extr(1,lmio),f_source_extr(2,lmio)
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione
               idx=myoffset+(i-f_source_extr(1,lmio))+(j-f_source_extr(3,lmio))*m1+(&
                k-f_source_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               
               rho(i,j,k)=f_source_buffmpi(idx)
             enddo
           enddo
         enddo
         !$acc end kernels
       
       
      end subroutine depackaging_float_buffmpi
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine exchange_isf_intpbc
       
       implicit none
       
       integer :: lmio,oi,oj,ok
       
       do lmio=1,nlinksmpi
         if(.not. lintpbc_dir(lmio))cycle
         do k=f_source_extr(5,lmio),f_source_extr(6,lmio)
           ok=k
           if(intpbc_dir(3,lmio))ok=mod(ok+nz-1,nz)+1
           do j=f_source_extr(3,lmio),f_source_extr(4,lmio)
             oj=j
             if(intpbc_dir(2,lmio))oj=mod(oj+ny-1,ny)+1
             do i=f_source_extr(1,lmio),f_source_extr(2,lmio)
               oi=i
               if(intpbc_dir(1,lmio))oi=mod(oi+nx-1,nx)+1
               isfluid(i,j,k)=isfluid(oi,oj,ok)
             enddo
           enddo
         enddo
       enddo
       
      end subroutine exchange_isf_intpbc
      
      subroutine exchange_isf_sendrecv
       
       implicit none
       
       integer :: l,ll,myoffset,tag,ierr
       
       do l=1,nlinksmpi
         if(ldest_dir(l))call packaging_isf_buffmpi(l)
       enddo
       !i_source_buffmpi=i_dest_buffmpi
#ifdef MPI       
       do l=1,nlinksmpi
         ll=(l+1)/2
         if(ldest_dir(l))then
           myoffset=i_nbuffmpi_dest(l)
           call mpi_isend(i_dest_buffmpi(myoffset),1,i_datampi(ll),dest_dir(l), &
            i_mpitag(l),lbecomm,i_reqs_dest(l),ierr)
         endif
         if(lsource_dir(l))then
           myoffset=i_nbuffmpi_source(l)
           call mpi_irecv(i_source_buffmpi(myoffset),1,i_datampi(ll),source_dir(l), &
            i_mpitag(l),lbecomm,i_reqs_source(l),ierr)
         endif
       enddo
#endif
       
      
      end subroutine exchange_isf_sendrecv
      
      subroutine exchange_isf_wait
       
       implicit none
       
       integer :: l,ll,myoffset,tag,ierr
       integer, dimension(nlinksmpi) :: ierr_dest,ierr_source
#ifdef MPI 
       integer, dimension(MPI_STATUS_SIZE) :: status_dest,status_source
#endif       
      
       ierr_dest=0
       ierr_source=0
#ifdef MPI 
       do l=1,nlinksmpi
         ll=(l+1)/2
         if(ldest_dir(l))then
           call mpi_wait(i_reqs_dest(l),status_dest,ierr_dest(l))
         endif
         if(lsource_dir(l))then
           call mpi_wait(i_reqs_source(l),status_source,ierr_source(l))
         endif
       enddo
#endif      
       if(any(ierr_dest.ne.0))call doerror(6,'ERROR in mpi_wait dest') 
       if(any(ierr_source.ne.0))call doerror(6,'ERROR in mpi_wait source')
       
       do l=1,nlinksmpi
         if(lsource_dir(l))call depackaging_isf_buffmpi(l)
       enddo
       
      end subroutine exchange_isf_wait
      
      subroutine packaging_isf_buffmpi(lmio)
      
       implicit none
       
       integer, intent(in) :: lmio
       integer :: myoffset
       
       integer :: i,j,k,l,ll,m1,m2,m3
       
       integer :: idx
       
       myoffset=i_nbuffmpi_dest(lmio)
       m1=f_dest_extr(2,lmio)-f_dest_extr(1,lmio)+1
       m2=f_dest_extr(4,lmio)-f_dest_extr(3,lmio)+1
       m3=f_dest_extr(6,lmio)-f_dest_extr(5,lmio)+1
       !scorro sul numero di campi da prendere (per scalare = 1)
       ll=1
         do k=f_dest_extr(5,lmio),f_dest_extr(6,lmio)
           do j=f_dest_extr(3,lmio),f_dest_extr(4,lmio)
             do i=f_dest_extr(1,lmio),f_dest_extr(2,lmio)
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione
               idx=myoffset+(i-f_dest_extr(1,lmio))+(j-f_dest_extr(3,lmio))*m1+(&
                k-f_dest_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               
               i_dest_buffmpi(idx)=isfluid(i,j,k)
             enddo
           enddo
         enddo
       
       
       
      end subroutine packaging_isf_buffmpi
      
      subroutine depackaging_isf_buffmpi(lmio)
      
       implicit none
       
       integer, intent(in) :: lmio
       integer :: myoffset
       
       integer :: i,j,k,l,ll,m1,m2,m3
       
       integer :: idx
       
       myoffset=i_nbuffmpi_source(lmio)
       m1=f_source_extr(2,lmio)-f_source_extr(1,lmio)+1
       m2=f_source_extr(4,lmio)-f_source_extr(3,lmio)+1
       m3=f_source_extr(6,lmio)-f_source_extr(5,lmio)+1
       !scorro sul numero di campi da prendere (per scalare = 1)
       ll=1
         do k=f_source_extr(5,lmio),f_source_extr(6,lmio)
           do j=f_source_extr(3,lmio),f_source_extr(4,lmio)
             do i=f_source_extr(1,lmio),f_source_extr(2,lmio)
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione
               idx=myoffset+(i-f_source_extr(1,lmio))+(j-f_source_extr(3,lmio))*m1+(&
                k-f_source_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               
               isfluid(i,j,k)=i_source_buffmpi(idx)
             enddo
           enddo
         enddo
       
       
       
      end subroutine depackaging_isf_buffmpi
 
      subroutine write_file_vtk_par(iframe,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
       implicit none
       
       integer, intent(in) ::iframe
       integer, intent(out) :: e_io
#ifdef MPI         
       integer(kind=MPI_OFFSET_KIND) :: tempoffset
#endif
       integer :: nns
       character(len=500) :: sheadervtk
       
       integer :: ioffset
       character(1), parameter :: end_rec = char(10)
       integer, parameter :: bytechar=kind(end_rec)
       integer, parameter :: byteint = 4
       integer, parameter :: byter4  = 4
       integer, parameter :: byter8  = 8
       integer, parameter :: nbuffsub = 0
       integer :: filetypesub,imemtype,filetypesubv,fdens,fvel,ierr
       
       integer, dimension(3) :: memDims,memOffs
       integer, dimension(4) :: velglobalDims,velldims,velmystarts, &
        velmemDims,velmemOffs

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!density!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       
       sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
        '_'//trim(write_fmtnumb(iframe)) // '.vti'
       
       
#ifdef MPI                 
       
       call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sevt1), &
			MPI_MODE_CREATE + MPI_MODE_WRONLY, &
			MPI_INFO_NULL,fdens,e_io)
                       
      tempoffset=int(0,kind=MPI_OFFSET_KIND)
      
      sheadervtk=repeat(' ',500)
      sheadervtk=headervtk(1)
      nns=nheadervtk(1)
      
      if(myrank==0)call MPI_File_write_at(fdens,tempoffset,sheadervtk(1:nns),nns, &
       MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)

      ioffset=vtkoffset(1)
      tempoffset=int(ioffset,kind=MPI_OFFSET_KIND)
      
      if(myrank==0)call MPI_File_write_at(fdens,tempoffset,int(ndatavtk(1),kind=4),1, &
        MPI_INTEGER,MPI_STATUS_IGNORE,e_io)
 
      call MPI_Type_create_subarray(3,gsizes,lsizes,myoffset, &
       MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)
        
      call MPI_Type_commit(filetypesub, e_io)
      
      ioffset=vtkoffset(1)+byteint
      tempoffset=int(ioffset,kind=MPI_OFFSET_KIND)
   
      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
       "native",MPI_INFO_NULL,e_io)
      ! We need full local sizes: memDims
      memDims = lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]
  

      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes,memOffs, &
       MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fdens,rhoprint,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      ioffset=vtkoffset(1)+byteint+ndatavtk(1)
      tempoffset=int(ioffset,kind=MPI_OFFSET_KIND)
     
      if(myrank==0)call MPI_File_write_at(fdens,tempoffset,footervtk(1),30, &
       MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)
      
      call MPI_FILE_CLOSE(fdens,e_io)
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!velocity!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
       '_'//trim(write_fmtnumb(iframe)) // '.vti'
      
      call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(sevt2), &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL,fvel,e_io)
                       
                       
      tempoffset=int(0,kind=MPI_OFFSET_KIND)
      
      sheadervtk=repeat(' ',500)
      sheadervtk=headervtk(2)
      nns=nheadervtk(2)
      
      if(myrank==0)call MPI_File_write_at(fvel,tempoffset,sheadervtk(1:nns),nns, &
       MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)
      
      velglobalDims(1)=3
      velglobalDims(2:4)=gsizes(1:3)
      velldims(1)=3
      velldims(2:4)=lsizes(1:3)
      velmystarts(1) = 0
      velmystarts(2:4) = myoffset(1:3)
      
      ioffset=vtkoffset(2)
      tempoffset=int(ioffset,kind=MPI_OFFSET_KIND)
      
      if(myrank==0)call MPI_File_write_at(fvel,tempoffset,int(ndatavtk(2),kind=4),1, &
        MPI_INTEGER,MPI_STATUS_IGNORE,e_io)
 
  
      call MPI_Type_create_subarray(4,velglobalDims,velldims,velmystarts, &
       MPI_ORDER_FORTRAN,MPI_REAL,filetypesubv,e_io)
        
      call MPI_Type_commit(filetypesubv, e_io)
  
      ioffset=vtkoffset(2)+byteint
      tempoffset=int(ioffset,kind=MPI_OFFSET_KIND)
   
      call MPI_File_Set_View(fvel,tempoffset,MPI_REAL,filetypesubv, &
        "native",MPI_INFO_NULL,e_io)
      ! We need full local sizes: memDims
      velmemDims(1) = vellDims(1)
      velmemDims(2:4) = vellDims(2:4) + 2*nbuffsub
      velmemOffs = [ 0, nbuffsub, nbuffsub, nbuffsub ]

      call MPI_TYPE_CREATE_SUBARRAY(4,velmemDims,velldims,velmemOffs, &
       MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fvel,velprint,1,imemtype,MPI_STATUS_IGNORE,e_io)
  
      ioffset=vtkoffset(2)+byteint+ndatavtk(2)
      tempoffset=int(ioffset,kind=MPI_OFFSET_KIND)
      
      if(myrank==0)call MPI_File_write_at(fvel,tempoffset,footervtk(2),30, &
       MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)
      
      call MPI_FILE_CLOSE(fvel, e_io)
      
#endif                       
      return
  
      end subroutine write_file_vtk_par
      
      subroutine write_file_raw_par(iframe,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
       implicit none
       
       integer, intent(in) ::iframe
       integer, intent(out) :: e_io
#ifdef MPI         
       integer(kind=MPI_OFFSET_KIND) :: tempoffset
#endif
       
       integer :: ioffset
       character(1), parameter :: end_rec = char(10)
       integer, parameter :: bytechar=kind(end_rec)
       integer, parameter :: byteint = 4
       integer, parameter :: byter4  = 4
       integer, parameter :: byter8  = 8
       integer, parameter :: nbuffsub = 0
       integer :: filetypesub,imemtype,filetypesubv,fdens,fvel,ierr
       
       integer, dimension(3) :: memDims,memOffs
       integer, dimension(4) :: velglobalDims,velldims,velmystarts, &
        velmemDims,velmemOffs

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!density!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       
       sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
        '_'//trim(write_fmtnumb(iframe)) // '.raw'
       
       
#ifdef MPI                 
       
       call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sevt1), &
			MPI_MODE_CREATE + MPI_MODE_WRONLY, &
			MPI_INFO_NULL,fdens,e_io)
                       
      tempoffset=int(0,kind=MPI_OFFSET_KIND)
      
      
      call MPI_Type_create_subarray(3,gsizes,lsizes,myoffset, &
       MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)
        
      call MPI_Type_commit(filetypesub, e_io)
      
      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
       "native",MPI_INFO_NULL,e_io)
      ! We need full local sizes: memDims
      memDims = lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]
  
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes,memOffs, &
       MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fdens,rhoprint,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      call MPI_FILE_CLOSE(fdens,e_io)
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!velocity!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
       '_'//trim(write_fmtnumb(iframe)) // '.raw'
      
      call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(sevt2), &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL,fvel,e_io)
                       
                       
      tempoffset=int(0,kind=MPI_OFFSET_KIND)
      
      
      velglobalDims(1)=3
      velglobalDims(2:4)=gsizes(1:3)
      velldims(1)=3
      velldims(2:4)=lsizes(1:3)
      velmystarts(1) = 0
      velmystarts(2:4) = myoffset(1:3)
  
      call MPI_Type_create_subarray(4,velglobalDims,velldims,velmystarts, &
       MPI_ORDER_FORTRAN,MPI_REAL,filetypesubv,e_io)
        
      call MPI_Type_commit(filetypesubv, e_io)
   
      call MPI_File_Set_View(fvel,tempoffset,MPI_REAL,filetypesubv, &
        "native",MPI_INFO_NULL,e_io)
      ! We need full local sizes: memDims
      velmemDims(1) = vellDims(1)
      velmemDims(2:4) = vellDims(2:4) + 2*nbuffsub
      velmemOffs = [ 0, nbuffsub, nbuffsub, nbuffsub ]

      call MPI_TYPE_CREATE_SUBARRAY(4,velmemDims,velldims,velmemOffs, &
       MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fvel,velprint,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      call MPI_FILE_CLOSE(fvel, e_io)
      
#endif                       
      return
  
      end subroutine write_file_raw_par      
      
      subroutine write_file_raw_par2D(iframe,mydir,mypoint,service1,service3, &
       myoffset_plane,lsizes_plane,gsizes_plane,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
      implicit none
       
      integer, intent(in) ::iframe,mydir,mypoint
      real(4), dimension(:,:,:), allocatable :: service1
      real(4), dimension(:,:,:,:), allocatable :: service3
      integer, dimension(mpid), intent(in) :: myoffset_plane,lsizes_plane,gsizes_plane
      integer, intent(out) :: e_io
#ifdef MPI         
      integer(kind=MPI_OFFSET_KIND) :: tempoffset
#endif
       
      integer :: ioffset
      character(1), parameter :: end_rec = char(10)
      integer, parameter :: bytechar=kind(end_rec)
      integer, parameter :: byteint = 4
      integer, parameter :: byter4  = 4
      integer, parameter :: byter8  = 8
      integer, parameter :: nbuffsub = 0
      integer :: filetypesub,imemtype,filetypesubv,fdens,fvel,ierr
       
      integer, dimension(3) :: memDims,memOffs
      integer, dimension(4) :: velglobalDims,velldims,velmystarts, &
       velmemDims,velmemOffs
       
       
#ifdef MPI       
      
      !qui ci arrivano solo i processi buoni, quelli pupu li ho esclusi prima di chiamare la sub
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!density!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      
      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
       '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(1))// &
       '_'//trim(write_fmtnumb(iframe)) // '.raw'           
       
      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sevt1), &
			MPI_MODE_CREATE + MPI_MODE_WRONLY, &
			MPI_INFO_NULL,fdens,e_io)
                       
      tempoffset=int(0,kind=MPI_OFFSET_KIND)
      
      
      call MPI_Type_create_subarray(3,gsizes_plane,lsizes_plane,myoffset_plane, &
       MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)
        
      call MPI_Type_commit(filetypesub, e_io)
      
      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
       "native",MPI_INFO_NULL,e_io)
      ! We need full local sizes: memDims
      memDims = lsizes_plane + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]
  
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes_plane,memOffs, &
       MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fdens,rhoprint,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      call MPI_FILE_CLOSE(fdens,e_io)
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!velocity!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
       '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(2))// &
       '_'//trim(write_fmtnumb(iframe)) // '.raw'
      
      
      call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(sevt2), &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL,fvel,e_io)
                       
                       
      tempoffset=int(0,kind=MPI_OFFSET_KIND)
      
      
      velglobalDims(1)=3
      velglobalDims(2:4)=gsizes_plane(1:3)
      velldims(1)=3
      velldims(2:4)=lsizes_plane(1:3)
      velmystarts(1) = 0
      velmystarts(2:4) = myoffset_plane(1:3)
  
      call MPI_Type_create_subarray(4,velglobalDims,velldims,velmystarts, &
       MPI_ORDER_FORTRAN,MPI_REAL,filetypesubv,e_io)
        
      call MPI_Type_commit(filetypesubv, e_io)
   
      call MPI_File_Set_View(fvel,tempoffset,MPI_REAL,filetypesubv, &
        "native",MPI_INFO_NULL,e_io)
      ! We need full local sizes: memDims
      velmemDims(1) = vellDims(1)
      velmemDims(2:4) = vellDims(2:4) + 2*nbuffsub
      velmemOffs = [ 0, nbuffsub, nbuffsub, nbuffsub ]

      call MPI_TYPE_CREATE_SUBARRAY(4,velmemDims,velldims,velmemOffs, &
       MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fvel,velprint,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      call MPI_FILE_CLOSE(fvel, e_io)
      
#endif                       
      return
  
      end subroutine write_file_raw_par2D      
      
      function GET_COORD_POINT(ii,jj,kk)
     
      implicit none

      integer, intent(in) :: ii,jj,kk
     
      integer :: i,j,k
      integer, dimension(mpid) :: GET_COORD_POINT
        
        do i=0,proc_x-1
          if(ii<=xfindom(i))then
            GET_COORD_POINT(1)=i
            exit
          endif
        enddo  
        
        do j=0,proc_y-1
          if(jj<=yfindom(j))then
            GET_COORD_POINT(2)=j
            exit
          endif
        enddo  
        
        do k=0,proc_z-1
          if(kk<=zfindom(k))then
            GET_COORD_POINT(3)=k
            exit
          endif
        enddo   
        
      
    end function GET_COORD_POINT
    
    function GET_RANK_POINT(ii,jj,kk,ierr)
    
      implicit none

      integer, intent(in) :: ii,jj,kk
      integer, dimension(mpid) :: temp
      integer :: ierr
      
      integer :: GET_RANK_POINT
      
      temp=GET_COORD_POINT(ii,jj,kk)
#ifdef MPI       
      call MPI_Cart_rank(lbecomm, temp, GET_RANK_POINT,ierr)
#else
      GET_RANK_POINT=0
#endif
           
    end function GET_RANK_POINT
    
    subroutine dostop(mystring)
    
     implicit none
     
     integer :: ierr
     
     character(len=*), optional :: mystring
     
     if(present(mystring))then
       if(myrank==0)then
         write(6,'(a)')mystring
         call flush(6)
       endif
     endif
     
#ifdef MPI
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_finalize(ierr)
#endif
     stop
    
    end subroutine dostop
    
    subroutine doerror(errcode,mystring)
    
     implicit none
     
     integer, value :: errcode
     integer :: ierr
     
     character(len=*), optional :: mystring
     
     if(present(mystring))then
         write(6,'(a)')mystring
         call flush(6)
     endif
     
#ifdef MPI
     call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
#endif
     stop
    
    end subroutine doerror
    
     subroutine or_world_larr(argument)
 
!***********************************************************************
!     
!     LBsoft global 'logical or' subroutine for a logical array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  
  
  logical, intent(inout) :: argument
  logical, dimension(1) :: buffersub,lbuffer
  
  integer ierr
  
#ifdef MPI  
    buffersub(1)=argument
    
    call MPI_ALLREDUCE(buffersub,lbuffer,1,MPI_LOGICAL, &
      MPI_LOR,MPI_COMM_WORLD,ierr)
      
    argument=lbuffer(1)
    
  
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
  
  return
  
 end subroutine or_world_larr
    
  end module mpi_template
