 module prints
    
    use vars
    
    implicit none
        
    contains
  
    subroutine header_vtk(nx,ny,nz,mystring500,namevar,extent,ncomps,iinisub,iend,myoffset, &
      new_myoffset,indent)
    
      implicit none
    
      integer, intent(in) :: nx,ny,nz
      character(len=8),intent(in) :: namevar
      character(len=120),intent(in) :: extent
      integer, intent(in) :: ncomps,iinisub,myoffset
      integer, intent(out) :: iend,new_myoffset
      integer, intent(inout) :: indent
    
      !namevar='density1'
    
      character(len=500), intent(out) :: mystring500
      ! End-character for binary-record finalize.
      character(1), parameter:: end_rec = char(10) 
      character(1) :: string1
      character(len=*),parameter :: topology='ImageData' 
      integer :: ioffset,nele,bytechar,byteint,byter4,byter8,iini
      
      iini=iinisub
      bytechar=kind(end_rec)
      byteint=kind(iini)
      byter4  = 4
      byter8  = 8
      
      mystring500=repeat(' ',500)
      
      iend=iini
      
      iini=iend+1
      nele=22
      iend=iend+nele
      mystring500(iini:iend)='<?xml version="1.0"?>'//end_rec
      
      new_myoffset=myoffset
      new_myoffset = new_myoffset + nele * bytechar
    
      
      iini=iend+1
      nele=67
      iend=iend+nele
      if(lelittle)then  
        mystring500(iini:iend) = '<VTKFile type="'//trim(topology)// &
        '" version="0.1" byte_order="LittleEndian">'//end_rec
      else
        mystring500(iini:iend) = '<VTKFile type="'//trim(topology)// &
        '" version="0.1" byte_order="BigEndian">   '//end_rec
      endif
    
      new_myoffset = new_myoffset + 67 * bytechar
  
    
      indent = indent + 2
      iini=iend+1
      nele=70
      iend=iend+nele
      mystring500(iini:iend) = repeat(' ',indent)//'<'//trim(topology)//' WholeExtent="'//&
                    trim(extent)//'">'//end_rec
    

      new_myoffset = new_myoffset + 70 * bytechar
    
      
      indent = indent + 2
      iini=iend+1
      nele=63
      iend=iend+nele
      mystring500(iini:iend) = repeat(' ',indent)//'<Piece Extent="'//trim(extent)//'">'//end_rec
      
      new_myoffset = new_myoffset + 63 * bytechar
    
      
      ! initializing offset pointer
      ioffset = 0 
      
      indent = indent + 2
      iini=iend+1
      nele=18
      iend=iend+nele
      mystring500(iini:iend)=repeat(' ',indent)//'<PointData>'//end_rec
      
      new_myoffset = new_myoffset + 18 * bytechar
      
      indent = indent + 2
      iini=iend+1
      nele=115
      iend=iend+nele
      
      if(ncomps/=1 .and. ncomps/=3)then
        write(6,'(a)')'ERROR in header_vtk'
        stop
      endif
      write(string1,'(i1)')ncomps
      mystring500(iini:iend)=repeat(' ',indent)//'<DataArray type="Float32" Name="'// &
      namevar//'" NumberOfComponents="'//string1// '" '//&
      'format="appended" offset="'//space_fmtnumb12(ioffset)//'"/>'//end_rec
      
      new_myoffset = new_myoffset + 115 * bytechar
    
      
      indent = indent - 2
      iini=iend+1
      nele=19
      iend=iend+nele
      mystring500(iini:iend)=repeat(' ',indent)//'</PointData>'//end_rec
      
      new_myoffset = new_myoffset + 19 * bytechar
      
      
      indent = indent - 2
      iini=iend+1
      nele=13
      iend=iend+nele
      mystring500(iini:iend)=repeat(' ',indent)//'</Piece>'//end_rec
      
      
      new_myoffset = new_myoffset + 13 * bytechar
    
      
      indent = indent - 2
      iini=iend+1
      nele=15
      iend=iend+nele
      mystring500(iini:iend)=repeat(' ',indent)//'</'//trim(topology)//'>'//end_rec

      new_myoffset = new_myoffset + 15 * bytechar
    

      iini=iend+1
      nele=32
      iend=iend+nele
      mystring500(iini:iend)=repeat(' ',indent)//'<AppendedData encoding="raw">'//end_rec
      
      new_myoffset = new_myoffset + 32 * bytechar
      
      iini=iend+1
      nele=1
      iend=iend+nele
      mystring500(iini:iend)='_'
      
      new_myoffset = new_myoffset + 1 * bytechar
      
      return
    
    end subroutine header_vtk
  
    subroutine footer_vtk(nx,ny,nz,mystring30,iinisub,iend,myoffset, &
      new_myoffset,indent)
  
      implicit none
      
      integer, intent(in) :: nx,ny,nz
      integer, intent(in) :: iinisub,myoffset
      integer, intent(out) :: iend,new_myoffset
      integer, intent(inout) :: indent
      
      
      character(len=30), intent(out) :: mystring30
      ! End-character for binary-record finalize.
      character(1), parameter:: end_rec = char(10) 
      character(1) :: string1
      character(len=*),parameter :: topology='ImageData' 
      integer :: ioffset,nele,bytechar,byteint,byter4,byter8,iini
      
      iini=iinisub
      bytechar=kind(end_rec)
      byteint=kind(iini)
      byter4  = 4
      byter8  = 8
      
      mystring30=repeat(' ',30)
      
      iend=iini
      
      iini=iend+1
      nele=1
      iend=iend+nele
      mystring30(iini:iend)=end_rec
      
      new_myoffset = myoffset
      new_myoffset = new_myoffset + 1 * bytechar
    
      
      
      iini=iend+1
      nele=18
      iend=iend+nele
      mystring30(iini:iend)=repeat(' ',indent)//'</AppendedData>'//end_rec
      
      new_myoffset = new_myoffset + 18 * bytechar
      
      iini=iend+1
      nele=11
      iend=iend+nele
      mystring30(iini:iend)='</VTKFile>'//end_rec
      
      if(iend/=30)then
        write(6,'(a)')'ERROR in footer_vtk'
        stop
      endif
      
      return
    
    end subroutine footer_vtk
    
    subroutine test_little_endian(ltest)
    
          !***********************************************************************
          !     
          !     LBsoft subroutine for checking if the computing architecture
          !     is working in little-endian or big-endian
          !     
          !     licensed under Open Software License v. 3.0 (OSL-3.0)
          !     author: M. Lauricella
          !     last modification October 2019
          !     
          !***********************************************************************
    
        implicit none 
        integer, parameter :: ik1 = selected_int_kind(2) 
        integer, parameter :: ik4 = selected_int_kind(9) 
        
        logical, intent(out) :: ltest
        
        if(btest(transfer(int((/1,0,0,0/),ik1),1_ik4),0)) then 
          !it is little endian
          ltest=.true.
        else 
          !it is big endian
          ltest=.false.
        end if 
        
        return
      
    end subroutine test_little_endian 
  
    subroutine init_output(ncomp,lvtk,lraw)
  
      !***********************************************************************
      !     
      !     LBsoft subroutine for creating the folders containing the files
      !     in image VTK legacy binary format in parallel IO
      !     
      !     licensed under Open Software License v. 3.0 (OSL-3.0)
      !     author: M. Lauricella
      !     last modification October 2018
      !     
      !***********************************************************************

    
      implicit none
      
      integer, intent(in) :: ncomp
      logical, intent(in) :: lvtk,lraw
      character(len=255) :: path,makedirectory
      logical :: lexist
      
      integer :: i,j,k,nn,indent,myoffset,new_myoffset,iend
      integer, parameter :: byter4=4
      integer, parameter :: byteint=4
      integer, allocatable :: printlistvtk(:)
      integer, parameter :: ioxyz=54
      character(len=*), parameter :: filexyz='isfluid.xyz'
      character(len=120) :: mystring120
      
      call test_little_endian(lelittle)
      
      sevt1=repeat(' ',mxln)
      sevt2=repeat(' ',mxln)
      
      path = repeat(' ',255)
      call getcwd(path)
      
      !call get_environment_variable('DELIMITER',delimiter)
      path = trim(path)
      delimiter = path(1:1)
      if (delimiter==' ') delimiter='/'


      
      makedirectory=repeat(' ',255)
      makedirectory = 'output'//delimiter
      dir_out=trim(makedirectory)
#ifdef _INTEL
        inquire(directory=trim(makedirectory),exist=lexist)
#else
        inquire(file=trim(makedirectory),exist=lexist)
#endif
      
      if(.not. lexist)then
        makedirectory=repeat(' ',255)
        makedirectory = 'mkdir output'
        call system(makedirectory)
      endif
      mystring120=repeat(' ',120)
      
      
      makedirectory=repeat(' ',255)
      makedirectory=trim(path)//delimiter//'output'//delimiter
      
      extentvtk =  space_fmtnumb(1) // ' ' // space_fmtnumb(nx) // ' ' &
            // space_fmtnumb(1) // ' ' // space_fmtnumb(ny) // ' ' &
            // space_fmtnumb(1) // ' ' // space_fmtnumb(nz)
      
      if(ncomp==1)then
        nfilevtk=2
      elseif(ncomp==2)then
        nfilevtk=3
      endif
      
      allocate(printlistvtk(nfilevtk))
      do i=1,nfilevtk
        printlistvtk(i)=i
      enddo
      
      allocate(varlistvtk(nfilevtk))
      allocate(namevarvtk(nfilevtk))
      allocate(ndimvtk(nfilevtk))
      allocate(headervtk(nfilevtk))
      allocate(footervtk(nfilevtk))
      allocate(nheadervtk(nfilevtk))
      allocate(vtkoffset(nfilevtk))
      allocate(ndatavtk(nfilevtk))
      varlistvtk(1:nfilevtk)=printlistvtk(1:nfilevtk)
      
      if(ncomp==1)then
        do i=1,nfilevtk
          select case(printlistvtk(i))
          case(1)
            namevarvtk(i)='rho     '
            ndimvtk(i)=1
          case(2)
            namevarvtk(i)='vel     '
            ndimvtk(i)=3
          case default
            write(6,'(a)')'ERROR in init_output'
            stop
          end select
        enddo
      elseif(ncomp==2)then
        do i=1,nfilevtk
          select case(printlistvtk(i))
          case(1)
            namevarvtk(i)='rho1    '
            ndimvtk(i)=1
          case(2)
            namevarvtk(i)='rho2    '
            ndimvtk(i)=1
          case(3)
            namevarvtk(i)='vel     '
            ndimvtk(i)=3
          case default
            write(6,'(a)')'ERROR in init_output'
            stop
          end select
        enddo
      endif
      nn=nx*ny*nz
      
      do i=1,nfilevtk
        myoffset=0
        indent=0
        call header_vtk(nx,ny,nz,headervtk(i),namevarvtk(i),extentvtk,ndimvtk(i),0,iend,myoffset, &
        new_myoffset,indent)
        vtkoffset(i)=new_myoffset
        myoffset=new_myoffset+byteint+ndimvtk(i)*nn*byter4
        ndatavtk(i)=ndimvtk(i)*nn*byter4
        nheadervtk(i)=iend
        call footer_vtk(nx,ny,nz,footervtk(i),0,iend,myoffset, &
        new_myoffset,indent)
      enddo
      
      return

    end subroutine init_output
  
    subroutine string_char(mychar,nstring,mystring)
    
      implicit none
      
      integer :: i
      character(1), allocatable, dimension(:) :: mychar
      integer, intent(in) :: nstring
      character(len=*), intent(in) :: mystring
      
      allocate(mychar(nstring))
      
      do i=1,nstring
        mychar(i)=mystring(i:i)
      enddo
      
    end subroutine string_char
  
    function space_fmtnumb(inum)

      !***********************************************************************
      !     
      !     LBsoft function for returning the string of six characters 
      !     with integer digits and leading spaces to the left
      !     originally written in JETSPIN by M. Lauricella et al.
      !     
      !     licensed under Open Software License v. 3.0 (OSL-3.0)
      !     author: M. Lauricella
      !     last modification October 2019
      !     
      !***********************************************************************

      implicit none

      integer,intent(in) :: inum
      character(len=6) :: space_fmtnumb
      integer :: numdigit,irest
      real(kind=8) :: tmp
      character(len=22) :: cnumberlabel

      numdigit=dimenumb(inum)
      irest=6-numdigit
      if(irest>0)then
        write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
        write(space_fmtnumb,fmt=cnumberlabel)repeat(' ',irest),inum
      else
        write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
        write(space_fmtnumb,fmt=cnumberlabel)inum
      endif
      
      return

    end function space_fmtnumb
 
    function space_fmtnumb12(inum)
  
      !***********************************************************************
      !     
      !     LBsoft function for returning the string of six characters 
      !     with integer digits and leading TWELVE spaces to the left
      !     originally written in JETSPIN by M. Lauricella et al.
      !     
      !     licensed under Open Software License v. 3.0 (OSL-3.0)
      !     author: M. Lauricella
      !     last modification October 2019
      !     
      !***********************************************************************
  
        implicit none

        integer,intent(in) :: inum
        character(len=12) :: space_fmtnumb12
        integer :: numdigit,irest
        real(kind=8) :: tmp
        character(len=22) :: cnumberlabel

        numdigit=dimenumb(inum)
        irest=12-numdigit
        if(irest>0)then
          write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
          write(space_fmtnumb12,fmt=cnumberlabel)repeat(' ',irest),inum
        else
          write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
          write(space_fmtnumb12,fmt=cnumberlabel)inum
        endif
        
        return

    end function space_fmtnumb12
  
    function dimenumb(inum)
  
        !***********************************************************************
        !    
        !     LBsoft function for returning the number of digits
        !     of an integer number
        !     originally written in JETSPIN by M. Lauricella et al.
        !    
        !     licensed under the 3-Clause BSD License (BSD-3-Clause)
        !     author: M. Lauricella
        !     last modification July 2018
        !    
        !***********************************************************************

        implicit none

        integer,intent(in) :: inum
        integer :: dimenumb
        integer :: i
        real(kind=db) :: tmp

        i=1
        tmp=real(inum,kind=db)
        do
        if(tmp< 10.0_db )exit
          i=i+1
          tmp=tmp/ 10.0_db
        enddo

        dimenumb=i

        return

    end function dimenumb

    function write_fmtnumb(inum)
  
          !***********************************************************************
          !    
          !     LBsoft function for returning the string of six characters
          !     with integer digits and leading zeros to the left
          !     originally written in JETSPIN by M. Lauricella et al.
          !    
          !     licensed under the 3-Clause BSD License (BSD-3-Clause)
          !     author: M. Lauricella
          !     last modification July 2018
          !    
          !***********************************************************************
      
          implicit none

          integer,intent(in) :: inum
          character(len=6) :: write_fmtnumb
          integer :: numdigit,irest
          !real*8 :: tmp
          character(len=22) :: cnumberlabel
          
          numdigit=dimenumb(inum)
          irest=6-numdigit
          if(irest>0)then
              write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
              write(write_fmtnumb,fmt=cnumberlabel)repeat('0',irest),inum
          else
              write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
              write(write_fmtnumb,fmt=cnumberlabel)inum
          endif
      
          return
    end function write_fmtnumb   
    
    subroutine get_memory_gpu(fout,fout2)

      !***********************************************************************
      !     
      !     LBsoft subroutine for register the memory usage
      !     
      !     licensed under the 3-Clause BSD License (BSD-3-Clause)
      !     modified by: M. Lauricella
      !     last modification July 2018
      !     
      !***********************************************************************  
#ifdef _OPENACC
        use openacc
        use accel_lib
        use iso_c_binding 
#elif defined _CUDA  
        use cudafor
#endif
        
        implicit none
        
        real(kind=db), intent(out) :: fout,fout2
        real(kind=db) :: myd(2),myd2(2)
        integer :: istat
#ifdef _OPENACC  
        integer(c_size_t) :: myfree, total
#elif defined _CUDA  
        integer(kind=cuda_count_kind) :: myfree, total
#else
        integer :: myfree, total
#endif  
        
#ifdef _OPENACC
        myfree=acc_get_free_memory()
        total=acc_get_memory() 
#elif defined _CUDA
        istat = cudaMemGetInfo( myfree, total )
#else
        myfree=0
        total=0
#endif  
        fout = real(total-myfree,kind=4)/(1024.0**3.0)
        fout2 = real(total,kind=4)/(1024.0**3.0)
        
        return
  
    end subroutine get_memory_gpu
    
    subroutine print_memory_registration_gpu(iu,mybanner,mybanner2,&
      mymemory,totmem)
 
        !***********************************************************************
        !     
        !     LBcuda subroutine for printing the memory registration
        !     
        !     licensed under the 3-Clause BSD License (BSD-3-Clause)
        !     author: M. Lauricella
        !     last modification April 2022
        !     
        !***********************************************************************
          
          implicit none
          
          integer, intent(in) :: iu
          character(len=*), intent(in) :: mybanner,mybanner2
          real(kind=db), intent(in) :: mymemory,totmem
          
          character(len=12) :: r_char,r_char2
          
          character(len=*),parameter :: of='(a)'
          
          
          
        
          write (r_char,'(f12.4)')mymemory
          write (r_char2,'(f12.4)')totmem
          write(iu,of)"                                                                               "
          write(iu,of)"******************************GPU MEMORY MONITOR*******************************"
          write(iu,of)"                                                                               "
          write(iu,'(4a)')trim(mybanner)," = ",trim(adjustl(r_char))," (GB)"
          write(iu,'(4a)')trim(mybanner2)," = ",trim(adjustl(r_char2))," (GB)"
          write(iu,of)"                                                                               "
          write(iu,of)"*******************************************************************************"
          write(iu,of)"                                                                               "
          
          return
  
    end subroutine print_memory_registration_gpu
    !******************************************************************************************************!
#ifdef _OPENACC
    subroutine printDeviceProperties(ngpus,dev_Num,dev_Type,iu)
    
    
    use openacc
    use iso_c_binding 
    
    integer :: ngpus,dev_Num
    integer(acc_device_kind) :: dev_Type
  
    integer,intent(in) :: iu 
    integer(c_size_t) :: tot_mem
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
    write (iu,903) "Total Global Memory: ",real(tot_mem)/(1024.0**3.0)," Gbytes"
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
#endif  
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

  subroutine dump_distros_1c_3d
  
   implicit none
   
   sevt1 = 'dump/' // 'out'//'_'//'f'//'_'// '.raw'
   open(unit=1045,file=trim(sevt1),status='replace',form='unformatted')
   do ll=0,26
      write(1045) f(:,:,:,ll)
   enddo
   close(1045)
   
  end subroutine dump_distros_1c_3d

  subroutine read_distros_1c_3d

      implicit none

      sevt1 = 'dump/' // 'out'//'_'//'f'// '_'// '.raw'
      open(unit=1045,file=trim(sevt1), form='unformatted', status='old')
      do ll=0,26
          read(1045) f(:,:,:,ll)
      enddo
      close(1045)
  end subroutine read_distros_1c_3d

  subroutine print_raw_slice_sync(iframe)
  
   implicit none
   
   integer, intent(in) :: iframe
   !rho
   sevt1 = trim(dir_out) // 'out'//'_'//'rhoxy'// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   open(unit=745,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted')
   write(745) rho(:,:,nz/2)
   close(745)
   sevt1 = trim(dir_out) // 'out'//'_'//'rhoxz'// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   open(unit=746,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted')
   write(746) rho(:,ny/2,:)
   close(746)
   sevt1 = trim(dir_out) // 'out'//'_'//'rhoyz'// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   open(unit=747,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted')
   write(747) rho(nx/2,:,:)
   close(747)
   
   !u
   sevt1 = trim(dir_out) // 'out'//'_'//'uxy'// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   open(unit=845,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted')
   write(845) u(:,:,nz/2)
   close(845)
   sevt1 = trim(dir_out) // 'out'//'_'//'uxz'// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   open(unit=846,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted')
   write(846) u(:,ny/2,:)
   close(846)
     sevt1 = trim(dir_out) // 'out'//'_'//'uyz'// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   open(unit=847,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted')
   write(847) u(nx/2,:,:)
   close(847)

   !v
   sevt1 = trim(dir_out) // 'out'//'_'//'vxy'// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   open(unit=848,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted')
   write(848) v(:,:,nz/2)
   close(848)
   sevt1 = trim(dir_out) // 'out'//'_'//'vyz'// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   open(unit=849,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted')
   write(849) v(nx/2,:,:)
   close(849)
   sevt1 = trim(dir_out) // 'out'//'_'//'vxz'// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   open(unit=850,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted')
   write(850) v(:,ny/2,:)
   close(850)

   !w
   sevt1 = trim(dir_out) // 'out'//'_'//'wxz'// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   open(unit=851,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted')
   write(851) w(:,ny/2,:)
   close(851)
   sevt1 = trim(dir_out) // 'out'//'_'//'wyz'// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   open(unit=852,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted')
   write(852) w(nx/2,:,:)
   close(852)
   sevt1 = trim(dir_out) // 'out'//'_'//'wxy'// &
    '_'//trim(write_fmtnumb(iframe)) // '.raw'
   open(unit=853,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted')
   write(853) w(:,:,nz/2)
   close(853)
   
  end subroutine print_raw_slice_sync

  subroutine print_raw_slice_2c_sync(iframe)
  
        implicit none
   
      integer, intent(in) :: iframe
      !phi
      sevt1 = trim(dir_out) // 'out'//'_'//'phixy'// &
        '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=745,file=trim(sevt1), &
        status='replace',action='write',access='stream',form='unformatted')
      write(745) phi(:,:,nz/2)
      close(745)
      sevt1 = trim(dir_out) // 'out'//'_'//'phixz'// &
        '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=746,file=trim(sevt1), &
        status='replace',action='write',access='stream',form='unformatted')
      write(746) phi(:,ny/2,:)
      close(746)
      sevt1 = trim(dir_out) // 'out'//'_'//'phiyz'// &
        '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=747,file=trim(sevt1), &
        status='replace',action='write',access='stream',form='unformatted')
      write(747) phi(nx/2,:,:)
      close(747)
      
      !u
      sevt1 = trim(dir_out) // 'out'//'_'//'uxy'// &
        '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=845,file=trim(sevt1), &
        status='replace',action='write',access='stream',form='unformatted')
      write(845) u(:,:,nz/2)
      close(845)
      sevt1 = trim(dir_out) // 'out'//'_'//'uxz'// &
        '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=846,file=trim(sevt1), &
        status='replace',action='write',access='stream',form='unformatted')
      write(846) u(:,ny/2,:)
      close(846)
        sevt1 = trim(dir_out) // 'out'//'_'//'uyz'// &
        '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=847,file=trim(sevt1), &
        status='replace',action='write',access='stream',form='unformatted')
      write(847) u(nx/2,:,:)
      close(847)

      !v
      sevt1 = trim(dir_out) // 'out'//'_'//'vxy'// &
        '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=848,file=trim(sevt1), &
        status='replace',action='write',access='stream',form='unformatted')
      write(848) v(:,:,nz/2)
      close(848)
      sevt1 = trim(dir_out) // 'out'//'_'//'vyz'// &
        '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=849,file=trim(sevt1), &
        status='replace',action='write',access='stream',form='unformatted')
      write(849) v(nx/2,:,:)
      close(849)
      sevt1 = trim(dir_out) // 'out'//'_'//'vxz'// &
        '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=850,file=trim(sevt1), &
        status='replace',action='write',access='stream',form='unformatted')
      write(850) v(:,ny/2,:)
      close(850)

      !w
      sevt1 = trim(dir_out) // 'out'//'_'//'wxz'// &
        '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=851,file=trim(sevt1), &
        status='replace',action='write',access='stream',form='unformatted')
      write(851) w(:,ny/2,:)
      close(851)
      sevt1 = trim(dir_out) // 'out'//'_'//'wyz'// &
        '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=852,file=trim(sevt1), &
        status='replace',action='write',access='stream',form='unformatted')
      write(852) w(nx/2,:,:)
      close(852)
      sevt1 = trim(dir_out) // 'out'//'_'//'wxy'// &
        '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=853,file=trim(sevt1), &
        status='replace',action='write',access='stream',form='unformatted')
      write(853) w(:,:,nz/2)
      close(853)
   
  end subroutine print_raw_slice_2c_sync
  
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

  subroutine print_vtk_phionly_sync(iframe)
    implicit none

    integer, intent(in) :: iframe
    sevt1 = trim(dir_out) // 'out'//'_'//'phi'// &
    '_'//trim(write_fmtnumb(iframe)) // '.vti'
    open(unit=345,file=trim(sevt1), &
    status='replace',action='write',access='stream',form='unformatted')
    write(345)head1,ndatavtk(1),phi,footervtk(1)
    close(345)
  endsubroutine 
  
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
  
  subroutine close_print_async(lvtk)
   
   implicit none
   logical, intent(in) :: lvtk
   
   wait(345)
   if(lvtk)write(345)footervtk(1)
   close(345)
   
   
   wait(780)
   if(lvtk)write(780)footervtk(2)
   close(780) 
   
  end subroutine close_print_async
 endmodule
