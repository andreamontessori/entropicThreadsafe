module bcs3D
    
    use vars
    use mpi_template
    !$if _OPENACC
    use openacc
    !$endif
    implicit none
        
    contains
    subroutine bcs_all_bback_2c
        !$acc kernels
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
							!gggggggggggggggggggggggggg
							g(i,j,k-1,6)=g(i,j,k,5)!gpc 
                            g(i,j,k+1,5)=g(i,j,k,6)!hpc 

                            g(i,j-1,k,4)=g(i,j,k,3)!gpc 
                            g(i,j+1,k,3)=g(i,j,k,4)!hpc 

                            g(i-1,j,k,2)=g(i,j,k,1)!gpc 
                            g(i+1,j,k,1)=g(i,j,k,2)!hpc 
                        endif
                    enddo
                enddo
            enddo    
		!$acc end kernels      
    endsubroutine
    !***************************************************
    subroutine bcs_poiseuille_w_bback
        !$acc kernels
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
        !*********************************** call other bcs:PERIODIC ************************     
            !periodic along x 
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
            !periodic along y
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
    endsubroutine
    !*****************************************************
    subroutine bcs_pois_w_regularized_non_eq_extrapolation
        !$acc kernels
        !$acc loop independent 
        do k=1,nz
            !$acc loop independent 
            do j=1,ny
                !$acc loop independent 
                do i=1,nx
                    if(isfluid(i,j,k).eq.0)then
                        if(k.lt.nz .and. k.gt.1 )then
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
                        if(k.eq.nz)then
                            
                            !6
                            feq=(rho(i,j,k-1)*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                            fneq1=(3*(-pyy(i,j,k-1) + 2*pzz(i,j,k-1) + 6*pxz(i,j,k-1)*u(i,j,k-1) + 6*pyz(i,j,k-1)*v(i,j,k-1) + 3*pyy(i,j,k-1)*w(i,j,k-1) + pxx(i,j,k-1)*(-1 + 3*w(i,j,k-1))))/2.
                            f(i,j,k-1,6)=feq + (1-omega)*fneq1*p1
                            
                            !16
                            feq=(rho(i-1,j,k-1)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*(-pyy(i-1,j,k-1) + 2*pzz(i-1,j,k-1) + 3*pyy(i-1,j,k-1)*u(i-1,j,k-1) - 6*pzz(i-1,j,k-1)*u(i-1,j,k-1) + 6*pxy(i-1,j,k-1)*v(i-1,j,k-1) + 6*pyz(i-1,j,k-1)*v(i-1,j,k-1) + pxx(i-1,j,k-1)*(2 - 6*w(i-1,j,k-1)) + 3*pyy(i-1,j,k-1)*w(i-1,j,k-1) - 6*pxz(i-1,j,k-1)*(-1 + 2*u(i-1,j,k-1) + 2*w(i-1,j,k-1))))/2.
                            f(i-1,j,k-1,16)=feq + (1-omega)*fneq1*p2
                            
                            !18
                            feq=(rho(i+1,j,k-1)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(-3*(pyy(i+1,j,k-1) - 2*pzz(i+1,j,k-1) + 3*pyy(i+1,j,k-1)*u(i+1,j,k-1) - 6*pzz(i+1,j,k-1)*u(i+1,j,k-1) + 6*pxy(i+1,j,k-1)*v(i+1,j,k-1) - 6*pyz(i+1,j,k-1)*v(i+1,j,k-1) + 6*pxz(i+1,j,k-1)*(1 + 2*u(i+1,j,k-1) - 2*w(i+1,j,k-1)) - 3*pyy(i+1,j,k-1)*w(i+1,j,k-1) + pxx(i+1,j,k-1)*(-2 + 6*w(i+1,j,k-1))))/2.
                            f(i+1,j,k-1,18)=feq + (1-omega)*fneq1*p2
                            
                            !12
                            feq=(rho(i,j-1,k-1)*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*pxx(i,j-1,k-1)*(-1 + 3*v(i,j-1,k-1) + 3*w(i,j-1,k-1)))/2. + 3*(pyy(i,j-1,k-1) + pzz(i,j-1,k-1) + 3*pxy(i,j-1,k-1)*u(i,j-1,k-1) + 3*pxz(i,j-1,k-1)*u(i,j-1,k-1) - 3*pzz(i,j-1,k-1)*v(i,j-1,k-1) + pyz(i,j-1,k-1)*(3 - 6*v(i,j-1,k-1) - 6*w(i,j-1,k-1)) - 3*pyy(i,j-1,k-1)*w(i,j-1,k-1))
                            f(i,j-1,k-1,12)=feq + (1-omega)*p2*fneq1

                            !13
                            
                            feq=(rho(i,j+1,k-1)*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                            fneq1=(-3*(pxx(i,j+1,k-1) - 2*pyy(i,j+1,k-1) + 6*pyz(i,j+1,k-1) - 2*pzz(i,j+1,k-1) + 6*pxy(i,j+1,k-1)*u(i,j+1,k-1) - 6*pxz(i,j+1,k-1)*u(i,j+1,k-1) + 3*pxx(i,j+1,k-1)*v(i,j+1,k-1) + 12*pyz(i,j+1,k-1)*v(i,j+1,k-1) - 6*pzz(i,j+1,k-1)*v(i,j+1,k-1) - 3*pxx(i,j+1,k-1)*w(i,j+1,k-1) + 6*pyy(i,j+1,k-1)*w(i,j+1,k-1) - 12*pyz(i,j+1,k-1)*w(i,j+1,k-1)))/2.
                            f(i,j+1,k-1,13)=feq + (1-omega)*p2*fneq1
                            
                            !20
                            feq=(rho(i-1,j-1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*((pyy(i-1,j-1,k-1) + 3*pyz(i-1,j-1,k-1) + pzz(i-1,j-1,k-1))*(-1 + 3*u(i-1,j-1,k-1)) + pxy(i-1,j-1,k-1)*(-3 + 6*u(i-1,j-1,k-1)) + pxz(i-1,j-1,k-1)*(-3 + 6*u(i-1,j-1,k-1)) + 3*(2*pxy(i-1,j-1,k-1) + 3*pxz(i-1,j-1,k-1) + 2*pyz(i-1,j-1,k-1) + pzz(i-1,j-1,k-1))*v(i-1,j-1,k-1) + 3*(3*pxy(i-1,j-1,k-1) + &
                            2*pxz(i-1,j-1,k-1) + pyy(i-1,j-1,k-1) + 2*pyz(i-1,j-1,k-1))*w(i-1,j-1,k-1) + pxx(i-1,j-1,k-1)*(-1 + 3*v(i-1,j-1,k-1) + 3*w(i-1,j-1,k-1)))
                            f(i-1,j-1,k-1,20)=feq + (1-omega)*p3*fneq1
                            
                            !22
                            feq=(rho(i-1,j+1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i-1,j+1,k-1) + 3*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 3*pyz(i-1,j+1,k-1) + pzz(i-1,j+1,k-1) - 3*(2*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 3*pyz(i-1,j+1,k-1) + pzz(i-1,j+1,k-1))*u(i-1,j+1,k-1) + pxy(i-1,j+1,k-1)*(-3 + 6*u(i-1,j+1,k-1)) + 3*pxx(i-1,j+1,k-1)*v(i-1,j+1,k-1) &
                            - 6*pxy(i-1,j+1,k-1)*v(i-1,j+1,k-1) + 9*pxz(i-1,j+1,k-1)*v(i-1,j+1,k-1) - 6*pyz(i-1,j+1,k-1)*v(i-1,j+1,k-1) + 3*pzz(i-1,j+1,k-1)*v(i-1,j+1,k-1) - 3*(pxx(i-1,j+1,k-1) - 3*pxy(i-1,j+1,k-1) + 2*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 2*pyz(i-1,j+1,k-1))*w(i-1,j+1,k-1))

                            f(i-1,j+1,k-1,22)=feq + (1-omega)*p3*fneq1

                            !24
                            feq=(rho(i+1,j+1,k-1)*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                            fneq1=3*(pxx(i+1,j+1,k-1) - 3*pxz(i+1,j+1,k-1)*(1 + 2*u(i+1,j+1,k-1)) + (pyy(i+1,j+1,k-1) - 3*pyz(i+1,j+1,k-1) + pzz(i+1,j+1,k-1))*(1 + 3*u(i+1,j+1,k-1)) + pxy(i+1,j+1,k-1)*(3 + 6*u(i+1,j+1,k-1)) + 3*pxx(i+1,j+1,k-1)*v(i+1,j+1,k-1) + 6*pxy(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 9*pxz(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 6*pyz(i+1,j+1,k-1)*v(i+1,j+1,k-1) &
                            + 3*pzz(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 3*(pxx(i+1,j+1,k-1) + 3*pxy(i+1,j+1,k-1) - 2*pxz(i+1,j+1,k-1) + pyy(i+1,j+1,k-1) - 2*pyz(i+1,j+1,k-1))*w(i+1,j+1,k-1))
                            f(i+1,j+1,k-1,24)=feq + (1-omega)*p3*fneq1

                            !25
                            feq=(rho(i+1,j-1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*(-pxx(i+1,j-1,k-1) - (pyy(i+1,j-1,k-1) + 3*pyz(i+1,j-1,k-1) + pzz(i+1,j-1,k-1))*(1 + 3*u(i+1,j-1,k-1)) + pxy(i+1,j-1,k-1)*(3 + 6*u(i+1,j-1,k-1)) + pxz(i+1,j-1,k-1)*(3 + 6*u(i+1,j-1,k-1)) + 3*pxx(i+1,j-1,k-1)*v(i+1,j-1,k-1) &
                            - 6*pxy(i+1,j-1,k-1)*v(i+1,j-1,k-1) - 9*pxz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 6*pyz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 3*pzz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 3*(pxx(i+1,j-1,k-1) - 3*pxy(i+1,j-1,k-1) - 2*pxz(i+1,j-1,k-1) + pyy(i+1,j-1,k-1) + 2*pyz(i+1,j-1,k-1))*w(i+1,j-1,k-1))
                            f(i+1,j-1,k-1,25)=feq + (1-omega)*p3*fneq1

                                
                        elseif(k.eq.1)then


                            !5
                            
                            feq=(rho(i,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k))))/27.
                            fneq1=(-3*(pxx(i,j,k+1) + pyy(i,j,k+1) - 2*pzz(i,j,k+1) + 6*pxz(i,j,k+1)*u(i,j,k+1) + 6*pyz(i,j,k+1)*v(i,j,k+1) + 3*pxx(i,j,k+1)*w(i,j,k+1) + 3*pyy(i,j,k+1)*w(i,j,k+1)))/2.
                            f(i,j,k+1,5)=feq + (1-omega)*p1*fneq1

                            !15
                            
                            feq=(rho(i+1,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i+1,j,k+1) + 2*pzz(i+1,j,k+1) - 3*pyy(i+1,j,k+1)*u(i+1,j,k+1) + 6*pzz(i+1,j,k+1)*u(i+1,j,k+1) - 6*pxy(i+1,j,k+1)*v(i+1,j,k+1) - 6*pyz(i+1,j,k+1)*v(i+1,j,k+1) - 3*pyy(i+1,j,k+1)*w(i+1,j,k+1) + 6*pxz(i+1,j,k+1)*(1 + 2*u(i+1,j,k+1) + 2*w(i+1,j,k+1)) + pxx(i+1,j,k+1)*(2 + 6*w(i+1,j,k+1))))/2.
                            f(i+1,j,k+1,15)=feq + (1-omega)*p2*fneq1

                            !17
                            
                            feq=(rho(i-1,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i-1,j,k+1) + 2*pzz(i-1,j,k+1) + 3*pyy(i-1,j,k+1)*u(i-1,j,k+1) - 6*pzz(i-1,j,k+1)*u(i-1,j,k+1) + 6*pxy(i-1,j,k+1)*v(i-1,j,k+1) - 6*pyz(i-1,j,k+1)*v(i-1,j,k+1) + 6*pxz(i-1,j,k+1)*(-1 + 2*u(i-1,j,k+1) - 2*w(i-1,j,k+1)) - 3*pyy(i-1,j,k+1)*w(i-1,j,k+1) + pxx(i-1,j,k+1)*(2 + 6*w(i-1,j,k+1))))/2.
                            f(i-1,j,k+1,17)=feq + (1-omega)*p2*fneq1

                            !11
                            
                            feq=(rho(i,j+1,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(-3*pxx(i,j+1,k+1)*(1 + 3*v(i,j+1,k+1) + 3*w(i,j+1,k+1)))/2. + 3*(pyy(i,j+1,k+1) + pzz(i,j+1,k+1) - 3*pxy(i,j+1,k+1)*u(i,j+1,k+1) - 3*pxz(i,j+1,k+1)*u(i,j+1,k+1) + 3*pzz(i,j+1,k+1)*v(i,j+1,k+1) + 3*pyy(i,j+1,k+1)*w(i,j+1,k+1) + pyz(i,j+1,k+1)*(3 + 6*v(i,j+1,k+1) + 6*w(i,j+1,k+1)))
                            f(i,j+1,k+1,11)=feq +(1-omega)*p2*fneq1
                            
                            !14
                            feq=(rho(i,j-1,k+1)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*pxx(i,j-1,k+1)*(-1 + 3*v(i,j-1,k+1) - 3*w(i,j-1,k+1)))/2. + 3*(pyy(i,j-1,k+1) + pzz(i,j-1,k+1) + 3*pxy(i,j-1,k+1)*u(i,j-1,k+1) - 3*pxz(i,j-1,k+1)*u(i,j-1,k+1) - 3*pzz(i,j-1,k+1)*v(i,j-1,k+1) + pyz(i,j-1,k+1)*(-3 + 6*v(i,j-1,k+1) - 6*w(i,j-1,k+1)) + 3*pyy(i,j-1,k+1)*w(i,j-1,k+1))
                            f(i,j-1,k+1,14)=feq + (1-omega)*p2*fneq1

                            !19
                            feq=(rho(i+1,j+1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i+1,j+1,k+1) + (pyy(i+1,j+1,k+1) + 3*pyz(i+1,j+1,k+1) + pzz(i+1,j+1,k+1))*(1 + 3*u(i+1,j+1,k+1)) + pxy(i+1,j+1,k+1)*(3 + 6*u(i+1,j+1,k+1)) + pxz(i+1,j+1,k+1)*(3 + 6*u(i+1,j+1,k+1)) + 3*pxx(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 6*pxy(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 9*pxz(i+1,j+1,k+1)*v(i+1,j+1,k+1) &
                            + 6*pyz(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 3*pzz(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 3*(pxx(i+1,j+1,k+1) + 3*pxy(i+1,j+1,k+1) + 2*pxz(i+1,j+1,k+1) + pyy(i+1,j+1,k+1) + 2*pyz(i+1,j+1,k+1))*w(i+1,j+1,k+1))
                            f(i+1,j+1,k+1,19)=feq + (1-omega)*fneq1*p3

                            !21
                            feq=(rho(i+1,j-1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i+1,j-1,k+1) - 3*pxy(i+1,j-1,k+1)*(1 + 2*u(i+1,j-1,k+1)) + (pyy(i+1,j-1,k+1) - 3*pyz(i+1,j-1,k+1) + pzz(i+1,j-1,k+1))*(1 + 3*u(i+1,j-1,k+1)) + pxz(i+1,j-1,k+1)*(3 + 6*u(i+1,j-1,k+1)) - 3*pxx(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 6*pxy(i+1,j-1,k+1)*v(i+1,j-1,k+1) &
                            - 9*pxz(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 6*pyz(i+1,j-1,k+1)*v(i+1,j-1,k+1) - 3*pzz(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 3*(pxx(i+1,j-1,k+1) - 3*pxy(i+1,j-1,k+1) + 2*pxz(i+1,j-1,k+1) + pyy(i+1,j-1,k+1) - 2*pyz(i+1,j-1,k+1))*w(i+1,j-1,k+1))
                            f(i+1,j-1,k+1,21)=feq + (1-omega)*fneq1*p3

                            !23
                            feq=(rho(i-1,j-1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i-1,j-1,k+1) + 3*pxy(i-1,j-1,k+1) - 3*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 3*pyz(i-1,j-1,k+1) + pzz(i-1,j-1,k+1) - 3*(2*pxy(i-1,j-1,k+1) - 2*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 3*pyz(i-1,j-1,k+1) + pzz(i-1,j-1,k+1))*u(i-1,j-1,k+1) &
                            - 3*pxx(i-1,j-1,k+1)*v(i-1,j-1,k+1) - 6*pxy(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 9*pxz(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 6*pyz(i-1,j-1,k+1)*v(i-1,j-1,k+1) - 3*pzz(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 3*(pxx(i-1,j-1,k+1) + 3*pxy(i-1,j-1,k+1) - 2*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 2*pyz(i-1,j-1,k+1))*w(i-1,j-1,k+1))
                            f(i-1,j-1,k+1,23)=feq + (1-omega)*fneq1*p3

                            !26
                            feq=(rho(i-1,j+1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i-1,j+1,k+1) - (pyy(i-1,j+1,k+1) + 3*pyz(i-1,j+1,k+1) + pzz(i-1,j+1,k+1))*(-1 + 3*u(i-1,j+1,k+1)) + pxy(i-1,j+1,k+1)*(-3 + 6*u(i-1,j+1,k+1)) + pxz(i-1,j+1,k+1)*(-3 + 6*u(i-1,j+1,k+1)) + 3*pxx(i-1,j+1,k+1)*v(i-1,j+1,k+1) - 6*pxy(i-1,j+1,k+1)*v(i-1,j+1,k+1) &
                            - 9*pxz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 6*pyz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 3*pzz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 3*(pxx(i-1,j+1,k+1) - 3*pxy(i-1,j+1,k+1) - 2*pxz(i-1,j+1,k+1) + pyy(i-1,j+1,k+1) + 2*pyz(i-1,j+1,k+1))*w(i-1,j+1,k+1))
                            f(i-1,j+1,k+1,26)=feq + (1-omega)*fneq1*p3
                        endif
                    endif
                enddo
            enddo
        enddo
        !*********************************** call other bcs:PERIODIC ************************     
        !periodic along x 
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
        !periodic along y
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
    endsubroutine
    !*****************************************************
    subroutine bcs_couette_w_regularized_non_eq_extrapolation
        !$acc kernels
        !$acc loop independent 
        do k=1,nz
            !$acc loop independent 
            do j=1,ny
                !$acc loop independent 
                do i=1,nx
                    if(isfluid(i,j,k).eq.0)then
                        if(k.lt.nz .and. k.gt.1 )then
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
                        if(k.eq.nz)then
                            
                            u(i,j,k)=uwall
                            v(i,j,k)=0
                            w(i,j,k)=0
                            
                            !6
                            feq=(rho(i,j,k-1)*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                            fneq1=(3*(-pyy(i,j,k-1) + 2*pzz(i,j,k-1) + 6*pxz(i,j,k-1)*u(i,j,k-1) + 6*pyz(i,j,k-1)*v(i,j,k-1) + 3*pyy(i,j,k-1)*w(i,j,k-1) + pxx(i,j,k-1)*(-1 + 3*w(i,j,k-1))))/2.
                            f(i,j,k-1,6)=feq + (1-omega)*fneq1*p1
                            
                            !16
                            feq=(rho(i-1,j,k-1)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*(-pyy(i-1,j,k-1) + 2*pzz(i-1,j,k-1) + 3*pyy(i-1,j,k-1)*u(i-1,j,k-1) - 6*pzz(i-1,j,k-1)*u(i-1,j,k-1) + 6*pxy(i-1,j,k-1)*v(i-1,j,k-1) + 6*pyz(i-1,j,k-1)*v(i-1,j,k-1) + pxx(i-1,j,k-1)*(2 - 6*w(i-1,j,k-1)) + 3*pyy(i-1,j,k-1)*w(i-1,j,k-1) - 6*pxz(i-1,j,k-1)*(-1 + 2*u(i-1,j,k-1) + 2*w(i-1,j,k-1))))/2.
                            f(i-1,j,k-1,16)=feq + (1-omega)*fneq1*p2
                            
                            !18
                            feq=(rho(i+1,j,k-1)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(-3*(pyy(i+1,j,k-1) - 2*pzz(i+1,j,k-1) + 3*pyy(i+1,j,k-1)*u(i+1,j,k-1) - 6*pzz(i+1,j,k-1)*u(i+1,j,k-1) + 6*pxy(i+1,j,k-1)*v(i+1,j,k-1) - 6*pyz(i+1,j,k-1)*v(i+1,j,k-1) + 6*pxz(i+1,j,k-1)*(1 + 2*u(i+1,j,k-1) - 2*w(i+1,j,k-1)) - 3*pyy(i+1,j,k-1)*w(i+1,j,k-1) + pxx(i+1,j,k-1)*(-2 + 6*w(i+1,j,k-1))))/2.
                            f(i+1,j,k-1,18)=feq + (1-omega)*fneq1*p2
                            
                            !12
                            feq=(rho(i,j-1,k-1)*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*pxx(i,j-1,k-1)*(-1 + 3*v(i,j-1,k-1) + 3*w(i,j-1,k-1)))/2. + 3*(pyy(i,j-1,k-1) + pzz(i,j-1,k-1) + 3*pxy(i,j-1,k-1)*u(i,j-1,k-1) + 3*pxz(i,j-1,k-1)*u(i,j-1,k-1) - 3*pzz(i,j-1,k-1)*v(i,j-1,k-1) + pyz(i,j-1,k-1)*(3 - 6*v(i,j-1,k-1) - 6*w(i,j-1,k-1)) - 3*pyy(i,j-1,k-1)*w(i,j-1,k-1))
                            f(i,j-1,k-1,12)=feq + (1-omega)*p2*fneq1

                            !13
                            
                            feq=(rho(i,j+1,k-1)*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                            fneq1=(-3*(pxx(i,j+1,k-1) - 2*pyy(i,j+1,k-1) + 6*pyz(i,j+1,k-1) - 2*pzz(i,j+1,k-1) + 6*pxy(i,j+1,k-1)*u(i,j+1,k-1) - 6*pxz(i,j+1,k-1)*u(i,j+1,k-1) + 3*pxx(i,j+1,k-1)*v(i,j+1,k-1) + 12*pyz(i,j+1,k-1)*v(i,j+1,k-1) - 6*pzz(i,j+1,k-1)*v(i,j+1,k-1) - 3*pxx(i,j+1,k-1)*w(i,j+1,k-1) + 6*pyy(i,j+1,k-1)*w(i,j+1,k-1) - 12*pyz(i,j+1,k-1)*w(i,j+1,k-1)))/2.
                            f(i,j+1,k-1,13)=feq + (1-omega)*p2*fneq1
                            
                            !20
                            feq=(rho(i-1,j-1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*((pyy(i-1,j-1,k-1) + 3*pyz(i-1,j-1,k-1) + pzz(i-1,j-1,k-1))*(-1 + 3*u(i-1,j-1,k-1)) + pxy(i-1,j-1,k-1)*(-3 + 6*u(i-1,j-1,k-1)) + pxz(i-1,j-1,k-1)*(-3 + 6*u(i-1,j-1,k-1)) + 3*(2*pxy(i-1,j-1,k-1) + 3*pxz(i-1,j-1,k-1) + 2*pyz(i-1,j-1,k-1) + pzz(i-1,j-1,k-1))*v(i-1,j-1,k-1) + 3*(3*pxy(i-1,j-1,k-1) + &
                            2*pxz(i-1,j-1,k-1) + pyy(i-1,j-1,k-1) + 2*pyz(i-1,j-1,k-1))*w(i-1,j-1,k-1) + pxx(i-1,j-1,k-1)*(-1 + 3*v(i-1,j-1,k-1) + 3*w(i-1,j-1,k-1)))
                            f(i-1,j-1,k-1,20)=feq + (1-omega)*p3*fneq1
                            
                            !22
                            feq=(rho(i-1,j+1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i-1,j+1,k-1) + 3*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 3*pyz(i-1,j+1,k-1) + pzz(i-1,j+1,k-1) - 3*(2*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 3*pyz(i-1,j+1,k-1) + pzz(i-1,j+1,k-1))*u(i-1,j+1,k-1) + pxy(i-1,j+1,k-1)*(-3 + 6*u(i-1,j+1,k-1)) + 3*pxx(i-1,j+1,k-1)*v(i-1,j+1,k-1) &
                            - 6*pxy(i-1,j+1,k-1)*v(i-1,j+1,k-1) + 9*pxz(i-1,j+1,k-1)*v(i-1,j+1,k-1) - 6*pyz(i-1,j+1,k-1)*v(i-1,j+1,k-1) + 3*pzz(i-1,j+1,k-1)*v(i-1,j+1,k-1) - 3*(pxx(i-1,j+1,k-1) - 3*pxy(i-1,j+1,k-1) + 2*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 2*pyz(i-1,j+1,k-1))*w(i-1,j+1,k-1))

                            f(i-1,j+1,k-1,22)=feq + (1-omega)*p3*fneq1

                            !24
                            feq=(rho(i+1,j+1,k-1)*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                            fneq1=3*(pxx(i+1,j+1,k-1) - 3*pxz(i+1,j+1,k-1)*(1 + 2*u(i+1,j+1,k-1)) + (pyy(i+1,j+1,k-1) - 3*pyz(i+1,j+1,k-1) + pzz(i+1,j+1,k-1))*(1 + 3*u(i+1,j+1,k-1)) + pxy(i+1,j+1,k-1)*(3 + 6*u(i+1,j+1,k-1)) + 3*pxx(i+1,j+1,k-1)*v(i+1,j+1,k-1) + 6*pxy(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 9*pxz(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 6*pyz(i+1,j+1,k-1)*v(i+1,j+1,k-1) &
                            + 3*pzz(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 3*(pxx(i+1,j+1,k-1) + 3*pxy(i+1,j+1,k-1) - 2*pxz(i+1,j+1,k-1) + pyy(i+1,j+1,k-1) - 2*pyz(i+1,j+1,k-1))*w(i+1,j+1,k-1))
                            f(i+1,j+1,k-1,24)=feq + (1-omega)*p3*fneq1

                            !25
                            feq=(rho(i+1,j-1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*(-pxx(i+1,j-1,k-1) - (pyy(i+1,j-1,k-1) + 3*pyz(i+1,j-1,k-1) + pzz(i+1,j-1,k-1))*(1 + 3*u(i+1,j-1,k-1)) + pxy(i+1,j-1,k-1)*(3 + 6*u(i+1,j-1,k-1)) + pxz(i+1,j-1,k-1)*(3 + 6*u(i+1,j-1,k-1)) + 3*pxx(i+1,j-1,k-1)*v(i+1,j-1,k-1) &
                            - 6*pxy(i+1,j-1,k-1)*v(i+1,j-1,k-1) - 9*pxz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 6*pyz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 3*pzz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 3*(pxx(i+1,j-1,k-1) - 3*pxy(i+1,j-1,k-1) - 2*pxz(i+1,j-1,k-1) + pyy(i+1,j-1,k-1) + 2*pyz(i+1,j-1,k-1))*w(i+1,j-1,k-1))
                            f(i+1,j-1,k-1,25)=feq + (1-omega)*p3*fneq1

                                
                        elseif(k.eq.1)then
                                
                            u(i,j,k)=-uwall
                            v(i,j,k)=0
                            w(i,j,k)=0

                            !5
                            
                            feq=(rho(i,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k))))/27.
                            fneq1=(-3*(pxx(i,j,k+1) + pyy(i,j,k+1) - 2*pzz(i,j,k+1) + 6*pxz(i,j,k+1)*u(i,j,k+1) + 6*pyz(i,j,k+1)*v(i,j,k+1) + 3*pxx(i,j,k+1)*w(i,j,k+1) + 3*pyy(i,j,k+1)*w(i,j,k+1)))/2.
                            f(i,j,k+1,5)=feq + (1-omega)*p1*fneq1

                            !15
                            
                            feq=(rho(i+1,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i+1,j,k+1) + 2*pzz(i+1,j,k+1) - 3*pyy(i+1,j,k+1)*u(i+1,j,k+1) + 6*pzz(i+1,j,k+1)*u(i+1,j,k+1) - 6*pxy(i+1,j,k+1)*v(i+1,j,k+1) - 6*pyz(i+1,j,k+1)*v(i+1,j,k+1) - 3*pyy(i+1,j,k+1)*w(i+1,j,k+1) + 6*pxz(i+1,j,k+1)*(1 + 2*u(i+1,j,k+1) + 2*w(i+1,j,k+1)) + pxx(i+1,j,k+1)*(2 + 6*w(i+1,j,k+1))))/2.
                            f(i+1,j,k+1,15)=feq + (1-omega)*p2*fneq1

                            !17
                            
                            feq=(rho(i-1,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i-1,j,k+1) + 2*pzz(i-1,j,k+1) + 3*pyy(i-1,j,k+1)*u(i-1,j,k+1) - 6*pzz(i-1,j,k+1)*u(i-1,j,k+1) + 6*pxy(i-1,j,k+1)*v(i-1,j,k+1) - 6*pyz(i-1,j,k+1)*v(i-1,j,k+1) + 6*pxz(i-1,j,k+1)*(-1 + 2*u(i-1,j,k+1) - 2*w(i-1,j,k+1)) - 3*pyy(i-1,j,k+1)*w(i-1,j,k+1) + pxx(i-1,j,k+1)*(2 + 6*w(i-1,j,k+1))))/2.
                            f(i-1,j,k+1,17)=feq + (1-omega)*p2*fneq1

                            !11
                            
                            feq=(rho(i,j+1,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(-3*pxx(i,j+1,k+1)*(1 + 3*v(i,j+1,k+1) + 3*w(i,j+1,k+1)))/2. + 3*(pyy(i,j+1,k+1) + pzz(i,j+1,k+1) - 3*pxy(i,j+1,k+1)*u(i,j+1,k+1) - 3*pxz(i,j+1,k+1)*u(i,j+1,k+1) + 3*pzz(i,j+1,k+1)*v(i,j+1,k+1) + 3*pyy(i,j+1,k+1)*w(i,j+1,k+1) + pyz(i,j+1,k+1)*(3 + 6*v(i,j+1,k+1) + 6*w(i,j+1,k+1)))
                            f(i,j+1,k+1,11)=feq +(1-omega)*p2*fneq1
                            
                            !14
                            feq=(rho(i,j-1,k+1)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*pxx(i,j-1,k+1)*(-1 + 3*v(i,j-1,k+1) - 3*w(i,j-1,k+1)))/2. + 3*(pyy(i,j-1,k+1) + pzz(i,j-1,k+1) + 3*pxy(i,j-1,k+1)*u(i,j-1,k+1) - 3*pxz(i,j-1,k+1)*u(i,j-1,k+1) - 3*pzz(i,j-1,k+1)*v(i,j-1,k+1) + pyz(i,j-1,k+1)*(-3 + 6*v(i,j-1,k+1) - 6*w(i,j-1,k+1)) + 3*pyy(i,j-1,k+1)*w(i,j-1,k+1))
                            f(i,j-1,k+1,14)=feq + (1-omega)*p2*fneq1

                            !19
                            feq=(rho(i+1,j+1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i+1,j+1,k+1) + (pyy(i+1,j+1,k+1) + 3*pyz(i+1,j+1,k+1) + pzz(i+1,j+1,k+1))*(1 + 3*u(i+1,j+1,k+1)) + pxy(i+1,j+1,k+1)*(3 + 6*u(i+1,j+1,k+1)) + pxz(i+1,j+1,k+1)*(3 + 6*u(i+1,j+1,k+1)) + 3*pxx(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 6*pxy(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 9*pxz(i+1,j+1,k+1)*v(i+1,j+1,k+1) &
                            + 6*pyz(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 3*pzz(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 3*(pxx(i+1,j+1,k+1) + 3*pxy(i+1,j+1,k+1) + 2*pxz(i+1,j+1,k+1) + pyy(i+1,j+1,k+1) + 2*pyz(i+1,j+1,k+1))*w(i+1,j+1,k+1))
                            f(i+1,j+1,k+1,19)=feq + (1-omega)*fneq1*p3

                            !21
                            feq=(rho(i+1,j-1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i+1,j-1,k+1) - 3*pxy(i+1,j-1,k+1)*(1 + 2*u(i+1,j-1,k+1)) + (pyy(i+1,j-1,k+1) - 3*pyz(i+1,j-1,k+1) + pzz(i+1,j-1,k+1))*(1 + 3*u(i+1,j-1,k+1)) + pxz(i+1,j-1,k+1)*(3 + 6*u(i+1,j-1,k+1)) - 3*pxx(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 6*pxy(i+1,j-1,k+1)*v(i+1,j-1,k+1) &
                            - 9*pxz(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 6*pyz(i+1,j-1,k+1)*v(i+1,j-1,k+1) - 3*pzz(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 3*(pxx(i+1,j-1,k+1) - 3*pxy(i+1,j-1,k+1) + 2*pxz(i+1,j-1,k+1) + pyy(i+1,j-1,k+1) - 2*pyz(i+1,j-1,k+1))*w(i+1,j-1,k+1))
                            f(i+1,j-1,k+1,21)=feq + (1-omega)*fneq1*p3

                            !23
                            feq=(rho(i-1,j-1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i-1,j-1,k+1) + 3*pxy(i-1,j-1,k+1) - 3*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 3*pyz(i-1,j-1,k+1) + pzz(i-1,j-1,k+1) - 3*(2*pxy(i-1,j-1,k+1) - 2*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 3*pyz(i-1,j-1,k+1) + pzz(i-1,j-1,k+1))*u(i-1,j-1,k+1) &
                            - 3*pxx(i-1,j-1,k+1)*v(i-1,j-1,k+1) - 6*pxy(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 9*pxz(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 6*pyz(i-1,j-1,k+1)*v(i-1,j-1,k+1) - 3*pzz(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 3*(pxx(i-1,j-1,k+1) + 3*pxy(i-1,j-1,k+1) - 2*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 2*pyz(i-1,j-1,k+1))*w(i-1,j-1,k+1))
                            f(i-1,j-1,k+1,23)=feq + (1-omega)*fneq1*p3

                            !26
                            feq=(rho(i-1,j+1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i-1,j+1,k+1) - (pyy(i-1,j+1,k+1) + 3*pyz(i-1,j+1,k+1) + pzz(i-1,j+1,k+1))*(-1 + 3*u(i-1,j+1,k+1)) + pxy(i-1,j+1,k+1)*(-3 + 6*u(i-1,j+1,k+1)) + pxz(i-1,j+1,k+1)*(-3 + 6*u(i-1,j+1,k+1)) + 3*pxx(i-1,j+1,k+1)*v(i-1,j+1,k+1) - 6*pxy(i-1,j+1,k+1)*v(i-1,j+1,k+1) &
                            - 9*pxz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 6*pyz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 3*pzz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 3*(pxx(i-1,j+1,k+1) - 3*pxy(i-1,j+1,k+1) - 2*pxz(i-1,j+1,k+1) + pyy(i-1,j+1,k+1) + 2*pyz(i-1,j+1,k+1))*w(i-1,j+1,k+1))
                            f(i-1,j+1,k+1,26)=feq + (1-omega)*fneq1*p3
                        endif
                    endif
                enddo
            enddo
        enddo
        !*********************************** call other bcs:PERIODIC ************************     
        !periodic along x 
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
        !periodic along y
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
    endsubroutine
	!*****************************************************
	subroutine bcs_turbulent_jet
	
	    implicit none
		
		if(mod(step,10).eq.0)then
			!$acc update host(w(nx/2-20:nx/2+20,ny/2-20:ny/2+20,1))
			do i=nx/2-20,nx/2+20
				do j=ny/2-20,ny/2+20
					if((float(i)-nx/2.0)**2 + (float(j)-ny/2.0)**2<=10**2)then
						call random_number(rrx)
						w(i,j,1)=uwall + 0.004*sqrt(-2.0*log(rrx))*cos(2*3.1415926535897932384626433832795028841971*rrx)
					endif
				enddo
			enddo
		   !$acc update device(w(nx/2-20:nx/2+20,ny/2-20:ny/2+20,1))
		endif
		!$acc kernels
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
						if(k.eq.1)then
							!5
							!w(i,j,k)=uwall + 0.02*sqrt(-2.0*log(ranx(ttt)))*cos(2*3.1415926535897932384626433832795028841971*rany(ttt))
                            feq=(rho(i,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k))))/27.
                            fneq1=(-3*(pxx(i,j,k+1) + pyy(i,j,k+1) - 2*pzz(i,j,k+1) + 6*pxz(i,j,k+1)*u(i,j,k+1) + 6*pyz(i,j,k+1)*v(i,j,k+1) + 3*pxx(i,j,k+1)*w(i,j,k+1) + 3*pyy(i,j,k+1)*w(i,j,k+1)))/2.
                            f(i,j,k+1,5)=feq + (1-omega)*p1*fneq1

                            !15
                            
                            feq=(rho(i+1,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i+1,j,k+1) + 2*pzz(i+1,j,k+1) - 3*pyy(i+1,j,k+1)*u(i+1,j,k+1) + 6*pzz(i+1,j,k+1)*u(i+1,j,k+1) - 6*pxy(i+1,j,k+1)*v(i+1,j,k+1) - 6*pyz(i+1,j,k+1)*v(i+1,j,k+1) - 3*pyy(i+1,j,k+1)*w(i+1,j,k+1) + 6*pxz(i+1,j,k+1)*(1 + 2*u(i+1,j,k+1) + 2*w(i+1,j,k+1)) + pxx(i+1,j,k+1)*(2 + 6*w(i+1,j,k+1))))/2.
                            f(i+1,j,k+1,15)=feq + (1-omega)*p2*fneq1

                            !17
                            
                            feq=(rho(i-1,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i-1,j,k+1) + 2*pzz(i-1,j,k+1) + 3*pyy(i-1,j,k+1)*u(i-1,j,k+1) - 6*pzz(i-1,j,k+1)*u(i-1,j,k+1) + 6*pxy(i-1,j,k+1)*v(i-1,j,k+1) - 6*pyz(i-1,j,k+1)*v(i-1,j,k+1) + 6*pxz(i-1,j,k+1)*(-1 + 2*u(i-1,j,k+1) - 2*w(i-1,j,k+1)) - 3*pyy(i-1,j,k+1)*w(i-1,j,k+1) + pxx(i-1,j,k+1)*(2 + 6*w(i-1,j,k+1))))/2.
                            f(i-1,j,k+1,17)=feq + (1-omega)*p2*fneq1

                            !11
                            
                            feq=(rho(i,j+1,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(-3*pxx(i,j+1,k+1)*(1 + 3*v(i,j+1,k+1) + 3*w(i,j+1,k+1)))/2. + 3*(pyy(i,j+1,k+1) + pzz(i,j+1,k+1) - 3*pxy(i,j+1,k+1)*u(i,j+1,k+1) - 3*pxz(i,j+1,k+1)*u(i,j+1,k+1) + 3*pzz(i,j+1,k+1)*v(i,j+1,k+1) + 3*pyy(i,j+1,k+1)*w(i,j+1,k+1) + pyz(i,j+1,k+1)*(3 + 6*v(i,j+1,k+1) + 6*w(i,j+1,k+1)))
                            f(i,j+1,k+1,11)=feq +(1-omega)*p2*fneq1
                            
                            !14
                            feq=(rho(i,j-1,k+1)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*pxx(i,j-1,k+1)*(-1 + 3*v(i,j-1,k+1) - 3*w(i,j-1,k+1)))/2. + 3*(pyy(i,j-1,k+1) + pzz(i,j-1,k+1) + 3*pxy(i,j-1,k+1)*u(i,j-1,k+1) - 3*pxz(i,j-1,k+1)*u(i,j-1,k+1) - 3*pzz(i,j-1,k+1)*v(i,j-1,k+1) + pyz(i,j-1,k+1)*(-3 + 6*v(i,j-1,k+1) - 6*w(i,j-1,k+1)) + 3*pyy(i,j-1,k+1)*w(i,j-1,k+1))
                            f(i,j-1,k+1,14)=feq + (1-omega)*p2*fneq1

                            !19
                            feq=(rho(i+1,j+1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i+1,j+1,k+1) + (pyy(i+1,j+1,k+1) + 3*pyz(i+1,j+1,k+1) + pzz(i+1,j+1,k+1))*(1 + 3*u(i+1,j+1,k+1)) + pxy(i+1,j+1,k+1)*(3 + 6*u(i+1,j+1,k+1)) + pxz(i+1,j+1,k+1)*(3 + 6*u(i+1,j+1,k+1)) + 3*pxx(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 6*pxy(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 9*pxz(i+1,j+1,k+1)*v(i+1,j+1,k+1) &
                            + 6*pyz(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 3*pzz(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 3*(pxx(i+1,j+1,k+1) + 3*pxy(i+1,j+1,k+1) + 2*pxz(i+1,j+1,k+1) + pyy(i+1,j+1,k+1) + 2*pyz(i+1,j+1,k+1))*w(i+1,j+1,k+1))
                            f(i+1,j+1,k+1,19)=feq + (1-omega)*fneq1*p3

                            !21
                            feq=(rho(i+1,j-1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i+1,j-1,k+1) - 3*pxy(i+1,j-1,k+1)*(1 + 2*u(i+1,j-1,k+1)) + (pyy(i+1,j-1,k+1) - 3*pyz(i+1,j-1,k+1) + pzz(i+1,j-1,k+1))*(1 + 3*u(i+1,j-1,k+1)) + pxz(i+1,j-1,k+1)*(3 + 6*u(i+1,j-1,k+1)) - 3*pxx(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 6*pxy(i+1,j-1,k+1)*v(i+1,j-1,k+1) &
                            - 9*pxz(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 6*pyz(i+1,j-1,k+1)*v(i+1,j-1,k+1) - 3*pzz(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 3*(pxx(i+1,j-1,k+1) - 3*pxy(i+1,j-1,k+1) + 2*pxz(i+1,j-1,k+1) + pyy(i+1,j-1,k+1) - 2*pyz(i+1,j-1,k+1))*w(i+1,j-1,k+1))
                            f(i+1,j-1,k+1,21)=feq + (1-omega)*fneq1*p3

                            !23
                            feq=(rho(i-1,j-1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i-1,j-1,k+1) + 3*pxy(i-1,j-1,k+1) - 3*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 3*pyz(i-1,j-1,k+1) + pzz(i-1,j-1,k+1) - 3*(2*pxy(i-1,j-1,k+1) - 2*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 3*pyz(i-1,j-1,k+1) + pzz(i-1,j-1,k+1))*u(i-1,j-1,k+1) &
                            - 3*pxx(i-1,j-1,k+1)*v(i-1,j-1,k+1) - 6*pxy(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 9*pxz(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 6*pyz(i-1,j-1,k+1)*v(i-1,j-1,k+1) - 3*pzz(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 3*(pxx(i-1,j-1,k+1) + 3*pxy(i-1,j-1,k+1) - 2*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 2*pyz(i-1,j-1,k+1))*w(i-1,j-1,k+1))
                            f(i-1,j-1,k+1,23)=feq + (1-omega)*fneq1*p3

                            !26
                            feq=(rho(i-1,j+1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i-1,j+1,k+1) - (pyy(i-1,j+1,k+1) + 3*pyz(i-1,j+1,k+1) + pzz(i-1,j+1,k+1))*(-1 + 3*u(i-1,j+1,k+1)) + pxy(i-1,j+1,k+1)*(-3 + 6*u(i-1,j+1,k+1)) + pxz(i-1,j+1,k+1)*(-3 + 6*u(i-1,j+1,k+1)) + 3*pxx(i-1,j+1,k+1)*v(i-1,j+1,k+1) - 6*pxy(i-1,j+1,k+1)*v(i-1,j+1,k+1) &
                            - 9*pxz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 6*pyz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 3*pzz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 3*(pxx(i-1,j+1,k+1) - 3*pxy(i-1,j+1,k+1) - 2*pxz(i-1,j+1,k+1) + pyy(i-1,j+1,k+1) + 2*pyz(i-1,j+1,k+1))*w(i-1,j+1,k+1))
                            f(i-1,j+1,k+1,26)=feq + (1-omega)*fneq1*p3
						elseif(k.eq.nz)then
							!6
							w(i,j,k)=w(i,j,nz-1)
                            feq=(rho(i,j,k-1)*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                            fneq1=(3*(-pyy(i,j,k-1) + 2*pzz(i,j,k-1) + 6*pxz(i,j,k-1)*u(i,j,k-1) + 6*pyz(i,j,k-1)*v(i,j,k-1) + 3*pyy(i,j,k-1)*w(i,j,k-1) + pxx(i,j,k-1)*(-1 + 3*w(i,j,k-1))))/2.
                            f(i,j,k-1,6)=feq + (1-omega)*fneq1*p1
                            
                            !16
                            feq=(rho(i-1,j,k-1)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*(-pyy(i-1,j,k-1) + 2*pzz(i-1,j,k-1) + 3*pyy(i-1,j,k-1)*u(i-1,j,k-1) - 6*pzz(i-1,j,k-1)*u(i-1,j,k-1) + 6*pxy(i-1,j,k-1)*v(i-1,j,k-1) + 6*pyz(i-1,j,k-1)*v(i-1,j,k-1) + pxx(i-1,j,k-1)*(2 - 6*w(i-1,j,k-1)) + 3*pyy(i-1,j,k-1)*w(i-1,j,k-1) - 6*pxz(i-1,j,k-1)*(-1 + 2*u(i-1,j,k-1) + 2*w(i-1,j,k-1))))/2.
                            f(i-1,j,k-1,16)=feq + (1-omega)*fneq1*p2
                            
                            !18
                            feq=(rho(i+1,j,k-1)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(-3*(pyy(i+1,j,k-1) - 2*pzz(i+1,j,k-1) + 3*pyy(i+1,j,k-1)*u(i+1,j,k-1) - 6*pzz(i+1,j,k-1)*u(i+1,j,k-1) + 6*pxy(i+1,j,k-1)*v(i+1,j,k-1) - 6*pyz(i+1,j,k-1)*v(i+1,j,k-1) + 6*pxz(i+1,j,k-1)*(1 + 2*u(i+1,j,k-1) - 2*w(i+1,j,k-1)) - 3*pyy(i+1,j,k-1)*w(i+1,j,k-1) + pxx(i+1,j,k-1)*(-2 + 6*w(i+1,j,k-1))))/2.
                            f(i+1,j,k-1,18)=feq + (1-omega)*fneq1*p2
                            
                            !12
                            feq=(rho(i,j-1,k-1)*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*pxx(i,j-1,k-1)*(-1 + 3*v(i,j-1,k-1) + 3*w(i,j-1,k-1)))/2. + 3*(pyy(i,j-1,k-1) + pzz(i,j-1,k-1) + 3*pxy(i,j-1,k-1)*u(i,j-1,k-1) + 3*pxz(i,j-1,k-1)*u(i,j-1,k-1) - 3*pzz(i,j-1,k-1)*v(i,j-1,k-1) + pyz(i,j-1,k-1)*(3 - 6*v(i,j-1,k-1) - 6*w(i,j-1,k-1)) - 3*pyy(i,j-1,k-1)*w(i,j-1,k-1))
                            f(i,j-1,k-1,12)=feq + (1-omega)*p2*fneq1

                            !13
                            
                            feq=(rho(i,j+1,k-1)*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                            fneq1=(-3*(pxx(i,j+1,k-1) - 2*pyy(i,j+1,k-1) + 6*pyz(i,j+1,k-1) - 2*pzz(i,j+1,k-1) + 6*pxy(i,j+1,k-1)*u(i,j+1,k-1) - 6*pxz(i,j+1,k-1)*u(i,j+1,k-1) + 3*pxx(i,j+1,k-1)*v(i,j+1,k-1) + 12*pyz(i,j+1,k-1)*v(i,j+1,k-1) - 6*pzz(i,j+1,k-1)*v(i,j+1,k-1) - 3*pxx(i,j+1,k-1)*w(i,j+1,k-1) + 6*pyy(i,j+1,k-1)*w(i,j+1,k-1) - 12*pyz(i,j+1,k-1)*w(i,j+1,k-1)))/2.
                            f(i,j+1,k-1,13)=feq + (1-omega)*p2*fneq1
                            
                            !20
                            feq=(rho(i-1,j-1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*((pyy(i-1,j-1,k-1) + 3*pyz(i-1,j-1,k-1) + pzz(i-1,j-1,k-1))*(-1 + 3*u(i-1,j-1,k-1)) + pxy(i-1,j-1,k-1)*(-3 + 6*u(i-1,j-1,k-1)) + pxz(i-1,j-1,k-1)*(-3 + 6*u(i-1,j-1,k-1)) + 3*(2*pxy(i-1,j-1,k-1) + 3*pxz(i-1,j-1,k-1) + 2*pyz(i-1,j-1,k-1) + pzz(i-1,j-1,k-1))*v(i-1,j-1,k-1) + 3*(3*pxy(i-1,j-1,k-1) + &
                            2*pxz(i-1,j-1,k-1) + pyy(i-1,j-1,k-1) + 2*pyz(i-1,j-1,k-1))*w(i-1,j-1,k-1) + pxx(i-1,j-1,k-1)*(-1 + 3*v(i-1,j-1,k-1) + 3*w(i-1,j-1,k-1)))
                            f(i-1,j-1,k-1,20)=feq + (1-omega)*p3*fneq1
                            
                            !22
                            feq=(rho(i-1,j+1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i-1,j+1,k-1) + 3*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 3*pyz(i-1,j+1,k-1) + pzz(i-1,j+1,k-1) - 3*(2*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 3*pyz(i-1,j+1,k-1) + pzz(i-1,j+1,k-1))*u(i-1,j+1,k-1) + pxy(i-1,j+1,k-1)*(-3 + 6*u(i-1,j+1,k-1)) + 3*pxx(i-1,j+1,k-1)*v(i-1,j+1,k-1) &
                            - 6*pxy(i-1,j+1,k-1)*v(i-1,j+1,k-1) + 9*pxz(i-1,j+1,k-1)*v(i-1,j+1,k-1) - 6*pyz(i-1,j+1,k-1)*v(i-1,j+1,k-1) + 3*pzz(i-1,j+1,k-1)*v(i-1,j+1,k-1) - 3*(pxx(i-1,j+1,k-1) - 3*pxy(i-1,j+1,k-1) + 2*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 2*pyz(i-1,j+1,k-1))*w(i-1,j+1,k-1))

                            f(i-1,j+1,k-1,22)=feq + (1-omega)*p3*fneq1

                            !24
                            feq=(rho(i+1,j+1,k-1)*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                            fneq1=3*(pxx(i+1,j+1,k-1) - 3*pxz(i+1,j+1,k-1)*(1 + 2*u(i+1,j+1,k-1)) + (pyy(i+1,j+1,k-1) - 3*pyz(i+1,j+1,k-1) + pzz(i+1,j+1,k-1))*(1 + 3*u(i+1,j+1,k-1)) + pxy(i+1,j+1,k-1)*(3 + 6*u(i+1,j+1,k-1)) + 3*pxx(i+1,j+1,k-1)*v(i+1,j+1,k-1) + 6*pxy(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 9*pxz(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 6*pyz(i+1,j+1,k-1)*v(i+1,j+1,k-1) &
                            + 3*pzz(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 3*(pxx(i+1,j+1,k-1) + 3*pxy(i+1,j+1,k-1) - 2*pxz(i+1,j+1,k-1) + pyy(i+1,j+1,k-1) - 2*pyz(i+1,j+1,k-1))*w(i+1,j+1,k-1))
                            f(i+1,j+1,k-1,24)=feq + (1-omega)*p3*fneq1

                            !25
                            feq=(rho(i+1,j-1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*(-pxx(i+1,j-1,k-1) - (pyy(i+1,j-1,k-1) + 3*pyz(i+1,j-1,k-1) + pzz(i+1,j-1,k-1))*(1 + 3*u(i+1,j-1,k-1)) + pxy(i+1,j-1,k-1)*(3 + 6*u(i+1,j-1,k-1)) + pxz(i+1,j-1,k-1)*(3 + 6*u(i+1,j-1,k-1)) + 3*pxx(i+1,j-1,k-1)*v(i+1,j-1,k-1) &
                            - 6*pxy(i+1,j-1,k-1)*v(i+1,j-1,k-1) - 9*pxz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 6*pyz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 3*pzz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 3*(pxx(i+1,j-1,k-1) - 3*pxy(i+1,j-1,k-1) - 2*pxz(i+1,j-1,k-1) + pyy(i+1,j-1,k-1) + 2*pyz(i+1,j-1,k-1))*w(i+1,j-1,k-1))
                            f(i+1,j-1,k-1,25)=feq + (1-omega)*p3*fneq1
						endif
                    enddo
                enddo
            enddo    
        !*********************************** call other bcs:PERIODIC ************************     
            !periodic along x 
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
            !periodic along y
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
		
	endsubroutine
    !*****************************************************
	subroutine bcs_2c_turbulent_jet
	
	    implicit none
		
		if(mod(step,10).eq.0)then
			!$acc update host(w(nx/2-20:nx/2+20,ny/2-20:ny/2+20,1))
			do i=nx/2-20,nx/2+20
				do j=ny/2-20,ny/2+20
					if((float(i)-nx/2.0)**2 + (float(j)-ny/2.0)**2<=10**2)then
						call random_number(rrx)
						w(i,j,1)=uwall + 0.004*sqrt(-2.0*log(rrx))*cos(2*3.1415926535897932384626433832795028841971*rrx)
					endif
				enddo
			enddo
		   !$acc update device(w(nx/2-20:nx/2+20,ny/2-20:ny/2+20,1))
		endif
		!$acc kernels
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
						if(k.eq.1)then
							!5
							!w(i,j,k)=uwall + 0.02*sqrt(-2.0*log(ranx(ttt)))*cos(2*3.1415926535897932384626433832795028841971*rany(ttt))
                            feq=(rho(i,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k))))/27.
                            fneq1=(-3*(pxx(i,j,k+1) + pyy(i,j,k+1) - 2*pzz(i,j,k+1) + 6*pxz(i,j,k+1)*u(i,j,k+1) + 6*pyz(i,j,k+1)*v(i,j,k+1) + 3*pxx(i,j,k+1)*w(i,j,k+1) + 3*pyy(i,j,k+1)*w(i,j,k+1)))/2.
                            f(i,j,k+1,5)=feq + (1-omega)*p1*fneq1

                            !15
                            
                            feq=(rho(i+1,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i+1,j,k+1) + 2*pzz(i+1,j,k+1) - 3*pyy(i+1,j,k+1)*u(i+1,j,k+1) + 6*pzz(i+1,j,k+1)*u(i+1,j,k+1) - 6*pxy(i+1,j,k+1)*v(i+1,j,k+1) - 6*pyz(i+1,j,k+1)*v(i+1,j,k+1) - 3*pyy(i+1,j,k+1)*w(i+1,j,k+1) + 6*pxz(i+1,j,k+1)*(1 + 2*u(i+1,j,k+1) + 2*w(i+1,j,k+1)) + pxx(i+1,j,k+1)*(2 + 6*w(i+1,j,k+1))))/2.
                            f(i+1,j,k+1,15)=feq + (1-omega)*p2*fneq1

                            !17
                            
                            feq=(rho(i-1,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i-1,j,k+1) + 2*pzz(i-1,j,k+1) + 3*pyy(i-1,j,k+1)*u(i-1,j,k+1) - 6*pzz(i-1,j,k+1)*u(i-1,j,k+1) + 6*pxy(i-1,j,k+1)*v(i-1,j,k+1) - 6*pyz(i-1,j,k+1)*v(i-1,j,k+1) + 6*pxz(i-1,j,k+1)*(-1 + 2*u(i-1,j,k+1) - 2*w(i-1,j,k+1)) - 3*pyy(i-1,j,k+1)*w(i-1,j,k+1) + pxx(i-1,j,k+1)*(2 + 6*w(i-1,j,k+1))))/2.
                            f(i-1,j,k+1,17)=feq + (1-omega)*p2*fneq1

                            !11
                            
                            feq=(rho(i,j+1,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(-3*pxx(i,j+1,k+1)*(1 + 3*v(i,j+1,k+1) + 3*w(i,j+1,k+1)))/2. + 3*(pyy(i,j+1,k+1) + pzz(i,j+1,k+1) - 3*pxy(i,j+1,k+1)*u(i,j+1,k+1) - 3*pxz(i,j+1,k+1)*u(i,j+1,k+1) + 3*pzz(i,j+1,k+1)*v(i,j+1,k+1) + 3*pyy(i,j+1,k+1)*w(i,j+1,k+1) + pyz(i,j+1,k+1)*(3 + 6*v(i,j+1,k+1) + 6*w(i,j+1,k+1)))
                            f(i,j+1,k+1,11)=feq +(1-omega)*p2*fneq1
                            
                            !14
                            feq=(rho(i,j-1,k+1)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*pxx(i,j-1,k+1)*(-1 + 3*v(i,j-1,k+1) - 3*w(i,j-1,k+1)))/2. + 3*(pyy(i,j-1,k+1) + pzz(i,j-1,k+1) + 3*pxy(i,j-1,k+1)*u(i,j-1,k+1) - 3*pxz(i,j-1,k+1)*u(i,j-1,k+1) - 3*pzz(i,j-1,k+1)*v(i,j-1,k+1) + pyz(i,j-1,k+1)*(-3 + 6*v(i,j-1,k+1) - 6*w(i,j-1,k+1)) + 3*pyy(i,j-1,k+1)*w(i,j-1,k+1))
                            f(i,j-1,k+1,14)=feq + (1-omega)*p2*fneq1

                            !19
                            feq=(rho(i+1,j+1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i+1,j+1,k+1) + (pyy(i+1,j+1,k+1) + 3*pyz(i+1,j+1,k+1) + pzz(i+1,j+1,k+1))*(1 + 3*u(i+1,j+1,k+1)) + pxy(i+1,j+1,k+1)*(3 + 6*u(i+1,j+1,k+1)) + pxz(i+1,j+1,k+1)*(3 + 6*u(i+1,j+1,k+1)) + 3*pxx(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 6*pxy(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 9*pxz(i+1,j+1,k+1)*v(i+1,j+1,k+1) &
                            + 6*pyz(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 3*pzz(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 3*(pxx(i+1,j+1,k+1) + 3*pxy(i+1,j+1,k+1) + 2*pxz(i+1,j+1,k+1) + pyy(i+1,j+1,k+1) + 2*pyz(i+1,j+1,k+1))*w(i+1,j+1,k+1))
                            f(i+1,j+1,k+1,19)=feq + (1-omega)*fneq1*p3

                            !21
                            feq=(rho(i+1,j-1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i+1,j-1,k+1) - 3*pxy(i+1,j-1,k+1)*(1 + 2*u(i+1,j-1,k+1)) + (pyy(i+1,j-1,k+1) - 3*pyz(i+1,j-1,k+1) + pzz(i+1,j-1,k+1))*(1 + 3*u(i+1,j-1,k+1)) + pxz(i+1,j-1,k+1)*(3 + 6*u(i+1,j-1,k+1)) - 3*pxx(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 6*pxy(i+1,j-1,k+1)*v(i+1,j-1,k+1) &
                            - 9*pxz(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 6*pyz(i+1,j-1,k+1)*v(i+1,j-1,k+1) - 3*pzz(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 3*(pxx(i+1,j-1,k+1) - 3*pxy(i+1,j-1,k+1) + 2*pxz(i+1,j-1,k+1) + pyy(i+1,j-1,k+1) - 2*pyz(i+1,j-1,k+1))*w(i+1,j-1,k+1))
                            f(i+1,j-1,k+1,21)=feq + (1-omega)*fneq1*p3

                            !23
                            feq=(rho(i-1,j-1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i-1,j-1,k+1) + 3*pxy(i-1,j-1,k+1) - 3*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 3*pyz(i-1,j-1,k+1) + pzz(i-1,j-1,k+1) - 3*(2*pxy(i-1,j-1,k+1) - 2*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 3*pyz(i-1,j-1,k+1) + pzz(i-1,j-1,k+1))*u(i-1,j-1,k+1) &
                            - 3*pxx(i-1,j-1,k+1)*v(i-1,j-1,k+1) - 6*pxy(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 9*pxz(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 6*pyz(i-1,j-1,k+1)*v(i-1,j-1,k+1) - 3*pzz(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 3*(pxx(i-1,j-1,k+1) + 3*pxy(i-1,j-1,k+1) - 2*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 2*pyz(i-1,j-1,k+1))*w(i-1,j-1,k+1))
                            f(i-1,j-1,k+1,23)=feq + (1-omega)*fneq1*p3

                            !26
                            feq=(rho(i-1,j+1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i-1,j+1,k+1) - (pyy(i-1,j+1,k+1) + 3*pyz(i-1,j+1,k+1) + pzz(i-1,j+1,k+1))*(-1 + 3*u(i-1,j+1,k+1)) + pxy(i-1,j+1,k+1)*(-3 + 6*u(i-1,j+1,k+1)) + pxz(i-1,j+1,k+1)*(-3 + 6*u(i-1,j+1,k+1)) + 3*pxx(i-1,j+1,k+1)*v(i-1,j+1,k+1) - 6*pxy(i-1,j+1,k+1)*v(i-1,j+1,k+1) &
                            - 9*pxz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 6*pyz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 3*pzz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 3*(pxx(i-1,j+1,k+1) - 3*pxy(i-1,j+1,k+1) - 2*pxz(i-1,j+1,k+1) + pyy(i-1,j+1,k+1) + 2*pyz(i-1,j+1,k+1))*w(i-1,j+1,k+1))
                            f(i-1,j+1,k+1,26)=feq + (1-omega)*fneq1*p3

                            !******************************************!
                            
                            feq=p1g*phi(i,j,k)*(1 + 3*w(i,j,k))
                            g(i,j,k+1,5)=feq
                            
						elseif(k.eq.nz)then
							!6
							w(i,j,k)=w(i,j,nz-1)
                            feq=(rho(i,j,k-1)*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                            fneq1=(3*(-pyy(i,j,k-1) + 2*pzz(i,j,k-1) + 6*pxz(i,j,k-1)*u(i,j,k-1) + 6*pyz(i,j,k-1)*v(i,j,k-1) + 3*pyy(i,j,k-1)*w(i,j,k-1) + pxx(i,j,k-1)*(-1 + 3*w(i,j,k-1))))/2.
                            f(i,j,k-1,6)=feq + (1-omega)*fneq1*p1
                            
                            !16
                            feq=(rho(i-1,j,k-1)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*(-pyy(i-1,j,k-1) + 2*pzz(i-1,j,k-1) + 3*pyy(i-1,j,k-1)*u(i-1,j,k-1) - 6*pzz(i-1,j,k-1)*u(i-1,j,k-1) + 6*pxy(i-1,j,k-1)*v(i-1,j,k-1) + 6*pyz(i-1,j,k-1)*v(i-1,j,k-1) + pxx(i-1,j,k-1)*(2 - 6*w(i-1,j,k-1)) + 3*pyy(i-1,j,k-1)*w(i-1,j,k-1) - 6*pxz(i-1,j,k-1)*(-1 + 2*u(i-1,j,k-1) + 2*w(i-1,j,k-1))))/2.
                            f(i-1,j,k-1,16)=feq + (1-omega)*fneq1*p2
                            
                            !18
                            feq=(rho(i+1,j,k-1)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(-3*(pyy(i+1,j,k-1) - 2*pzz(i+1,j,k-1) + 3*pyy(i+1,j,k-1)*u(i+1,j,k-1) - 6*pzz(i+1,j,k-1)*u(i+1,j,k-1) + 6*pxy(i+1,j,k-1)*v(i+1,j,k-1) - 6*pyz(i+1,j,k-1)*v(i+1,j,k-1) + 6*pxz(i+1,j,k-1)*(1 + 2*u(i+1,j,k-1) - 2*w(i+1,j,k-1)) - 3*pyy(i+1,j,k-1)*w(i+1,j,k-1) + pxx(i+1,j,k-1)*(-2 + 6*w(i+1,j,k-1))))/2.
                            f(i+1,j,k-1,18)=feq + (1-omega)*fneq1*p2
                            
                            !12
                            feq=(rho(i,j-1,k-1)*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*pxx(i,j-1,k-1)*(-1 + 3*v(i,j-1,k-1) + 3*w(i,j-1,k-1)))/2. + 3*(pyy(i,j-1,k-1) + pzz(i,j-1,k-1) + 3*pxy(i,j-1,k-1)*u(i,j-1,k-1) + 3*pxz(i,j-1,k-1)*u(i,j-1,k-1) - 3*pzz(i,j-1,k-1)*v(i,j-1,k-1) + pyz(i,j-1,k-1)*(3 - 6*v(i,j-1,k-1) - 6*w(i,j-1,k-1)) - 3*pyy(i,j-1,k-1)*w(i,j-1,k-1))
                            f(i,j-1,k-1,12)=feq + (1-omega)*p2*fneq1

                            !13
                            
                            feq=(rho(i,j+1,k-1)*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                            fneq1=(-3*(pxx(i,j+1,k-1) - 2*pyy(i,j+1,k-1) + 6*pyz(i,j+1,k-1) - 2*pzz(i,j+1,k-1) + 6*pxy(i,j+1,k-1)*u(i,j+1,k-1) - 6*pxz(i,j+1,k-1)*u(i,j+1,k-1) + 3*pxx(i,j+1,k-1)*v(i,j+1,k-1) + 12*pyz(i,j+1,k-1)*v(i,j+1,k-1) - 6*pzz(i,j+1,k-1)*v(i,j+1,k-1) - 3*pxx(i,j+1,k-1)*w(i,j+1,k-1) + 6*pyy(i,j+1,k-1)*w(i,j+1,k-1) - 12*pyz(i,j+1,k-1)*w(i,j+1,k-1)))/2.
                            f(i,j+1,k-1,13)=feq + (1-omega)*p2*fneq1
                            
                            !20
                            feq=(rho(i-1,j-1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*((pyy(i-1,j-1,k-1) + 3*pyz(i-1,j-1,k-1) + pzz(i-1,j-1,k-1))*(-1 + 3*u(i-1,j-1,k-1)) + pxy(i-1,j-1,k-1)*(-3 + 6*u(i-1,j-1,k-1)) + pxz(i-1,j-1,k-1)*(-3 + 6*u(i-1,j-1,k-1)) + 3*(2*pxy(i-1,j-1,k-1) + 3*pxz(i-1,j-1,k-1) + 2*pyz(i-1,j-1,k-1) + pzz(i-1,j-1,k-1))*v(i-1,j-1,k-1) + 3*(3*pxy(i-1,j-1,k-1) + &
                            2*pxz(i-1,j-1,k-1) + pyy(i-1,j-1,k-1) + 2*pyz(i-1,j-1,k-1))*w(i-1,j-1,k-1) + pxx(i-1,j-1,k-1)*(-1 + 3*v(i-1,j-1,k-1) + 3*w(i-1,j-1,k-1)))
                            f(i-1,j-1,k-1,20)=feq + (1-omega)*p3*fneq1
                            
                            !22
                            feq=(rho(i-1,j+1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i-1,j+1,k-1) + 3*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 3*pyz(i-1,j+1,k-1) + pzz(i-1,j+1,k-1) - 3*(2*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 3*pyz(i-1,j+1,k-1) + pzz(i-1,j+1,k-1))*u(i-1,j+1,k-1) + pxy(i-1,j+1,k-1)*(-3 + 6*u(i-1,j+1,k-1)) + 3*pxx(i-1,j+1,k-1)*v(i-1,j+1,k-1) &
                            - 6*pxy(i-1,j+1,k-1)*v(i-1,j+1,k-1) + 9*pxz(i-1,j+1,k-1)*v(i-1,j+1,k-1) - 6*pyz(i-1,j+1,k-1)*v(i-1,j+1,k-1) + 3*pzz(i-1,j+1,k-1)*v(i-1,j+1,k-1) - 3*(pxx(i-1,j+1,k-1) - 3*pxy(i-1,j+1,k-1) + 2*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 2*pyz(i-1,j+1,k-1))*w(i-1,j+1,k-1))

                            f(i-1,j+1,k-1,22)=feq + (1-omega)*p3*fneq1

                            !24
                            feq=(rho(i+1,j+1,k-1)*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                            fneq1=3*(pxx(i+1,j+1,k-1) - 3*pxz(i+1,j+1,k-1)*(1 + 2*u(i+1,j+1,k-1)) + (pyy(i+1,j+1,k-1) - 3*pyz(i+1,j+1,k-1) + pzz(i+1,j+1,k-1))*(1 + 3*u(i+1,j+1,k-1)) + pxy(i+1,j+1,k-1)*(3 + 6*u(i+1,j+1,k-1)) + 3*pxx(i+1,j+1,k-1)*v(i+1,j+1,k-1) + 6*pxy(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 9*pxz(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 6*pyz(i+1,j+1,k-1)*v(i+1,j+1,k-1) &
                            + 3*pzz(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 3*(pxx(i+1,j+1,k-1) + 3*pxy(i+1,j+1,k-1) - 2*pxz(i+1,j+1,k-1) + pyy(i+1,j+1,k-1) - 2*pyz(i+1,j+1,k-1))*w(i+1,j+1,k-1))
                            f(i+1,j+1,k-1,24)=feq + (1-omega)*p3*fneq1

                            !25
                            feq=(rho(i+1,j-1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*(-pxx(i+1,j-1,k-1) - (pyy(i+1,j-1,k-1) + 3*pyz(i+1,j-1,k-1) + pzz(i+1,j-1,k-1))*(1 + 3*u(i+1,j-1,k-1)) + pxy(i+1,j-1,k-1)*(3 + 6*u(i+1,j-1,k-1)) + pxz(i+1,j-1,k-1)*(3 + 6*u(i+1,j-1,k-1)) + 3*pxx(i+1,j-1,k-1)*v(i+1,j-1,k-1) &
                            - 6*pxy(i+1,j-1,k-1)*v(i+1,j-1,k-1) - 9*pxz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 6*pyz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 3*pzz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 3*(pxx(i+1,j-1,k-1) - 3*pxy(i+1,j-1,k-1) - 2*pxz(i+1,j-1,k-1) + pyy(i+1,j-1,k-1) + 2*pyz(i+1,j-1,k-1))*w(i+1,j-1,k-1))
                            f(i+1,j-1,k-1,25)=feq + (1-omega)*p3*fneq1

                            !******************************************!
                            
                            
                            !g6
                            feq=p1g*phi(i,j,k-1)*(1 - 3*0.05) !3*w(i,j,k))
                            g(i,j,k-1,6)=feq

						endif
                        
                    enddo
                enddo
            enddo    
        !*********************************** call other bcs:PERIODIC ************************     
            !periodic along x 
                ! 0  1   2  3   4   5   6   7    8   9   10  11   12  13   14  15   16   17   18  19  20  21  22  23  24  25  26

            !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1,  1, -1,  1, -1, -1,  1,  1, -1/)
            !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0,  1  -1, -1,  1, -1,  1, -1,  1/)
            !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1,  1, -1,  1, -1,  1, -1, -1,  1/)
            !f
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

            ! periodic g and phi x-y
            g(2,2:ny-1,2:nz-1,1)=g(nx-1,2:ny-1,2:nz-1,1)
            g(nx-1,2:ny-1,2:nz-1,2)=g(2,2:ny-1,2:nz-1,2)
            g(2:nx-1,2,2:nz-1,3)=g(2:nx-1,ny-1,2:nz-1,3)
            g(2:nx-1,ny-1,2:nz-1,4)=g(2:nx-1,2,2:nz-1,4)
            phi(1,2:ny-1,2:nz-1)=phi(nx-1,2:ny-1,2:nz-1)
            phi(nx,2:ny-1,2:nz-1)=phi(2,2:ny-1,2:nz-1)
            phi(2:nx-1,2:ny-1,nz)=phi(2:nx-1,2:ny-1,nz-1)

            !periodic along y
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
		
	endsubroutine
    !*****************************************************
    subroutine bcs_TSLB_turbojet
        
        implicit none 

        if(mod(step,10).eq.0)then
        !$acc update host(w(nx/2-20:nx/2+20,ny/2-20:ny/2+20,1))
        do i=nx/2-20,nx/2+20
          do j=ny/2-20,ny/2+20
            if((float(i)-nx/2.0)**2 + (float(j)-ny/2.0)**2<=10**2)then
              call random_number(rrx)
              w(i,j,1)=uwall + 0.004*sqrt(-2.0*log(rrx))*cos(2*3.1415926535897932384626433832795028841971*rrx)
            endif
          enddo
        enddo
        !$acc update device(w(nx/2-20:nx/2+20,ny/2-20:ny/2+20,1))
		    endif
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !$acc kernels
        !$acc loop independent 
        do k=1,nz
            !$acc loop independent 
            do j=1,ny
                !east and west
                rho(1,j,k)=rho(2,j,k)
                u(1,j,k)=u(2,j,k)
                v(1,j,k)=v(2,j,k)
                w(1,j,k)=w(2,j,k)
                pxx(1,j,k)=pxx(2,j,k)
                pxy(1,j,k)=pxy(2,j,k)
                pxz(1,j,k)=pxz(2,j,k)
                pyy(1,j,k)=pyy(2,j,k)
                pyz(1,j,k)=pyz(2,j,k)
                pzz(1,j,k)=pzz(2,j,k)
                !
                rho(nx,j,k)=rho(nx-1,j,k)
                u(nx,j,k)=u(nx-1,j,k)
                v(nx,j,k)=v(nx-1,j,k)
                v(nx,j,k)=v(nx-1,j,k)
                pxx(nx,j,k)=pxx(nx-1,j,k)
                pxy(nx,j,k)=pxy(nx-1,j,k)
                pxz(nx,j,k)=pxz(nx-1,j,k)
                pyy(nx,j,k)=pyy(nx-1,j,k)
                pyz(nx,j,k)=pyz(nx-1,j,k)
                pzz(nx,j,k)=pzz(nx-1,j,k)
            enddo
        enddo
        !$acc end kernels
        !$acc kernels
        !$acc loop independent 
        do k=1,nz
            !$acc loop independent 
            do i=1,nx
                !front and rear
                rho(i,1,k)=rho(i,2,k)
                u(i,1,k)=u(i,2,k)
                v(i,1,k)=v(i,2,k)
                w(i,1,k)=w(i,2,k)
                pxx(i,1,k)=pxx(i,2,k)
                pxy(i,1,k)=pxy(i,2,k)
                pxz(i,1,k)=pxz(i,2,k)
                pyy(i,1,k)=pyy(i,2,k)
                pyz(i,1,k)=pyz(i,2,k)
                pzz(i,1,k)=pzz(i,2,k)
                !
                rho(i,ny,k)=rho(i,ny-1,k)
                u(i,ny,k)=u(i,ny-1,k)
                v(i,ny,k)=v(i,ny-1,k)
                w(i,ny,k)=w(i,ny-1,k)
                pxx(i,ny,k)=pxx(i,ny-1,k)
                pxy(i,ny,k)=pxy(i,ny-1,k)
                pxz(i,ny,k)=pxz(i,ny-1,k)
                pyy(i,ny,k)=pyy(i,ny-1,k)
                pyz(i,ny,k)=pyz(i,ny-1,k)
                pzz(i,ny,k)=pzz(i,ny-1,k)
            enddo
        enddo
        !$acc end kernels
        
        !$acc kernels
        !$acc loop independent 
        do j=1,ny
            !$acc loop independent 
            do i=1,nx
                !north south
                rho(i,j,1)=rho(i,j,2)
                u(i,j,1)=0.0 !u(i,j,2)
                v(i,j,1)=0.0 !v(i,j,2)
                !w(i,j,1)=w(i,j,1)
                pxx(i,j,1)=pxx(i,j,2)
                pxy(i,j,1)=pxy(i,j,2)
                pxz(i,j,1)=pxz(i,j,2)
                pyy(i,j,1)=pyy(i,j,2)
                pyz(i,j,1)=pyz(i,j,2)
                pzz(i,j,1)=pzz(i,j,2)
                !
                rho(i,j,nz)=rho(i,j,nz-1)
                u(i,j,nz)=u(i,j,nz-1)
                v(i,j,nz)=v(i,j,nz-1)
                w(i,j,nz)=w(i,j,nz-1)
                pxx(i,j,nz)=pxx(i,j,nz-1)
                pxy(i,j,nz)=pxy(i,j,nz-1)
                pxz(i,j,nz)=pxz(i,j,nz-1)
                pyy(i,j,nz)=pyy(i,j,nz-1)
                pyz(i,j,nz)=pyz(i,j,nz-1)
                pzz(i,j,nz)=pzz(i,j,nz-1)
            enddo
        enddo
        !$acc end kernels
    endsubroutine
    !*****************************************************
    !*****************************************************
    subroutine bcs_TSLB_only_z_turbojet
        
        implicit none 
        
        integer :: subchords(3)
        
      !devo trovare quali processi hanno in carico i nodi lungo il piano gk=1
      gk=1
      !subchords(1)=(gi-1)/nx
      !subchords(2)=(gj-1)/ny
      subchords(3)=(gk-1)/nz
      !sono buoni tutti i processi che hanno subchords(3)==coords(3)
      if(subchords(3)==coords(3))then
        if(mod(step,10).eq.0)then
            !devo fare update di step che mi serve come seed 
            !$acc update device(step)
            !$acc kernels present(w,step,myoffset,coords,nx,ny) 
            !$acc loop independent collapse(2)  private(i,j,k,gi,gj,gk,rrx)
            do j=1,ny
              do i=1,nx
                gi=nx*coords(1)+i
                gj=ny*coords(2)+j
                if((float(gi)-lx/2.0)**2 + (float(gj)-ly/2.0)**2<=10**2)then
                              !call random_number(rrx)
                              !sto sul piano gk=1
                              gk=1
                              !myoffset(3) è il mio offset lungo z del mio sottodominio MPI e mi ridà il valore di k nel sottodominio
                              k=gk-myoffset(3)
                              !è uno pseudo generatore che da un numero randomico partendo da 4 integer come seed
                              !devi fare in modo che ogni lattice point ad ogni time step abbia seed diversi
                              !quindi uso come seed la posizione i j k e il timestep come quarto seed 
                              !il fatto che tutti i seed siano diversi è perchè può essere chiamata da più threads contemporaneamente
                              !invece se tu hai un unico seed lo devi mettere in save per tutti i threads e poi dipende da chi chiama prima (dipende dall'ordine di chiamata)
                    rrx=rand_noseeded(gi,gj,gk,step)		
                    w(i,j,k)=uwall + 0.004*sqrt(-2.0*log(rrx))*cos(2*3.1415926535897932384626433832795028841971*rrx)
                endif
            enddo
			    enddo
          !$acc end kernels
          
		    endif
      
      
        !$acc kernels present(rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz,myoffset) async
        !$acc loop independent collapse(2) private(i,j,k,gk)
        do j=1,ny
            do i=1,nx
                !sto sul piano gk=1
                 gk=1
                 !myoffset(3) è il mio offset lungo z del mio sottodominio MPI e mi ridà il valore di k nel sottodominio
                 k=gk-myoffset(3)
                !south
                !davanti a me è k+1
                rho(i,j,k)=rho(i,j,k+1)
                u(i,j,k)=0.0 !u(i,j,k+1)
                v(i,j,k)=0.0 !v(i,j,k+1)
                !w(i,j,k)=w(i,j,k)
                pxx(i,j,k)=pxx(i,j,k+1)
                pxy(i,j,k)=pxy(i,j,k+1)
                pxz(i,j,k)=pxz(i,j,k+1)
                pyy(i,j,k)=pyy(i,j,k+1)
                pyz(i,j,k)=pyz(i,j,k+1)
                pzz(i,j,k)=pzz(i,j,k+1)
            enddo
        enddo
        !$acc end kernels
      
      endif
      
      !devo trovare quali processi hanno in carico i nodi lungo il piano gk=1
      gk=lz
      !subchords(1)=(oi-1)/nx
      !subchords(2)=(oj-1)/ny
      subchords(3)=(gk-1)/nz
      !sono buoni tutti i processi che hanno subchords(3)==coords(3)
      if(subchords(3)==coords(3))then     
        !$acc kernels present(rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz,myoffset,lz) async
        !$acc loop independent collapse(2) private(i,j,k)
        do j=1,ny
            do i=1,nx
                !north
                !sto sul piano gk=lz
                !gk=lz
                !myoffset(3) è il mio offset lungo z del mio sottodominio MPI e mi ridà il valore di k nel sottodominio
                k=lz-myoffset(3)!gk-myoffset(3)
                !dietro a me è k-1
                rho(i,j,k)=rho(i,j,k-1)
                u(i,j,k)=u(i,j,k-1)
                v(i,j,k)=v(i,j,k-1)
                w(i,j,k)=w(i,j,k-1)
                pxx(i,j,k)=pxx(i,j,k-1)
                pxy(i,j,k)=pxy(i,j,k-1)
                pxz(i,j,k)=pxz(i,j,k-1)
                pyy(i,j,k)=pyy(i,j,k-1)
                pyz(i,j,k)=pyz(i,j,k-1)
                pzz(i,j,k)=pzz(i,j,k-1)
            enddo
        enddo
        !$acc end kernels
      endif
      
      !$acc wait
      
    endsubroutine
    
    !****************************************************
    subroutine pbcs
        
        implicit none
        integer :: oi,oj,ok

        !periodic along x
                ! 0  1   2  3   4   5   6   7    8   9   10  11   12  13   14  15   16   17   18  19  20  21  22  23  24  25  26
            !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1,  1, -1,  1, -1, -1,  1,  1, -1/)
            !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0,  1, -1, -1,  1, -1,  1, -1,  1/)
            !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1,  1, -1,  1, -1,  1, -1, -1,  1/)
       
       
        !dir l   1 disp    1   0   0
        !$acc kernels
        !$acc loop independent
        do k=1,nz
          !$acc loop independent
          do j=1,ny
            !$acc loop independent
            do i=nx+1,nx+1
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,1)=f(i,j,k,1)
              f(oi,oj,ok,7)=f(i,j,k,7)
              f(oi,oj,ok,9)=f(i,j,k,9)
              f(oi,oj,ok,15)=f(i,j,k,15)
              f(oi,oj,ok,18)=f(i,j,k,18)
              f(oi,oj,ok,19)=f(i,j,k,19)
              f(oi,oj,ok,21)=f(i,j,k,21)
              f(oi,oj,ok,24)=f(i,j,k,24)
              f(oi,oj,ok,25)=f(i,j,k,25)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l   2 disp   -1   0   0
        !$acc kernels
        !$acc loop independent
        do k=1,nz
          !$acc loop independent
          do j=1,ny
            !$acc loop independent
            do i=0,0
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,2)=f(i,j,k,2)
              f(oi,oj,ok,8)=f(i,j,k,8)
              f(oi,oj,ok,10)=f(i,j,k,10)
              f(oi,oj,ok,16)=f(i,j,k,16)
              f(oi,oj,ok,17)=f(i,j,k,17)
              f(oi,oj,ok,20)=f(i,j,k,20)
              f(oi,oj,ok,22)=f(i,j,k,22)
              f(oi,oj,ok,23)=f(i,j,k,23)
              f(oi,oj,ok,26)=f(i,j,k,26)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l   3 disp    0   1   0
        !$acc kernels
        !$acc loop independent
        do k=1,nz
          !$acc loop independent
          do j=ny+1,ny+1
            !$acc loop independent
            do i=1,nx
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,3)=f(i,j,k,3)
              f(oi,oj,ok,7)=f(i,j,k,7)
              f(oi,oj,ok,10)=f(i,j,k,10)
              f(oi,oj,ok,11)=f(i,j,k,11)
              f(oi,oj,ok,13)=f(i,j,k,13)
              f(oi,oj,ok,19)=f(i,j,k,19)
              f(oi,oj,ok,22)=f(i,j,k,22)
              f(oi,oj,ok,24)=f(i,j,k,24)
              f(oi,oj,ok,26)=f(i,j,k,26)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l   4 disp    0  -1   0
        !$acc kernels
        !$acc loop independent
        do k=1,nz
          !$acc loop independent
          do j=0,0
            !$acc loop independent
            do i=1,nx
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,4)=f(i,j,k,4)
              f(oi,oj,ok,8)=f(i,j,k,8)
              f(oi,oj,ok,9)=f(i,j,k,9)
              f(oi,oj,ok,12)=f(i,j,k,12)
              f(oi,oj,ok,14)=f(i,j,k,14)
              f(oi,oj,ok,20)=f(i,j,k,20)
              f(oi,oj,ok,21)=f(i,j,k,21)
              f(oi,oj,ok,23)=f(i,j,k,23)
              f(oi,oj,ok,25)=f(i,j,k,25)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l   5 disp    0   0   1
        !$acc kernels
        !$acc loop independent
        do k=nz+1,nz+1
          !$acc loop independent
          do j=1,ny
            !$acc loop independent
            do i=1,nx
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,5)=f(i,j,k,5)
              f(oi,oj,ok,11)=f(i,j,k,11)
              f(oi,oj,ok,14)=f(i,j,k,14)
              f(oi,oj,ok,15)=f(i,j,k,15)
              f(oi,oj,ok,17)=f(i,j,k,17)
              f(oi,oj,ok,19)=f(i,j,k,19)
              f(oi,oj,ok,21)=f(i,j,k,21)
              f(oi,oj,ok,23)=f(i,j,k,23)
              f(oi,oj,ok,26)=f(i,j,k,26)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l   6 disp    0   0  -1
        !$acc kernels
        !$acc loop independent
        do k=0,0
          !$acc loop independent
          do j=1,ny
            !$acc loop independent
            do i=1,nx
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,6)=f(i,j,k,6)
              f(oi,oj,ok,12)=f(i,j,k,12)
              f(oi,oj,ok,13)=f(i,j,k,13)
              f(oi,oj,ok,16)=f(i,j,k,16)
              f(oi,oj,ok,18)=f(i,j,k,18)
              f(oi,oj,ok,20)=f(i,j,k,20)
              f(oi,oj,ok,22)=f(i,j,k,22)
              f(oi,oj,ok,24)=f(i,j,k,24)
              f(oi,oj,ok,25)=f(i,j,k,25)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l   7 disp    1   1   0
        !$acc kernels
        !$acc loop independent
        do k=1,nz
          !$acc loop independent
          do j=ny+1,ny+1
            !$acc loop independent
            do i=nx+1,nx+1
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,7)=f(i,j,k,7)
              f(oi,oj,ok,19)=f(i,j,k,19)
              f(oi,oj,ok,24)=f(i,j,k,24)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l   8 disp   -1  -1   0
        !$acc kernels
        !$acc loop independent
        do k=1,nz
          !$acc loop independent
          do j=0,0
            !$acc loop independent
            do i=0,0
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,8)=f(i,j,k,8)
              f(oi,oj,ok,20)=f(i,j,k,20)
              f(oi,oj,ok,23)=f(i,j,k,23)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l   9 disp    1  -1   0
        !$acc kernels
        !$acc loop independent
        do k=1,nz
          !$acc loop independent
          do j=0,0
            !$acc loop independent
            do i=nx+1,nx+1
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,9)=f(i,j,k,9)
              f(oi,oj,ok,21)=f(i,j,k,21)
              f(oi,oj,ok,25)=f(i,j,k,25)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  10 disp   -1   1   0
        !$acc kernels
        !$acc loop independent
        do k=1,nz
          !$acc loop independent
          do j=ny+1,ny+1
            !$acc loop independent
            do i=0,0
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,10)=f(i,j,k,10)
              f(oi,oj,ok,22)=f(i,j,k,22)
              f(oi,oj,ok,26)=f(i,j,k,26)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  11 disp    0   1   1
        !$acc kernels
        !$acc loop independent
        do k=nz+1,nz+1
          !$acc loop independent
          do j=ny+1,ny+1
            !$acc loop independent
            do i=1,nx
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,11)=f(i,j,k,11)
              f(oi,oj,ok,19)=f(i,j,k,19)
              f(oi,oj,ok,26)=f(i,j,k,26)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  12 disp    0  -1  -1
        !$acc kernels
        !$acc loop independent
        do k=0,0
          !$acc loop independent
          do j=0,0
            !$acc loop independent
            do i=1,nx
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,12)=f(i,j,k,12)
              f(oi,oj,ok,20)=f(i,j,k,20)
              f(oi,oj,ok,25)=f(i,j,k,25)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  13 disp    0   1  -1
        !$acc kernels
        !$acc loop independent
        do k=0,0
          !$acc loop independent
          do j=ny+1,ny+1
            !$acc loop independent
            do i=1,nx
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,13)=f(i,j,k,13)
              f(oi,oj,ok,22)=f(i,j,k,22)
              f(oi,oj,ok,24)=f(i,j,k,24)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  14 disp    0  -1   1
        !$acc kernels
        !$acc loop independent
        do k=nz+1,nz+1
          !$acc loop independent
          do j=0,0
            !$acc loop independent
            do i=1,nx
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,14)=f(i,j,k,14)
              f(oi,oj,ok,21)=f(i,j,k,21)
              f(oi,oj,ok,23)=f(i,j,k,23)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  15 disp    1   0   1
        !$acc kernels
        !$acc loop independent
        do k=nz+1,nz+1
          !$acc loop independent
          do j=1,ny
            !$acc loop independent
            do i=nx+1,nx+1
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,15)=f(i,j,k,15)
              f(oi,oj,ok,19)=f(i,j,k,19)
              f(oi,oj,ok,21)=f(i,j,k,21)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  16 disp   -1   0  -1
        !$acc kernels
        !$acc loop independent
        do k=0,0
          !$acc loop independent
          do j=1,ny
            !$acc loop independent
            do i=0,0
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,16)=f(i,j,k,16)
              f(oi,oj,ok,20)=f(i,j,k,20)
              f(oi,oj,ok,22)=f(i,j,k,22)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  17 disp   -1   0   1
        !$acc kernels
        !$acc loop independent
        do k=nz+1,nz+1
          !$acc loop independent
          do j=1,ny
            !$acc loop independent
            do i=0,0
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,17)=f(i,j,k,17)
              f(oi,oj,ok,23)=f(i,j,k,23)
              f(oi,oj,ok,26)=f(i,j,k,26)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  18 disp    1   0  -1
        !$acc kernels
        !$acc loop independent
        do k=0,0
          !$acc loop independent
          do j=1,ny
            !$acc loop independent
            do i=nx+1,nx+1
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,18)=f(i,j,k,18)
              f(oi,oj,ok,24)=f(i,j,k,24)
              f(oi,oj,ok,25)=f(i,j,k,25)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  19 disp    1   1   1
        !$acc kernels
        !$acc loop independent
        do k=nz+1,nz+1
          !$acc loop independent
          do j=ny+1,ny+1
            !$acc loop independent
            do i=nx+1,nx+1
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,19)=f(i,j,k,19)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  20 disp   -1  -1  -1
        !$acc kernels
        !$acc loop independent
        do k=0,0
          !$acc loop independent
          do j=0,0
            !$acc loop independent
            do i=0,0
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,20)=f(i,j,k,20)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  21 disp    1  -1   1
        !$acc kernels
        !$acc loop independent
        do k=nz+1,nz+1
          !$acc loop independent
          do j=0,0
            !$acc loop independent
            do i=nx+1,nx+1
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,21)=f(i,j,k,21)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  22 disp   -1   1  -1
        !$acc kernels
        !$acc loop independent
        do k=0,0
          !$acc loop independent
          do j=ny+1,ny+1
            !$acc loop independent
            do i=0,0
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,22)=f(i,j,k,22)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  23 disp   -1  -1   1
        !$acc kernels
        !$acc loop independent
        do k=nz+1,nz+1
          !$acc loop independent
          do j=0,0
            !$acc loop independent
            do i=0,0
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,23)=f(i,j,k,23)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  24 disp    1   1  -1
        !$acc kernels
        !$acc loop independent
        do k=0,0
          !$acc loop independent
          do j=ny+1,ny+1
            !$acc loop independent
            do i=nx+1,nx+1
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,24)=f(i,j,k,24)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  25 disp    1  -1  -1
        !$acc kernels
        !$acc loop independent
        do k=0,0
          !$acc loop independent
          do j=0,0
            !$acc loop independent
            do i=nx+1,nx+1
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,25)=f(i,j,k,25)
            enddo
          enddo
        enddo
        !$acc end kernels
        !dir l  26 disp   -1   1   1
        !$acc kernels
        !$acc loop independent
        do k=nz+1,nz+1
          !$acc loop independent
          do j=ny+1,ny+1
            !$acc loop independent
            do i=0,0
              oi=i
              oj=j
              ok=k
              if(periodic(1))oi=mod(oi+nx-1,nx)+1
              if(periodic(2))oj=mod(oj+ny-1,ny)+1
              if(periodic(3))ok=mod(ok+nz-1,nz)+1
              f(oi,oj,ok,26)=f(i,j,k,26)
            enddo
          enddo
        enddo
        !$acc end kernels
		

		
        

    endsubroutine

endmodule
