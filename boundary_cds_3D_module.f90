module bcs3D
    
    use vars
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
                            f(i,j,k+1,5)=feq + (1-omega(k))*p1*fneq1

                            !15
                            
                            feq=(rho(i+1,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i+1,j,k+1) + 2*pzz(i+1,j,k+1) - 3*pyy(i+1,j,k+1)*u(i+1,j,k+1) + 6*pzz(i+1,j,k+1)*u(i+1,j,k+1) - 6*pxy(i+1,j,k+1)*v(i+1,j,k+1) - 6*pyz(i+1,j,k+1)*v(i+1,j,k+1) - 3*pyy(i+1,j,k+1)*w(i+1,j,k+1) + 6*pxz(i+1,j,k+1)*(1 + 2*u(i+1,j,k+1) + 2*w(i+1,j,k+1)) + pxx(i+1,j,k+1)*(2 + 6*w(i+1,j,k+1))))/2.
                            f(i+1,j,k+1,15)=feq + (1-omega(k))*p2*fneq1

                            !17
                            
                            feq=(rho(i-1,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i-1,j,k+1) + 2*pzz(i-1,j,k+1) + 3*pyy(i-1,j,k+1)*u(i-1,j,k+1) - 6*pzz(i-1,j,k+1)*u(i-1,j,k+1) + 6*pxy(i-1,j,k+1)*v(i-1,j,k+1) - 6*pyz(i-1,j,k+1)*v(i-1,j,k+1) + 6*pxz(i-1,j,k+1)*(-1 + 2*u(i-1,j,k+1) - 2*w(i-1,j,k+1)) - 3*pyy(i-1,j,k+1)*w(i-1,j,k+1) + pxx(i-1,j,k+1)*(2 + 6*w(i-1,j,k+1))))/2.
                            f(i-1,j,k+1,17)=feq + (1-omega(k))*p2*fneq1

                            !11
                            
                            feq=(rho(i,j+1,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(-3*pxx(i,j+1,k+1)*(1 + 3*v(i,j+1,k+1) + 3*w(i,j+1,k+1)))/2. + 3*(pyy(i,j+1,k+1) + pzz(i,j+1,k+1) - 3*pxy(i,j+1,k+1)*u(i,j+1,k+1) - 3*pxz(i,j+1,k+1)*u(i,j+1,k+1) + 3*pzz(i,j+1,k+1)*v(i,j+1,k+1) + 3*pyy(i,j+1,k+1)*w(i,j+1,k+1) + pyz(i,j+1,k+1)*(3 + 6*v(i,j+1,k+1) + 6*w(i,j+1,k+1)))
                            f(i,j+1,k+1,11)=feq +(1-omega(k))*p2*fneq1
                            
                            !14
                            feq=(rho(i,j-1,k+1)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*pxx(i,j-1,k+1)*(-1 + 3*v(i,j-1,k+1) - 3*w(i,j-1,k+1)))/2. + 3*(pyy(i,j-1,k+1) + pzz(i,j-1,k+1) + 3*pxy(i,j-1,k+1)*u(i,j-1,k+1) - 3*pxz(i,j-1,k+1)*u(i,j-1,k+1) - 3*pzz(i,j-1,k+1)*v(i,j-1,k+1) + pyz(i,j-1,k+1)*(-3 + 6*v(i,j-1,k+1) - 6*w(i,j-1,k+1)) + 3*pyy(i,j-1,k+1)*w(i,j-1,k+1))
                            f(i,j-1,k+1,14)=feq + (1-omega(k))*p2*fneq1

                            !19
                            feq=(rho(i+1,j+1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i+1,j+1,k+1) + (pyy(i+1,j+1,k+1) + 3*pyz(i+1,j+1,k+1) + pzz(i+1,j+1,k+1))*(1 + 3*u(i+1,j+1,k+1)) + pxy(i+1,j+1,k+1)*(3 + 6*u(i+1,j+1,k+1)) + pxz(i+1,j+1,k+1)*(3 + 6*u(i+1,j+1,k+1)) + 3*pxx(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 6*pxy(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 9*pxz(i+1,j+1,k+1)*v(i+1,j+1,k+1) &
                            + 6*pyz(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 3*pzz(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 3*(pxx(i+1,j+1,k+1) + 3*pxy(i+1,j+1,k+1) + 2*pxz(i+1,j+1,k+1) + pyy(i+1,j+1,k+1) + 2*pyz(i+1,j+1,k+1))*w(i+1,j+1,k+1))
                            f(i+1,j+1,k+1,19)=feq + (1-omega(k))*fneq1*p3

                            !21
                            feq=(rho(i+1,j-1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i+1,j-1,k+1) - 3*pxy(i+1,j-1,k+1)*(1 + 2*u(i+1,j-1,k+1)) + (pyy(i+1,j-1,k+1) - 3*pyz(i+1,j-1,k+1) + pzz(i+1,j-1,k+1))*(1 + 3*u(i+1,j-1,k+1)) + pxz(i+1,j-1,k+1)*(3 + 6*u(i+1,j-1,k+1)) - 3*pxx(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 6*pxy(i+1,j-1,k+1)*v(i+1,j-1,k+1) &
                            - 9*pxz(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 6*pyz(i+1,j-1,k+1)*v(i+1,j-1,k+1) - 3*pzz(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 3*(pxx(i+1,j-1,k+1) - 3*pxy(i+1,j-1,k+1) + 2*pxz(i+1,j-1,k+1) + pyy(i+1,j-1,k+1) - 2*pyz(i+1,j-1,k+1))*w(i+1,j-1,k+1))
                            f(i+1,j-1,k+1,21)=feq + (1-omega(k))*fneq1*p3

                            !23
                            feq=(rho(i-1,j-1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i-1,j-1,k+1) + 3*pxy(i-1,j-1,k+1) - 3*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 3*pyz(i-1,j-1,k+1) + pzz(i-1,j-1,k+1) - 3*(2*pxy(i-1,j-1,k+1) - 2*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 3*pyz(i-1,j-1,k+1) + pzz(i-1,j-1,k+1))*u(i-1,j-1,k+1) &
                            - 3*pxx(i-1,j-1,k+1)*v(i-1,j-1,k+1) - 6*pxy(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 9*pxz(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 6*pyz(i-1,j-1,k+1)*v(i-1,j-1,k+1) - 3*pzz(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 3*(pxx(i-1,j-1,k+1) + 3*pxy(i-1,j-1,k+1) - 2*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 2*pyz(i-1,j-1,k+1))*w(i-1,j-1,k+1))
                            f(i-1,j-1,k+1,23)=feq + (1-omega(k))*fneq1*p3

                            !26
                            feq=(rho(i-1,j+1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i-1,j+1,k+1) - (pyy(i-1,j+1,k+1) + 3*pyz(i-1,j+1,k+1) + pzz(i-1,j+1,k+1))*(-1 + 3*u(i-1,j+1,k+1)) + pxy(i-1,j+1,k+1)*(-3 + 6*u(i-1,j+1,k+1)) + pxz(i-1,j+1,k+1)*(-3 + 6*u(i-1,j+1,k+1)) + 3*pxx(i-1,j+1,k+1)*v(i-1,j+1,k+1) - 6*pxy(i-1,j+1,k+1)*v(i-1,j+1,k+1) &
                            - 9*pxz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 6*pyz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 3*pzz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 3*(pxx(i-1,j+1,k+1) - 3*pxy(i-1,j+1,k+1) - 2*pxz(i-1,j+1,k+1) + pyy(i-1,j+1,k+1) + 2*pyz(i-1,j+1,k+1))*w(i-1,j+1,k+1))
                            f(i-1,j+1,k+1,26)=feq + (1-omega(k))*fneq1*p3
						elseif(k.eq.nz)then
							!6
							w(i,j,k)=w(i,j,nz-1)
                            feq=(rho(i,j,k-1)*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                            fneq1=(3*(-pyy(i,j,k-1) + 2*pzz(i,j,k-1) + 6*pxz(i,j,k-1)*u(i,j,k-1) + 6*pyz(i,j,k-1)*v(i,j,k-1) + 3*pyy(i,j,k-1)*w(i,j,k-1) + pxx(i,j,k-1)*(-1 + 3*w(i,j,k-1))))/2.
                            f(i,j,k-1,6)=feq + (1-omega(k))*fneq1*p1
                            
                            !16
                            feq=(rho(i-1,j,k-1)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*(-pyy(i-1,j,k-1) + 2*pzz(i-1,j,k-1) + 3*pyy(i-1,j,k-1)*u(i-1,j,k-1) - 6*pzz(i-1,j,k-1)*u(i-1,j,k-1) + 6*pxy(i-1,j,k-1)*v(i-1,j,k-1) + 6*pyz(i-1,j,k-1)*v(i-1,j,k-1) + pxx(i-1,j,k-1)*(2 - 6*w(i-1,j,k-1)) + 3*pyy(i-1,j,k-1)*w(i-1,j,k-1) - 6*pxz(i-1,j,k-1)*(-1 + 2*u(i-1,j,k-1) + 2*w(i-1,j,k-1))))/2.
                            f(i-1,j,k-1,16)=feq + (1-omega(k))*fneq1*p2
                            
                            !18
                            feq=(rho(i+1,j,k-1)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(-3*(pyy(i+1,j,k-1) - 2*pzz(i+1,j,k-1) + 3*pyy(i+1,j,k-1)*u(i+1,j,k-1) - 6*pzz(i+1,j,k-1)*u(i+1,j,k-1) + 6*pxy(i+1,j,k-1)*v(i+1,j,k-1) - 6*pyz(i+1,j,k-1)*v(i+1,j,k-1) + 6*pxz(i+1,j,k-1)*(1 + 2*u(i+1,j,k-1) - 2*w(i+1,j,k-1)) - 3*pyy(i+1,j,k-1)*w(i+1,j,k-1) + pxx(i+1,j,k-1)*(-2 + 6*w(i+1,j,k-1))))/2.
                            f(i+1,j,k-1,18)=feq + (1-omega(k))*fneq1*p2
                            
                            !12
                            feq=(rho(i,j-1,k-1)*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*pxx(i,j-1,k-1)*(-1 + 3*v(i,j-1,k-1) + 3*w(i,j-1,k-1)))/2. + 3*(pyy(i,j-1,k-1) + pzz(i,j-1,k-1) + 3*pxy(i,j-1,k-1)*u(i,j-1,k-1) + 3*pxz(i,j-1,k-1)*u(i,j-1,k-1) - 3*pzz(i,j-1,k-1)*v(i,j-1,k-1) + pyz(i,j-1,k-1)*(3 - 6*v(i,j-1,k-1) - 6*w(i,j-1,k-1)) - 3*pyy(i,j-1,k-1)*w(i,j-1,k-1))
                            f(i,j-1,k-1,12)=feq + (1-omega(k))*p2*fneq1

                            !13
                            
                            feq=(rho(i,j+1,k-1)*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                            fneq1=(-3*(pxx(i,j+1,k-1) - 2*pyy(i,j+1,k-1) + 6*pyz(i,j+1,k-1) - 2*pzz(i,j+1,k-1) + 6*pxy(i,j+1,k-1)*u(i,j+1,k-1) - 6*pxz(i,j+1,k-1)*u(i,j+1,k-1) + 3*pxx(i,j+1,k-1)*v(i,j+1,k-1) + 12*pyz(i,j+1,k-1)*v(i,j+1,k-1) - 6*pzz(i,j+1,k-1)*v(i,j+1,k-1) - 3*pxx(i,j+1,k-1)*w(i,j+1,k-1) + 6*pyy(i,j+1,k-1)*w(i,j+1,k-1) - 12*pyz(i,j+1,k-1)*w(i,j+1,k-1)))/2.
                            f(i,j+1,k-1,13)=feq + (1-omega(k))*p2*fneq1
                            
                            !20
                            feq=(rho(i-1,j-1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*((pyy(i-1,j-1,k-1) + 3*pyz(i-1,j-1,k-1) + pzz(i-1,j-1,k-1))*(-1 + 3*u(i-1,j-1,k-1)) + pxy(i-1,j-1,k-1)*(-3 + 6*u(i-1,j-1,k-1)) + pxz(i-1,j-1,k-1)*(-3 + 6*u(i-1,j-1,k-1)) + 3*(2*pxy(i-1,j-1,k-1) + 3*pxz(i-1,j-1,k-1) + 2*pyz(i-1,j-1,k-1) + pzz(i-1,j-1,k-1))*v(i-1,j-1,k-1) + 3*(3*pxy(i-1,j-1,k-1) + &
                            2*pxz(i-1,j-1,k-1) + pyy(i-1,j-1,k-1) + 2*pyz(i-1,j-1,k-1))*w(i-1,j-1,k-1) + pxx(i-1,j-1,k-1)*(-1 + 3*v(i-1,j-1,k-1) + 3*w(i-1,j-1,k-1)))
                            f(i-1,j-1,k-1,20)=feq + (1-omega(k))*p3*fneq1
                            
                            !22
                            feq=(rho(i-1,j+1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i-1,j+1,k-1) + 3*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 3*pyz(i-1,j+1,k-1) + pzz(i-1,j+1,k-1) - 3*(2*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 3*pyz(i-1,j+1,k-1) + pzz(i-1,j+1,k-1))*u(i-1,j+1,k-1) + pxy(i-1,j+1,k-1)*(-3 + 6*u(i-1,j+1,k-1)) + 3*pxx(i-1,j+1,k-1)*v(i-1,j+1,k-1) &
                            - 6*pxy(i-1,j+1,k-1)*v(i-1,j+1,k-1) + 9*pxz(i-1,j+1,k-1)*v(i-1,j+1,k-1) - 6*pyz(i-1,j+1,k-1)*v(i-1,j+1,k-1) + 3*pzz(i-1,j+1,k-1)*v(i-1,j+1,k-1) - 3*(pxx(i-1,j+1,k-1) - 3*pxy(i-1,j+1,k-1) + 2*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 2*pyz(i-1,j+1,k-1))*w(i-1,j+1,k-1))

                            f(i-1,j+1,k-1,22)=feq + (1-omega(k))*p3*fneq1

                            !24
                            feq=(rho(i+1,j+1,k-1)*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                            fneq1=3*(pxx(i+1,j+1,k-1) - 3*pxz(i+1,j+1,k-1)*(1 + 2*u(i+1,j+1,k-1)) + (pyy(i+1,j+1,k-1) - 3*pyz(i+1,j+1,k-1) + pzz(i+1,j+1,k-1))*(1 + 3*u(i+1,j+1,k-1)) + pxy(i+1,j+1,k-1)*(3 + 6*u(i+1,j+1,k-1)) + 3*pxx(i+1,j+1,k-1)*v(i+1,j+1,k-1) + 6*pxy(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 9*pxz(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 6*pyz(i+1,j+1,k-1)*v(i+1,j+1,k-1) &
                            + 3*pzz(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 3*(pxx(i+1,j+1,k-1) + 3*pxy(i+1,j+1,k-1) - 2*pxz(i+1,j+1,k-1) + pyy(i+1,j+1,k-1) - 2*pyz(i+1,j+1,k-1))*w(i+1,j+1,k-1))
                            f(i+1,j+1,k-1,24)=feq + (1-omega(k))*p3*fneq1

                            !25
                            feq=(rho(i+1,j-1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*(-pxx(i+1,j-1,k-1) - (pyy(i+1,j-1,k-1) + 3*pyz(i+1,j-1,k-1) + pzz(i+1,j-1,k-1))*(1 + 3*u(i+1,j-1,k-1)) + pxy(i+1,j-1,k-1)*(3 + 6*u(i+1,j-1,k-1)) + pxz(i+1,j-1,k-1)*(3 + 6*u(i+1,j-1,k-1)) + 3*pxx(i+1,j-1,k-1)*v(i+1,j-1,k-1) &
                            - 6*pxy(i+1,j-1,k-1)*v(i+1,j-1,k-1) - 9*pxz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 6*pyz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 3*pzz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 3*(pxx(i+1,j-1,k-1) - 3*pxy(i+1,j-1,k-1) - 2*pxz(i+1,j-1,k-1) + pyy(i+1,j-1,k-1) + 2*pyz(i+1,j-1,k-1))*w(i+1,j-1,k-1))
                            f(i+1,j-1,k-1,25)=feq + (1-omega(k))*p3*fneq1
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
						if(k.eq.1)then
							!5
							!w(i,j,k)=uwall + 0.02*sqrt(-2.0*log(ranx(ttt)))*cos(2*3.1415926535897932384626433832795028841971*rany(ttt))
                            feq=(rho(i,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k))))/27.
                            fneq1=(-3*(pxx(i,j,k+1) + pyy(i,j,k+1) - 2*pzz(i,j,k+1) + 6*pxz(i,j,k+1)*u(i,j,k+1) + 6*pyz(i,j,k+1)*v(i,j,k+1) + 3*pxx(i,j,k+1)*w(i,j,k+1) + 3*pyy(i,j,k+1)*w(i,j,k+1)))/2.
                            f(i,j,k+1,5)=feq + (1-omega(k))*p1*fneq1

                            !15
                            
                            feq=(rho(i+1,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i+1,j,k+1) + 2*pzz(i+1,j,k+1) - 3*pyy(i+1,j,k+1)*u(i+1,j,k+1) + 6*pzz(i+1,j,k+1)*u(i+1,j,k+1) - 6*pxy(i+1,j,k+1)*v(i+1,j,k+1) - 6*pyz(i+1,j,k+1)*v(i+1,j,k+1) - 3*pyy(i+1,j,k+1)*w(i+1,j,k+1) + 6*pxz(i+1,j,k+1)*(1 + 2*u(i+1,j,k+1) + 2*w(i+1,j,k+1)) + pxx(i+1,j,k+1)*(2 + 6*w(i+1,j,k+1))))/2.
                            f(i+1,j,k+1,15)=feq + (1-omega(k))*p2*fneq1

                            !17
                            
                            feq=(rho(i-1,j,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*u(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*v(i,j,k)**2*(1 + 3*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*(-pyy(i-1,j,k+1) + 2*pzz(i-1,j,k+1) + 3*pyy(i-1,j,k+1)*u(i-1,j,k+1) - 6*pzz(i-1,j,k+1)*u(i-1,j,k+1) + 6*pxy(i-1,j,k+1)*v(i-1,j,k+1) - 6*pyz(i-1,j,k+1)*v(i-1,j,k+1) + 6*pxz(i-1,j,k+1)*(-1 + 2*u(i-1,j,k+1) - 2*w(i-1,j,k+1)) - 3*pyy(i-1,j,k+1)*w(i-1,j,k+1) + pxx(i-1,j,k+1)*(2 + 6*w(i-1,j,k+1))))/2.
                            f(i-1,j,k+1,17)=feq + (1-omega(k))*p2*fneq1

                            !11
                            
                            feq=(rho(i,j+1,k+1)*(2 + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 3*u(i,j,k)**2*(1 + 3*v(i,j,k) + 3*w(i,j,k)) + 2*v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(-3*pxx(i,j+1,k+1)*(1 + 3*v(i,j+1,k+1) + 3*w(i,j+1,k+1)))/2. + 3*(pyy(i,j+1,k+1) + pzz(i,j+1,k+1) - 3*pxy(i,j+1,k+1)*u(i,j+1,k+1) - 3*pxz(i,j+1,k+1)*u(i,j+1,k+1) + 3*pzz(i,j+1,k+1)*v(i,j+1,k+1) + 3*pyy(i,j+1,k+1)*w(i,j+1,k+1) + pyz(i,j+1,k+1)*(3 + 6*v(i,j+1,k+1) + 6*w(i,j+1,k+1)))
                            f(i,j+1,k+1,11)=feq +(1-omega(k))*p2*fneq1
                            
                            !14
                            feq=(rho(i,j-1,k+1)*(2 + u(i,j,k)**2*(-3 + 9*v(i,j,k) - 9*w(i,j,k)) + 6*w(i,j,k)*(1 + w(i,j,k)) + 6*v(i,j,k)**2*(1 + 3*w(i,j,k)) - 6*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)))))/108.
                            fneq1=(3*pxx(i,j-1,k+1)*(-1 + 3*v(i,j-1,k+1) - 3*w(i,j-1,k+1)))/2. + 3*(pyy(i,j-1,k+1) + pzz(i,j-1,k+1) + 3*pxy(i,j-1,k+1)*u(i,j-1,k+1) - 3*pxz(i,j-1,k+1)*u(i,j-1,k+1) - 3*pzz(i,j-1,k+1)*v(i,j-1,k+1) + pyz(i,j-1,k+1)*(-3 + 6*v(i,j-1,k+1) - 6*w(i,j-1,k+1)) + 3*pyy(i,j-1,k+1)*w(i,j-1,k+1))
                            f(i,j-1,k+1,14)=feq + (1-omega(k))*p2*fneq1

                            !19
                            feq=(rho(i+1,j+1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i+1,j+1,k+1) + (pyy(i+1,j+1,k+1) + 3*pyz(i+1,j+1,k+1) + pzz(i+1,j+1,k+1))*(1 + 3*u(i+1,j+1,k+1)) + pxy(i+1,j+1,k+1)*(3 + 6*u(i+1,j+1,k+1)) + pxz(i+1,j+1,k+1)*(3 + 6*u(i+1,j+1,k+1)) + 3*pxx(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 6*pxy(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 9*pxz(i+1,j+1,k+1)*v(i+1,j+1,k+1) &
                            + 6*pyz(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 3*pzz(i+1,j+1,k+1)*v(i+1,j+1,k+1) + 3*(pxx(i+1,j+1,k+1) + 3*pxy(i+1,j+1,k+1) + 2*pxz(i+1,j+1,k+1) + pyy(i+1,j+1,k+1) + 2*pyz(i+1,j+1,k+1))*w(i+1,j+1,k+1))
                            f(i+1,j+1,k+1,19)=feq + (1-omega(k))*fneq1*p3

                            !21
                            feq=(rho(i+1,j-1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i+1,j-1,k+1) - 3*pxy(i+1,j-1,k+1)*(1 + 2*u(i+1,j-1,k+1)) + (pyy(i+1,j-1,k+1) - 3*pyz(i+1,j-1,k+1) + pzz(i+1,j-1,k+1))*(1 + 3*u(i+1,j-1,k+1)) + pxz(i+1,j-1,k+1)*(3 + 6*u(i+1,j-1,k+1)) - 3*pxx(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 6*pxy(i+1,j-1,k+1)*v(i+1,j-1,k+1) &
                            - 9*pxz(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 6*pyz(i+1,j-1,k+1)*v(i+1,j-1,k+1) - 3*pzz(i+1,j-1,k+1)*v(i+1,j-1,k+1) + 3*(pxx(i+1,j-1,k+1) - 3*pxy(i+1,j-1,k+1) + 2*pxz(i+1,j-1,k+1) + pyy(i+1,j-1,k+1) - 2*pyz(i+1,j-1,k+1))*w(i+1,j-1,k+1))
                            f(i+1,j-1,k+1,21)=feq + (1-omega(k))*fneq1*p3

                            !23
                            feq=(rho(i-1,j-1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) + 9*w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*w(i,j,k)*(1 + w(i,j,k)) - 3*v(i,j,k)*(1 + 3*w(i,j,k)))))/216.
                            fneq1=3*(pxx(i-1,j-1,k+1) + 3*pxy(i-1,j-1,k+1) - 3*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 3*pyz(i-1,j-1,k+1) + pzz(i-1,j-1,k+1) - 3*(2*pxy(i-1,j-1,k+1) - 2*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 3*pyz(i-1,j-1,k+1) + pzz(i-1,j-1,k+1))*u(i-1,j-1,k+1) &
                            - 3*pxx(i-1,j-1,k+1)*v(i-1,j-1,k+1) - 6*pxy(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 9*pxz(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 6*pyz(i-1,j-1,k+1)*v(i-1,j-1,k+1) - 3*pzz(i-1,j-1,k+1)*v(i-1,j-1,k+1) + 3*(pxx(i-1,j-1,k+1) + 3*pxy(i-1,j-1,k+1) - 2*pxz(i-1,j-1,k+1) + pyy(i-1,j-1,k+1) - 2*pyz(i-1,j-1,k+1))*w(i-1,j-1,k+1))
                            f(i-1,j-1,k+1,23)=feq + (1-omega(k))*fneq1*p3

                            !26
                            feq=(rho(i-1,j+1,k+1)*(1 + 3*w(i,j,k)*(1 + w(i,j,k)) + v(i,j,k)**2*(3 + 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) + 9*w(i,j,k)) + v(i,j,k)*(3 + 9*w(i,j,k)*(1 + w(i,j,k))) - 3*u(i,j,k)*(1 + 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 + 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i-1,j+1,k+1) - (pyy(i-1,j+1,k+1) + 3*pyz(i-1,j+1,k+1) + pzz(i-1,j+1,k+1))*(-1 + 3*u(i-1,j+1,k+1)) + pxy(i-1,j+1,k+1)*(-3 + 6*u(i-1,j+1,k+1)) + pxz(i-1,j+1,k+1)*(-3 + 6*u(i-1,j+1,k+1)) + 3*pxx(i-1,j+1,k+1)*v(i-1,j+1,k+1) - 6*pxy(i-1,j+1,k+1)*v(i-1,j+1,k+1) &
                            - 9*pxz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 6*pyz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 3*pzz(i-1,j+1,k+1)*v(i-1,j+1,k+1) + 3*(pxx(i-1,j+1,k+1) - 3*pxy(i-1,j+1,k+1) - 2*pxz(i-1,j+1,k+1) + pyy(i-1,j+1,k+1) + 2*pyz(i-1,j+1,k+1))*w(i-1,j+1,k+1))
                            f(i-1,j+1,k+1,26)=feq + (1-omega(k))*fneq1*p3
						elseif(k.eq.nz)then
							!6
							w(i,j,k)=w(i,j,nz-1)
                            feq=(rho(i,j,k-1)*(2 + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*w(i,j,k)) + v(i,j,k)**2*(-3 + 9*w(i,j,k))))/27.
                            fneq1=(3*(-pyy(i,j,k-1) + 2*pzz(i,j,k-1) + 6*pxz(i,j,k-1)*u(i,j,k-1) + 6*pyz(i,j,k-1)*v(i,j,k-1) + 3*pyy(i,j,k-1)*w(i,j,k-1) + pxx(i,j,k-1)*(-1 + 3*w(i,j,k-1))))/2.
                            f(i,j,k-1,6)=feq + (1-omega(k))*fneq1*p1
                            
                            !16
                            feq=(rho(i-1,j,k-1)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(-2 + 3*v(i,j,k)**2 - 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*(-pyy(i-1,j,k-1) + 2*pzz(i-1,j,k-1) + 3*pyy(i-1,j,k-1)*u(i-1,j,k-1) - 6*pzz(i-1,j,k-1)*u(i-1,j,k-1) + 6*pxy(i-1,j,k-1)*v(i-1,j,k-1) + 6*pyz(i-1,j,k-1)*v(i-1,j,k-1) + pxx(i-1,j,k-1)*(2 - 6*w(i-1,j,k-1)) + 3*pyy(i-1,j,k-1)*w(i-1,j,k-1) - 6*pxz(i-1,j,k-1)*(-1 + 2*u(i-1,j,k-1) + 2*w(i-1,j,k-1))))/2.
                            f(i-1,j,k-1,16)=feq + (1-omega(k))*fneq1*p2
                            
                            !18
                            feq=(rho(i+1,j,k-1)*(2 + u(i,j,k)**2*(6 - 18*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)**2*(-3 + 9*w(i,j,k)) + 3*u(i,j,k)*(2 - 3*v(i,j,k)**2 + 6*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(-3*(pyy(i+1,j,k-1) - 2*pzz(i+1,j,k-1) + 3*pyy(i+1,j,k-1)*u(i+1,j,k-1) - 6*pzz(i+1,j,k-1)*u(i+1,j,k-1) + 6*pxy(i+1,j,k-1)*v(i+1,j,k-1) - 6*pyz(i+1,j,k-1)*v(i+1,j,k-1) + 6*pxz(i+1,j,k-1)*(1 + 2*u(i+1,j,k-1) - 2*w(i+1,j,k-1)) - 3*pyy(i+1,j,k-1)*w(i+1,j,k-1) + pxx(i+1,j,k-1)*(-2 + 6*w(i+1,j,k-1))))/2.
                            f(i+1,j,k-1,18)=feq + (1-omega(k))*fneq1*p2
                            
                            !12
                            feq=(rho(i,j-1,k-1)*(2 + 2*v(i,j,k)**2*(3 - 9*w(i,j,k)) + 6*(-1 + w(i,j,k))*w(i,j,k) + u(i,j,k)**2*(-3 + 9*v(i,j,k) + 9*w(i,j,k)) + 2*v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k))))/108.
                            fneq1=(3*pxx(i,j-1,k-1)*(-1 + 3*v(i,j-1,k-1) + 3*w(i,j-1,k-1)))/2. + 3*(pyy(i,j-1,k-1) + pzz(i,j-1,k-1) + 3*pxy(i,j-1,k-1)*u(i,j-1,k-1) + 3*pxz(i,j-1,k-1)*u(i,j-1,k-1) - 3*pzz(i,j-1,k-1)*v(i,j-1,k-1) + pyz(i,j-1,k-1)*(3 - 6*v(i,j-1,k-1) - 6*w(i,j-1,k-1)) - 3*pyy(i,j-1,k-1)*w(i,j-1,k-1))
                            f(i,j-1,k-1,12)=feq + (1-omega(k))*p2*fneq1

                            !13
                            
                            feq=(rho(i,j+1,k-1)*(2 - 6*w(i,j,k) + u(i,j,k)**2*(-3 - 9*v(i,j,k) + 9*w(i,j,k)) + 6*(v(i,j,k) + v(i,j,k)**2*(1 - 3*w(i,j,k)) + 3*v(i,j,k)*(-1 + w(i,j,k))*w(i,j,k) + w(i,j,k)**2)))/108.
                            fneq1=(-3*(pxx(i,j+1,k-1) - 2*pyy(i,j+1,k-1) + 6*pyz(i,j+1,k-1) - 2*pzz(i,j+1,k-1) + 6*pxy(i,j+1,k-1)*u(i,j+1,k-1) - 6*pxz(i,j+1,k-1)*u(i,j+1,k-1) + 3*pxx(i,j+1,k-1)*v(i,j+1,k-1) + 12*pyz(i,j+1,k-1)*v(i,j+1,k-1) - 6*pzz(i,j+1,k-1)*v(i,j+1,k-1) - 3*pxx(i,j+1,k-1)*w(i,j+1,k-1) + 6*pyy(i,j+1,k-1)*w(i,j+1,k-1) - 12*pyz(i,j+1,k-1)*w(i,j+1,k-1)))/2.
                            f(i,j+1,k-1,13)=feq + (1-omega(k))*p2*fneq1
                            
                            !20
                            feq=(rho(i-1,j-1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*((pyy(i-1,j-1,k-1) + 3*pyz(i-1,j-1,k-1) + pzz(i-1,j-1,k-1))*(-1 + 3*u(i-1,j-1,k-1)) + pxy(i-1,j-1,k-1)*(-3 + 6*u(i-1,j-1,k-1)) + pxz(i-1,j-1,k-1)*(-3 + 6*u(i-1,j-1,k-1)) + 3*(2*pxy(i-1,j-1,k-1) + 3*pxz(i-1,j-1,k-1) + 2*pyz(i-1,j-1,k-1) + pzz(i-1,j-1,k-1))*v(i-1,j-1,k-1) + 3*(3*pxy(i-1,j-1,k-1) + &
                            2*pxz(i-1,j-1,k-1) + pyy(i-1,j-1,k-1) + 2*pyz(i-1,j-1,k-1))*w(i-1,j-1,k-1) + pxx(i-1,j-1,k-1)*(-1 + 3*v(i-1,j-1,k-1) + 3*w(i-1,j-1,k-1)))
                            f(i-1,j-1,k-1,20)=feq + (1-omega(k))*p3*fneq1
                            
                            !22
                            feq=(rho(i-1,j+1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 + 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3 + 9*(-1 + w(i,j,k))*w(i,j,k)) - 3*u(i,j,k)*(1 - 3*w(i,j,k) + 3*(v(i,j,k) + v(i,j,k)**2 - 3*v(i,j,k)*w(i,j,k) + w(i,j,k)**2))))/216.
                            fneq1=3*(pxx(i-1,j+1,k-1) + 3*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 3*pyz(i-1,j+1,k-1) + pzz(i-1,j+1,k-1) - 3*(2*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 3*pyz(i-1,j+1,k-1) + pzz(i-1,j+1,k-1))*u(i-1,j+1,k-1) + pxy(i-1,j+1,k-1)*(-3 + 6*u(i-1,j+1,k-1)) + 3*pxx(i-1,j+1,k-1)*v(i-1,j+1,k-1) &
                            - 6*pxy(i-1,j+1,k-1)*v(i-1,j+1,k-1) + 9*pxz(i-1,j+1,k-1)*v(i-1,j+1,k-1) - 6*pyz(i-1,j+1,k-1)*v(i-1,j+1,k-1) + 3*pzz(i-1,j+1,k-1)*v(i-1,j+1,k-1) - 3*(pxx(i-1,j+1,k-1) - 3*pxy(i-1,j+1,k-1) + 2*pxz(i-1,j+1,k-1) + pyy(i-1,j+1,k-1) - 2*pyz(i-1,j+1,k-1))*w(i-1,j+1,k-1))

                            f(i-1,j+1,k-1,22)=feq + (1-omega(k))*p3*fneq1

                            !24
                            feq=(rho(i+1,j+1,k-1)*(1 + 3*v(i,j,k) - 3*w(i,j,k)+ 3*(u(i,j,k) + u(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k) + 3*u(i,j,k)**2*v(i,j,k) + v(i,j,k)**2 + 3*u(i,j,k)*v(i,j,k)**2 - 3*(u(i,j,k) + u(i,j,k)**2 + v(i,j,k) + 3*u(i,j,k)*v(i,j,k) + v(i,j,k)**2)*w(i,j,k) + (1 + 3*u(i,j,k) + 3*v(i,j,k))*w(i,j,k)**2)))/216.
                            fneq1=3*(pxx(i+1,j+1,k-1) - 3*pxz(i+1,j+1,k-1)*(1 + 2*u(i+1,j+1,k-1)) + (pyy(i+1,j+1,k-1) - 3*pyz(i+1,j+1,k-1) + pzz(i+1,j+1,k-1))*(1 + 3*u(i+1,j+1,k-1)) + pxy(i+1,j+1,k-1)*(3 + 6*u(i+1,j+1,k-1)) + 3*pxx(i+1,j+1,k-1)*v(i+1,j+1,k-1) + 6*pxy(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 9*pxz(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 6*pyz(i+1,j+1,k-1)*v(i+1,j+1,k-1) &
                            + 3*pzz(i+1,j+1,k-1)*v(i+1,j+1,k-1) - 3*(pxx(i+1,j+1,k-1) + 3*pxy(i+1,j+1,k-1) - 2*pxz(i+1,j+1,k-1) + pyy(i+1,j+1,k-1) - 2*pyz(i+1,j+1,k-1))*w(i+1,j+1,k-1))
                            f(i+1,j+1,k-1,24)=feq + (1-omega(k))*p3*fneq1

                            !25
                            feq=(rho(i+1,j-1,k-1)*(1 + v(i,j,k)**2*(3 - 9*w(i,j,k)) + u(i,j,k)**2*(3 - 9*v(i,j,k) - 9*w(i,j,k)) + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 - 9*(-1 + w(i,j,k))*w(i,j,k)) + 3*u(i,j,k)*(1 + 3*v(i,j,k)**2 + 3*(-1 + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3 + 9*w(i,j,k)))))/216.
                            fneq1=-3*(-pxx(i+1,j-1,k-1) - (pyy(i+1,j-1,k-1) + 3*pyz(i+1,j-1,k-1) + pzz(i+1,j-1,k-1))*(1 + 3*u(i+1,j-1,k-1)) + pxy(i+1,j-1,k-1)*(3 + 6*u(i+1,j-1,k-1)) + pxz(i+1,j-1,k-1)*(3 + 6*u(i+1,j-1,k-1)) + 3*pxx(i+1,j-1,k-1)*v(i+1,j-1,k-1) &
                            - 6*pxy(i+1,j-1,k-1)*v(i+1,j-1,k-1) - 9*pxz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 6*pyz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 3*pzz(i+1,j-1,k-1)*v(i+1,j-1,k-1) + 3*(pxx(i+1,j-1,k-1) - 3*pxy(i+1,j-1,k-1) - 2*pxz(i+1,j-1,k-1) + pyy(i+1,j-1,k-1) + 2*pyz(i+1,j-1,k-1))*w(i+1,j-1,k-1))
                            f(i+1,j-1,k-1,25)=feq + (1-omega(k))*p3*fneq1
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

            ! periodic g and psi x-y
            g(2,2:ny-1,2:nz-1,1)=g(nx-1,2:ny-1,2:nz-1,1)
            g(nx-1,2:ny-1,2:nz-1,2)=g(2,2:ny-1,2:nz-1,2)
            g(2:nx-1,2,2:nz-1,3)=g(2:nx-1,ny-1,2:nz-1,3)
            g(2:nx-1,ny-1,2:nz-1,4)=g(2:nx-1,2,2:nz-1,4)
            psi(1,2:ny-1,2:nz-1)=psi(nx-1,2:ny-1,2:nz-1)
            psi(nx,2:ny-1,2:nz-1)=psi(2,2:ny-1,2:nz-1)

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
endmodule