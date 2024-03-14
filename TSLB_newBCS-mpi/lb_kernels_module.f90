module lb_kernels

   use vars
   use mpi_template
   !$if _OPENACC
   use openacc
   !$endif
   implicit none

contains

   subroutine moments_TSLB

      implicit none

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
   endsubroutine
   !****************************************************************************!
   subroutine moments_TSLB_multiphase
	  implicit none

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
				  phi(i,j,k) = g(i,j,k,0)+g(i,j,k,1)+g(i,j,k,2)+g(i,j,k,3)+g(i,j,k,4)+g(i,j,k,5)+g(i,j,k,6)
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
   endsubroutine
   !****************************************************************************!
   subroutine fused_TSLB_1c

      implicit none

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
   endsubroutine
   !************************************************************************************************************************************************************************************************************************************!
   subroutine fused_TSLB_force_1c

      implicit none
	   
      !$acc kernels
      !$acc loop collapse(3) private(feq,uu,temp,udotc,us,vs,ws)
      do k=1,nz
         do j=1,ny
            do i=1,nx
               !if(isfluid(i,j,k).eq.1)then
               
			   us=u(i,j,k)+fx
			   vs=v(i,j,k)+fy
			   ws=w(i,j,k)+fz
			   !0
               feq=(-4*rho(i,j,k)*(-2 + 3*us**2 + 3*vs**2 + 3*ws**2))/27.
               fneq1=(-3*(pxx(i,j,k) + pyy(i,j,k) + pzz(i,j,k)))/2.
               f(i,j,k,0)=feq + (1-omega)*fneq1*p0

               !1

               feq=(rho(i,j,k)*(2 + 6*us + 6*us**2 - 3*vs**2 - 9*us*vs**2 - 3*(1 + 3*us)*ws**2))/27.
               fneq1=(3*(2*pxx(i,j,k) - pzz(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - pyy(i,j,k)*(1 + 3*u(i,j,k)) - 6*pxy(i,j,k)*v(i,j,k) - 6*pxz(i,j,k)*w(i,j,k)))/2.
               f(i+1,j,k,1)=feq + (1-omega)*fneq1*p1

               !2
               feq=(rho(i,j,k)*(2 - 3*vs**2 - 3*ws**2 + 3*us*(-2 + 2*us + 3*vs**2 + 3*ws**2)))/27.
               fneq1=(3*(2*pxx(i,j,k) - pzz(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + pyy(i,j,k)*(-1 + 3*u(i,j,k)) + 6*pxy(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*w(i,j,k)))/2.
               f(i-1,j,k,2)=feq + (1-omega)*fneq1*p1

               !3

               feq=(rho(i,j,k)*(2 - 3*us**2*(1 + 3*vs) - 3*ws**2 + 3*vs*(2 + 2*vs - 3*ws**2)))/27.
               fneq1=(-3*(pxx(i,j,k) - 2*pyy(i,j,k) + pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) + 3*pxx(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*w(i,j,k)))/2.
               f(i,j+1,k,3)=feq+ (1-omega)*fneq1*p1

               !4
               feq=(rho(i,j,k)*(2 + us**2*(-3 + 9*vs) - 3*ws**2 + 3*vs*(-2 + 2*vs + 3*ws**2)))/27.
               fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(-1 + 3*v(i,j,k)) + 6*pyz(i,j,k)*w(i,j,k)))/2.
               f(i,j-1,k,4)=feq+ (1-omega)*fneq1*p1

               !7

               feq=(2*rho(i,j,k)*(1 + 3*vs*(1 + vs) + us**2*(3 + 9*vs) + us*(3 + 9*vs*(1 + vs))) - 3*rho(i,j,k)*(1 + 3*us + 3*vs)*ws**2)/108.
               fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) + 6*pyy(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*(1 + 2*u(i,j,k) + 2*v(i,j,k)) + pxx(i,j,k)*(2 + 6*v(i,j,k)) - 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
               f(i+1,j+1,k,7)=feq + (1-omega)*fneq1*p2

               !8
               feq=(2*rho(i,j,k)*(1 + us**2*(3 - 9*vs) + 3*(-1 + vs)*vs + us*(-3 - 9*(-1 + vs)*vs)) + 3*rho(i,j,k)*(-1 + 3*us + 3*vs)*ws**2)/108.
               fneq1=(-3*(-2*pyy(i,j,k) + pzz(i,j,k) + 6*pyy(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*(-1 + 2*u(i,j,k) + 2*v(i,j,k)) + pxx(i,j,k)*(-2 + 6*v(i,j,k)) - 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
               f(i-1,j-1,k,8)=feq + (1-omega)*fneq1*p2

               !10

               feq=(2*rho(i,j,k)*(1 + 3*vs*(1 + vs) + us**2*(3 + 9*vs) - 3*us*(1 + 3*vs*(1 + vs))) + 3*rho(i,j,k)*(-1 + 3*us - 3*vs)*ws**2)/108.
               fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) - 6*pyy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*(-1 + 2*u(i,j,k) - 2*v(i,j,k)) - 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(2 + 6*v(i,j,k)) + 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
               f(i-1,j+1,k,10)=feq+ (1-omega)*fneq1*p2

               !9
               feq=(2*rho(i,j,k)*(1 + us**2*(3 - 9*vs) + 3*(-1 + vs)*vs + us*(3 + 9*(-1 + vs)*vs)) - 3*rho(i,j,k)*(1 + 3*us - 3*vs)*ws**2)/108.
               fneq1=(-3*(-2*pyy(i,j,k) + pzz(i,j,k) - 6*pyy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*(1 + 2*u(i,j,k) - 2*v(i,j,k)) - 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(-2 + 6*v(i,j,k)) + 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
               f(i+1,j-1,k,9)=feq+ (1-omega)*fneq1*p2

               !5

               feq=(rho(i,j,k)*(2 + 6*ws*(1 + ws) - 3*us**2*(1 + 3*ws) - 3*vs**2*(1 + 3*ws)))/27.
               fneq1=(-3*(pxx(i,j,k) + pyy(i,j,k) - 2*pzz(i,j,k) + 6*pxz(i,j,k)*u(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pxx(i,j,k)*w(i,j,k) + 3*pyy(i,j,k)*w(i,j,k)))/2.
               f(i,j,k+1,5)=feq+ (1-omega)*fneq1*p1

               !6
               feq=(rho(i,j,k)*(2 + 6*(-1 + ws)*ws + us**2*(-3 + 9*ws) + vs**2*(-3 + 9*ws)))/27.
               fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 6*pxz(i,j,k)*u(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(-1 + 3*w(i,j,k))))/2.
               f(i,j,k-1,6)=feq+ (1-omega)*fneq1*p1

               !15

               feq=(rho(i,j,k)*(2 + 6*ws*(1 + ws) + 6*us**2*(1 + 3*ws) - 3*vs**2*(1 + 3*ws) + 3*us*(2 - 3*vs**2 + 6*ws*(1 + ws))))/108.
               fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) - 3*pyy(i,j,k)*u(i,j,k) + 6*pzz(i,j,k)*u(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) - 3*pyy(i,j,k)*w(i,j,k) + 6*pxz(i,j,k)*(1 + 2*u(i,j,k) + 2*w(i,j,k)) + pxx(i,j,k)*(2 + 6*w(i,j,k))))/2.
               f(i+1,j,k+1,15)=feq+ (1-omega)*fneq1*p2

               !16
               feq=(rho(i,j,k)*(2 + us**2*(6 - 18*ws) + 6*(-1 + ws)*ws + vs**2*(-3 + 9*ws) + 3*us*(-2 + 3*vs**2 - 6*(-1 + ws)*ws)))/108.
               fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(2 - 6*w(i,j,k)) + 3*pyy(i,j,k)*w(i,j,k) - 6*pxz(i,j,k)*(-1 + 2*u(i,j,k) + 2*w(i,j,k))))/2.
               f(i-1,j,k-1,16)=feq+ (1-omega)*fneq1*p2

               !17

               feq=(rho(i,j,k)*(2 + 6*ws*(1 + ws) + 6*us**2*(1 + 3*ws) - 3*vs**2*(1 + 3*ws) + 3*us*(-2 + 3*vs**2 - 6*ws*(1 + ws))))/108.
               fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*(-1 + 2*u(i,j,k) - 2*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(2 + 6*w(i,j,k))))/2.
               f(i-1,j,k+1,17)=feq+ (1-omega)*fneq1*p2

               !18
               feq=(rho(i,j,k)*(2 + us**2*(6 - 18*ws) + 6*(-1 + ws)*ws + vs**2*(-3 + 9*ws) + 3*us*(2 - 3*vs**2 + 6*(-1 + ws)*ws)))/108.
               fneq1=(-3*(pyy(i,j,k) - 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*(1 + 2*u(i,j,k) - 2*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(-2 + 6*w(i,j,k))))/2.
               f(i+1,j,k-1,18)=feq+ (1-omega)*fneq1*p2

               !11

               feq=(rho(i,j,k)*(2 + 6*ws*(1 + ws) + 6*vs**2*(1 + 3*ws) - 3*us**2*(1 + 3*vs + 3*ws) + 2*vs*(3 + 9*ws*(1 + ws))))/108.
               fneq1=(-3*pxx(i,j,k)*(1 + 3*v(i,j,k) + 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) - 3*pxy(i,j,k)*u(i,j,k) - 3*pxz(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*pyy(i,j,k)*w(i,j,k) + pyz(i,j,k)*(3 + 6*v(i,j,k) + 6*w(i,j,k)))
               f(i,j+1,k+1,11)=feq+ (1-omega)*fneq1*p2

               !12
               feq=(rho(i,j,k)*(2 + 2*vs**2*(3 - 9*ws) + 6*(-1 + ws)*ws + us**2*(-3 + 9*vs + 9*ws) + 2*vs*(-3 - 9*(-1 + ws)*ws)))/108.
               fneq1=(3*pxx(i,j,k)*(-1 + 3*v(i,j,k) + 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) + 3*pxy(i,j,k)*u(i,j,k) + 3*pxz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + pyz(i,j,k)*(3 - 6*v(i,j,k) - 6*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k))
               f(i,j-1,k-1,12)=feq+ (1-omega)*fneq1*p2

               !13

               feq=(rho(i,j,k)*(2 - 6*ws + us**2*(-3 - 9*vs + 9*ws) + 6*(vs + vs**2*(1 - 3*ws) + 3*vs*(-1 + ws)*ws + ws**2)))/108.
               fneq1=(-3*(pxx(i,j,k) - 2*pyy(i,j,k) + 6*pyz(i,j,k) - 2*pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) - 6*pxz(i,j,k)*u(i,j,k) + 3*pxx(i,j,k)*v(i,j,k) + 12*pyz(i,j,k)*v(i,j,k) - 6*pzz(i,j,k)*v(i,j,k) - 3*pxx(i,j,k)*w(i,j,k) + 6*pyy(i,j,k)*w(i,j,k) - 12*pyz(i,j,k)*w(i,j,k)))/2.
               f(i,j+1,k-1,13)=feq+ (1-omega)*fneq1*p2

               !14
               feq=(rho(i,j,k)*(2 + us**2*(-3 + 9*vs - 9*ws) + 6*ws*(1 + ws) + 6*vs**2*(1 + 3*ws) - 6*vs*(1 + 3*ws*(1 + ws))))/108.
               fneq1=(3*pxx(i,j,k)*(-1 + 3*v(i,j,k) - 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) + 3*pxy(i,j,k)*u(i,j,k) - 3*pxz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + pyz(i,j,k)*(-3 + 6*v(i,j,k) - 6*w(i,j,k)) + 3*pyy(i,j,k)*w(i,j,k))
               f(i,j-1,k+1,14)=feq+ (1-omega)*fneq1*p2

               !19
               feq=(rho(i,j,k)*(1 + 3*ws*(1 + ws) + vs**2*(3 + 9*ws) + us**2*(3 + 9*vs + 9*ws) + vs*(3 + 9*ws*(1 + ws)) + 3*us*(1 + 3*ws + 3*(vs + vs**2 + 3*vs*ws+ ws**2))))/216.
               fneq1=3*(pxx(i,j,k) + (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) + 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
               f(i+1,j+1,k+1,19)=feq + (1-omega)*fneq1*p3

               !20
               feq=(rho(i,j,k)*(1 + vs**2*(3 - 9*ws) + us**2*(3 - 9*vs - 9*ws) + 3*(-1 + ws)*ws + vs*(-3 - 9*(-1 + ws)*ws) - 3*us*(1 + 3*vs**2 + 3*(-1 + ws)*ws + vs*(-3 + 9*ws))))/216.
               fneq1=-3*((pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(-1 + 3*u(i,j,k)) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + pxz(i,j,k)*(-3 + 6*u(i,j,k)) + 3*(2*pxy(i,j,k) + 3*pxz(i,j,k) + 2*pyz(i,j,k) + pzz(i,j,k))*v(i,j,k) + 3*(3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k) + pxx(i,j,k)*(-1 + 3*v(i,j,k) + 3*w(i,j,k)))
               f(i-1,j-1,k-1,20)=feq+ (1-omega)*fneq1*p3

               !21
               feq=(rho(i,j,k)*(1 + 3*ws*(1 + ws) + vs**2*(3 + 9*ws) + us**2*(3 - 9*vs + 9*ws) - 3*vs*(1 + 3*ws*(1 + ws)) + 3*us*(1 + 3*vs**2 + 3*ws*(1 + ws) - 3*vs*(1 + 3*ws))))/216.
               fneq1=3*(pxx(i,j,k) - 3*pxy(i,j,k)*(1 + 2*u(i,j,k)) + (pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) - 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
               f(i+1,j-1,k+1,21)=feq+ (1-omega)*fneq1*p3

               !22
               feq=(rho(i,j,k)*(1 + vs**2*(3 - 9*ws) + us**2*(3 + 9*vs - 9*ws) + 3*(-1 + ws)*ws + vs*(3 + 9*(-1 + ws)*ws) - 3*us*(1 - 3*ws + 3*(vs + vs**2 - 3*vs*ws + ws**2))))/216.
               fneq1=3*(pxx(i,j,k) + 3*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k) - 3*(2*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) - 3*(pxx(i,j,k) - 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
               f(i-1,j+1,k-1,22)=feq+ (1-omega)*fneq1*p3

               !23
               feq=(rho(i,j,k)*(1 + 3*ws*(1 + ws) + vs**2*(3 + 9*ws) + us**2*(3 - 9*vs + 9*ws) - 3*vs*(1 + 3*ws*(1 + ws)) - 3*us*(1 + 3*vs**2 + 3*ws*(1 + ws) - 3*vs*(1 + 3*ws))))/216.
               fneq1=3*(pxx(i,j,k) + 3*pxy(i,j,k) - 3*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k) - 3*(2*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k) - 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) + 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
               f(i-1,j-1,k+1,23)=feq+ (1-omega)*fneq1*p3

               !24
               feq=(rho(i,j,k)*(1 + 3*vs - 3*ws+ 3*(us + us**2 + 3*us*vs + 3*us**2*vs + vs**2 + 3*us*vs**2 - 3*(us + us**2 + vs + 3*us*vs + vs**2)*ws + (1 + 3*us + 3*vs)*ws**2)))/216.
               fneq1=3*(pxx(i,j,k) - 3*pxz(i,j,k)*(1 + 2*u(i,j,k)) + (pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) - 3*(pxx(i,j,k) + 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
               f(i+1,j+1,k-1,24)=feq+ (1-omega)*fneq1*p3

               !25
               feq=(rho(i,j,k)*(1 + vs**2*(3 - 9*ws) + us**2*(3 - 9*vs - 9*ws) + 3*(-1 + ws)*ws + vs*(-3 - 9*(-1 + ws)*ws) + 3*us*(1 + 3*vs**2 + 3*(-1 + ws)*ws + vs*(-3 + 9*ws))))/216.
               fneq1=-3*(-pxx(i,j,k) - (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
               f(i+1,j-1,k-1,25)=feq + (1-omega)*fneq1*p3

               !26
               feq=(rho(i,j,k)*(1 + 3*ws*(1 + ws) + vs**2*(3 + 9*ws) + us**2*(3 + 9*vs + 9*ws) + vs*(3 + 9*ws*(1 + ws)) - 3*us*(1 + 3*ws + 3*(vs + vs**2 + 3*vs*ws + ws**2))))/216.
               fneq1=3*(pxx(i,j,k) - (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(-1 + 3*u(i,j,k)) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + pxz(i,j,k)*(-3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
               f(i-1,j+1,k+1,26)=feq + (1-omega)*fneq1*p3
               !endif
            enddo
         enddo
      enddo
      !$acc end kernels
   endsubroutine
   !************************************************************************************************************************************************************************************************************************************!
   subroutine fused_TSLB_multiphase
	
	implicit none
	real*4 :: grad_fix,grad_fiy,grad_fiz,mod_grad,curvature,sigma,sharp_c,ffx_s,ffy_s,ffz_s,indicator,feqs
	!$acc kernels
          !$acc loop collapse(3) private(feq,feqs,curvature,ffx_s,ffy_s,ffz_s,grad_fix,grad_fiy,grad_fiz,mod_grad,curvature,sigma,sharp_c,ffx_s,ffy_s,ffz_s,indicator)
          do k=1,nz
              do j=1,ny
                  do i=1,nx
                      if(isfluid(i,j,k).eq.1)then
                          ! surface tension caluclation
                          curvature=0.0
                          indicator=(phi(i,j,k)*(1-phi(i,j,k))*(normx(i,j,k)**2+normy(i,j,k)**2+normz(i,j,k)**2))
                          if(indicator .gt. 0.2)then
                              curvature=-(3.0*(p1*(normx(i+1,j,k)-normx(i-1,j,k)) + &
                              p2*(normx(i+1,j+1,k)+normx(i+1,j-1,k)+normx(i+1,j,k+1)+normx(i+1,j,k-1)-(normx(i-1,j+1,k)+normx(i-1,j-1,k)+normx(i-1,j,k+1)+normx(i-1,j,k-1)) ) + &
                              p3*(normx(i+1,j+1,k+1)+normx(i+1,j-1,k-1)+normx(i+1,j-1,k+1)+normx(i+1,j+1,k-1) - (normx(i-1,j+1,k+1)+normx(i-1,j-1,k-1)+normx(i-1,j-1,k+1)+normx(i-1,j+1,k-1))) ) + &
                                            
											   3.0*(p1*(normy(i,j+1,k)-normy(i,j-1,k)) + &
                              p2*(normy(i+1,j+1,k)+normy(i-1,j+1,k)+normy(i,j+1,k+1)+normy(i,j+1,k-1)-(normy(i+1,j-1,k)+normy(i-1,j-1,k)+normy(i,j-1,k+1)+normy(i,j-1,k-1)) ) + &
                              p3*(normy(i+1,j+1,k+1)+normy(i-1,j+1,k-1)+normy(i+1,j+1,k-1)+normy(i-1,j+1,k+1)-(normy(i+1,j-1,k+1)+normy(i-1,j-1,k-1)+normy(i+1,j-1,k-1)+normy(i-1,j-1,k+1)))) + &
                                            
											   3.0*(p1*(normz(i,j,k+1)-normz(i,j,k-1)) + &
                              p2*(normz(i+1,j,k+1)+normz(i-1,j,k+1)+normz(i,j+1,k+1)+normz(i,j-1,k+1)-(normz(i+1,j,k-1)+normz(i-1,j,k-1)+normz(i,j+1,k-1)+normz(i,j-1,k-1)) ) + &
                              p3*(normz(i+1,j+1,k+1)+normz(i-1,j-1,k+1)+normz(i+1,j-1,k+1)+normz(i-1,j+1,k+1)-(normz(i+1,j+1,k-1)+normz(i-1,j-1,k-1)+normz(i+1,j-1,k-1)+normz(i-1,j+1,k-1)))))
             
                          endif
                          ffx_s=sigma*curvature*normx(i,j,k)
                          ffy_s=sigma*curvature*normy(i,j,k)
                          ffz_s=sigma*curvature*normz(i,j,k)
                          
						  us=u(i,j,k)+ffx_s+fx
                          vs=v(i,j,k)+ffy_s+fy
                          ws=w(i,j,k)+ffz_s+fz
                          !0
						  feq=(-4*rho(i,j,k)*(-2 + 3*us**2 + 3*vs**2 + 3*ws**2))/27.
                          fneq1=(-3*(pxx(i,j,k) + pyy(i,j,k) + pzz(i,j,k)))/2.
                          f(i,j,k,0)=feq + (1-omega)*fneq1*p0 
                          
                          !1
                          feq=(rho(i,j,k)*(2 + 6*us + 6*us**2 - 3*vs**2 - 9*us*vs**2 - 3*(1 + 3*us)*ws**2))/27.
                          fneq1=(3*(2*pxx(i,j,k) - pzz(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - pyy(i,j,k)*(1 + 3*u(i,j,k)) - 6*pxy(i,j,k)*v(i,j,k) - 6*pxz(i,j,k)*w(i,j,k)))/2.
                          f(i+1,j,k,1)=feq + (1-omega)*fneq1*p1 !+ (feqs-feq)
                          
                          !2
                          feq=(rho(i,j,k)*(2 - 3*vs**2 - 3*ws**2 + 3*us*(-2 + 2*us + 3*vs**2 + 3*ws**2)))/27.
                          fneq1=(3*(2*pxx(i,j,k) - pzz(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + pyy(i,j,k)*(-1 + 3*u(i,j,k)) + 6*pxy(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*w(i,j,k)))/2.
                          f(i-1,j,k,2)=feq + (1-omega)*fneq1*p1 !- (feqs-feq)
                          
                          !3
                          feq=(rho(i,j,k)*(2 - 3*us**2*(1 + 3*vs) - 3*ws**2 + 3*vs*(2 + 2*vs - 3*ws**2)))/27.
                          fneq1=(-3*(pxx(i,j,k) - 2*pyy(i,j,k) + pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) + 3*pxx(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i,j+1,k,3)=feq + (1-omega)*fneq1*p1 !+ (feqs-feq)
                          
                          !4
                          feq=(rho(i,j,k)*(2 + us**2*(-3 + 9*vs) - 3*ws**2 + 3*vs*(-2 + 2*vs + 3*ws**2)))/27.
                          fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(-1 + 3*v(i,j,k)) + 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i,j-1,k,4)=feq + (1-omega)*fneq1*p1 !- (feqs-feq)
                          
                          !7
                          feq=(2*rho(i,j,k)*(1 + 3*vs*(1 + vs) + us**2*(3 + 9*vs) + us*(3 + 9*vs*(1 + vs))) - 3*rho(i,j,k)*(1 + 3*us + 3*vs)*ws**2)/108.
                          fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) + 6*pyy(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*(1 + 2*u(i,j,k) + 2*v(i,j,k)) + pxx(i,j,k)*(2 + 6*v(i,j,k)) - 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i+1,j+1,k,7)=feq + (1-omega)*fneq1*p2 !+ (feqs-feq)
                          
                          !8
                          feq=(2*rho(i,j,k)*(1 + us**2*(3 - 9*vs) + 3*(-1 + vs)*vs + us*(-3 - 9*(-1 + vs)*vs)) + 3*rho(i,j,k)*(-1 + 3*us + 3*vs)*ws**2)/108.
                          fneq1=(-3*(-2*pyy(i,j,k) + pzz(i,j,k) + 6*pyy(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*(-1 + 2*u(i,j,k) + 2*v(i,j,k)) + pxx(i,j,k)*(-2 + 6*v(i,j,k)) - 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i-1,j-1,k,8)=feq + (1-omega)*fneq1*p2 !- (feqs-feq)
                          
                          !10
                          feq=(2*rho(i,j,k)*(1 + 3*vs*(1 + vs) + us**2*(3 + 9*vs) - 3*us*(1 + 3*vs*(1 + vs))) + 3*rho(i,j,k)*(-1 + 3*us - 3*vs)*ws**2)/108.
                          fneq1=(3*(2*pyy(i,j,k) - pzz(i,j,k) - 6*pyy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*(-1 + 2*u(i,j,k) - 2*v(i,j,k)) - 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(2 + 6*v(i,j,k)) + 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i-1,j+1,k,10)=feq + (1-omega)*fneq1*p2 !+ (feqs-feq)
                          
                          !9
                          feq=(2*rho(i,j,k)*(1 + us**2*(3 - 9*vs) + 3*(-1 + vs)*vs + us*(3 + 9*(-1 + vs)*vs)) - 3*rho(i,j,k)*(1 + 3*us - 3*vs)*ws**2)/108.
                          fneq1=(-3*(-2*pyy(i,j,k) + pzz(i,j,k) - 6*pyy(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*(1 + 2*u(i,j,k) - 2*v(i,j,k)) - 3*pzz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(-2 + 6*v(i,j,k)) + 6*pxz(i,j,k)*w(i,j,k) - 6*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i+1,j-1,k,9)=feq + (1-omega)*fneq1*p2 !+ (feqs-feq)

                          !5
                          feq=(rho(i,j,k)*(2 + 6*ws*(1 + ws) - 3*us**2*(1 + 3*ws) - 3*vs**2*(1 + 3*ws)))/27.
                          fneq1=(-3*(pxx(i,j,k) + pyy(i,j,k) - 2*pzz(i,j,k) + 6*pxz(i,j,k)*u(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pxx(i,j,k)*w(i,j,k) + 3*pyy(i,j,k)*w(i,j,k)))/2.
                          f(i,j,k+1,5)=feq + (1-omega)*fneq1*p1 !+ (feqs-feq)
                          
                          !6
                          feq=(rho(i,j,k)*(2 + 6*(-1 + ws)*ws + us**2*(-3 + 9*ws) + vs**2*(-3 + 9*ws)))/27.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 6*pxz(i,j,k)*u(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(-1 + 3*w(i,j,k))))/2.
                          f(i,j,k-1,6)=feq + (1-omega)*fneq1*p1 !- (feqs-feq)

                          !15
                          feq=(rho(i,j,k)*(2 + 6*ws*(1 + ws) + 6*us**2*(1 + 3*ws) - 3*vs**2*(1 + 3*ws) + 3*us*(2 - 3*vs**2 + 6*ws*(1 + ws))))/108.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) - 3*pyy(i,j,k)*u(i,j,k) + 6*pzz(i,j,k)*u(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) - 3*pyy(i,j,k)*w(i,j,k) + 6*pxz(i,j,k)*(1 + 2*u(i,j,k) + 2*w(i,j,k)) + pxx(i,j,k)*(2 + 6*w(i,j,k))))/2.
                          f(i+1,j,k+1,15)=feq + (1-omega)*fneq1*p2 !+ (feqs-feq) 
                          
                          !16
                          feq=(rho(i,j,k)*(2 + us**2*(6 - 18*ws) + 6*(-1 + ws)*ws + vs**2*(-3 + 9*ws) + 3*us*(-2 + 3*vs**2 - 6*(-1 + ws)*ws)))/108.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + pxx(i,j,k)*(2 - 6*w(i,j,k)) + 3*pyy(i,j,k)*w(i,j,k) - 6*pxz(i,j,k)*(-1 + 2*u(i,j,k) + 2*w(i,j,k))))/2.
                          f(i-1,j,k-1,16)=feq + (1-omega)*fneq1*p2 !- (feqs-feq)

                          !17
                          feq=(rho(i,j,k)*(2 + 6*ws*(1 + ws) + 6*us**2*(1 + 3*ws) - 3*vs**2*(1 + 3*ws) + 3*us*(-2 + 3*vs**2 - 6*ws*(1 + ws))))/108.
                          fneq1=(3*(-pyy(i,j,k) + 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*(-1 + 2*u(i,j,k) - 2*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(2 + 6*w(i,j,k))))/2.
                          f(i-1,j,k+1,17)=feq + (1-omega)*fneq1*p2 !+ (feqs-feq)
                          
                          !18
                          feq=(rho(i,j,k)*(2 + us**2*(6 - 18*ws) + 6*(-1 + ws)*ws + vs**2*(-3 + 9*ws) + 3*us*(2 - 3*vs**2 + 6*(-1 + ws)*ws)))/108.
                          fneq1=(-3*(pyy(i,j,k) - 2*pzz(i,j,k) + 3*pyy(i,j,k)*u(i,j,k) - 6*pzz(i,j,k)*u(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 6*pxz(i,j,k)*(1 + 2*u(i,j,k) - 2*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k) + pxx(i,j,k)*(-2 + 6*w(i,j,k))))/2.
                          f(i+1,j,k-1,18)=feq + (1-omega)*fneq1*p2 !+ (feqs-feq)

                          !11
                          feq=(rho(i,j,k)*(2 + 6*ws*(1 + ws) + 6*vs**2*(1 + 3*ws) - 3*us**2*(1 + 3*vs + 3*ws) + 2*vs*(3 + 9*ws*(1 + ws))))/108.
                          fneq1=(-3*pxx(i,j,k)*(1 + 3*v(i,j,k) + 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) - 3*pxy(i,j,k)*u(i,j,k) - 3*pxz(i,j,k)*u(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*pyy(i,j,k)*w(i,j,k) + pyz(i,j,k)*(3 + 6*v(i,j,k) + 6*w(i,j,k)))
                          f(i,j+1,k+1,11)=feq + (1-omega)*fneq1*p2 !+ (feqs-feq)
                          
                          !12
                          feq=(rho(i,j,k)*(2 + 2*vs**2*(3 - 9*ws) + 6*(-1 + ws)*ws + us**2*(-3 + 9*vs + 9*ws) + 2*vs*(-3 - 9*(-1 + ws)*ws)))/108.
                          fneq1=(3*pxx(i,j,k)*(-1 + 3*v(i,j,k) + 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) + 3*pxy(i,j,k)*u(i,j,k) + 3*pxz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + pyz(i,j,k)*(3 - 6*v(i,j,k) - 6*w(i,j,k)) - 3*pyy(i,j,k)*w(i,j,k))
                          f(i,j-1,k-1,12)=feq + (1-omega)*fneq1*p2 !- (feqs-feq)

                          !13
                          feq=(rho(i,j,k)*(2 - 6*ws + us**2*(-3 - 9*vs + 9*ws) + 6*(vs + vs**2*(1 - 3*ws) + 3*vs*(-1 + ws)*ws + ws**2)))/108.
                          fneq1=(-3*(pxx(i,j,k) - 2*pyy(i,j,k) + 6*pyz(i,j,k) - 2*pzz(i,j,k) + 6*pxy(i,j,k)*u(i,j,k) - 6*pxz(i,j,k)*u(i,j,k) + 3*pxx(i,j,k)*v(i,j,k) + 12*pyz(i,j,k)*v(i,j,k) - 6*pzz(i,j,k)*v(i,j,k) - 3*pxx(i,j,k)*w(i,j,k) + 6*pyy(i,j,k)*w(i,j,k) - 12*pyz(i,j,k)*w(i,j,k)))/2.
                          f(i,j+1,k-1,13)=feq + (1-omega)*fneq1*p2 !+ (feqs-feq)
                          
                          !14
                          feq=(rho(i,j,k)*(2 + us**2*(-3 + 9*vs - 9*ws) + 6*ws*(1 + ws) + 6*vs**2*(1 + 3*ws) - 6*vs*(1 + 3*ws*(1 + ws))))/108.
                          fneq1=(3*pxx(i,j,k)*(-1 + 3*v(i,j,k) - 3*w(i,j,k)))/2. + 3*(pyy(i,j,k) + pzz(i,j,k) + 3*pxy(i,j,k)*u(i,j,k) - 3*pxz(i,j,k)*u(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + pyz(i,j,k)*(-3 + 6*v(i,j,k) - 6*w(i,j,k)) + 3*pyy(i,j,k)*w(i,j,k))
                          f(i,j-1,k+1,14)=feq + (1-omega)*fneq1*p2 !+ (feqs-feq)

                          !19
                          feq=(rho(i,j,k)*(1 + 3*ws*(1 + ws) + vs**2*(3 + 9*ws) + us**2*(3 + 9*vs + 9*ws) + vs*(3 + 9*ws*(1 + ws)) + 3*us*(1 + 3*ws + 3*(vs + vs**2 + 3*vs*ws + ws**2))))/216.
                          fneq1=3*(pxx(i,j,k) + (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) + 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j+1,k+1,19)=feq + (1-omega)*fneq1*p3 !+ (feqs-feq)
                          
                          !20
                          feq=(rho(i,j,k)*(1 + vs**2*(3 - 9*ws) + us**2*(3 - 9*vs - 9*ws) + 3*(-1 + ws)*ws + vs*(-3 - 9*(-1 + ws)*ws) - 3*us*(1 + 3*vs**2 + 3*(-1 + ws)*ws + vs*(-3 + 9*ws))))/216.
                          fneq1=-3*((pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(-1 + 3*u(i,j,k)) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + pxz(i,j,k)*(-3 + 6*u(i,j,k)) + 3*(2*pxy(i,j,k) + 3*pxz(i,j,k) + 2*pyz(i,j,k) + pzz(i,j,k))*v(i,j,k) + 3*(3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k) + pxx(i,j,k)*(-1 + 3*v(i,j,k) + 3*w(i,j,k)))
                          f(i-1,j-1,k-1,20)=feq + (1-omega)*fneq1*p3 !- (feqs-feq)

                          !21
                          feq=(rho(i,j,k)*(1 + 3*ws*(1 + ws) + vs**2*(3 + 9*ws) + us**2*(3 - 9*vs + 9*ws) - 3*vs*(1 + 3*ws*(1 + ws)) + 3*us*(1 + 3*vs**2 + 3*ws*(1 + ws) - 3*vs*(1 + 3*ws))))/216.
						  fneq1=3*(pxx(i,j,k) - 3*pxy(i,j,k)*(1 + 2*u(i,j,k)) + (pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) - 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j-1,k+1,21)=feq + (1-omega)*fneq1*p3 !+ (feqs-feq)
                          
                          !22
                          feq=(rho(i,j,k)*(1 + vs**2*(3 - 9*ws) + us**2*(3 + 9*vs - 9*ws) + 3*(-1 + ws)*ws + vs*(3 + 9*(-1 + ws)*ws) - 3*us*(1 - 3*ws + 3*(vs + vs**2 - 3*vs*ws + ws**2))))/216.
                          fneq1=3*(pxx(i,j,k) + 3*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k) - 3*(2*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) - 3*(pxx(i,j,k) - 3*pxy(i,j,k) + 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i-1,j+1,k-1,22)=feq+ (1-omega)*fneq1*p3 !- (feqs-feq)

                          !23
                          feq=(rho(i,j,k)*(1 + 3*ws*(1 + ws) + vs**2*(3 + 9*ws) + us**2*(3 - 9*vs + 9*ws) - 3*vs*(1 + 3*ws*(1 + ws)) - 3*us*(1 + 3*vs**2 + 3*ws*(1 + ws) - 3*vs*(1 + 3*ws))))/216.
                          fneq1=3*(pxx(i,j,k) + 3*pxy(i,j,k) - 3*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k) - 3*(2*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k) - 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) + 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) - 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) + 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i-1,j-1,k+1,23)=feq + (1-omega)*fneq1*p3 !+ (feqs-feq)

                          !24
                          feq=(rho(i,j,k)*(1 + 3*vs - 3*ws+ 3*(us + us**2 + 3*us*vs + 3*us**2*vs + vs**2 + 3*us*vs**2 - 3*(us + us**2 + vs + 3*us*vs + vs**2)*ws + (1 + 3*us + 3*vs)*ws**2)))/216.
                          fneq1=3*(pxx(i,j,k) - 3*pxz(i,j,k)*(1 + 2*u(i,j,k)) + (pyy(i,j,k) - 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) + 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) - 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) - 3*(pxx(i,j,k) + 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) - 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j+1,k-1,24)=feq + (1-omega)*fneq1*p3 !+ (feqs-feq)

                          !25
                          feq=(rho(i,j,k)*(1 + vs**2*(3 - 9*ws) + us**2*(3 - 9*vs - 9*ws) + 3*(-1 + ws)*ws + vs*(-3 - 9*(-1 + ws)*ws) + 3*us*(1 + 3*vs**2 + 3*(-1 + ws)*ws + vs*(-3 + 9*ws))))/216.
                          fneq1=-3*(-pxx(i,j,k) - (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(1 + 3*u(i,j,k)) + pxy(i,j,k)*(3 + 6*u(i,j,k)) + pxz(i,j,k)*(3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                          f(i+1,j-1,k-1,25)=feq + (1-omega)*fneq1*p3 !+ (feqs-feq)

                          !26
                          feq=(rho(i,j,k)*(1 + 3*ws*(1 + ws) + vs**2*(3 + 9*ws) + us**2*(3 + 9*vs + 9*ws) + vs*(3 + 9*ws*(1 + ws)) - 3*us*(1 + 3*ws + 3*(vs + vs**2 + 3*vs*ws + ws**2))))/216.
                          fneq1=3*(pxx(i,j,k) - (pyy(i,j,k) + 3*pyz(i,j,k) + pzz(i,j,k))*(-1 + 3*u(i,j,k)) + pxy(i,j,k)*(-3 + 6*u(i,j,k)) + pxz(i,j,k)*(-3 + 6*u(i,j,k)) + 3*pxx(i,j,k)*v(i,j,k) - 6*pxy(i,j,k)*v(i,j,k) - 9*pxz(i,j,k)*v(i,j,k) + 6*pyz(i,j,k)*v(i,j,k) + 3*pzz(i,j,k)*v(i,j,k) + 3*(pxx(i,j,k) - 3*pxy(i,j,k) - 2*pxz(i,j,k) + pyy(i,j,k) + 2*pyz(i,j,k))*w(i,j,k))
                          f(i-1,j+1,k+1,26)=feq + (1-omega)*fneq1*p3 !+ (feqs-feq)
							
						              !g0
                          feq=p0g*phi(i,j,k)
                          g(i,j,k,0)=feq
                          !g1
                          feqs=p1g*sharp_c*phi(i,j,k)*(1-phi(i,j,k))*(normx(i,j,k))
                          feq=p1g*phi(i,j,k)*(1 + 3*u(i,j,k))
                          g(i+1,j,k,1)=feq+feqs
                          !g2
                          feq=p1g*phi(i,j,k)*(1 - 3*u(i,j,k))
                          g(i-1,j,k,2)=feq-feqs
                          !g3
                          feqs=p1g*sharp_c*phi(i,j,k)*(1-phi(i,j,k))*(normy(i,j,k))
                          feq=p1g*phi(i,j,k)*(1 + 3*v(i,j,k))
                          g(i,j+1,k,3)=feq+feqs
                          !g4
                          feq=p1g*phi(i,j,k)*(1 - 3*v(i,j,k))
                          g(i,j-1,k,4)=feq-feqs
                          !g5
                          feqs=p1g*sharp_c*phi(i,j,k)*(1-phi(i,j,k))*(normz(i,j,k))
                          feq=p1g*phi(i,j,k)*(1 + 3*w(i,j,k))
                          g(i,j,k+1,5)=feq+feqs
                          !g6
                          feq=p1g*phi(i,j,k)*(1 - 3*w(i,j,k))
                          g(i,j,k-1,6)=feq-feqs
                      endif
                  enddo
              enddo
          enddo
          !$acc end kernels
   endsubroutine
   !************************************************************************************************************************************************************************************************************************************!
   subroutine initial_conditions

      implicit none
      integer:: subchords(3)
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
               if((float(gi)-lx/2.0)**2.0 + (float(gj)-ly/2.0)**2.0<=radius**2.0)then
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

   endsubroutine
   !************************************************************************************************************************************************************************************************************************************!
   subroutine initial_conditions_multiphase
   
   endsubroutine
endmodule
