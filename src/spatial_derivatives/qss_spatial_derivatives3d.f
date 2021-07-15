c***********************************************************************
      subroutine qss3dHJENO2testNew(
     &  phi,
     &  lse_rhs,
     &  vel_n,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  dx, dy, dz)
c***********************************************************************
c { begin subroutine
      implicit none
      
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      real phi_x_plus
      real phi_y_plus
      real phi_z_plus
      real phi_x_minus
      real phi_y_minus
      real phi_z_minus
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real vel_n(ilo_phi_gb:ihi_phi_gb,
     &           jlo_phi_gb:jhi_phi_gb,
     &           klo_phi_gb:khi_phi_gb)
      real lse_rhs(ilo_phi_gb:ihi_phi_gb,
     &             jlo_phi_gb:jhi_phi_gb,
     &             klo_phi_gb:khi_phi_gb)
      real D1_x, D1_x_plus, D2_x, D2_x_plus, D2_x_minus
      real D1_y, D1_y_plus, D2_y, D2_y_plus, D2_y_minus
      real D1_z, D1_z_plus, D2_z, D2_z_plus, D2_z_minus
      real temp
      real dx, dy, dz
      real inv_dx, inv_dy, inv_dz
      integer i, j, k
      real half
      parameter (half=0.5d0)
      integer order_1, order_2
      parameter (order_1=1,order_2=2)
      integer x_dir, y_dir, z_dir
      parameter (x_dir=1,y_dir=2,z_dir=3)
      real vel_n_cur
      real norm_grad_phi_sq
      real zero_tol
      parameter (zero_tol=1.e-5)
      

c     compute inv_dx, inv_dy, and inv_dz
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy
      inv_dz = 1.0d0/dz
      
c     { begin loop over grid 
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

c           phi_x_plus
            D1_x = phi(i,j,k) - phi(i-1,j,k)
            D1_x_plus = phi(i+1,j,k) - phi(i,j,k)
            D2_x = phi(i+1,j,k) + phi(i-1,j,k) - 2*phi(i,j,k)
            D2_x_plus = phi(i+2,j,k) + phi(i,j,k) - 2*phi(i+1,j,k)
            D2_x_minus = phi(i-2,j,k) + phi(i,j,k) - 2*phi(i-1,j,k)
            
            if (abs(D2_x).lt.abs(D2_x_plus)) then
              phi_x_plus = (D1_x_plus 
     &                          - half*D2_x)*inv_dx
            else
              phi_x_plus = (D1_x_plus 
     &                          - half*D2_x_plus)*inv_dx
            endif

c           phi_x_minus
            if (abs(D2_x_minus).lt.abs(D2_x)) then
              phi_x_minus = (D1_x 
     &                           + half*D2_x_minus)*inv_dx
            else
              phi_x_minus = (D1_x
     &                           + half*D2_x)*inv_dx
            endif

c           phi_y_plus
            D1_y = phi(i,j,k) - phi(i,j-1,k)
            D1_y_plus = phi(i,j+1,k) - phi(i,j,k)
            D2_y = phi(i,j+1,k) + phi(i,j-1,k) - 2*phi(i,j,k)
            D2_y_plus = phi(i,j+2,k) + phi(i,j,k) - 2*phi(i,j+1,k)
            D2_y_minus = phi(i,j-2,k) + phi(i,j,k) - 2*phi(i,j-1,k)
            
            if (abs(D2_y).lt.abs(D2_y_plus)) then
              phi_y_plus = (D1_y_plus 
     &                          - half*D2_y)*inv_dy
            else
              phi_y_plus = (D1_y_plus 
     &                          - half*D2_y_plus)*inv_dy
            endif

c           phi_y_minus
            if (abs(D2_y_minus).lt.abs(D2_y)) then
              phi_y_minus = (D1_y 
     &                           + half*D2_y_minus)*inv_dy
            else
              phi_y_minus = (D1_y
     &                           + half*D2_y)*inv_dy
            endif
            
c           phi_z_plus
            D1_z = phi(i,j,k) - phi(i,j,k-1)
            D1_z_plus = phi(i,j,k+1) - phi(i,j,k)
            D2_z = phi(i,j,k+1) + phi(i,j,k-1) - 2*phi(i,j,k)
            D2_z_plus = phi(i,j,k+2) + phi(i,j,k) - 2*phi(i,j,k+1)
            D2_z_minus = phi(i,j,k-2) + phi(i,j,k) - 2*phi(i,j,k-1)
            
            if (abs(D2_z).lt.abs(D2_z_plus)) then
              phi_z_plus = (D1_z_plus 
     &                          - half*D2_z)*inv_dz
            else
              phi_z_plus = (D1_z_plus 
     &                          - half*D2_z_plus)*inv_dz
            endif

c           phi_z_minus
            if (abs(D2_z_minus).lt.abs(D2_z)) then
              phi_z_minus = (D1_z 
     &                           + half*D2_z_minus)*inv_dz
            else
              phi_z_minus = (D1_z
     &                           + half*D2_z)*inv_dz
            endif
            
            
            
            vel_n_cur = vel_n(i,j,k)
            if (abs(vel_n_cur) .ge. zero_tol) then

c             { begin Godunov selection of grad_phi

              if (vel_n_cur .gt. 0.d0) then
                norm_grad_phi_sq = max(max(phi_x_minus,0.d0)**2,
     &                                 min(phi_x_plus,0.d0)**2 )
     &                           + max(max(phi_y_minus,0.d0)**2,
     &                                 min(phi_y_plus,0.d0)**2 )
     &                           + max(max(phi_z_minus,0.d0)**2,
     &                                 min(phi_z_plus,0.d0)**2 )
              else
                norm_grad_phi_sq = max(min(phi_x_minus,0.d0)**2,
     &                                 max(phi_x_plus,0.d0)**2 )
     &                           + max(min(phi_y_minus,0.d0)**2,
     &                                 max(phi_y_plus,0.d0)**2 )
     &                           + max(min(phi_z_minus,0.d0)**2,
     &                                 max(phi_z_plus,0.d0)**2 )
              endif

c             } end Godunov selection of grad_phi

c            lse_rhs(i,j,k) = - vel_n_cur
c             compute contribution to lse_rhs(i,j,k) 
              temp = lse_rhs(i,j,k) - vel_n_cur*sqrt(norm_grad_phi_sq)
              lse_rhs(i,j,k) = temp
     
            endif
          enddo
        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm3dSetVarNorm() 
c  
c  Arguments:
c    *_gb (in):         index range for ghostbox
c    *_fb (in):         index range for fillbox
c    index_[xyz](in):  [xyz] coordinates of local (narrow band) points
c    n*_index(in):    index range of points to loop over in index_*
c    narrow_band(in): array that marks voxels outside desired fillbox
c    mark_fb(in):     upper limit narrow band value for voxels in 
c                     fillbox
c
c***********************************************************************
      subroutine qss3dSetVarNorm(
     &  mask,
     &  mask_x,
     &  mask_y,
     &  mask_z,
     &  var_a,
     &  a0,
     &  ca,
     &  ilo_gb, 
     &  ihi_gb,
     &  jlo_gb, 
     &  jhi_gb,
     &  klo_gb, 
     &  khi_gb,
     &  ilo_fb,
     &  ihi_fb,
     &  jlo_fb,
     &  jhi_fb,
     &  klo_fb,
     &  khi_fb,
     &  dx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      integer klo_fb, khi_fb
      real mask(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)    
      real mask_x(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)
      real mask_y(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)
      real mask_z(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)
      real var_a(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)      
      real dx
      real ca
      real a0
      integer i, j, k
      real h1, h2, s, C2, C1, eps_var, absgrad, x
      eps_var = 1.5 * dx      
      C1 = 0.04
      C2 = 1
      
c       loop over included cells {
        do k=klo_fb,khi_fb
          do j=jlo_fb,jhi_fb
            do i=ilo_fb,ihi_fb
                call qssHeaviside(h1, -mask(i,j,k), eps_var)
                call qssHeaviside(h2, mask(i,j,k), eps_var)
                absgrad = mask_x(i,j,k) * mask_x(i,j,k) +
     &                     mask_y(i,j,k) * mask_y(i,j,k) + 
     &                     mask_z(i,j,k) * mask_z(i,j,k)
                absgrad = SQRT(absgrad)
                x = mask(i,j,k)*mask(i,j,k) + absgrad*dx*absgrad*dx
                x = SQRT(x)
                s = mask(i,j,k)/x
                var_a(i,j,k) = h1 * a0 - s * h2 * C2 * COS(ca)*absgrad
            enddo
          enddo
        enddo
c       } end loop over grid
          return
      end
c***********************************************************************

c***********************************************************************
c
c  lsm3dSetVarCurvAdv() 
c  
c  Arguments:
c    *_gb (in):         index range for ghostbox
c    *_fb (in):         index range for fillbox
c    index_[xyz](in):  [xyz] coordinates of local (narrow band) points
c    n*_index(in):    index range of points to loop over in index_*
c    narrow_band(in): array that marks voxels outside desired fillbox
c    mark_fb(in):     upper limit narrow band value for voxels in 
c                     fillbox
c
c***********************************************************************
      subroutine qss3dSetVarCurvAdv(
     &  mask,
     &  mask_x,
     &  mask_y,
     &  mask_z,
     &  var_b,
     &  vel_x,
     &  vel_y,
     &  vel_z,
     &  b_max_over_dx,
     &  max_U_over_dx,
     &  const_b,
     &  ilo_gb, 
     &  ihi_gb,
     &  jlo_gb, 
     &  jhi_gb,
     &  klo_gb, 
     &  khi_gb,
     &  ilo_fb,
     &  ihi_fb,
     &  jlo_fb,
     &  jhi_fb,
     &  klo_fb,
     &  khi_fb,
     &  dx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      integer klo_fb, khi_fb
      real mask(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)      
      real mask_x(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)
      real mask_y(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)
      real mask_z(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)
      real var_b(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)      
      real vel_x(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)
      real vel_y(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb) 
      real vel_z(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb) 
      real dx
      real const_b
      integer i, j, k
      real h1, h2, s, C2, C1, eps_var, x, absgrad
      real b_max_over_dx, max_U_over_dx, b_max, U_over_dX_cur
      real inv_dx, inv_dy, inv_dz

      eps_var = 1.5 * dx      
      C1 = const_b
      C2 = 1
      
c     initialize max_U_over_dX to -1
      max_U_over_dX = -1.0d0
      
c     initialize b_max to -1
      b_max = -1.0d0  

c     compute inv_dx, inv_dy, and inv_dz
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dx
      inv_dz = 1.0d0/dx
      
c       loop over included cells {
        do k=klo_fb,khi_fb
          do j=jlo_fb,jhi_fb
            do i=ilo_fb,ihi_fb
                call qssHeaviside(h1, -mask(i,j,k), eps_var)
                call qssHeaviside(h2, mask(i,j,k), eps_var)
                absgrad = mask_x(i,j,k) * mask_x(i,j,k) +
     &                     mask_y(i,j,k) * mask_y(i,j,k) + 
     &                     mask_z(i,j,k) * mask_z(i,j,k)
                absgrad = SQRT(absgrad)
                x = mask(i,j,k)*mask(i,j,k) + absgrad*dx*absgrad*dx
                x = SQRT(x)
                s = mask(i,j,k)/x
                var_b(i,j,k) = h1*C1
c                h2*C1
                
                if (var_b(i,j,k) > b_max) then
                    b_max = var_b(i,j,k)
                endif
                
                vel_x(i,j,k) = s * C2 * mask_x(i,j,k) * h2
                vel_y(i,j,k) = s * C2 * mask_y(i,j,k) * h2
                vel_z(i,j,k) = s * C2 * mask_z(i,j,k) * h2
                
                U_over_dX_cur = abs(vel_x(i,j,k))*inv_dx
     &                        + abs(vel_y(i,j,k))*inv_dy
     &                        + abs(vel_z(i,j,k))*inv_dz

                if (U_over_dX_cur .gt. max_U_over_dX) then
                  max_U_over_dX = U_over_dX_cur
                endif  
            enddo
          enddo
        enddo
c       } end loop over grid

          b_max_over_dx = 2 * b_max * 
     &               (inv_dx*inv_dx + inv_dy*inv_dy + inv_dz*inv_dz)
          return
      end
c***********************************************************************

c***********************************************************************
      subroutine qss3dGetRHS2nd(
     &  phi,
     &  vel_n,
     &  b,
     &  vel_x, vel_y, vel_z,
     &  max_H_over_dX,
     &  lse_rhs,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  dx, dy, dz)
c***********************************************************************
c { begin subroutine
      implicit none
      
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      
c     plus and minus upwind derivatives
      real phi_x_plus, phi_y_plus, phi_z_plus
      real phi_x_minus, phi_y_minus, phi_z_minus
      
c     upwinded derivatives
      real phi_x_adv, phi_y_adv, phi_z_adv
    
c     central difference derivatives
      real phi_x, phi_y, phi_z
      
c     second order central difference derivatives
      real phi_xx, phi_yy, phi_zz, phi_xy, phi_xz, phi_yz
      
      real curv, grad_mag2
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real vel_n(ilo_phi_gb:ihi_phi_gb,
     &           jlo_phi_gb:jhi_phi_gb,
     &           klo_phi_gb:khi_phi_gb)
      real b(ilo_phi_gb:ihi_phi_gb,
     &           jlo_phi_gb:jhi_phi_gb,
     &           klo_phi_gb:khi_phi_gb)
      real vel_x(ilo_phi_gb:ihi_phi_gb,
     &           jlo_phi_gb:jhi_phi_gb,
     &           klo_phi_gb:khi_phi_gb)
      real vel_y(ilo_phi_gb:ihi_phi_gb,
     &           jlo_phi_gb:jhi_phi_gb,
     &           klo_phi_gb:khi_phi_gb)
      real vel_z(ilo_phi_gb:ihi_phi_gb,
     &           jlo_phi_gb:jhi_phi_gb,
     &           klo_phi_gb:khi_phi_gb)
      real lse_rhs(ilo_phi_gb:ihi_phi_gb,
     &             jlo_phi_gb:jhi_phi_gb,
     &             klo_phi_gb:khi_phi_gb)
      real D1_x, D1_x_plus, D2_x, D2_x_plus, D2_x_minus
      real D1_y, D1_y_plus, D2_y, D2_y_plus, D2_y_minus
      real D1_z, D1_z_plus, D2_z, D2_z_plus, D2_z_minus
      
      real dx, dy, dz
      real inv_dx, inv_dy, inv_dz
      integer i, j, k
      real half
      parameter (half=0.5d0)
      integer order_1, order_2
      parameter (order_1=1,order_2=2)
      integer x_dir, y_dir, z_dir
      parameter (x_dir=1,y_dir=2,z_dir=3)
      real vel_n_cur
      real norm_grad_phi_sq
      real zero_tol
      parameter (zero_tol=1.e-5)
      real max_H_over_dX
      real H_over_dX_cur
      real max_dx_sq
      real norm_grad_phi, phi_x_cur, phi_y_cur, phi_z_cur
      real temp_norm, temp_curv, temp_vel

      real dx_factor, dy_factor, dz_factor

c     compute denominator values, for curvature term
      dx_factor = 0.5d0/dx
      dy_factor = 0.5d0/dy
      dz_factor = 0.5d0/dz

c     compute inv_dx, inv_dy, and inv_dz
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy
      inv_dz = 1.0d0/dz
      
c     { begin loop over grid 
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

c           phi_x_plus
            D1_x = phi(i,j,k) - phi(i-1,j,k)
            D1_x_plus = phi(i+1,j,k) - phi(i,j,k)
            D2_x = phi(i+1,j,k) + phi(i-1,j,k) - 2*phi(i,j,k)
            D2_x_plus = phi(i+2,j,k) + phi(i,j,k) - 2*phi(i+1,j,k)
            D2_x_minus = phi(i-2,j,k) + phi(i,j,k) - 2*phi(i-1,j,k)
            
            if (abs(D2_x).lt.abs(D2_x_plus)) then
              phi_x_plus = (D1_x_plus 
     &                          - half*D2_x)*inv_dx
            else
              phi_x_plus = (D1_x_plus 
     &                          - half*D2_x_plus)*inv_dx
            endif

c           phi_x_minus
            if (abs(D2_x_minus).lt.abs(D2_x)) then
              phi_x_minus = (D1_x 
     &                           + half*D2_x_minus)*inv_dx
            else
              phi_x_minus = (D1_x
     &                           + half*D2_x)*inv_dx
            endif

c           phi_y_plus
            D1_y = phi(i,j,k) - phi(i,j-1,k)
            D1_y_plus = phi(i,j+1,k) - phi(i,j,k)
            D2_y = phi(i,j+1,k) + phi(i,j-1,k) - 2*phi(i,j,k)
            D2_y_plus = phi(i,j+2,k) + phi(i,j,k) - 2*phi(i,j+1,k)
            D2_y_minus = phi(i,j-2,k) + phi(i,j,k) - 2*phi(i,j-1,k)
            
            if (abs(D2_y).lt.abs(D2_y_plus)) then
              phi_y_plus = (D1_y_plus 
     &                          - half*D2_y)*inv_dy
            else
              phi_y_plus = (D1_y_plus 
     &                          - half*D2_y_plus)*inv_dy
            endif

c           phi_y_minus
            if (abs(D2_y_minus).lt.abs(D2_y)) then
              phi_y_minus = (D1_y 
     &                           + half*D2_y_minus)*inv_dy
            else
              phi_y_minus = (D1_y
     &                           + half*D2_y)*inv_dy
            endif
            
c           phi_z_plus
            D1_z = phi(i,j,k) - phi(i,j,k-1)
            D1_z_plus = phi(i,j,k+1) - phi(i,j,k)
            D2_z = phi(i,j,k+1) + phi(i,j,k-1) - 2*phi(i,j,k)
            D2_z_plus = phi(i,j,k+2) + phi(i,j,k) - 2*phi(i,j,k+1)
            D2_z_minus = phi(i,j,k-2) + phi(i,j,k) - 2*phi(i,j,k-1)
            
            if (abs(D2_z).lt.abs(D2_z_plus)) then
              phi_z_plus = (D1_z_plus 
     &                          - half*D2_z)*inv_dz
            else
              phi_z_plus = (D1_z_plus 
     &                          - half*D2_z_plus)*inv_dz
            endif

c           phi_z_minus
            if (abs(D2_z_minus).lt.abs(D2_z)) then
              phi_z_minus = (D1_z 
     &                           + half*D2_z_minus)*inv_dz
            else
              phi_z_minus = (D1_z
     &                           + half*D2_z)*inv_dz
            endif
            
            
            vel_n_cur = vel_n(i,j,k)
            if (abs(vel_n_cur) .ge. zero_tol) then

c             { begin Godunov selection of grad_phi

              if (vel_n_cur .gt. 0.d0) then
                norm_grad_phi_sq = max(max(phi_x_minus,0.d0)**2,
     &                                 min(phi_x_plus,0.d0)**2 )
     &                           + max(max(phi_y_minus,0.d0)**2,
     &                                 min(phi_y_plus,0.d0)**2 )
     &                           + max(max(phi_z_minus,0.d0)**2,
     &                                 min(phi_z_plus,0.d0)**2 )
              else
                norm_grad_phi_sq = max(min(phi_x_minus,0.d0)**2,
     &                                 max(phi_x_plus,0.d0)**2 )
     &                           + max(min(phi_y_minus,0.d0)**2,
     &                                 max(phi_y_plus,0.d0)**2 )
     &                           + max(min(phi_z_minus,0.d0)**2,
     &                                 max(phi_z_plus,0.d0)**2 )
              endif

c             } end Godunov selection of grad_phi


c             compute contribution to lse_rhs(i,j,k) 
              temp_norm = - vel_n_cur*sqrt(norm_grad_phi_sq)

            endif
        
c           get max_H_over_dX for calculating stable dt
c           initialize max_H_over_dX to -1
            max_H_over_dX = -1.0d0

c           compute max_dx_sq
            max_dx_sq = max(dx,dy,dz)
            max_dx_sq = max(dx,dy,dz) * max(dx,dy,dz)
      
            phi_x_cur = max(abs(phi_x_plus),
     &                          abs(phi_x_minus))
            phi_y_cur = max(abs(phi_y_plus),
     &                          abs(phi_y_minus))
            phi_z_cur = max(abs(phi_z_plus),
     &                       abs(phi_z_minus))
            norm_grad_phi = sqrt( phi_x_cur*phi_x_cur 
     &                            + phi_y_cur*phi_y_cur 
     &                            + phi_z_cur*phi_z_cur + max_dx_sq )

            H_over_dX_cur = abs(vel_n(i,j,k)) / norm_grad_phi
     &                        * ( phi_x_cur*inv_dx 
     &                        + phi_y_cur*inv_dy 
     &                        + phi_z_cur*inv_dz )

            if (H_over_dX_cur .gt. max_H_over_dX) then
                  max_H_over_dX = H_over_dX_cur  
            endif
        
            
    
c           compute advective term contributions
c              phi_x
              if (vel_x(i,j,k).ge.0) then
                phi_x_adv = phi_x_minus
              else
                phi_x_adv = phi_x_plus
              endif

c              phi_y
              if (vel_y(i,j,k).ge.0) then
                phi_y_adv = phi_y_minus
              else
                phi_y_adv = phi_y_plus
              endif    
              
c              phi_z
              if (vel_z(i,j,k).ge.0) then
                phi_z_adv = phi_z_minus
              else
                phi_z_adv= phi_z_plus
              endif   
              
            temp_vel = - ( vel_x(i,j,k)*phi_x_adv
     &                       + vel_y(i,j,k)*phi_y_adv 
     &                       + vel_z(i,j,k)*phi_z_adv )
     
     
c         compute curvature term
            phi_x = (phi(i+1,j,k) - phi(i-1,j,k))*dx_factor
            phi_y = (phi(i,j+1,k) - phi(i,j-1,k))*dy_factor
            phi_z = (phi(i,j,k+1) - phi(i,j,k-1))*dz_factor
   
            phi_xx = (phi(i+2,j,k) + phi(i-2,j,k)
     &               - 2*phi(i,j,k))*dx_factor*dx_factor
            phi_zz = (phi(i,j,k+2) + phi(i,j,k-2)
     &               - 2*phi(i,j,k))*dz_factor*dz_factor
            phi_yy = (phi(i,j+2,k) + phi(i,j-2,k)
     &               - 2*phi(i,j,k))*dy_factor*dy_factor
        
            phi_xy = (phi(i+1,j+1,k) + phi(i-1,j-1,k)
     &               - phi(i-1,j+1,k) - phi(i+1,j-1,k))
     &               *dx_factor*dy_factor
            phi_xz = (phi(i+1,j,k+1) + phi(i-1,j,k-1)
     &               - phi(i-1,j,k+1) - phi(i+1,j,k-1))
     &               *dx_factor*dz_factor
            phi_yz = (phi(i,j+1,k+1) + phi(i,j-1,k-1)
     &               - phi(i,j-1,k+1) - phi(i,j+1,k-1))
     &               *dy_factor*dz_factor
    
c           compute squared magnitude of gradient
            grad_mag2 = phi_x * phi_x
     &                + phi_y * phi_y
     &                + phi_z * phi_z
            if (grad_mag2 .lt. zero_tol) then
              curv = 0.d0
            else
              curv = phi_xx*phi_y*phi_y 
     &             +   phi_yy*phi_x*phi_x  
     &             - 2*phi_xy*phi_x*phi_y
     &             +   phi_xx*phi_z*phi_z 
     &             +   phi_zz*phi_x*phi_x 
     &             - 2*phi_xz*phi_x*phi_z
     &             +   phi_yy*phi_z*phi_z  
     &             +   phi_zz*phi_y*phi_y 
     &             - 2*phi_yz*phi_y*phi_z
              curv = curv / grad_mag2 
            endif

            temp_curv = b(i,j,k)*curv

	    lse_rhs(i,j,k) = temp_norm + temp_curv + temp_vel    
          enddo
        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************
c***********************************************************************
c
c  lsm3dCentralGradOrder2() computes the second-order, central,
c  finite difference approximation to the gradient of phi.
c
c  Arguments:
c    phi_* (out):      components of grad(phi) 
c    phi (in):         phi
c    dx, dy, dz (in):  grid spacing
c    *_gb (in):        index range for ghostbox
c    *_fb (in):        index range for fillbox
c
c***********************************************************************
      subroutine qss3dCentralGradOrder2(
     &  phi_x, phi_y, phi_z,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb, 
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  klo_grad_phi_gb, khi_grad_phi_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  dx, dy, dz)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_gb refers to ghostbox for grad_phi data
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer klo_grad_phi_gb, khi_grad_phi_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb,
     &           klo_grad_phi_gb:khi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb,
     &           klo_grad_phi_gb:khi_grad_phi_gb)
      real phi_z(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb,
     &           klo_grad_phi_gb:khi_grad_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real dx, dy, dz
      integer i,j,k
      real dx_factor, dy_factor, dz_factor

c     compute denominator values
      dx_factor = 0.5d0/dx
      dy_factor = 0.5d0/dy
      dz_factor = 0.5d0/dz

c     { begin loop over grid 
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

            phi_x(i,j,k) = (phi(i+1,j,k) - phi(i-1,j,k))*dx_factor
            phi_y(i,j,k) = (phi(i,j+1,k) - phi(i,j-1,k))*dy_factor
            phi_z(i,j,k) = (phi(i,j,k+1) - phi(i,j,k-1))*dz_factor
   
          enddo
        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
      subroutine qss3dSignedLinearExtrapolation(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      
c     local variables       
      integer i,j,k
      real s, abs_diff, dist, slope
      real one
      parameter (one = 1.0d0)

      if (bdry_location_idx .eq. 0) then
c     { extrapolate data in x-direction at lower end

c       { begin k,j loop
        do k = klo_gb, khi_gb
          do j = jlo_gb, jhi_gb
            s = sign(one,phi(ilo_fb,j,k))
            abs_diff = abs(phi(ilo_fb,j,k) - phi(ilo_fb+1,j,k))
            slope = s*abs_diff
            do i = ilo_gb, ilo_fb-1
              dist = ilo_fb - i
              phi(i,j,k) = phi(ilo_fb,j,k) + slope*dist
            enddo
          enddo
        enddo
c       } end k,j loop

c     } end extrapolate data in x-direction at lower end

      elseif (bdry_location_idx .eq. 1) then
c     { extrapolate data in x-direction at upper end

c       { begin k,j loop
        do k = klo_gb, khi_gb
          do j = jlo_gb, jhi_gb
            s = sign(one,phi(ihi_fb,j,k))
            abs_diff = abs(phi(ihi_fb,j,k) - phi(ihi_fb-1,j,k))
            slope = s*abs_diff
            do i = ihi_fb+1, ihi_gb
              dist = i - ihi_fb
              phi(i,j,k) = phi(ihi_fb,j,k) + slope*dist
            enddo 
          enddo
        enddo
c       } end k,j loop

c     } end extrapolate data in x-direction at upper end

      elseif (bdry_location_idx .eq. 2) then
c     { extrapolate data in y-direction at lower end

c       { begin i,k loop
        do i = ilo_gb, ihi_gb
          do k = klo_gb, khi_gb
            s = sign(one,phi(i,jlo_fb,k))
            abs_diff = abs(phi(i,jlo_fb,k) - phi(i,jlo_fb+1,k))
            slope = s*abs_diff
            do j = jlo_gb, jlo_fb-1
              dist = jlo_fb - j
              phi(i,j,k) = phi(i,jlo_fb,k) + slope*dist
            enddo
          enddo
        enddo
c       } end i,k loop

c     } end extrapolate data in y-direction at lower end

      elseif (bdry_location_idx .eq. 3) then
c     { extrapolate data in y-direction at upper end

c       { begin i,k loop
        do i = ilo_gb, ihi_gb
          do k = klo_gb, khi_gb
            s = sign(one,phi(i,jhi_fb,k))
            abs_diff = abs(phi(i,jhi_fb,k) - phi(i,jhi_fb-1,k))
            slope = s*abs_diff
            do j = jhi_fb+1, jhi_gb
              dist = j - jhi_fb
              phi(i,j,k) = phi(i,jhi_fb,k) + slope*dist
            enddo 
          enddo
        enddo
c       } end i,k loop

c     } end extrapolate data in y-direction at upper end

      elseif (bdry_location_idx .eq. 4) then
c     { extrapolate data in z-direction at lower end

c       { begin i,j loop
        do i = ilo_gb, ihi_gb
          do j = jlo_gb, jhi_gb
            s = sign(one,phi(i,j,klo_fb))
            abs_diff = abs(phi(i,j,klo_fb) - phi(i,j,klo_fb+1))
            slope = s*abs_diff
            do k = klo_gb, klo_fb-1
              dist = klo_fb - k
              phi(i,j,k) = phi(i,j,klo_fb) + slope*dist
            enddo
          enddo
        enddo
c       } end i,j loop

c     } end extrapolate data in z-direction at lower end

      elseif (bdry_location_idx .eq. 5) then
c     { extrapolate data in z-direction at upper end

c       { begin i,j loop
        do i = ilo_gb, ihi_gb
          do j = jlo_gb, jhi_gb
            s = sign(one,phi(i,j,khi_fb))
            abs_diff = abs(phi(i,j,khi_fb) - phi(i,j,khi_fb-1))
            slope = s*abs_diff
            do k = khi_fb+1, khi_gb
              dist = k - khi_fb
              phi(i,j,k) = phi(i,j,khi_fb) + slope*dist
            enddo 
          enddo
        enddo
c       } end i,j loop

c     } end extrapolate data in z-direction at upper end

      endif
      return
      end
c } end subroutine
c***********************************************************************
