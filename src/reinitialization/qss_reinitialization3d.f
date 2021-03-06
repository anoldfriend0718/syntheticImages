c***********************************************************************
c
c  File:        lsm_reinitialization3d.f
c***********************************************************************
c The algorithms and notation in these subroutines closely follows
c the discussion in Osher & Fedkiw (2003).
c***********************************************************************

c***********************************************************************
c
c  lsm3dComputeReinitializationEqnRHS() computes the right-hand side of 
c  the reinitialization equation using a Godunov scheme to select the 
c  numerical discretization of the sgn(phi) |grad(phi)| term.  
c  Forward (plus) and backward (minus) spatial derivatives used in 
c  the Godunov calculation must be supplied by the user.
c
c  Arguments:
c    reinit_rhs (out):       right-hand side of reinitialization 
c                            equation
c    phi (in):               level set function at current iteration
c                            of reinitialization process
c    phi0 (in):              level set function at initial iteration
c                            iteration of reinitialization process
c    phi_*_plus (in):        forward spatial derivatives for grad(phi)
c    phi_*_minus (in):       backward spatial derivatives for grad(phi)
c    use_phi0_for_sgn (in):  flag to specify whether phi0 should be
c                            used in the computation of sgn(phi).
c                              0 = use phi (do NOT use phi0)
c                              1 = use phi0
c    *_gb (in):              index range for ghostbox
c    *_fb (in):              index range for fillbox
c
c  NOTES:
c   (1) if use_phi0_for_sgn is not equal to 0 or 1, the default
c       behavior of lsm3dComputeReinitializationEqnRHS() is to use 
c       phi (i.e. equivalent to setting use_phi0_for_sgn to 0)
c
c***********************************************************************
      subroutine qss3dComputeReinitializationEqnRHS(
     &  reinit_rhs,
     &  ilo_rhs_gb, ihi_rhs_gb, 
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  klo_rhs_gb, khi_rhs_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  phi0,
     &  ilo_phi0_gb, ihi_phi0_gb, 
     &  jlo_phi0_gb, jhi_phi0_gb,
     &  klo_phi0_gb, khi_phi0_gb,
     &  phi_x_plus, phi_y_plus, phi_z_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  klo_grad_phi_plus_gb, khi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus, phi_z_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  klo_grad_phi_minus_gb, khi_grad_phi_minus_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  dx, dy, dz,
     &  use_phi0_for_sgn)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer klo_rhs_gb, khi_rhs_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_phi0_gb, ihi_phi0_gb
      integer jlo_phi0_gb, jhi_phi0_gb
      integer klo_phi0_gb, khi_phi0_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer klo_grad_phi_plus_gb, khi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer klo_grad_phi_minus_gb, khi_grad_phi_minus_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      real reinit_rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                jlo_rhs_gb:jhi_rhs_gb,
     &                klo_rhs_gb:khi_rhs_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real phi0(ilo_phi0_gb:ihi_phi0_gb,
     &          jlo_phi0_gb:jhi_phi0_gb,
     &          klo_phi0_gb:khi_phi0_gb)
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_z_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real phi_z_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real dx, dy, dz
      integer use_phi0_for_sgn
      real phi_cur
      integer DIM
      parameter (DIM=3)
      real grad_phi_plus_cur(1:DIM)
      real grad_phi_minus_cur(1:DIM)
      real grad_phi_star(1:DIM)
      integer i,j,k
      integer dir
      real sgn_phi
      real norm_grad_phi_sq
      real dx_sq
      real zero_tol
      parameter (zero_tol=1.e-5)
      real one
      parameter (one=1.d0)

c     set value of dx_sq to be square of max{dx,dy,dz}
      dx_sq = max(dx,dy,dz)
      dx_sq = dx_sq*dx_sq

c----------------------------------------------------
c      compute RHS of reinitialization equation
c      using Godunov's method
c----------------------------------------------------
c     { begin condition on use_phi0_for_sgn
      if (use_phi0_for_sgn .ne. 1) then

c       -----------------------------------------------
c       use phi in computation of smoothed sgn function
c       -----------------------------------------------
c       { begin loop over grid
        do k=klo_fb,khi_fb
          do j=jlo_fb,jhi_fb
            do i=ilo_fb,ihi_fb

c             cache phi and spatial derivative approximations
              phi_cur = phi(i,j,k)
              grad_phi_plus_cur(1) = phi_x_plus(i,j,k)
              grad_phi_plus_cur(2) = phi_y_plus(i,j,k)
              grad_phi_plus_cur(3) = phi_z_plus(i,j,k)
              grad_phi_minus_cur(1) = phi_x_minus(i,j,k)
              grad_phi_minus_cur(2) = phi_y_minus(i,j,k)
              grad_phi_minus_cur(3) = phi_z_minus(i,j,k)
  
c             { begin Godunov selection of grad_phi
              do dir=1,DIM

                if (phi_cur .gt. 0.d0) then
                  grad_phi_plus_cur(dir) = max(-grad_phi_plus_cur(dir),
     &                                         0.d0)
                  grad_phi_minus_cur(dir) = max(grad_phi_minus_cur(dir),
     &                                          0.d0)
                else
                  grad_phi_plus_cur(dir) = 
     &              max(grad_phi_plus_cur(dir), 0.d0)
                  grad_phi_minus_cur(dir) = 
     &              max(-grad_phi_minus_cur(dir), 0.d0)
                endif

                grad_phi_star(dir) = max(grad_phi_plus_cur(dir),
     &                                   grad_phi_minus_cur(dir)) 

              enddo
c             } end Godunov selection of grad_phi

c             compute reinit_rhs(i,j) using smoothed sgn(phi)
              if (abs(phi_cur) .ge. zero_tol) then
                norm_grad_phi_sq = grad_phi_star(1)*grad_phi_star(1)
     &                           + grad_phi_star(2)*grad_phi_star(2)
     &                           + grad_phi_star(3)*grad_phi_star(3)
                sgn_phi = phi_cur
     &                  / sqrt(phi_cur*phi_cur + norm_grad_phi_sq*dx_sq)
                reinit_rhs(i,j,k) = 
     &            sgn_phi*(one - sqrt(norm_grad_phi_sq))
              else
                reinit_rhs(i,j,k) = 0.d0
              endif

            enddo
          enddo
        enddo
c       } end loop over grid

      else

c       ------------------------------------------------
c       use phi0 in computation of smoothed sgn function 
c       ------------------------------------------------
c       { begin loop over grid
        do k=klo_fb,khi_fb
          do j=jlo_fb,jhi_fb
            do i=ilo_fb,ihi_fb
  
c             cache phi and spatial derivative approximations
              phi_cur = phi0(i,j,k)
              grad_phi_plus_cur(1) = phi_x_plus(i,j,k)
              grad_phi_plus_cur(2) = phi_y_plus(i,j,k)
              grad_phi_plus_cur(3) = phi_z_plus(i,j,k)
              grad_phi_minus_cur(1) = phi_x_minus(i,j,k)
              grad_phi_minus_cur(2) = phi_y_minus(i,j,k)
              grad_phi_minus_cur(3) = phi_z_minus(i,j,k)
  
c             { begin Godunov selection of grad_phi
              do dir=1,DIM

                if (phi_cur .gt. 0.d0) then
                  grad_phi_plus_cur(dir) = max(-grad_phi_plus_cur(dir),
     &                                         0.d0)
                  grad_phi_minus_cur(dir) = max(grad_phi_minus_cur(dir),
     &                                          0.d0)
                else
                  grad_phi_plus_cur(dir) = 
     &              max(grad_phi_plus_cur(dir), 0.d0)
                  grad_phi_minus_cur(dir) = 
     &              max(-grad_phi_minus_cur(dir), 0.d0)
                endif

                grad_phi_star(dir) = max(grad_phi_plus_cur(dir),
     &                                   grad_phi_minus_cur(dir)) 

              enddo
c             } end Godunov selection of grad_phi

c             compute reinit_rhs(i,j) using smoothed sgn(phi)
              if (abs(phi_cur) .ge. zero_tol) then
                norm_grad_phi_sq = grad_phi_star(1)*grad_phi_star(1)
     &                           + grad_phi_star(2)*grad_phi_star(2)
     &                           + grad_phi_star(3)*grad_phi_star(3)
                sgn_phi = phi_cur / sqrt(phi_cur*phi_cur + dx_sq)
                reinit_rhs(i,j,k) = 
     &            sgn_phi*(one - sqrt(norm_grad_phi_sq))
              else
                reinit_rhs(i,j,k) = 0.d0
              endif

            enddo
          enddo
        enddo
c       } end loop over grid

      endif
c     } end condition on use_phi0_for_sgn

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm3dComputeOrthogonalizationEqnRHS() computes the right-hand side 
c  of the orthogonalization equation:
c
c    phi_t + grad(phi) dot { sgn(psi)/|grad(psi)| grad(psi) } = 0
c
c  Upwinding is used to select whether the forward (plus) or backward
c  (minus) spatial derivative should be used for grad(phi).  grad(psi)
c  is computed by averaging the forward and backward spatial derivatives
c  for grad(psi).  Forward and backward spatial derivatives used in the 
c  calculation must be supplied by the user.
c
c  Arguments:
c    othro_rhs (out):        right-hand side of orthogonalization
c                            equation
c    psi (in):               data array for psi
c    phi_*_plus (in):        forward spatial derivatives for grad(phi)
c    phi_*_minus (in):       backward spatial derivatives for grad(phi)
c    psi_*_plus (in):        forward spatial derivatives for grad(psi)
c    psi_*_minus (in):       backward spatial derivatives for grad(psi)
c    *_gb (in):              index range for ghostbox
c    *_fb (in):              index range for fillbox
c
c***********************************************************************
      subroutine qss3dComputeOrthogonalizationEqnRHS(
     &  ortho_rhs,
     &  ilo_rhs_gb, ihi_rhs_gb, 
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  klo_rhs_gb, khi_rhs_gb,
     &  phi_x_plus, phi_y_plus, phi_z_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  klo_grad_phi_plus_gb, khi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus, phi_z_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  klo_grad_phi_minus_gb, khi_grad_phi_minus_gb,
     &  psi,
     &  ilo_psi_gb, ihi_psi_gb, 
     &  jlo_psi_gb, jhi_psi_gb,
     &  klo_psi_gb, khi_psi_gb,
     &  psi_x_plus, psi_y_plus, psi_z_plus,
     &  ilo_grad_psi_plus_gb, ihi_grad_psi_plus_gb, 
     &  jlo_grad_psi_plus_gb, jhi_grad_psi_plus_gb,
     &  klo_grad_psi_plus_gb, khi_grad_psi_plus_gb,
     &  psi_x_minus, psi_y_minus, psi_z_minus,
     &  ilo_grad_psi_minus_gb, ihi_grad_psi_minus_gb, 
     &  jlo_grad_psi_minus_gb, jhi_grad_psi_minus_gb,
     &  klo_grad_psi_minus_gb, khi_grad_psi_minus_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  dx, dy, dz)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer klo_rhs_gb, khi_rhs_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer klo_grad_phi_plus_gb, khi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer klo_grad_phi_minus_gb, khi_grad_phi_minus_gb
      integer ilo_psi_gb, ihi_psi_gb
      integer jlo_psi_gb, jhi_psi_gb
      integer klo_psi_gb, khi_psi_gb
      integer ilo_grad_psi_plus_gb, ihi_grad_psi_plus_gb
      integer jlo_grad_psi_plus_gb, jhi_grad_psi_plus_gb
      integer klo_grad_psi_plus_gb, khi_grad_psi_plus_gb
      integer ilo_grad_psi_minus_gb, ihi_grad_psi_minus_gb
      integer jlo_grad_psi_minus_gb, jhi_grad_psi_minus_gb
      integer klo_grad_psi_minus_gb, khi_grad_psi_minus_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      real ortho_rhs(ilo_rhs_gb:ihi_rhs_gb,
     &               jlo_rhs_gb:jhi_rhs_gb,
     &               klo_rhs_gb:khi_rhs_gb)
      real psi(ilo_psi_gb:ihi_psi_gb,
     &         jlo_psi_gb:jhi_psi_gb,
     &         klo_psi_gb:khi_psi_gb)
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_z_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real phi_z_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real psi_x_plus(ilo_grad_psi_plus_gb:ihi_grad_psi_plus_gb,
     &                jlo_grad_psi_plus_gb:jhi_grad_psi_plus_gb,
     &                klo_grad_psi_plus_gb:khi_grad_psi_plus_gb)
      real psi_y_plus(ilo_grad_psi_plus_gb:ihi_grad_psi_plus_gb,
     &                jlo_grad_psi_plus_gb:jhi_grad_psi_plus_gb,
     &                klo_grad_psi_plus_gb:khi_grad_psi_plus_gb)
      real psi_z_plus(ilo_grad_psi_plus_gb:ihi_grad_psi_plus_gb,
     &                jlo_grad_psi_plus_gb:jhi_grad_psi_plus_gb,
     &                klo_grad_psi_plus_gb:khi_grad_psi_plus_gb)
      real psi_x_minus(ilo_grad_psi_minus_gb:ihi_grad_psi_minus_gb,
     &                 jlo_grad_psi_minus_gb:jhi_grad_psi_minus_gb,
     &                 klo_grad_psi_minus_gb:khi_grad_psi_minus_gb)
      real psi_y_minus(ilo_grad_psi_minus_gb:ihi_grad_psi_minus_gb,
     &                 jlo_grad_psi_minus_gb:jhi_grad_psi_minus_gb,
     &                 klo_grad_psi_minus_gb:khi_grad_psi_minus_gb)
      real psi_z_minus(ilo_grad_psi_minus_gb:ihi_grad_psi_minus_gb,
     &                 jlo_grad_psi_minus_gb:jhi_grad_psi_minus_gb,
     &                 klo_grad_psi_minus_gb:khi_grad_psi_minus_gb)
      real dx, dy, dz
      real dx_sq
      integer DIM
      parameter (DIM=3)
      real grad_psi_star(1:DIM)
      real norm_grad_psi
      real sgn_psi
      integer i,j,k
      real psi_zero_tol, grad_psi_zero_tol
      parameter (psi_zero_tol=1.e-5)
      parameter (grad_psi_zero_tol=1.e-5)
      real one, half
      parameter (one=1.d0,half=0.5d0)

c     set value of dx_sq to be square of max{dx,dy,dz}
      dx_sq = max(dx,dy,dz)
      dx_sq = dx_sq*dx_sq

c----------------------------------------------------
c      compute RHS of orthogonalization equation
c      using upwinding
c----------------------------------------------------

c     { begin loop over grid
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

c           average forward and backward derivatives of 
c           grad(psi) to get grad_psi_star
            grad_psi_star(1) = half 
     &                       * (psi_x_plus(i,j,k)+ psi_x_minus(i,j,k))
            grad_psi_star(2) = half 
     &                       * (psi_y_plus(i,j,k)+ psi_y_minus(i,j,k))
            grad_psi_star(3) = half 
     &                       * (psi_z_plus(i,j,k)+ psi_z_minus(i,j,k))
 
c           compute norm of grad(psi)
            norm_grad_psi = grad_psi_star(1)*grad_psi_star(1)
     &                    + grad_psi_star(2)*grad_psi_star(2)
     &                    + grad_psi_star(3)*grad_psi_star(3)
            norm_grad_psi = sqrt(norm_grad_psi)

c           compute ortho_rhs(i,j,k) using upwinding on sgn_psi*grad(psi) 
            if ( (abs(psi(i,j,k)) .gt. psi_zero_tol) .and.
     &           (norm_grad_psi .gt. grad_psi_zero_tol) ) then

c             CASE: nontrivial psi and grad(psi) 

              sgn_psi = psi(i,j,k)/sqrt( psi(i,j,k)*psi(i,j,k) 
     &                                 + dx_sq*norm_grad_psi**2)

              ortho_rhs(i,j,k) = -1.d0/norm_grad_psi 
     &          * ( max(sgn_psi*grad_psi_star(1),0.d0)
     &                         *phi_x_minus(i,j,k)
     &            + min(sgn_psi*grad_psi_star(1),0.d0)
     &                          *phi_x_plus(i,j,k)
     &            + max(sgn_psi*grad_psi_star(2),0.d0)
     &                          *phi_y_minus(i,j,k)
     &            + min(sgn_psi*grad_psi_star(2),0.d0)
     &                          *phi_y_plus(i,j,k)
     &            + max(sgn_psi*grad_psi_star(3),0.d0)
     &                          *phi_z_minus(i,j,k)
     &            + min(sgn_psi*grad_psi_star(3),0.d0)
     &                          *phi_z_plus(i,j,k) )
            else

c             CASE: grad(psi) = 0 or psi = 0 

              ortho_rhs(i,j,k) = 0.d0

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
      subroutine qss3dComputeDistanceForSubcellFix(
     &  distance0,
     &  phi0,
     &  ilo_phi0_gb, ihi_phi0_gb, jlo_phi0_gb, jhi_phi0_gb,
     &  klo_phi0_gb, khi_phi0_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  dx, dy, dz)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_phi0_gb, ihi_phi0_gb, jlo_phi0_gb, jhi_phi0_gb
      integer klo_phi0_gb, khi_phi0_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      real distance0(ilo_phi0_gb:ihi_phi0_gb,
     &               jlo_phi0_gb:jhi_phi0_gb,
     &               klo_phi0_gb:khi_phi0_gb)
      real phi0(ilo_phi0_gb:ihi_phi0_gb,
     &          jlo_phi0_gb:jhi_phi0_gb,
     &          klo_phi0_gb:khi_phi0_gb)
      real dx, dy, dz, max_dx
      real zero_tol, zero
      parameter (zero_tol=1.e-5)
      parameter (zero=0.d0)
      real large_distance_flag
      real near_plus_x, near_minus_x
      real near_plus_y, near_minus_y
      real near_plus_z, near_minus_z
      logical near
      real d1, d2, d3, delta
      integer i,j,k
      
      max_dx = max(dx,dy,dz);
      large_distance_flag = -1000.d0*max_dx;
      
c       { begin loop over grid
        do k=klo_fb,khi_fb
         do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

            near = ( phi0(i+1,j,k) .gt. zero .and.
     &               phi0(i,j,k)   .le. zero ) .or. 
     &             ( phi0(i+1,j,k) .le. zero .and.
     &	             phi0(i,j,k)   .gt. zero ) .or.
     &             ( phi0(i-1,j,k) .gt. zero .and.
     &               phi0(i,j,k)   .le. zero ) .or. 
     &             ( phi0(i-1,j,k) .le. zero .and.
     &	             phi0(i,j,k)   .gt. zero ) .or.
     	     
     &             ( phi0(i,j+1,k) .gt. zero .and.
     &               phi0(i,j,k)   .le. zero ) .or. 
     &             ( phi0(i,j+1,k) .le. zero .and.
     &	             phi0(i,j,k)   .gt. zero ) .or.
     &             ( phi0(i,j-1,k) .gt. zero .and.
     &               phi0(i,j,k)   .le. zero ) .or. 
     &             ( phi0(i,j-1,k) .le. zero .and.
     &	             phi0(i,j,k)   .gt. zero ) .or.	  	     
     
     &             ( phi0(i,j,k+1) .gt. zero .and.
     &               phi0(i,j,k)   .le. zero ) .or. 
     &             ( phi0(i,j,k+1) .le. zero .and.
     &	             phi0(i,j,k)   .gt. zero ) .or.
     &             ( phi0(i,j,k-1) .gt. zero .and.
     &               phi0(i,j,k)   .le. zero ) .or. 
     &             ( phi0(i,j,k-1) .le. zero .and.
     &	             phi0(i,j,k)   .gt. zero );
     	 	     
	   if( near .eqv. .false. ) then
	     distance0(i,j,k) = large_distance_flag;
	   else
c version 1
c	      d1 = phi0(i+1,j,k) - phi0(i-1,j,k);
c	      d2 = phi0(i,j+1,k) - phi0(i,j-1,k);
c	      d3 = phi0(i,j,k+1) - phi0(i,j,k-1);
c	      delta  = sqrt( d1*d1 + d2*d2 + d3*d3 +zero_tol);

c version 2 - (Russo/Smereka)
              d1 = abs(phi0(i+1,j,k) - phi0(i-1,j,k))/2.d0;
	      d2 = abs(phi0(i+1,j,k) - phi0(i,j,k));
	      d3 = abs(phi0(i,j,k)   - phi0(i-1,j,k));
	      delta = max(d1,d2,d3,zero_tol);
	      
	      d1 = abs(phi0(i,j+1,k) - phi0(i,j-1,k))/2.d0;
	      d2 = abs(phi0(i,j+1,k) - phi0(i,j,k));
	      d3 = abs(phi0(i,j,k)   - phi0(i,j-1,k));
	      delta = max(d1,d2,d3,delta);
	      
	      d1 = abs(phi0(i,j,k+1) - phi0(i,j,k-1))/2.d0;
	      d2 = abs(phi0(i,j,k+1) - phi0(i,j,k));
	      d3 = abs(phi0(i,j,k)   - phi0(i,j,k-1));
	      delta = max(d1,d2,d3,delta);
	      
              distance0(i,j,k) = phi0(i,j,k)/delta;     
	   endif
	   
	  enddo 
	 enddo
        enddo
c       } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine qss3dComputeReinitializationEqnRHSSubcellFixOrder1(
     &  reinit_rhs,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  phi0, distance0,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,klo_fb, khi_fb,
     &  dx, dy, dz)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      real reinit_rhs(ilo_phi_gb:ihi_phi_gb,
     &                jlo_phi_gb:jhi_phi_gb,
     &                klo_phi_gb:khi_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real phi0(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real distance0(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real phi_x_plus
      real phi_y_plus
      real phi_z_plus
      real phi_x_minus
      real phi_y_minus
      real phi_z_minus
      real dx, dy, dz, max_dx
      real phi_cur
      integer DIM
      parameter (DIM=3)
      real grad_phi_plus_cur(1:DIM)
      real grad_phi_minus_cur(1:DIM)
      real grad_phi_star(1:DIM)
      integer i,j,k
      integer dir
      real sgn_phi0
      real norm_grad_phi_sq
      real dx_sq
      real zero_tol
      parameter (zero_tol=1.e-5)
      real one, zero, large_distance_flag
      parameter (one=1.d0,  zero=0.d0)

      real D1_x, D1_x_plus, D2_x, D2_x_plus, D2_x_minus
      real D1_y, D1_y_plus, D2_y, D2_y_plus, D2_y_minus
      real D1_z, D1_z_plus, D2_z, D2_z_plus, D2_z_minus
      
      real inv_dx, inv_dy, inv_dz
      
c     set value of dx_sq to be square of max{dx,dy}
      max_dx = max(dx,dy,dz)
      large_distance_flag = -1000.d0*max_dx;
      
c     compute inv_dx, inv_dy, and inv_dz
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy
      inv_dz = 1.0d0/dz
      
      dx_sq = max_dx*max_dx

c----------------------------------------------------
c      compute RHS of reinitialization equation
c      using Godunov's method with subcell fix near interface
c----------------------------------------------------

c       { begin loop over grid
      do k=klo_fb,khi_fb
       do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c           phi_x_plus, minus
            D1_x = phi(i,j,k) - phi(i-1,j,k)
            D1_x_plus = phi(i+1,j,k) - phi(i,j,k)
            phi_x_plus = (D1_x_plus)*inv_dx
            phi_x_minus = (D1_x)*inv_dx            
            
c           phi_y_plus, minus
            D1_y = phi(i,j,k) - phi(i,j-1,k)
            D1_y_plus = phi(i,j+1,k) - phi(i,j,k)
            phi_y_plus = (D1_y_plus)*inv_dy
            phi_y_minus = (D1_y)*inv_dy
            
c           phi_z_plus, minus
            D1_z = phi(i,j,k) - phi(i,j,k-1)
            D1_z_plus = phi(i,j,k+1) - phi(i,j,k)
            phi_z_plus = (D1_z_plus)*inv_dz
            phi_z_minus = (D1_z)*inv_dz
            
                       
	      if (abs(phi0(i,j,k)) .gt. zero_tol) then
	        sgn_phi0 = sign(one,phi0(i,j,k));
	      else
	        sgn_phi0 = zero;
	      endif
     
          if( distance0(i,j,k) .eq. large_distance_flag ) then   
c           Away from interface we do standard Godunov gradient selection
c             cache phi and spatial derivative approximations
             phi_cur = phi(i,j,k)
             grad_phi_plus_cur(1) = phi_x_plus
             grad_phi_plus_cur(2) = phi_y_plus
             grad_phi_plus_cur(3) = phi_z_plus
             grad_phi_minus_cur(1) = phi_x_minus
             grad_phi_minus_cur(2) = phi_y_minus
             grad_phi_minus_cur(3) = phi_z_minus

c             { begin Godunov selection of grad_phi
             do dir=1,DIM
               if (phi_cur .gt. 0.d0) then
                grad_phi_plus_cur(dir)=max(-grad_phi_plus_cur(dir),
     &                                         0.d0)
                grad_phi_minus_cur(dir)=max(grad_phi_minus_cur(dir),
     &                                          0.d0)
               else
                 grad_phi_plus_cur(dir) = 
     &              max(grad_phi_plus_cur(dir), 0.d0)
                 grad_phi_minus_cur(dir) = 
     &              max(-grad_phi_minus_cur(dir), 0.d0)
               endif

               grad_phi_star(dir) = max(grad_phi_plus_cur(dir),
     &                                  grad_phi_minus_cur(dir)) 

             enddo
c             } end Godunov selection of grad_phi    

c           compute reinit_rhs(i,j,k) - no need for smoothing since
c           we're away from the interface
            if (abs(phi_cur) .gt. zero_tol) then 
               norm_grad_phi_sq = grad_phi_star(1)*grad_phi_star(1)
     &                          + grad_phi_star(2)*grad_phi_star(2)
     &                          + grad_phi_star(3)*grad_phi_star(3)
               
	           reinit_rhs(i,j,k) = sgn_phi0 - sgn_phi0*
     &                  	            sqrt(norm_grad_phi_sq)
            else
               reinit_rhs(i,j,k) = 0.d0
            endif
          else

c        Close to the interface make sure you don't use information from
c        the other side of the interface (Russo/Smereka 2000. subcell fix)
	      reinit_rhs(i,j,k) = distance0(i,j,k) - 
     &                       sgn_phi0*abs(phi(i,j,k))/max_dx;
	      endif
	      
        enddo
       enddo
      enddo
c       } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************
