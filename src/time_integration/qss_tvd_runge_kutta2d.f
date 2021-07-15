c***********************************************************************
      subroutine qss2dRK1Step(
     &  u_next,
     &  ilo_u_next_gb, ihi_u_next_gb,
     &  jlo_u_next_gb, jhi_u_next_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dt)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_u_next_gb, ihi_u_next_gb
      integer jlo_u_next_gb, jhi_u_next_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real u_next(ilo_u_next_gb:ihi_u_next_gb,
     &                        jlo_u_next_gb:jhi_u_next_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb)
      integer i, j
      real dt

c     { begin loop over grid
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

            u_next(i,j) = u_cur(i,j) + dt*rhs(i,j)

          enddo
        enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************
