!int_edge_commute_z_filter.f90
!      module int_edge_commute_z_filter
!
!      Written by H. Matsui
!
!
!!      subroutine int_edge_commutative_filter(numnod, numele,          &
!!     &          zz, ie_edge, gauss, g_int)
!!        type(gauss_points), intent(in) :: gauss
!!        type(gauss_integrations), intent(inout) :: g_int
!
      module int_edge_commute_z_filter
!
      use m_precision
      use t_gauss_points
!
      implicit none
!
!   --------------------------------------------------------------------
!
      contains
!
!   --------------------------------------------------------------------
!
      subroutine int_edge_commutative_filter(numnod, numele,            &
     &          zz, ie_edge, gauss, g_int)
!
      use m_constants
      use m_commute_filter_z
      use m_int_commtative_filter
      use m_z_filter_values
      use m_work_4_integration
      use m_int_edge_data
      use set_filter_moments
!
      integer(kind = kint), intent(in) :: numnod, numele
      integer(kind = kint), intent(in) :: ie_edge(numele,2)
      real(kind = kreal), intent(in) :: zz(numnod)
      type(gauss_points), intent(in) :: gauss
!
      type(gauss_integrations), intent(inout) :: g_int
!
      integer (kind = kint) :: inod, jnod, kf
      integer (kind = kint) :: i, j, je, jele
      integer (kind = kint) :: jnod1, jnod2, j0
!
      real(kind = kreal) :: zz0, zz1, zz2, zs, ze
      real(kind = kreal) :: filter_0(gauss%n_point)
!
!
      do inod = 1, numnod
        zz0 = zz(inod)
        do je = 1, ncomp_mat-1
          j0 = je - nneib_nod2(inod,1) - 1
          jele = inod + j0
          zs = dble(2*(j0)  )
          ze = dble(2*(j0+1))
          jnod1 = ie_edge(jele,1)
          jnod2 = ie_edge(jele,2)
          zz1 = zz(jnod1)
          zz2 = zz(jnod2)
!
          call set_points_4_integration(zs, ze, gauss, g_int)
!
          do j = 1, 2
            jnod = je + j - 1
!
            if ( iflag_filter .eq. 0) then
              call filter_moment_tophat(izero, gauss%n_point,           &
     &            f_width, filter_0, g_int%x_point)
            else if (iflag_filter .eq. 1) then
              call filter_moment_linear(izero, gauss%n_point,           &
     &            f_width, filter_0, g_int%x_point)
            else
              call filter_moment_gaussian(izero, gauss%n_point,         &
     &            f_width, filter_0, g_int%x_point)
            end if
!
            do i = 1, gauss%n_point
              g_int%f_point(1,i) = half * filter_0(i)                   &
     &             * (one + (-1)**j * (g_int%x_point(i)-dble(2*j0+1))) 
             do kf = 2, nfilter6_1+1
               g_int%f_point(kf,i) = half * g_int%f_point(kf-1,i)       &
     &             * ( zz2 + zz1 - two*zz0                              &
     &              + ( (zz2-zz1)*(g_int%x_point(i)-dble(2*j0+1)) ))
             end do
            end do
!
            call cal_gauss_integrals(gauss, g_int, sk_norm_n(0))
!
            do kf = 0, 2
              xmom_int_org(inod,jnod,kf) = xmom_int_org(inod,jnod,kf)   &
     &                         + sk_norm_n(kf)
            end do
!
!
            do i = 1, gauss%n_point
              g_int%f_point(1,i) = filter_0(i)  * quad * dz(jele)       &
     &            * (one + (-1)**j * (g_int%x_point(i)-dble(2*j0+1)))   &
     &             * ( c_filter(je,inod)                                &
     &             * ( one - (g_int%x_point(i)-dble(2*j0+1)) )          &
     &              + c_filter(je+1,inod)                               &
     &             * ( one + (g_int%x_point(i)-dble(2*j0+1)) ) )
             do kf = 2, nfilter6_1+1
               g_int%f_point(kf,i) = half * g_int%f_point(kf-1,i)       &
     &            * ( zz2 + zz1 - two*zz0                               &
     &             + ( (zz2-zz1)*(g_int%x_point(i)-dble(2*j0+1)) ))
             end do
            end do
!
            call cal_gauss_integrals(gauss, g_int, sk_norm_n(0))
!
            do kf = 0, 2
              xmom_int(inod,jnod,kf) = xmom_int(inod,jnod,kf)           &
     &                         + sk_norm_n(kf)
            end do
!
          end do
        end do
      end do
!
!
      do kf = 0, 2
       do inod = 1, numnod
        do jnod = 1, ncomp_mat
          xmom_int_t(inod,kf) = xmom_int_t(inod,kf)                     &
     &                          + xmom_int(inod,jnod,kf)
          xmom_int_to(inod,kf) = xmom_int_to(inod,kf)                   &
     &                          + xmom_int_org(inod,jnod,kf)
        end do
       end do
      end do
!
      end subroutine int_edge_commutative_filter
!
!   --------------------------------------------------------------------
!
      end module int_edge_commute_z_filter
