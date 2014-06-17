!>@file   legendre_fwd_trans_org.f90
!!@brief  module legendre_fwd_trans_org
!!
!!@author H. Matsui
!!@date Programmed in Aug., 2007
!!@n    Modified in Apr. 2013
!
!>@brief  forward Legendre transform
!!       (Original version)
!!
!!@verbatim
!!      subroutine legendre_f_trans_vector_org(ncomp, nvector)
!!        Input:  vr_rtm   (Order: radius,theta,phi)
!!        Output: sp_rlm   (Order: poloidal,diff_poloidal,toroidal)
!!      subroutine legendre_f_trans_scalar_org(ncomp, nvector, nscalar)
!!        Input:  vr_rtm
!!        Output: sp_rlm
!!@endverbatim
!!
!!@param   ncomp    Total number of components for spherical transform
!!@param   nvector  Number of vector for spherical transform
!!@param   nscalar  Number of scalar (including tensor components)
!!                  for spherical transform
!
      module legendre_fwd_trans_org
!
      use m_precision
!
      use m_machine_parameter
      use m_spheric_parameter
      use m_spheric_param_smp
      use m_schmidt_poly_on_rtm
      use m_work_4_sph_trans
!
      implicit none
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine legendre_f_trans_vector_org(ncomp, nvector)
!
      integer(kind = kint), intent(in) :: ncomp, nvector
!
      integer(kind = kint) :: i_rlm, k_rlm, j_rlm
      integer(kind = kint) :: l_rtm
      integer(kind = kint) :: ip_rtm, in_rtm
      integer(kind = kint) :: nd
      real(kind = kreal) :: r2_1d_rlm_r
!
!
!$omp parallel do private(l_rtm,j_rlm,k_rlm,nd,i_rlm,ip_rtm,in_rtm,     &
!$omp&                    r2_1d_rlm_r)
      do k_rlm = 1, nidx_rlm(1)
        r2_1d_rlm_r = radius_1d_rlm_r(k_rlm)*radius_1d_rlm_r(k_rlm)
!
        do j_rlm = 1, nidx_rlm(2)
          do nd = 1, nvector
            i_rlm = 3*nd + (j_rlm-1) * ncomp                            &
     &                 + (k_rlm-1) * ncomp*nidx_rlm(2)
!
            do l_rtm = 1, nidx_rtm(2)
              ip_rtm = 3*nd + (l_rtm-1)  * ncomp                        &
     &                 + (k_rlm-1)  * ncomp * nidx_rtm(2)               &
     &                 + (mdx_p_rlm_rtm(j_rlm)-1)                       &
     &                  * ncomp * nidx_rtm(1) * nidx_rtm(2)
              in_rtm = 3*nd + (l_rtm-1)  * ncomp                        &
     &                 + (k_rlm-1)  * ncomp * nidx_rtm(2)               &
     &                 + (mdx_n_rlm_rtm(j_rlm)-1)                       &
     &                  * ncomp * nidx_rtm(1) * nidx_rtm(2)
!
              sp_rlm(i_rlm-2) = sp_rlm(i_rlm-2)                         &
     &                     + vr_rtm(ip_rtm-2) * Pvw_lj(l_rtm,j_rlm)
              sp_rlm(i_rlm-1) = sp_rlm(i_rlm-1)                         &
     &                 + ( vr_rtm(ip_rtm-1) * dPvw_lj(l_rtm,j_rlm)      &
     &                   - vr_rtm(in_rtm  ) * Pgvw_lj(l_rtm,j_rlm))
              sp_rlm(i_rlm  ) = sp_rlm(i_rlm  )                         &
     &                 - ( vr_rtm(in_rtm-1) * Pgvw_lj(l_rtm,j_rlm)      &
     &                   + vr_rtm(ip_rtm  ) * dPvw_lj(l_rtm,j_rlm))
            end do
!
            sp_rlm(i_rlm-2) = sp_rlm(i_rlm-2) * r2_1d_rlm_r
            sp_rlm(i_rlm-1) = sp_rlm(i_rlm-1) * radius_1d_rlm_r(k_rlm)
            sp_rlm(i_rlm  ) = sp_rlm(i_rlm  ) * radius_1d_rlm_r(k_rlm)
          end do
        end do
      end do
!$omp end parallel do
!
      end subroutine legendre_f_trans_vector_org
!
! -----------------------------------------------------------------------
!
      subroutine legendre_f_trans_scalar_org(ncomp, nvector, nscalar)
!
      integer(kind = kint), intent(in) :: ncomp, nvector, nscalar
!
      integer(kind = kint) :: i_rlm, k_rlm, j_rlm
      integer(kind = kint) :: l_rtm
      integer(kind = kint) :: ip_rtm
      integer(kind = kint) :: nd
!
!
!$omp parallel do private(j_rlm,k_rlm,nd,i_rlm,ip_rtm,l_rtm)
      do k_rlm = 1, nidx_rlm(1)
        do j_rlm = 1, nidx_rlm(2)
!
          do nd = 1, nscalar
            i_rlm = nd + 3*nvector + (j_rlm-1) * ncomp                  &
     &                             + (k_rlm-1) * ncomp*nidx_rlm(2)
!
            do l_rtm = 1, nidx_rtm(2)
              ip_rtm = nd + 3*nvector + (l_rtm-1)  * ncomp              &
     &                + (k_rlm-1)  * ncomp * nidx_rtm(2)                &
     &                + (mdx_p_rlm_rtm(j_rlm)-1)                        &
     &                  * ncomp * nidx_rtm(1) * nidx_rtm(2)
!
              sp_rlm(i_rlm) = sp_rlm(i_rlm)                             &
     &                       + vr_rtm(ip_rtm) * Pws_lj(l_rtm,j_rlm)
            end do
!
          end do
        end do
      end do
!$omp end parallel do
!
      end subroutine legendre_f_trans_scalar_org
!
! -----------------------------------------------------------------------
!
      end module legendre_fwd_trans_org
