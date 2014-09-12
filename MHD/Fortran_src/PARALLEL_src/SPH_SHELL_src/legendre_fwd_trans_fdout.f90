!>@file   legendre_fwd_trans_fdout.f90
!!@brief  module legendre_fwd_trans_fdout
!!
!!@author H. Matsui
!!@date Programmed in Aug., 2007
!!@n    Modified in Apr. 2013
!
!>@brief  forward Legendre transform
!!       (outmost field loop version)
!!
!!@verbatim
!!      subroutine legendre_f_trans_vector_fdout(nvector,               &
!!     &          vr_rtm_fdout, sp_rlm_fdout)
!!      subroutine legendre_f_trans_scalar_fdout(nscalar,               &
!!     &          vr_rtm_fdout, sp_rlm_fdout)
!!        Input:  vr_rtm_fdout
!!        Output: sp_rlm_fdout
!!
!!      field data for Legendre transform  @f$ f(r,\theta,m) @f$ 
!!       Order: vr_rtm_fdout(l_rtm,k_rtm,m_rtm,i_comp,i_fld)
!!       size:  vr_rtm_fdout(nidx_rtm(2)*nidx_rtm(1)*nidx_rtm(3),3*nb)
!!      real(kind = kreal), allocatable :: vr_rtm_fdout(:,:)
!!
!!      Spectr data for Legendre transform  @f$ f(r,l,m) @f$ 
!!        Order: sp_rlm(i_comp,i_fld,j_rlm,k_rtm)
!!        size: sp_rlm(nidx_rlm(2)*nidx_rtm(1),3*nb)
!!      real(kind = kreal), allocatable :: sp_rlm_fdout(:,:)
!!@endverbatim
!!
!!@param   nvector  Number of vector for spherical transform
!!@param   nscalar  Number of scalar (including tensor components)
!!                  for spherical transform
!
      module legendre_fwd_trans_fdout
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
      subroutine legendre_f_trans_vector_fdout(nvector,                 &
     &          vr_rtm_fdout, sp_rlm_fdout)
!
      integer(kind = kint), intent(in) :: nvector
      real(kind = kreal), intent(in)                                    &
     &      :: vr_rtm_fdout(nnod_rtm,3*nvector)
      real(kind = kreal), intent(inout)                                 &
     &      :: sp_rlm_fdout(nnod_rlm,3*nvector)
!
      integer(kind = kint) :: i_rlm, k_rlm, j_rlm
      integer(kind = kint) :: l_rtm
      integer(kind = kint) :: ip_rtm, in_rtm
      integer(kind = kint) :: nd, ip, kst, ked, lp, lst, led
      real(kind = kreal) :: sp1, sp2, sp3
      real(kind = kreal) :: Pvw_l(nidx_rtm(2))
      real(kind = kreal) :: dPvw_l(nidx_rtm(2))
      real(kind = kreal) :: Pgvw_l(nidx_rtm(2))
!
!
!$omp parallel do schedule(static)                                      &
!$omp&            private(ip,kst,ked,lp,lst,led,nd,l_rtm,j_rlm,         &
!$omp&                    k_rlm,i_rlm,ip_rtm,in_rtm,                    &
!$omp&                    sp1,sp2,sp3,Pvw_l,dPvw_l,Pgvw_l)
      do ip = 1, np_smp
        kst = idx_rlm_smp_stack(ip-1,1) + 1
        ked = idx_rlm_smp_stack(ip,  1)
        do nd = 1, nvector
          do k_rlm = kst, ked
            do j_rlm = 1, nidx_rlm(2)
              i_rlm = j_rlm + (k_rlm-1) * nidx_rlm(2)
              sp_rlm_fdout(i_rlm,3*nd-2) = 0.0d0
              sp_rlm_fdout(i_rlm,3*nd-1) = 0.0d0
              sp_rlm_fdout(i_rlm,3*nd  ) = 0.0d0
            end do
          end do
        end do
!
        do lp = 1, nblock_l_rtm
          lst = lstack_block_rtm(lp-1) + 1
          led = lstack_block_rtm(lp  )
!
          do nd = 1, nvector
            do j_rlm = 1, nidx_rlm(2)
              do l_rtm = lst, led
                Pvw_l(l_rtm) = P_rtm(l_rtm,j_rlm)                       &
     &                     * g_sph_rlm(j_rlm,7)* weight_rtm(l_rtm)
                dPvw_l(l_rtm) = dPdt_rtm(l_rtm,j_rlm)                   &
     &                     * g_sph_rlm(j_rlm,7)* weight_rtm(l_rtm)
                Pgvw_l(l_rtm) = P_rtm(l_rtm,j_rlm)                      &
     &                     * dble(idx_gl_1d_rlm_j(j_rlm,3))             &
     &                      * asin_theta_1d_rtm(l_rtm)                  &
     &                      * g_sph_rlm(j_rlm,7)* weight_rtm(l_rtm)
              end do
!
              do k_rlm = kst, ked
                sp1 = 0.0d0
                sp2 = 0.0d0
                sp3 = 0.0d0
                do l_rtm = lst, led
                  ip_rtm = 1 + (l_rtm-1) * istep_rtm(2)                 &
     &                       + (k_rlm-1) * istep_rtm(1)                 &
     &                       + (mdx_p_rlm_rtm(j_rlm)-1) * istep_rtm(3)
                  in_rtm = 1 + (l_rtm-1) * istep_rtm(2)                 &
     &                       + (k_rlm-1) * istep_rtm(1)                 &
     &                       + (mdx_n_rlm_rtm(j_rlm)-1) * istep_rtm(3)
!
                  sp1 = sp1 + vr_rtm_fdout(ip_rtm,3*nd-2)               &
     &                        * Pvw_l(l_rtm)
                  sp2 = sp2 + (vr_rtm_fdout(ip_rtm,3*nd-1)              &
     &                        * dPvw_l(l_rtm)                           &
     &                       - vr_rtm_fdout(in_rtm,3*nd  )              &
     &                        * Pgvw_l(l_rtm))
                  sp3 = sp3 - (vr_rtm_fdout(in_rtm,3*nd-1)              &
     &                        * Pgvw_l(l_rtm)                           &
     &                       + vr_rtm_fdout(ip_rtm,3*nd  )              &
     &                        * dPvw_l(l_rtm))
                end do
!
                i_rlm = j_rlm + (k_rlm-1) * nidx_rlm(2)
                sp_rlm_fdout(i_rlm,3*nd-2)                              &
     &                = sp_rlm_fdout(i_rlm,3*nd-2)                      &
     &                 + sp1 * radius_1d_rlm_r(k_rlm)                   &
     &                        *radius_1d_rlm_r(k_rlm)
                sp_rlm_fdout(i_rlm,3*nd-1)                              &
     &                = sp_rlm_fdout(i_rlm,3*nd-1)                      &
     &                 + sp2 * radius_1d_rlm_r(k_rlm)
                sp_rlm_fdout(i_rlm,3*nd  )                              &
     &                = sp_rlm_fdout(i_rlm,3*nd  )                      &
     &                 + sp3 * radius_1d_rlm_r(k_rlm)
              end do
            end do
          end do
!
        end do
      end do
!$omp end parallel do
!
      end subroutine legendre_f_trans_vector_fdout
!
! -----------------------------------------------------------------------
!
      subroutine legendre_f_trans_scalar_fdout(nscalar,                 &
     &          vr_rtm_fdout, sp_rlm_fdout)
!
      integer(kind = kint), intent(in) :: nscalar
      real(kind = kreal), intent(in)                                    &
     &      :: vr_rtm_fdout(nnod_rtm,nscalar)
      real(kind = kreal), intent(inout)                                 &
     &      :: sp_rlm_fdout(nnod_rlm,nscalar)
!
      integer(kind = kint) :: i_rlm, k_rlm, j_rlm
      integer(kind = kint) :: l_rtm
      integer(kind = kint) :: ip_rtm
      integer(kind = kint) :: nd, ip, kst, ked, lp, lst, led
      real(kind = kreal) :: sp1
      real(kind = kreal) :: Pws_l(nidx_rtm(2))
!
!
!$omp parallel do schedule(static)                                      &
!$omp&            private(nd,ip,kst,ked,lp,lst,led,j_rlm,k_rlm,         &
!$omp&                    i_rlm,l_rtm,ip_rtm,sp1,Pws_l)
      do ip = 1, np_smp
        kst = idx_rlm_smp_stack(ip-1,1) + 1
        ked = idx_rlm_smp_stack(ip,  1)
        do nd = 1, nscalar
          do k_rlm = kst, ked
            do j_rlm = 1, nidx_rlm(2)
              i_rlm = j_rlm + (k_rlm-1) * nidx_rlm(2)
              sp_rlm_fdout(i_rlm,nd) = 0.0d0
            end do
          end do
        end do
!
        do lp = 1, nblock_l_rtm
          lst = lstack_block_rtm(lp-1) + 1
          led = lstack_block_rtm(lp  )
!
          do nd = 1, nscalar
            do j_rlm = 1, nidx_rlm(2)
              do l_rtm = lst, led
                Pws_l(l_rtm) = P_rtm(l_rtm,j_rlm)                       &
     &                        * g_sph_rlm(j_rlm,6)*weight_rtm(l_rtm)
              end do
!
              do k_rlm = kst, ked
                sp1 = 0.0d0
                do l_rtm = lst, led
                  ip_rtm = 1 + (l_rtm-1) * istep_rtm(2)                 &
     &                       + (k_rlm-1) * istep_rtm(1)                 &
     &                       + (mdx_p_rlm_rtm(j_rlm)-1) * istep_rtm(3)
!
                  sp1 = sp1 + vr_rtm_fdout(ip_rtm,nd) * Pws_l(l_rtm)
                end do
                i_rlm = j_rlm + (k_rlm-1) * nidx_rlm(2)
                sp_rlm_fdout(i_rlm,nd) = sp_rlm_fdout(i_rlm,nd) + sp1
              end do
            end do
          end do
!
        end do
      end do
!$omp end parallel do
!
      end subroutine legendre_f_trans_scalar_fdout
!
! -----------------------------------------------------------------------
!
      end module legendre_fwd_trans_fdout
