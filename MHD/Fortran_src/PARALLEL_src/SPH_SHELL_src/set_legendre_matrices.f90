!>@file   set_legendre_matrices.f90
!!@brief  module set_legendre_matrices
!!
!!@author H. Matsui
!!@date Programmed in Dec., 2014
!
!>@brief  set legendre polynomials into matrices
!!
!!@verbatim
!!      subroutine set_trans_legendre_rtm(nth_rtm, jmax_rlm,            &
!!     &          P_rtm, dPdt_rtm, P_jl, dPdt_jl)
!!
!!      subroutine set_sym_legendre_stack                               &
!!     &         (mphi_rtm, lstack_rlm, lstack_even_rlm)
!!      subroutine set_symmetric_legendre_lj(nth_rtm, mphi_rtm,         &
!!     &           jmax_rlm, nth_hemi_rtm, lstack_rlm, lstack_even_rlm, &
!!     &           P_rtm, dPdt_rtm, Ps_lj, dPsdt_lj)
!!      subroutine set_symmetric_legendre_jl(nth_rtm, mphi_rtm,         &
!!     &           jmax_rlm, nth_hemi_rtm, lstack_rlm, lstack_even_rlm, &
!!     &           P_rtm, dPdt_rtm, Ps_jl, dPsdt_jl)
!!
!!      subroutine symmetric_legendre_4_tensor_lj                       &
!!     &         (nth_rtm, mphi_rtm, jmax_rlm, nth_hemi_rtm, lstack_rlm,&
!!     &          lstack_even_rlm, idx_gl_1d_rlm_j, g_sph_rlm,          &
!!     &          asin_theta_1d_rtm, cos_theta_1d_rtm, Ps_lj, dPsdt_lj, &
!!     &          P2s_lj, dP2sdt_lj)
!!      subroutine symmetric_legendre_4_tensor_jl                       &
!!     &         (nth_rtm, mphi_rtm, jmax_rlm, nth_hemi_rtm, lstack_rlm,&
!!     &          lstack_even_rlm, idx_gl_1d_rlm_j, g_sph_rlm,          &
!!     &          asin_theta_1d_rtm, cos_theta_1d_rtm, Ps_jl, dPsdt_jl, &
!!     &          P2s_jl, dP2sdt_jl)
!!@endverbatim
!
      module set_legendre_matrices
!
      use m_precision
      use m_constants
!
      implicit none
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine set_trans_legendre_rtm(nth_rtm, jmax_rlm,              &
     &          P_rtm, dPdt_rtm, P_jl, dPdt_jl)
!
      integer(kind = kint), intent(in) :: nth_rtm, jmax_rlm
      real(kind= kreal), intent(in) :: P_rtm(nth_rtm,jmax_rlm)
      real(kind= kreal), intent(in) :: dPdt_rtm(nth_rtm,jmax_rlm)
!
      real(kind= kreal), intent(inout) :: P_jl(jmax_rlm,nth_rtm)
      real(kind= kreal), intent(inout) :: dPdt_jl(jmax_rlm,nth_rtm)
!
!
      integer(kind = kint) :: l_rtm, j_rlm
!
!
!$omp parallel do private(j_rlm,l_rtm)
      do j_rlm = 1, jmax_rlm
        do l_rtm = 1, nth_rtm
          P_jl(j_rlm,l_rtm) =     P_rtm(l_rtm,j_rlm)
          dPdt_jl(j_rlm,l_rtm) =  dPdt_rtm(l_rtm,j_rlm)
        end do
      end do
!$omp end parallel do
!
      end subroutine set_trans_legendre_rtm
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine set_sym_legendre_stack                                 &
     &         (mphi_rtm, lstack_rlm, lstack_even_rlm)
!
      integer(kind = kint), intent(in) :: mphi_rtm
      integer(kind = kint), intent(in) :: lstack_rlm(0:mphi_rtm)
!
      integer(kind = kint), intent(inout)                               &
     &                     :: lstack_even_rlm(0:mphi_rtm)
!
      integer(kind = kint) :: mp_rlm, jst, nj_rlm
!
!
      do mp_rlm = 1, mphi_rtm
        jst = lstack_rlm(mp_rlm-1)
        nj_rlm = lstack_rlm(mp_rlm) - lstack_rlm(mp_rlm-1)
        lstack_even_rlm(mp_rlm) = jst + (nj_rlm+1) / 2
      end do
!
      end subroutine set_sym_legendre_stack
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine set_symmetric_legendre_lj(nth_rtm, mphi_rtm,           &
     &           jmax_rlm, nth_hemi_rtm, lstack_rlm, lstack_even_rlm,   &
     &           P_rtm, dPdt_rtm, Ps_lj, dPsdt_lj)
!
      integer(kind = kint), intent(in) :: nth_rtm, mphi_rtm, jmax_rlm
      integer(kind = kint), intent(in) :: nth_hemi_rtm
      integer(kind = kint), intent(in) :: lstack_rlm(0:mphi_rtm)
      integer(kind = kint), intent(in) :: lstack_even_rlm(0:mphi_rtm)
!
      real(kind= kreal), intent(in) :: P_rtm(nth_rtm,jmax_rlm)
      real(kind= kreal), intent(in) :: dPdt_rtm(nth_rtm,jmax_rlm)
!
      real(kind= kreal), intent(inout) :: Ps_lj(nth_hemi_rtm,jmax_rlm)
      real(kind= kreal), intent(inout)                                  &
     &                  :: dPsdt_lj(nth_hemi_rtm,jmax_rlm)
!
      integer(kind = kint) :: l_rtm, j_rlm
      integer(kind = kint) :: mp_rlm, jst, n_jk_e, n_jk_o, jj
!
!
!$omp parallel do private(jst,j_rlm,l_rtm,jj,n_jk_e,n_jk_o)
      do mp_rlm = 1, mphi_rtm
        jst = lstack_rlm(mp_rlm-1)
        n_jk_e = lstack_even_rlm(mp_rlm) - lstack_rlm(mp_rlm-1)
        n_jk_o = lstack_rlm(mp_rlm) - lstack_even_rlm(mp_rlm)
!
        do jj = 1, n_jk_e
          j_rlm = 2*jj + jst - 1
          do l_rtm = 1, nth_hemi_rtm
            Ps_lj(l_rtm,jj+jst) =     P_rtm(l_rtm,j_rlm)
            dPsdt_lj(l_rtm,jj+jst) =  dPdt_rtm(l_rtm,j_rlm)
          end do
        end do
!
        do jj = 1, n_jk_o
          j_rlm = 2*jj + jst
          do l_rtm = 1, nth_hemi_rtm
            Ps_lj(l_rtm,jj+jst+n_jk_e) =     P_rtm(l_rtm,j_rlm)
            dPsdt_lj(l_rtm,jj+jst+n_jk_e) =  dPdt_rtm(l_rtm,j_rlm)
          end do
        end do
      end do
!$omp end parallel do
!
      end subroutine set_symmetric_legendre_lj
!
! -----------------------------------------------------------------------
!
      subroutine set_symmetric_legendre_jl(nth_rtm, mphi_rtm,           &
     &           jmax_rlm, nth_hemi_rtm, lstack_rlm, lstack_even_rlm,   &
     &           P_rtm, dPdt_rtm, Ps_jl, dPsdt_jl)
!
      integer(kind = kint), intent(in) :: nth_rtm, mphi_rtm, jmax_rlm
      integer(kind = kint), intent(in) :: nth_hemi_rtm
      integer(kind = kint), intent(in) :: lstack_rlm(0:mphi_rtm)
      integer(kind = kint), intent(in) :: lstack_even_rlm(0:mphi_rtm)
!
      real(kind= kreal), intent(in) :: P_rtm(nth_rtm,jmax_rlm)
      real(kind= kreal), intent(in) :: dPdt_rtm(nth_rtm,jmax_rlm)
!
      real(kind= kreal), intent(inout) :: Ps_jl(jmax_rlm,nth_hemi_rtm)
      real(kind= kreal), intent(inout)                                  &
     &                  :: dPsdt_jl(jmax_rlm,nth_hemi_rtm)
!
      integer(kind = kint) :: l_rtm, j_rlm
      integer(kind = kint) :: mp_rlm, jst, n_jk_e, n_jk_o, jj
!
!
!$omp parallel do private(jst,j_rlm,l_rtm,jj,n_jk_e,n_jk_o)
      do mp_rlm = 1, mphi_rtm
        jst = lstack_rlm(mp_rlm-1)
        n_jk_e = lstack_even_rlm(mp_rlm) - lstack_rlm(mp_rlm-1)
        n_jk_o = lstack_rlm(mp_rlm) - lstack_even_rlm(mp_rlm)
!
        do jj = 1, n_jk_e
          j_rlm = 2*jj + jst - 1
          do l_rtm = 1, nth_hemi_rtm
            Ps_jl(jj+jst,l_rtm) =     P_rtm(l_rtm,j_rlm)
            dPsdt_jl(jj+jst,l_rtm) =  dPdt_rtm(l_rtm,j_rlm)
          end do
        end do
!
        do jj = 1, n_jk_o
          j_rlm = 2*jj + jst
          do l_rtm = 1, nth_hemi_rtm
            Ps_jl(jj+jst+n_jk_e,l_rtm) =     P_rtm(l_rtm,j_rlm)
            dPsdt_jl(jj+jst+n_jk_e,l_rtm) =  dPdt_rtm(l_rtm,j_rlm)
          end do
        end do
      end do
!$omp end parallel do
!
      end subroutine set_symmetric_legendre_jl
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine symmetric_legendre_4_tensor_lj                         &
     &         (nth_rtm, mphi_rtm, jmax_rlm, nth_hemi_rtm, lstack_rlm,  &
     &          lstack_even_rlm, idx_gl_1d_rlm_j, g_sph_rlm,            &
     &          asin_theta_1d_rtm, cos_theta_1d_rtm, Ps_lj, dPsdt_lj,   &
     &          P2s_lj, dP2sdt_lj)
!
      integer(kind = kint), intent(in) :: nth_rtm, mphi_rtm, jmax_rlm
      integer(kind = kint), intent(in) :: nth_hemi_rtm
      integer(kind = kint), intent(in) :: lstack_rlm(0:mphi_rtm)
      integer(kind = kint), intent(in) :: lstack_even_rlm(0:mphi_rtm)
!
      integer(kind = kint), intent(in) :: idx_gl_1d_rlm_j(jmax_rlm,3)
      real(kind = kreal), intent(in) :: g_sph_rlm(jmax_rlm,17)
      real(kind= kreal), intent(in) :: asin_theta_1d_rtm(nth_rtm)
      real(kind= kreal), intent(in) :: cos_theta_1d_rtm(nth_rtm)
!
      real(kind= kreal), intent(in) :: Ps_lj(nth_hemi_rtm,jmax_rlm)
      real(kind= kreal), intent(in)                                  &
     &                  :: dPsdt_lj(nth_hemi_rtm,jmax_rlm)
!
      real(kind= kreal), intent(inout) :: P2s_lj(nth_hemi_rtm,jmax_rlm)
      real(kind= kreal), intent(inout)                                  &
     &                  :: dP2sdt_lj(nth_hemi_rtm,jmax_rlm)
!
      integer(kind = kint) :: l_rtm, j_rlm
      integer(kind = kint) :: mp_rlm, jst, n_jk_e, n_jk_o, jj
      real(kind = kreal) :: g3, gm, ast, cst
!
!
!$omp parallel do                                                       &
!$omp&  private(jst,j_rlm,l_rtm,jj,n_jk_e,n_jk_o,g3,gm,ast,cst)
      do mp_rlm = 1, mphi_rtm
        jst = lstack_rlm(mp_rlm-1)
        n_jk_e = lstack_even_rlm(mp_rlm) - lstack_rlm(mp_rlm-1)
        n_jk_o = lstack_rlm(mp_rlm) - lstack_even_rlm(mp_rlm)
!
        do jj = 1, n_jk_e
          j_rlm = 2*jj + jst - 1
          g3 = g_sph_rlm(j_rlm,3)
          gm = dble(idx_gl_1d_rlm_j(j_rlm,3))
          do l_rtm = 1, nth_hemi_rtm
            ast = asin_theta_1d_rtm(l_rtm)
            cst = cos_theta_1d_rtm(l_rtm)
            P2s_lj(l_rtm,jj+jst)                                        &
     &       = (-g3 + two * gm**2 * ast**2) * Ps_lj(l_rtm,jj+jst)       &
     &        - two * ast * cst * dPsdt_lj(l_rtm,jj+jst)
            dP2sdt_lj(l_rtm,jj+jst)                                     &
     &       = - two * gm  * (ast * dPsdt_lj(l_rtm,jj+jst)              &
     &        - ast**2 * cst * Ps_lj(l_rtm,jj+jst))
          end do
        end do
!
        do jj = 1, n_jk_o
          j_rlm = 2*jj + jst
          g3 = g_sph_rlm(j_rlm,3)
          gm = dble(idx_gl_1d_rlm_j(j_rlm,3))
          do l_rtm = 1, nth_hemi_rtm
            ast = asin_theta_1d_rtm(l_rtm)
            cst = cos_theta_1d_rtm(l_rtm)
            P2s_lj(l_rtm,jj+jst+n_jk_e)                                 &
     &       =(-g3 + two * gm**2 * ast**2) * Ps_lj(l_rtm,jj+jst+n_jk_e) &
     &            - two * ast * cst * dPsdt_lj(l_rtm,jj+jst+n_jk_e)
            dP2sdt_lj(l_rtm,jj+jst+n_jk_e)                              &
     &       = - two * gm * (ast * dPsdt_lj(l_rtm,jj+jst+n_jk_e)        &
     &        - ast**2 * cst * Ps_lj(l_rtm,jj+jst+n_jk_e))
          end do
        end do
      end do
!$omp end parallel do
!
      end subroutine symmetric_legendre_4_tensor_lj
!
! -----------------------------------------------------------------------
!
      subroutine symmetric_legendre_4_tensor_jl                         &
     &         (nth_rtm, mphi_rtm, jmax_rlm, nth_hemi_rtm, lstack_rlm,  &
     &          lstack_even_rlm, idx_gl_1d_rlm_j, g_sph_rlm,            &
     &          asin_theta_1d_rtm, cos_theta_1d_rtm, Ps_jl, dPsdt_jl,   &
     &          P2s_jl, dP2sdt_jl)
!
      integer(kind = kint), intent(in) :: nth_rtm, mphi_rtm, jmax_rlm
      integer(kind = kint), intent(in) :: nth_hemi_rtm
      integer(kind = kint), intent(in) :: lstack_rlm(0:mphi_rtm)
      integer(kind = kint), intent(in) :: lstack_even_rlm(0:mphi_rtm)
!
      integer(kind = kint), intent(in) :: idx_gl_1d_rlm_j(jmax_rlm,3)
      real(kind = kreal), intent(in) :: g_sph_rlm(jmax_rlm,17)
      real(kind= kreal), intent(in) :: asin_theta_1d_rtm(nth_rtm)
      real(kind= kreal), intent(in) :: cos_theta_1d_rtm(nth_rtm)
!
      real(kind= kreal), intent(in) :: Ps_jl(jmax_rlm,nth_hemi_rtm)
      real(kind= kreal), intent(in)                                  &
     &                  :: dPsdt_jl(jmax_rlm,nth_hemi_rtm)
!
      real(kind= kreal), intent(inout) :: P2s_jl(jmax_rlm,nth_hemi_rtm)
      real(kind= kreal), intent(inout)                                  &
     &                  :: dP2sdt_jl(jmax_rlm,nth_hemi_rtm)
!
      integer(kind = kint) :: l_rtm, j_rlm
      integer(kind = kint) :: mp_rlm, jst, n_jk_e, n_jk_o, jj
      real(kind = kreal) :: g3, gm, ast, cst
!
!
!$omp parallel do                                                       &
!$omp& private(jst,j_rlm,l_rtm,jj,n_jk_e,n_jk_o,g3,gm,ast,cst)
      do mp_rlm = 1, mphi_rtm
        jst = lstack_rlm(mp_rlm-1)
        n_jk_e = lstack_even_rlm(mp_rlm) - lstack_rlm(mp_rlm-1)
        n_jk_o = lstack_rlm(mp_rlm) - lstack_even_rlm(mp_rlm)
!
        do jj = 1, n_jk_e
          j_rlm = 2*jj + jst - 1
          g3 = g_sph_rlm(j_rlm,3)
          gm = dble(idx_gl_1d_rlm_j(j_rlm,3))
          do l_rtm = 1, nth_hemi_rtm
            ast = asin_theta_1d_rtm(l_rtm)
            cst = cos_theta_1d_rtm(l_rtm)
            P2s_jl(jj+jst,l_rtm)                                        &
     &       = (-g3 + two * gm**2 * ast**2) * Ps_jl(jj+jst,l_rtm)       &
     &        - two * ast * cst * dPsdt_jl(jj+jst,l_rtm)
            dP2sdt_jl(jj+jst,l_rtm)                                     &
     &       = - two * gm  * (ast * dPsdt_jl(jj+jst,l_rtm)              &
     &        - ast**2 * cst * Ps_jl(jj+jst,l_rtm))
          end do
        end do
!
        do jj = 1, n_jk_o
          j_rlm = 2*jj + jst
          g3 = g_sph_rlm(j_rlm,3)
          gm = dble(idx_gl_1d_rlm_j(j_rlm,3))
          do l_rtm = 1, nth_hemi_rtm
            ast = asin_theta_1d_rtm(l_rtm)
            cst = cos_theta_1d_rtm(l_rtm)
            P2s_jl(jj+jst+n_jk_e,l_rtm)                                 &
     &       =(-g3 + two * gm**2 * ast**2) * Ps_jl(jj+jst+n_jk_e,l_rtm) &
     &        - two * ast * cst * dPsdt_jl(jj+jst+n_jk_e,l_rtm)
            dP2sdt_jl(jj+jst+n_jk_e,l_rtm)                              &
     &       = - two * gm  * (ast * dPsdt_jl(jj+jst+n_jk_e,l_rtm)       &
     &        - ast**2 * cst *Ps_jl(jj+jst+n_jk_e,l_rtm))
          end do
        end do
      end do
!$omp end parallel do
!
      end subroutine symmetric_legendre_4_tensor_jl
!
! -----------------------------------------------------------------------
!
      end module set_legendre_matrices
