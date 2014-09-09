!>@file   legendre_fwd_trans_testloop.f90
!!@brief  module legendre_fwd_trans_testloop
!!
!!@author H. Matsui
!!@date Programmed in Aug., 2013
!
!>@brief  forward Legendre transform for testing
!!
!!@verbatim
!!      subroutine alloc_vec_fleg_mat_test(nvector)
!!      subroutine alloc_scl_fleg_mat_test(nscalar)
!!      subroutine dealloc_vec_fleg_mat_test
!!      subroutine dealloc_scl_fleg_mat_test
!!
!!      subroutine legendre_f_trans_vector_test(ncomp, nvector)
!!        Input:  vr_rtm   (Order: radius,theta,phi)
!!        Output: sp_rlm   (Order: poloidal,diff_poloidal,toroidal)
!!      subroutine legendre_f_trans_scalar_test(ncomp, nvector, nscalar)
!!        Input:  vr_rtm
!!        Output: sp_rlm
!!@endverbatim
!!
!!@param   ncomp    Total number of components for spherical transform
!!@param   nvector  Number of vector for spherical transform
!!@param   nscalar  Number of scalar (including tensor components)
!!                  for spherical transform
!
      module legendre_fwd_trans_testloop
!
      use m_precision
      use m_constants
!
      use m_machine_parameter
      use m_spheric_parameter
      use m_spheric_param_smp
      use m_schmidt_poly_on_rtm
      use m_work_4_sph_trans
!
      implicit none
!
      integer(kind = kint), private :: num_jl
      real(kind = kreal), allocatable, private :: Pvw_le(:,:)
      real(kind = kreal), allocatable, private :: dPvw_le(:,:)
      real(kind = kreal), allocatable, private :: Pgvw_le(:,:)
!
      real(kind = kreal), allocatable, private :: Pws_le(:,:)
!
      integer(kind = kint), private :: nvec_jk
      real(kind = kreal), allocatable, private :: pol_e(:,:)
      real(kind = kreal), allocatable, private :: dpl_e(:,:)
      real(kind = kreal), allocatable, private :: tor_e(:,:)
      real(kind = kreal), allocatable, private :: dpl_o(:,:)
      real(kind = kreal), allocatable, private :: tor_o(:,:)
!
      integer(kind = kint), private :: nvec_lk
      real(kind = kreal), allocatable, private :: symp_r(:,:)
      real(kind = kreal), allocatable, private :: asmp_t(:,:)
      real(kind = kreal), allocatable, private :: asmp_p(:,:)
      real(kind = kreal), allocatable, private :: symn_t(:,:)
      real(kind = kreal), allocatable, private :: symn_p(:,:)
!
      integer(kind = kint), private :: nscl_jk
      real(kind = kreal), allocatable, private :: scl_e(:,:)
!
      integer(kind = kint), private :: nscl_lk
      real(kind = kreal), allocatable, private :: symp(:,:)
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine alloc_vec_fleg_mat_test(nvector)
!
      integer(kind = kint), intent(in) ::nvector
!
!
      num_jl = maxdegree_rlm * lmax_block_rtm
      allocate(Pvw_le(num_jl,np_smp))
      allocate(dPvw_le(num_jl,np_smp))
      allocate(Pgvw_le(num_jl,np_smp))
!
      nvec_jk = lmax_block_rtm * maxidx_rlm_smp(1)*nvector
      allocate(pol_e(nvec_jk,np_smp))
      allocate(dpl_e(nvec_jk,np_smp))
      allocate(tor_e(nvec_jk,np_smp))
      allocate(dpl_o(nvec_jk,np_smp))
      allocate(tor_o(nvec_jk,np_smp))
!
      nvec_lk = lmax_block_rtm * maxidx_rlm_smp(1)*nvector
      allocate(symp_r(nvec_lk,np_smp))
      allocate(symn_t(nvec_lk,np_smp))
      allocate(symn_p(nvec_lk,np_smp))
!
      allocate(asmp_t(nvec_lk,np_smp))
      allocate(asmp_p(nvec_lk,np_smp))
!
      end subroutine alloc_vec_fleg_mat_test
!
! -----------------------------------------------------------------------
!
      subroutine alloc_scl_fleg_mat_test(nscalar)
!
      integer(kind = kint), intent(in) :: nscalar
!
!
      num_jl = maxdegree_rlm * lmax_block_rtm
      allocate(Pws_le(num_jl,np_smp))
!
      nscl_jk = lmax_block_rtm * maxidx_rlm_smp(1)*nscalar
      allocate(scl_e(nscl_jk,np_smp))
!
      nscl_lk = lmax_block_rtm * maxidx_rlm_smp(1)*nscalar
      allocate(symp(nscl_lk,np_smp))
!
      end subroutine alloc_scl_fleg_mat_test
!
! -----------------------------------------------------------------------
!
      subroutine dealloc_vec_fleg_mat_test
!
!
      deallocate(Pvw_le, dPvw_le, Pgvw_le)
!
      deallocate(pol_e, dpl_e, tor_e, dpl_o, tor_o)
!
      deallocate(symp_r, symn_t, symn_p)
      deallocate(asmp_t, asmp_p)
!
      end subroutine dealloc_vec_fleg_mat_test
!
! -----------------------------------------------------------------------
!
      subroutine dealloc_scl_fleg_mat_test
!
!
      deallocate(Pws_le)
      deallocate(scl_e)
      deallocate(symp)
!
      end subroutine dealloc_scl_fleg_mat_test
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine legendre_f_trans_vector_test(ncomp, nvector,           &
     &          vr_rtm_spin, sp_rlm_spin)
!
      integer(kind = kint), intent(in) :: ncomp, nvector
      real(kind = kreal), intent(in)                                    &
     &      :: vr_rtm_spin(nidx_rtm(2),nidx_rtm(1)*ncomp,nidx_rtm(3))
      real(kind = kreal), intent(inout)                                 &
     &      :: sp_rlm_spin(nidx_rlm(2),nidx_rtm(1)*ncomp)
!
      integer(kind = kint) :: ip, nb_nri, kr_nd, kk, k_rlm
      integer(kind = kint) :: l_rtm, ls_rtm, i_lk, i_jl, i_jk
      integer(kind = kint) :: nd, ll, mp_rlm, mn_rlm, j_rlm, jj
      integer(kind = kint) :: kst(np_smp), nkr(np_smp)
      integer(kind = kint) :: jst(np_smp), nj_rlm(np_smp)
      real(kind = kreal) :: r2_1d_rlm_r
!
!
      call alloc_vec_fleg_mat_test(nvector)
!
      nb_nri = nvector*nidx_rlm(1)
!$omp parallel do schedule(static)                                      &
!$omp             private(ip,kr_nd,ll,l_rtm,ls_rtm,jj,j_rlm,nd,         &
!$omp&                    k_rlm,kk,mp_rlm,mn_rlm,i_lk,i_jl,i_jk)
      do ip = 1, np_smp
        kst(ip) = nvector*idx_rlm_smp_stack(ip-1,1)
        nkr(ip) = nvector                                               &
     &       * (idx_rtm_smp_stack(ip,  1) - idx_rtm_smp_stack(ip-1,1))
!
        do mp_rlm = 1, nidx_rtm(3)
          mn_rlm = nidx_rtm(3) - mp_rlm + 1
          jst(ip) = lstack_rlm(mp_rlm-1)
          nj_rlm(ip) = lstack_rlm(mp_rlm) - lstack_rlm(mp_rlm-1)
!    even l-m
          do l_rtm = 1, nidx_rtm(2)
            do jj = 1, nj_rlm(ip)
              j_rlm = jj + jst(ip)
              i_jl = jj + (l_rtm-1) * nj_rlm(ip)
              Pvw_le(i_jl,ip) = P_rtm(l_rtm,j_rlm)                      &
     &               * g_sph_rlm(j_rlm,7)* weight_rtm(l_rtm)
              dPvw_le(i_jl,ip) = dPdt_rtm(l_rtm,j_rlm)                  &
     &               * g_sph_rlm(j_rlm,7)* weight_rtm(l_rtm)
              Pgvw_le(i_jl,ip) = P_rtm(l_rtm,j_rlm)                     &
     &               * dble(idx_gl_1d_rlm_j(j_rlm,3))                   &
     &                * asin_theta_1d_rtm(l_rtm)                        &
     &                * g_sph_rlm(j_rlm,7)* weight_rtm(l_rtm)
            end do
          end do
!
          do kk = 1, nkr(ip)
            kr_nd = kk + kst(ip)
            do l_rtm = 1, nidx_rtm(2)
              ls_rtm = nidx_rtm(2) - l_rtm + 1
              i_lk = l_rtm + (kk-1) * nidx_rtm(2)
!
              symp_r(i_lk,ip)                                           &
     &               =  vr_rtm_spin(l_rtm, kr_nd,         mp_rlm)
!
              asmp_t(i_lk,ip)                                           &
     &               =  vr_rtm_spin(l_rtm, kr_nd+nb_nri,  mp_rlm)
              asmp_p(i_lk,ip)                                           &
     &               =  vr_rtm_spin(l_rtm, kr_nd+2*nb_nri,mp_rlm)
!
              symn_t(i_lk,ip)                                           &
     &               =  vr_rtm_spin(l_rtm, kr_nd+nb_nri,  mn_rlm)
              symn_p(i_lk,ip)                                           &
     &               =  vr_rtm_spin(l_rtm, kr_nd+2*nb_nri,mn_rlm)
            end do
          end do
!
          call matmul_fwd_leg_trans(nj_rlm(ip), nkr(ip), nidx_rtm(2),   &
     &        Pvw_le(1,ip), symp_r(1,ip), pol_e(1,ip))
          call matmul_fwd_leg_trans(nj_rlm(ip), nkr(ip), nidx_rtm(2),   &
     &        dPvw_le(1,ip), asmp_t(1,ip), dpl_o(1,ip))
          call matmul_fwd_leg_trans(nj_rlm(ip), nkr(ip), nidx_rtm(2),   &
     &        Pgvw_le(1,ip), symn_p(1,ip), dpl_e(1,ip))
          call matmul_fwd_leg_trans(nj_rlm(ip), nkr(ip), nidx_rtm(2),   &
     &        Pgvw_le(1,ip), symn_t(1,ip), tor_e(1,ip))
          call matmul_fwd_leg_trans(nj_rlm(ip), nkr(ip), nidx_rtm(2),   &
     &        dPvw_le(1,ip), asmp_p(1,ip), tor_o(1,ip))
!
          do kk = 1, nkr(ip)
            kr_nd = kk + kst(ip)
            k_rlm = 1 + mod((kr_nd-1),nidx_rlm(1))
            do jj = 1, nj_rlm(ip)
              j_rlm = jj + jst(ip)
              i_jk = jj + (kk-1) * nj_rlm(ip)
              sp_rlm_spin(j_rlm,kr_nd         )                         &
     &               = sp_rlm_spin(j_rlm,kr_nd         )                &
     &                + pol_e(i_jk,ip)
              sp_rlm_spin(j_rlm,kr_nd+nb_nri  )                         &
     &               = sp_rlm_spin(j_rlm,kr_nd+nb_nri  )                &
     &                - dpl_e(i_jk,ip) + dpl_o(i_jk,ip)
              sp_rlm_spin(j_rlm,kr_nd+2*nb_nri)                         &
     &               = sp_rlm_spin(j_rlm,kr_nd+2*nb_nri)                &
     &                - tor_e(i_jk,ip) - tor_o(i_jk,ip)
            end do
          end do
!
        end do
      end do
!$omp end parallel do
!
!$omp parallel do schedule(static)                                      &
!$omp&            private(ip,kk,kr_nd,j_rlm,k_rlm,r2_1d_rlm_r)
      do ip = 1, np_smp
        kst(ip) = nvector*idx_rtm_smp_stack(ip-1,1)
        nkr(ip) = nvector                                               &
     &       * (idx_rtm_smp_stack(ip,  1) - idx_rtm_smp_stack(ip-1,1))
        do kk = 1, nkr(ip)
          kr_nd = kk + kst(ip)
          k_rlm = 1 + mod((kr_nd-1),nidx_rlm(1))
          r2_1d_rlm_r = radius_1d_rlm_r(k_rlm) * radius_1d_rlm_r(k_rlm)
          do j_rlm = 1, nidx_rlm(2)
!
            sp_rlm_spin(j_rlm,kr_nd         )                           &
     &        = sp_rlm_spin(j_rlm,kr_nd         ) * r2_1d_rlm_r
            sp_rlm_spin(j_rlm,kr_nd+nb_nri  )                           &
     &        = sp_rlm_spin(j_rlm,kr_nd+nb_nri  )                       &
     &         * radius_1d_rlm_r(k_rlm)
            sp_rlm_spin(j_rlm,kr_nd+2*nb_nri)                           &
     &        = sp_rlm_spin(j_rlm,kr_nd+2*nb_nri)                       &
     &         * radius_1d_rlm_r(k_rlm)
            end do
        end do
      end do
!$omp end parallel do
!
      call dealloc_vec_fleg_mat_test
!
      end subroutine legendre_f_trans_vector_test
!
! -----------------------------------------------------------------------
!
      subroutine legendre_f_trans_scalar_test(ncomp, nvector, nscalar,  &
     &          vr_rtm_spin, sp_rlm_spin)
!
      integer(kind = kint), intent(in) :: ncomp, nvector, nscalar
      real(kind = kreal), intent(in)                                    &
     &      :: vr_rtm_spin(nidx_rtm(2),nidx_rtm(1)*ncomp,nidx_rtm(3))
      real(kind = kreal), intent(inout)                                 &
     &      :: sp_rlm_spin(nidx_rlm(2),nidx_rtm(1)*ncomp)
!
      integer(kind = kint) :: ip, kr_nd, kk
      integer(kind = kint) :: mp_rlm, j_rlm, jj
      integer(kind = kint) :: ll, l_rtm, ls_rtm, i_lk, i_jl, i_jk
      integer(kind = kint) :: kst(np_smp), nkr(np_smp)
      integer(kind = kint) :: jst(np_smp), nj_rlm(np_smp)
!
!
      call alloc_scl_fleg_mat_test(nscalar)
!
!$omp parallel do schedule(static)                                      &
!$omp&            private(ip,kr_nd,kk,jj,j_rlm,mp_rlm,                  &
!$omp&                    ll,l_rtm,ls_rtm,i_lk,i_jl,i_jk)
      do ip = 1, np_smp
        kst(ip) = 3*nvector*nidx_rlm(1)                                 &
     &           + nscalar*idx_rlm_smp_stack(ip-1,1)
        nkr(ip) = nscalar                                               &
     &       * (idx_rtm_smp_stack(ip,  1) - idx_rtm_smp_stack(ip-1,1))
        do mp_rlm = 1, nidx_rtm(3)
          jst(ip) = lstack_rlm(mp_rlm-1)
          nj_rlm(ip) = lstack_rlm(mp_rlm) - lstack_rlm(mp_rlm-1)
!
!    even l-m
          do l_rtm = 1, nidx_rtm(2)
            do jj = 1, nj_rlm(ip)
              j_rlm = jj + jst(ip)
              i_jl = jj + (l_rtm-1) * nj_rlm(ip)
              Pws_le(i_jl,ip) = P_rtm(l_rtm,j_rlm)                      &
     &                         * g_sph_rlm(j_rlm,6) * weight_rtm(l_rtm)
            end do
          end do
!
          do kk = 1, nkr(ip)
            kr_nd = kk + kst(ip)
            do l_rtm = 1, nidx_rtm(2)
              ls_rtm = nidx_rtm(2) - l_rtm + 1
              i_lk = l_rtm + (kk-1) * nidx_rtm(2)
              symp(i_lk,ip) =  vr_rtm_spin(l_rtm, kr_nd,mp_rlm)
            end do
          end do
!
          call matmul_fwd_leg_trans(nj_rlm(ip), nkr(ip), nidx_rtm(2),   &
     &        Pws_le(1,ip), symp(1,ip), scl_e(1,ip))
!
          do kk = 1, nkr(ip)
              kr_nd = kk + kst(ip)
            do jj = 1, nj_rlm(ip)
              j_rlm = jj + jst(ip)
              i_jk = jj + (kk-1) * nj_rlm(ip)
              sp_rlm_spin(j_rlm,kr_nd)                                  &
     &            = sp_rlm_spin(j_rlm,kr_nd) + scl_e(i_jk,ip)
            end do
          end do
!
        end do
      end do
!$omp end parallel do
!
!
      call dealloc_scl_fleg_mat_test
!
      end subroutine legendre_f_trans_scalar_test
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine matmat_fwd_leg_trans(n_jk, nkr, nl_rtm,                &
     &          P_jl, V_lk, S_jk)
!
      integer(kind = kint), intent(in) :: n_jk, nkr, nl_rtm
      real(kind = kreal), intent(in) :: P_jl(n_jk,nl_rtm)
      real(kind = kreal), intent(in) :: V_lk(nl_rtm,nkr)
!
      real(kind = kreal), intent(inout) :: S_jk(n_jk,nkr)
!
      integer(kind = kint) :: jj, kk, ll
      real(kind = kreal) :: s
!
!
      do kk = 1, nkr
        do jj = 1, n_jk
          s = 0.0d0
          do ll = 1, nl_rtm
            s = s + P_jl(jj,ll) * V_lk(ll,kk)
          end do
          S_jk(jj,kk) = s
        end do
      end do
!
      end subroutine matmat_fwd_leg_trans
!
! ----------------------------------------------------------------------
!
      subroutine matmul_fwd_leg_trans(n_jk, nkr, nl_rtm,                &
     &          P_jl, V_lk, S_jk)
!
      integer(kind = kint), intent(in) :: n_jk, nkr, nl_rtm
      real(kind = kreal), intent(in) :: P_jl(n_jk,nl_rtm)
      real(kind = kreal), intent(in) :: V_lk(nl_rtm,nkr)
!
      real(kind = kreal), intent(inout) :: S_jk(n_jk,nkr)
!
!
      S_jk = matmul(P_jl,V_lk)
!
      end subroutine matmul_fwd_leg_trans
!
! ----------------------------------------------------------------------
!
!      subroutine dgemm_fwd_leg_trans(n_jk, nkr, nl_rtm,                &
!     &          P_jl, V_lk, S_jk)
!
!      integer(kind = kint), intent(in) :: n_jk, nkr, nl_rtm
!      real(kind = kreal), intent(in) :: P_jl(n_jk,nl_rtm)
!      real(kind = kreal), intent(in) :: V_lk(nl_rtm,nkr)
!
!      real(kind = kreal), intent(inout) :: S_jk(n_jk,nkr)
!
!
!      call DGEMM('N', 'N', n_jk, nkr, nl_rtm, one,                     &
!     &    P_jl, n_jk, V_lk, nl_rtm, zero, S_jk, n_jk)
!
!      end subroutine dgemm_fwd_leg_trans
!
! ----------------------------------------------------------------------
!
      end module legendre_fwd_trans_testloop
