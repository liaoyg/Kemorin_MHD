!>@file   legendre_fwd_trans_sym_spin.f90
!!@brief  module legendre_fwd_trans_sym_spin
!!
!!@author H. Matsui
!!@date Programmed in Aug., 2013
!
!>@brief  forward Legendre transform considering symmetry
!!
!!@verbatim
!!      subroutine leg_fwd_trans_vector_sym_spin(ncomp, nvector,        &
!!     &          irev_sr_rtm, irev_sr_rlm, n_WR, n_WS, WR, WS)
!!        Input:  vr_rtm   (Order: radius,theta,phi)
!!        Output: sp_rlm   (Order: poloidal,diff_poloidal,toroidal)
!!      subroutine leg_fwd_trans_scalar_sym_spin                        &
!!     &         (ncomp, nvector, nscalar, irev_sr_rtm, irev_sr_rlm,    &
!!     &          n_WR, n_WS, WR, WS)
!!        Input:  vr_rtm
!!        Output: sp_rlm
!!@endverbatim
!!
!!@param   ncomp    Total number of components for spherical transform
!!@param   nvector  Number of vector for spherical transform
!!@param   nscalar  Number of scalar (including tensor components)
!!                  for spherical transform
!
      module legendre_fwd_trans_sym_spin
!
      use m_precision
!
      use m_machine_parameter
      use m_spheric_parameter
      use m_spheric_param_smp
      use m_schmidt_poly_on_rtm
      use m_work_4_sph_trans
      use m_legendre_work_sym_matmul
!
      implicit none
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine leg_fwd_trans_vector_sym_spin(ncomp, nvector,          &
     &          irev_sr_rtm, irev_sr_rlm, n_WR, n_WS, WR, WS)
!
      use set_vr_rtm_for_leg_vecprod
      use cal_sp_rlm_by_vecprod
!
      integer(kind = kint), intent(in) :: ncomp, nvector
      integer(kind = kint), intent(in) :: n_WR, n_WS
      integer(kind = kint), intent(in) :: irev_sr_rtm(nnod_rtm)
      integer(kind = kint), intent(in) :: irev_sr_rlm(nnod_rlm)
      real (kind=kreal), intent(inout):: WR(n_WR)
      real (kind=kreal), intent(inout):: WS(n_WS)
!
!
      integer(kind = kint) :: ip, kst, ked, k_rlm, nd, ie_rlm, io_rlm
      integer(kind = kint) :: mp_rlm, mn_rlm, ie_send, io_send
      integer(kind = kint) :: jst, nj_rlm, jj, je_rlm, jo_rlm
      integer(kind = kint) :: nle_rtm, nlo_rtm, n_jk_e, n_jk_o
      real(kind = kreal) :: r1_1d_rlm_r, r2_1d_rlm_r, gme, gmo
!
!
      nle_rtm = (nidx_rtm(2) + 1)/2
      nlo_rtm = nidx_rtm(2) / 2
!$omp parallel do schedule(static)                                      &
!$omp             private(ip,kst,ked,jj,k_rlm,nd,je_rlm,jo_rlm,         &
!$omp&                    mp_rlm,mn_rlm,jst,nj_rlm,n_jk_e,n_jk_o,       &
!$omp&                    r1_1d_rlm_r,r2_1d_rlm_r,                      &
!$omp&                    ie_rlm,io_rlm,ie_send,io_send,gme,gmo)
      do ip = 1, np_smp
        kst = idx_rlm_smp_stack(ip-1,1) + 1
        ked = idx_rlm_smp_stack(ip,  1)
        do k_rlm = kst, ked
          r1_1d_rlm_r = sph_rlm1%radius_1d_rlm_r(k_rlm)
          r2_1d_rlm_r = r1_1d_rlm_r*r1_1d_rlm_r
          do nd = 1, nvector
!
            do mp_rlm = 1, nidx_rtm(3)
              mn_rlm = nidx_rtm(3) - mp_rlm + 1
              jst = lstack_rlm(mp_rlm-1)
              nj_rlm = lstack_rlm(mp_rlm) - lstack_rlm(mp_rlm-1)
              n_jk_e = (nj_rlm+1) / 2
              n_jk_o =  nj_rlm - n_jk_e
!    even l-m
!    odd  l-m
              do jj = 1, nj_rlm/2
                je_rlm = 2*jj + jst - 1
                jo_rlm = 2*jj + jst
                gme = dble(sph_rlm1%idx_gl_1d_rlm_j(je_rlm,3))
                gmo = dble(sph_rlm1%idx_gl_1d_rlm_j(jo_rlm,3))
!
                ie_rlm = 1 + (je_rlm-1) * sph_rlm1%istep_rlm(2)         &
     &                     + (k_rlm-1) *  sph_rlm1%istep_rlm(1)
                io_rlm = 1 + (jo_rlm-1) * sph_rlm1%istep_rlm(2)         &
     &                     + (k_rlm-1)  * sph_rlm1%istep_rlm(1)
                ie_send = 3*nd-2 + (irev_sr_rlm(ie_rlm) - 1) * ncomp
                io_send = 3*nd-2 + (irev_sr_rlm(io_rlm) - 1) * ncomp
!
                call set_vr_rtm_vector_symmetry                         &
     &             (nd, k_rlm, mp_rlm, mn_rlm,                          &
     &              asin_theta_1d_rtm(1), izero, nle_rtm, nlo_rtm,      &
     &              ncomp, irev_sr_rtm, n_WR, WR, symp_r(1,ip),         &
     &              asmp_t(1,ip), asmp_p(1,ip), symn_t(1,ip),           &
     &              symn_p(1,ip), asmp_r(1,ip), symp_t(1,ip),           &
     &              symp_p(1,ip), asmn_t(1,ip), asmn_p(1,ip))
!
                call cal_vector_sp_rlm_dotprod(nth_hemi_rtm,            &
     &              g_sph_rlm(je_rlm,7), gme, r1_1d_rlm_r, r2_1d_rlm_r, &
     &              Ps_rtm(1,jj+jst), dPsdt_rtm(1,jj+jst),              &
     &              symp_r(1,ip), asmp_t(1,ip), asmp_p(1,ip),           &
     &              symn_t(1,ip), symn_p(1,ip), WS(ie_send))
                call cal_vector_sp_rlm_dotprod(nth_hemi_rtm,            &
     &             g_sph_rlm(jo_rlm,7), gmo, r1_1d_rlm_r, r2_1d_rlm_r,  &
     &             Ps_rtm(1,jj+jst+n_jk_e), dPsdt_rtm(1,jj+jst+n_jk_e), &
     &             asmp_r(1,ip), symp_t(1,ip), symp_p(1,ip),            &
     &             asmn_t(1,ip), asmn_p(1,ip), WS(io_send))
              end do
!
!   the last even l-m
              do jj = nj_rlm/2+1, (nj_rlm+1)/2
                je_rlm = 2*jj + jst - 1
                gme = dble(sph_rlm1%idx_gl_1d_rlm_j(je_rlm,3))
!
                ie_rlm = 1 + (je_rlm-1) * sph_rlm1%istep_rlm(2)         &
     &                     + (k_rlm-1) *  sph_rlm1%istep_rlm(1)
                ie_send = 3*nd-2 + (irev_sr_rlm(ie_rlm) - 1) * ncomp
!
                call set_vr_rtm_vector_symmetry                         &
     &             (nd, k_rlm, mp_rlm, mn_rlm,                          &
     &              asin_theta_1d_rtm(1), izero, nle_rtm, nlo_rtm,      &
     &              ncomp, irev_sr_rtm, n_WR, WR, symp_r(1,ip),         &
     &              asmp_t(1,ip), asmp_p(1,ip), symn_t(1,ip),           &
     &              symn_p(1,ip), asmp_r(1,ip), symp_t(1,ip),           &
     &              symp_p(1,ip), asmn_t(1,ip), asmn_p(1,ip))
                call cal_vector_sp_rlm_dotprod(nth_hemi_rtm,            &
     &              g_sph_rlm(je_rlm,7), gme, r1_1d_rlm_r, r2_1d_rlm_r, &
     &              Ps_rtm(1,jj+jst), dPsdt_rtm(1,jj+jst),              &
     &              symp_r(1,ip), asmp_t(1,ip), asmp_p(1,ip),           &
     &              symn_t(1,ip), symn_p(1,ip), WS(ie_send))
              end do
!
            end do
          end do
        end do
      end do
!$omp end parallel do
!
      end subroutine leg_fwd_trans_vector_sym_spin
!
! -----------------------------------------------------------------------
!
      subroutine leg_fwd_trans_scalar_sym_spin                          &
     &         (ncomp, nvector, nscalar,  irev_sr_rtm, irev_sr_rlm,     &
     &          n_WR, n_WS, WR, WS)
!
      use set_vr_rtm_for_leg_vecprod
      use cal_sp_rlm_by_vecprod
!
      integer(kind = kint), intent(in) :: ncomp, nvector, nscalar
      integer(kind = kint), intent(in) :: n_WR, n_WS
      integer(kind = kint), intent(in) :: irev_sr_rtm(nnod_rtm)
      integer(kind = kint), intent(in) :: irev_sr_rlm(nnod_rlm)
      real (kind=kreal), intent(inout):: WR(n_WR)
      real (kind=kreal), intent(inout):: WS(n_WS)
!
!
      integer(kind = kint) :: ip, kst, ked, k_rlm, nd
      integer(kind = kint) :: nle_rtm, nlo_rtm, je_rlm, jo_rlm
      integer(kind = kint) :: ie_rlm, io_rlm, ie_send, io_send
      integer(kind = kint) :: mp_rlm, jst, nj_rlm, jj, n_jk_e, n_jk_o
!
!
      nle_rtm = (nidx_rtm(2) + 1)/2
      nlo_rtm = nidx_rtm(2) / 2
!$omp parallel do schedule(static)                                      &
!$omp&            private(ip,kst,ked,jj,k_rlm,mp_rlm,n_jk_e,n_jk_o,     &
!$omp&                    nd,jst,nj_rlm,ie_rlm,io_rlm,je_rlm,jo_rlm,    &
!$omp&                    ie_send,io_send)
      do ip = 1, np_smp
        kst = idx_rlm_smp_stack(ip-1,1) + 1
        ked = idx_rlm_smp_stack(ip,  1)
        do k_rlm = kst, ked
          do nd = 1, nscalar
            do mp_rlm = 1, nidx_rtm(3)
              jst = lstack_rlm(mp_rlm-1)
              nj_rlm = lstack_rlm(mp_rlm) - lstack_rlm(mp_rlm-1)
              n_jk_e = (nj_rlm+1) / 2
              n_jk_o =  nj_rlm - n_jk_e
!    even l-m
!    odd  l-m
              do jj = 1, nj_rlm/2
                je_rlm = 2*jj + jst - 1
                jo_rlm = 2*jj + jst
                ie_rlm = 1 + (je_rlm-1) * sph_rlm1%istep_rlm(2)         &
     &                     + (k_rlm-1) *  sph_rlm1%istep_rlm(1)
                io_rlm = 1 + (jo_rlm-1) * sph_rlm1%istep_rlm(2)         &
     &                     + (k_rlm-1) *  sph_rlm1%istep_rlm(1)
                ie_send = nd + 3*nvector                                &
     &                       + (irev_sr_rlm(ie_rlm) - 1) * ncomp
                io_send = nd + 3*nvector                                &
     &                       + (irev_sr_rlm(io_rlm) - 1) * ncomp
                call set_vr_rtm_scalar_symmetry(nd, k_rlm, mp_rlm,      &
     &              izero, nle_rtm, nlo_rtm, ncomp, nvector,            &
     &              irev_sr_rtm, n_WR, WR, symp(1,ip), asmp(1,ip))
!
                call cal_scalar_sp_rlm_dotprod                          &
     &             (nle_rtm, g_sph_rlm(je_rlm,6),                       &
     &              Ps_rtm(1,jj+jst), symp(1,ip), WS(ie_send))
!
                call cal_scalar_sp_rlm_dotprod                          &
     &             (nlo_rtm, g_sph_rlm(jo_rlm,6),                       &
     &              Ps_rtm(1,jj+jst+n_jk_e), asmp(1,ip), WS(io_send))
              end do
!
!   the last even l-m
              do jj = nj_rlm/2+1, (nj_rlm+1)/2
                je_rlm = 2*jj + jst - 1
                call set_vr_rtm_scalar_symmetry(nd, k_rlm, mp_rlm,      &
     &              izero, nle_rtm, nlo_rtm, ncomp, nvector,            &
     &              irev_sr_rtm, n_WR, WR, symp(1,ip), asmp(1,ip))
!
                ie_rlm = 1 + (je_rlm-1) * sph_rlm1%istep_rlm(2)         &
     &                     + (k_rlm-1) *  sph_rlm1%istep_rlm(1)
                ie_send = nd + 3*nvector                                &
     &                       + (irev_sr_rlm(ie_rlm) - 1) * ncomp
                call cal_scalar_sp_rlm_dotprod                          &
     &             (nle_rtm, g_sph_rlm(je_rlm,6),                       &
     &              Ps_rtm(1,jj+jst), symp(1,ip), WS(ie_send))
              end do
!
            end do
          end do
        end do
      end do
!$omp end parallel do
!
      end subroutine leg_fwd_trans_scalar_sym_spin
!
! -----------------------------------------------------------------------
!
      end module legendre_fwd_trans_sym_spin
