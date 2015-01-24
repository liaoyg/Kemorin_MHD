!>@file   legendre_bwd_trans_testloop.f90
!!@brief  module legendre_bwd_trans_testloop
!!
!!@author H. Matsui
!!@date Programmed in Aug., 2007
!!@n    Modified in Apr. 2013
!
!>@brief  backward Legendre transform
!!       (Blocked loop version)
!!
!!@verbatim
!!      subroutine legendre_b_trans_vector_test(ncomp, nvector,         &
!!     &          irev_sr_rlm, n_WR, WR, WS)
!!        Input:  vr_rtm   (Order: radius,theta,phi)
!!        Output: sp_rlm   (Order: poloidal,diff_poloidal,toroidal)
!!      subroutine legendre_b_trans_scalar_test(ncomp, nvector, nscalar,&
!!     &          irev_sr_rlm, n_WR, WR, WS)
!!        Input:  vr_rtm
!!        Output: sp_rlm
!!@endverbatim
!!
!!@param   ncomp    Total number of components for spherical transform
!!@param   nvector  Number of vector for spherical transform
!!@param   nscalar  Number of scalar (including tensor components)
!!                  for spherical transform
!
      module legendre_bwd_trans_testloop
!
      use m_precision
!
      use m_constants
      use m_work_time
      use calypso_mpi
!
      use m_machine_parameter
      use m_spheric_parameter
      use m_spheric_param_smp
      use m_schmidt_poly_on_rtm
      use m_work_4_sph_trans
      use m_legendre_work_testlooop
      use matmul_for_legendre_trans
!
      implicit none
!
      real(kind = kreal), private :: st_elapsed
      real(kind = kreal), private :: elaps(4)
      integer, external :: omp_get_max_threads
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine legendre_b_trans_vector_test(ncomp, nvector, nscalar,  &
     &          irev_sr_rlm, irev_sr_rtm, n_WR, n_WS, WR, WS)
!
      use set_sp_rlm_for_leg_matmul
      use cal_vr_rtm_by_matmul
!
      integer(kind = kint), intent(in) :: ncomp, nvector, nscalar
      integer(kind = kint), intent(in) :: n_WR, n_WS
      integer(kind = kint), intent(in) :: irev_sr_rlm(nnod_rlm)
      integer(kind = kint), intent(in) :: irev_sr_rtm(nnod_rtm)
      real (kind=kreal), intent(inout):: WR(n_WR)
      real (kind=kreal), intent(inout):: WS(n_WS)
!
      integer(kind = kint) :: ip
      integer(kind = kint) :: nl_rtm, mp_rlm, mn_rlm
      integer(kind = kint) :: kstv(np_smp), nkrv(np_smp)
      integer(kind = kint) :: ksts(np_smp), nkrs(np_smp)
      integer(kind = kint) :: jst(np_smp), jst_h(np_smp)
      integer(kind = kint) :: n_jk_e(np_smp), n_jk_o(np_smp)
!
!
      elaps(1:4) = 0
      nl_rtm = (nidx_rtm(2) + 1)/2
!$omp parallel do schedule(static)                                      &
!$omp&            private(ip,mp_rlm,mn_rlm,st_elapsed)                  &
!$omp& reduction(+:elaps)
      do ip = 1, np_smp
        kstv(ip) = nvector*idx_rtm_smp_stack(ip-1,1)
        nkrv(ip) = nvector                                              &
     &       * (idx_rtm_smp_stack(ip,  1) - idx_rtm_smp_stack(ip-1,1))
        ksts(ip) = nscalar*idx_rtm_smp_stack(ip-1,1)
        nkrs(ip) = nscalar                                              &
     &       * (idx_rtm_smp_stack(ip,  1) - idx_rtm_smp_stack(ip-1,1))
!
        do mp_rlm = 1, nidx_rtm(3)
          mn_rlm = nidx_rtm(3) - mp_rlm + 1
          jst(ip) = lstack_rlm(mp_rlm-1)
          jst_h(ip) = lstack_even_rlm(mp_rlm) + 1
          n_jk_e(ip) = lstack_even_rlm(mp_rlm) - lstack_rlm(mp_rlm-1)
          n_jk_o(ip) = lstack_rlm(mp_rlm) - lstack_even_rlm(mp_rlm)
!
!          st_elapsed = MPI_WTIME()
          call set_sp_rlm_vector_sym_matmul                             &
     &       (kstv(ip), nkrv(ip), jst(ip), n_jk_e(ip), n_jk_o(ip),      &
     &        ncomp, irev_sr_rlm, n_WR, WR, nvec_jk, pol_e(1,ip),       &
     &        dpoldt_e(1,ip), dpoldp_e(1,ip), dtordt_e(1,ip),           &
     &        dtordp_e(1,ip), pol_o(1,ip), dpoldt_o(1,ip),              &
     &        dpoldp_o(1,ip), dtordt_o(1,ip), dtordp_o(1,ip))
          call set_sp_rlm_scalar_sym_matmul                             &
     &       (ksts(ip), nkrs(ip), jst(ip), n_jk_e(ip), n_jk_o(ip),      &
     &        ncomp, nvector, irev_sr_rlm, n_WR, WR,                    &
     &        nscl_jk, scl_e(1,ip), scl_o(1,ip) )
!          elaps(2) = MPI_WTIME() - st_elapsed + elaps(2)
!
!   even l-m
!          st_elapsed = MPI_WTIME()
         if(nkrv(ip) .gt. 0) then
          call matmul_bwd_leg_trans(nl_rtm, nkrv(ip), n_jk_e(ip),       &
     &        Ps_rtm(1,jst(ip)+1), pol_e(1,ip), symp_r(1,ip))
          call matmul_bwd_leg_trans(nl_rtm, nkrv(ip), n_jk_e(ip),       &
     &        dPsdt_rtm(1,jst(ip)+1), dpoldt_e(1,ip), asmp_t(1,ip))
          call matmul_bwd_leg_trans(nl_rtm, nkrv(ip), n_jk_e(ip),       &
     &        dPsdt_rtm(1,jst(ip)+1), dtordt_e(1,ip), asmp_p(1,ip))
!
          call matmul_bwd_leg_trans(nl_rtm, nkrv(ip), n_jk_e(ip),       &
     &        Ps_rtm(1,jst(ip)+1), dtordp_e(1,ip), symn_t(1,ip))
          call matmul_bwd_leg_trans(nl_rtm, nkrv(ip), n_jk_e(ip),       &
     &        Ps_rtm(1,jst(ip)+1), dpoldp_e(1,ip), symn_p(1,ip))
!   odd l-m
          call matmul_bwd_leg_trans(nl_rtm, nkrv(ip), n_jk_o(ip),       &
     &        Ps_rtm(1,jst_h(ip)), pol_o(1,ip), asmp_r(1,ip))
          call matmul_bwd_leg_trans(nl_rtm, nkrv(ip), n_jk_o(ip),       &
     &        dPsdt_rtm(1,jst_h(ip)), dpoldt_o(1,ip), symp_t(1,ip))
          call matmul_bwd_leg_trans(nl_rtm, nkrv(ip), n_jk_o(ip),       &
     &        dPsdt_rtm(1,jst_h(ip)), dtordt_o(1,ip), symp_p(1,ip))
!
          call matmul_bwd_leg_trans(nl_rtm, nkrv(ip), n_jk_o(ip),       &
     &        Ps_rtm(1,jst_h(ip)), dtordp_o(1,ip), asmn_t(1,ip))
          call matmul_bwd_leg_trans(nl_rtm, nkrv(ip), n_jk_o(ip),       &
     &        Ps_rtm(1,jst_h(ip)), dpoldp_o(1,ip), asmn_p(1,ip))
         end if
!
!   even l-m
         if(nkrs(ip) .gt. 0) then
          call matmul_bwd_leg_trans(nl_rtm, nkrs(ip), n_jk_e(ip),       &
     &        Ps_rtm(1,jst(ip)+1), scl_e(1,ip), symp(1,ip))
!   odd l-m
          call matmul_bwd_leg_trans(nl_rtm, nkrs(ip), n_jk_o(ip),       &
     &        Ps_rtm(1,jst_h(ip)), scl_o(1,ip), asmp(1,ip))
         end if
!          elaps(3) = MPI_WTIME() - st_elapsed + elaps(3)
!
!          st_elapsed = MPI_WTIME()
          call cal_vr_rtm_vector_sym_matmul                             &
     &       (kstv(ip), nkrv(ip), mp_rlm, mn_rlm, nl_rtm, nvec_lk,      &
     &        symp_r(1,ip), asmp_t(1,ip), asmp_p(1,ip), symn_t(1,ip),   &
     &        symn_p(1,ip), asmp_r(1,ip), symp_t(1,ip), symp_p(1,ip),   &
     &        asmn_t(1,ip), asmn_p(1,ip), ncomp, irev_sr_rtm, n_WS, WS)
          call cal_vr_rtm_scalar_sym_matmul(ksts(ip), nkrs(ip),         &
     &        mp_rlm, nl_rtm, nscl_lk, symp(1,ip), asmp(1,ip),          &
     &        ncomp, nvector, irev_sr_rtm, n_WS, WS)
!          elaps(4) = MPI_WTIME() - st_elapsed + elaps(4)
!
        end do
      end do
!$omp end parallel do
!
!      elapsed(41:44)                                                   &
!     &     = elaps(1:4)/ dble(omp_get_max_threads()) + elapsed(41:44)
!
      end subroutine legendre_b_trans_vector_test
!
! -----------------------------------------------------------------------
!
!      subroutine legendre_b_trans_scalar_test(ncomp, nvector, nscalar,  &
!     &          irev_sr_rlm, irev_sr_rtm, n_WR, n_WS, WR, WS)
!
!
!      end subroutine legendre_b_trans_scalar_test
!
! -----------------------------------------------------------------------
!
      end module legendre_bwd_trans_testloop
