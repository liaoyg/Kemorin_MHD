!>@file   legendre_transform_testloop.f90
!!@brief  module legendre_transform_testloop
!!
!!@author H. Matsui
!!@date Programmed in Aug., 2007
!!@n    Modified in Apr. 2013
!
!>@brief  Legendre transforms
!!       (longest loop version)
!!
!!
!!@verbatim
!!      subroutine leg_backward_trans_test                              &
!!     &         (ncomp, nvector, nscalar, n_WR, n_WS, WR, WS)
!!        Input:  sp_rlm   (Order: poloidal,diff_poloidal,toroidal)
!!        Output: vr_rtm   (Order: radius,theta,phi)
!!
!!    Forward transforms
!!      subroutine leg_forward_trans_test                               &
!!     &         (ncomp, nvector, nscalar, n_WR, n_WS, WR, WS)
!!        Input:  vr_rtm   (Order: radius,theta,phi)
!!        Output: sp_rlm   (Order: poloidal,diff_poloidal,toroidal)
!!@endverbatim
!!
!!@param   ncomp    Total number of components for spherical transform
!!@param   nvector  Number of vector for spherical transform
!!@param   nscalar  Number of scalar (including tensor components)
!!                  for spherical transform
!
      module legendre_transform_testloop
!
      use m_precision
!
      implicit none
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine leg_backward_trans_test                                &
     &         (ncomp, nvector, nscalar, n_WR, n_WS, WR, WS)
!
      use m_work_4_sph_trans_spin
      use legendre_bwd_trans_testloop
      use ordering_schmidt_trans_spin
      use ordering_schmidt_trans_krin
      use merge_polidal_toroidal_v
      use spherical_SRs_N
!
      integer(kind = kint), intent(in) :: ncomp, nvector, nscalar
      integer(kind = kint), intent(in) :: n_WR, n_WS
      real (kind=kreal), intent(inout):: WR(n_WR)
      real (kind=kreal), intent(inout):: WS(n_WS)
!
!
      call calypso_rlm_from_recv_N(ncomp, n_WR, WR, sp_rlm)
      call order_b_trans_fields_spin(ncomp, nvector, nscalar,           &
     &    sp_rlm(1), sp_rlm_wk(1))
      call clear_bwd_legendre_work(ncomp)
      call clear_bwd_legendre_trans(ncomp)
      if(nscalar .gt. 0) then
        call legendre_b_trans_scalar_test(ncomp, nvector, nscalar,      &
     &      sp_rlm_wk(1), vr_rtm_wk(1))
      end if
      if(nvector .gt. 0) then
        call legendre_b_trans_vector_test(ncomp, nvector,               &
     &      sp_rlm_wk(1), vr_rtm_wk(1))
      end if
!
      call finish_send_recv_rj_2_rlm
      call back_b_trans_fields_spin(ncomp, nvector, nscalar,            &
     &    vr_rtm_wk(1), vr_rtm(1))
      call calypso_rtm_to_send_N(ncomp, n_WS, vr_rtm, WS)
!
      end subroutine leg_backward_trans_test
!
! -----------------------------------------------------------------------
!
      subroutine leg_forward_trans_test                                 &
     &         (ncomp, nvector, nscalar, n_WR, n_WS, WR, WS)
!
      use m_work_4_sph_trans_spin
      use legendre_fwd_trans_testloop
      use ordering_schmidt_trans_spin
      use spherical_SRs_N
!
      integer(kind = kint), intent(in) :: ncomp, nvector, nscalar
      integer(kind = kint), intent(in) :: n_WR, n_WS
      real (kind=kreal), intent(inout):: WR(n_WR)
      real (kind=kreal), intent(inout):: WS(n_WS)
!
!
      call calypso_rtm_from_recv_N(ncomp, n_WR, WR, vr_rtm)
      call order_f_trans_fields_spin(ncomp, nvector, nscalar,           &
     &    vr_rtm(1), vr_rtm_wk(1))
      call clear_fwd_legendre_work(ncomp)
!
      if(nvector .gt. 0) then
        call legendre_f_trans_vector_test(ncomp, nvector,               &
     &      vr_rtm_wk(1), sp_rlm_wk(1))
      end if
      if(nscalar .gt. 0) then
        call legendre_f_trans_scalar_test(ncomp, nvector, nscalar,      &
     &      vr_rtm_wk(1), sp_rlm_wk(1))
      end if
!
      call finish_send_recv_rtp_2_rtm
      call back_f_trans_fields_spin(ncomp, nvector, nscalar,            &
     &    sp_rlm_wk(1), sp_rlm(1))
      call calypso_rlm_to_send_N(ncomp, n_WS, sp_rlm, WS)
!
      end subroutine leg_forward_trans_test
!
! -----------------------------------------------------------------------
!
      end module legendre_transform_testloop

