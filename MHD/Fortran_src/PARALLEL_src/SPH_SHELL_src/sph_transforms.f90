!>@file   sph_transforms.f90
!!@brief  module sph_transforms
!!
!!@author H. Matsui
!!@date Programmed in Aug., 2007
!!@n    Modified in Apr. 2013
!
!>@brief Spherical harmonics transform for vector
!!       and gradient of scalar
!!
!!@verbatim
!!      subroutine sph_backward_transforms                              &
!!     &         (ncomp_trans, nvector, nscalar, n_WS, n_WR, WS, WR)
!!      subroutine sph_forward_transforms                               &
!!     &         (ncomp_trans, nvector, nscalar, n_WS, n_WR, WS, WR)
!!
!!   input /outpt arrays for single field
!!
!!      radial component:      vr_rtp(i_rtp,1)
!!      elevetional component: vr_rtp(i_rtp,2)
!!      azimuthal component:   vr_rtp(i_rtp,3)
!!
!!      Poloidal component:          sp_rj(3*i_rj-2)
!!      diff. of Poloidal component: sp_rj(3*i_rj-1)
!!      Toroidal component:          sp_rj(3*i_rj  )
!!@endverbatim
!!
!!@param ncomp_trans Number of components for transform
!
      module sph_transforms
!
      use m_precision
!
      use calypso_mpi
      use m_work_time
      use m_machine_parameter
      use sph_FFT_selector
      use legendre_transform_select
      use spherical_SRs_N
!
      implicit none
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine sph_backward_transforms                                &
     &         (ncomp_trans, nvector, nscalar, n_WS, n_WR, WS, WR)
!
      integer(kind = kint), intent(in) :: ncomp_trans
      integer(kind = kint), intent(in) :: nvector, nscalar
      integer(kind = kint), intent(in) :: n_WS, n_WR
      real(kind = kreal), intent(inout) :: WS(n_WS), WR(n_WR)
!
!
      START_SRtime= MPI_WTIME()
      call start_eleps_time(18)
      call calypso_sph_comm_rj_2_rlm_N(ncomp_trans)
      call end_eleps_time(18)
      SendRecvtime = MPI_WTIME() - START_SRtime + SendRecvtime
!
      call start_eleps_time(22)
      if(iflag_debug .gt. 0) write(*,*) 'sel_backward_legendre_trans'
      call sel_backward_legendre_trans                                  &
     &   (ncomp_trans, nvector, nscalar, n_WR, n_WS, WR, WS)
      call end_eleps_time(22)
!
      START_SRtime= MPI_WTIME()
      call start_eleps_time(19)
      call calypso_sph_comm_rtm_2_rtp_N(ncomp_trans)
      call end_eleps_time(19)
      SendRecvtime = MPI_WTIME() - START_SRtime + SendRecvtime
!
!      call check_vr_rtp(my_rank, ncomp_trans)
!
      call start_eleps_time(24)
      call back_FFT_select_from_recv(ncomp_trans, n_WR, WR)
      call end_eleps_time(24)
!
      call finish_send_recv_rtm_2_rtp
!
!      call check_vr_rtp(my_rank, ncomp_trans)
!
      end subroutine sph_backward_transforms
!
! -----------------------------------------------------------------------
!
      subroutine sph_forward_transforms                                 &
     &         (ncomp_trans, nvector, nscalar, n_WS, n_WR, WS, WR)
!
      integer(kind = kint), intent(in) :: ncomp_trans
      integer(kind = kint), intent(in) :: nvector, nscalar
      integer(kind = kint), intent(in) :: n_WS, n_WR
      real(kind = kreal), intent(inout) :: WS(n_WS), WR(n_WR)
!
!      call check_vr_rtp(my_rank, ncomp_trans)
      call start_eleps_time(24)
      call fwd_FFT_select_to_send(ncomp_trans, n_WS, WS)
      call end_eleps_time(24)
!      call check_vr_rtp(my_rank, ncomp_trans)
!
      START_SRtime= MPI_WTIME()
      call start_eleps_time(20)
      call calypso_sph_comm_rtp_2_rtm_N(ncomp_trans)
      call end_eleps_time(20)
      SendRecvtime = MPI_WTIME() - START_SRtime + SendRecvtime
!
      call start_eleps_time(23)
      if(iflag_debug .gt. 0) write(*,*) 'sel_forward_legendre_trans'
      call sel_forward_legendre_trans                                   &
     &   (ncomp_trans, nvector, nscalar, n_WR, n_WS, WR, WS)
      call end_eleps_time(23)
!
      START_SRtime= MPI_WTIME()
      call start_eleps_time(21)
      call calypso_sph_comm_rlm_2_rj_N(ncomp_trans)
      call finish_send_recv_rlm_2_rj
      call end_eleps_time(21)
      SendRecvtime = MPI_WTIME() - START_SRtime + SendRecvtime
!
      end subroutine sph_forward_transforms
!
! -----------------------------------------------------------------------
!
      end module sph_transforms
