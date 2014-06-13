!>@file   spherical_SRs_N.f90
!!@brief  module spherical_SRs_N
!!
!!@author H. Matsui
!!@date Programmed in Aug., 2007
!
!>@brief  Arbitrary components data communications 
!!@n      for spherical harmonics transform
!!
!!@verbatim
!!      subroutine init_sph_send_recv_N(NB, X_rtp, X_rtm, X_rlm, X_rj)
!!
!!      subroutine send_recv_rtp_2_rtm_N(NB, X_rtp, X_rtm)
!!      subroutine send_recv_rtm_2_rtp_N(NB, X_rtm, X_rtp)
!!      subroutine send_recv_rj_2_rlm_N(NB, X_rj, X_rlm)
!!      subroutine send_recv_rlm_2_rj_N(NB, X_rlm, X_rj)
!!
!!      subroutine finish_send_recv_rtp_2_rtm
!!      subroutine finish_send_recv_rtm_2_rtp
!!      subroutine finish_send_recv_rj_2_rlm
!!      subroutine finish_send_recv_rlm_2_rj
!!@endverbatim
!!
!!
!!@n @param  NB    Number of components for communication
!!@n @param  X_rtp(NB*nnod_rtp)  @f$ f(r,\theta,\phi) @f$
!!@n               (Order, X_rtp(i_comp,inod))
!!@n @param  X_rtm(NB*nnod_rtm)  @f$ f(r,\theta,m) @f$
!!@n               (Order, X_rtm(i_comp,inod))
!!@n @param  X_rlm(NB*nnod_rlm)  @f$ f(r,l,m) @f$
!!@n               (Order, X_rlm(i_comp,inod))
!!@n @param  X_rj(NB*nnod_rj)    @f$ f(r,j) @f$
!!@n               (Order, X_rj(i_comp,inod))
!
!
      module spherical_SRs_N
!
      use m_precision
!
      use m_constants
      use m_spheric_parameter
      use m_sph_trans_comm_table
!
      use select_calypso_SR
!
      implicit none
!
!>      Data communication mode for arbitrary size data
      integer(kind = kint) :: iflag_sph_SRN = iflag_import_item
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine init_sph_send_recv_N(NB, X_rtp, X_rtm, X_rlm, X_rj)
!
      use calypso_mpi
!
      use m_spheric_parameter
      use m_sph_trans_comm_table
!
      use m_solver_SR
!
      integer (kind=kint), intent(in) :: NB
      real (kind=kreal), intent(inout) :: X_rtp(NB*nnod_rtp)
      real (kind=kreal), intent(inout) :: X_rtm(NB*nnod_rtm)
      real (kind=kreal), intent(inout) :: X_rlm(NB*nnod_rlm)
      real (kind=kreal), intent(inout)::  X_rj(NB*nnod_rj)
!
      integer (kind=kint) :: nneib_max_send, nneib_max_recv
      integer (kind=kint) :: nnod_max_send,  nnod_max_recv
!
      real(kind = kreal) :: stime, etime
      real(kind = kreal) :: etime_item_import, etime_irev_import
!
!
      nneib_max_send = nneib_domain_rtp
      nneib_max_recv = nneib_domain_rtm
      nnod_max_send =  ntot_item_sr_rtp
      nnod_max_recv =  ntot_item_sr_rtm
!
      nneib_max_send = max(nneib_max_send,nneib_domain_rtm)
      nneib_max_recv = max(nneib_max_recv,nneib_domain_rtp)
      nnod_max_send =  max(nnod_max_send,ntot_item_sr_rtm)
      nnod_max_recv =  max(nnod_max_recv,ntot_item_sr_rtp)
!
      nneib_max_send = max(nneib_max_send,nneib_domain_rj)
      nneib_max_recv = max(nneib_max_recv,nneib_domain_rlm)
      nnod_max_send =  max(nnod_max_send,ntot_item_sr_rj)
      nnod_max_recv =  max(nnod_max_recv,ntot_item_sr_rlm)
!
      nneib_max_send = max(nneib_max_send,nneib_domain_rlm)
      nneib_max_recv = max(nneib_max_recv,nneib_domain_rj)
      nnod_max_send =  max(nnod_max_send,ntot_item_sr_rlm)
      nnod_max_recv =  max(nnod_max_recv,ntot_item_sr_rj)
!
      call resize_work_sph_SR(NB, nneib_max_send, nneib_max_recv,       &
     &    ntot_item_sr_rtp, ntot_item_sr_rtm)
!
!
!
      iflag_sph_SRN = iflag_import_item
      stime = MPI_WTIME()
      call send_recv_rtp_2_rtm_N(NB, X_rtp, X_rtm)
      call send_recv_rtm_2_rtp_N(NB, X_rtm, X_rtp)

      call send_recv_rj_2_rlm_N(NB, X_rj, X_rlm)
      call send_recv_rlm_2_rj_N(NB, X_rlm, X_rj)
!
      etime = MPI_WTIME() - stime
      call MPI_allREDUCE (etime, etime_item_import, ione,               &
     &    CALYPSO_REAL, MPI_SUM, CALYPSO_COMM, ierr_MPI)
      etime_item_import = etime_item_import / dble(nprocs)
!
      iflag_sph_SRN = iflag_import_rev
      stime = MPI_WTIME()
      call send_recv_rtp_2_rtm_N(NB, X_rtp, X_rtm)
      call send_recv_rtm_2_rtp_N(NB, X_rtm, X_rtp)
      call send_recv_rj_2_rlm_N(NB, X_rj, X_rlm)
      call send_recv_rlm_2_rj_N(NB, X_rlm, X_rj)
!
      etime = MPI_WTIME() - stime
      call MPI_allREDUCE (etime, etime_irev_import, ione,               &
     &    CALYPSO_REAL, MPI_SUM, CALYPSO_COMM, ierr_MPI)
      etime_irev_import = etime_irev_import / dble(nprocs)
!
      if(etime_irev_import .le. etime_item_import) then
        iflag_sph_SRN = iflag_import_rev
      end if
!
      if(my_rank .eq. 0) then
        write(*,*) 'Comm. mode for sph. trans.: ', iflag_sph_SRN
        write(*,*) '0: Time by reg. import list: ', etime_item_import
        write(*,*) '1: Time by rev. import list: ', etime_irev_import
      end if
!
      end subroutine init_sph_send_recv_N
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine send_recv_rtp_2_rtm_N(NB, X_rtp, X_rtm)
!
      use m_spheric_parameter
      use m_sph_trans_comm_table
!
      integer (kind=kint), intent(in) :: NB
      real (kind=kreal), intent(in)::    X_rtp(NB*nnod_rtp)
      real (kind=kreal), intent(inout):: X_rtm(NB*nnod_rtm)
!
!
      call sel_calypso_send_recv_N                                      &
     &             (iflag_sph_SRN, NB, nnod_rtp, nnod_rtm,              &
     &              nneib_domain_rtp, iflag_self_rtp,                   &
     &              id_domain_rtp, istack_sr_rtp, item_sr_rtp,          &
     &              nneib_domain_rtm, iflag_self_rtm,                   &
     &              id_domain_rtm, istack_sr_rtm, item_sr_rtm,          &
     &              irev_sr_rtm, X_rtp, X_rtm)
!
      end subroutine send_recv_rtp_2_rtm_N
!
! ----------------------------------------------------------------------
!
      subroutine send_recv_rtm_2_rtp_N(NB, X_rtm, X_rtp)
!
      use m_spheric_parameter
      use m_sph_trans_comm_table
!
      integer (kind=kint), intent(in) :: NB
      real (kind=kreal), intent(in)::    X_rtm(NB*nnod_rtm)
      real (kind=kreal), intent(inout):: X_rtp(NB*nnod_rtp)
!
!
      call sel_calypso_send_recv_N                                      &
     &             (iflag_sph_SRN, NB, nnod_rtm, nnod_rtp,              &
     &              nneib_domain_rtm, iflag_self_rtm,                   &
     &              id_domain_rtm, istack_sr_rtm, item_sr_rtm,          &
     &              nneib_domain_rtp, iflag_self_rtp,                   &
     &              id_domain_rtp, istack_sr_rtp, item_sr_rtp,          &
     &              irev_sr_rtp, X_rtm, X_rtp)
!
      end subroutine send_recv_rtm_2_rtp_N
!
! ----------------------------------------------------------------------
!
      subroutine send_recv_rj_2_rlm_N(NB, X_rj, X_rlm)
!
      use m_spheric_parameter
      use m_sph_trans_comm_table
!
      integer (kind=kint), intent(in) :: NB
      real (kind=kreal), intent(in)::    X_rj(NB*nnod_rj)
      real (kind=kreal), intent(inout):: X_rlm(NB*nnod_rlm)
!
!
      call sel_calypso_send_recv_N                                      &
     &             (iflag_sph_SRN, NB, nnod_rj, nnod_rlm,               &
     &              nneib_domain_rj, iflag_self_rj,                     &
     &              id_domain_rj, istack_sr_rj, item_sr_rj,             &
     &              nneib_domain_rlm, iflag_self_rlm,                   &
     &              id_domain_rlm, istack_sr_rlm, item_sr_rlm,          &
     &              irev_sr_rlm, X_rj, X_rlm)
!
      end subroutine send_recv_rj_2_rlm_N
!
! ----------------------------------------------------------------------
!
      subroutine send_recv_rlm_2_rj_N(NB, X_rlm, X_rj)
!
      use m_spheric_parameter
      use m_sph_trans_comm_table
!
      integer (kind=kint), intent(in) :: NB
      real (kind=kreal), intent(in)::    X_rlm(NB*nnod_rlm)
      real (kind=kreal), intent(inout):: X_rj(NB*nnod_rj)
!
!
      call sel_calypso_send_recv_N                                      &
     &             (iflag_sph_SRN, NB, nnod_rlm, nnod_rj,               &
     &              nneib_domain_rlm, iflag_self_rlm,                   &
     &              id_domain_rlm, istack_sr_rlm, item_sr_rlm,          &
     &              nneib_domain_rj, iflag_self_rj,                     &
     &              id_domain_rj, istack_sr_rj, item_sr_rj,             &
     &              irev_sr_rj, X_rlm, X_rj)
!
      end subroutine send_recv_rlm_2_rj_N
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine finish_send_recv_rtp_2_rtm
!
      use m_sph_trans_comm_table
!
      call finish_calypso_send_recv(nneib_domain_rtp, iflag_self_rtp)
!
      end subroutine finish_send_recv_rtp_2_rtm
!
! ----------------------------------------------------------------------
!
      subroutine finish_send_recv_rtm_2_rtp
!
      use m_sph_trans_comm_table
!
      call finish_calypso_send_recv(nneib_domain_rtm, iflag_self_rtm)
!
      end subroutine finish_send_recv_rtm_2_rtp
!
! ----------------------------------------------------------------------
!
      subroutine finish_send_recv_rj_2_rlm
!
      use m_sph_trans_comm_table
!
      call finish_calypso_send_recv(nneib_domain_rj, iflag_self_rj)
!
      end subroutine finish_send_recv_rj_2_rlm
!
! ----------------------------------------------------------------------
!
      subroutine finish_send_recv_rlm_2_rj
!
      use m_sph_trans_comm_table
!
      call finish_calypso_send_recv(nneib_domain_rlm, iflag_self_rlm)
!
      end subroutine finish_send_recv_rlm_2_rj
!
! ----------------------------------------------------------------------
!
      end module spherical_SRs_N
