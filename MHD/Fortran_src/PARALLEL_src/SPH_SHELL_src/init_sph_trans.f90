!>@file   init_sph_trans.f90
!!@brief  module init_sph_trans
!!
!!@author H. Matsui
!!@date Programmed in Aug., 2007
!
!>@brief  Initialize spherical harmonics transform
!!
!!@verbatim
!!      subroutine copy_sph_trans_nums_from_rtp
!!      subroutine initialize_legendre_trans
!!      subroutine initialize_sph_trans
!!@endverbatim
!
      module init_sph_trans
!
      use m_precision
      use m_constants
!
      implicit none
!
      private :: set_blocks_4_leg_trans
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine initialize_sph_trans
!
      use init_FFT_4_sph
      use m_work_4_sph_trans
!
!
      call initialize_legendre_trans
      call init_fourier_transform_4_sph(ncomp_sph_trans)
!
      end subroutine initialize_sph_trans
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine initialize_legendre_trans
!
      use m_spheric_parameter
      use m_schmidt_poly_on_rtm
      use m_work_4_sph_trans
      use m_FFT_selector
      use schmidt_poly_on_rtm_grid
      use set_legendre_matrices
      use set_params_sph_trans
!
!
      call allocate_work_4_sph_trans
!
      call radial_4_sph_trans(sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1)
      call set_mdx_rlm_rtm(l_truncation, nidx_rtm, nidx_rlm,            &
     &    sph_rtm1%idx_gl_1d_rtm_m, sph_rlm1%idx_gl_1d_rlm_j)
!
      call s_cal_schmidt_poly_rtm                                       &
     &   (l_truncation, sph_rj1, sph_rtm1, sph_rlm1)
!
      call set_sin_theta_rtm(sph_rtm1%nidx_rtm(2))
      call set_sin_theta_rtp                                            &
     &   (sph_rtp1%nidx_rtp(2), sph_rtp1%idx_gl_1d_rtp_t)
!
      call allocate_trans_schmidt_rtm                                   &
     &   (sph_rtm1%nidx_rtm(2), sph_rlm1%nidx_rlm(2))
      call set_trans_legendre_rtm
!
      call allocate_hemi_schmidt_rtm                                    &
     &   (sph_rtm1%nidx_rtm(2), sph_rlm1%nidx_rlm(2))
      call set_legendre_hemispher_rtm
!
      call set_blocks_4_leg_trans
!
      end subroutine initialize_legendre_trans
!
! -----------------------------------------------------------------------
!
      subroutine set_blocks_4_leg_trans
!
      use calypso_mpi
      use m_machine_parameter
      use m_sph_communicators
      use m_sph_trans_comm_table
      use m_spheric_parameter
      use m_work_4_sph_trans
      use init_spherical_SRs
      use cal_minmax_and_stacks
      use legendre_transform_select
!
!
      if(nvector_legendre .le. 0                                        &
     &     .or. nvector_legendre .gt. nidx_rtm(2)) then
        nblock_l_rtm =  1
      else
        nblock_l_rtm =  nidx_rtm(2) / nvector_legendre
      end if
      if(nvector_legendre .le. 0                                        &
     &     .or. nvector_legendre .gt. nidx_rlm(2)) then
        nblock_j_rlm =  1
      else
        nblock_j_rlm =  nidx_rlm(2) / nvector_legendre
      end if
!
      call allocate_l_rtm_block
      call count_number_4_smp(nblock_l_rtm, ione, nidx_rtm(2),          &
     &    lstack_block_rtm, lmax_block_rtm)
      call count_number_4_smp(nblock_j_rlm, ione, nidx_rlm(2),          &
     &    jstack_block_rlm, jmax_block_rlm)
!
!
      call split_rtp_comms(comm_rtp1%nneib_domain, comm_rtp1%id_domain, &
     &    nneib_domain_rj) 
      call init_sph_send_recv_N(ncomp_sph_trans)
!
!      if(iflag_sph_commN .eq. iflag_alltoall) then
!        call set_rev_all2all_import_tbl(nnod_rtp, nmax_sr_rtp,         &
!     &      comm_rtp1%nneib_domain, comm_rtp1%istack_sr,               &
!     &      comm_rtp1%item_sr, comm_rtp1%irev_sr)
!        call set_rev_all2all_import_tbl(nnod_rtm, nmax_sr_rtp,         &
!     &      comm_rtm1%nneib_domain, comm_rtm1%istack_sr,               &
!     &      comm_rtm1%item_sr, comm_rtm1%irev_sr)
!        call set_rev_all2all_import_tbl(nnod_rlm, nmax_sr_rj,          &
!     &      nneib_domain_rlm, istack_sr_rlm, item_sr_rlm, irev_sr_rlm)
!        call set_rev_all2all_import_tbl(nnod_rj, nmax_sr_rj,           &
!     &      nneib_domain_rj,  istack_sr_rj,  item_sr_rj,  irev_sr_rj)
!      end if
!
      if(my_rank .ne. 0) return
      write(*,*) 'Vector length for Legendre transform:',               &
     &          nvector_legendre
      write(*,*) 'Block number for meridinal grid: ', nblock_l_rtm
      write(*,*) 'Block number for Legendre transform: ', nblock_j_rlm
!
      end subroutine set_blocks_4_leg_trans
!
! -----------------------------------------------------------------------
!
      end module init_sph_trans
