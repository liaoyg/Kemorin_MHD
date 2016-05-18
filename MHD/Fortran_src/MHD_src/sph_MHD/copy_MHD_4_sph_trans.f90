!>@file   copy_MHD_4_sph_trans.f90
!!@brief  module copy_MHD_4_sph_trans
!!
!!@author H. Matsui
!!@date    programmed by H.Matsui in Oct., 2009
!
!>@brief Copy spectrum data and field data to spherical transform buffer
!!       for dynamo simulation
!!
!!@verbatim
!!  routines for backward transform
!!      subroutine select_mhd_field_from_trans                          &
!!     &         (sph_rtp, f_trns, ncomp_rtp_2_rj, frc_rtp, frm_rtp)
!!        type(sph_rtp_grid), intent(in) :: sph_rtp
!!      subroutine copy_forces_to_snapshot_rtp                          &
!!     &          (m_folding, sph_rtp, f_trns, ncomp_rtp_2_rj, node,    &
!!     &           iphys, nod_fld)
!!        type(node_data), intent(in) :: node
!!        type(phys_address), intent(in) :: iphys
!!        type(phys_data), intent(inout) :: nod_fld
!!@endverbatim
!
      module copy_MHD_4_sph_trans
!
      use m_precision
      use m_machine_parameter
!
      use t_geometry_data
      use t_phys_address
      use t_phys_data
      use t_spheric_rtp_data
!
      implicit  none
!
      private :: sel_force_from_MHD_trans, copy_force_from_MHD_trans
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine select_mhd_field_from_trans                            &
     &         (sph_rtp, f_trns, ncomp_rtp_2_rj, frc_rtp, frm_rtp)
!
      use m_node_phys_data
!
      type(sph_rtp_grid), intent(in) :: sph_rtp
      type(phys_address), intent(in) :: f_trns
!
      integer(kind = kint), intent(in) :: ncomp_rtp_2_rj
      real(kind = kreal), intent(in)                                 &
     &                   :: frc_rtp(sph_rtp%nnod_rtp,ncomp_rtp_2_rj)
      real(kind = kreal), intent(inout)                                 &
     &                   :: frm_rtp(sph_rtp%nnod_rtp,ncomp_rtp_2_rj)
!
!
!   advection flag
      call sel_force_from_MHD_trans                                     &
     &   (sph_rtp, f_trns%i_m_advect, ncomp_rtp_2_rj, frc_rtp, frm_rtp)
!   Coriolis flag
      call sel_force_from_MHD_trans                                     &
     &   (sph_rtp, f_trns%i_coriolis, ncomp_rtp_2_rj, frc_rtp, frm_rtp)
!   Lorentz flag
      call sel_force_from_MHD_trans                                     &
     &   (sph_rtp, f_trns%i_lorentz, ncomp_rtp_2_rj, frc_rtp, frm_rtp)
!
!   induction flag
      call sel_force_from_MHD_trans                                     &
     &   (sph_rtp, f_trns%i_vp_induct, ncomp_rtp_2_rj,                  &
     &    frc_rtp, frm_rtp)
!   divergence of heat flux flag
      call sel_force_from_MHD_trans                                     &
     &   (sph_rtp, f_trns%i_h_flux, ncomp_rtp_2_rj, frc_rtp, frm_rtp)
!
!   divergence of composition flux flag
      call sel_force_from_MHD_trans                                     &
     &   (sph_rtp, f_trns%i_c_flux, ncomp_rtp_2_rj, frc_rtp, frm_rtp)
!
      end subroutine select_mhd_field_from_trans
!
!-----------------------------------------------------------------------
!
      subroutine copy_forces_to_snapshot_rtp                            &
     &          (m_folding, sph_rtp, f_trns, ncomp_rtp_2_rj, node,      &
     &           iphys, nod_fld)
!
      use m_addresses_trans_sph_snap
!
      type(sph_rtp_grid), intent(in) :: sph_rtp
      type(phys_address), intent(in) :: f_trns
      type(node_data), intent(in) :: node
      type(phys_address), intent(in) :: iphys
!
      integer(kind = kint), intent(in) :: ncomp_rtp_2_rj
!
      integer(kind = kint), intent(in) :: m_folding
!
      type(phys_data), intent(inout) :: nod_fld
!
!
!   advection flag
      call copy_force_from_MHD_trans                                    &
     &   (f_trns%i_m_advect, iphys%i_m_advect,                          &
     &    m_folding, sph_rtp, node, ncomp_rtp_2_rj, frm_rtp, nod_fld)
!   Coriolis flag
      call copy_force_from_MHD_trans                                    &
     &   (f_trns%i_coriolis, iphys%i_coriolis,                          &
     &    m_folding, sph_rtp, node, ncomp_rtp_2_rj, frm_rtp, nod_fld)
!   Lorentz flag
      call copy_force_from_MHD_trans                                    &
     &   (f_trns%i_lorentz, iphys%i_lorentz,                            &
     &    m_folding, sph_rtp, node, ncomp_rtp_2_rj, frm_rtp, nod_fld)
!
!   induction flag
      call copy_force_from_MHD_trans                                    &
     &   (f_trns%i_vp_induct, iphys%i_vp_induct,                        &
     &    m_folding, sph_rtp, node, ncomp_rtp_2_rj, frm_rtp, nod_fld)
!   divergence of heat flux flag
      call copy_force_from_MHD_trans                                    &
     &   (f_trns%i_h_flux, iphys%i_h_flux,                              &
     &    m_folding, sph_rtp, node, ncomp_rtp_2_rj, frm_rtp, nod_fld)
!
!   divergence of composition flux flag
      call copy_force_from_MHD_trans                                    &
     &   (f_trns%i_c_flux, iphys%i_c_flux,                              &
     &    m_folding, sph_rtp, node, ncomp_rtp_2_rj, frm_rtp, nod_fld)
!
      end subroutine copy_forces_to_snapshot_rtp
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine sel_force_from_MHD_trans                               &
     &         (sph_rtp, i_trns, ncomp_rtp_2_rj, frc_rtp, frm_rtp)
!
      use sel_fld_copy_4_sph_trans
!
      type(sph_rtp_grid), intent(in) :: sph_rtp
      integer(kind = kint), intent(in) :: i_trns
      integer(kind = kint), intent(in) :: ncomp_rtp_2_rj
!
      real(kind = kreal), intent(in)                                    &
     &                   :: frc_rtp(sph_rtp%nnod_rtp,ncomp_rtp_2_rj)
      real(kind = kreal), intent(inout)                                 &
     &                   :: frm_rtp(sph_rtp%nnod_rtp,ncomp_rtp_2_rj)
!
!
      if(i_trns .le. 0) return
      call sel_vector_from_trans(sph_rtp%nnod_rtp, sph_rtp%nidx_rtp,    &
     &    ione, sph_rtp%istack_inod_rtp_smp, sph_rtp%nnod_rtp,          &
     &    frc_rtp(1,i_trns), frm_rtp(1,i_trns) )
!
      end subroutine sel_force_from_MHD_trans
!
!-----------------------------------------------------------------------
!
      subroutine copy_force_from_MHD_trans                              &
     &         (i_trns, i_field, m_folding, sph_rtp, node,              &
     &          ncomp_rtp_2_rj, frm_rtp, nod_fld)
!
      use copy_nodal_fld_4_sph_trans
!
      type(sph_rtp_grid), intent(in) :: sph_rtp
      type(node_data), intent(in) :: node
      integer(kind = kint), intent(in) :: i_field, i_trns
      integer(kind = kint), intent(in) :: m_folding
      integer(kind = kint), intent(in) :: ncomp_rtp_2_rj
      real(kind = kreal), intent(in)                                    &
     &                   :: frm_rtp(sph_rtp%nnod_rtp,ncomp_rtp_2_rj)
!
      type(phys_data), intent(inout) :: nod_fld
!
!
      if( (i_field*i_trns) .le. 0) return
      call copy_nodal_vector_from_trans                                 &
     &   (sph_rtp, m_folding, ncomp_rtp_2_rj, i_trns, frm_rtp,          &
     &    node%numnod, nod_fld%ntot_phys, i_field, nod_fld%d_fld)
!
      end subroutine copy_force_from_MHD_trans
!
!-----------------------------------------------------------------------
!
      end module copy_MHD_4_sph_trans
