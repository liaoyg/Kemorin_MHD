!>@file   SGS_diff_field_from_trans.f90
!!@brief  module SGS_diff_field_from_trans
!!
!!@author H. Matsui
!!@date Programmed in Oct., 2009
!
!>@brief  Copy data from/to sphrical transform buffer for snapshots
!!
!!@verbatim
!!      subroutine copy_SGS_diff_field_from_trans                       &
!!     &         (sph_params, sph_rtp, trns_snap, node, iphys, nod_fld)
!!        type(sph_shell_parameters), intent(in) :: sph_params
!!        type(sph_rtp_grid), intent(in) :: sph_rtp
!!        type(address_4_sph_trans), intent(in) :: trns_snap
!!        type(node_data), intent(in) :: node
!!        type(phys_address), intent(in) :: iphys
!!        type(phys_data), intent(inout) :: nod_fld
!!@endverbatim
!
      module SGS_diff_field_from_trans
!
      use m_precision
      use m_machine_parameter
!
      use t_geometry_data
      use t_phys_address
      use t_phys_data
      use t_spheric_parameter
      use t_addresses_sph_transform
!
      implicit  none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine copy_SGS_diff_field_from_trans                         &
     &         (sph_params, sph_rtp, trns_snap, node, iphys, nod_fld)
!
      use copy_fields_from_sph_trans
!
      type(sph_shell_parameters), intent(in) :: sph_params
      type(sph_rtp_grid), intent(in) :: sph_rtp
      type(address_4_sph_trans), intent(in) :: trns_snap
      type(node_data), intent(in) :: node
      type(phys_address), intent(in) :: iphys
      type(phys_data), intent(inout) :: nod_fld
!
!
      call copy_vector_from_snap_trans                                  &
     &   (trns_snap%b_trns%i_SGS_rot_inertia, iphys%i_SGS_rot_inertia,  &
     &    sph_params%m_folding, sph_rtp, trns_snap, node, nod_fld)
!
      call copy_vector_from_snap_trans                                  &
     &   (trns_snap%b_trns%i_SGS_rot_Lorentz, iphys%i_SGS_rot_Lorentz,  &
     &    sph_params%m_folding, sph_rtp, trns_snap, node, nod_fld)
!
      call copy_vector_from_snap_trans                                  &
     &   (trns_snap%b_trns%i_SGS_induction, iphys%i_SGS_induction,      &
     &    sph_params%m_folding, sph_rtp, trns_snap, node, nod_fld)
!
!
      call copy_scalar_from_snap_trans                                  &
     &   (trns_snap%b_trns%i_SGS_div_inertia, iphys%i_SGS_div_inertia,  &
     &    sph_params%m_folding, sph_rtp, trns_snap, node, nod_fld)
!
      call copy_scalar_from_snap_trans                                  &
     &   (trns_snap%b_trns%i_SGS_div_Lorentz, iphys%i_SGS_div_Lorentz,  &
     &    sph_params%m_folding, sph_rtp, trns_snap, node, nod_fld)
!
      call copy_scalar_from_snap_trans                                  &
     &   (trns_snap%b_trns%i_SGS_div_h_flux, iphys%i_SGS_div_h_flux,    &
     &    sph_params%m_folding, sph_rtp, trns_snap, node, nod_fld)
!
      call copy_scalar_from_snap_trans                                  &
     &   (trns_snap%b_trns%i_SGS_div_c_flux, iphys%i_SGS_div_c_flux,    &
     &    sph_params%m_folding, sph_rtp, trns_snap, node, nod_fld)
!
      end  subroutine copy_SGS_diff_field_from_trans
!
!-----------------------------------------------------------------------
!
      end module SGS_diff_field_from_trans
