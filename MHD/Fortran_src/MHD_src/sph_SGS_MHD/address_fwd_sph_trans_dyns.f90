!>@file   address_fwd_sph_trans_dyns.f90
!!@brief  module address_fwd_sph_trans_dyns
!!
!!@author H. Matsui
!!@date Programmed in Jan., 2010
!
!>@brief Field addresses for spherical harmonics transform
!!       in MHD dynamo simulation
!!
!!@verbatim
!!      subroutine f_trans_address_vector_DYNS(trns_DYNS)
!!      subroutine f_trans_address_scalar_DYNS(ipol, trns_DYNS)
!!        type(phys_address), intent(in) :: ipol
!!        type(address_4_sph_trans), intent(inout) :: trns_DYNS
!!
!!      subroutine set_f_trans_vector_field_DYNS                        &
!!     &         (icou, ipol, itor, iphys, trns_DYNS)
!!      subroutine set_f_trans_scalar_field_DYNS                        &
!!     &         (icou, ipol, itor, iphys, trns_DYNS)
!!        type(phys_address), intent(in) :: ipol, itor, iphys
!!        type(address_4_sph_trans), intent(inout) :: trns_DYNS
!!@endverbatim
!
      module address_fwd_sph_trans_dyns
!
      use m_precision
!
      use m_phys_labels
      use t_phys_address
      use t_addresses_sph_transform
      use t_control_parameter
      use t_physical_property
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine f_trans_address_vector_DYNS(trns_DYNS)
!
!      type(phys_address), intent(in) :: ipol
      type(address_4_sph_trans), intent(inout) :: trns_DYNS
!
!
      trns_DYNS%forward%num_vector = 0
!   SGS advection flag
!      call add_vector_trans_flag(ipol%i_SGS_inertia,                   &
!     &    trns_DYNS%forward%num_vector, trns_DYNS%f_trns%i_SGS_inertia)
!   SGS Lorentz force flag
!      call add_vector_trans_flag(ipol%i_SGS_Lorentz,                   &
!     &    trns_DYNS%forward%num_vector, trns_DYNS%f_trns%i_SGS_Lorentz)
!   SGS induction flag
!      call add_vector_trans_flag(ipol%i_SGS_vp_induct,                 &
!     &    trns_DYNS%forward%num_vector, trns_DYNS%f_trns%i_SGS_vp_induct)
!   SGS heat flux flag
!      call add_vector_trans_flag(ipol%i_SGS_h_flux,                    &
!     &    trns_DYNS%forward%num_vector, trns_DYNS%f_trns%i_SGS_h_flux)
!   SGS composition flux flag
!      call add_vector_trans_flag(ipol%i_SGS_c_flux,                    &
!     &    trns_DYNS%forward%num_vector, trns_DYNS%f_trns%i_SGS_c_flux)
!
      end subroutine f_trans_address_vector_DYNS
!
!-----------------------------------------------------------------------
!
      subroutine f_trans_address_scalar_DYNS(ipol, trns_DYNS)
!
      type(phys_address), intent(in) :: ipol
      type(address_4_sph_trans), intent(inout) :: trns_DYNS
!
!
      trns_DYNS%forward%num_scalar = 0
!   work of Reynolds stress
      call add_scalar_trans_flag(ipol%i_reynolds_wk,                    &
     &    trns_DYNS%forward%num_vector, trns_DYNS%forward%num_scalar,       &
     &    trns_DYNS%f_trns%i_reynolds_wk)
!   work of SGS buoyancy
      call add_scalar_trans_flag(ipol%i_SGS_buo_wk,                     &
     &    trns_DYNS%forward%num_vector, trns_DYNS%forward%num_scalar,       &
     &    trns_DYNS%f_trns%i_SGS_buo_wk)
!   work of SGS compositional buoyancy
      call add_scalar_trans_flag(ipol%i_SGS_comp_buo_wk,                &
     &    trns_DYNS%forward%num_vector, trns_DYNS%forward%num_scalar,       &
     &    trns_DYNS%f_trns%i_SGS_comp_buo_wk)
!
      end subroutine f_trans_address_scalar_DYNS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine set_f_trans_vector_field_DYNS                          &
     &         (icou, ipol, itor, iphys, trns_DYNS)
!
      type(phys_address), intent(in) :: ipol, itor, iphys
      type(address_4_sph_trans), intent(inout) :: trns_DYNS
      integer(kind = kint), intent(inout) :: icou
!
!
!   SGS advection flag
      call set_field_name_4_fwd_trns                                    &
     &   (fhd_SGS_inertia, trns_DYNS%f_trns%i_SGS_inertia,              &
     &    ipol%i_SGS_inertia, itor%i_SGS_inertia, iphys%i_SGS_inertia,  &
     &    icou, trns_DYNS)
      call set_field_name_4_fwd_trns                                    &
     &   (fhd_SGS_Lorentz, trns_DYNS%f_trns%i_SGS_Lorentz,              &
     &    ipol%i_SGS_Lorentz, itor%i_SGS_Lorentz, iphys%i_SGS_Lorentz,  &
     &    icou, trns_DYNS)
      call set_field_name_4_fwd_trns                                    &
     &   (fhd_SGS_vp_induct, trns_DYNS%f_trns%i_SGS_vp_induct,          &
     &    ipol%i_SGS_vp_induct, itor%i_SGS_vp_induct,                   &
     &    iphys%i_SGS_vp_induct, icou, trns_DYNS)
      call set_field_name_4_fwd_trns                                    &
     &   (fhd_SGS_h_flux, trns_DYNS%f_trns%i_SGS_h_flux,                &
     &    ipol%i_SGS_h_flux, itor%i_SGS_h_flux, iphys%i_SGS_h_flux,     &
     &    icou, trns_DYNS)
      call set_field_name_4_fwd_trns                                    &
     &   (fhd_SGS_c_flux, trns_DYNS%f_trns%i_SGS_c_flux,                &
     &    ipol%i_SGS_c_flux, itor%i_SGS_c_flux, iphys%i_SGS_c_flux,     &
     &    icou, trns_DYNS)
!
      end subroutine set_f_trans_vector_field_DYNS
!
!-----------------------------------------------------------------------
!
      subroutine set_f_trans_scalar_field_DYNS                          &
     &         (icou, ipol, itor, iphys, trns_DYNS)
!
      type(phys_address), intent(in) :: ipol, itor, iphys
      type(address_4_sph_trans), intent(inout) :: trns_DYNS
      integer(kind = kint), intent(inout) :: icou
!
!
      call set_field_name_4_fwd_trns                                    &
     &   (fhd_Reynolds_work, trns_DYNS%f_trns%i_reynolds_wk,            &
     &    ipol%i_reynolds_wk, itor%i_reynolds_wk, iphys%i_reynolds_wk,  &
     &    icou, trns_DYNS)
      call set_field_name_4_fwd_trns                                    &
     &   (fhd_SGS_buo_flux, trns_DYNS%f_trns%i_SGS_buo_wk,              &
     &    ipol%i_SGS_buo_wk, itor%i_SGS_buo_wk, iphys%i_SGS_buo_wk,     &
     &    icou, trns_DYNS)
      call set_field_name_4_fwd_trns                                    &
     &   (fhd_SGS_comp_buo_flux, trns_DYNS%f_trns%i_SGS_comp_buo_wk,    &
     &    ipol%i_SGS_comp_buo_wk, itor%i_SGS_comp_buo_wk,               &
     &    iphys%i_SGS_comp_buo_wk, icou, trns_DYNS)
!
!   SGS buoyancy
      call set_field_name_4_fwd_trns                                    &
     &   (fhd_Csim_SGS_buoyancy, trns_DYNS%f_trns%i_Csim_SGS_buoyancy,  &
     &    ipol%i_Csim_SGS_buoyancy, itor%i_Csim_SGS_buoyancy,           &
     &    iphys%i_Csim_SGS_buoyancy, icou, trns_DYNS)
      call set_field_name_4_fwd_trns                                    &
     &   (fhd_Csim_SGS_comp_buo, trns_DYNS%f_trns%i_Csim_SGS_comp_buo,  &
     &    ipol%i_Csim_SGS_comp_buo, itor%i_Csim_SGS_comp_buo,           &
     &    iphys%i_Csim_SGS_comp_buo, icou, trns_DYNS)
!
      end subroutine set_f_trans_scalar_field_DYNS
!
!-----------------------------------------------------------------------
!
      end module address_fwd_sph_trans_dyns
