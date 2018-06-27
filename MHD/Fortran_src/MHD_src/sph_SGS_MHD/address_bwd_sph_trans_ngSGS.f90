!>@file   address_bwd_sph_trans_ngSGS.f90
!!@brief  module address_bwd_sph_trans_ngSGS
!!
!!@author H. Matsui
!!@date Programmed in Jan., 2010
!
!>@brief Field addresses for spherical harmonics transform
!!       in MHD dynamo simulation
!!
!!@verbatim
!!      subroutine b_trans_address_vector_ngSGS                         &
!!     &         (ipol, itor, iphys, b_trns, trns_back)
!!      subroutine b_trans_filter_vector_grads                          &
!!     &         (ipol, itor, iphys, b_trns, trns_back)
!!        type(phys_address), intent(in) :: ipol, itor, iphys
!!        type(address_each_sph_trans), intent(inout) :: trns_back
!!        type(phys_address), intent(inout) :: b_trns
!!@endverbatim
!
      module address_bwd_sph_trans_ngSGS
!
      use m_precision
!
      use m_phys_labels
      use m_phys_constants
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
      subroutine b_trans_vector_gradients                            &
     &         (ipol, itor, iphys, b_trns, trns_back)
!
      type(phys_address), intent(in) :: ipol, itor, iphys
      type(address_each_sph_trans), intent(inout) :: trns_back
      type(phys_address), intent(inout) :: b_trns
!
!
!   Gradient of Radial velocity
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_vx, fhd_grad_v_1, n_vector, ipol%i_grad_vx,       &
     &    itor%i_grad_vx, iphys%i_grad_vx, b_trns%i_grad_vx, trns_back)
!   Gradient of meridional velocity
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_vy, fhd_grad_v_2, n_vector, ipol%i_grad_vy,       &
     &    itor%i_grad_vy, iphys%i_grad_vy, b_trns%i_grad_vy, trns_back)
!   Gradient of zonal velocity
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_vz, fhd_grad_v_3, n_vector, ipol%i_grad_vz,       &
     &    itor%i_grad_vz, iphys%i_grad_vz, b_trns%i_grad_vz, trns_back)
!
!   Gradient of Radial vorticity
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_wx, fhd_grad_w_1, n_vector, ipol%i_grad_wx,       &
     &    itor%i_grad_wx, iphys%i_grad_wx, b_trns%i_grad_wx, trns_back)
!   Gradient of meridional vorticity
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_wy, fhd_grad_w_2, n_vector, ipol%i_grad_wy,       &
     &    itor%i_grad_wy, iphys%i_grad_wy, b_trns%i_grad_wy, trns_back)
!   Gradient of zonal vorticity
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_wz, fhd_grad_w_3, n_vector, ipol%i_grad_wz,       &
     &    itor%i_grad_wz, iphys%i_grad_wz, b_trns%i_grad_wz, trns_back)
!
!   Gradient of Radial magnetic field
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_bx, fhd_grad_b_1, n_vector, ipol%i_grad_bx,       &
     &    itor%i_grad_bx, iphys%i_grad_bx, b_trns%i_grad_bx, trns_back)
!   Gradient of meridional magnetic field
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_by, fhd_grad_b_2, n_vector, ipol%i_grad_by,       &
     &    itor%i_grad_by, iphys%i_grad_by, b_trns%i_grad_by, trns_back)
!   Gradient of zonal magnetic field
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_bz, fhd_grad_b_3, n_vector, ipol%i_grad_bz,       &
     &    itor%i_grad_bz, iphys%i_grad_bz, b_trns%i_grad_bz, trns_back)
!
!   Gradient of Radial current density
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_jx, fhd_grad_j_1, n_vector, ipol%i_grad_jx,       &
     &    itor%i_grad_jx, iphys%i_grad_jx, b_trns%i_grad_jx, trns_back)
!   Gradient of meridional current density
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_jy, fhd_grad_j_2, n_vector, ipol%i_grad_jy,       &
     &    itor%i_grad_jy, iphys%i_grad_jy, b_trns%i_grad_jy, trns_back)
!   Gradient of zonal current density
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_jz, fhd_grad_j_3, n_vector, ipol%i_grad_jz,       &
     &    itor%i_grad_jz, iphys%i_grad_jz, b_trns%i_grad_jz, trns_back)
!
!   Gradient of temperature
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_t, fhd_grad_temp, n_vector, ipol%i_grad_t,        &
     &    itor%i_grad_t, iphys%i_grad_t, b_trns%i_grad_t, trns_back)
!
!   Gradient of composition
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_composit, fhd_grad_composit, n_vector,            &
     &    ipol%i_grad_composit, itor%i_grad_composit,                   &
     &    iphys%i_grad_composit, b_trns%i_grad_composit, trns_back)
!
      end subroutine b_trans_vector_gradients
!
!-----------------------------------------------------------------------
!
      subroutine b_trans_filter_vector_grads                            &
     &         (ipol, itor, iphys, b_trns, trns_back)
!
      type(phys_address), intent(in) :: ipol, itor, iphys
      type(address_each_sph_trans), intent(inout) :: trns_back
      type(phys_address), intent(inout) :: b_trns
!
!
!   Gradient of Radial velocity
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_filter_vx, fhd_grad_filter_v_1, n_vector,         &
     &    ipol%i_grad_filter_vx, itor%i_grad_filter_vx,                 &
     &    iphys%i_grad_filter_vx, b_trns%i_grad_filter_vx, trns_back)
!   Gradient of meridional velocity
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_filter_vy, fhd_grad_filter_v_2, n_vector,         &
     &    ipol%i_grad_filter_vy, itor%i_grad_filter_vy,                 &
     &    iphys%i_grad_filter_vy, b_trns%i_grad_filter_vy, trns_back)
!   Gradient of zonal velocity
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_filter_vz, fhd_grad_filter_v_3, n_vector,         &
     &    ipol%i_grad_filter_vz, itor%i_grad_filter_vz,                 &
     &    iphys%i_grad_filter_vz, b_trns%i_grad_filter_vz, trns_back)
!
!   Gradient of Radial vorticity
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_filter_wx, fhd_grad_filter_w_1, n_vector,         &
     &    ipol%i_grad_filter_wx, itor%i_grad_filter_wx,                 &
     &    iphys%i_grad_filter_wx, b_trns%i_grad_filter_wx, trns_back)
!   Gradient of meridional vorticity
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_filter_wy, fhd_grad_filter_w_2, n_vector,         &
     &    ipol%i_grad_filter_wy, itor%i_grad_filter_wy,                 &
     &    iphys%i_grad_filter_wy, b_trns%i_grad_filter_wy, trns_back)
!   Gradient of zonal vorticity
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_filter_wz, fhd_grad_filter_w_3, n_vector,         &
     &    ipol%i_grad_filter_wz, itor%i_grad_filter_wz,                 &
     &    iphys%i_grad_filter_wz, b_trns%i_grad_filter_wz, trns_back)
!
!   Gradient of Radial magnetic field
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_filter_bx, fhd_grad_filter_b_1, n_vector,         &
     &    ipol%i_grad_filter_bx, itor%i_grad_filter_bx,                 &
     &    iphys%i_grad_filter_bx, b_trns%i_grad_filter_bx, trns_back)
!   Gradient of meridional magnetic field
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_filter_by, fhd_grad_filter_b_2, n_vector,         &
     &    ipol%i_grad_filter_by, itor%i_grad_filter_by,                 &
     &    iphys%i_grad_filter_by, b_trns%i_grad_filter_by, trns_back)
!   Gradient of zonal magnetic field
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_filter_bz, fhd_grad_filter_b_3, n_vector,         &
     &    ipol%i_grad_filter_bz, itor%i_grad_filter_bz,                 &
     &    iphys%i_grad_filter_bz, b_trns%i_grad_filter_bz, trns_back)
!
!   Gradient of Radial current density
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_filter_jx, fhd_grad_filter_j_1, n_vector,         &
     &    ipol%i_grad_filter_jx, itor%i_grad_filter_jx,                 &
     &    iphys%i_grad_filter_jx, b_trns%i_grad_filter_jx, trns_back)
!   Gradient of meridional current density
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_filter_jy, fhd_grad_filter_j_2, n_vector,         &
     &    ipol%i_grad_filter_jy, itor%i_grad_filter_jy,                 &
     &    iphys%i_grad_filter_jy, b_trns%i_grad_filter_jy, trns_back)
!   Gradient of zonal current density
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_filter_jz, fhd_grad_filter_j_3, n_vector,         &
     &    ipol%i_grad_filter_jz, itor%i_grad_filter_jz,                 &
     &    iphys%i_grad_filter_jz, b_trns%i_grad_filter_jz, trns_back)
!
!   Gradient of temperature
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_filter_temp, fhd_grad_filter_temp, n_vector,      &
     &    ipol%i_grad_filter_temp, itor%i_grad_filter_temp,             &
     &    iphys%i_grad_filter_temp, b_trns%i_grad_filter_temp,          &
     &    trns_back)
!
!   Gradient of composition
      call add_field_name_4_sph_trns                                    &
     &   (ipol%i_grad_filter_comp, fhd_grad_filter_comp, n_vector,      &
     &    ipol%i_grad_filter_comp, itor%i_grad_filter_comp,             &
     &    iphys%i_grad_filter_comp, b_trns%i_grad_filter_comp,          &
     &    trns_back)
!
      end subroutine b_trans_filter_vector_grads
!
!-----------------------------------------------------------------------
!
      end module address_bwd_sph_trans_ngSGS