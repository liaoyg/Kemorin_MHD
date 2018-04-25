!>@file   dynamic_model_sph_MHD.f90
!!@brief  module dynamic_model_sph_MHD
!!
!!@author H. Matsui
!!@date Programmed in Oct., 2009
!
!>@brief  Evaluate nonlinear terms in spherical coordinate grid
!!
!!@verbatim
!!      subroutine const_model_coefs_4_sph                              &
!!     &         (SGS_param, sph_rtp, trns_SGS, dynamic_SPH)
!!        type(SGS_model_control_params), intent(in) :: SGS_param
!!        type(sph_rtp_grid), intent(in) :: sph_rtp
!!        type(address_4_sph_trans), intent(inout) :: trns_SGS
!!        type(dynamic_SGS_data_4_sph), intent(inout) :: dynamic_SPH
!!
!!      subroutine const_dynamic_SGS_4_buo_sph(stab_weight, sph_rtp,    &
!!     &          fl_prop, trns_MHD, trns_SGS, trns_DYNS, dynamic_SPH)
!!        type(SGS_model_control_params), intent(in) :: SGS_param
!!        type(sph_rtp_grid), intent(in) :: sph_rtp
!!        type(sph_dynamic_model_group), intent(in) :: sph_d_grp
!!        type(fluid_property), intent(in) :: fl_prop
!!        type(address_4_sph_trans), intent(in) :: trns_MHD
!!        type(address_4_sph_trans), intent(inout) :: trns_snap
!!        type(address_4_sph_trans), intent(inout) :: trns_SGS
!!        type(dynamic_SGS_data_4_sph), intent(inout) :: dynamic_SPH
!!@endverbatim
!
      module dynamic_model_sph_MHD
!
      use m_precision
!
      use calypso_mpi
      use m_constants
      use m_machine_parameter
      use m_phys_constants
!
      use t_SGS_control_parameter
      use t_spheric_rtp_data
      use t_groups_sph_dynamic
      use t_phys_address
      use t_addresses_sph_transform
      use t_ele_info_4_dynamic
      use t_addresses_sph_transform
      use t_SGS_model_coefs
      use t_sph_filtering
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine const_model_coefs_4_sph                                &
     &         (SGS_param, sph_rtp, trns_SGS, trns_DYNS, dynamic_SPH)
!
      use zonal_lsq_4_model_coefs
!
      type(SGS_model_control_params), intent(in) :: SGS_param
      type(sph_rtp_grid), intent(in) :: sph_rtp
!
      type(address_4_sph_trans), intent(in) :: trns_SGS, trns_DYNS
      type(dynamic_SGS_data_4_sph), intent(inout) :: dynamic_SPH
!
!
      if(dynamic_SPH%ifld_sgs%i_mom_flux .gt. 0) then
        if (iflag_debug.eq.1) write(*,*) 'cal_dynamic_SGS_4_sph_MHD MF'
        call cal_dynamic_SGS_4_sph_MHD                                  &
     &     (sph_rtp, dynamic_SPH%sph_d_grp, SGS_param%stab_weight,      &
     &      n_vector, dynamic_SPH%ifld_sgs%i_mom_flux,                  &
     &      trns_SGS%frc_rtp(1,trns_SGS%f_trns%i_SGS_Lorentz),          &
     &      trns_DYNS%fld_rtp(1,trns_DYNS%b_trns%i_wide_SGS_inertia),   &
     &      trns_DYNS%fld_rtp(1,trns_DYNS%b_trns%i_dbl_SGS_inertia),    &
     &      dynamic_SPH%wk_sgs)
      end if
!
      if(dynamic_SPH%ifld_sgs%i_lorentz .gt. 0) then
        if (iflag_debug.eq.1) write(*,*) 'cal_dynamic_SGS_4_sph_MHD LZ'
        call cal_dynamic_SGS_4_sph_MHD                                  &
     &     (sph_rtp, dynamic_SPH%sph_d_grp, SGS_param%stab_weight,      &
     &      n_vector, dynamic_SPH%ifld_sgs%i_lorentz,                   &
     &      trns_SGS%frc_rtp(1,trns_SGS%f_trns%i_SGS_Lorentz),          &
     &      trns_DYNS%fld_rtp(1,trns_DYNS%b_trns%i_wide_SGS_Lorentz),   &
     &      trns_DYNS%fld_rtp(1,trns_DYNS%b_trns%i_dbl_SGS_Lorentz),    &
     &      dynamic_SPH%wk_sgs)
      end if
!
      if(dynamic_SPH%ifld_sgs%i_induction .gt. 0) then
        if (iflag_debug.eq.1) write(*,*) 'cal_dynamic_SGS_4_sph_MHD ID'
        call cal_dynamic_SGS_4_sph_MHD                                  &
     &     (sph_rtp, dynamic_SPH%sph_d_grp, SGS_param%stab_weight,      &
     &      n_vector, dynamic_SPH%ifld_sgs%i_induction,                 &
     &      trns_SGS%frc_rtp(1,trns_SGS%f_trns%i_SGS_vp_induct),        &
     &      trns_DYNS%fld_rtp(1,trns_DYNS%b_trns%i_wide_SGS_vp_induct), &
     &      trns_DYNS%fld_rtp(1,trns_DYNS%b_trns%i_dbl_SGS_vp_induct),  &
     &      dynamic_SPH%wk_sgs)
      end if
!
      if(dynamic_SPH%ifld_sgs%i_heat_flux .gt. 0) then
        if (iflag_debug.eq.1) write(*,*) 'cal_dynamic_SGS_4_sph_MHD HF'
        call cal_dynamic_SGS_4_sph_MHD                                  &
     &     (sph_rtp, dynamic_SPH%sph_d_grp, SGS_param%stab_weight,      &
     &      n_vector, dynamic_SPH%ifld_sgs%i_heat_flux,                 &
     &      trns_SGS%frc_rtp(1,trns_SGS%f_trns%i_SGS_h_flux),           &
     &      trns_DYNS%fld_rtp(1,trns_DYNS%b_trns%i_wide_SGS_h_flux),    &
     &      trns_DYNS%fld_rtp(1,trns_DYNS%b_trns%i_dbl_SGS_h_flux),     &
     &      dynamic_SPH%wk_sgs)
      end if
!
      if(dynamic_SPH%ifld_sgs%i_comp_flux .gt. 0) then
        if (iflag_debug.eq.1) write(*,*) 'cal_dynamic_SGS_4_sph_MHD CF'
        call cal_dynamic_SGS_4_sph_MHD                                  &
     &     (sph_rtp, dynamic_SPH%sph_d_grp, SGS_param%stab_weight,      &
     &      n_vector, dynamic_SPH%ifld_sgs%i_comp_flux,                 &
     &      trns_SGS%frc_rtp(1,trns_SGS%f_trns%i_SGS_c_flux),           &
     &      trns_DYNS%fld_rtp(1,trns_DYNS%b_trns%i_wide_SGS_c_flux),    &
     &      trns_DYNS%fld_rtp(1,trns_DYNS%b_trns%i_dbl_SGS_c_flux),     &
     &      dynamic_SPH%wk_sgs)
      end if
!
      end subroutine const_model_coefs_4_sph
!
! ----------------------------------------------------------------------
!
      subroutine const_dynamic_SGS_4_buo_sph(stab_weight, sph_rtp,      &
     &          fl_prop, trns_MHD, trns_SGS, trns_DYNS, dynamic_SPH)
!
      use SGS_buo_coefs_sph_MHD
      use cal_SGS_buo_flux_sph_MHD
!
      real(kind = kreal), intent(in) :: stab_weight
      type(sph_rtp_grid), intent(in) :: sph_rtp
      type(fluid_property), intent(in) :: fl_prop
      type(address_4_sph_trans), intent(in) :: trns_MHD
      type(address_4_sph_trans), intent(in) :: trns_SGS
!
      type(address_4_sph_trans), intent(inout) :: trns_DYNS
      type(dynamic_SGS_data_4_sph), intent(inout) :: dynamic_SPH
!
!
      call SGS_fluxes_for_buo_coefs(sph_rtp, fl_prop,                   &
     &    trns_MHD%b_trns, trns_SGS%f_trns, trns_DYNS%f_trns,           &
     &    trns_MHD%ncomp_rj_2_rtp, trns_SGS%ncomp_rtp_2_rj,             &
     &    trns_DYNS%ncomp_rtp_2_rj, trns_MHD%fld_rtp, trns_SGS%frc_rtp, &
     &    trns_DYNS%frc_rtp)
!
      if(dynamic_SPH%ifld_sgs%i_buoyancy .gt. 0) then
        call cal_SGS_buo_coefs_sph_MHD                                  &
     &     (sph_rtp, dynamic_SPH%sph_d_grp, stab_weight,                &
     &      trns_DYNS%frc_rtp, trns_DYNS%ncomp_rtp_2_rj,                &
     &      trns_DYNS%f_trns%i_reynolds_wk,                             &
     &      trns_DYNS%f_trns%i_SGS_buo_wk,                              &
     &      dynamic_SPH%ifld_sgs%i_buoyancy,                            &
     &      dynamic_SPH%icomp_sgs%i_buoyancy, dynamic_SPH%wk_sgs)
      end if
!
      if(dynamic_SPH%ifld_sgs%i_comp_buoyancy .gt. 0) then
        call cal_SGS_buo_coefs_sph_MHD                                  &
     &     (sph_rtp, dynamic_SPH%sph_d_grp, stab_weight,                &
     &      trns_DYNS%frc_rtp, trns_DYNS%ncomp_rtp_2_rj,                &
     &      trns_DYNS%f_trns%i_reynolds_wk,                             &
     &      trns_DYNS%f_trns%i_SGS_comp_buo_wk,                         &
     &      dynamic_SPH%ifld_sgs%i_comp_buoyancy,                       &
     &      dynamic_SPH%icomp_sgs%i_comp_buoyancy, dynamic_SPH%wk_sgs)
      end if
!
      end subroutine const_dynamic_SGS_4_buo_sph 
!
! ----------------------------------------------------------------------
!
      end module dynamic_model_sph_MHD
