!>@file   add_sph_SGS_MHD_fld_2_ctl.f90
!!@brief  module add_sph_SGS_MHD_fld_2_ctl
!!
!!@author H. Matsui
!!@date    programmed by H.Matsui in Sep., 2009
!
!>@brief Add fields in control list for MHD dynamo simulation
!!
!!@verbatim
!!      subroutine add_field_name_4_SGS(SGS_param, field_ctl)
!!      subroutine add_field_name_dynamic_SGS                           &
!!     &         (SGS_param, fl_prop, field_ctl)
!!        type(SGS_model_control_params), intent(in) :: SGS_param
!!        type(fluid_property), intent(in) :: fl_prop
!!        type(conductive_property), intent(in) :: cd_prop
!!        type(scalar_property), intent(in) :: ht_prop, cp_prop
!!        type(ctl_array_c3), intent(inout) :: field_ctl
!!@endverbatim
!
      module add_sph_SGS_MHD_fld_2_ctl
!
      use m_precision
!
      use m_phys_labels
      use t_SGS_control_parameter
      use t_read_control_arrays
      use t_physical_property
!
      implicit  none
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine add_field_name_4_SGS(SGS_param, field_ctl)
!
      use add_nodal_fields_ctl
!
      type(SGS_model_control_params), intent(in) :: SGS_param
      type(ctl_array_c3), intent(inout) :: field_ctl
!
!
      if(SGS_param%iflag_SGS_h_flux .gt. id_SGS_none) then
        call add_phys_name_ctl(fhd_filter_velo, field_ctl)
        call add_phys_name_ctl(fhd_filter_temp, field_ctl)
!
        call add_phys_name_ctl(fhd_SGS_h_flux, field_ctl)
        call add_phys_name_ctl(fhd_div_SGS_h_flux, field_ctl)
      end if
!
      if(SGS_param%iflag_SGS_c_flux .gt. id_SGS_none) then
        call add_phys_name_ctl(fhd_filter_velo, field_ctl)
        call add_phys_name_ctl(fhd_filter_comp, field_ctl)
!
        call add_phys_name_ctl(fhd_SGS_c_flux, field_ctl)
        call add_phys_name_ctl(fhd_div_SGS_c_flux, field_ctl)
      end if
!
      if(SGS_param%iflag_SGS_m_flux .gt. id_SGS_none) then
        call add_phys_name_ctl(fhd_filter_velo, field_ctl)
        call add_phys_name_ctl(fhd_filter_vort, field_ctl)
!
        call add_phys_name_ctl(fhd_SGS_inertia, field_ctl)
        call add_phys_name_ctl(fhd_SGS_rot_inertia, field_ctl)
        call add_phys_name_ctl(fhd_SGS_div_inertia, field_ctl)
      end if
!
      if(SGS_param%iflag_SGS_lorentz .gt. id_SGS_none) then
        call add_phys_name_ctl(fhd_filter_magne, field_ctl)
        call add_phys_name_ctl(fhd_filter_current, field_ctl)
!
        call add_phys_name_ctl(fhd_SGS_Lorentz, field_ctl)
        call add_phys_name_ctl(fhd_SGS_rot_Lorentz, field_ctl)
        call add_phys_name_ctl(fhd_SGS_div_Lorentz, field_ctl)
      end if
!
      if(SGS_param%iflag_SGS_uxb .gt. id_SGS_none) then
        call add_phys_name_ctl(fhd_filter_velo, field_ctl)
        call add_phys_name_ctl(fhd_filter_magne, field_ctl)
!
        call add_phys_name_ctl(fhd_SGS_vp_induct, field_ctl)
        call add_phys_name_ctl(fhd_SGS_induction, field_ctl)
      end if
!
      end subroutine add_field_name_4_SGS
!
! -----------------------------------------------------------------------
!
      subroutine add_field_name_dynamic_SGS                             &
     &         (SGS_param, fl_prop, field_ctl)
!
      use add_nodal_fields_ctl
!
      type(SGS_model_control_params), intent(in) :: SGS_param
      type(fluid_property), intent(in) :: fl_prop
      type(ctl_array_c3), intent(inout) :: field_ctl
!
!
      if(SGS_param%iflag_dynamic .eq. id_SGS_DYNAMIC_OFF) return
!
      if(SGS_param%iflag_SGS_h_flux .gt. id_SGS_none) then
        call add_phys_name_ctl(fhd_w_filter_velo, field_ctl)
        call add_phys_name_ctl(fhd_w_filter_temp, field_ctl)
!
        call add_phys_name_ctl(fhd_wide_SGS_h_flux, field_ctl)
        call add_phys_name_ctl(fhd_Csim_SGS_h_flux, field_ctl)
      end if
!
      if(SGS_param%iflag_SGS_c_flux .gt. id_SGS_none) then
        call add_phys_name_ctl(fhd_w_filter_velo, field_ctl)
        call add_phys_name_ctl(fhd_w_filter_comp, field_ctl)
!
        call add_phys_name_ctl(fhd_wide_SGS_c_flux, field_ctl)
        call add_phys_name_ctl(fhd_Csim_SGS_c_flux, field_ctl)
      end if
!
      if(SGS_param%iflag_SGS_m_flux .gt. id_SGS_none) then
        call add_phys_name_ctl(fhd_w_filter_velo, field_ctl)
        call add_phys_name_ctl(fhd_w_filter_vort, field_ctl)
!
        call add_phys_name_ctl(fhd_wide_SGS_inertia, field_ctl)
        call add_phys_name_ctl(fhd_Csim_SGS_m_flux, field_ctl)
      end if
!
      if(SGS_param%iflag_SGS_lorentz .gt. id_SGS_none) then
        call add_phys_name_ctl(fhd_w_filter_magne, field_ctl)
        call add_phys_name_ctl(fhd_w_filter_current, field_ctl)
!
        call add_phys_name_ctl(fhd_wide_SGS_Lorentz, field_ctl)
        call add_phys_name_ctl(fhd_Csim_SGS_Lorentz, field_ctl)
      end if
!
      if(SGS_param%iflag_SGS_uxb .gt. id_SGS_none) then
        call add_phys_name_ctl(fhd_w_filter_velo, field_ctl)
        call add_phys_name_ctl(fhd_w_filter_magne, field_ctl)
!
        call add_phys_name_ctl(fhd_wide_SGS_vp_induct, field_ctl)
        call add_phys_name_ctl(fhd_Csim_SGS_induction, field_ctl)
      end if
!
      if(SGS_param%iflag_SGS_gravity .gt. id_SGS_none) then
        call add_phys_name_ctl(fhd_Reynolds_work, field_ctl)
!
        if(fl_prop%iflag_4_gravity .gt. id_turn_OFF) then
          call add_phys_name_ctl(fhd_SGS_buo_flux, field_ctl)
          call add_phys_name_ctl(fhd_Csim_SGS_buoyancy, field_ctl)
        end if
!
        if(fl_prop%iflag_4_composit_buo .gt. id_turn_OFF) then
          call add_phys_name_ctl(fhd_SGS_comp_buo_flux, field_ctl)
          call add_phys_name_ctl(fhd_Csim_SGS_comp_buo, field_ctl)
        end if
      end if
!
!      if(SGS_param%iflag_SGS_gravity .gt. id_SGS_none) then
!        call add_phys_name_ctl(fhd_filter_velo, field_ctl)
!        call add_phys_name_ctl(fhd_filter_temp, field_ctl)
!
!        call add_phys_name_ctl(fhd_SGS_comp_buo, field_ctl)
!      end if
!
      end subroutine add_field_name_dynamic_SGS
!
! -----------------------------------------------------------------------
!
      end module add_sph_SGS_MHD_fld_2_ctl
