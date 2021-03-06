!>@file   set_mean_square_array.f90
!!        module set_mean_square_array
!!
!! @author H. Matsui
!! @date   Programmed in 2002
!! @n      Modified  on Jan., 2013
!!
!
!> @brief addresses for volume integrated data
!!
!!@verbatim
!!      subroutine count_mean_square_values(nod_fld, fem_msq)
!!      subroutine set_mean_square_values                               &
!!     &         (nod_fld, i_rms, j_ave, ifld_msq)
!!        type(phys_data), intent(in) :: nod_fld
!!        type(phys_address), intent(inout) :: i_rms, j_ave
!!        type(mean_square_address), intent(inout) :: ifld_msq
!!        type(mean_square_values), intent(inout) :: fem_msq
!!@endverbatim
!
      module set_mean_square_array
!
      use m_precision
!
      use t_phys_address
      use t_phys_data
      use t_mean_square_values
!
      implicit  none
!
      private :: set_rms_address
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine count_mean_square_values(nod_fld, fem_msq)
!
      use m_phys_labels
      use m_phys_constants
!
      type(phys_data), intent(in) :: nod_fld
      type(mean_square_values), intent(inout) :: fem_msq
!
      integer (kind = kint) :: i, i0, j0
      character(len = kchara) :: field_name
!
      i0 = 0
      j0 = 0
      do i = 1, nod_fld%num_phys
        field_name = nod_fld%phys_name(i)
        if (nod_fld%iflag_monitor(i) .eq. 1) then
          if     (field_name .eq. fhd_velo                              &
     &       .or. field_name .eq. fhd_filter_velo                       &
     &      ) then
         i0 = i0 + 2
         j0 = j0 + 7
        else if ( field_name .eq. fhd_magne                             &
     &       .or. field_name .eq. fhd_filter_magne                      &
     &      ) then
         i0 = i0 + 3
         j0 = j0 + 7
        else if ( field_name .eq. fhd_current ) then
         i0 = i0 + 4
         j0 = j0 + 6
        else if ( field_name .eq. fhd_vort ) then
         i0 = i0 + 2
         j0 = j0 + 3
        else if ( field_name .eq. fhd_vecp                              &
     &       .or. field_name .eq. fhd_filter_vecp                       &
     &       .or. field_name .eq. fhd_temp                              &
     &       .or. field_name .eq. fhd_part_temp                         &
     &       .or. field_name .eq. fhd_filter_temp                       &
     &       .or. field_name .eq. fhd_filter_comp                       &
     &       .or. field_name .eq. fhd_filter_part_temp                  &
     &       .or. field_name .eq. fhd_press                             &
     &       .or. field_name .eq. fhd_mag_potential                     &
     &       .or. field_name .eq. fhd_light                             &
     &       .or. field_name .eq. fhd_part_light                        &
     &       .or. field_name .eq. fhd_entropy                           &
     &       .or. field_name .eq. fhd_per_entropy                       &
     &       .or. field_name .eq. fhd_heat_source                       &
     &       .or. field_name .eq. fhd_light_source                      &
     &       .or. field_name .eq. fhd_entropy_source                    &
     &       .or. field_name .eq. fhd_density                           &
     &       .or. field_name .eq. fhd_per_density                       &
     &       .or. field_name .eq. fhd_mag_ene_gen                       &
     &       .or. field_name .eq. fhd_work_agst_Lorentz                 &
     &       .or. field_name .eq. fhd_Lorentz_work                      &
     &       .or. field_name .eq. fhd_mag_tension_work                  &
     &       .or. field_name .eq. fhd_buoyancy_flux                     &
     &       .or. field_name .eq. fhd_buoyancy_work                     &
     &       .or. field_name .eq. fhd_comp_buo_flux                     &
     &       .or. field_name .eq. fhd_filter_buo_flux                   &
     &       .or. field_name .eq. fhd_heat_advect                       &
     &       .or. field_name .eq. fhd_part_h_advect                     &
     &       .or. field_name .eq. fhd_part_c_advect                     &
     &       .or. field_name .eq. fhd_div_h_flux                        &
     &       .or. field_name .eq. fhd_div_ph_flux                       &
     &       .or. field_name .eq. fhd_div_c_flux                        &
     &       .or. field_name .eq. fhd_div_pc_flux                       &
     &       .or. field_name .eq. fhd_temp_generation                   &
     &       .or. field_name .eq. fhd_part_temp_gen                     &
     &       .or. field_name .eq. fhd_part_comp_gen                     &
     &       .or. field_name .eq. fhd_div_SGS_h_flux                    &
     &       .or. field_name .eq. fhd_SGS_temp_gen                      &
     &       .or. field_name .eq. fhd_SGS_m_ene_gen                     &
     &       .or. field_name .eq. fhd_SGS_Lorentz_work                  &
     &       .or. field_name .eq. fhd_Reynolds_work                     &
     &       .or. field_name .eq. fhd_SGS_buo_flux                      &
     &       .or. field_name .eq. fhd_SGS_comp_buo_flux                 &
     &       .or. field_name .eq. fhd_c_diffuse                         &
     &       .or. field_name .eq. fhd_thermal_diffusion                 &
     &       .or. field_name .eq. fhd_vis_ene_diffuse                   &
     &       .or. field_name .eq. fhd_mag_ene_diffuse                   &
     &       .or. field_name .eq. fhd_SGS_div_h_flux_true               &
     &       .or. field_name .eq. fhd_SGS_div_c_flux_true               &
     &       .or. field_name .eq. fhd_SGS_Lorentz_wk_true               &
     &       .or. field_name .eq. fhd_Reynolds_work_true                &
     &       .or. field_name .eq. fhd_SGS_temp_gen_true                 &
     &       .or. field_name .eq. fhd_SGS_comp_gen_true                 &
     &       .or. field_name .eq. fhd_SGS_m_ene_gen_true                &
     &       .or. field_name .eq. fhd_div_SGS_h_flux                    &
     &       .or. field_name .eq. fhd_div_SGS_c_flux                    &
     &       .or. field_name .eq. fhd_SGS_div_inertia                   &
     &       .or. field_name .eq. fhd_SGS_div_Lorentz                   &
     &       .or. field_name .eq. fhd_w_filter_temp                     &
     &       .or. field_name .eq. fhd_w_filter_comp                     &
     &      ) then
         i0 = i0 + 1
         j0 = j0 + n_scalar
        else if ( field_name .eq. fhd_mag_tension                       &
     &       .or. field_name .eq. fhd_inertia                           &
     &       .or. field_name .eq. fhd_div_m_flux                        &
     &       .or. field_name .eq. fhd_div_maxwell_t                     &
     &       .or. field_name .eq. fhd_div_induct_t                      &
     &       .or. field_name .eq. fhd_mag_induct                        &
     &       .or. field_name .eq. fhd_vp_induct                         &
     &       .or. field_name .eq. fhd_press_grad                        &
     &       .or. field_name .eq. fhd_mag_stretch                       &
     &       .or. field_name .eq. fhd_Lorentz                           &
     &       .or. field_name .eq. fhd_Coriolis                          &
     &       .or. field_name .eq. fhd_buoyancy                          &
     &       .or. field_name .eq. fhd_comp_buo                          &
     &       .or. field_name .eq. fhd_filter_buo                        &
     &       .or. field_name .eq. fhd_h_flux                            &
     &       .or. field_name .eq. fhd_ph_flux                           &
     &       .or. field_name .eq. fhd_c_flux                            &
     &       .or. field_name .eq. fhd_pc_flux                           &
     &       .or. field_name .eq. fhd_e_field                           &
     &       .or. field_name .eq. fhd_poynting                          &
     &       .or. field_name .eq. fhd_grad_v_1                          &
     &       .or. field_name .eq. fhd_grad_v_2                          &
     &       .or. field_name .eq. fhd_grad_v_3                          &
     &       .or. field_name .eq. fhd_filter_vort                       &
     &       .or. field_name .eq. fhd_w_filter_velo                     &
     &       .or. field_name .eq. fhd_w_filter_vort                     &
     &       .or. field_name .eq. fhd_w_filter_magne                    &
     &       .or. field_name .eq. fhd_w_filter_current                  &
     &      ) then
         i0 = i0 + 1
         j0 = j0 + n_vector
        else if(  field_name .eq. fhd_SGS_h_flux                        &
     &       .or. field_name .eq. fhd_SGS_c_flux                        &
     &       .or. field_name .eq. fhd_div_SGS_m_flux                    &
     &       .or. field_name .eq. fhd_SGS_Lorentz                       &
     &       .or. field_name .eq. fhd_SGS_inertia                       &
     &       .or. field_name .eq. fhd_SGS_induction                     &
     &       .or. field_name .eq. fhd_SGS_vp_induct                     &
     &       .or. field_name .eq. fhd_viscous                           &
     &       .or. field_name .eq. fhd_vecp_diffuse                      &
     &       .or. field_name .eq. fhd_mag_diffuse                       &
     &       .or. field_name .eq. fhd_SGS_div_m_flux_true               &
     &       .or. field_name .eq. fhd_SGS_Lorentz_true                  &
     &       .or. field_name .eq. fhd_SGS_mag_induct_true               &
     &       .or. field_name .eq. fhd_SGS_buoyancy                      &
     &       .or. field_name .eq. fhd_SGS_comp_buo                      &
     &       .or. field_name .eq. fhd_wide_SGS_h_flux                   &
     &       .or. field_name .eq. fhd_wide_SGS_c_flux                   &
     &       .or. field_name .eq. fhd_wide_SGS_inertia                  &
     &       .or. field_name .eq. fhd_wide_SGS_Lorentz                  &
     &       .or. field_name .eq. fhd_wide_SGS_vp_induct                &
     &       .or. field_name .eq. fhd_SGS_rot_inertia                   &
     &       .or. field_name .eq. fhd_SGS_rot_Lorentz                   &
     &       .or. field_name .eq. fhd_geostrophic                       &
     &       .or. field_name .eq. fhd_h_flux_w_sgs                      &
     &       .or. field_name .eq. fhd_c_flux_w_sgs                      &
     &       .or. field_name .eq. fhd_inertia_w_sgs                     &
     &       .or. field_name .eq. fhd_Lorentz_w_sgs                     &
     &       .or. field_name .eq. fhd_vp_induct_w_sgs                   &
     &       .or. field_name .eq. fhd_mag_induct_w_sgs                  &
     &      ) then
         i0 = i0 + 1
         j0 = j0 + n_vector
        else if ( field_name .eq. fhd_mom_flux                          &
     &       .or. field_name .eq. fhd_maxwell_t                         &
     &       .or. field_name .eq. fhd_induct_t                          &
     &       .or. field_name .eq. fhd_SGS_m_flux                        &
     &       .or. field_name .eq. fhd_SGS_maxwell_t                     &
     &       .or. field_name .eq. fhd_SGS_induct_t                      &
     &       .or. field_name .eq. fhd_mom_flux_w_sgs                    &
     &       .or. field_name .eq. fhd_maxwell_t_w_sgs                   &
     &      ) then
         i0 = i0 + 1
         j0 = j0 + n_sym_tensor
!
        else if ( field_name .eq. fhd_velocity_scale                    &
     &       .or. field_name .eq. fhd_magnetic_scale                    &
     &       .or. field_name .eq. fhd_temp_scale                        &
     &       .or. field_name .eq. fhd_composition_scale                 &
     &       .or. field_name .eq. fhd_Csim_SGS_h_flux                   &
     &       .or. field_name .eq. fhd_Csim_SGS_c_flux                   &
     &       .or. field_name .eq. fhd_Csim_SGS_m_flux                   &
     &       .or. field_name .eq. fhd_Csim_SGS_Lorentz                  &
     &       .or. field_name .eq. fhd_Csim_SGS_induction                &
     &       .or. field_name .eq. fhd_Csim_SGS_buoyancy                 &
     &       .or. field_name .eq. fhd_Csim_SGS_comp_buo                 &
     &      ) then
         i0 = i0 + 1
         j0 = j0 + n_scalar
        end if
!
        else
          if    ( field_name .eq. fhd_velo                              &
     &       .or. field_name .eq. fhd_filter_velo                       &
     &       .or. field_name .eq. fhd_magne                             &
     &       .or. field_name .eq. fhd_filter_magne                      &
     &       .or. field_name .eq. fhd_vecp                              &
     &       .or. field_name .eq. fhd_filter_vecp                       &
     &       .or. field_name .eq. fhd_mag_potential                     &
     &      ) then
            i0 = i0 + 1
            j0 = j0 + n_scalar
          end if
        end if
!
      end do
!
      fem_msq%num_rms = i0 + 1
      fem_msq%num_ave = j0
!
       return
       end subroutine count_mean_square_values
!
!-----------------------------------------------------------------------
!
      subroutine set_mean_square_values                                 &
     &         (nod_fld, i_rms, j_ave, ifld_msq)
!
      use m_phys_labels
      use m_phys_constants
!
      type(phys_data), intent(in) :: nod_fld
!
      type(phys_address), intent(inout) :: i_rms, j_ave
      type(mean_square_address), intent(inout) :: ifld_msq
!
      integer (kind = kint) :: i, i0, j0, num_comps
      character(len = kchara) :: field_name
!
!
      i0 = 0
      j0 = 0
      do i = 1, nod_fld%num_phys
        field_name = nod_fld%phys_name(i)
        num_comps =  nod_fld%num_component(i)
        if (nod_fld%iflag_monitor(i) .eq. 1) then
          if ( field_name .eq. fhd_velo) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_velo, j_ave%i_velo)
            call set_rms_address(n_scalar, i0, j0,                      &
     &          i_rms%i_div_v, j_ave%i_div_v)
!
            ifld_msq%ja_amom = j0 + 1
            j0 = j0 + 3
          end if
!
          if ( field_name .eq. fhd_magne ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_magne, j_ave%i_magne)
            call set_rms_address(num_comps, i0, j0,                     &
     &          ifld_msq%ir_me_ic, ifld_msq%ja_mag_ic)
            call set_rms_address(n_scalar, i0, j0,                      &
     &          i_rms%i_div_b, j_ave%i_div_b)
          end if
!
          if ( field_name .eq. fhd_vecp ) then
            call set_rms_address(n_scalar, i0, j0,                      &
     &          i_rms%i_div_a, j_ave%i_div_a)
          end if

          if ( field_name .eq. fhd_vort ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_vort, j_ave%i_vort)
            ifld_msq%ir_rms_w   = i0 + 1
            i0 = i0 + 1
          end if
!
          if ( field_name .eq. fhd_current ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_current, j_ave%i_current)
            call set_rms_address(num_comps, i0, j0,                     &
     &          ifld_msq%ir_sqj_ic, ifld_msq%ja_j_ic)
!
            ifld_msq%ir_rms_j    = i0 + 1
            ifld_msq%ir_rms_j_ic = i0 + 2
            i0 = i0 + 2
          end if
!
          if ( field_name .eq. fhd_e_field ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_electric, j_ave%i_electric)
          else if ( field_name .eq. fhd_poynting ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_poynting, j_ave%i_poynting)
          else if ( field_name .eq. fhd_temp ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_temp, j_ave%i_temp)
          else if ( field_name .eq. fhd_press ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_press, j_ave%i_press)
          else if ( field_name .eq. fhd_mag_potential ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_mag_p, j_ave%i_mag_p)
          end if
!
          if ( field_name .eq. fhd_part_temp ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_par_temp, j_ave%i_par_temp)
          else if ( field_name .eq. fhd_light ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_light, j_ave%i_light)
          else if ( field_name .eq. fhd_part_light ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_par_light, j_ave%i_par_light)
          else if ( field_name .eq. fhd_entropy ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_entropy, j_ave%i_entropy)
          else if ( field_name .eq. fhd_per_entropy ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_par_entropy, j_ave%i_par_entropy)
          else if ( field_name .eq. fhd_density ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_density, j_ave%i_density)
          else if ( field_name .eq. fhd_per_density ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_par_density, j_ave%i_par_density)
!
          else if ( field_name .eq. fhd_heat_source ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_heat_source, j_ave%i_heat_source)
          else if ( field_name .eq. fhd_light_source ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_light_source, j_ave%i_light_source)
          else if ( field_name .eq. fhd_entropy_source ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_entropy_source, j_ave%i_entropy_source)
          end if
!
          if ( field_name .eq. fhd_press_grad ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_press_grad, j_ave%i_press_grad)
          end if
!
          if ( field_name .eq. fhd_mag_tension ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_m_tension, j_ave%i_m_tension)
          end if
!
          if ( field_name .eq. fhd_filter_velo ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_filter_velo, j_ave%i_filter_velo)
            call set_rms_address(n_scalar, i0, j0,                      &
     &          i_rms%i_div_filter_v, j_ave%i_div_filter_v)
            ifld_msq%jr_amom_f = i0 + 1
            j0 = j0 + 3
          end if
!
          if ( field_name .eq. fhd_filter_magne ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_filter_magne, j_ave%i_filter_magne)
            call set_rms_address(num_comps, i0, j0,                     &
     &          ifld_msq%ir_me_f_ic, ifld_msq%ja_mag_f_ic)
            call set_rms_address(n_scalar, i0, j0,                      &
     &          i_rms%i_div_filter_b, j_ave%i_div_filter_b)
          end if
!
          if ( field_name .eq. fhd_filter_vecp ) then
            call set_rms_address(n_scalar, i0, j0,                      &
     &          i_rms%i_div_filter_a, j_ave%i_div_filter_a)
          else if ( field_name .eq. fhd_filter_temp ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_filter_temp, j_ave%i_filter_temp)
          else if ( field_name .eq. fhd_filter_comp ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_filter_comp, j_ave%i_filter_comp)
          else if ( field_name .eq. fhd_w_filter_temp ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_wide_fil_temp, j_ave%i_wide_fil_temp)
          else if ( field_name .eq. fhd_w_filter_comp ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_wide_fil_comp, j_ave%i_wide_fil_comp)
          end if
!
          if ( field_name .eq. fhd_grad_v_1 ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_grad_vx, j_ave%i_grad_vx)
          else if ( field_name .eq. fhd_grad_v_2 ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_grad_vy, j_ave%i_grad_vy)
          else if ( field_name .eq. fhd_grad_v_3 ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_grad_vz, j_ave%i_grad_vz)
          end if
!
          if ( field_name .eq. fhd_mom_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_m_flux, j_ave%i_m_flux)
          else if ( field_name .eq. fhd_maxwell_t ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_maxwell, j_ave%i_maxwell)
          else if ( field_name .eq. fhd_induct_t ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_induct_t, j_ave%i_induct_t)
          else if ( field_name .eq. fhd_inertia ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_m_advect, j_ave%i_m_advect)
          else if ( field_name .eq. fhd_div_m_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_m_flux_div, j_ave%i_m_flux_div)
          else if ( field_name .eq. fhd_div_maxwell_t ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_maxwell_div, j_ave%i_maxwell_div)
          else if ( field_name .eq. fhd_div_induct_t ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_induct_div, j_ave%i_induct_div)
          else if ( field_name .eq. fhd_mag_induct ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_induction, j_ave%i_induction)
          else if ( field_name .eq. fhd_vp_induct ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_vp_induct, j_ave%i_vp_induct)
          else if ( field_name .eq. fhd_mag_stretch ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_mag_stretch, j_ave%i_mag_stretch)
          else if ( field_name .eq. fhd_Lorentz ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_lorentz, j_ave%i_lorentz)
          else if ( field_name .eq. fhd_Coriolis ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_coriolis, j_ave%i_coriolis)
          else if ( field_name .eq. fhd_buoyancy ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_buoyancy, j_ave%i_buoyancy)
          else if ( field_name .eq. fhd_comp_buo ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_comp_buo, j_ave%i_comp_buo)
          else if ( field_name .eq. fhd_filter_buo ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_filter_buo, j_ave%i_filter_buo)
          end if
!
          if ( field_name .eq. fhd_viscous ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_v_diffuse, j_ave%i_v_diffuse)
          else if ( field_name .eq. fhd_vecp_diffuse ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_vp_diffuse, j_ave%i_vp_diffuse)
          else if ( field_name .eq. fhd_mag_diffuse ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_b_diffuse, j_ave%i_b_diffuse)
          else if ( field_name .eq. fhd_thermal_diffusion ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_t_diffuse, j_ave%i_t_diffuse)
          else if ( field_name .eq. fhd_c_diffuse) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_c_diffuse, j_ave%i_c_diffuse)
          end if
!
          if ( field_name .eq. fhd_SGS_m_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_m_flux, j_ave%i_SGS_m_flux)
          else if ( field_name .eq. fhd_SGS_maxwell_t ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_maxwell, j_ave%i_SGS_maxwell)
          else if ( field_name .eq. fhd_SGS_induct_t ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_induct_t, j_ave%i_SGS_induct_t)
          else if ( field_name .eq. fhd_div_SGS_m_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_div_m_flux, j_ave%i_SGS_div_m_flux)
          else if ( field_name .eq. fhd_mom_flux_w_sgs ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_mom_flux_w_sgs, j_ave%i_mom_flux_w_sgs)
          else if ( field_name .eq. fhd_maxwell_t_w_sgs ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_maxwell_t_w_sgs, j_ave%i_maxwell_t_w_sgs)
          else if ( field_name .eq. fhd_SGS_Lorentz ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_Lorentz, j_ave%i_SGS_Lorentz)
          else if ( field_name .eq. fhd_SGS_inertia ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_inertia, j_ave%i_SGS_inertia)
          else if ( field_name .eq. fhd_SGS_induction ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_induction, j_ave%i_SGS_induction)
          else if ( field_name .eq. fhd_SGS_vp_induct ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_vp_induct, j_ave%i_SGS_vp_induct)
          else if ( field_name .eq. fhd_SGS_buoyancy ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_buoyancy, j_ave%i_SGS_buoyancy)
          else if ( field_name .eq. fhd_SGS_comp_buo ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_comp_buo, j_ave%i_SGS_comp_buo)
          end if
!
          if ( field_name .eq. fhd_wide_SGS_h_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_wide_SGS_h_flux, j_ave%i_wide_SGS_h_flux)
          else if ( field_name .eq. fhd_wide_SGS_c_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_wide_SGS_c_flux, j_ave%i_wide_SGS_c_flux)
          else if ( field_name .eq. fhd_wide_SGS_inertia ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_wide_SGS_inertia, j_ave%i_wide_SGS_inertia)
          else if ( field_name .eq. fhd_wide_SGS_Lorentz ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_wide_SGS_Lorentz, j_ave%i_wide_SGS_Lorentz)
          else if ( field_name .eq. fhd_wide_SGS_vp_induct ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_wide_SGS_vp_induct, j_ave%i_wide_SGS_vp_induct)
          else if ( field_name .eq. fhd_SGS_rot_inertia ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_rot_inertia, j_ave%i_SGS_rot_inertia)
          else if ( field_name .eq. fhd_SGS_rot_Lorentz ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_rot_Lorentz, j_ave%i_SGS_rot_Lorentz)
          else if ( field_name .eq. fhd_geostrophic ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_geostrophic, j_ave%i_geostrophic)
          else if ( field_name .eq. fhd_h_flux_w_sgs ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_h_flux_w_sgs, j_ave%i_h_flux_w_sgs)
          else if ( field_name .eq. fhd_c_flux_w_sgs ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_c_flux_w_sgs, j_ave%i_c_flux_w_sgs)
          else if ( field_name .eq. fhd_inertia_w_sgs ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_inertia_w_sgs, j_ave%i_inertia_w_sgs)
          else if ( field_name .eq. fhd_Lorentz_w_sgs ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_Lorentz_w_sgs, j_ave%i_Lorentz_w_sgs)
          else if ( field_name .eq. fhd_vp_induct_w_sgs ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_vp_induct_w_sgs, j_ave%i_vp_induct_w_sgs)
          else if ( field_name .eq. fhd_mag_induct_w_sgs ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_mag_induct_w_sgs, j_ave%i_mag_induct_w_sgs)
          end if
!
          if ( field_name .eq. fhd_mag_ene_gen ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_me_gen, j_ave%i_me_gen)
          else if ( field_name .eq. fhd_Lorentz_work ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_ujb, j_ave%i_ujb)
          else if ( field_name .eq. fhd_work_agst_Lorentz ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_nega_ujb, j_ave%i_nega_ujb)
          else if ( field_name .eq. fhd_mag_tension_work ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_m_tension_wk, j_ave%i_m_tension_wk)
          else if ( field_name .eq. fhd_buoyancy_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_buo_gen, j_ave%i_buo_gen)
          else if ( field_name .eq. fhd_comp_buo_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_c_buo_gen, j_ave%i_c_buo_gen)
          else if ( field_name .eq. fhd_filter_buo_flux) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_f_buo_gen, j_ave%i_f_buo_gen)
          end if
!
          if ( field_name .eq. fhd_vis_ene_diffuse ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_vis_e_diffuse, j_ave%i_vis_e_diffuse)
          else if ( field_name .eq. fhd_mag_ene_diffuse ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_mag_e_diffuse, j_ave%i_mag_e_diffuse)
          else if ( field_name .eq. fhd_h_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_h_flux, j_ave%i_h_flux)
          else if ( field_name .eq. fhd_ph_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_ph_flux, j_ave%i_ph_flux)
          else if ( field_name .eq. fhd_c_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_c_flux, j_ave%i_c_flux)
          else if ( field_name .eq. fhd_pc_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_pc_flux, j_ave%i_pc_flux)
          else if ( field_name .eq. fhd_heat_advect ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_h_advect, j_ave%i_h_advect)
          else if ( field_name .eq. fhd_part_h_advect ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_ph_advect, j_ave%i_ph_advect)
          else if ( field_name .eq. fhd_part_c_advect ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_pc_advect, j_ave%i_pc_advect)
          else if ( field_name .eq. fhd_div_h_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_h_flux_div, j_ave%i_h_flux_div)
          else if ( field_name .eq. fhd_div_ph_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_ph_flux_div, j_ave%i_ph_flux_div)
          else if ( field_name .eq. fhd_div_c_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_c_flux_div, j_ave%i_c_flux_div)
          else if ( field_name .eq. fhd_div_pc_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_pc_flux_div, j_ave%i_pc_flux_div)
          else if ( field_name .eq. fhd_temp_generation ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_temp_gen, j_ave%i_temp_gen)
          else if ( field_name .eq. fhd_part_temp_gen ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_par_t_gen, j_ave%i_par_t_gen)
          else if ( field_name .eq. fhd_part_comp_gen ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_par_c_gen, j_ave%i_par_c_gen)
          else if ( field_name .eq. fhd_SGS_h_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_h_flux, j_ave%i_SGS_h_flux)
          else if ( field_name .eq. fhd_div_SGS_h_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_div_h_flux, j_ave%i_SGS_div_h_flux)
          else if ( field_name .eq. fhd_SGS_c_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_c_flux, j_ave%i_SGS_c_flux)
          else if ( field_name .eq. fhd_SGS_temp_gen ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_temp_gen, j_ave%i_SGS_temp_gen)
          else if ( field_name .eq. fhd_SGS_m_ene_gen ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_me_gen, j_ave%i_SGS_me_gen)
          else if ( field_name .eq. fhd_SGS_Lorentz_work ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_Lor_wk, j_ave%i_SGS_Lor_wk)
          else if ( field_name .eq. fhd_Reynolds_work ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_reynolds_wk, j_ave%i_reynolds_wk)
          else if ( field_name .eq. fhd_SGS_buo_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_buo_wk, j_ave%i_SGS_buo_wk)
          else if ( field_name .eq. fhd_SGS_comp_buo_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_comp_buo_wk, j_ave%i_SGS_comp_buo_wk)
          end if
!
          if (field_name .eq. fhd_SGS_div_h_flux_true) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_div_hf_true, j_ave%i_SGS_div_hf_true)
          else if (field_name .eq. fhd_SGS_div_c_flux_true) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_div_cf_true, j_ave%i_SGS_div_cf_true)
          else if (field_name .eq. fhd_SGS_div_m_flux_true) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_div_mf_true, j_ave%i_SGS_div_mf_true)
          else if ( field_name .eq. fhd_SGS_Lorentz_true ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_Lor_true, j_ave%i_SGS_Lor_true)
          else if ( field_name .eq. fhd_SGS_mag_induct_true ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_idct_true, j_ave%i_SGS_idct_true)
          end if
!
          if ( field_name .eq. fhd_SGS_Lorentz_wk_true ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_Lor_wk_tr, j_ave%i_SGS_Lor_wk_tr)
          else if ( field_name .eq. fhd_Reynolds_work_true ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_reynolds_wk_tr, j_ave%i_reynolds_wk_tr)
          else if ( field_name .eq. fhd_SGS_temp_gen_true ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_t_gen_tr, j_ave%i_SGS_t_gen_tr)
          else if ( field_name .eq. fhd_SGS_comp_gen_true ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_c_gen_tr, j_ave%i_SGS_c_gen_tr)
          else if ( field_name .eq. fhd_SGS_m_ene_gen_true ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_me_gen_tr, j_ave%i_SGS_me_gen_tr)
          else if ( field_name .eq. fhd_div_SGS_h_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_div_h_flux, j_ave%i_SGS_div_h_flux)
          else if ( field_name .eq. fhd_div_SGS_c_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_div_c_flux, j_ave%i_SGS_div_c_flux)
          else if ( field_name .eq. fhd_SGS_div_inertia ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_div_inertia, j_ave%i_SGS_div_inertia)
          else if ( field_name .eq. fhd_SGS_div_Lorentz ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_SGS_div_Lorentz, j_ave%i_SGS_div_Lorentz)
          end if
!
          if ( field_name .eq. fhd_w_filter_velo ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_wide_fil_velo, j_ave%i_wide_fil_velo)
          else if ( field_name .eq. fhd_filter_vort ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_wide_fil_vort, j_ave%i_wide_fil_vort)
          else if ( field_name .eq. fhd_w_filter_vort ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_wide_fil_vort, j_ave%i_wide_fil_vort)
          else if ( field_name .eq. fhd_w_filter_magne ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_wide_fil_magne, j_ave%i_wide_fil_magne)
          else if ( field_name .eq. fhd_w_filter_current ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_wide_fil_current, j_ave%i_wide_fil_current)
          end if
!
          if ( field_name .eq. fhd_velocity_scale ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_velo_scale, j_ave%i_velo_scale)
          else if ( field_name .eq. fhd_magnetic_scale ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_magne_scale, j_ave%i_magne_scale)
          else if ( field_name .eq. fhd_temp_scale ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_temp_scale, j_ave%i_temp_scale)
          else if ( field_name .eq. fhd_composition_scale ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_comp_scale, j_ave%i_comp_scale)
          end if
!
          if ( field_name .eq. fhd_Csim_SGS_h_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_Csim_SGS_h_flux, j_ave%i_Csim_SGS_h_flux)
          else if ( field_name .eq. fhd_Csim_SGS_c_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_Csim_SGS_c_flux, j_ave%i_Csim_SGS_c_flux)
          else if ( field_name .eq. fhd_Csim_SGS_m_flux ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_Csim_SGS_m_flux, j_ave%i_Csim_SGS_m_flux)
          else if ( field_name .eq. fhd_Csim_SGS_Lorentz ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_Csim_SGS_Lorentz, j_ave%i_Csim_SGS_Lorentz)
          else if ( field_name .eq. fhd_Csim_SGS_induction ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_Csim_SGS_induction, j_ave%i_Csim_SGS_induction)
          else if ( field_name .eq. fhd_Csim_SGS_buoyancy ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_Csim_SGS_buoyancy, j_ave%i_Csim_SGS_buoyancy)
          else if ( field_name .eq. fhd_Csim_SGS_comp_buo ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_Csim_SGS_comp_buo, j_ave%i_Csim_SGS_comp_buo)
          end if
!
!   Old field label... Should be deleted later!!
          if ( field_name .eq. fhd_buoyancy_work ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_buo_gen, j_ave%i_buo_gen)
          end if
!
        else
          if ( field_name .eq. fhd_velo) then
            call set_rms_address(n_scalar, i0, j0,                      &
     &          i_rms%i_div_v, j_ave%i_div_v)
          else if ( field_name .eq. fhd_magne ) then
            call set_rms_address(n_scalar, i0, j0,                      &
     &          i_rms%i_div_b, j_ave%i_div_b)
          else if ( field_name .eq. fhd_vecp ) then
            call set_rms_address(n_scalar, i0, j0,                      &
     &          i_rms%i_div_a, j_ave%i_div_a)
          else if ( field_name .eq. fhd_filter_velo ) then
            call set_rms_address(n_scalar, i0, j0,                      &
     &          i_rms%i_div_filter_v, j_ave%i_div_filter_v)
          else if ( field_name .eq. fhd_filter_magne ) then
            call set_rms_address(n_scalar, i0, j0,                      &
     &          i_rms%i_div_filter_b, j_ave%i_div_filter_b)
          else if ( field_name .eq. fhd_filter_vecp ) then
            call set_rms_address(n_scalar, i0, j0,                      &
     &          i_rms%i_div_filter_a, j_ave%i_div_filter_a)
          else if ( field_name .eq. fhd_mag_potential ) then
            call set_rms_address(num_comps, i0, j0,                     &
     &          i_rms%i_mag_p, j_ave%i_mag_p)
          end if

        end if
      end do
!
      ifld_msq%ivol = i0 + 1
      i0 = i0 + 1
!
      end subroutine set_mean_square_values
!
! ----------------------------------------------------------------------
!
      subroutine set_rms_address(numdir, i0, j0, ir_rms, ja_ave)
!
      integer(kind = kint), intent(in) :: numdir
      integer(kind = kint), intent(inout) :: ir_rms, ja_ave, i0, j0
!
      ir_rms = i0 + 1
      ja_ave = j0 + 1
!
      i0 = i0 + 1
      j0 = j0 + numdir
!
      end subroutine set_rms_address
!
! ----------------------------------------------------------------------
!
      end module set_mean_square_array
