!>@file   cal_momentum_eq_explicit.f90
!!@brief  module cal_momentum_eq_explicit
!!
!!@author H. Matsui
!!@date    programmed by H.Matsui in Oct., 2010
!
!>@brief Time integration for momentum equation by explicit scheme
!!
!!@verbatim
!!      subroutine cal_expricit_sph_adams(evo_V, evo_B, evo_T, evo_C,   &
!!     &          sph_rj, fl_prop,cd_prop,  ht_prop, cp_prop,           &
!!     &         ipol, itor, rj_fld)
!!     &          sph_rj, ht_prop, cp_prop, ipol, itor, rj_fld)
!!      subroutine cal_expricit_sph_euler                               &
!!     &        (i_step, evo_V, evo_B, evo_T, evo_C, sph_rj,            &
!!     &         fl_prop, cd_prop, ht_prop, cp_prop, ipol, itor, rj_fld)
!!        type(time_evolution_params), intent(in) :: evo_V, evo_B
!!        type(time_evolution_params), intent(in) :: evo_T, evo_C
!!        type(sph_rj_grid), intent(in) ::  sph_rj
!!        type(fdm_matrices), intent(in) :: r_2nd
!!        type(scalar_property), intent(in) :: fl_prop
!!        type(conductive_property), intent(in) :: cd_prop
!!        type(scalar_property), intent(in) :: ht_prop, cp_prop
!!        type(legendre_4_sph_trans), intent(in) :: leg
!!        type(phys_address), intent(in) :: ipol, itor
!!        type(phys_data), intent(inout) :: rj_fld
!!@endverbatim
!!
!!@param i_step  time step
!
      module cal_momentum_eq_explicit
!
      use m_precision
!
      use t_time_stepping_parameter
      use t_physical_property
      use t_spheric_rj_data
      use t_phys_address
      use t_phys_data
      use t_fdm_coefs
      use t_schmidt_poly_on_rtm
!
      implicit  none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine cal_expricit_sph_adams(evo_V, evo_B, evo_T, evo_C,     &
     &          sph_rj, fl_prop,cd_prop,  ht_prop, cp_prop,             &
     &         ipol, itor, rj_fld)
!
      use m_boundary_params_sph_MHD
      use cal_explicit_terms
      use cal_vorticity_terms_adams
      use cal_nonlinear_sph_MHD
      use select_diff_adv_source
!
      type(time_evolution_params), intent(in) :: evo_V, evo_B
      type(time_evolution_params), intent(in) :: evo_T, evo_C
      type(sph_rj_grid), intent(in) ::  sph_rj
      type(fluid_property), intent(in) :: fl_prop
      type(conductive_property), intent(in) :: cd_prop
      type(scalar_property), intent(in) :: ht_prop, cp_prop
      type(phys_address), intent(in) :: ipol, itor
      type(phys_data), intent(inout) :: rj_fld
!
!
!$omp parallel
      if(evo_V%iflag_scheme .gt.     id_no_evolution) then
        call cal_vorticity_eq_adams(ipol, itor,                         &
     &      sph_bc_U%kr_in, sph_bc_U%kr_out, fl_prop%coef_exp,          &
     &      rj_fld%n_point,sph_rj%nidx_rj(2), rj_fld%ntot_phys,         &
     &      rj_fld%d_fld)
      end if
!
      if(evo_B%iflag_scheme .gt.    id_no_evolution) then
        call cal_diff_induction_MHD_adams(cd_prop%coef_exp,             &
     &      ipol, itor, rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
      end if
      if(evo_T%iflag_scheme .gt.     id_no_evolution) then
        call sel_scalar_diff_adv_src_adams                              &
     &     (sph_bc_T%kr_in, sph_bc_T%kr_out,                            &
     &      ipol%i_t_diffuse, ipol%i_h_advect, ipol%i_heat_source,      &
     &      ipol%i_temp, ipol%i_pre_heat, ht_prop%coef_exp,             &
     &      ht_prop%coef_source, sph_rj, rj_fld)
      end if
      if(evo_C%iflag_scheme .gt. id_no_evolution) then
        call sel_scalar_diff_adv_src_adams                              &
     &     (sph_bc_C%kr_in, sph_bc_C%kr_out,                            &
     &      ipol%i_c_diffuse, ipol%i_c_advect, ipol%i_light_source,     &
     &      ipol%i_light, ipol%i_pre_composit, cp_prop%coef_exp,        &
     &      cp_prop%coef_source, sph_rj, rj_fld)
      end if
!$omp end parallel
!
      end subroutine cal_expricit_sph_adams
!
! ----------------------------------------------------------------------
!
      subroutine cal_expricit_sph_euler                                 &
     &        (i_step, evo_V, evo_B, evo_T, evo_C, sph_rj,              &
     &         fl_prop, cd_prop, ht_prop, cp_prop, ipol, itor, rj_fld)
!
      use m_boundary_params_sph_MHD
      use cal_explicit_terms
      use cal_vorticity_terms_adams
      use select_diff_adv_source
!
      integer(kind = kint), intent(in) :: i_step
      type(time_evolution_params), intent(in) :: evo_V, evo_B
      type(time_evolution_params), intent(in) :: evo_T, evo_C
      type(sph_rj_grid), intent(in) ::  sph_rj
      type(fluid_property), intent(in) :: fl_prop
      type(conductive_property), intent(in) :: cd_prop
      type(scalar_property), intent(in) :: ht_prop, cp_prop
      type(phys_address), intent(in) :: ipol, itor
      type(phys_data), intent(inout) :: rj_fld
!
!$omp parallel
      if(evo_V%iflag_scheme .gt.     id_no_evolution) then
        call cal_vorticity_eq_euler(ipol, itor,                         &
     &      sph_bc_U%kr_in, sph_bc_U%kr_out, fl_prop%coef_exp,          &
     &      rj_fld%n_point, sph_rj%nidx_rj(2), rj_fld%ntot_phys,        &
     &      rj_fld%d_fld)
      end if
!
      if(evo_T%iflag_scheme .gt.     id_no_evolution) then
        call sel_scalar_diff_adv_src_euler                              &
     &     (sph_bc_T%kr_in, sph_bc_T%kr_out,                            &
     &      ipol%i_t_diffuse, ipol%i_h_advect, ipol%i_heat_source,      &
     &      ipol%i_temp, ht_prop%coef_exp, ht_prop%coef_advect,         &
     &      ht_prop%coef_source, sph_rj, rj_fld)
      end if
      if(evo_B%iflag_scheme .gt.    id_no_evolution) then
        call cal_diff_induction_MHD_euler(cd_prop%coef_exp,             &
     &      ipol, itor, rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
      end if
      if(evo_C%iflag_scheme .gt. id_no_evolution) then
        call sel_scalar_diff_adv_src_euler                              &
     &     (sph_bc_C%kr_in, sph_bc_C%kr_out,                            &
     &      ipol%i_c_diffuse, ipol%i_c_advect, ipol%i_light_source,     &
     &      ipol%i_light, cp_prop%coef_exp, cp_prop%coef_advect,        &
     &      cp_prop%coef_source, sph_rj, rj_fld)
      end if
!
      if (i_step .eq. 1) then
        if(evo_V%iflag_scheme .gt.     id_no_evolution) then
          call set_ini_adams_inertia(ipol, itor,                        &
     &        rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
        end if
        if(evo_T%iflag_scheme .gt.     id_no_evolution) then
          call sel_ini_adams_scalar_w_src                               &
     &       (sph_bc_T%kr_in, sph_bc_T%kr_out, ipol%i_h_advect,         &
     &        ipol%i_heat_source, ipol%i_pre_heat,                      &
     &        ht_prop%coef_source, sph_rj, rj_fld)
        end if
        if(evo_B%iflag_scheme .gt.    id_no_evolution) then
          call set_ini_adams_mag_induct(ipol, itor,                     &
     &        rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
        end if
        if(evo_C%iflag_scheme .gt. id_no_evolution) then
          call sel_ini_adams_scalar_w_src                               &
     &       (sph_bc_C%kr_in, sph_bc_C%kr_out, ipol%i_c_advect,         &
     &        ipol%i_light_source, ipol%i_pre_composit,                 &
     &        cp_prop%coef_source, sph_rj, rj_fld)
        end if
      end if
!$omp end parallel
!
      end subroutine cal_expricit_sph_euler
!
! ----------------------------------------------------------------------
!
      end module cal_momentum_eq_explicit
