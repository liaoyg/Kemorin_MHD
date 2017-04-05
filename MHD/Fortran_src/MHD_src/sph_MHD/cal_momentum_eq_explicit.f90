!>@file   cal_momentum_eq_explicit.f90
!!@brief  module cal_momentum_eq_explicit
!!
!!@author H. Matsui
!!@date    programmed by H.Matsui in Oct., 2010
!
!>@brief Time integration for momentum equation by explicit scheme
!!
!!@verbatim
!!      subroutine cal_explicit_sph_adams                               &
!!     &         (dt, sph_rj, fl_prop, cd_prop, ht_prop, cp_prop,       &
!!     &          ipol, itor, rj_fld)
!!      subroutine cal_explicit_sph_euler(i_step, dt, SGS_param, sph_rj,&
!!     &         fl_prop, cd_prop, ht_prop, cp_prop, ipol, itor, rj_fld)
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
      use t_SGS_control_parameter
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
      subroutine cal_explicit_sph_adams(dt, SGS_param,                  &
     &          sph_rj, fl_prop, cd_prop, ht_prop, cp_prop,             &
     &          ipol, itor, rj_fld)
!
      use m_boundary_params_sph_MHD
      use cal_explicit_terms
      use cal_vorticity_terms_adams
      use cal_nonlinear_sph_MHD
      use select_diff_adv_source
!
      real(kind = kreal), intent(in) :: dt
!
      type(SGS_model_control_params), intent(in) :: SGS_param
      type(sph_rj_grid), intent(in) ::  sph_rj
      type(fluid_property), intent(in) :: fl_prop
      type(conductive_property), intent(in) :: cd_prop
      type(scalar_property), intent(in) :: ht_prop, cp_prop
      type(phys_address), intent(in) :: ipol, itor
      type(phys_data), intent(inout) :: rj_fld
!
!
!$omp parallel
      if(fl_prop%iflag_scheme .gt.     id_no_evolution) then
        call cal_vorticity_eq_adams(ipol, itor,                         &
     &      sph_bc_U%kr_in, sph_bc_U%kr_out, dt, fl_prop%coef_exp,      &
     &      rj_fld%n_point,sph_rj%nidx_rj(2), rj_fld%ntot_phys,         &
     &      rj_fld%d_fld)
      end if
!
      if(cd_prop%iflag_Bevo_scheme .gt.    id_no_evolution) then
        if(SGS_param%iflag_SGS_uxb .gt. id_SGS_none) then
          if(iflag_debug .gt. 0) write(*,*)                             &
     &                'cal_diff_induction_wSGS_adams'
          call cal_diff_induction_wSGS_adams                            &
     &       (dt, cd_prop%coef_exp, ipol, itor,                         &
     &        rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
        else
          if(iflag_debug .gt. 0) write(*,*)                             &
     &                'cal_diff_induction_MHD_adams'
          call cal_diff_induction_MHD_adams                             &
     &       (dt, cd_prop%coef_exp, ipol, itor,                         &
     &        rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
        end if
      end if
!
      if(ht_prop%iflag_scheme .gt.     id_no_evolution) then
          if(iflag_debug .gt. 0) write(*,*)                             &
     &                'sel_scalar_diff_adv_src_adams temperature'
        call sel_scalar_diff_adv_src_adams(SGS_param%iflag_SGS_h_flux,  &
     &      sph_bc_T%kr_in, sph_bc_T%kr_out, ipol%i_t_diffuse,          &
     &      ipol%i_h_advect, ipol%i_SGS_div_h_flux,                     &
     &      ipol%i_heat_source, ipol%i_temp, ipol%i_pre_heat,           &
     &      dt, ht_prop%coef_exp, ht_prop%coef_source, sph_rj, rj_fld)
      end if
!
      if(cp_prop%iflag_scheme .gt. id_no_evolution) then
          if(iflag_debug .gt. 0) write(*,*)                             &
     &                'sel_scalar_diff_adv_src_adams composition'
        call sel_scalar_diff_adv_src_adams(SGS_param%iflag_SGS_c_flux,  &
     &      sph_bc_C%kr_in, sph_bc_C%kr_out, ipol%i_c_diffuse,          &
     &      ipol%i_c_advect, ipol%i_SGS_div_c_flux,                     &
     &      ipol%i_light_source, ipol%i_light, ipol%i_pre_composit,     &
     &      dt, cp_prop%coef_exp, cp_prop%coef_source, sph_rj, rj_fld)
      end if
!$omp end parallel
!
!  Center evolution
!
      if(ht_prop%iflag_scheme .gt.     id_no_evolution) then
        if(iflag_debug .gt. 0) write(*,*)                               &
     &                'sel_ct_scl_diff_adv_src_adams temperature'
        call sel_ct_scl_diff_adv_src_adams                              &
     &     (SGS_param%iflag_SGS_h_flux, ipol%i_t_diffuse,               &
     &      ipol%i_h_advect, ipol%i_SGS_div_h_flux,                     &
     &      ipol%i_heat_source, ipol%i_temp, ipol%i_pre_heat,           &
     &      dt, ht_prop%coef_exp, ht_prop%coef_source, sph_rj, rj_fld)
      end if
!
      if(cp_prop%iflag_scheme .gt. id_no_evolution) then
          if(iflag_debug .gt. 0) write(*,*)                             &
     &                'sel_ct_scl_diff_adv_src_adams composition'
        call sel_ct_scl_diff_adv_src_adams                              &
     &     (SGS_param%iflag_SGS_c_flux, ipol%i_c_diffuse,               &
     &      ipol%i_c_advect, ipol%i_SGS_div_c_flux,                     &
     &      ipol%i_light_source, ipol%i_light, ipol%i_pre_composit,     &
     &      dt, cp_prop%coef_exp, cp_prop%coef_source, sph_rj, rj_fld)
      end if
!
      end subroutine cal_explicit_sph_adams
!
! ----------------------------------------------------------------------
!
      subroutine cal_explicit_sph_euler(i_step, dt, SGS_param, sph_rj,  &
     &         fl_prop, cd_prop, ht_prop, cp_prop, ipol, itor, rj_fld)
!
      use m_boundary_params_sph_MHD
      use cal_explicit_terms
      use cal_vorticity_terms_adams
      use select_diff_adv_source
!
      integer(kind = kint), intent(in) :: i_step
      real(kind = kreal), intent(in) :: dt
!
      type(SGS_model_control_params), intent(in) :: SGS_param
      type(sph_rj_grid), intent(in) ::  sph_rj
      type(fluid_property), intent(in) :: fl_prop
      type(conductive_property), intent(in) :: cd_prop
      type(scalar_property), intent(in) :: ht_prop, cp_prop
      type(phys_address), intent(in) :: ipol, itor
      type(phys_data), intent(inout) :: rj_fld
!
!$omp parallel
      if(fl_prop%iflag_scheme .gt.     id_no_evolution) then
        call cal_vorticity_eq_euler(ipol, itor,                         &
     &      sph_bc_U%kr_in, sph_bc_U%kr_out, dt, fl_prop%coef_exp,      &
     &      rj_fld%n_point, sph_rj%nidx_rj(2), rj_fld%ntot_phys,        &
     &      rj_fld%d_fld)
      end if
!
      if(cd_prop%iflag_Bevo_scheme .gt.    id_no_evolution) then
        if(SGS_param%iflag_SGS_uxb .gt. id_SGS_none) then
          if(iflag_debug .gt. 0) write(*,*)                             &
     &                'cal_diff_induction_wSGS_euler'
          call cal_diff_induction_wSGS_euler                            &
     &       (dt, cd_prop%coef_exp, ipol, itor,                         &
     &        rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
        else
          if(iflag_debug .gt. 0) write(*,*)                             &
     &                'cal_diff_induction_MHD_euler'
          call cal_diff_induction_MHD_euler                             &
     &       (dt, cd_prop%coef_exp, ipol, itor,                         &
     &        rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
        end if
      end if
!
      if(ht_prop%iflag_scheme .gt.     id_no_evolution) then
        if(iflag_debug .gt. 0) write(*,*)                               &
     &                'sel_scalar_diff_adv_src_euler temperature'
        call sel_scalar_diff_adv_src_euler(SGS_param%iflag_SGS_h_flux,  &
     &      sph_bc_T%kr_in, sph_bc_T%kr_out, ipol%i_t_diffuse,          &
     &      ipol%i_h_advect, ipol%i_SGS_div_h_flux,                     &
     &      ipol%i_heat_source, ipol%i_temp, dt, ht_prop%coef_exp,      &
     &      ht_prop%coef_advect, ht_prop%coef_source, sph_rj, rj_fld)
      end if
!
      if(cp_prop%iflag_scheme .gt. id_no_evolution) then
        if(iflag_debug .gt. 0) write(*,*)                               &
     &                'sel_scalar_diff_adv_src_euler composition'
        call sel_scalar_diff_adv_src_euler(SGS_param%iflag_SGS_c_flux,  &
     &      sph_bc_C%kr_in, sph_bc_C%kr_out, ipol%i_c_diffuse,          &
     &      ipol%i_c_advect, ipol%i_SGS_div_c_flux,                     &
     &      ipol%i_light_source, ipol%i_light, dt, cp_prop%coef_exp,    &
     &      cp_prop%coef_advect, cp_prop%coef_source, sph_rj, rj_fld)
      end if
!
      if (i_step .eq. 1) then
        if(fl_prop%iflag_scheme .gt.     id_no_evolution) then
          call set_ini_adams_inertia(ipol, itor,                        &
     &        rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
        end if
!
        if(cd_prop%iflag_Bevo_scheme .gt.    id_no_evolution) then
          call set_ini_adams_mag_induct(ipol, itor,                     &
     &        rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
        end if
!
        if(ht_prop%iflag_scheme .gt.     id_no_evolution) then
          call sel_ini_adams_scalar_w_src(SGS_param%iflag_SGS_h_flux,   &
     &        sph_bc_T%kr_in, sph_bc_T%kr_out, ipol%i_h_advect,         &
     &        ipol%i_SGS_div_h_flux, ipol%i_heat_source,                &
     &        ipol%i_pre_heat, ht_prop%coef_source, sph_rj, rj_fld)
        end if
!
        if(cp_prop%iflag_scheme .gt. id_no_evolution) then
          call sel_ini_adams_scalar_w_src(SGS_param%iflag_SGS_c_flux,   &
     &        sph_bc_C%kr_in, sph_bc_C%kr_out, ipol%i_c_advect,         &
     &        ipol%i_SGS_div_c_flux, ipol%i_light_source,               &
     &        ipol%i_pre_composit, cp_prop%coef_source, sph_rj, rj_fld)
        end if
      end if
!$omp end parallel
!
!   Center evolution
!
      if(ht_prop%iflag_scheme .gt.     id_no_evolution) then
        if(iflag_debug .gt. 0) write(*,*)                               &
     &                'sel_ctr_scl_diff_adv_src_euler temperature'
        call sel_ctr_scl_diff_adv_src_euler                             &
     &     (SGS_param%iflag_SGS_h_flux, ipol%i_t_diffuse,               &
     &      ipol%i_h_advect, ipol%i_SGS_div_h_flux,                     &
     &      ipol%i_heat_source, ipol%i_temp, dt, ht_prop%coef_exp,      &
     &      ht_prop%coef_advect, ht_prop%coef_source, sph_rj, rj_fld)
!
        if (i_step .eq. 1) then
          call sel_ctr_ini_adams_scl_w_src                              &
     &       (SGS_param%iflag_SGS_h_flux, ipol%i_h_advect,              &
     &        ipol%i_SGS_div_h_flux, ipol%i_heat_source,                &
     &        ipol%i_pre_heat, ht_prop%coef_source, sph_rj, rj_fld)
        end if
      end if
!
      if(cp_prop%iflag_scheme .gt. id_no_evolution) then
        if(iflag_debug .gt. 0) write(*,*)                               &
     &                'sel_ctr_scl_diff_adv_src_euler composition'
        call sel_ctr_scl_diff_adv_src_euler                             &
     &     (SGS_param%iflag_SGS_c_flux, ipol%i_c_diffuse,               &
     &      ipol%i_c_advect, ipol%i_SGS_div_c_flux,                     &
     &      ipol%i_light_source, ipol%i_light, dt, cp_prop%coef_exp,    &
     &      cp_prop%coef_advect, cp_prop%coef_source, sph_rj, rj_fld)
!
        if (i_step .eq. 1) then
          call sel_ctr_ini_adams_scl_w_src                              &
     &       (SGS_param%iflag_SGS_c_flux, ipol%i_c_advect,              &
     &        ipol%i_SGS_div_c_flux, ipol%i_light_source,               &
     &        ipol%i_pre_composit, cp_prop%coef_source, sph_rj, rj_fld)
        end if
      end if
!
      end subroutine cal_explicit_sph_euler
!
! ----------------------------------------------------------------------
!
      end module cal_momentum_eq_explicit
