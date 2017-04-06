!>@file   init_radial_infos_sph_mhd.f90
!!@brief  module init_radial_infos_sph_mhd
!!
!!@author H. Matsui
!!@date Programmed in June., 1994
!!@n    Modified in Jan, 2010
!
!>@brief  Coefficients to obtain radial derivatives
!!        by finite difference method
!!
!!@verbatim
!!      subroutine init_r_infos_sph_mhd_evo(fl_prop, sph_grps, ipol,    &
!!     &          fl_prop, cd_prop, ht_prop, cp_prop,                   &
!!     &          ref_param_T1, ref_param_C1, takepito_T1, takepito_C1, &
!!     &          r_2nd, rj_fld)
!!      subroutine init_r_infos_sph_mhd(sph_grps, ipol,                 &
!!     &          fl_prop, cd_prop, ht_prop, cp_prop,                   &
!!     &          sph, omega_sph, ref_temp, ref_comp, rj_fld,           &
!!     &          ref_param_T1, ref_param_C1, takepito_T1, takepito_C1)
!!        type(fluid_property), intent(in) :: fl_prop
!!        type(sph_group_data), intent(in) :: sph_grps
!!        type(phys_address), intent(in) :: ipol
!!        type(sph_grids), intent(inout) :: sph
!!        type(sph_rotation), intent(inout) :: omega_sph
!!        type(reference_temperature), intent(inout) :: ref_temp
!!        type(fluid_property), intent(inout) :: fl_prop
!!        type(conductive_property), intent(inout) :: cd_prop
!!        type(scalar_property), intent(inout) :: ht_prop, cp_prop
!!        type(fdm_matrices), intent(inout) :: r_2nd
!!        type(phys_data), intent(inout) :: rj_fld
!!        type(reference_scalar_param), intent(inout)  :: ref_param_T1
!!        type(reference_scalar_param), intent(inout)  :: ref_param_C1
!!        type(takepiro_model_param), intent(inout)  :: takepito_T1
!!        type(takepiro_model_param), intent(inout)  :: takepito_C1
!!@endverbatim
!!
!!@n @param r_hot        radius at highest temperature point
!!@n @param r_cold       radius at lowest temperature point
!!@n @param temp_hot     temperature at highest temperature point
!!@n @param temp_cold    temperature at lowest temperature point
!!@n @param rotate(3)    rotation vector
!
      module init_radial_infos_sph_mhd
!
      use m_precision
      use calypso_mpi
!
      use m_constants
      use m_machine_parameter
      use m_spheric_constants
!
      use t_physical_property
      use t_reference_scalar_param
      use t_spheric_parameter
      use t_spheric_mesh
      use t_group_data
      use t_poloidal_rotation
      use t_radial_reference_temp
      use t_phys_address
      use t_phys_data
      use t_fdm_coefs
!
      implicit none
!
      private :: set_radius_rot_reft_dat_4_sph, init_reference_temps
!
!  -------------------------------------------------------------------
!
      contains
!
!  -------------------------------------------------------------------
!
      subroutine init_r_infos_sph_mhd_evo(sph_grps, ipol,               &
     &          sph, omega_sph, ref_temp, ref_comp,                     &
     &          fl_prop, cd_prop, ht_prop, cp_prop,                     &
     &          ref_param_T1, ref_param_C1, takepito_T1, takepito_C1,   &
     &          r_2nd, rj_fld)
!
      use calypso_mpi
      use const_fdm_coefs
      use material_property
!
      type(sph_group_data), intent(in) :: sph_grps
      type(phys_address), intent(in) :: ipol
!
      type(sph_grids), intent(inout) :: sph
      type(sph_rotation), intent(inout) :: omega_sph
      type(reference_temperature), intent(inout) :: ref_temp, ref_comp
      type(fluid_property), intent(inout) :: fl_prop
      type(conductive_property), intent(inout) :: cd_prop
      type(scalar_property), intent(inout) :: ht_prop, cp_prop
      type(reference_scalar_param), intent(inout)  :: ref_param_T1
      type(reference_scalar_param), intent(inout)  :: ref_param_C1
      type(takepiro_model_param), intent(inout)  :: takepito_T1
      type(takepiro_model_param), intent(inout)  :: takepito_C1
      type(fdm_matrices), intent(inout) :: r_2nd
      type(phys_data), intent(inout) :: rj_fld
!
!
      call init_r_infos_sph_mhd(sph_grps, ipol,                         &
     &    fl_prop, cd_prop, ht_prop, cp_prop,                           &
     &    sph, omega_sph, ref_temp, ref_comp, rj_fld,                   &
     &    ref_param_T1, ref_param_C1, takepito_T1, takepito_C1)
!
      if (iflag_debug.gt.0) write(*,*) 'const_2nd_fdm_matrices'
      call const_2nd_fdm_matrices(sph%sph_params, sph%sph_rj, r_2nd)
!
      if(iflag_debug.gt.0) write(*,*)' set_material_property'
      call set_material_property                                        &
     &   (ipol, sph%sph_params%radius_CMB, sph%sph_params%radius_ICB,   &
     &    fl_prop, cd_prop, ht_prop, cp_prop)
!
      end subroutine init_r_infos_sph_mhd_evo
!
!  -------------------------------------------------------------------
!
      subroutine init_r_infos_sph_mhd(sph_grps, ipol,                   &
     &          fl_prop, cd_prop, ht_prop, cp_prop,                     &
     &          sph, omega_sph, ref_temp, ref_comp, rj_fld,             &
     &          ref_param_T1, ref_param_C1, takepito_T1, takepito_C1)
!
      use m_boundary_params_sph_MHD
!
      use set_bc_sph_mhd
!
      type(fluid_property), intent(in) :: fl_prop
      type(conductive_property), intent(in) :: cd_prop
      type(scalar_property), intent(in) :: ht_prop, cp_prop
      type(sph_group_data), intent(in) :: sph_grps
      type(phys_address), intent(in) :: ipol
!
      type(sph_grids), intent(inout) :: sph
      type(sph_rotation), intent(inout) :: omega_sph
      type(reference_temperature), intent(inout) :: ref_temp, ref_comp
      type(reference_scalar_param), intent(inout)  :: ref_param_T1
      type(reference_scalar_param), intent(inout)  :: ref_param_C1
      type(takepiro_model_param), intent(inout)  :: takepito_T1
      type(takepiro_model_param), intent(inout)  :: takepito_C1
      type(phys_data), intent(inout) :: rj_fld
!
!
      if (iflag_debug.gt.0) write(*,*) 'set_radius_rot_reft_dat_4_sph'
      call set_radius_rot_reft_dat_4_sph                                &
     &   (sph_grps%radial_rj_grp, sph%sph_params, sph%sph_rj)
!
!*  ----------  rotation of earth  ---------------
!
      if (iflag_debug .ge. iflag_routine_msg)                           &
     &                write(*,*) 'set_rot_earth_4_sph'
      call set_rot_earth_4_sph(sph%sph_rlm, sph%sph_rj,                 &
     &    fl_prop, omega_sph)
!
!*  ---------- boudary conditions  ---------------
      if(iflag_debug.gt.0) write(*,*) 's_set_bc_sph_mhd'
      call s_set_bc_sph_mhd                                             &
     &   (sph%sph_params, sph%sph_rj, sph_grps%radial_rj_grp,           &
     &    fl_prop, cd_prop, ht_prop, cp_prop,                           &
     &    CTR_nod_grp_name, CTR_sf_grp_name)
!
      call init_reference_temps(ref_param_T1, takepito_T1,              &
     &    sph%sph_params, sph%sph_rj, ipol%i_ref_t, ipol%i_gref_t,      &
     &    ref_temp, rj_fld, sph_bc_T)
      call init_reference_temps(ref_param_C1, takepito_C1,              &
     &    sph%sph_params, sph%sph_rj, ipol%i_ref_c, ipol%i_gref_c,      &
     &    ref_comp, rj_fld, sph_bc_C)
!
      end subroutine init_r_infos_sph_mhd
!
!  -------------------------------------------------------------------
!  -------------------------------------------------------------------
!
      subroutine set_radius_rot_reft_dat_4_sph                          &
     &         (radial_rj_grp, sph_params, sph_rj)
!
      use set_radius_func_noequi
      use set_radius_4_sph_dynamo
      use set_reference_temp_sph
!
      type(group_data), intent(in) :: radial_rj_grp
!
      type(sph_rj_grid), intent(inout) :: sph_rj
      type(sph_shell_parameters), intent(inout) :: sph_params
!
!* --------  radius  --------------
!
      if (iflag_debug .ge. iflag_routine_msg)                           &
     &      write(*,*) 'set_radius_dat_4_sph_dynamo'
      call set_radius_dat_4_sph_dynamo                                  &
     &   (sph_rj%nidx_rj(1), sph_rj%radius_1d_rj_r, radial_rj_grp,      &
     &    sph_params%iflag_radial_grid, sph_params%nlayer_ICB,          &
     &    sph_params%nlayer_CMB, sph_params%nlayer_2_center,            &
     &    sph_rj%ar_1d_rj, sph_rj%r_ele_rj, sph_rj%ar_ele_rj,           &
     &    sph_params%radius_ICB, sph_params%radius_CMB,                 &
     &    sph_params%R_earth)
!
!   Choose radial grid mode
      if (iflag_debug .ge. iflag_routine_msg)                           &
     &      write(*,*) 'set_dr_for_nonequi'
      call allocate_dr_rj_noequi(sph_rj%nidx_rj(1))
      call set_dr_for_nonequi(sph_params%nlayer_CMB,                    &
     &    sph_rj%nidx_rj(1), sph_rj%radius_1d_rj_r)
!*
      end subroutine set_radius_rot_reft_dat_4_sph
!
!  -------------------------------------------------------------------
!
      subroutine init_reference_temps(ref_param, takepiro,              &
     &          sph_params, sph_rj, i_ref, i_gref,                      &
     &          ref_fld, rj_fld, sph_bc_S)
!
      use t_boundary_params_sph_MHD
      use t_reference_scalar_param
      use set_reference_sph_mhd
      use set_reference_temp_sph
!
      integer(kind = kint), intent(in) :: i_ref, i_gref
      type(reference_scalar_param), intent(inout) :: ref_param
      type(takepiro_model_param), intent(in) :: takepiro
      type(sph_shell_parameters), intent(in) :: sph_params
      type(sph_rj_grid), intent(in) ::  sph_rj
!
      type(reference_temperature), intent(inout) :: ref_fld
      type(phys_data), intent(inout) :: rj_fld
      type(sph_boundary_type), intent(inout) :: sph_bc_S
!
!      Set reference temperature and adjust boundary conditions
!
      if(iflag_debug .gt. 0) write(*,*) 'alloc_reft_rj_data'
      call alloc_reft_rj_data(sph_rj%nidx_rj(1), ref_fld)
!
      if (ref_param%iflag_reference .eq. id_sphere_ref_temp) then
        if(iflag_debug .gt. 0) write(*,*) 'set_ref_temp_sph_mhd'
        call set_ref_temp_sph_mhd                                       &
     &   (ref_param%low_value, ref_param%depth_top,                     &
     &    ref_param%high_value, ref_param%depth_bottom,                 &
     &    sph_rj%nidx_rj, sph_rj%radius_1d_rj_r, sph_rj%ar_1d_rj,       &
     &    ref_fld%t_rj)
        call adjust_sph_temp_bc_by_reftemp                              &
     &     (sph_rj%idx_rj_degree_zero, sph_rj%nidx_rj(1),               &
     &      ref_fld%t_rj, sph_bc_S)
!
      else if(ref_param%iflag_reference .eq. id_takepiro_temp) then
        call set_stratified_sph_mhd(takepiro%stratified_sigma,          &
     &    takepiro%stratified_width, takepiro%stratified_outer_r,       &
     &    sph_rj%nidx_rj, sph_params%radius_ICB, sph_params%radius_CMB, &
     &    sph_params%nlayer_ICB, sph_params%nlayer_CMB,                 &
     &    sph_rj%radius_1d_rj_r, ref_fld%t_rj)
!!        call adjust_sph_temp_bc_by_reftemp                            &
!!     &     (sph_rj%idx_rj_degree_zero, sph_rj%nidx_rj(1),             &
!!     &      ref_fld%t_rj, sph_bc_S)
!
      else
        call no_ref_temp_sph_mhd(ref_param%depth_top,                   &
     &      ref_param%depth_bottom, sph_rj%nidx_rj(1),                  &
     &      sph_params%radius_ICB, sph_params%radius_CMB, ref_fld%t_rj)
      end if
!
!      call set_reftemp_4_sph(sph_rj%idx_rj_degree_zero, sph_rj%nidx_rj,&
!     &    ref_fld%t_rj, i_ref, i_gref,                                 &
!     &    rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
!
      end subroutine init_reference_temps
!
! -----------------------------------------------------------------------
!
      end module init_radial_infos_sph_mhd
