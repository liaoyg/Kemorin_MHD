!>@file   lead_physical_values.f90
!!        module lead_physical_values
!!
!! @author H. Matsui
!! @date ...when???
!!
!> @brief Evaluate many kind of field data
!!
!!@verbatim
!!      subroutine lead_fields_by_FEM(time_d, FEM_prm, SGS_par,         &
!!     &          mesh, group, ele_mesh, MHD_mesh, fl_prop, cd_prop,    &
!!     &          ht_prop, cp_prop, nod_bcs, surf_bcs,                  &
!!     &          iphys, iphys_ele, ak_MHD, fem_int, FEM_elens,         &
!!     &          icomp_sgs, icomp_diff, ifld_diff, iphys_elediff,      &
!!     &          sgs_coefs, sgs_coefs_nod, filtering, wide_filtering,  &
!!     &          layer_tbl, wk_cor, wk_lsq, wk_diff, wk_filter,        &
!!     &          mhd_fem_wk, rhs_mat, nod_fld, ele_fld, diff_coefs)
!!        type(FEM_MHD_paremeters), intent(in) :: FEM_prm
!!        type(SGS_paremeters), intent(in) :: SGS_par
!!        type(time_data), intent(in) :: time_d
!!        type(mesh_geometry), intent(in) :: mesh
!!        type(mesh_groups), intent(in) ::   group
!!        type(element_geometry), intent(in) :: ele_mesh
!!        type(mesh_data_MHD), intent(in) :: MHD_mesh
!!        type(fluid_property), intent(in) :: fl_prop
!!        type(conductive_property), intent(in) :: cd_prop
!!        type(scalar_property), intent(in) :: ht_prop, cp_prop
!!        type(nodal_boundarty_conditions), intent(in) :: nod_bcs
!!        type(surface_boundarty_conditions), intent(in) :: surf_bcs
!!        type(phys_address), intent(in) :: iphys
!!        type(phys_address), intent(in) :: iphys_ele
!!        type(coefs_4_MHD_type), intent(in) :: ak_MHD
!!        type(finite_element_integration), intent(in) :: fem_int
!!        type(gradient_model_data_type), intent(in) :: FEM_elens
!!        type(SGS_terms_address), intent(in) :: icomp_sgs
!!        type(SGS_terms_address), intent(in) :: ifld_diff
!!        type(SGS_terms_address), intent(in) :: icomp_diff
!!        type(SGS_terms_address), intent(in) :: iphys_elediff
!!        type(SGS_coefficients_type), intent(in) :: sgs_coefs
!!        type(SGS_coefficients_type), intent(in) :: sgs_coefs_nod
!!        type(filtering_data_type), intent(in) :: filtering
!!        type(filtering_data_type), intent(in) :: wide_filtering
!!        type(layering_tbl), intent(in) :: layer_tbl
!!        type(filtering_work_type), intent(inout) :: wk_filter
!!        type(dynamis_correlation_data), intent(inout) :: wk_cor
!!        type(dynamis_least_suare_data), intent(inout) :: wk_lsq
!!        type(dynamic_model_data), intent(inout) :: wk_diff
!!        type(work_MHD_fe_mat), intent(inout) :: mhd_fem_wk
!!        type(arrays_finite_element_mat), intent(inout) :: rhs_mat
!!        type(phys_data), intent(inout) :: nod_fld
!!        type(phys_data), intent(inout) :: ele_fld
!!        type(SGS_coefficients_type), intent(inout) :: diff_coefs
!!@endverbatim
!
      module lead_physical_values
!
      use m_precision
      use m_machine_parameter
!
      use t_FEM_control_parameter
      use t_SGS_control_parameter
      use t_physical_property
      use t_reference_scalar_param
      use t_time_data
      use t_mesh_data
      use t_comm_table
      use t_geometry_data_MHD
      use t_geometry_data
      use t_surface_data
      use t_edge_data
      use t_group_data
      use t_phys_data
      use t_phys_address
      use t_jacobians
      use t_table_FEM_const
      use t_finite_element_mat
      use t_int_surface_data
      use t_MHD_finite_element_mat
      use t_filter_elength
      use t_filtering_data
      use t_layering_ele_list
      use t_bc_data_MHD
      use t_MHD_boundary_data
      use t_material_property
      use t_SGS_model_coefs
      use t_ele_info_4_dynamic
      use t_work_4_dynamic_model
      use t_work_layer_correlate
      use t_work_FEM_integration
!
      implicit none
!
      private :: cal_energy_fluxes
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine lead_fields_by_FEM(time_d, FEM_prm, SGS_par,           &
     &          mesh, group, ele_mesh, MHD_mesh, fl_prop, cd_prop,      &
     &          ht_prop, cp_prop, nod_bcs, surf_bcs,                    &
     &          iphys, iphys_ele, ak_MHD, fem_int, FEM_elens,           &
     &          icomp_sgs, icomp_diff, ifld_diff, iphys_elediff,        &
     &          sgs_coefs, sgs_coefs_nod, filtering, wide_filtering,    &
     &          layer_tbl, wk_cor, wk_lsq, wk_diff, wk_filter,          &
     &          mhd_fem_wk, rhs_mat, nod_fld, ele_fld, diff_coefs)
!
      use update_after_evolution
      use itp_potential_on_edge
      use MHD_field_by_rotation
      use cal_helicities
!
      type(time_data), intent(in) :: time_d
      type(FEM_MHD_paremeters), intent(in) :: FEM_prm
      type(SGS_paremeters), intent(in) :: SGS_par
      type(mesh_geometry), intent(in) :: mesh
      type(mesh_groups), intent(in) ::   group
      type(element_geometry), intent(in) :: ele_mesh
      type(mesh_data_MHD), intent(in) :: MHD_mesh
      type(fluid_property), intent(in) :: fl_prop
      type(conductive_property), intent(in) :: cd_prop
      type(scalar_property), intent(in) :: ht_prop, cp_prop
      type(nodal_boundarty_conditions), intent(in) :: nod_bcs
      type(surface_boundarty_conditions), intent(in) :: surf_bcs
      type(phys_address), intent(in) :: iphys
      type(phys_address), intent(in) :: iphys_ele
      type(coefs_4_MHD_type), intent(in) :: ak_MHD
      type(finite_element_integration), intent(in) :: fem_int
      type(gradient_model_data_type), intent(in) :: FEM_elens
      type(SGS_terms_address), intent(in) :: icomp_sgs
      type(SGS_terms_address), intent(in) :: ifld_diff
      type(SGS_terms_address), intent(in) :: icomp_diff
      type(SGS_terms_address), intent(in) :: iphys_elediff
      type(SGS_coefficients_type), intent(in) :: sgs_coefs
      type(SGS_coefficients_type), intent(in) :: sgs_coefs_nod
      type(filtering_data_type), intent(in) :: filtering
      type(filtering_data_type), intent(in) :: wide_filtering
      type(layering_tbl), intent(in) :: layer_tbl
!
      type(dynamis_correlation_data), intent(inout) :: wk_cor
      type(dynamis_least_suare_data), intent(inout) :: wk_lsq
      type(dynamic_model_data), intent(inout) :: wk_diff
      type(filtering_work_type), intent(inout) :: wk_filter
      type(work_MHD_fe_mat), intent(inout) :: mhd_fem_wk
      type(arrays_finite_element_mat), intent(inout) :: rhs_mat
      type(phys_data), intent(inout) :: nod_fld
      type(phys_data), intent(inout) :: ele_fld
      type(SGS_coefficients_type), intent(inout) :: diff_coefs
!
!
      if (iflag_debug.gt.0) write(*,*) 'cal_potential_on_edge'
      call cal_potential_on_edge                                        &
     &   (mesh%node, mesh%ele, ele_mesh%edge, iphys, nod_fld)
!
      if (iflag_debug.gt.0) write(*,*) 'update_fields'
      call update_fields(time_d, FEM_prm, SGS_par, mesh, group,         &
     &    ele_mesh, MHD_mesh, nod_bcs, surf_bcs, iphys, iphys_ele,      &
     &    fem_int, FEM_elens, ifld_diff, icomp_diff, iphys_elediff,     &
     &    filtering, wide_filtering, layer_tbl, wk_cor, wk_lsq,         &
     &    wk_diff, wk_filter, mhd_fem_wk, rhs_mat,                      &
     &    nod_fld, ele_fld, diff_coefs)
!
      call cal_field_by_rotation                                        &
     &   (time_d%dt, FEM_prm, SGS_par%model_p, SGS_par%commute_p,       &
     &    mesh%nod_comm, mesh%node, mesh%ele, ele_mesh%surf,            &
     &    MHD_mesh%fluid, MHD_mesh%conduct, group%surf_grp,             &
     &    cd_prop, nod_bcs, surf_bcs, iphys, iphys_ele, ele_fld,        &
     &    fem_int%jacobians%jac_3d, fem_int%jacobians%jac_sf_grp,       &
     &    fem_int%rhs_tbl, FEM_elens, ifld_diff, diff_coefs,            &
     &    fem_int%m_lump, mhd_fem_wk, rhs_mat%fem_wk, rhs_mat%surf_wk,  &
     &    rhs_mat%f_l, rhs_mat%f_nl, nod_fld)
!
      if (iflag_debug.gt.0) write(*,*) 'cal_helicity'
      call cal_helicity(iphys, nod_fld)
!
      if (iflag_debug.gt.0) write(*,*) 'cal_energy_fluxes'
      call cal_energy_fluxes                                            &
     &   (time_d%dt, FEM_prm, SGS_par, mesh, group, ele_mesh, MHD_mesh, &
     &    fl_prop, cd_prop, ht_prop, cp_prop,                           &
     &    nod_bcs, surf_bcs, iphys, iphys_ele, ak_MHD, fem_int,         &
     &    FEM_elens, icomp_sgs, ifld_diff, iphys_elediff,               &
     &    sgs_coefs, sgs_coefs_nod, diff_coefs, filtering,              &
     &    wk_filter, mhd_fem_wk, rhs_mat, nod_fld, ele_fld)
!
      end subroutine lead_fields_by_FEM
!
! ----------------------------------------------------------------------
!
      subroutine cal_energy_fluxes                                      &
     &        (dt, FEM_prm, SGS_par, mesh, group, ele_mesh, MHD_mesh,   &
     &         fl_prop, cd_prop, ht_prop, cp_prop,                      &
     &         nod_bcs, surf_bcs, iphys, iphys_ele, ak_MHD, fem_int,    &
     &         FEM_elens, icomp_sgs, ifld_diff, iphys_elediff,          &
     &         sgs_coefs, sgs_coefs_nod, diff_coefs, filtering,         &
     &         wk_filter, mhd_fem_wk, rhs_mat, nod_fld, ele_fld)
!
      use cal_MHD_forces_4_monitor
      use cal_sgs_4_monitor
      use cal_true_sgs_terms
!
      real(kind = kreal), intent(in) :: dt
      type(FEM_MHD_paremeters), intent(in) :: FEM_prm
      type(SGS_paremeters), intent(in) :: SGS_par
      type(mesh_geometry), intent(in) :: mesh
      type(mesh_groups), intent(in) ::   group
      type(element_geometry), intent(in) :: ele_mesh
      type(mesh_data_MHD), intent(in) :: MHD_mesh
      type(fluid_property), intent(in) :: fl_prop
      type(conductive_property), intent(in) :: cd_prop
      type(scalar_property), intent(in) :: ht_prop, cp_prop
      type(nodal_boundarty_conditions), intent(in) :: nod_bcs
      type(surface_boundarty_conditions), intent(in) :: surf_bcs
      type(phys_address), intent(in) :: iphys
      type(phys_address), intent(in) :: iphys_ele
      type(coefs_4_MHD_type), intent(in) :: ak_MHD
      type(finite_element_integration), intent(in) :: fem_int
      type(gradient_model_data_type), intent(in) :: FEM_elens
      type(SGS_terms_address), intent(in) :: icomp_sgs
      type(SGS_terms_address), intent(in) :: ifld_diff
      type(SGS_terms_address), intent(in) :: iphys_elediff
      type(SGS_coefficients_type), intent(in) :: sgs_coefs
      type(SGS_coefficients_type), intent(in) :: sgs_coefs_nod
      type(SGS_coefficients_type), intent(in) :: diff_coefs
      type(filtering_data_type), intent(in) :: filtering
!
      type(filtering_work_type), intent(inout) :: wk_filter
      type(work_MHD_fe_mat), intent(inout) :: mhd_fem_wk
      type(arrays_finite_element_mat), intent(inout) :: rhs_mat
      type(phys_data), intent(inout) :: nod_fld
      type(phys_data), intent(inout) :: ele_fld
!
!
      call cal_true_sgs_terms_pre(dt, FEM_prm, SGS_par,                 &
     &    mesh%nod_comm, mesh%node, mesh%ele, ele_mesh%surf,            &
     &    group%surf_grp, MHD_mesh%fluid, MHD_mesh%conduct,             &
     &    fl_prop, cd_prop, ht_prop, cp_prop,                           &
     &    nod_bcs, surf_bcs, iphys, iphys_ele, ak_MHD,                  &
     &    fem_int%jacobians, fem_int%rhs_tbl, FEM_elens,                &
     &    ifld_diff, diff_coefs, mhd_fem_wk, rhs_mat%fem_wk,            &
     &    rhs_mat%surf_wk, rhs_mat%f_l, rhs_mat%f_nl, nod_fld, ele_fld)
!
      call cal_sgs_terms_4_monitor                                      &
     &   (dt, FEM_prm, SGS_par%model_p, SGS_par%filter_p,               &
     &    mesh%nod_comm, mesh%node, mesh%ele,                           &
     &    MHD_mesh%fluid, MHD_mesh%conduct, cd_prop, iphys,             &
     &    iphys_ele, ele_fld, fem_int%jacobians, fem_int%rhs_tbl,       &
     &    FEM_elens, icomp_sgs, iphys_elediff,                          &
     &    sgs_coefs, sgs_coefs_nod, filtering, wk_filter, mhd_fem_wk,   &
     &    rhs_mat%fem_wk, rhs_mat%f_l, rhs_mat%f_nl, nod_fld)
!
      call cal_fluxes_4_monitor                                         &
     &   (mesh%node, fl_prop, cd_prop, iphys, nod_fld)
!
      call cal_forces_4_monitor(dt, FEM_prm, SGS_par,                   &
     &    mesh%nod_comm, mesh%node, mesh%ele, ele_mesh%surf,            &
     &    MHD_mesh%fluid, MHD_mesh%conduct, group%surf_grp,             &
     &    fl_prop, cd_prop, ht_prop, cp_prop, nod_bcs, surf_bcs,        &
     &    iphys, iphys_ele, ak_MHD, fem_int%jacobians, fem_int%rhs_tbl, &
     &    FEM_elens, ifld_diff, diff_coefs, fem_int%m_lump, mhd_fem_wk, &
     &    rhs_mat%fem_wk, rhs_mat%surf_wk, rhs_mat%f_l, rhs_mat%f_nl,   &
     &    nod_fld, ele_fld)
      call cal_diff_of_sgs_terms                                        &
     &   (dt, FEM_prm, SGS_par%model_p, SGS_par%commute_p,              &
     &    mesh%nod_comm, mesh%node, mesh%ele, ele_mesh%surf,            &
     &    group%surf_grp, MHD_mesh%fluid, MHD_mesh%conduct,             &
     &    fl_prop, cd_prop, ht_prop, cp_prop, nod_bcs, surf_bcs,        &
     &    iphys, iphys_ele, ak_MHD, fem_int%jacobians, fem_int%rhs_tbl, &
     &    FEM_elens, ifld_diff, diff_coefs, mhd_fem_wk,                 &
     &    rhs_mat%fem_wk, rhs_mat%surf_wk, rhs_mat%f_l, rhs_mat%f_nl,   &
     &    nod_fld, ele_fld)
!
      call cal_true_sgs_terms_post                                      &
     &   (SGS_par%filter_p, mesh%nod_comm, mesh%node, iphys,            &
     &    filtering, wk_filter, nod_fld)
!
      call cal_work_4_forces                                            &
     &   (FEM_prm, mesh%nod_comm, mesh%node, mesh%ele, fl_prop,         &
     &    cd_prop, iphys, fem_int%jacobians, fem_int%rhs_tbl,           &
     &    mhd_fem_wk, rhs_mat%fem_wk, rhs_mat%f_nl, nod_fld)
!
      call cal_work_4_sgs_terms(FEM_prm,                                &
     &   mesh%nod_comm, mesh%node, mesh%ele, MHD_mesh%conduct,          &
     &   fl_prop, cd_prop, iphys, fem_int%jacobians, fem_int%rhs_tbl,   &
     &   mhd_fem_wk, rhs_mat%fem_wk, rhs_mat%f_nl, nod_fld)
! 
      end subroutine cal_energy_fluxes
!
!  ---------------------------------------------------------------------
!
      end module lead_physical_values
