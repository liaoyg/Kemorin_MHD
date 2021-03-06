!cal_vector_potential_pre.f90
!      module cal_vector_potential_pre
!
!        programmed by H.Matsui and H.Okuda
!                                    on July 2000 (ver 1.1)
!        modieied by H. Matsui on Sep., 2005
!
!!      subroutine cal_vector_p_pre                                     &
!!     &         (iak_diff_b, icomp_sgs_uxb, ie_dvx, ak_d_magne, dt,    &
!!     &          FEM_prm, SGS_param, cmt_param, filter_param,          &
!!     &          nod_comm, node, ele, surf, conduct, sf_grp,           &
!!     &          cd_prop, Bnod_bcs, Asf_bcs, iphys, iphys_ele, ele_fld,&
!!     &          jacs, rhs_tbl, FEM_elens, sgs_coefs, diff_coefs,      &
!!     &          filtering, mlump_cd, Bmatrix, MG_vector,              &
!!     &          wk_filter, mhd_fem_wk, fem_wk, f_l, f_nl, nod_fld)
!!      subroutine cal_vector_p_co(iak_diff_b, ak_d_magne, dt,          &
!!     &          FEM_prm, SGS_param, cmt_param,                        &
!!     &          nod_comm, node, ele, surf, conduct, sf_grp, cd_prop,  &
!!     &          Bnod_bcs, Fsf_bcs, iphys, iphys_ele, ele_fld,         &
!!     &          jacs, rhs_tbl, FEM_elens, diff_coefs, m_lump,         &
!!     &          Bmatrix, MG_vector, mhd_fem_wk, fem_wk, surf_wk,      &
!!     &          f_l, f_nl, nod_fld)
!!        type(FEM_MHD_paremeters), intent(in) :: FEM_prm
!!        type(SGS_model_control_params), intent(in) :: SGS_param
!!        type(commutation_control_params), intent(in) :: cmt_param
!!        type(SGS_filtering_params), intent(in) :: filter_param
!!        type(communication_table), intent(in) :: nod_comm
!!        type(node_data), intent(in) :: node
!!        type(element_data), intent(in) :: ele
!!        type(surface_data), intent(in) :: surf
!!        type(surface_group_data), intent(in) :: sf_grp
!!        type(field_geometry_data), intent(in) :: conduct
!!        type(conductive_property), intent(in) :: cd_prop
!!        type(nodal_bcs_4_induction_type), intent(in) :: Bnod_bcs
!!        type(velocity_surf_bc_type), intent(in) :: Asf_bcs
!!        type(potential_surf_bc_type), intent(in) :: Fsf_bcs
!!        type(phys_address), intent(in) :: iphys
!!        type(phys_address), intent(in) :: iphys_ele
!!        type(phys_data), intent(in) :: ele_fld
!!        type(jacobians_type), intent(in) :: jacs
!!        type(tables_4_FEM_assembles), intent(in) :: rhs_tbl
!!        type(lumped_mass_matrices), intent(in) :: m_lump
!!        type(gradient_model_data_type), intent(in) :: FEM_elens
!!        type(SGS_coefficients_type), intent(in) :: sgs_coefs
!!        type(SGS_coefficients_type), intent(in) :: diff_coefs
!!        type(filtering_data_type), intent(in) :: filtering
!!        type(lumped_mass_matrices), intent(in) :: mlump_cd
!!        type(MHD_MG_matrix), intent(in) :: Bmatrix
!!        type(vectors_4_solver), intent(inout)                         &
!!       &           :: MG_vector(0:num_MG_level)
!!        type(filtering_work_type), intent(inout) :: wk_filter
!!        type(work_MHD_fe_mat), intent(inout) :: mhd_fem_wk
!!        type(work_finite_element_mat), intent(inout) :: fem_wk
!!        type(work_surface_element_mat), intent(inout) :: surf_wk
!!        type(finite_ele_mat_node), intent(inout) :: f_l, f_nl
!!        type(phys_data), intent(inout) :: nod_fld
!
      module cal_vector_potential_pre
!
      use m_precision
!
      use m_machine_parameter
      use m_phys_constants
!
      use t_FEM_control_parameter
      use t_SGS_control_parameter
      use t_comm_table
      use t_geometry_data_MHD
      use t_geometry_data
      use t_surface_data
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
      use t_solver_djds
      use t_solver_djds_MHD
      use t_interpolate_table
      use t_vector_for_solver
      use t_material_property
      use t_bc_data_magne
      use t_surface_bc_data
      use t_physical_property
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine cal_vector_p_pre                                       &
     &         (iak_diff_b, icomp_sgs_uxb, ie_dvx, ak_d_magne, dt,      &
     &          FEM_prm, SGS_param, cmt_param, filter_param,            &
     &          nod_comm, node, ele, surf, conduct, sf_grp,             &
     &          cd_prop, Bnod_bcs, Asf_bcs, iphys, iphys_ele, ele_fld,  &
     &          jacs, rhs_tbl, FEM_elens, sgs_coefs, diff_coefs,        &
     &          filtering, mlump_cd, Bmatrix, MG_vector,                &
     &          wk_filter, mhd_fem_wk, fem_wk, f_l, f_nl, nod_fld)
!
      use calypso_mpi
!
      use set_boundary_scalars
      use nod_phys_send_recv
      use cal_sgs_fluxes
      use int_vol_diffusion_ele
      use int_vol_vect_p_pre
      use int_surf_fixed_gradients
      use evolve_by_1st_euler
      use evolve_by_adams_bashforth
      use evolve_by_lumped_crank
      use evolve_by_consist_crank
      use copy_nodal_fields
!
      type(FEM_MHD_paremeters), intent(in) :: FEM_prm
      type(SGS_model_control_params), intent(in) :: SGS_param
      type(commutation_control_params), intent(in) :: cmt_param
      type(SGS_filtering_params), intent(in) :: filter_param
      type(communication_table), intent(in) :: nod_comm
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(surface_data), intent(in) :: surf
      type(surface_group_data), intent(in) :: sf_grp
      type(field_geometry_data), intent(in) :: conduct
      type(conductive_property), intent(in) :: cd_prop
      type(nodal_bcs_4_induction_type), intent(in) :: Bnod_bcs
      type(velocity_surf_bc_type), intent(in) :: Asf_bcs
      type(phys_address), intent(in) :: iphys
      type(phys_address), intent(in) :: iphys_ele
      type(phys_data), intent(in) :: ele_fld
      type(jacobians_type), intent(in) :: jacs
      type(tables_4_FEM_assembles), intent(in) :: rhs_tbl
      type(gradient_model_data_type), intent(in) :: FEM_elens
      type(SGS_coefficients_type), intent(in) :: sgs_coefs
      type(SGS_coefficients_type), intent(in) :: diff_coefs
      type(filtering_data_type), intent(in) :: filtering
      type(lumped_mass_matrices), intent(in) :: mlump_cd
      type(MHD_MG_matrix), intent(in) :: Bmatrix
!
      integer(kind = kint), intent(in) :: iak_diff_b, icomp_sgs_uxb
      integer(kind = kint), intent(in) :: ie_dvx
      real(kind = kreal), intent(in) :: ak_d_magne(ele%numele)
      real(kind = kreal), intent(in) :: dt
!
      type(vectors_4_solver), intent(inout)                             &
     &           :: MG_vector(0:Bmatrix%nlevel_MG)
      type(filtering_work_type), intent(inout) :: wk_filter
      type(work_MHD_fe_mat), intent(inout) :: mhd_fem_wk
      type(work_finite_element_mat), intent(inout) :: fem_wk
      type(finite_ele_mat_node), intent(inout) :: f_l, f_nl
      type(phys_data), intent(inout) :: nod_fld
!
!
      call reset_ff_smps(node%max_nod_smp, f_l, f_nl)
!
!   lead diffusion term
!
      if (cd_prop%coef_magne .gt. zero                                  &
     &     .and. cd_prop%coef_exp .gt. zero) then
        call int_vol_vector_diffuse_ele(SGS_param%ifilter_final,        &
     &      ele%istack_ele_smp, FEM_prm%npoint_t_evo_int,               &
     &      node, ele, nod_fld, jacs%g_FEM, jacs%jac_3d, rhs_tbl,       &
     &      FEM_elens, diff_coefs, iak_diff_b, cd_prop%coef_exp,        &
     &      ak_d_magne, iphys%i_vecp, fem_wk, f_l)
      end if
!
!  lead induction terms
!
      if ( SGS_param%iflag_SGS_uxb .ne. id_SGS_none) then
        call cal_sgs_uxb_2_evo(icomp_sgs_uxb, ie_dvx, dt,               &
     &      FEM_prm, SGS_param, filter_param, nod_comm, node, ele,      &
     &      conduct, cd_prop, iphys, iphys_ele, ele_fld,                &
     &      jacs, rhs_tbl, FEM_elens, filtering, sgs_coefs,             &
     &      wk_filter, mhd_fem_wk, fem_wk, f_nl, nod_fld)
      end if
!
      if (FEM_prm%iflag_magne_supg .gt. id_turn_OFF) then
        call int_vol_vect_p_pre_ele_upm(FEM_prm%npoint_t_evo_int, dt,   &
     &      node, ele, conduct, cd_prop, iphys, nod_fld,                &
     &      ele_fld%ntot_phys, iphys_ele%i_magne, ele_fld%d_fld,        &
     &      jacs%g_FEM, jacs%jac_3d, rhs_tbl, mhd_fem_wk,               &
     &      fem_wk, f_nl)
      else
        call int_vol_vect_p_pre_ele(FEM_prm%npoint_t_evo_int,           &
     &      node, ele, conduct, cd_prop, iphys, nod_fld,                &
     &      ele_fld%ntot_phys, iphys_ele%i_magne, ele_fld%d_fld,        &
     &      jacs%g_FEM, jacs%jac_3d, rhs_tbl, mhd_fem_wk,               &
     &      fem_wk, f_nl)
      end if
!
      call int_sf_grad_velocity                                         &
     &   (node, ele, surf, sf_grp, jacs%g_FEM, jacs%jac_sf_grp,         &
     &    rhs_tbl, Asf_bcs%grad, FEM_prm%npoint_t_evo_int, ak_d_magne,  &
     &    fem_wk, f_l)
!
!      call check_nodal_data                                            &
!     &   ((50+my_rank), nod_fld, n_vector, iphys%i_velo)
!      call check_nodal_data                                            &
!     &   ((50+my_rank), ele_fld, n_vector, iphys_ele%i_magne)
!      call check_ff_smp(my_rank, n_vector, node%max_nod_smp, f_l)
!      call check_ff_smp(my_rank, n_vector, node%max_nod_smp, f_nl)
!
      if (iflag_debug.eq.1) write(*,*) 'coefs_4_time_evolution_end'
!
!  -----for explicit euler
      if (cd_prop%iflag_Aevo_scheme .eq. id_explicit_euler) then
        call cal_magne_pre_euler(iphys%i_vecp, dt,                      &
     &      FEM_prm, nod_comm, node, ele, conduct, iphys_ele, ele_fld,  &
     &      jacs%g_FEM, jacs%jac_3d, rhs_tbl, mlump_cd, mhd_fem_wk,     &
     &      fem_wk,  f_l, f_nl, nod_fld)
!
!  -----for Adams_Bashforth
      else if (cd_prop%iflag_Aevo_scheme .eq. id_explicit_adams2) then
        call cal_magne_pre_adams(iphys%i_vecp, iphys%i_pre_uxb, dt,     &
     &      FEM_prm, nod_comm, node, ele, conduct, iphys_ele, ele_fld,  &
     &      jacs%g_FEM, jacs%jac_3d, rhs_tbl, mlump_cd,                 &
     &      mhd_fem_wk, fem_wk, f_l, f_nl, nod_fld)
!
!  -----for Ceank-nicolson
      else if (cd_prop%iflag_Aevo_scheme .eq. id_Crank_nicolson) then
        call cal_vect_p_pre_lumped_crank                                &
     &     (cmt_param%iflag_c_magne, SGS_param%ifilter_final,           &
     &      iphys%i_vecp, iphys%i_pre_uxb, iak_diff_b, ak_d_magne,      &
     &      Bnod_bcs%nod_bc_a, dt, FEM_prm, nod_comm, node, ele,        &
     &      conduct, cd_prop, iphys_ele, ele_fld, jacs%g_FEM,           &
     &      jacs%jac_3d, rhs_tbl, FEM_elens, diff_coefs, mlump_cd,      &
     &      Bmatrix, MG_vector, mhd_fem_wk, fem_wk, f_l, f_nl, nod_fld)
      else if(cd_prop%iflag_Aevo_scheme .eq. id_Crank_nicolson_cmass)   &
     & then
        call cal_vect_p_pre_consist_crank                               &
     &     (cmt_param%iflag_c_magne, SGS_param%ifilter_final,           &
     &      iphys%i_vecp, iphys%i_pre_uxb, iak_diff_b, ak_d_magne,      &
     &      Bnod_bcs%nod_bc_a, dt, FEM_prm, node, ele, conduct,         &
     &      cd_prop, jacs%g_FEM, jacs%jac_3d, rhs_tbl, FEM_elens,       &
     &      diff_coefs, Bmatrix, MG_vector, mhd_fem_wk, fem_wk,         &
     &      f_l, f_nl, nod_fld)
      end if
!
      call set_boundary_vect(Bnod_bcs%nod_bc_a, iphys%i_vecp, nod_fld)
!
      call vector_send_recv(iphys%i_vecp, nod_comm, nod_fld)
      call clear_field_data(nod_fld, n_scalar, iphys%i_m_phi)
!
!      call check_nodal_data                                            &
!     &   ((50+my_rank), nod_fld, n_vector, iphys%i_vecp)
!
      end subroutine cal_vector_p_pre
!
! ----------------------------------------------------------------------
!
      subroutine cal_vector_p_co(iak_diff_b, ak_d_magne, dt,            &
     &          FEM_prm, SGS_param, cmt_param,                          &
     &          nod_comm, node, ele, surf, conduct, sf_grp, cd_prop,    &
     &          Bnod_bcs, Fsf_bcs, iphys, iphys_ele, ele_fld,           &
     &          jacs, rhs_tbl, FEM_elens, diff_coefs, m_lump,           &
     &          Bmatrix, MG_vector, mhd_fem_wk, fem_wk, surf_wk,        &
     &          f_l, f_nl, nod_fld)
!
      use set_boundary_scalars
      use nod_phys_send_recv
      use int_vol_solenoid_correct
      use int_surf_grad_sgs
      use implicit_vector_correct
      use copy_nodal_fields
      use cal_multi_pass
      use cal_sol_vector_co_crank
!
      type(FEM_MHD_paremeters), intent(in) :: FEM_prm
      type(SGS_model_control_params), intent(in) :: SGS_param
      type(commutation_control_params), intent(in) :: cmt_param
      type(communication_table), intent(in) :: nod_comm
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(surface_data), intent(in) :: surf
      type(surface_group_data), intent(in) :: sf_grp
      type(field_geometry_data), intent(in) :: conduct
      type(conductive_property), intent(in) :: cd_prop
      type(nodal_bcs_4_induction_type), intent(in) :: Bnod_bcs
      type(potential_surf_bc_type), intent(in) :: Fsf_bcs
      type(phys_address), intent(in) :: iphys
      type(phys_address), intent(in) :: iphys_ele
      type(phys_data), intent(in) :: ele_fld
      type(jacobians_type), intent(in) :: jacs
      type(tables_4_FEM_assembles), intent(in) :: rhs_tbl
      type(lumped_mass_matrices), intent(in) :: m_lump
      type(gradient_model_data_type), intent(in) :: FEM_elens
      type(SGS_coefficients_type), intent(in) :: diff_coefs
      type(MHD_MG_matrix), intent(in) :: Bmatrix
!
      integer(kind = kint), intent(in) :: iak_diff_b
      real(kind = kreal), intent(in) :: ak_d_magne(ele%numele)
      real(kind = kreal), intent(in) :: dt
!
      type(vectors_4_solver), intent(inout)                             &
     &           :: MG_vector(0:Bmatrix%nlevel_MG)
      type(work_MHD_fe_mat), intent(inout) :: mhd_fem_wk
      type(work_finite_element_mat), intent(inout) :: fem_wk
      type(work_surface_element_mat), intent(inout) :: surf_wk
      type(finite_ele_mat_node), intent(inout) :: f_l, f_nl
      type(phys_data), intent(inout) :: nod_fld
!
!
      call reset_ff_smps(node%max_nod_smp, f_l, f_nl)
!
      if (iflag_debug.eq.1) write(*,*) 'int_vol_magne_co'
      call int_vol_solenoid_co                                          &
     &   (FEM_prm%npoint_poisson_int, SGS_param%ifilter_final,          &
     &    ele%istack_ele_smp, iphys%i_m_phi, iak_diff_b,                &
     &    node, ele, nod_fld, jacs%g_FEM, jacs%jac_3d, jacs%jac_3d_l,   &
     &    rhs_tbl, FEM_elens, diff_coefs, fem_wk, f_nl)
!
      if (cmt_param%iflag_c_magne .eq. id_SGS_commute_ON                &
     &     .and. Fsf_bcs%sgs%ngrp_sf_dat .gt. 0) then
        if (iflag_debug.eq.1) write(*,*)                                &
                             'int_surf_sgs_velo_co_ele', iphys%i_m_phi
         call int_surf_sgs_velo_co_ele                                  &
     &      (node, ele, surf, sf_grp, nod_fld,                          &
     &       jacs%g_FEM, jacs%jac_sf_grp, jacs%jac_sf_grp_l,            &
     &       rhs_tbl, FEM_elens, FEM_prm%npoint_poisson_int,            &
     &       Fsf_bcs%sgs%ngrp_sf_dat, Fsf_bcs%sgs%id_grp_sf_dat,        &
     &       SGS_param%ifilter_final, diff_coefs%num_field,             &
     &       iak_diff_b, diff_coefs%ak, iphys%i_m_phi,                  &
     &       fem_wk, surf_wk, f_nl)
      end if
!
!
      if (   FEM_prm%iflag_imp_correct .eq. id_Crank_nicolson           &
     &  .or. FEM_prm%iflag_imp_correct .eq. id_Crank_nicolson_cmass)    &
     & then
        call cal_vector_p_co_imp(iphys%i_vecp, iak_diff_b, ak_d_magne,  &
     &      dt, FEM_prm, SGS_param, cmt_param, nod_comm, node, ele,     &
     &      conduct, cd_prop, Bnod_bcs, iphys_ele, ele_fld, jacs%g_FEM, &
     &      jacs%jac_3d, rhs_tbl, FEM_elens, diff_coefs, m_lump,        &
     &      Bmatrix, MG_vector, mhd_fem_wk, fem_wk, f_l, f_nl, nod_fld)
        call clear_field_data(nod_fld, n_scalar, iphys%i_m_phi)
      else
        call cal_vector_p_co_exp                                        &
     &     (iphys%i_vecp, FEM_prm, nod_comm, node, ele, jacs%g_FEM,     &
     &      jacs%jac_3d, rhs_tbl, m_lump, mhd_fem_wk, fem_wk,           &
     &      f_l, f_nl, nod_fld)
      end if
!
      if (iflag_debug.eq.1) write(*,*) 'set_boundary_vect vect_p'
      call set_boundary_vect(Bnod_bcs%nod_bc_a, iphys%i_vecp, nod_fld)
!
      if (iflag_debug.eq.1) write(*,*) 'vector_send_recv for vector_p'
      call vector_send_recv(iphys%i_vecp, nod_comm, nod_fld)
      if (iflag_debug.eq.1) write(*,*) 'scalar_send_recv for potential'
      call scalar_send_recv(iphys%i_mag_p, nod_comm, nod_fld)
!
      end subroutine cal_vector_p_co
!
! ----------------------------------------------------------------------
!
      end module cal_vector_potential_pre
