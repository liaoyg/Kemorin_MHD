!
!      module cal_velocity
!
!        programmed by H.Matsui and H.Okuda
!                                    on July 2000 (ver 1.1)
!        modified by H.Matsui on July, 2006
!
!!      subroutine velocity_evolution(nod_comm, node, ele, surf, fluid, &
!!     &          sf_grp, sf_grp_nod, Vnod_bcs, Vsf_bcs, Bsf_bcs,       &
!!     &          Psf_bcs, iphys, iphys_ele, ele_fld,                   &
!!     &          jac_3d_q, jac_3d_l, jac_sf_grp_q, jac_sf_grp_l,       &
!!     &          rhs_tbl, FEM_elens, filtering, layer_tbl,             &
!!     &          mhd_fem_wk, fem_wk, f_l, f_nl, nod_fld)
!!        type(communication_table), intent(in) :: nod_comm
!!        type(node_data), intent(in) :: node
!!        type(element_data), intent(in) :: ele
!!        type(surface_data), intent(in) :: surf
!!        type(surface_group_data), intent(in) :: sf_grp
!!        type(surface_node_grp_data), intent(in) :: sf_grp_nod
!!        type(field_geometry_data), intent(in) :: fluid
!!        type(nodal_bcs_4_momentum_type), intent(in) :: Vnod_bcs
!!        type(velocity_surf_bc_type), intent(in) :: Vsf_bcs
!!        type(vector_surf_bc_type), intent(in) :: Bsf_bcs
!!        type(phys_address), intent(in) :: iphys
!!        type(phys_address), intent(in) :: iphys_ele
!!        type(phys_data), intent(in) :: ele_fld
!!        type(jacobians_3d), intent(in) :: jac_3d_q, jac_3d_l
!!        type(jacobians_2d), intent(in) :: jac_sf_grp_q, jac_sf_grp_l
!!        type(tables_4_FEM_assembles), intent(in) :: rhs_tbl
!!        type(gradient_model_data_type), intent(in) :: FEM_elens
!!        type(filtering_data_type), intent(in) :: filtering
!!        type(layering_tbl), intent(in) :: layer_tbl
!!        type(work_MHD_fe_mat), intent(inout) :: mhd_fem_wk
!!        type(work_finite_element_mat), intent(inout) :: fem_wk
!!        type(finite_ele_mat_node), intent(inout) :: f_l, f_nl
!!        type(phys_data), intent(inout) :: nod_fld
!
      module cal_velocity
!
      use m_precision
!
      use t_comm_table
      use t_geometry_data_MHD
      use t_geometry_data
      use t_surface_data
      use t_group_data
      use t_surface_group_connect
      use t_phys_data
      use t_phys_address
      use t_jacobian_3d
      use t_jacobian_2d
      use t_table_FEM_const
      use t_finite_element_mat
      use t_MHD_finite_element_mat
      use t_filter_elength
      use t_filtering_data
      use t_layering_ele_list
      use t_bc_data_velo
      use t_surface_bc_data
!
      implicit none
!
      real(kind = kreal) :: ave_pr0, rms_pr0
      private :: ave_pr0, rms_pr0
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine velocity_evolution(nod_comm, node, ele, surf, fluid,   &
     &          sf_grp, sf_grp_nod, Vnod_bcs, Vsf_bcs, Bsf_bcs,         &
     &          Psf_bcs, iphys, iphys_ele, ele_fld,                     &
     &          jac_3d_q, jac_3d_l, jac_sf_grp_q, jac_sf_grp_l,         &
     &          rhs_tbl, FEM_elens, filtering, layer_tbl,               &
     &          mhd_fem_wk, fem_wk, f_l, f_nl, nod_fld)
!
      use m_control_parameter
      use m_machine_parameter
      use m_physical_property
      use m_type_AMG_data
      use m_SGS_model_coefs
      use m_SGS_address
      use m_solver_djds_MHD
!
      use cal_velocity_pre
      use cal_mod_vel_potential
      use cal_sol_pressure_MHD
      use init_4_sol_potentials
      use int_rms_div_MHD
      use int_norm_div_MHD
      use cal_rms_potentials
!
      type(communication_table), intent(in) :: nod_comm
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(surface_data), intent(in) :: surf
      type(surface_group_data), intent(in) :: sf_grp
      type(surface_node_grp_data), intent(in) :: sf_grp_nod
      type(field_geometry_data), intent(in) :: fluid
      type(nodal_bcs_4_momentum_type), intent(in) :: Vnod_bcs
      type(velocity_surf_bc_type), intent(in) :: Vsf_bcs
      type(vector_surf_bc_type), intent(in) :: Bsf_bcs
      type(potential_surf_bc_type), intent(in) :: Psf_bcs
      type(phys_address), intent(in) :: iphys
      type(phys_address), intent(in) :: iphys_ele
      type(phys_data), intent(in) :: ele_fld
      type(jacobians_3d), intent(in) :: jac_3d_q, jac_3d_l
      type(jacobians_2d), intent(in) :: jac_sf_grp_q, jac_sf_grp_l
      type(tables_4_FEM_assembles), intent(in) :: rhs_tbl
      type(gradient_model_data_type), intent(in) :: FEM_elens
      type(filtering_data_type), intent(in) :: filtering
      type(layering_tbl), intent(in) :: layer_tbl
!
      type(work_MHD_fe_mat), intent(inout) :: mhd_fem_wk
      type(work_finite_element_mat), intent(inout) :: fem_wk
      type(finite_ele_mat_node), intent(inout) :: f_l, f_nl
      type(phys_data), intent(inout) :: nod_fld
!
      integer(kind=kint) :: iloop
      real(kind = kreal) :: rel_correct
!
!
      if (iflag_4_lorentz .eq. id_turn_ON) then
        if (iflag_4_rotate .eq. id_turn_OFF) then
          call cal_sol_pressure_w_mag_ene                               &
     &       (node%numnod, node%istack_internal_smp,                    &
     &        nod_fld%ntot_phys, iphys%i_p_phi, iphys%i_magne,          &
     &        iphys%i_press, nod_fld%d_fld)
        else if (iflag_magneto_cv .eq. id_turn_ON                       &
     &     .and. iflag_4_rotate .eq. id_turn_OFF) then
          call cal_sol_pressure_mcv                                     &
     &       (node%numnod, node%istack_internal_smp,                    &
     &        nod_fld%ntot_phys, iphys%i_p_phi, iphys%i_magne,          &
     &        iphys%i_press, nod_fld%d_fld)
        else
          call init_sol_potential(node%numnod, node%istack_nod_smp,     &
     &        coef_press, nod_fld%ntot_phys, iphys%i_p_phi,             &
     &        iphys%i_press, nod_fld%d_fld)
        end if
      else
        call init_sol_potential(node%numnod, node%istack_nod_smp,       &
     &      coef_press, nod_fld%ntot_phys, iphys%i_p_phi,               &
     &      iphys%i_press, nod_fld%d_fld)
      end if
!
!     --------------------- 
!
      if (iflag_debug.eq.1)  write(*,*) 's_cal_velocity_pre'
      call s_cal_velocity_pre                                           &
     &   (nod_comm, node, ele, surf, fluid, sf_grp, sf_grp_nod,         &
     &    Vnod_bcs, Vsf_bcs, Bsf_bcs, iphys, iphys_ele, ele_fld,        &
     &    jac_3d_q, jac_3d_l, jac_sf_grp_q, rhs_tbl, FEM_elens,         &
     &    filtering, layer_tbl, num_MG_level,                           &
          MHD1_matrices%MG_interpolate, MHD1_matrices%MG_comm_fluid,    &
     &    MHD1_matrices%MG_DJDS_fluid, MHD1_matrices%Vmat_MG_DJDS,      &
     &     MG_vector, mhd_fem_wk, fem_wk, f_l, f_nl, nod_fld)
!
!     --------------------- 
!
      iloop = -1
      call int_norm_div_v_monitor(iloop, node, ele, fluid,              &
     &    iphys, nod_fld, jac_3d_q, fem_wk, rel_correct)
!      call int_rms_div_v_monitor(iloop, node, ele, fluid,              &
!     &    iphys, nod_fld, jac_3d_q, fem_wk, rel_correct)
!
      do iloop = 0, maxiter
        call cal_mod_potential(iak_diff_v,                              &
     &      node, ele, surf, fluid, sf_grp, Vnod_bcs, Vsf_bcs, Psf_bcs, &
     &      iphys, jac_3d_q, jac_3d_l, jac_sf_grp_l, rhs_tbl,           &
     &      FEM_elens, num_MG_level, MHD1_matrices%MG_interpolate,      &
     &      MHD1_matrices%MG_comm_fluid, MHD1_matrices%MG_DJDS_lin_fl,  &
     &      MHD1_matrices%Pmat_MG_DJDS, MG_vector,                      &
     &      fem_wk, f_l, f_nl, nod_fld)
!
        call cal_sol_pressure                                           &
     &     (node%numnod, node%istack_internal_smp, nod_fld%ntot_phys,   &
     &      iphys%i_p_phi, iphys%i_press,  nod_fld%d_fld)
!
        call cal_velocity_co                                            &
     &     (nod_comm, node, ele, surf, fluid,  sf_grp, sf_grp_nod,      &
     &      Vnod_bcs, Vsf_bcs, Psf_bcs, iphys, iphys_ele, ele_fld,      &
     &      jac_3d_q, jac_3d_l, jac_sf_grp_q, jac_sf_grp_l, rhs_tbl,    &
     &      FEM_elens, num_MG_level, MHD1_matrices%MG_interpolate,      &
     &      MHD1_matrices%MG_comm_fluid, MHD1_matrices%MG_DJDS_fluid,   &
     &      MHD1_matrices%Vmat_MG_DJDS, MG_vector,                      &
     &      mhd_fem_wk, fem_wk, f_l, f_nl, nod_fld)
!
!
        call cal_rms_scalar_potential(iloop, fluid%istack_ele_fld_smp,  &
     &      iphys%i_press, i_rms%i_press, j_ave%i_press,                &
     &      node, ele, nod_fld, jac_3d_q, jac_3d_l, fem_wk,             &
     &      rel_correct, ave_pr0, rms_pr0)
!
!
        if (iflag_debug.eq.1)                                           &
     &         write(12,*) 'average and RMS of presssur correction: ',  &
     &         iloop, ave_pr0, rms_pr0
!
!
        call int_norm_div_v_monitor(iloop, node, ele, fluid,            &
     &      iphys, nod_fld, jac_3d_q, fem_wk, rel_correct)
!        call int_rms_div_v_monitor(iloop, node, ele, fluid,            &
!     &      iphys, nod_fld, jac_3d_q, fem_wk, rel_correct)
!
        if ( abs(rel_correct) .lt. eps_4_velo ) go to 10
!
      end do
 10   continue
!
      if (iflag_4_rotate .eq. id_turn_ON) then
        call cal_sol_pressure_rotate                                    &
     &     (node%numnod, node%istack_internal_smp, nod_fld%ntot_phys,   &
     &      iphys%i_velo, iphys%i_press, nod_fld%d_fld)
      end if
!
      end subroutine velocity_evolution
!
!-----------------------------------------------------------------------
!
      end module cal_velocity
