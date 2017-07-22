!initialize_4_snapshot.f90
!     module initialize_4_snapshot
!
!      Written by H. Matsui
!
!!      subroutine init_analyzer_snap(fst_file_IO, FEM_prm, SGS_par,    &
!!     &          IO_bc, MHD_step, mesh, group, ele_mesh, MHD_mesh,     &
!!     &          layer_tbl, MHD_prop, ak_MHD, Csims_FEM_MHD,           &
!!     &          iphys, nod_fld, t_IO, rst_step, label_sim)
!!        type(field_IO_params), intent(in) :: fst_file_IO
!!        type(FEM_MHD_paremeters), intent(inout) :: FEM_prm
!!        type(SGS_paremeters), intent(in) :: SGS_par
!!        type(IO_boundary), intent(in) :: IO_bc
!!        type(MHD_step_param), intent(inout) :: MHD_step
!!        type(mesh_geometry), intent(inout) :: mesh
!!        type(mesh_groups), intent(inout) ::   group
!!        type(element_geometry), intent(inout) :: ele_mesh
!!        type(mesh_data_MHD), intent(inout) :: MHD_mesh
!!        type(layering_tbl), intent(inout) :: layer_tbl
!!        type(MHD_evolution_param), intent(inout) :: MHD_prop
!!        type(SGS_coefficients_data), intent(inout) :: Csims_FEM_MHD
!!        type(phys_address), intent(inout) :: iphys
!!        type(phys_data), intent(inout) :: nod_fld
!!        type(time_data), intent(inout) :: t_IO
!!        type(IO_step_param), intent(inout) :: rst_step
!
      module initialize_4_snapshot
!
      use m_precision
      use m_machine_parameter
      use calypso_mpi
!
      use t_control_parameter
      use t_FEM_control_parameter
      use t_SGS_control_parameter
      use t_MHD_step_parameter
      use t_time_data
      use t_phys_data
      use t_phys_address
!
      use t_mesh_data
      use t_geometry_data_MHD
      use t_layering_ele_list
      use t_work_layer_correlate
      use t_time_data
      use t_boundary_field_IO
      use t_IO_step_parameter
      use t_file_IO_parameter
      use t_FEM_SGS_model_coefs
      use t_material_property
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine init_analyzer_snap(fst_file_IO, FEM_prm, SGS_par,      &
     &          IO_bc, MHD_step, mesh, group, ele_mesh, MHD_mesh,       &
     &          layer_tbl, MHD_prop, ak_MHD, Csims_FEM_MHD,             &
     &          iphys, nod_fld, t_IO, rst_step, label_sim)
!
      use m_fem_mhd_restart
!
      use m_work_4_dynamic_model
      use m_boundary_condition_IDs
      use m_array_for_send_recv
      use m_finite_element_matrix
      use m_bc_data_velo
      use m_3d_filter_coef_MHD
      use m_bc_data_list
!
      use count_whole_num_element
!
      use cal_volume_node_MHD
      use int_MHD_mass_matrices
      use int_surface_params_MHD
      use set_nodal_bc_id_data
      use allocate_array_MHD
      use ordering_line_filter_smp
      use const_ele_layering_table
      use const_comm_table_fluid
      use const_bc_infty_surf_type
      use set_reference_value
      use material_property
      use reordering_by_layers
      use set_layers_4_MHD
      use set_istart_3d_filtering
      use count_sgs_components
      use init_sgs_diff_coefs
      use set_layer_list_by_table
      use set_normal_vectors
      use set_table_type_RHS_assemble
      use const_jacobians_sf_grp
      use const_element_comm_tables
      use const_mesh_information
      use init_ele_material_property
!
      use nod_phys_send_recv
!
      type(field_IO_params), intent(in) :: fst_file_IO
      type(SGS_paremeters), intent(in) :: SGS_par
      type(IO_boundary), intent(in) :: IO_bc
      type(MHD_step_param), intent(in) :: MHD_step
!
      type(FEM_MHD_paremeters), intent(inout) :: FEM_prm
      type(mesh_geometry), intent(inout) :: mesh
      type(mesh_groups), intent(inout) ::   group
      type(element_geometry), intent(inout) :: ele_mesh
      type(mesh_data_MHD), intent(inout) :: MHD_mesh
      type(layering_tbl), intent(inout) :: layer_tbl
      type(MHD_evolution_param), intent(inout) :: MHD_prop
      type(coefs_4_MHD_type), intent(inout) :: ak_MHD
      type(SGS_coefficients_data), intent(inout) :: Csims_FEM_MHD
      type(phys_address), intent(inout) :: iphys
      type(phys_data), intent(inout) :: nod_fld
      type(time_data), intent(inout) :: t_IO
      type(IO_step_param), intent(inout) :: rst_step
      character(len=kchara), intent(inout)   :: label_sim
!
      integer(kind = kint) :: iflag
!
!     --------------------- 
!
      if (iflag_debug.eq.1) write(*,*)' reordering_by_layers_snap'
      call reordering_by_layers_snap                                    &
     &   (FEM_prm, SGS_par, mesh%ele, group, MHD_mesh)
!
      if (iflag_debug.eq.1) write(*,*)' set_layers'
      call set_layers                                                   &
     &   (FEM_prm, mesh%node, mesh%ele, group%ele_grp, MHD_mesh)
!
      if (SGS_par%model_p%iflag_dynamic  .ne. id_SGS_DYNAMIC_OFF) then
        if (iflag_debug.eq.1) write(*,*)' const_layers_4_dynamic'
        call const_layers_4_dynamic(group%ele_grp, layer_tbl)
        call alloc_work_4_dynamic(layer_tbl%e_grp%num_grp, wk_lsq1)
        call alloc_work_layer_correlate                                 &
     &     (layer_tbl%e_grp%num_grp, inine, wk_cor1)
      end if
!
!
!     ---------------------
!
      if (iflag_debug.ge.1 ) write(*,*) 'allocate_vector_for_solver'
      call allocate_vector_for_solver(n_sym_tensor, mesh%node%numnod)
!
      call init_send_recv(mesh%nod_comm)
!
      if (iflag_debug .gt. 0) write(*,*) 'const_mesh_infos'
      call const_mesh_infos(my_rank, mesh, group, ele_mesh)
!
      if(iflag_debug.gt.0) write(*,*)' const_element_comm_tbls'
      call const_element_comm_tbls(mesh, ele_mesh)
!
      if(i_debug .eq. iflag_full_msg) then
        call check_whole_num_of_elements(mesh%ele)
      end if
!
!     ---------------------
!
      iflag = SGS_par%filter_p%iflag_SGS_filter
      if(     (SGS_par%model_p%iflag_dynamic .ne. id_SGS_DYNAMIC_OFF    &
     &    .or. SGS_par%model_p%iflag_SGS.eq.id_SGS_similarity)) then
!
        if   (iflag .eq. id_SGS_3D_FILTERING                            &
     &   .or. iflag .eq. id_SGS_3D_EZ_FILTERING) then
          if (iflag_debug .gt. 0)                                       &
     &      write(*,*)' s_set_istart_3d_filtering'
          call s_set_istart_3d_filtering(filtering1%filter)
!
        else if (iflag.eq.id_SGS_3D_SMP_FILTERING                       &
     &     .or. iflag.eq.id_SGS_3D_EZ_SMP_FILTERING) then
          if (iflag_debug .gt. 0)                                       &
     &      write(*,*) ' const_tbl_3d_filtering_smp'
          call const_tbl_3d_filtering_smp(filtering1)
!
        else if (iflag .eq. id_SGS_LINE_FILTERING) then
          if (iflag_debug.gt.0) write(*,*)' ordering_l_filter_smp'
          call ordering_l_filter_smp(mesh%node%istack_nod_smp,          &
     &        filtering1%fil_l, filtering1%fil_l_smp)
        end if
      end if
!
      if(      (SGS_par%model_p%iflag_dynamic .ne. id_SGS_DYNAMIC_OFF   &
     &    .and. SGS_par%model_p%iflag_SGS.eq.id_SGS_similarity) ) then
        if    (iflag .eq. id_SGS_3D_FILTERING                           &
     &    .or. iflag .eq. id_SGS_3D_EZ_FILTERING) then
          if (iflag_debug .gt. 0) write(*,*)' s_set_istart_w_filtering'
          call s_set_istart_3d_filtering(wide_filtering%filter)
!
        else if (iflag.eq.id_SGS_3D_SMP_FILTERING                       &
     &     .or. iflag.eq.id_SGS_3D_EZ_SMP_FILTERING) then
!
          if (iflag_debug .gt. 0)                                       &
     &      write(*,*) 'const_tbl_3d_filtering_smp'
          call const_tbl_3d_filtering_smp(wide_filtering)
        end if
      end if
!
!     ---------------------
!
      if (iflag_debug.eq.1) write(*,*)' allocate_array'
      call allocate_array(SGS_par, mesh, MHD_prop,                      &
     &    iphys, nod_fld, Csims_FEM_MHD%iphys_elediff,                  &
     &    mhd_fem1_wk, rhs_mat1, fem_int1, label_sim)
!
      if (iflag_debug.eq.1) write(*,*)' set_reference_temp'
      call set_reference_temp                                           &
     &   (MHD_prop%ref_param_T, MHD_prop%takepito_T, mesh%node,         &
     &    MHD_mesh%fluid, iphys%i_ref_t, iphys%i_gref_t, nod_fld)
      call set_reference_temp                                           &
     &   (MHD_prop%ref_param_C, MHD_prop%takepito_C, mesh%node,         &
     &    MHD_mesh%fluid, iphys%i_ref_c, iphys%i_gref_c, nod_fld)
!
      if (iflag_debug.eq.1) write(*,*)' set_material_property'
      call set_material_property                                        &
     &   (iphys, MHD_prop%ref_param_T%depth_top,                        &
     &    MHD_prop%ref_param_T%depth_bottom, MHD_prop)
      call s_init_ele_material_property                                 &
     &   (mesh%ele%numele, MHD_prop, ak_MHD)
!
      call define_sgs_components                                        &
     &   (mesh%node%numnod, mesh%ele%numele,                            &
     &    SGS_par%model_p, layer_tbl, MHD_prop,                         &
     &    Csims_FEM_MHD%ifld_sgs, Csims_FEM_MHD%icomp_sgs, wk_sgs1,     &
     &    Csims_FEM_MHD%sgs_coefs, Csims_FEM_MHD%sgs_coefs_nod)
      call define_sgs_diff_coefs(mesh%ele%numele,                       &
     &    SGS_par%model_p, SGS_par%commute_p, layer_tbl, MHD_prop,      &
     &    Csims_FEM_MHD%ifld_diff, Csims_FEM_MHD%icomp_diff,            &
     &    wk_diff1, Csims_FEM_MHD%diff_coefs)
!
      call deallocate_surface_geom_type(ele_mesh%surf)
      call deallocate_edge_geom_type(ele_mesh%edge)
!
!     --------------------- 
!
      if (rst_step%increment .gt. 0) then
        if (iflag_debug.eq.1) write(*,*)' init_restart_4_snapshot'
        call init_restart_4_snapshot(MHD_step%init_d%i_time_step,       &
     &      fst_file_IO, mesh%node, t_IO, rst_step)
      end if
!
!     ---------------------
!
      call const_bc_infinity_surf_grp                                   &
     &   (iflag_surf_infty, group%surf_grp, group%infty_grp)
!
!     --------------------- 
!
      if (iflag_debug.eq.1) write(*,*)' const_MHD_jacobian_and_volumes'
      call const_MHD_jacobian_and_volumes(SGS_par%model_p,              &
     &    ele_mesh, group, mesh, layer_tbl, fem_int1%jcs,               &
     &    MHD_mesh)
!
!     --------------------- 
!
      if (iflag_debug.eq.1) write(*,*)' set_connect_RHS_assemble'
      call s_set_RHS_assemble_table                                     &
     &   (mesh%node, mesh%ele, fem_int1%next_tbl, fem_int1%rhs_tbl)
!
!     ---------------------
!
      if (iflag_debug.eq.1) write(*,*)  'const_normal_vector'
      call const_normal_vector                                          &
     &   (my_rank, nprocs, mesh%node, ele_mesh%surf,                    &
     &    fem_int1%jcs)
!
      if (iflag_debug.eq.1) write(*,*)' int_surface_parameters'
      call int_surface_parameters                                       &
     &   (mesh, ele_mesh%surf, group, rhs_mat1%surf_wk)
!
!     --------------------- 
!
      if (iflag_debug.eq.1) write(*,*)' set_boundary_data'
      call set_boundary_data                                            &
     &   (MHD_step%time_d, IO_bc, mesh, ele_mesh, MHD_mesh, group,      &
     &    MHD_prop, MHD_BC1, iphys, nod_fld)
!
!     ---------------------
!
      call int_RHS_mass_matrices(FEM_prm%npoint_t_evo_int,              &
     &     mesh, MHD_mesh, fem_int1%jcs, fem_int1%rhs_tbl,              &
     &     mhd_fem1_wk, rhs_mat1%fem_wk, rhs_mat1%f_l, fem_int1%m_lump)
!
!     ---------------------
!
      call deallocate_surf_bc_lists(MHD_prop, MHD_BC1)
!
      end subroutine init_analyzer_snap
!
! ----------------------------------------------------------------------
!
      end module initialize_4_snapshot
