!
!     module initialize_4_MHD_AMG
!
!        programmed H.Matsui on Dec., 2008
!
!!      subroutine s_initialize_4_MHD_AMG                               &
!!     &         (dt, FEM_prm, node_1st, ele_1st, ifld_diff, diff_coefs,&
!!     &          DJDS_param, MGCG_WK, MGCG_FEM, MGCG_MHD_FEM,          &
!!     &          MHD_matrices)
!!        type(FEM_MHD_paremeters), intent(in) :: FEM_prm
!!        type(node_data), intent(inout) :: node_1st
!!        type(element_data), intent(inout) :: ele_1st
!!      subroutine const_MGCG_MHD_matrices                              &
!!     &         (dt, FEM_prm, SGS_param, cmt_param, ifld_diff,         &
!!     &          MGCG_WK, MGCG_FEM, MGCG_MHD_FEM, MHD_matrices)
!!        type(FEM_MHD_paremeters), intent(in) :: FEM_prm
!!        type(SGS_model_control_params), intent(in) :: SGS_param
!!        type(commutation_control_params), intent(in) :: cmt_param
!!        type(SGS_terms_address), intent(in) :: ifld_diff
!!        type(SGS_coefficients_type), intent(in) :: diff_coefs
!!        type(DJDS_poarameter), intent(in) :: DJDS_param
!!        type(MGCG_data), intent(inout) :: MGCG_WK
!!        type(mesh_4_MGCG), intent(inout) :: MGCG_FEM
!!        type(MGCG_MHD_data), intent(inout) :: MGCG_MHD_FEM
!!        type(MHD_MG_matrices), intent(inout) :: MHD_matrices
!
      module initialize_4_MHD_AMG
!
      use m_precision
!
      use m_machine_parameter
!
      use m_physical_property
      use t_FEM_control_parameter
      use t_iccg_parameter
      use t_solver_djds_MHD
      use t_material_property
      use t_SGS_model_coefs
      use t_next_node_ele_4_node
      use t_MGCG_data
      use t_MGCG_data_4_MHD
!
      use calypso_mpi
!
      implicit none
!
! ---------------------------------------------------------------------
!
      contains
!
! ---------------------------------------------------------------------
!
      subroutine s_initialize_4_MHD_AMG                                 &
     &         (dt, FEM_prm, node_1st, ele_1st, ifld_diff, diff_coefs,  &
     &          DJDS_param, MGCG_WK, MGCG_FEM, MGCG_MHD_FEM,            &
     &          MHD_matrices)
!
      use t_geometry_data
      use t_edge_data
      use t_surface_data
      use t_bc_data_MHD
      use t_jacobians
!
      use m_boundary_condition_IDs
      use set_layers_4_MHD_AMG
      use const_mesh_information
      use allocate_MHD_AMG_array
      use set_diffusivities_MHD_AMG
      use const_comm_tbl_type_fluid
      use const_bc_infty_surf_type
      use set_table_type_RHS_assemble
      use set_djds_connectivity_type
      use set_djds_connect_type_MHD
      use set_normal_vectors
      use int_surface_params_MHD
      use set_surface_id_MHD
      use int_type_mass_matrices
      use set_MHD_idx_4_mat_type
      use link_MG_MHD_mesh_data
      use const_element_comm_tables
!
      real(kind = kreal), intent(in) :: dt
      type(FEM_MHD_paremeters), intent(in) :: FEM_prm
      type(SGS_terms_address), intent(in) :: ifld_diff
      type(SGS_coefficients_type), intent(in) :: diff_coefs
      type(DJDS_poarameter), intent(in) :: DJDS_param
!
      type(node_data), intent(inout) :: node_1st
      type(element_data), intent(inout) :: ele_1st
      type(MGCG_data), intent(inout) :: MGCG_WK
      type(mesh_4_MGCG), intent(inout) :: MGCG_FEM
      type(MGCG_MHD_data), intent(inout) :: MGCG_MHD_FEM
      type(MHD_MG_matrices), intent(inout) :: MHD_matrices
!
      integer(kind = kint) :: i_level
!
!
      call split_multigrid_comms(MGCG_WK)
!
      if (iflag_debug .gt. 0) write(*,*) 'alloc_iccgN_vec_type'
      MGCG_WK%MG_vector(0)%isize_solver_vect = -1
      call alloc_iccgN_vec_type                                         &
     &   (isix, node_1st%numnod,  MGCG_WK%MG_vector(0))
!
!     --------------------- 
!
!
      do i_level = 1, MGCG_WK%num_MG_level
        if(my_rank .lt. MGCG_WK%MG_mpi(i_level)%nprocs ) then
          if(iflag_debug .gt. 0) write(*,*)                             &
     &            'set_layers_type_4_MHD', i_level
          call set_layers_type_4_MHD                                    &
     &       (FEM_prm, MGCG_FEM%MG_mesh(i_level)%mesh,                  &
     &        MGCG_FEM%MG_mesh(i_level)%group,                          &
     &        MGCG_MHD_FEM%MG_MHD_mesh(i_level) )
          if(iflag_debug .gt. 0) write(*,*)                             &
     &            'const_mesh_infos', i_level
          call const_mesh_infos                                         &
     &       (my_rank, MGCG_FEM%MG_mesh(i_level)%mesh,                  &
     &        MGCG_FEM%MG_mesh(i_level)%group,                          &
     &        MGCG_FEM%MG_ele_mesh(i_level))
        else
          call set_empty_layers_type_4_MHD                              &
     &       (MGCG_MHD_FEM%MG_MHD_mesh(i_level) )
          call empty_mesh_info                                          &
     &       (MGCG_FEM%MG_mesh(i_level)%mesh,                           &
     &        MGCG_FEM%MG_mesh(i_level)%group,                          &
     &        MGCG_FEM%MG_ele_mesh(i_level))
        end if
!
        call deallocate_edge_geom_type                                  &
     &     (MGCG_FEM%MG_ele_mesh(i_level)%edge)
      end do
!
!     ---------------------
!
      do i_level = 1, MGCG_WK%num_MG_level
        if(iflag_debug .gt. 0) write(*,*)                               &
     &            's_allocate_MHD_AMG_array', i_level
        call s_allocate_MHD_AMG_array                                   &
     &     (MGCG_FEM%MG_mesh(i_level), MGCG_FEM%MG_ele_mesh(i_level),   &
     &      MGCG_WK%MG_vector(i_level), MGCG_FEM%MG_FEM_mat(i_level),   &
     &      MGCG_FEM%MG_FEM_int(i_level),                               &
     &      MGCG_MHD_FEM%MG_mk_MHD(i_level) )
      end do
!
!     --------------------- 
!
      do i_level = 1, MGCG_WK%num_MG_level
        if(iflag_debug .gt. 0) write(*,*)                               &
     &            's_set_diffusivities_MHD_AMG', i_level
        call s_set_diffusivities_MHD_AMG                                &
     &    (MGCG_FEM%MG_mesh(i_level)%mesh%ele,                          &
     &     fl_prop1, cd_prop1, ht_prop1, cp_prop1,                      &
     &     MGCG_MHD_FEM%ak_MHD_AMG(i_level))
        if(iflag_debug .gt. 0) write(*,*)                               &
     &            's_set_sgs_diff_array_MHD_AMG', i_level
!
        call copy_SGS_num_coefs                                         &
     &     (diff_coefs, MGCG_MHD_FEM%MG_diff_coefs(i_level))
        call alloc_SGS_coefs                                            &
     &     (MGCG_FEM%MG_mesh(i_level)%mesh%ele%numele,                  &
     &      MGCG_MHD_FEM%MG_diff_coefs(i_level))
      end do
!
!     --------------------- 
!
      do i_level = 1, MGCG_WK%num_MG_level
        if(iflag_debug .gt. 0) write(*,*)                               &
     &            's_const_comm_tbl_type_fluid', i_level
        call s_const_comm_tbl_type_fluid(MGCG_WK%MG_mpi(i_level),       &
     &      MGCG_FEM%MG_mesh(i_level)%mesh,                             &
     &      MGCG_MHD_FEM%MG_MHD_mesh(i_level) )
!
        call const_element_comm_tbls                                    &
     &     (MGCG_FEM%MG_mesh(i_level)%mesh,                             &
     &      MGCG_FEM%MG_ele_mesh(i_level))
      end do
!
!     ---------------------
!
      do i_level = 1, MGCG_WK%num_MG_level
        if(my_rank .lt. MGCG_WK%MG_mpi(i_level)%nprocs ) then
          if(iflag_debug .gt. 0) write(*,*)                             &
     &            'const_bc_infinity_surf_grp', i_level
          call const_bc_infinity_surf_grp(iflag_surf_infty,             &
     &        MGCG_FEM%MG_mesh(i_level)%group%surf_grp,                 &
     &        MGCG_FEM%MG_mesh(i_level)%group%infty_grp)
        else
          call empty_infty_surf_type(MGCG_FEM%MG_mesh(i_level)%group)
        end if
      end do
!
!     --------------------- 
!
      do i_level = 1, MGCG_WK%num_MG_level
        call const_jacobians_element                                    &
     &     (my_rank, MGCG_WK%MG_mpi(i_level)%nprocs,                    &
     &      MGCG_FEM%MG_mesh(i_level)%mesh%node,                        &
     &      MGCG_FEM%MG_mesh(i_level)%mesh%ele,                         &
     &      MGCG_FEM%MG_mesh(i_level)%group%surf_grp,                   &
     &      MGCG_FEM%MG_mesh(i_level)%group%infty_grp,                  &
     &      MGCG_FEM%MG_FEM_int(i_level)%jacobians)
!
        call const_jacobians_surf_group                                 &
     &     (my_rank, MGCG_WK%MG_mpi(i_level)%nprocs,                    &
     &      MGCG_FEM%MG_mesh(i_level)%mesh%node,                        &
     &      MGCG_FEM%MG_mesh(i_level)%mesh%ele,                         &
     &      MGCG_FEM%MG_ele_mesh(i_level)%surf,                         &
     &      MGCG_FEM%MG_mesh(i_level)%group%surf_grp,                   &
     &      MGCG_FEM%MG_FEM_int(i_level)%jacobians)
      end do
!
!
!     -----  set DJDS matrix connectivity
!
      if(iflag_debug .gt. 0) write(*,*) 'set_MG_djds_connect_type'
      call set_MG_djds_connect_type(DJDS_param,                         &
     &    MGCG_WK, MGCG_MHD_FEM, MGCG_FEM, MHD_matrices)
!
!     --------------------- 
!
      if(iflag_debug .gt. 0) write(*,*) 's_link_MG_MHD_mesh_data'
      call s_link_MG_MHD_mesh_data                                      &
     &   (MGCG_WK, MGCG_FEM%MG_mesh, MGCG_MHD_FEM%MG_MHD_mesh,          &
     &    ele_1st, MHD_matrices)
!
!     --------------------- 
!
      do i_level = 1, MGCG_WK%num_MG_level
        if(iflag_debug .gt. 0) write(*,*)                               &
     &         'int_normal_4_all_surface', i_level
        call int_normal_4_all_surface                                   &
     &     (MGCG_FEM%MG_ele_mesh(i_level)%surf,                         &
     &      MGCG_FEM%MG_FEM_int(i_level)%jacobians%jac_2d)
        call int_surface_parameters(MGCG_FEM%MG_mesh(i_level)%mesh,     &
     &      MGCG_FEM%MG_ele_mesh(i_level)%surf,                         &
     &      MGCG_FEM%MG_mesh(i_level)%group,                            &
     &      MGCG_FEM%MG_FEM_mat(i_level)%surf_wk)
      end do
!
!     --------------------- 
!
      do i_level = 1, MGCG_WK%num_MG_level
        call set_bc_id_data                                             &
     &     (dt, MGCG_MHD_FEM%IO_MG_bc(i_level),                         &
     &      MGCG_FEM%MG_mesh(i_level)%mesh,                             &
     &      MGCG_FEM%MG_mesh(i_level)%group,                            &
     &      MGCG_MHD_FEM%MG_MHD_mesh(i_level),                          &
     &      fl_prop1, cd_prop1, ht_prop1, cp_prop1,                     &
     &      MGCG_MHD_FEM%MG_node_bc(i_level))
!
        call set_bc_surface_data(MGCG_MHD_FEM%IO_MG_bc(i_level),        &
     &      MGCG_FEM%MG_mesh(i_level)%mesh%node,                        &
     &      MGCG_FEM%MG_mesh(i_level)%mesh%ele,                         &
     &      MGCG_FEM%MG_ele_mesh(i_level)%surf,                         &
     &      MGCG_FEM%MG_mesh(i_level)%group%surf_grp,                   &
     &      MGCG_FEM%MG_mesh(i_level)%group%surf_nod_grp,               &
     &      MGCG_FEM%MG_mesh(i_level)%group%surf_grp_geom,              &
     &      fl_prop1, cd_prop1, ht_prop1, cp_prop1,                     &
     &      MGCG_MHD_FEM%MG_surf_bc(i_level) )
      end do
!
!     --------------------- 
!
      do i_level = 1, MGCG_WK%num_MG_level
        if(iflag_debug .gt. 0) write(*,*) 's_int_type_mass_matrices'
        call s_int_type_mass_matrices                                   &
     &     (FEM_prm, MGCG_FEM%MG_mesh(i_level)%mesh,                    &
     &      MGCG_MHD_FEM%MG_MHD_mesh(i_level),                          &
     &      MGCG_FEM%MG_FEM_int(i_level)%jacobians,                     &
     &      MGCG_FEM%MG_FEM_int(i_level)%rhs_tbl,                       &
     &      MGCG_FEM%MG_FEM_mat(i_level),                               &
     &      MGCG_FEM%MG_FEM_int(i_level),                               &
     &      MGCG_MHD_FEM%MG_mk_MHD(i_level) )
      end do
!
!     ---------------------
!
      do i_level = 1, MGCG_WK%num_MG_level
        if(iflag_debug .gt. 0) write(*,*) 's_set_MHD_idx_4_mat_type'
        call s_set_MHD_idx_4_mat_type                                   &
     &     (MGCG_FEM%MG_mesh(i_level)%mesh,                             &
     &      MGCG_MHD_FEM%MG_MHD_mesh(i_level),                          &
     &      fl_prop1, cd_prop1, ht_prop1, cp_prop1,                     &
     &      MGCG_FEM%MG_FEM_int(i_level)%rhs_tbl,                       &
     &      MHD_matrices%MG_DJDS_table(i_level),                        &
     &      MHD_matrices%MG_DJDS_fluid(i_level),                        &
     &      MHD_matrices%MG_DJDS_linear(i_level),                       &
     &      MHD_matrices%MG_DJDS_lin_fl(i_level),                       &
     &      MHD_matrices%MG_mat_tbls(i_level)%base,                     &
     &      MHD_matrices%MG_mat_tbls(i_level)%fluid_q,                  &
     &      MHD_matrices%MG_mat_tbls(i_level)%full_conduct_q,           &
     &      MHD_matrices%MG_mat_tbls(i_level)%linear,                   &
     &      MHD_matrices%MG_mat_tbls(i_level)%fluid_l)
      end do
!
!     ---------------------
!
      do i_level = 1, MGCG_WK%num_MG_level
        if(my_rank .lt. MGCG_WK%MG_mpi(i_level)%nprocs) then
          if(iflag_debug .gt. 0) write(*,*) 'alloc_aiccg_matrices'
          call alloc_aiccg_matrices                                     &
     &       (MGCG_FEM%MG_mesh(i_level)%mesh%node,                      &
     &        fl_prop1, cd_prop1, ht_prop1, cp_prop1,                   &
     &        MHD_matrices%MG_DJDS_table(i_level),                      &
     &        MHD_matrices%MG_DJDS_fluid(i_level),                      &
     &        MHD_matrices%MG_DJDS_linear(i_level),                     &
     &        MHD_matrices%MG_DJDS_lin_fl(i_level),                     &
     &        MHD_matrices%Vmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Bmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Tmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Cmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Pmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Fmat_MG_DJDS(i_level) )
        else
          if(iflag_debug .gt. 0) write(*,*) 'alloc_MG_zero_matrices'
          call alloc_MG_zero_matrices                                   &
     &       (fl_prop1, cd_prop1, ht_prop1, cp_prop1,                   &
     &        MHD_matrices%Vmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Bmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Tmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Cmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Pmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Fmat_MG_DJDS(i_level) )
        end if
      end do
!
      end subroutine s_initialize_4_MHD_AMG
!
! ---------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine set_MG_djds_connect_type(DJDS_param,                   &
     &          MGCG_WK, MGCG_MHD_FEM, MGCG_FEM, MHD_matrices)
!
      use t_mesh_data
      use t_geometry_data_MHD
      use t_next_node_ele_4_node
      use set_djds_connect_type_MHD
!
      type(DJDS_poarameter), intent(in) :: DJDS_param
      type(MGCG_data), intent(in) :: MGCG_WK
      type(MGCG_MHD_data), intent(in) :: MGCG_MHD_FEM
      type(mesh_4_MGCG), intent(inout) :: MGCG_FEM
      type(MHD_MG_matrices), intent(inout) :: MHD_matrices
!
      integer(kind = kint) :: i_level
!
!
      do i_level = 1, MGCG_WK%num_MG_level
        if(my_rank .ge. MGCG_WK%MG_mpi(i_level)%nprocs ) then
          call empty_whole_djds_connectivity                            &
     &       (MGCG_FEM%MG_mesh(i_level)%mesh,                           &
     &        MGCG_FEM%MG_FEM_int(i_level)%next_tbl,                    &
     &        MGCG_FEM%MG_FEM_int(i_level)%rhs_tbl,                     &
     &        MHD_matrices%MG_DJDS_table(i_level))
          call empty_MHD_djds_connectivities                            &
     &       (MGCG_FEM%MG_mesh(i_level)%mesh,                           &
     &        MHD_matrices%MG_DJDS_fluid(i_level),                      &
     &        MHD_matrices%MG_DJDS_linear(i_level),                     &
     &        MHD_matrices%MG_DJDS_lin_fl(i_level))
        else
          call set_MHD_whole_connectivity(DJDS_param,                   &
     &        MGCG_FEM%MG_mesh(i_level)%mesh, MGCG_WK%MG_mpi(i_level),  &
     &        MGCG_FEM%MG_FEM_int(i_level)%next_tbl,                    &
     &        MGCG_FEM%MG_FEM_int(i_level)%rhs_tbl,                     &
     &        MHD_matrices%MG_DJDS_table(i_level),                      &
     &        MHD_matrices%MG_comm_table(i_level))
!
          call set_MHD_djds_connectivities(DJDS_param,                  &
     &        MGCG_FEM%MG_mesh(i_level)%mesh,                           &
     &        MGCG_MHD_FEM%MG_MHD_mesh(i_level)%fluid,                  &
     &        MGCG_MHD_FEM%MG_MHD_mesh(i_level)%nod_fl_comm,            &
     &        MGCG_WK%MG_mpi(i_level),                                  &
     &        MHD_matrices%MG_DJDS_table(i_level),                      &
     &        MHD_matrices%MG_DJDS_fluid(i_level),                      &
     &        MHD_matrices%MG_DJDS_linear(i_level),                     &
     &        MHD_matrices%MG_DJDS_lin_fl(i_level))
        end if
!
        call link_comm_tbl_types                                        &
     &     (MGCG_MHD_FEM%MG_MHD_mesh(i_level)%nod_fl_comm,              &
     &      MHD_matrices%MG_comm_fluid(i_level))
!
        call dealloc_iele_belonged                                      &
     &    (MGCG_FEM%MG_FEM_int(i_level)%next_tbl%neib_ele)
        call dealloc_inod_next_node                                     &
     &    (MGCG_FEM%MG_FEM_int(i_level)%next_tbl%neib_nod)
      end do
!
      end subroutine set_MG_djds_connect_type
!
!-----------------------------------------------------------------------
!
      subroutine const_MGCG_MHD_matrices                                &
     &         (dt, FEM_prm, SGS_param, cmt_param, ifld_diff,           &
     &          MGCG_WK, MGCG_FEM, MGCG_MHD_FEM, MHD_matrices)
!
      use t_FEM_control_parameter
      use t_SGS_control_parameter
      use set_aiccg_matrices_type
      use precond_djds_MHD
!
      real(kind = kreal), intent(in) :: dt
      type(FEM_MHD_paremeters), intent(in) :: FEM_prm
      type(SGS_model_control_params), intent(in) :: SGS_param
      type(commutation_control_params), intent(in) :: cmt_param
      type(SGS_terms_address), intent(in) :: ifld_diff
      type(MGCG_data), intent(in) :: MGCG_WK
!
      type(mesh_4_MGCG), intent(inout) :: MGCG_FEM
      type(MGCG_MHD_data), intent(inout) :: MGCG_MHD_FEM
      type(MHD_MG_matrices), intent(inout) :: MHD_matrices
!
      integer(kind = kint) :: i_level
!
!
      do i_level = 1, MGCG_WK%num_MG_level
        if(my_rank .lt. MGCG_WK%MG_mpi(i_level)%nprocs) then
          if (iflag_debug.eq.1) write(*,*) 'set MG matrices', i_level
          call s_set_aiccg_matrices                                     &
     &       (dt, FEM_prm, SGS_param, cmt_param,                        &
     &        MGCG_FEM%MG_mesh(i_level)%mesh,                           &
     &        MGCG_FEM%MG_mesh(i_level)%group,                          &
     &        MGCG_FEM%MG_ele_mesh(i_level),                            &
     &        MGCG_MHD_FEM%MG_MHD_mesh(i_level),                        &
     &        MGCG_MHD_FEM%MG_node_bc(i_level),                         &
     &        MGCG_MHD_FEM%MG_surf_bc(i_level),                         &
     &        fl_prop1, cd_prop1, ht_prop1, cp_prop1,                   &
     &        MGCG_MHD_FEM%ak_MHD_AMG(i_level),                         &
     &        MGCG_FEM%MG_FEM_int(i_level)%jacobians,                   &
     &        MGCG_MHD_FEM%MG_filter_MHD(i_level), ifld_diff,           &
     &        MGCG_MHD_FEM%MG_diff_coefs(i_level),                      &
     &        MGCG_FEM%MG_FEM_int(i_level)%rhs_tbl,                     &
     &        MHD_matrices%MG_DJDS_table(i_level),                      &
     &        MHD_matrices%MG_DJDS_fluid(i_level),                      &
     &        MHD_matrices%MG_DJDS_linear(i_level),                     &
     &        MHD_matrices%MG_DJDS_lin_fl(i_level),                     &
     &        MHD_matrices%MG_mat_tbls(i_level)%base,                   &
     &        MHD_matrices%MG_mat_tbls(i_level)%fluid_q,                &
     &        MHD_matrices%MG_mat_tbls(i_level)%full_conduct_q,         &
     &        MHD_matrices%MG_mat_tbls(i_level)%linear,                 &
     &        MHD_matrices%MG_mat_tbls(i_level)%fluid_l,                &
     &        MGCG_MHD_FEM%MG_mk_MHD(i_level)%fluid,                    &
     &        MGCG_MHD_FEM%MG_mk_MHD(i_level)%conduct,                  &
     &        MGCG_FEM%MG_FEM_mat(i_level)%surf_wk,                     &
     &        MGCG_FEM%MG_FEM_mat(i_level)%fem_wk,                      &
     &        MHD_matrices%Vmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Bmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Tmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Cmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Pmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Fmat_MG_DJDS(i_level))
        end if
      end do
!
      do i_level = 1, MGCG_WK%num_MG_level
        if(my_rank .lt. MGCG_WK%MG_mpi(i_level)%nprocs) then
          if (iflag_debug.gt.0) write(*,*) 'preconditioning', i_level
          call matrix_precondition                                      &
     &       (FEM_prm%MG_param%PRECOND_MG, FEM_prm%MG_param%PRECOND_MG, &
     &        FEM_prm%CG11_param%sigma_diag,                            &
     &        fl_prop1, cd_prop1, ht_prop1, cp_prop1,                   &
     &        MHD_matrices%MG_DJDS_table(i_level),                      &
     &        MHD_matrices%MG_DJDS_fluid(i_level),                      &
     &        MHD_matrices%MG_DJDS_linear(i_level),                     &
     &        MHD_matrices%MG_DJDS_lin_fl(i_level),                     &
     &        MHD_matrices%Vmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Bmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Tmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Cmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Pmat_MG_DJDS(i_level),                       &
     &        MHD_matrices%Fmat_MG_DJDS(i_level))
        end if
      end do
!
      end subroutine const_MGCG_MHD_matrices
!
! ---------------------------------------------------------------------
!
      end module initialize_4_MHD_AMG
