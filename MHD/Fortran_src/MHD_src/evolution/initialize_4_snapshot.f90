!initialize_4_snapshot.f90
!     module initialize_4_snapshot
!
!      Written by H. Matsui
!
!------- subroutine init_analyzer_snap ---------------------
!
      module initialize_4_snapshot
!
      use m_precision
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine init_analyzer_snap
!
      use m_parallel_var_dof
      use m_machine_parameter
      use m_control_parameter
      use m_t_step_parameter
!
      use m_geometry_data
      use m_surface_geometry_data
      use m_surf_data_infinity
      use m_edge_geometry_data
      use m_node_phys_address
      use m_ele_material_property
      use m_bulk_values
      use m_jacobians
      use m_work_4_dynamic_model
      use m_work_layer_correlate
!
      use m_check_subroutines
!
      use const_mesh_info
      use cal_mesh_position
      use count_whole_num_element
!
      use cal_jacobian
      use cal_volume_node_MHD
      use int_MHD_mass_matrices
      use int_surface_params_MHD
      use set_nodal_bc_id_data
      use set_surface_bc_data
      use allocate_array_MHD
      use const_RHS_assemble_list
      use ordering_line_filter_smp
      use const_ele_layering_table
      use const_comm_table_fluid
      use set_reference_value
      use material_property
      use reordering_by_layers
      use set_layers_4_MHD
      use const_tbl_3d_filtering_smp
      use const_tbl_w_filtering_smp
      use set_istart_3d_filtering
      use count_sgs_components
      use set_layer_list_by_table
      use set_normal_vectors
      use mhd_restart_file_IO_control
!
      use nodal_vector_send_recv
!
!     --------------------- 
!
      if (iflag_debug.eq.1) write(*,*)' s_reordering_by_layers_snap'
      call s_reordering_by_layers_snap
!
      if (iflag_debug.eq.1) write(*,*)' set_layers'
      call set_layers
!      call check_numbers_of_nodes(my_rank)
!      call check_nodes_4_layers(my_rank)
!
      if (iflag_dynamic_SGS  .ne. id_SGS_DYNAMIC_OFF) then
        ncomp_correlate = 9
        if (iflag_debug.eq.1) write(*,*)' const_layers_4_dynamic'
        call const_layers_4_dynamic
        call allocate_work_4_dynamic
        call allocate_work_layer_correlate
      end if
!
!
      if (iflag_debug.eq.1) write(*,*)' const_mesh_informations'
! 
      call const_mesh_informations(my_rank)
!
      call deallocate_surface_geometry
      call deallocate_edge_geometry
!
      if(i_debug .eq. iflag_full_msg) call check_whole_num_of_elements
!
!     ---------------------
!
      if( (iflag_dynamic_SGS .ne. id_SGS_DYNAMIC_OFF                    &
     &       .or. iflag_SGS_model.eq.id_SGS_similarity)) then
!
        if (iflag_SGS_filter.eq.1 .or. iflag_SGS_filter.eq.11) then
          if (iflag_debug .gt. 0)                                       &
     &      write(*,*)' s_set_istart_3d_filtering'
          call s_set_istart_3d_filtering
!
        else if (iflag_SGS_filter.eq.21                                 &
     &     .or. iflag_SGS_filter.eq.31) then
          if (iflag_debug .gt. 0)                                       &
     &      write(*,*)' s_const_tbl_3d_filtering_smp'
          call s_const_tbl_3d_filtering_smp
!
        else if (iflag_SGS_filter.eq.2) then
          if (iflag_debug.gt.0) write(*,*)' ordering_l_filter_smp'
          call ordering_l_filter_smp
        end if
      end if
!
      if( (iflag_dynamic_SGS .ne. id_SGS_DYNAMIC_OFF                    &
     &      .and. iflag_SGS_model.eq.id_SGS_similarity) ) then
        if (iflag_SGS_filter.eq.1 .or. iflag_SGS_filter.eq.11) then
          if (iflag_debug .gt. 0) write(*,*)' s_set_istart_w_filtering'
          call s_set_istart_w_filtering
!
        else if (iflag_SGS_filter.eq.21                                 &
     &     .or. iflag_SGS_filter.eq.31) then
!
          if (iflag_debug .gt. 0)                                       &
     &      write(*,*)' s_const_tbl_w_filtering_smp'
          call s_const_tbl_w_filtering_smp
        end if
      end if
!
!     ---------------------
!
      if (iflag_debug.eq.1) write(*,*)' allocate_array'
      call allocate_array
!
      if (iflag_debug.eq.1) write(*,*)' set_reference_temp'
      call set_reference_temp
!
      if (iflag_debug.eq.1) write(*,*)' set_material_property'
      call set_material_property
      call init_ele_material_property
      call s_count_sgs_components
!
      call init_send_recv
!
      if (iflag_debug.gt.0)  write(*,*)' make comm. table for fluid'
      call s_const_comm_table_fluid
!
!     --------------------- 
!
      if (i_step_output_rst .gt. 0) then
        if (iflag_debug.eq.1) write(*,*)' init_restart_4_snapshot'
        call init_restart_4_snapshot
      end if
!
!     ---------------------
!
      call const_bc_infinity_surf_grp
!
!     --------------------- 
!
      call set_max_int_point_by_etype
      call cal_jacobian_element
      call cal_jacobian_surface
!
      call deallocate_dxi_dx_quad
      call deallocate_dxi_dx_linear
!
      call cal_jacobian_surf_grp
!
!     --------------------- 
!
      if (iflag_debug.eq.1) write(*,*)' set_connect_RHS_assemble'
      call set_connect_RHS_assemble
!
!     --------------------- 
!
      if (iflag_debug.eq.1) write(*,*)' cal_volume_node'
      call cal_volume_node
!
!     --------------------- 
!
      if (iflag_debug.gt.0) write(*,*) 's_cal_normal_vector'
      call s_cal_normal_vector
!
      if (iflag_debug.eq.1) write(*,*)' int_surface_parameters'
      call int_surface_parameters
!
!     --------------------- 
!
      if (iflag_debug.eq.1) write(*,*)' set_bc_id_data'
      call set_bc_id_data 
!
      if (iflag_debug.eq.1) write(*,*)' set_surf_bc_data'
      call set_surf_bc_data 
      call deallocate_surf_bc_lists
!
!     --------------------- 
!
      call int_mass_matrices
!
!     --------------------- 
!
      call time_prog_barrier
!
!     --------------------- 
!
       if (my_rank.eq.0) write(*,*)' end init_analyzer_snap'
!
!
      end subroutine init_analyzer_snap
!
! ----------------------------------------------------------------------
!
      end module initialize_4_snapshot
