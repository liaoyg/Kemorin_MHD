!analyzer_udt_correlation.f90
!      module analyzer_udt_correlation
!
!      subroutine initialize_udt_correlate
!      subroutine analyze_udt_correlate
!
!..................................................
!
!      modified by H. Matsui on Nov., 2006 
!
      module analyzer_udt_correlation
!
      use m_precision
      use m_constants
      use m_machine_parameter
      use calypso_mpi
!
      use t_mesh_data
      use t_phys_data
!
      use transfer_correlate_field
!
      implicit none
!
      type(mesh_geometry), save :: mesh_ref
      type(mesh_groups), save :: group_ref
      type(phys_data), save :: phys_ref
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine initialize_udt_correlate
!
      use m_array_for_send_recv
      use m_nod_comm_table
      use m_geometry_data
      use m_group_data
      use m_jacobians
      use m_layering_ele_list
      use m_node_phys_address
      use m_node_phys_data
      use m_ele_sf_eg_comm_tables
!
      use copy_mesh_structures
      use input_control_udt_diff
      use const_mesh_info
      use nod_phys_send_recv
!
      use m_2nd_pallalel_vector
      use const_ele_layering_table
      use int_volume_of_domain
      use correlation_all_layerd_data
!
      use m_jacobians
!
!
      if (my_rank.eq.0) then
        write(*,*) 'diff. udt files'
        write(*,*) 'Input file: mesh data, udt data'
      end if
!
!     --------------------- 
!
      if (iflag_debug.eq.1) write(*,*) 's_input_control_corr_udt'
      call s_input_control_corr_udt
      if (iflag_debug.eq.1) write(*,*) 's_input_mesh_udt_diff'
      call s_input_mesh_udt_diff
!
!     --------------------- 
!
      if (iflag_debug.eq.1) write(*,*) 'const_layers_4_dynamic'
      call const_layers_4_dynamic(ele_grp1, layer_tbl1)
!
!     --------------------- 
!
      if (iflag_debug.eq.1) write(*,*) 'const_mesh_informations'
      call const_mesh_informations(my_rank)
!
!     --------------------- 
!
      if (iflag_debug.eq.1) write(*,*) 'initialize_nod_field_data'
      call initialize_nod_field_data
!
!     --------------------- 
!
      call copy_num_processes_to_2nd
      call copy_comm_tbl_types(nod_comm, mesh_ref%nod_comm)
      call link_new_nod_geometry_type(node1, mesh_ref%node)
      call link_new_ele_connect_type(ele1, mesh_ref%ele)
      call link_new_overlaped_ele_type(ele1, mesh_ref%ele)
!
      call link_group_type(nod_grp1, group_ref%nod_grp)
      call link_group_type(ele_grp1, group_ref%ele_grp)
      call link_surf_group_type(sf_grp1, group_ref%surf_grp)
!
      call link_field_name_type(nod_fld1, phys_ref)
      call alloc_phys_data_type(mesh_ref%node%numnod, phys_ref)
      call allocate_vec_transfer(node1%numnod)
!
!     ---------------------
!
      if (iflag_debug.eq.1) write(*,*) 'allocate_vector_for_solver'
      call allocate_vector_for_solver(isix, node1%numnod)
      call allocate_2nd_iccg_matrix(isix, mesh_ref%node%numnod)
!
      call init_send_recv(nod_comm)
!
      if(iflag_debug.gt.0) write(*,*)' const_element_comm_tables_1st'
      call const_element_comm_tables_1st
!
!     --------------------- 
!
      if (iflag_debug.eq.1) write(*,*)  'cal_jacobian_element'
      call set_max_int_point_by_etype
      call cal_jacobian_element
!
      call dealloc_dxi_dx_type(jac1_3d_q)
      call dealloc_dxi_dx_type(jac1_3d_l)
!
!     --------------------- 
!
      if (iflag_debug.eq.1) write(*,*)  's_int_whole_volume_w_layer'
      call s_int_whole_volume_w_layer(ele1, jac1_3d_q, layer_tbl1)
!
      end subroutine initialize_udt_correlate
!
! ----------------------------------------------------------------------
!
      subroutine analyze_udt_correlate
!
      use m_nod_comm_table
      use m_geometry_data
      use m_node_phys_data
      use m_geometry_constants
      use m_layering_ele_list
      use m_node_phys_data
      use m_t_step_parameter
      use m_control_params_2nd_files
      use m_ctl_params_4_diff_udt
      use m_ucd_data
      use m_ucd_input_data
      use set_ucd_data_to_type
      use set_ucd_data
      use m_work_layer_correlate
      use ucd_IO_select
      use nod_phys_send_recv
!
      use fields_type_send_recv
      use correlation_all_layerd_data
!
      integer(kind=kint) :: istep, istep_ucd
!
!
      call link_fem_num_field_2_ucd_out
!
!     ---------------------
!
      ntot_correlate =  nod_fld1%ntot_phys
      ncomp_correlate = nod_fld1%ntot_phys
      nlayer_correlate = layer_tbl1%e_grp%num_grp
      call allocate_name_layer_correlate
      call allocate_all_layer_correlate
      call allocate_work_layer_correlate(layer_tbl1%e_grp%num_grp)
!
      call set_correlate_data_names
!
!     ---------------------
!
      do istep = i_step_init, i_step_number
        if ( mod(istep,i_step_output_ucd) .eq. izero) then
!
          istep_ucd = istep / i_step_output_ucd
!
          call set_data_by_read_ucd_once(my_rank, istep_ucd,            &
     &        ifmt_org_ucd, ref_udt_file_head)
!
          fem_ucd%ifmt_file = ifmt_org_ucd
          fem_ucd%file_prefix = tgt_udt_file_head
          call set_ucd_data_type_from_IO_once(my_rank, istep_ucd,       &
     &        mesh_ref%node%numnod, fem_ucd, phys_ref)
          fem_ucd%nnod = mesh_ref%node%numnod
!
          call nod_fields_send_recv(node1, nod_comm, nod_fld1)
          call phys_type_send_recv_all(mesh_ref, phys_ref)
!
!    output udt data
!
          call coord_transfer_4_1st_field
          call coord_transfer_4_2nd_field(mesh_ref%node%numnod,         &
     &        phys_ref)

!
          if (iflag_debug .gt. 0) write(*,*)                            &
     &          's_correlation_all_layerd_data'
          call s_correlation_all_layerd_data(mesh_ref%node%numnod,      &
     &        phys_ref)
!
          if (iflag_debug .gt. 0) write(*,*)                            &
     &          ' write_layerd_correlate_data', istep_ucd
          call write_layerd_correlate_data(my_rank, istep_ucd)
!
        end if
      end do
!
      end subroutine analyze_udt_correlate
!
! ----------------------------------------------------------------------
!
      end module analyzer_udt_correlation

