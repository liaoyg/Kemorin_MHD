!FEM_analyzer_sph_trans.f90
!
!      subroutine FEM_initialize_sph_trans
!      subroutine FEM_analyze_sph_trans(i_step, visval, ucd)
!
!      subroutine SPH_to_FEM_bridge_sph_trans
!      subroutine FEM_finalize_sph_trans(ucd, m_ucd)
!        type(ucd_data), intent(inout) :: ucd
!        type(merged_ucd_data), intent(inout) :: m_ucd
!
      module FEM_analyzer_sph_trans
!
      use m_precision
      use m_constants
      use m_machine_parameter
      use calypso_mpi
!
      use t_ucd_data
!
      implicit none
!
!>        Instance for FEM field data IO
      type(ucd_data), save :: input_ucd
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine FEM_initialize_sph_trans
!
      use m_nod_comm_table
      use m_geometry_data
      use m_control_params_2nd_files
      use m_array_for_send_recv
      use m_t_step_parameter
      use m_node_phys_data
      use m_ele_sf_eg_comm_tables
!
      use nod_phys_send_recv
      use const_mesh_info
      use int_volume_of_domain
      use set_normal_vectors
      use set_surf_grp_vectors
      use sum_normal_4_surf_group
      use output_parallel_ucd_file
      use ucd_IO_select
!
      use copy_all_field_4_sph_trans
!
!  -----    construct geometry informations
!
      if (iflag_debug.gt.0) write(*,*) 'allocate_vector_for_solver'
      call allocate_vector_for_solver(isix, node1%numnod)
!
      if(iflag_debug.gt.0) write(*,*)' init_send_recv'
      call init_send_recv(nod_comm)
!
      if (iflag_debug.gt.0) write(*,*) 'const_mesh_informations'
      call const_mesh_informations(my_rank)
!
      if(iflag_debug.gt.0) write(*,*)' const_element_comm_tables_1st'
      call const_element_comm_tables_1st
!
      call deallocate_edge_geom_type(edge1)
!
!  -------------------------------
!
      if (iflag_debug.gt.0) write(*,*) 'alloc_phys_data_type'
      call alloc_phys_data_type(node1%numnod, nod_fld1)
!
!  -------------------------------
!
      call calypso_MPI_barrier
      call set_ucd_file_format(ifmt_org_ucd, input_ucd)
      call set_ucd_file_prefix(org_ucd_header, input_ucd)
!
      input_ucd%nnod = ione
      call sel_read_udt_param(my_rank, i_step_init, input_ucd)
!
      end subroutine FEM_initialize_sph_trans
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine FEM_analyze_sph_trans(i_step, visval)
!
      use m_control_params_2nd_files
      use m_t_step_parameter
      use m_node_phys_data
      use set_ucd_data_to_type
!
      integer (kind =kint), intent(in) :: i_step
      integer (kind =kint), intent(inout) :: visval
!
!
!*  ----------   Count steps for visualization
!*
      visval =  mod(i_step,i_step_output_ucd)
!
!*  -----------  Output volume data --------------
!*
      if(visval .eq. 0) then
        call set_ucd_file_prefix(org_ucd_header, input_ucd)
        call set_data_by_read_ucd(my_rank, i_step, input_ucd, nod_fld1)
      end if
!
      end subroutine FEM_analyze_sph_trans
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine SPH_to_FEM_bridge_sph_trans(fld_IO)
!
      use t_field_data_IO
      use copy_rj_phys_data_4_IO
!
      type(field_IO), intent(inout) :: fld_IO
!
!
      if (iflag_debug.gt.0) write(*,*) 'copy_rj_all_phys_name_to_IO'
      call copy_rj_all_phys_name_to_IO(fld_IO)
      call alloc_phys_data_IO(fld_IO)
      call alloc_merged_field_stack(nprocs, fld_IO)
!
      end subroutine SPH_to_FEM_bridge_sph_trans
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine FEM_finalize_sph_trans(ucd, m_ucd)
!
      use m_t_step_parameter
      use t_ucd_data
      use output_parallel_ucd_file
!
      type(ucd_data), intent(inout) :: ucd
      type(merged_ucd_data), intent(inout) :: m_ucd
!
!
      if(i_step_output_ucd .gt. 0) then
        call finalize_ucd_file_output(ucd, m_ucd)
      end if
!
      end subroutine FEM_finalize_sph_trans
!
!-----------------------------------------------------------------------
!
      end module FEM_analyzer_sph_trans
