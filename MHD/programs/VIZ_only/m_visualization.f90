!>@file   m_visualization.f90
!!@brief  module m_visualization
!!
!!@author H. Matsui
!!@date Programmed in June, 2006
!
!>@brief Arrays for Field data IO for visualizers
!!
!!@verbatim
!!      subroutine mesh_setup_4_VIZ
!!      subroutine jacobian_4_VIZ
!!      subroutine element_normals_4_VIZ
!!@endverbatim
!
      module m_visualization
!
      use m_precision
      use m_machine_parameter
!
      use t_mesh_data
      use t_phys_data
      use t_ucd_data
      use t_next_node_ele_4_node
      use t_jacobian_3d
!
      implicit none
!
!>     Structure for mesh data
!>        (position, connectivity, group, and communication)
      type(mesh_data), save :: femmesh_VIZ
!
!>     Structure for element data (communication)
      type(element_comms), save :: elemesh_VIZ
!>     Structure for surface data
!>        (position, connectivity, and communication)
      type(surface_geometry), save :: surfmesh_VIZ
!>     Structure for edge data
!>        (position, connectivity, and communication)
      type(edge_geometry), save :: edgemesh_VIZ
!
!
!>       Structure for nodal field data
      type(phys_data), save :: field_VIZ
!
!
!>        Instance for FEM field data IO
     type(ucd_data), save :: ucd_VIZ
!>        Instance for numbers of FEM mesh for merged IO
!      type(merged_ucd_data), save :: m_ucd_SPH_TRNS
!
!>   Structure of included element list for each node
      type(element_around_node), save :: ele_4_nod_VIZ
!
!>     Stracture for Jacobians for linear element
      type(jacobians_3d), save :: jac_VIZ_l
!>     Stracture for Jacobians for quad element
      type(jacobians_3d), save :: jac_VIZ_q
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine mesh_setup_4_VIZ
!
      use calypso_mpi
      use m_t_step_parameter
      use m_array_for_send_recv
      use m_read_mesh_data
      use load_mesh_data
      use nod_phys_send_recv
      use const_mesh_information
      use set_parallel_file_name
      use const_element_comm_tables
      use set_ucd_data_to_type
      use ucd_IO_select
!
!   --------------------------------
!       setup mesh information
!   --------------------------------
!
!       load mesh informations
      if (iflag_debug.gt.0) write(*,*) 'input_mesh', mesh_file_head
      call input_mesh_data_type(my_rank, femmesh_VIZ,                   &
     &    surfmesh_VIZ%surf%nnod_4_surf, edgemesh_VIZ%edge%nnod_4_edge)
!
      if (iflag_debug.eq.1) write(*,*) 'const_mesh_infos'
      call s_const_mesh_types_info(my_rank, femmesh_VIZ,                &
     &    surfmesh_VIZ, edgemesh_VIZ)
!
      call allocate_vector_for_solver                                   &
     &   (isix, femmesh_VIZ%mesh%node%numnod)
!
      if(iflag_debug.gt.0) write(*,*)' init_send_recv'
      call init_send_recv(femmesh_VIZ%mesh%nod_comm)
!
      if(iflag_debug.gt.0) write(*,*)' const_element_comm_tbls'
      call const_ele_comm_tbl_global_id(femmesh_VIZ%mesh, elemesh_VIZ,  &
     &                                  surfmesh_VIZ, edgemesh_VIZ)
!
!     ---------------------
!
      ucd_VIZ%nnod =      femmesh_VIZ%mesh%node%numnod
      call sel_read_udt_param(my_rank, i_step_init, ucd_VIZ)
      call alloc_phys_data_type_by_output                               &
     &   (ucd_VIZ, femmesh_VIZ%mesh%node, field_VIZ)
!
      end subroutine mesh_setup_4_VIZ
!
! ----------------------------------------------------------------------
!
      subroutine jacobian_4_VIZ
!
      use m_fem_gauss_int_coefs
      use int_volume_of_domain
!
!
      if (iflag_debug.gt.0) write(*,*) 'const_jacobian_and_volume'
      call max_int_point_by_etype(femmesh_VIZ%mesh%ele%nnod_4_ele)
      call const_jacobian_and_volume(femmesh_VIZ%mesh%node,             &
     &    femmesh_VIZ%group%surf_grp, femmesh_VIZ%group%infty_grp,      &
     &    femmesh_VIZ%mesh%ele, jac_VIZ_l, jac_VIZ_q)
!
      end subroutine jacobian_4_VIZ
!
! ----------------------------------------------------------------------
!
      subroutine element_normals_4_VIZ
!
      use set_ele_id_4_node_type
      use set_normal_vectors
      use set_surf_grp_vectors
      use sum_normal_4_surf_group
!
!     --------------------- Connection information for PVR and fieldline
!     --------------------- init for fieldline and PVR
!
      if (iflag_debug.gt.0) write(*,*) 'set_ele_id_4_node'
      call set_ele_id_4_node                                            &
     &   (femmesh_VIZ%mesh%node, femmesh_VIZ%mesh%ele, ele_4_nod_VIZ)
!
      call jacobian_4_VIZ
!
!     --------------------- Surface jacobian for fieldline
!
      if (iflag_debug.eq.1) write(*,*)  'const_normal_vector'
        call const_normal_vector                                        &
     &     (femmesh_VIZ%mesh%node, surfmesh_VIZ%surf)
!
      if (iflag_debug.eq.1)  write(*,*) 'pick_normal_of_surf_group'
      call pick_normal_of_surf_group                                    &
     &   (surfmesh_VIZ%surf, femmesh_VIZ%group%surf_grp,                &
     &    femmesh_VIZ%group%tbls_surf_grp,                              &
     &    femmesh_VIZ%group%surf_grp_geom)
!
      if (iflag_debug.eq.1)  write(*,*) 's_sum_normal_4_surf_group'
      call s_sum_normal_4_surf_group(femmesh_VIZ%mesh%ele,              &
     &    femmesh_VIZ%group%surf_grp, femmesh_VIZ%group%surf_grp_geom)
!
      end subroutine element_normals_4_VIZ
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine set_field_data_4_VIZ(iflag, istep_ucd)
!
      use set_ucd_data_to_type
      use nod_phys_send_recv
!
      integer(kind = kint), intent(in) :: iflag, istep_ucd
!
!
      if(iflag .ne. 0) return
      call set_data_by_read_ucd(my_rank, istep_ucd, ucd_VIZ, field_VIZ)
!
      if (iflag_debug.gt.0)  write(*,*) 'phys_send_recv_all'
      call nod_fields_send_recv(femmesh_VIZ%mesh%node,                  &
     &    femmesh_VIZ%mesh%nod_comm, field_VIZ)
!
      end subroutine set_field_data_4_VIZ
!
! ----------------------------------------------------------------------
!
      end module m_visualization
