!analyzer_sph_bc_temp.f90
!
!      module analyzer_sph_bc_temp
!
!      modified by H. Matsui on Aug., 2006 
!
!      subroutine initilize_bc_temp
!      subroutine analyze_bc_temp
!
!..................................................
!
      module analyzer_sph_bc_temp
!
      use m_precision
      use m_machine_parameter
!
      use calypso_mpi
      use t_mesh_data
!
      implicit none
!
      type(mesh_data), save :: femmesh
      type(element_comms), save ::    ele_mesh
      type(surface_geometry), save :: surf_mesh
      type(edge_geometry), save ::    edge_mesh
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine initilize_bc_temp
!
      use m_ctl_data_test_bc_temp
      use m_ctl_params_test_bc_temp
      use load_mesh_data
      use const_mesh_information
!
!
!     ----- read control data
!
      if (iflag_debug.gt.0) write(*,*) 'read_control_4_bc_temp'
      call read_control_4_bc_temp
!
      if (iflag_debug.gt.0) write(*,*) 'set_ctl_params_sph_bc_temp'
      call set_ctl_params_sph_bc_temp
!
!  --  read geometry
!
      if (iflag_debug.gt.0) write(*,*) 'input_mesh'
      call input_mesh_data_type(my_rank, femmesh,                       &
     &    surf_mesh%surf%nnod_4_surf, edge_mesh%edge%nnod_4_edge)
!
      if (iflag_debug.eq.1) write(*,*) 'const_mesh_infos'
      call s_const_mesh_types_info                                      &
     &   (my_rank, femmesh, surf_mesh, edge_mesh)
!
       end subroutine initilize_bc_temp
!
! ----------------------------------------------------------------------
!
      subroutine analyze_bc_temp
!
      use const_sph_boundary_temp
!
      if (iflag_debug.gt.0) write(*,*) 'const_sph_temp_bc'
      call const_sph_temp_bc                                            &
     &   (femmesh%mesh%node, femmesh%group%nod_grp)
!
      end subroutine analyze_bc_temp
!
! ----------------------------------------------------------------------
!
      end module analyzer_sph_bc_temp
