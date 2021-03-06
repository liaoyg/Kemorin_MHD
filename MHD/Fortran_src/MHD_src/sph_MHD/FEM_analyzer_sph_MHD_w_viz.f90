!>@file   FEM_analyzer_sph_MHD_w_viz.f90
!!@brief  module FEM_analyzer_sph_MHD_w_viz
!!
!!@author H. Matsui
!!@date Programmed in Nov., 2014
!
!>@brief Top routines to transfer spherical harmonics grids data
!!       to FEM data for data visualization
!!
!!@verbatim
!!      subroutine FEM_initialize_w_viz                                 &
!!     &         (MHD_files, MHD_step, geofem, ele_mesh,                &
!!     &          iphys, nod_fld, next_tbl, jacobians, MHD_IO)
!!        type(MHD_file_IO_params), intent(in) :: MHD_files
!!        type(MHD_step_param), intent(in) :: MHD_step
!!        type(mesh_data), intent(inout) :: geofem
!!        type(element_geometry), intent(inout) :: ele_mesh
!!        type(phys_address), intent(in) :: iphys
!!        type(phys_data), intent(inout) :: nod_fld
!!        type(next_nod_ele_table), intent(inout) :: next_tbl
!!        type(jacobians_type), intent(inout) :: jacobians
!!        type(MHD_IO_data), intent(inout) :: MHD_IO
!!@endverbatim
!!
!!@n @param  i_step       Current time step
!!@n @param  visval       Return flag to call visualization routines
!
      module FEM_analyzer_sph_MHD_w_viz
!
      use m_precision
      use m_constants
!
      use calypso_mpi
      use m_MHD_step_parameter
      use m_machine_parameter
      use m_work_time
!
      use t_mesh_data
      use t_phys_data
      use t_phys_address
      use t_next_node_ele_4_node
      use t_shape_functions
      use t_jacobians
      use t_VIZ_step_parameter
      use t_MHD_step_parameter
      use t_MHD_file_parameter
      use t_cal_max_indices
      use t_MHD_IO_data
!
      implicit none
!
      type(shape_finctions_at_points), save, private :: spfs_M
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine FEM_initialize_w_viz                                   &
     &         (MHD_files, MHD_step, geofem, ele_mesh,                  &
     &          iphys, nod_fld, next_tbl, jacobians, MHD_IO)
!
      use set_table_4_RHS_assemble
      use FEM_analyzer_sph_MHD
      use int_volume_of_domain
      use set_normal_vectors
!
      type(MHD_file_IO_params), intent(in) :: MHD_files
      type(MHD_step_param), intent(in) :: MHD_step
      type(mesh_data), intent(inout) :: geofem
      type(element_geometry), intent(inout) :: ele_mesh
      type(phys_address), intent(inout) :: iphys
      type(phys_data), intent(inout) :: nod_fld
      type(next_nod_ele_table), intent(inout) :: next_tbl
      type(jacobians_type), intent(inout) :: jacobians
      type(MHD_IO_data), intent(inout) :: MHD_IO
!
!   --------------------------------
!       setup mesh information
!   --------------------------------
!
!  --  init FEM mesh data
!
      if(iflag_debug .gt. 0) write(*,*) 'FEM_initialize_sph_MHD'
      call FEM_initialize_sph_MHD(MHD_files, MHD_step,                  &
     &    geofem, ele_mesh, iphys, nod_fld, MHD_IO)
!
!  -------------------------------
!
      if(MHD_step%viz_step%FLINE_t%increment .gt. 0) then
        if (iflag_debug.gt.0) write(*,*) 'set_element_on_node_in_mesh'
        call set_element_on_node_in_mesh                                &
     &     (geofem%mesh, next_tbl%neib_ele)
      end if
!
      if(MHD_step%viz_step%PVR_t%increment .le. 0) Return
!
!  -----  If there is no volume rendering... return
!
      if (iflag_debug.eq.1) write(*,*)  'set_max_integration_points'
      allocate(jacobians%g_FEM)
      call set_max_integration_points(ione, jacobians%g_FEM)
      call const_jacobian_volume_normals(my_rank, nprocs,               &
     &    geofem%mesh, ele_mesh%surf, geofem%group, spfs_M, jacobians)
!
      end subroutine FEM_initialize_w_viz
!
!-----------------------------------------------------------------------
!
      end module FEM_analyzer_sph_MHD_w_viz
