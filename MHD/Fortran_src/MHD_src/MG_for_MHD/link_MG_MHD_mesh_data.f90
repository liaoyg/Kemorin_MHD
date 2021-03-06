!link_MG_MHD_mesh_data.f90
!     module link_MG_MHD_mesh_data
!
!        programmed H.Matsui on Dec., 2008
!
!!      subroutine s_link_MG_MHD_mesh_data                              &
!!     &        (MGCG_WK, MG_mesh, MG_MHD_mesh, ele_1st, MHD_matrices)
!!        type(MGCG_data), intent(in) :: MGCG_WK
!!        type(mesh_data), intent(in), target                           &
!!     &                    :: MG_mesh(MGCG_WK%num_MG_level)
!!        type(mesh_data_MHD), intent(in), target                       &
!!     &                    :: MG_MHD_mesh(MGCG_WK%num_MG_level)
!!        type(element_data), intent(in) :: ele_1st
!!        type(MHD_MG_matrices), intent(inout) :: MHD_matrices
!
      module link_MG_MHD_mesh_data
!
      use m_precision
!
!
      implicit none
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine s_link_MG_MHD_mesh_data                                &
     &         (MGCG_WK, MG_mesh, MG_MHD_mesh, ele_1st, MHD_matrices)
!
      use t_MGCG_data
      use t_mesh_data
      use t_geometry_data
      use t_solver_djds_MHD
      use t_interpolate_table
      use interpolate_by_module
!
      type(MGCG_data), intent(in) :: MGCG_WK
      type(mesh_data), intent(in), target                               &
     &                    :: MG_mesh(MGCG_WK%num_MG_level)
      type(mesh_data_MHD), intent(in), target                           &
     &                    :: MG_MHD_mesh(MGCG_WK%num_MG_level)
      type(element_data), intent(in) :: ele_1st
      type(MHD_MG_matrices), intent(inout) :: MHD_matrices
      integer(kind = kint) :: i_level
!
!
!
      call init_interpolate_mat_type                                    &
     &    (ele_1st, MHD_matrices%MG_interpolate(1)%f2c%tbl_org,         &
     &     MHD_matrices%MG_interpolate(1)%f2c%mat)
      call init_interpolate_mat_type(MG_mesh(1)%mesh%ele,               &
     &    MHD_matrices%MG_interpolate(1)%c2f%tbl_org,                   &
     &    MHD_matrices%MG_interpolate(1)%c2f%mat)
!
      do i_level = 2, MGCG_WK%num_MG_level
        call link_comm_tbl_types(MG_MHD_mesh(i_level)%nod_fl_comm,      &
     &      MHD_matrices%MG_comm_fluid(i_level))
!
        call init_interpolate_mat_type(MG_mesh(i_level-1)%mesh%ele,     &
     &      MHD_matrices%MG_interpolate(i_level)%f2c%tbl_org,           &
     &      MHD_matrices%MG_interpolate(i_level)%f2c%mat)
        call init_interpolate_mat_type(MG_mesh(i_level)%mesh%ele,       &
     &      MHD_matrices%MG_interpolate(i_level)%c2f%tbl_org,           &
     &      MHD_matrices%MG_interpolate(i_level)%c2f%mat)
      end do
!
      end subroutine s_link_MG_MHD_mesh_data
!
!  ---------------------------------------------------------------------
!
      end module link_MG_MHD_mesh_data
