!link_MG_MHD_mesh_data.f90
!     module link_MG_MHD_mesh_data
!
!        programmed H.Matsui on Dec., 2008
!
!      subroutine s_link_MG_MHD_mesh_data
!
      module link_MG_MHD_mesh_data
!
      use m_precision
!
      use t_comm_table
!
      implicit none
!
      private :: link_first_comm_to_type
      private :: link_first_ele_connect_type
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine s_link_MG_MHD_mesh_data
!
      use m_solver_djds_MHD
      use m_type_AMG_data_4_MHD
      use m_type_AMG_mesh
      use m_type_AMG_data
      use t_interpolate_table
      use t_geometry_data
      use interpolate_by_type
!
      integer(kind = kint) :: i_level
      type(element_data) :: original_ele
!
!
      call link_first_comm_to_type(MG_comm(0))
      call link_comm_tbl_types(DJDS_comm_fl, MG_comm_fl(0))
!
      call link_comm_tbl_types(MG_mesh(1)%mesh%nod_comm, MG_comm(1) )
      call link_comm_tbl_types(MG_MHD_mesh(1)%nod_fl_comm,              &
     &    MG_comm_fl(1))
!
      call link_first_ele_connect_type(original_ele)
      call init_interpolate_mat_type                                    &
     &   (original_ele, MG_itp(1)%f2c)
      call init_interpolate_mat_type                                    &
     &   (MG_mesh(1)%mesh%ele, MG_itp(1)%c2f)
!
      do i_level = 2, num_MG_level
        call link_comm_tbl_types(MG_mesh(i_level)%mesh%nod_comm,        &
     &      MG_comm(i_level) )
        call link_comm_tbl_types(MG_MHD_mesh(i_level)%nod_fl_comm,      &
     &      MG_comm_fl(i_level))
!
        call init_interpolate_mat_type                                  &
     &     (MG_mesh(i_level-1)%mesh%ele, MG_itp(i_level)%f2c)
        call init_interpolate_mat_type                                  &
     &     (MG_mesh(i_level)%mesh%ele, MG_itp(i_level)%c2f)
      end do
!
!
      end subroutine s_link_MG_MHD_mesh_data
!
!  ---------------------------------------------------------------------
!
      subroutine link_first_comm_to_type(MG_comm_0)
!
      use m_nod_comm_table
!
      type(communication_table), intent(inout) :: MG_comm_0
!
!
      MG_comm_0%num_neib =    nod_comm%num_neib
      MG_comm_0%ntot_import = nod_comm%ntot_import
      MG_comm_0%ntot_export = nod_comm%ntot_export
!
      MG_comm_0%id_neib =>       nod_comm%id_neib
      MG_comm_0%num_import =>    nod_comm%num_import
      MG_comm_0%istack_import => nod_comm%istack_import
      MG_comm_0%item_import =>   nod_comm%item_import
      MG_comm_0%num_export =>    nod_comm%num_export
      MG_comm_0%istack_export => nod_comm%istack_export
      MG_comm_0%item_export =>   nod_comm%item_export
!
      end subroutine link_first_comm_to_type
!
!  ---------------------------------------------------------------------
!
      subroutine link_first_ele_connect_type(ele)
!
      use m_geometry_data
      use t_geometry_data
!
      type(element_data), intent(inout) :: ele
!
!
      ele%numele =     ele1%numele
      ele%nnod_4_ele = ele1%nnod_4_ele
!
      ele%iele_global => ele1%iele_global
      ele%elmtyp =>      ele1%elmtyp
      ele%nodelm =>      nodelm
      ele%ie =>          ele1%ie
!
      end subroutine link_first_ele_connect_type
!
!-----------------------------------------------------------------------
!
      end module link_MG_MHD_mesh_data
