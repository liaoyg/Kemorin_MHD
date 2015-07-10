!link_group_type_2_1st_mesh.f90
!     module link_group_type_2_1st_mesh
!
!      Written by H. Matsui on Sep., 2006
!
!      subroutine link_node_group_to_type(nod_grp)
!      subroutine link_element_group_to_type(ele_grp)
!      subroutine link_surface_group_to_type(sf_grp)
!
!      subroutine link_1st_ele_grp_connect_type(tbls_ele_grp)
!
      module link_group_type_2_1st_mesh
!
      use m_precision
!
      implicit none
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine link_node_group_to_type(nod_grp)
!
      use t_group_data
      use m_node_group
!
      type(group_data), intent(inout) :: nod_grp
!
      nod_grp%num_grp = num_bc
      nod_grp%num_item = num_nod_bc
!
      nod_grp%grp_name =>  bc_name
      nod_grp%istack_grp => bc_istack
      nod_grp%item_grp =>   bc_item
!
      end subroutine link_node_group_to_type
!
!  ---------------------------------------------------------------------
!
      subroutine link_element_group_to_type(ele_grp)
!
      use t_group_data
      use m_element_group
!
      type(group_data), intent(inout) :: ele_grp
!
      ele_grp%num_grp =  num_mat
      ele_grp%num_item = num_mat_bc
!
      ele_grp%grp_name =>   mat_name
      ele_grp%istack_grp => mat_istack
      ele_grp%item_grp =>   mat_item
!
      end subroutine link_element_group_to_type
!
!  ---------------------------------------------------------------------
!
      subroutine link_surface_group_to_type(sf_grp)
!
      use t_group_data
      use m_surface_group
!
      type(surface_group_data), intent(inout) :: sf_grp
!
      sf_grp%num_grp = num_surf
      sf_grp%num_item = num_surf_bc
!
      sf_grp%grp_name =>   surf_name
      sf_grp%istack_grp => surf_istack
      sf_grp%item_sf_grp => surf_item
!
      end subroutine link_surface_group_to_type
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine link_1st_ele_grp_connect_type(tbls_ele_grp)
!
      use m_element_group_connect
      use t_group_connects
!
      type(element_group_table), intent(inout) :: tbls_ele_grp
!
!
      tbls_ele_grp%surf%ntot_e_grp = ele_grp_data1%surf%ntot_e_grp
      tbls_ele_grp%edge%ntot_e_grp = ele_grp_data1%edge%ntot_e_grp
      tbls_ele_grp%node%ntot_e_grp = ele_grp_data1%node%ntot_e_grp
!
      tbls_ele_grp%surf%nitem_e_grp =>  ele_grp_data1%surf%nitem_e_grp
      tbls_ele_grp%surf%istack_e_grp => ele_grp_data1%surf%istack_e_grp
      tbls_ele_grp%surf%item_e_grp =>   ele_grp_data1%surf%item_e_grp
!
      tbls_ele_grp%edge%nitem_e_grp =>  ele_grp_data1%edge%nitem_e_grp
      tbls_ele_grp%edge%istack_e_grp => ele_grp_data1%edge%istack_e_grp
      tbls_ele_grp%edge%item_e_grp =>   ele_grp_data1%edge%item_e_grp
!
      tbls_ele_grp%node%nitem_e_grp =>  ele_grp_data1%node%nitem_e_grp
      tbls_ele_grp%node%istack_e_grp => ele_grp_data1%node%istack_e_grp
      tbls_ele_grp%node%item_e_grp =>   ele_grp_data1%node%item_e_grp
!
      end subroutine link_1st_ele_grp_connect_type
!
!  ---------------------------------------------------------------------
!
      end module link_group_type_2_1st_mesh
