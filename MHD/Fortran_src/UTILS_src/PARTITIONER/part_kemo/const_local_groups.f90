!const_local_groups.f90
!      module const_local_groups
!
!      Written by H. Matsui on Aug., 2007
!
!      subroutine s_const_local_groups(orggroup, newgroup)
!
      module const_local_groups
!
      use m_precision
!
      use set_group_4_subdomain
!
      implicit none
!
!
      private :: const_local_nod_group, const_local_ele_group
      private :: const_local_surf_group
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine s_const_local_groups(orggroup, newgroup)
!
      use t_group_data
      use t_mesh_data
!
      type(mesh_groups), intent(in) :: orggroup
      type(mesh_groups), intent(inout) :: newgroup
!
!
      call const_local_nod_group(orggroup%nod_grp, newgroup%nod_grp)
      call const_local_ele_group(orggroup%ele_grp, newgroup%ele_grp)
      call const_local_surf_group(orggroup%surf_grp, newgroup%surf_grp)
!
      end subroutine s_const_local_groups
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine const_local_nod_group(nod_grp, new_nod_grp)
!
      use t_group_data
!
      type(group_data), intent(in) :: nod_grp
      type(group_data), intent(inout) :: new_nod_grp
!
!
      new_nod_grp%num_grp = nod_grp%num_grp
      call allocate_grp_type_num(new_nod_grp)
      call count_local_node_group(nod_grp, new_nod_grp)
!
      call allocate_grp_type_item(new_nod_grp)
      call set_local_node_group(nod_grp, new_nod_grp)
!
      end subroutine const_local_nod_group
!
!  ---------------------------------------------------------------------
!
      subroutine const_local_ele_group(ele_grp, new_ele_grp)
!
      use t_group_data
!
      type(group_data), intent(in) :: ele_grp
      type(group_data), intent(inout) :: new_ele_grp
!
!
      new_ele_grp%num_grp = ele_grp%num_grp
      call allocate_grp_type_num(new_ele_grp)
      call count_local_ele_group(ele_grp, new_ele_grp)
!
      call allocate_grp_type_item(new_ele_grp)
      call set_local_ele_group(ele_grp, new_ele_grp)
!
      end subroutine const_local_ele_group
!
!  ---------------------------------------------------------------------
!
      subroutine const_local_surf_group(sf_grp, new_sf_grp)
!
      use t_group_data
!
      type(surface_group_data), intent(in) :: sf_grp
      type(surface_group_data), intent(inout) :: new_sf_grp
!
!
      new_sf_grp%num_grp = sf_grp%num_grp
      call allocate_sf_grp_type_num(new_sf_grp)
      call count_local_surf_group(sf_grp, new_sf_grp)
!
      call allocate_sf_grp_type_item(new_sf_grp)
      call set_local_surf_group(sf_grp, new_sf_grp)
!
      end subroutine const_local_surf_group
!
!  ---------------------------------------------------------------------
!
      end module const_local_groups
