!cvt_quad_2_linear_mesh.f90
!      module cvt_quad_2_linear_mesh
!
!      Written by H. Matsui on Apr., 2006
!
!!      subroutine generate_linear_nod_by_quad                          &
!!     &         (node_q, ele_q, surf_q, mesh_l)
!!      subroutine connect_quad_mesh_2_linear                           &
!!     &         (node_q, ele_q, surf_q, mesh_l)
!!      subroutine connect_lag_mesh_2_linear(ele_q, mesh_l)
!!      subroutine gen_linear_group_info                                &
!!     &         (nod_grp, ele_grp, sf_grp, group_l)
!!
!!      subroutine set_internal_list_lin_20(node_q, ele_q, surf_q,      &
!!     &          mesh_l, surf_l, edge_l)
!!      subroutine set_internal_list_lin_27                             &
!!     &         (node_q, mesh_l, surf_l, edge_l)
!!
!!      subroutine init_linear_nod_phys(mesh_l, nod_fld_q, nod_fld_l)
!!      subroutine copy_nod_phys_2_linear                               &
!!     &         (node_q, nod_fld_q, mesh_l, nod_fld_l)
!!      subroutine generate_phys_on_surf(node_q, ele_q, surf_q,         &
!!     &          mesh_l, nod_fld_l)
!
      module cvt_quad_2_linear_mesh
!
      use m_precision
!
      implicit none
!
      integer(kind = kint), allocatable :: ie_4_333(:,:)
      private :: ie_4_333
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine generate_linear_nod_by_quad                            &
     &         (node_q, ele_q, surf_q, mesh_l)
!
      use t_mesh_data
      use t_geometry_data
      use t_surface_data
!
      use set_geometry_4_quad27
      use cal_mesh_position
!
      type(node_data), intent(in) ::    node_q
      type(element_data), intent(in) :: ele_q
      type(surface_data), intent(in) :: surf_q
!
      type(mesh_geometry), intent(inout) :: mesh_l
!
!
      mesh_l%node%numnod = node_q%numnod                                &
     &                    + surf_q%numsurf + ele_q%numele
!
      call allocate_node_geometry_type(mesh_l%node)
!
      call set_position_on_surf                                         &
     &   (node_q%numnod, surf_q%numsurf, ele_q%numele,                  &
     &    node_q%xx, ele_q%x_ele, surf_q%x_surf,                        &
     &    mesh_l%node%numnod,  mesh_l%node%xx)
!
      call set_spherical_position(mesh_l%node)
!
      end subroutine generate_linear_nod_by_quad
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine connect_quad_mesh_2_linear                             &
     &         (node_q, ele_q, surf_q, mesh_l)
!
      use m_geometry_constants
      use t_mesh_data
      use t_geometry_data
!
      use m_27quad_2_8x8linear
!
      type(node_data), intent(in) ::    node_q
      type(element_data), intent(in) :: ele_q
      type(surface_data), intent(in) :: surf_q
      type(mesh_geometry), intent(inout) :: mesh_l
!
!
      mesh_l%ele%numele = 8 * ele_q%numele
      mesh_l%ele%nnod_4_ele = num_t_linear
!
      call allocate_ele_connect_type(mesh_l%ele)
      call allocate_ele_geometry_type(mesh_l%ele)
!
      allocate(ie_4_333(ele_q%numele,27) )
!
      call gen_connect_quad27_from_quad20                               &
     &   (node_q%numnod, ele_q%numele, surf_q%numsurf, ele_q%ie,        &
     &    surf_q%isf_4_ele, ie_4_333)
!
      call set_27quad_2_8x8linear(ele_q%numele, ie_4_333,               &
     &    mesh_l%node%numnod, mesh_l%ele%ie)
!
      deallocate(ie_4_333)
!
      end subroutine connect_quad_mesh_2_linear
!
!  ---------------------------------------------------------------------
!
      subroutine connect_lag_mesh_2_linear(ele_q, mesh_l)
!
      use m_geometry_constants
      use t_mesh_data
      use t_geometry_data
!
      use m_27quad_2_8x8linear
!
      type(element_data), intent(in) :: ele_q
      type(mesh_geometry), intent(inout) :: mesh_l
!
!
      mesh_l%ele%numele = 8 * ele_q%numele
      mesh_l%ele%nnod_4_ele = num_t_linear
!
      call allocate_ele_connect_type(mesh_l%ele)
      call allocate_ele_geometry_type(mesh_l%ele)
!
      call set_27quad_2_8x8linear(ele_q%numele, ele_q%ie,               &
     &    mesh_l%node%numnod, mesh_l%ele%ie)
!
      end subroutine connect_lag_mesh_2_linear
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine gen_linear_group_info                                  &
     &         (nod_grp, ele_grp, sf_grp, group_l)
!
      use t_mesh_data
      use t_group_data
      use convert_group_2_linear
!
      type(group_data), intent(in) :: nod_grp
      type(group_data), intent(in) :: ele_grp
      type(surface_group_data), intent(in) :: sf_grp
      type(mesh_groups), intent(inout) :: group_l
!
!
      call link_group_type(nod_grp, group_l%nod_grp)
!
      group_l%ele_grp%num_grp = ele_grp%num_grp
      group_l%ele_grp%num_item = 8 * ele_grp%num_item
      call allocate_grp_type_num(group_l%ele_grp)
      call allocate_grp_type_item(group_l%ele_grp)
!
      call convert_ele_group_2_linear                                   &
     &   (ele_grp%num_grp, ele_grp%num_item, ele_grp%grp_name,          &
     &    ele_grp%istack_grp, ele_grp%item_grp,                         &
     &    group_l%ele_grp%num_item, group_l%ele_grp%grp_name,           &
     &    group_l%ele_grp%istack_grp, group_l%ele_grp%item_grp)
!
      group_l%surf_grp%num_grp = sf_grp%num_grp
      group_l%surf_grp%num_item = 4 * sf_grp%num_item
      call allocate_sf_grp_type_num(group_l%surf_grp)
      call allocate_sf_grp_type_item(group_l%surf_grp)
!
      call convert_surf_group_2_linear(sf_grp%num_grp, sf_grp%num_item, &
     &    sf_grp%grp_name, sf_grp%istack_grp, sf_grp%item_sf_grp,       &
     &    group_l%surf_grp%num_item, group_l%surf_grp%grp_name,         &
     &    group_l%surf_grp%istack_grp, group_l%surf_grp%item_sf_grp)
!
      end subroutine gen_linear_group_info
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine set_internal_list_lin_20(node_q, ele_q, surf_q,        &
     &          mesh_l, surf_l, edge_l)
!
      use m_machine_parameter
      use t_geometry_data
      use t_mesh_data
      use t_surface_data
      use t_edge_data
!
      use set_internal_list_4_linear
!
      type(node_data), intent(in) ::    node_q
      type(element_data), intent(in) :: ele_q
      type(surface_data), intent(in) :: surf_q
!
      type(mesh_geometry), intent(inout) :: mesh_l
      type(surface_data), intent(inout) :: surf_l
      type(edge_data),    intent(inout) :: edge_l
!
!
      call set_internal_list_4_linear_20                                &
     &   (node_q%numnod, node_q%internal_node, ele_q%numele,            &
     &    surf_q%numsurf, ele_q%interior_ele, surf_q%interior_surf,     &
     &    mesh_l%node%numnod, mesh_l%ele%numele, surf_l%numsurf,        &
     &    edge_l%numedge, mesh_l%ele%ie, surf_l%ie_surf,                &
     &    edge_l%ie_edge, mesh_l%ele%interior_ele,                      &
     &    surf_l%interior_surf, edge_l%interior_edge)
!
      end subroutine set_internal_list_lin_20
!
!  ---------------------------------------------------------------------
!
      subroutine set_internal_list_lin_27                               &
     &         (node_q, mesh_l, surf_l, edge_l)
!
      use m_machine_parameter
      use t_geometry_data
      use t_mesh_data
      use t_surface_data
      use t_edge_data
!
      use set_internal_list_4_linear
!
      type(node_data), intent(in) ::    node_q
!
      type(mesh_geometry), intent(inout) :: mesh_l
      type(surface_data), intent(inout) :: surf_l
      type(edge_data),    intent(inout) :: edge_l
!
!
      call set_internal_list_4_linear_27(node_q%internal_node,          &
     &    mesh_l%node%numnod, mesh_l%ele%numele, surf_l%numsurf,        &
     &    edge_l%numedge, mesh_l%ele%ie, surf_l%ie_surf,                &
     &    edge_l%ie_edge, mesh_l%ele%interior_ele,                      &
     &    surf_l%interior_surf, edge_l%interior_edge)
!
      end subroutine set_internal_list_lin_27
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine init_linear_nod_phys(mesh_l, nod_fld_q, nod_fld_l)
!
      use t_mesh_data
      use t_phys_data
!
      type(mesh_geometry), intent(in) :: mesh_l
      type(phys_data), intent(in) :: nod_fld_q
      type(phys_data), intent(inout) :: nod_fld_l
!
!
      call link_field_name_type(nod_fld_q, nod_fld_l)
      call alloc_phys_data_type(mesh_l%node%numnod, nod_fld_l)
!
      end subroutine init_linear_nod_phys
!
!  ---------------------------------------------------------------------
!
      subroutine copy_nod_phys_2_linear                                 &
     &         (node_q, nod_fld_q, mesh_l, nod_fld_l)
!
      use t_geometry_data
      use t_mesh_data
      use t_phys_data
!
      use set_data_4_quad27
!
      type(node_data), intent(in) ::    node_q
      type(phys_data), intent(in) :: nod_fld_q
!
      type(mesh_geometry), intent(in) :: mesh_l
      type(phys_data), intent(inout) :: nod_fld_l
!
!
      call copy_original_data(node_q%numnod, nod_fld_q%ntot_phys,       &
     &    nod_fld_q%d_fld, mesh_l%node%numnod, nod_fld_l%d_fld)
!
      end subroutine copy_nod_phys_2_linear
!
!  ---------------------------------------------------------------------
!
      subroutine generate_phys_on_surf(node_q, ele_q, surf_q,           &
     &          mesh_l, nod_fld_l)
!
      use t_geometry_data
      use t_mesh_data
      use t_phys_data
!
      use set_data_4_quad27
!
      type(node_data), intent(in) ::    node_q
      type(element_data), intent(in) :: ele_q
      type(surface_data), intent(in) :: surf_q
!
      type(mesh_geometry), intent(in) :: mesh_l
      type(phys_data), intent(inout) :: nod_fld_l
!
!
      call set_fields_on_surf                                           &
     &   (node_q%numnod, surf_q%numsurf, ele_q%numele,                  &
     &    ele_q%ie, surf_q%ie_surf, nod_fld_l%ntot_phys,                &
     &    mesh_l%node%numnod, nod_fld_l%d_fld)
!
      end subroutine generate_phys_on_surf
!
!  ---------------------------------------------------------------------
!
      end module cvt_quad_2_linear_mesh
