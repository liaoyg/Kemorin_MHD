!m_element_group_connect.f90
!     module m_element_group_connect
!
!> @brief connectivity data for element group items
!
!     Writteg by H.Matsui on Aug., 2006
!
!      subroutine const_ele_group_connect
!      subroutine deallocate_ele_group_connect
!
      module m_element_group_connect
!
      use m_precision
      use t_group_connects
!
      implicit  none
!
!>   Structure of connectivities for element group
      type(element_group_table), save :: ele_grp_data1
!
!>   Array for marking
      integer(kind=kint), allocatable :: imark_4_grp(:)
!
      private :: imark_4_grp
      private :: set_surf_4_ele_group, set_edge_4_ele_group
      private :: set_node_4_ele_group
      private :: allocate_imark_4_grp, deallocate_imark_4_grp
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine const_ele_group_connect
!
      use m_machine_parameter
!
!
       if (iflag_debug.eq.1) write(*,*) 'set_surf_4_ele_group'
      call set_surf_4_ele_group
!
       if (iflag_debug.eq.1) write(*,*) 'set_edge_4_ele_group'
      call set_edge_4_ele_group
!
       if (iflag_debug.eq.1) write(*,*) 'set_node_4_ele_group'
      call set_node_4_ele_group
!
      end subroutine const_ele_group_connect
!
! ----------------------------------------------------------------------
!
      subroutine deallocate_ele_group_connect
!
!
      call dealloc_grp_connect(ele_grp_data1%surf)
      call dealloc_grp_connect(ele_grp_data1%edge)
      call dealloc_grp_connect(ele_grp_data1%node)
!
      end subroutine deallocate_ele_group_connect
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine set_surf_4_ele_group
!
      use m_geometry_constants
      use m_geometry_parameter
      use m_geometry_data
      use m_element_group
!
      use set_node_4_group
!
!
      call allocate_imark_4_grp(numsurf)
      call alloc_num_other_grp(ele_grp1%num_grp, ele_grp_data1%surf)
!
      call count_nod_4_ele_grp(numsurf, numele, nsurf_4_ele, isf_4_ele, &
     &   ele_grp1%num_grp, ele_grp1%num_item, ele_grp1%istack_grp,      &
     &   ele_grp1%item_grp, ele_grp_data1%surf%ntot_e_grp,              &
     &   ele_grp_data1%surf%nitem_e_grp,                                &
     &   ele_grp_data1%surf%istack_e_grp, imark_4_grp)
!
!
      call alloc_item_other_grp(ele_grp_data1%surf)
!
      call set_nod_4_ele_grp(numsurf, numele, nsurf_4_ele, isf_4_ele,   &
     &    ele_grp1%num_grp, ele_grp1%num_item, ele_grp1%istack_grp,     &
     &    ele_grp1%item_grp, ele_grp_data1%surf%ntot_e_grp,             &
     &    ele_grp_data1%surf%nitem_e_grp,                               &
     &    ele_grp_data1%surf%istack_e_grp,                              &
     &    ele_grp_data1%surf%item_e_grp, imark_4_grp)
!
      call deallocate_imark_4_grp
!
      end subroutine set_surf_4_ele_group
!
!-----------------------------------------------------------------------
!
      subroutine set_edge_4_ele_group
!
      use m_geometry_constants
      use m_geometry_parameter
      use m_geometry_data
      use m_element_group
!
      use set_node_4_group
!
!
      call allocate_imark_4_grp(numedge)
      call alloc_num_other_grp(ele_grp1%num_grp, ele_grp_data1%edge)
!
      call count_nod_4_ele_grp                                          &
     &   (numedge, numele, nedge_4_ele, iedge_4_ele,                    &
     &    ele_grp1%num_grp, ele_grp1%num_item, ele_grp1%istack_grp,     &
     &    ele_grp1%item_grp, ele_grp_data1%edge%ntot_e_grp,             &
     &    ele_grp_data1%edge%nitem_e_grp,                               &
     &    ele_grp_data1%edge%istack_e_grp, imark_4_grp)
!
!
      call alloc_item_other_grp(ele_grp_data1%edge)
!
      call set_nod_4_ele_grp(numedge, numele, nedge_4_ele, iedge_4_ele, &
     &    ele_grp1%num_grp, ele_grp1%num_item, ele_grp1%istack_grp,     &
     &    ele_grp1%item_grp, ele_grp_data1%edge%ntot_e_grp,             &
     &    ele_grp_data1%edge%nitem_e_grp,                               &
     &    ele_grp_data1%edge%istack_e_grp,                              &
     &    ele_grp_data1%edge%item_e_grp, imark_4_grp)
!
      call deallocate_imark_4_grp
!
      end subroutine set_edge_4_ele_group
!
!-----------------------------------------------------------------------
!
      subroutine set_node_4_ele_group
!
      use m_geometry_parameter
      use m_geometry_data
      use m_element_group
!
      use set_node_4_group
!
!
      call allocate_imark_4_grp(numnod)
      call alloc_num_other_grp(ele_grp1%num_grp, ele_grp_data1%node)
!
      call count_nod_4_ele_grp(numnod, numele, nnod_4_ele, ie,          &
     &   ele_grp1%num_grp, ele_grp1%num_item, ele_grp1%istack_grp,      &
     &   ele_grp1%item_grp, ele_grp_data1%node%ntot_e_grp,              &
     &   ele_grp_data1%node%nitem_e_grp,                                &
     &   ele_grp_data1%node%istack_e_grp, imark_4_grp)
!
!
      call alloc_item_other_grp(ele_grp_data1%node)
!
      call set_nod_4_ele_grp(numnod, numele, nnod_4_ele, ie,            &
     &    ele_grp1%num_grp, ele_grp1%num_item, ele_grp1%istack_grp,     &
     &    ele_grp1%item_grp, ele_grp_data1%node%ntot_e_grp,             &
     &    ele_grp_data1%node%nitem_e_grp,                               &
     &    ele_grp_data1%node%istack_e_grp,                              &
     &    ele_grp_data1%node%item_e_grp, imark_4_grp)
!
      call deallocate_imark_4_grp
!
      end subroutine set_node_4_ele_group
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine allocate_imark_4_grp(num)
!
      integer(kind = kint), intent(in) :: num
!
      allocate(imark_4_grp(num))
      imark_4_grp = 0
!
      end subroutine allocate_imark_4_grp
!
!-----------------------------------------------------------------------
!
      subroutine deallocate_imark_4_grp
!
      deallocate(imark_4_grp)
!
      end subroutine deallocate_imark_4_grp
!
!-----------------------------------------------------------------------
!
      end module m_element_group_connect
