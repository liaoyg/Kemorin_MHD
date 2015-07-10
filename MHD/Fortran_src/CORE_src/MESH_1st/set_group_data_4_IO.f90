!set_group_data_4_IO.f90
!     module set_group_data_4_IO
!
!      written by H. Matsui on Dec., 2006
!
!      subroutine copy_group_data_from_IO
!      subroutine copy_group_data_to_IO
!
      module set_group_data_4_IO
!
      use m_precision
!
      implicit  none
!
      private :: copy_node_group_from_IO, set_node_group_to_IO
      private :: copy_element_group_from_IO, set_element_group_to_IO
      private :: copy_surface_group_from_IO, set_surface_group_to_IO
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine copy_group_data_from_IO
!
!
      call copy_node_group_from_IO
      call copy_element_group_from_IO
      call copy_surface_group_from_IO
!
      end subroutine copy_group_data_from_IO
!
!-----------------------------------------------------------------------
!
      subroutine copy_group_data_to_IO
!
!
      call set_node_group_to_IO
      call set_element_group_to_IO
      call set_surface_group_to_IO
!
      end subroutine copy_group_data_to_IO
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine copy_node_group_from_IO
!
      use m_read_boundary_data
      use m_node_group
!
!   set node group
!
      num_bc = num_bc_dummy
      if (num_bc/=0) then
!
        num_nod_bc = num_nod_bc_dummy
        call allocate_boundary_data
!
        bc_name(1:num_bc) =     bc_name_dummy(1:num_bc)
        bc_istack(0:num_bc) =   bc_istack_dummy(0:num_bc)
        bc_item(1:num_nod_bc) = bc_item_dummy(1:num_nod_bc)
      end if
      call deallocate_bc_item_dummy
!
      end subroutine copy_node_group_from_IO
!
!-----------------------------------------------------------------------
!
      subroutine copy_element_group_from_IO
!
      use m_read_boundary_data
      use m_element_group
!
!    set element group
!
      num_mat =     num_mat_dummy
      if (num_mat/=0) then
!
        num_mat_bc = num_mat_bc_dummy
        call allocate_material_data
!
        mat_name(1:num_mat) =    mat_name_dummy(1:num_mat)
        mat_istack(0:num_mat) =  mat_istack_dummy(0:num_mat)
        mat_item(1:num_mat_bc) = mat_item_dummy(1:num_mat_bc)
      end if
      call deallocate_bc_ele_item_dummy
!
      end subroutine copy_element_group_from_IO
!
!-----------------------------------------------------------------------
!
      subroutine copy_surface_group_from_IO
!
      use m_read_boundary_data
      use m_surface_group
!
!   set surface group
!
      sf_grp1%num_grp = num_surf_dummy
      if (sf_grp1%num_grp .gt. 0) then
!
        sf_grp1%num_item = num_surf_bc_dummy
        call allocate_surface_data
!
        surf_name(1:sf_grp1%num_grp)                                    &
     &     = surf_name_dummy(1:sf_grp1%num_grp)
        surf_istack(0:sf_grp1%num_grp)                                  &
     &     =  surf_istack_dummy(0:sf_grp1%num_grp)
        surf_item(1,1:sf_grp1%num_item)                                 &
     &     = surf_item_dummy(1:sf_grp1%num_item,1)
        surf_item(2,1:sf_grp1%num_item)                                 &
     &     = surf_item_dummy(1:sf_grp1%num_item,2)
      end if
      call deallocate_bc_sf_item_dummy
!
      end subroutine copy_surface_group_from_IO
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine set_node_group_to_IO
!
      use m_node_group
      use m_read_boundary_data
!
!
      num_bc_dummy = num_bc
      num_nod_bc_dummy = num_nod_bc
      call allocate_bc_stack_dummy
!
      bc_name_dummy(1:num_bc) = bc_name(1:num_bc)
      bc_istack_dummy(0:num_bc) = bc_istack(0:num_bc)
!
      call allocate_bc_item_dummy
      bc_item_dummy(1:num_nod_bc) = bc_item(1:num_nod_bc)
!
      call deallocate_boundary_data
!
      end subroutine set_node_group_to_IO
!
!-----------------------------------------------------------------------
!
      subroutine set_element_group_to_IO
!
      use m_element_group
      use m_read_boundary_data
!
!
      num_mat_dummy = num_mat
      num_mat_bc_dummy = num_mat_bc
      call allocate_bc_ele_stack_dummy
!
      mat_name_dummy(1:num_mat) = mat_name(1:num_mat)
      mat_istack_dummy(0:num_mat) = mat_istack(0:num_mat)
!
      call allocate_bc_ele_item_dummy
      mat_item_dummy(1:num_mat_bc) = mat_item(1:num_mat_bc)
!
      call deallocate_material_data
!
      end subroutine set_element_group_to_IO
!
!-----------------------------------------------------------------------
!
      subroutine set_surface_group_to_IO
!
      use m_surface_group
      use m_read_boundary_data
!
!
      num_surf_dummy = sf_grp1%num_grp
      num_surf_bc_dummy = sf_grp1%num_item
      call allocate_bc_sf_stack_dummy
!
      surf_name_dummy(1:sf_grp1%num_grp) = surf_name(1:sf_grp1%num_grp)
      surf_istack_dummy(0:sf_grp1%num_grp)                              &
     &     = surf_istack(0:sf_grp1%num_grp)
!
      call allocate_bc_sf_item_dummy
      surf_item_dummy(1:sf_grp1%num_item,1)                             &
     &     = surf_item(1,1:sf_grp1%num_item)
      surf_item_dummy(1:sf_grp1%num_item,2)                             &
     &     = surf_item(2,1:sf_grp1%num_item)
!
      call deallocate_surface_data
!
      end subroutine set_surface_group_to_IO
!
!-----------------------------------------------------------------------
!
      end module set_group_data_4_IO
