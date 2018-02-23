!
!      module set_merged_geometry
!
!      Written by H. Matsui
!
!!      subroutine set_merged_mesh_and_group(mesh_file)
!!      subroutine set_merged_node_and_element(mesh_file)
!!      subroutine set_overlapped_mesh_and_group(mesh_file, nnod_4_ele)
!!        type(field_IO_params), intent(in) :: mesh_file
!
      module set_merged_geometry
!
      use m_precision
!
      use t_file_IO_parameter
      use m_geometry_data_4_merge
      use set_geometry_to_merge
      use count_number_with_overlap
!
      implicit none
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine set_merged_mesh_and_group(mesh_file)
!
      use const_merged_groups
!
      type(field_IO_params), intent(in) :: mesh_file
!
      integer (kind = kint) :: nnod_4_ele
!
!
!      write(*,*) 'allocate_number_of_mesh'
      call allocate_number_of_mesh

!     count number of node for each domain
!
!       write(*,*) 'count_number_w_overlap'
      call alloc_subdomain_groups(mgd_mesh1)
      call count_number_w_overlap(mesh_file, nnod_4_ele)
!
!     array allocation
!
!      write(*,*) 'allocate_geometry_data_4_merge'
      call allocate_geometry_data_4_merge
!
!  set mesh_information
!
!      write(*,*) 'set_geometry_data_2_merge'
      call set_geometry_data_2_merge
!
!      write(*,*) 'count_num_group_w_overlap'
      call allocate_subdomain_grp_stack(mgd_mesh1%num_pe)
      call count_num_group_w_overlap(mgd_mesh1%num_pe,                  &
     &    mgd_mesh1%sub_nod_grp, mgd_mesh1%sub_ele_grp,                 &
     &    mgd_mesh1%sub_surf_grp)
!
!      write(*,*) 'count_merged_mesh_groups'
      call count_merged_mesh_groups(mgd_mesh1%num_pe,                   &
     &    mgd_mesh1%sub_nod_grp, mgd_mesh1%sub_ele_grp,                 &
     &    mgd_mesh1%sub_surf_grp, mgd_mesh1%merged_grp)
!      write(*,*) 'const_merged_mesh_groups'
      call const_merged_mesh_groups(mgd_mesh1%num_pe,                   &
     &    mgd_mesh1%subdomain, mgd_mesh1%merged, merge_tbl,             &
     &    mgd_mesh1%sub_nod_grp, mgd_mesh1%sub_ele_grp,                 &
     &    mgd_mesh1%sub_surf_grp, mgd_mesh1%merged_grp)
!
      call deallocate_subdomain_grp_stack
!
      end subroutine set_merged_mesh_and_group
!
!  ---------------------------------------------------------------------
!
      subroutine set_merged_node_and_element(mesh_file)
!
      type(field_IO_params), intent(in) :: mesh_file
      integer (kind = kint) :: nnod_4_ele
!
!
!      write(*,*) 'allocate_number_of_mesh'
      call allocate_number_of_mesh
!
!     count number of node for each domain
!
!       write(*,*) 'count_number_w_overlap'
      call alloc_subdomain_groups(mgd_mesh1)
      call count_number_w_overlap(mesh_file, nnod_4_ele)
!
!     array allocation
!
!      write(*,*) 'allocate_geometry_data_4_merge'
      call allocate_geometry_data_4_merge
!
!  set mesh_information
!
!       write(*,*) 'set_geometry_data_2_merge'
      call set_geometry_data_2_merge
!
!
      end subroutine set_merged_node_and_element
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine set_overlapped_mesh_and_group(mesh_file, nnod_4_ele)
!
      use const_merged_groups
!
      type(field_IO_params), intent(in) :: mesh_file
      integer (kind = kint), intent(inout) :: nnod_4_ele
!
!
!       write(*,*) 'allocate_number_of_mesh'
      call allocate_number_of_mesh
!
!     count number of node for each domain
!
       write(*,*) 'count_number_w_overlap'
      call alloc_subdomain_groups(mgd_mesh1)
      call count_number_w_overlap(mesh_file, nnod_4_ele)
!
!     array allocation
!
      write(*,*) 'allocate_geometry_data_4_merge'
      call allocate_geometry_data_4_merge
!
!  set mesh_information
!
      write(*,*) 'set_geometry_data_w_overlap'
      call set_geometry_data_w_overlap
!
!
      call allocate_subdomain_grp_stack(mgd_mesh1%num_pe)
      call count_num_group_w_overlap(mgd_mesh1%num_pe,                  &
     &    mgd_mesh1%sub_nod_grp, mgd_mesh1%sub_ele_grp,                 &
     &    mgd_mesh1%sub_surf_grp)
!
      write(*,*) 'count_merged_mesh_groups'
      call count_merged_mesh_groups(mgd_mesh1%num_pe,                   &
     &    mgd_mesh1%sub_nod_grp, mgd_mesh1%sub_ele_grp,                 &
     &    mgd_mesh1%sub_surf_grp, mgd_mesh1%merged_grp)
      write(*,*) 'const_merged_overlapped_groups'
      call const_merged_overlapped_groups(mgd_mesh1%num_pe,             &
     &    mgd_mesh1%sub_nod_grp, mgd_mesh1%sub_ele_grp,                 &
     &    mgd_mesh1%sub_surf_grp, mgd_mesh1%merged_grp)
!
      call deallocate_subdomain_grp_stack
      call dealloc_subdomain_groups(mgd_mesh1)
!
      end subroutine set_overlapped_mesh_and_group
!
!  ---------------------------------------------------------------------
!
      end module set_merged_geometry
