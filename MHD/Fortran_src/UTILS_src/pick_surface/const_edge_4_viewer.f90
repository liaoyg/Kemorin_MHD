!
!      module const_edge_4_viewer
!
!     Written by H. Matsui on Jan., 2007
!
!      subroutine construct_edge_4_viewer(surf, edge)
!
      module const_edge_4_viewer
!
      use m_precision
!
      use m_geometry_constants
      use m_surface_mesh_4_merge
!
      use t_sum_hash
!
      implicit    none
!
      type(sum_hash_tbl), save :: edge_sf_tbl
!
      private :: edge_sf_tbl
      private :: const_all_edge_4_viewer, construct_edge_4_domain
!
!------------------------------------------------------------------
!
      contains
!
!------------------------------------------------------------------
!
      subroutine construct_edge_4_viewer(surf, edge)
!
      use t_surface_data
      use t_edge_data
      use const_grp_edge_4_viewer
      use count_edge_domain_4_viewer
!
      type(surface_data), intent(in) :: surf
      type(edge_data), intent(in) :: edge
!
!
      call alloc_sum_hash(nodpetot_viewer, surfpetot_viewer,            &
     &    nedge_4_surf, edge%nnod_4_edge, edge_sf_tbl)
!
      call const_all_edge_4_viewer(surf, edge, edge_sf_tbl)
!
      call construct_edge_4_domain                                      &
     &   (surf%nnod_4_surf, edge%nnod_4_edge, edge_sf_tbl)
!
      call construct_edge_4_ele_grp                                     &
     &   (surf%nnod_4_surf, edge%nnod_4_edge, edge_sf_tbl)
      call construct_edge_4_surf_grp                                    &
     &   (surf%nnod_4_surf, edge%nnod_4_edge, edge_sf_tbl)
!
      call dealloc_sum_hash(edge_sf_tbl)
!
         write(*,*)  'count_nedge_4_each_domain'
      call count_nedge_4_each_domain
         write(*,*)  'count_nedge_domain_4_domain'
      call count_nedge_domain_4_domain
         write(*,*)  'count_nedge_ele_grp_4_domain'
      call count_nedge_ele_grp_4_domain
         write(*,*)  'count_nedge_surf_grp_4_domain'
      call count_nedge_surf_grp_4_domain
!
      end subroutine construct_edge_4_viewer
!
!------------------------------------------------------------------
!
      subroutine const_all_edge_4_viewer(surf, edge, ed_sf_tbl)
!
      use t_surface_data
      use t_edge_data
      use set_edge_hash_by_sf
      use set_edge_data_by_sf
!
      type(surface_data), intent(in) :: surf
      type(edge_data), intent(in) :: edge
!
      type(sum_hash_tbl), intent(inout) :: ed_sf_tbl
!
!   set hash data for edge elements using sum of local node ID
!
      call clear_sum_hash(ed_sf_tbl)
!
      write(*,*) 'const_edge_hash_4_sf'
      call const_edge_hash_4_sf(nodpetot_viewer, surfpetot_viewer,      &
     &    surf%nnod_4_surf, edge%nnod_4_edge, ie_sf_viewer,             &
     &    ed_sf_tbl%num_hash, ed_sf_tbl%istack_hash,                    &
     &    ed_sf_tbl%iend_hash, ed_sf_tbl%id_hash, ed_sf_tbl%iflag_hash)
!
      write(*,*) 'count_num_edges_by_sf'
      call count_num_edges_by_sf                                        &
     &   (nodpetot_viewer, surfpetot_viewer, edge%nnod_4_edge,          &
     &    ed_sf_tbl%istack_hash, ed_sf_tbl%iend_hash,                   &
     &    ed_sf_tbl%iflag_hash, edgepetot_viewer)
!
      call allocate_edge_data_4_sf(edge%nnod_4_edge)
!
      write(*,*) 'set_edges_connect_by_sf'
      call set_edges_connect_by_sf                                      &
     &   (nodpetot_viewer, surfpetot_viewer, edgepetot_viewer,          &
     &    surf%nnod_4_surf, edge%nnod_4_edge, ie_sf_viewer,             &
     &    ed_sf_tbl%istack_hash, ed_sf_tbl%iend_hash,                   &
     &    ed_sf_tbl%id_hash, ed_sf_tbl%iflag_hash,                      &
     &    ie_edge_viewer, iedge_sf_viewer, edge%node_on_edge_sf)
!
      end subroutine const_all_edge_4_viewer
!
!------------------------------------------------------------------
!
      subroutine construct_edge_4_domain                                &
     &         (nnod_4_surf, nnod_4_edge, ed_sf_tbl)
!
      use set_edge_hash_by_sf
      use set_edge_data_by_sf
!
      integer(kind = kint), intent(in) :: nnod_4_surf, nnod_4_edge
      type(sum_hash_tbl), intent(inout) :: ed_sf_tbl
!
!   set hash data for edge elements using sum of local node ID
!
      call clear_sum_hash(ed_sf_tbl)
!
!      write(*,*) 'const_part_edge_hash_4_sf'
      call const_part_edge_hash_4_sf                                    &
     &   (nodpetot_viewer, surfpetot_viewer, nsurf_domain_sf,           &
     &    nnod_4_surf, nnod_4_edge, ie_sf_viewer, isurf_domain_sf,      &
     &    ed_sf_tbl%num_hash, ed_sf_tbl%istack_hash,                    &
     &    ed_sf_tbl%iend_hash, ed_sf_tbl%id_hash, ed_sf_tbl%iflag_hash)
!
!
      call count_num_edges_by_sf                                        &
     &   (nodpetot_viewer, surfpetot_viewer, nnod_4_edge,               &
     &    ed_sf_tbl%istack_hash, ed_sf_tbl%iend_hash,                   &
     &    ed_sf_tbl%iflag_hash, nedge_domain_sf)
!
      call allocate_domain_edge_item_sf
!
!      write(*,*) 'set_part_edges_4_sf'
      call set_part_edges_4_sf(nodpetot_viewer, surfpetot_viewer,       &
     &    nnod_4_edge, nedge_domain_sf, iedge_sf_viewer,                &
     &    ed_sf_tbl%istack_hash, ed_sf_tbl%iend_hash,                   &
     &    ed_sf_tbl%id_hash, ed_sf_tbl%iflag_hash,                      &
     &    edge_item_domain_sf)
!
      end subroutine construct_edge_4_domain
!
!------------------------------------------------------------------
!
      end module const_edge_4_viewer
