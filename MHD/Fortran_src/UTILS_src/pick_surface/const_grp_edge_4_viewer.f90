!
!      module const_grp_edge_4_viewer
!
!     Written by H. Matsui on Jan., 2007
!
!!      subroutine construct_edge_4_ele_grp                             &
!!     &         (nnod_4_surf, nnod_4_edge, edge_sf_tbl)
!!      subroutine construct_edge_4_surf_grp                            &
!!     &         (nnod_4_surf, nnod_4_edge, edge_sf_tbl)
!        type(sum_hash_tbl), intent(inout) :: edge_sf_tbl
!
      module const_grp_edge_4_viewer
!
      use m_precision
!
      use m_surface_mesh_4_merge
      use m_pickup_table_4_viewer
!
      use t_sum_hash
!
      implicit    none
!
!------------------------------------------------------------------
!
      contains
!
!------------------------------------------------------------------
!
      subroutine construct_edge_4_ele_grp                               &
     &         (nnod_4_surf, nnod_4_edge, edge_sf_tbl)
!
      use set_edge_hash_by_sf
      use set_edge_data_by_sf
!
      integer(kind = kint), intent(in) :: nnod_4_surf, nnod_4_edge
      type(sum_hash_tbl), intent(inout) :: edge_sf_tbl
!
      integer(kind = kint) :: igrp, ngrp, ist, nedge_grp
!
!
      call allocate_ele_grp_edge_item_sf
!
      do igrp = 1, ngrp_ele_sf
        ngrp = ele_stack_sf(igrp*num_pe_sf)                             &
     &        - ele_stack_sf( (igrp-1)*num_pe_sf )
        ist = ele_stack_sf((igrp-1)*num_pe_sf) + 1
!
!   set hash data for edge elements using sum of local node ID
!
        call clear_sum_hash(edge_sf_tbl)
!
!        write(*,*) 'const_part_edge_hash_4_sf', igrp
        call const_part_edge_hash_4_sf                                  &
     &     (nodpetot_viewer, surfpetot_viewer, ngrp,                    &
     &      nnod_4_surf,nnod_4_edge, ie_sf_viewer, ele_item_sf(ist),    &
     &      edge_sf_tbl%num_hash, edge_sf_tbl%istack_hash,              &
     &      edge_sf_tbl%iend_hash, edge_sf_tbl%id_hash,                 &
     &      edge_sf_tbl%iflag_hash)
!
!
        ist = ele_edge_stack_sf( (igrp-1)*num_pe_sf )
!
        call allocate_ele_edge_item_tmp
        ele_edge_item_tmp(1:nedge_ele_sf)                               &
     &          = ele_edge_item_sf(1:nedge_ele_sf)
        call deallocate_ele_grp_edge_item_sf
!
        call count_num_edges_by_sf(nodpetot_viewer, surfpetot_viewer,   &
     &      nnod_4_edge, edge_sf_tbl%istack_hash,                       &
     &      edge_sf_tbl%iend_hash, edge_sf_tbl%iflag_hash, nedge_grp)
        ele_edge_stack_sf(igrp*num_pe_sf)                               &
     &        = ele_edge_stack_sf((igrp-1)*num_pe_sf) + nedge_grp
        nedge_ele_sf = ele_edge_stack_sf(igrp*num_pe_sf)
!
        call allocate_ele_grp_edge_item_sf
        ele_edge_item_sf(1:ist) = ele_edge_item_tmp(1:ist)
        call deallocate_ele_edge_item_tmp
!
!        write(*,*) 'set_part_edges_4_sf', igrp
        call set_part_edges_4_sf(nodpetot_viewer, surfpetot_viewer,     &
     &      nnod_4_edge, nedge_grp, iedge_sf_viewer,                    &
     &      edge_sf_tbl%istack_hash, edge_sf_tbl%iend_hash,             &
     &      edge_sf_tbl%id_hash, edge_sf_tbl%iflag_hash,                &
     &      ele_edge_item_sf(ist+1) )
!
      end do
!
!      write(50,*) 'ele_edge_item_sf', ele_edge_item_sf
!
      end subroutine construct_edge_4_ele_grp
!
!------------------------------------------------------------------
!
      subroutine construct_edge_4_surf_grp                              &
     &         (nnod_4_surf, nnod_4_edge, edge_sf_tbl)
!
      use set_edge_hash_by_sf
      use set_edge_data_by_sf
!
      integer(kind = kint), intent(in) :: nnod_4_surf, nnod_4_edge
      type(sum_hash_tbl), intent(inout) :: edge_sf_tbl
!
      integer(kind = kint) :: igrp, ngrp, ist, nedge_grp
!
!
      call allocate_sf_grp_edge_item_sf
!
      do igrp = 1, ngrp_surf_sf
        ngrp = surf_stack_sf( igrp*num_pe_sf )                          &
     &        - surf_stack_sf( (igrp-1)*num_pe_sf )
        ist = surf_stack_sf( (igrp-1)*num_pe_sf ) + 1
!
!   set hash data for edge elements using sum of local node ID
!
        call clear_sum_hash(edge_sf_tbl)
!
!        write(*,*) 'const_part_edge_hash_4_sf', igrp
        call const_part_edge_hash_4_sf                                  &
     &     (nodpetot_viewer, surfpetot_viewer, ngrp, nnod_4_surf,       &
     &      nnod_4_edge, ie_sf_viewer, surf_item_sf(ist),               &
     &      edge_sf_tbl%num_hash, edge_sf_tbl%istack_hash,              &
     &      edge_sf_tbl%iend_hash, edge_sf_tbl%id_hash,                 &
     &      edge_sf_tbl%iflag_hash)
!
!
        call allocate_sf_edge_item_tmp
        surf_edge_item_tmp(1:nedge_surf_sf)                             &
     &          = surf_edge_item_sf(1:nedge_surf_sf)
        call deallocate_sf_grp_edge_item_sf
!
        call count_num_edges_by_sf(nodpetot_viewer, surfpetot_viewer,   &
     &      nnod_4_edge, edge_sf_tbl%istack_hash,                       &
     &      edge_sf_tbl%iend_hash, edge_sf_tbl%iflag_hash, nedge_grp)
        surf_edge_stack_sf(igrp*num_pe_sf)                              &
     &        = surf_edge_stack_sf((igrp-1)*num_pe_sf) + nedge_grp
        nedge_surf_sf = surf_edge_stack_sf(igrp*num_pe_sf)
!
        call allocate_sf_grp_edge_item_sf
        ist = surf_edge_stack_sf( (igrp-1)*num_pe_sf )
        surf_edge_item_sf(1:ist) = surf_edge_item_tmp(1:ist)
        call deallocate_sf_edge_item_tmp
!
!        write(*,*) 'set_part_edges_4_sf', igrp
        call set_part_edges_4_sf(nodpetot_viewer, surfpetot_viewer,     &
     &      nnod_4_edge, nedge_grp, iedge_sf_viewer,                    &
     &      edge_sf_tbl%istack_hash, edge_sf_tbl%iend_hash,             &
     &      edge_sf_tbl%id_hash, edge_sf_tbl%iflag_hash,                &
     &      surf_edge_item_sf(ist+1) )
      end do
!
      end subroutine construct_edge_4_surf_grp
!
!------------------------------------------------------------------
!
      end module const_grp_edge_4_viewer
