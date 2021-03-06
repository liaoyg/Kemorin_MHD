!
!      module m_pickup_table_4_viewer
!
      module m_pickup_table_4_viewer
!
!      Written by Kemorin in Jan., 2007
!
      use m_precision
!
      implicit none
!
!
      integer(kind = kint), allocatable :: imark_surf(:)
      integer(kind = kint), allocatable :: imark_edge(:)
      integer(kind = kint), allocatable :: imark_node(:)
!
      integer(kind = kint), allocatable :: isf_merge2viewer(:)
      integer(kind = kint), allocatable :: isf_viewer2merge(:)
!
      integer(kind = kint), allocatable :: inod_merge2viewer(:)
      integer(kind = kint), allocatable :: inod_viewer2merge(:)
!
!
      integer(kind=kint ), allocatable :: ele_edge_item_tmp(:)
      integer(kind=kint ), allocatable :: surf_edge_item_tmp(:)
!
      integer(kind = kint), allocatable :: ele_nod_item_tmp(:)
      integer(kind = kint), allocatable :: surf_nod_item_tmp(:)
!
!
!      subroutine allocate_imark_surf
!      subroutine allocate_imark_node(numnod)
!      subroutine allocate_ele_edge_item_tmp
!      subroutine allocate_sf_edge_item_tmp
!      subroutine allocate_sf_cvt_table_viewer
!      subroutine allocate_nod_cvt_table_viewer
!
!
!      subroutine deallocate_imark_surf
!      subroutine deallocate_imark_node
!      subroutine deallocate_sf_cvt_table_viewer
!      subroutine deallocate_nod_cvt_table_viewer
!      subroutine deallocate_ele_edge_item_tmp
!      subroutine deallocate_sf_edge_item_tmp
!
!------------------------------------------------------------------
!
      contains
!
!------------------------------------------------------------------
!
      subroutine allocate_imark_surf
!
      use m_surf_geometry_4_merge
!
      allocate( imark_surf(merged_surf%numsurf) )
      imark_surf = 0
!
      end subroutine allocate_imark_surf
!
!------------------------------------------------------------------
!
      subroutine allocate_imark_node(numnod)
!
      integer(kind = kint) :: numnod
!
      allocate( imark_node(numnod) )
      imark_node = 0
!
      end subroutine allocate_imark_node
!
!------------------------------------------------------------------
!
      subroutine allocate_sf_cvt_table_viewer
!
      use m_surf_geometry_4_merge
      use m_surface_mesh_4_merge
!
      allocate( isf_merge2viewer(merged_surf%numsurf) )
      allocate( isf_viewer2merge(surfpetot_viewer) )
      isf_merge2viewer = 0
      isf_viewer2merge = 0
!
      end subroutine allocate_sf_cvt_table_viewer
!
!------------------------------------------------------------------
!
      subroutine allocate_nod_cvt_table_viewer
!
      use m_geometry_data_4_merge
      use m_surface_mesh_4_merge
!
      allocate( inod_merge2viewer(merged%node%numnod) )
      allocate( inod_viewer2merge(nodpetot_viewer) )
      inod_merge2viewer = 0
      inod_viewer2merge = 0
!
      end subroutine allocate_nod_cvt_table_viewer
!
!------------------------------------------------------------------
!
      subroutine allocate_ele_edge_item_tmp
!
      use m_surface_mesh_4_merge
!
      allocate( ele_edge_item_tmp(nedge_ele_sf) )
!
      end subroutine allocate_ele_edge_item_tmp
!
!------------------------------------------------------------------
!
      subroutine allocate_sf_edge_item_tmp
!
      use m_surface_mesh_4_merge
!
      allocate( surf_edge_item_tmp(nedge_surf_sf) )
!
      end subroutine allocate_sf_edge_item_tmp
!
!------------------------------------------------------------------
!
      subroutine allocate_ele_gp_nod_item_tmp
!
      use m_surface_mesh_4_merge
!
      allocate( ele_nod_item_tmp(nnod_ele_sf) )
      ele_nod_item_tmp = 0
!
      end subroutine allocate_ele_gp_nod_item_tmp
!
!------------------------------------------------------------------
!
      subroutine allocate_sf_gp_nod_item_tmp
!
      use m_surface_mesh_4_merge
!
      allocate( surf_nod_item_tmp(nnod_surf_sf) )
      surf_nod_item_tmp = 0
!
      end subroutine allocate_sf_gp_nod_item_tmp
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
      subroutine deallocate_imark_surf
!
      deallocate( imark_surf )
!
      end subroutine deallocate_imark_surf
!
!------------------------------------------------------------------
!
      subroutine deallocate_imark_node
!
      deallocate( imark_node )
!
      end subroutine deallocate_imark_node
!
!------------------------------------------------------------------
!
      subroutine deallocate_sf_cvt_table_viewer
!
      deallocate( isf_merge2viewer )
      deallocate( isf_viewer2merge )
!
      end subroutine deallocate_sf_cvt_table_viewer
!
!------------------------------------------------------------------
!
      subroutine deallocate_nod_cvt_table_viewer
!
      deallocate( inod_merge2viewer )
      deallocate( inod_viewer2merge )
!
      end subroutine deallocate_nod_cvt_table_viewer
!
!------------------------------------------------------------------
!
      subroutine deallocate_ele_edge_item_tmp
!
      deallocate( ele_edge_item_tmp )
!
      end subroutine deallocate_ele_edge_item_tmp
!
!------------------------------------------------------------------
!
      subroutine deallocate_sf_edge_item_tmp
!
      deallocate( surf_edge_item_tmp )
!
      end subroutine deallocate_sf_edge_item_tmp
!
!------------------------------------------------------------------
!
      subroutine deallocate_ele_gp_nod_item_tmp
!
      deallocate( ele_nod_item_tmp )
!
      end subroutine deallocate_ele_gp_nod_item_tmp
!
!------------------------------------------------------------------
!
      subroutine deallocate_sf_gp_nod_item_tmp
!
      deallocate( surf_nod_item_tmp )
!
      end subroutine deallocate_sf_gp_nod_item_tmp
!
!------------------------------------------------------------------
!
      end module m_pickup_table_4_viewer
