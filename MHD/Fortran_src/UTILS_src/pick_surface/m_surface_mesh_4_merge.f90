!
!      module m_surface_mesh_4_merge
!
!      Written by Kemorin
!
!      subroutine allocate_num_mesh_sf
!      subroutine allocate_edge_data_4_sf
!      subroutine allocate_domain_stack_4_surf
!      subroutine allocate_domain_surf_item_sf
!      subroutine allocate_domain_nod_item_sf
!      subroutine allocate_domain_edge_item_sf
!      subroutine allocate_nod_grp_item_4_surf
!      subroutine allocate_ele_grp_item_4_surf
!      subroutine allocate_ele_gp_nod_item_sf
!      subroutine allocate_ele_grp_edge_item_sf
!      subroutine allocate_surf_grp_item_4_surf
!      subroutine allocate_sf_gp_nod_item_sf
!      subroutine allocate_sf_grp_edge_item_sf
!
!      subroutine deallocate_ele_gp_nod_item_sf
!      subroutine deallocate_ele_grp_edge_item_sf
!      subroutine deallocate_sf_grp_edge_item_sf
!
!      subroutine check_edge_connent_viewer(nnod_4_edge)
!
      module m_surface_mesh_4_merge
!
      use m_precision
      use m_constants
      use t_surface_mesh_4_merge
!
      implicit none
!
!
!view_mesh%num_pe_sf
      integer(kind = kint)  :: num_pe_sf
!
      integer(kind=kint ), allocatable :: inod_sf_stack(:)
      integer(kind=kint ), allocatable :: iedge_sf_stack(:)
      integer(kind=kint ), allocatable :: isurf_sf_stack(:)
!
      integer(kind=kint ), allocatable  ::  ie_edge_viewer(:,:)
!
      type(viewer_mesh_data), save :: view_mesh
!
      type(viewer_surface_groups), save :: domain_grps
!
      type(viewer_node_groups), save :: view_nod_grps
      type(viewer_surface_groups), save :: view_ele_grps
      type(viewer_surface_groups), save :: view_sf_grps
!
!
      character (len = kchara) :: surface_file_head = 'in_surface'
      character (len = kchara) :: surface_file_name = 'in_surface.ksm'
!
!------------------------------------------------------------------
!
      contains
!
!------------------------------------------------------------------
!
      subroutine allocate_num_mesh_sf
!
      allocate( inod_sf_stack(0:num_pe_sf)  )
      allocate( isurf_sf_stack(0:num_pe_sf) )
      allocate( iedge_sf_stack(0:num_pe_sf) )
      inod_sf_stack  = 0
      isurf_sf_stack = 0
      iedge_sf_stack = 0
!
      end subroutine allocate_num_mesh_sf
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
      subroutine allocate_edge_data_4_sf(nnod_4_edge)
!
      use m_geometry_constants
!
      integer(kind = kint), intent(in) :: nnod_4_edge
!
      integer(kind = kint) :: num
!
!
      num = view_mesh%edgepetot_viewer
      allocate(ie_edge_viewer(num,nnod_4_edge))
      if(num .gt. 0) ie_edge_viewer = 0
!
      num = view_mesh%surfpetot_viewer
      allocate(view_mesh%iedge_sf_viewer(num,nedge_4_surf))
      if(num .gt. 0) view_mesh%iedge_sf_viewer = 0
!
      end subroutine allocate_edge_data_4_sf
!
!------------------------------------------------------------------
!
      subroutine allocate_domain_stack_4_surf
!
!
      domain_grps%num_grp = 1
!
      call alloc_viewer_surf_grps_stack(num_pe_sf, domain_grps)
!
      domain_grps%grp_name = 'subdomains'
!
      end subroutine allocate_domain_stack_4_surf
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
      subroutine check_edge_connent_viewer(nnod_4_edge)
!
      integer(kind = kint), intent(in) :: nnod_4_edge
      integer(kind = kint) :: iedge
!
!
      write(50,*) 'edgepetot_viewer', view_mesh%edgepetot_viewer
      write(50,*) 'iedge_sf_stack', iedge_sf_stack
      write(50,*) 'ie_edge_viewer'
      do iedge = 1, view_mesh%edgepetot_viewer
        write(50,*) iedge, ie_edge_viewer(iedge,1:nnod_4_edge)
      end do
!
      end subroutine check_edge_connent_viewer
!
!------------------------------------------------------------------
!
      end module m_surface_mesh_4_merge
