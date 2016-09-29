!>@file  t_rendering_vr_image.f90
!!       module t_rendering_vr_image
!!
!!@author H. Matsui
!!@date   Programmed in Aug., 2011
!
!> @brief Structures for position in the projection coordinate 
!!
!!@verbatim
!!      subroutine set_fixed_view_and_image(node, ele, surf, group,     &
!!     &          pvr_param, pvr_data)
!!      subroutine rendering_with_fixed_view(istep_pvr, node, ele, surf,&
!!     &          pvr_param, pvr_data)
!!      subroutine rendering_with_rotation(istep_pvr, node, ele, surf,  &
!!     &          group, pvr_param, pvr_data)
!!        type(node_data), intent(in) :: node
!!        type(element_data), intent(in) :: ele
!!        type(surface_data), intent(in) :: surf
!!        type(mesh_groups), intent(in) :: group
!!        type(PVR_control_params), intent(in) :: pvr_param
!!        type(PVR_image_generator), intent(inout) :: pvr_data
!!@endverbatim
!
      module t_rendering_vr_image
!
      use m_precision
      use m_machine_parameter
      use m_constants
!
      use calypso_mpi
!
      use t_mesh_data
      use t_geometry_data
      use t_surface_data
      use t_group_data
      use t_control_params_4_pvr
      use t_geometries_in_pvr_screen
      use t_surf_grp_4_pvr_domain
      use t_pvr_ray_startpoints
      use t_pvr_image_array
      use generate_vr_image
!
      implicit  none
!
!>      Structure of PVR control parameters
      type PVR_control_params
!>        Parameters for output files
        type(pvr_output_parameter) :: file
!>        Parameters for image pixels
        type(pvr_pixel_position_type) :: pixel
!>        Structure for field parameter for PVR
        type(pvr_field_parameter) :: field_def
!>        Structure for rough serch of subdomains
        type(pvr_domain_outline) :: outline
!>        Field data for volume rendering
        type(pvr_projected_field) :: field
!>        Structure for PVR colormap
        type(pvr_colorbar_parameter):: colorbar
      end type PVR_control_params
!
!
!>      Structure of PVR image generation
      type PVR_image_generator
!>        Viewer coordinate information
        type(pvr_view_parameter) :: view
!>        color paramter for volume rendering
        type(pvr_colormap_parameter) :: color
!>        Domain boundary information
        type(pvr_bounds_surf_ctl) :: bound
!>        Data on screen oordinate
        type(pvr_projected_data) :: screen
!>        Start point structure for volume rendering
        type(pvr_ray_start_type) :: start_pt
!>        Stored start point structure for volume rendering
        type(pvr_ray_start_type) :: start_pt_saved
!
!>        Pixel data structure for volume rendering
        type(pvr_image_type) :: image
      end type PVR_image_generator
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine set_fixed_view_and_image(node, ele, surf, group,       &
     &          pvr_param, pvr_data)
!
      use cal_pvr_modelview_mat
!
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(surface_data), intent(in) :: surf
      type(mesh_groups), intent(in) :: group
      type(PVR_control_params), intent(in) :: pvr_param
!
      type(PVR_image_generator), intent(inout) :: pvr_data
!
!
      call cal_pvr_modelview_matrix                                     &
     &   (izero, pvr_param%outline, pvr_data%view, pvr_data%color,      &
     &    pvr_data%screen)
!
      call transfer_to_screen(IFLAG_NORMAL,                             &
     &    node, ele, surf, group%surf_grp, group%surf_grp_geom,         &
     &    pvr_param%field, pvr_data%view, pvr_data%bound,               &
     &    pvr_param%pixel, pvr_data%screen, pvr_data%start_pt)
!
      call set_subimages(pvr_data%start_pt, pvr_data%image)
!
      pvr_data%start_pt_saved%num_pvr_ray                               &
     &               = pvr_data%start_pt%num_pvr_ray
      call allocate_item_pvr_ray_start(pvr_data%start_pt_saved)
      call copy_item_pvr_ray_start                                      &
     &   (pvr_data%start_pt, pvr_data%start_pt_saved)
!
      end subroutine set_fixed_view_and_image
!
!  ---------------------------------------------------------------------
!
      subroutine rendering_with_fixed_view(istep_pvr, node, ele, surf,  &
     &          pvr_param, pvr_data)
!
      use composite_pvr_images
      use set_pvr_ray_start_point
!
      integer(kind = kint), intent(in) :: istep_pvr
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(surface_data), intent(in) :: surf
      type(PVR_control_params), intent(in) :: pvr_param
!
      type(PVR_image_generator), intent(inout) :: pvr_data
!
      integer(kind = kint), parameter :: i_rot = -1
!
!
      call copy_item_pvr_ray_start                                      &
     &   (pvr_data%start_pt_saved, pvr_data%start_pt)
!
      call rendering_image(i_rot, istep_pvr, node, ele, surf,           &
     &    pvr_param%file, pvr_data%color, pvr_param%colorbar,           &
     &    pvr_data%view, pvr_param%field, pvr_data%screen,              &
     &    pvr_data%start_pt,pvr_data%image)
!
!      call dealloc_pvr_local_subimage(pvr_data%image)
!
      end subroutine rendering_with_fixed_view
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine rendering_with_rotation(istep_pvr, node, ele, surf,    &
     &          group, pvr_param, pvr_data)
!
      use cal_pvr_modelview_mat
!
      integer(kind = kint), intent(in) :: istep_pvr
!
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(surface_data), intent(in) :: surf
      type(mesh_groups), intent(in) :: group
      type(PVR_control_params), intent(in) :: pvr_param
!
      type(PVR_image_generator), intent(inout) :: pvr_data
!
!
      integer(kind = kint) :: i_rot, ist_rot, ied_rot
!
      ist_rot = pvr_data%view%istart_rot
      ied_rot = pvr_data%view%iend_rot
      do i_rot = ist_rot, ied_rot
        call cal_pvr_modelview_matrix                                   &
     &     (i_rot, pvr_param%outline, pvr_data%view, pvr_data%color,    &
     &      pvr_data%screen)
        call transfer_to_screen(IFLAG_NORMAL,                           &
     &      node, ele, surf, group%surf_grp, group%surf_grp_geom,       &
     &      pvr_param%field, pvr_data%view, pvr_data%bound,             &
     &      pvr_param%pixel, pvr_data%screen, pvr_data%start_pt)
        call set_subimages(pvr_data%start_pt, pvr_data%image)
!
        call rendering_image(i_rot, istep_pvr, node, ele, surf,         &
     &      pvr_param%file, pvr_data%color, pvr_param%colorbar,         &
     &      pvr_data%view, pvr_param%field, pvr_data%screen,            &
     &      pvr_data%start_pt, pvr_data%image)
!
        call dealloc_pvr_local_subimage(pvr_data%image)
        call deallocate_pvr_ray_start(pvr_data%start_pt)
      end do
!
      end subroutine rendering_with_rotation
!
!  ---------------------------------------------------------------------
!
      end module t_rendering_vr_image
