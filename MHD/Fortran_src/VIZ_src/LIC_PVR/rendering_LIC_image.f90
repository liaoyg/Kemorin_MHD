!>@file  rendering_LIC_image.f90
!!       module rendering_LIC_image
!!
!!@author  Y. Liao and H. Matsui
!!@date Programmed in Feb., 2018
!
!> @brief Structures for position in the projection coordinate 
!!
!!@verbatim
!!      subroutine lic_rendering_with_rotation(istep_pvr,               &
!!     &          node, ele, surf, group, lic_p, pvr_param, pvr_data)
!!      subroutine lic_rendering_with_fixed_view(istep_pvr,             &
!!     &          node, ele, surf, lic_p, pvr_param, pvr_data)
!!
!!      subroutine rendering_lic_at_once(isel_projection, istep_pvr,    &
!!     &           node, ele, surf, group, lic_p, pvr_param, pvr_data)
!!        type(node_data), intent(in) :: node
!!        type(element_data), intent(in) :: ele
!!        type(surface_data), intent(in) :: surf
!!        type(mesh_groups), intent(in) :: group
!!        type(lic_parameters), intent(in) :: lic_p
!!        type(PVR_control_params), intent(in) :: pvr_param
!!        type(PVR_image_generator), intent(inout) :: pvr_data
!!@endverbatim
!
      module rendering_LIC_image
!
      use m_precision
      use m_machine_parameter
      use m_constants
      use m_work_time
!
      use calypso_mpi
!
      use t_mesh_data
      use t_geometry_data
      use t_surface_data
      use t_group_data
      use t_control_params_4_pvr
      use t_control_param_LIC
      use t_geometries_in_pvr_screen
      use t_surf_grp_4_pvr_domain
      use t_pvr_ray_startpoints
      use t_pvr_image_array
      use t_rendering_vr_image
      use generate_vr_image
!
      implicit  none
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine lic_rendering_with_rotation(istep_pvr,                 &
     &          node, ele, surf, group, lic_p, pvr_param, pvr_data)
!
      use cal_pvr_modelview_mat
      use composite_pvr_images
      use write_PVR_image
!
      integer(kind = kint), intent(in) :: istep_pvr
!
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(surface_data), intent(in) :: surf
      type(mesh_groups), intent(in) :: group
      type(lic_parameters), intent(in) :: lic_p
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
!
        call rendering_lic_at_once(IFLAG_NORMAL, istep_pvr,             &
     &       node, ele, surf, group, lic_p, pvr_param, pvr_data)
!
        call end_elapsed_time(76)
        call start_elapsed_time(77)
        if(iflag_debug .gt. 0) write(*,*) 'sel_write_pvr_image_file'
        call sel_write_pvr_image_file                                   &
     &   (pvr_param%file, i_rot, istep_pvr, IFLAG_NORMAL, pvr_data%rgb)
        call calypso_mpi_barrier
        call end_elapsed_time(77)
        call start_elapsed_time(76)
!
        call dealloc_pvr_local_subimage(pvr_data%image)
        call deallocate_pvr_ray_start(pvr_data%start_pt)
      end do
!
      end subroutine lic_rendering_with_rotation
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine lic_rendering_with_fixed_view(istep_pvr,               &
     &          node, ele, surf, lic_p, pvr_param, pvr_data)
!
      use composite_pvr_images
      use write_LIC_image
      use write_PVR_image
!
      integer(kind = kint), intent(in) :: istep_pvr
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(surface_data), intent(in) :: surf
      type(lic_parameters), intent(in) :: lic_p
      type(PVR_control_params), intent(in) :: pvr_param
!
      type(PVR_image_generator), intent(inout) :: pvr_data
!
!
      call copy_item_pvr_ray_start                                      &
     &   (pvr_data%start_pt_1, pvr_data%start_pt)
!
      if(iflag_debug .gt. 0) write(*,*) 'rendering_image_4_lic'
      call rendering_image_4_lic(istep_pvr, pvr_param%file,             &
     &    node, ele, surf, lic_p, pvr_data%color, pvr_param%colorbar,   &
     &    pvr_param%field, pvr_data%screen, pvr_data%start_pt,          &
     &    pvr_data%image, pvr_data%rgb)
!
      call end_elapsed_time(76)
      call start_elapsed_time(77)
      if(iflag_debug .gt. 0) write(*,*) 'sel_write_pvr_image_file'
      call sel_write_pvr_image_file(pvr_param%file, iminus,             &
     &    istep_pvr, IFLAG_NORMAL, pvr_data%rgb)
!
      if(pvr_param%file%iflag_monitoring .gt. 0) then
        call sel_write_pvr_image_file(pvr_param%file, iminus,           &
     &      iminus, IFLAG_NORMAL, pvr_data%rgb)
      end if
      call calypso_mpi_barrier
      call end_elapsed_time(77)
      call start_elapsed_time(76)
!
      end subroutine lic_rendering_with_fixed_view
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine rendering_lic_at_once(isel_projection, istep_pvr,      &
     &           node, ele, surf, group, lic_p, pvr_param, pvr_data)
!
      use cal_pvr_modelview_mat
      use composite_pvr_images
      use write_LIC_image
!
      integer(kind = kint), intent(in) :: isel_projection
      integer(kind = kint), intent(in) :: istep_pvr
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(surface_data), intent(in) :: surf
      type(mesh_groups), intent(in) :: group
      type(lic_parameters), intent(in) :: lic_p
      type(PVR_control_params), intent(in) :: pvr_param
!
      type(PVR_image_generator), intent(inout) :: pvr_data
!
!
      call transfer_to_screen(isel_projection,                          &
     &      node, ele, surf, group%surf_grp, group%surf_grp_geom,       &
     &      pvr_param%field, pvr_data%view, pvr_param%pixel,            &
     &      pvr_data%bound, pvr_data%screen, pvr_data%start_pt)
      call set_subimages(pvr_data%rgb%num_pixel_xy,                     &
     &      pvr_data%start_pt, pvr_data%image)
!
      if(iflag_debug .gt. 0) write(*,*) 'rendering_image_4_lic'
      call rendering_image_4_lic(istep_pvr, pvr_param%file,             &
     &    node, ele, surf, lic_p, pvr_data%color, pvr_param%colorbar,   &
     &    pvr_param%field, pvr_data%screen, pvr_data%start_pt,          &
     &    pvr_data%image, pvr_data%rgb)
!
      end subroutine rendering_lic_at_once
!
!  ---------------------------------------------------------------------
!
      end module rendering_LIC_image