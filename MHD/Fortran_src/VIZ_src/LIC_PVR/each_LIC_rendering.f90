!>@file  each_LIC_rendering.f90
!!       module each_LIC_rendering
!!
!!@author  Y. Liao and H. Matsui
!!@date Programmed in Feb., 2018
!
!> @brief Structures for position in the projection coordinate 
!!
!!@verbatim
!!      subroutine s_each_LIC_rendering                                 &
!!     &         (istep_pvr, mesh, group, ele_mesh, jacs, nod_fld,      &
!!     &          pvr_param, pvr_data)
!!        type(mesh_geometry), intent(in) :: mesh
!!        type(mesh_groups), intent(in) :: group
!!        type(element_geometry), intent(in) :: ele_mesh
!!        type(node_data), intent(in) :: node
!!        type(element_data), intent(in) :: ele
!!        type(surface_data), intent(in) :: surf
!!        type(phys_data), intent(in) :: nod_fld
!!        type(jacobians_type), intent(in) :: jacs
!!        type(PVR_control_params), intent(inout) :: pvr_param
!!        type(PVR_image_generator), intent(inout) :: pvr_data
!!      subroutine dealloc_each_lic_data(pvr_param, pvr_data)
!!        type(PVR_control_params), intent(inout) :: pvr_param
!!        type(PVR_image_generator), intent(inout) :: pvr_data
!!@endverbatim
!
!
      module each_LIC_rendering
!
      use m_precision
      use calypso_mpi
!
      use m_constants
      use m_machine_parameter
      use m_geometry_constants
!
      use t_mesh_data
      use t_phys_data
      use t_jacobians
!
      use t_rendering_vr_image
      use t_control_params_4_pvr
      use t_surf_grp_4_pvr_domain
      use t_pvr_ray_startpoints
      use t_pvr_image_array
      use t_geometries_in_pvr_screen
!
      use set_default_pvr_params
      use set_position_pvr_screen
      use mesh_outline_4_pvr
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
      subroutine s_each_LIC_rendering                                   &
     &         (istep_pvr, mesh, group, ele_mesh, jacs, nod_fld,        &
     &          pvr_param, pvr_data)
!
      use cal_pvr_modelview_mat
      use field_data_4_LIC
      use rendering_LIC_image
      use rendering_streo_LIC_image
!
      integer(kind = kint), intent(in) :: istep_pvr
!
      type(mesh_geometry), intent(in) :: mesh
      type(mesh_groups), intent(in) :: group
      type(element_geometry), intent(in) :: ele_mesh
      type(phys_data), intent(in) :: nod_fld
      type(jacobians_type), intent(in) :: jacs
!
      type(PVR_control_params), intent(inout) :: pvr_param
      type(PVR_image_generator), intent(inout) :: pvr_data
!
      integer(kind = kint) :: ierr
!
!
      if(iflag_debug .gt. 0) write(*,*) 'cal_field_4_pvr'
      call cal_field_4_each_lic                                         &
     &   (mesh%node, mesh%ele, jacs%g_FEM, jacs%jac_3d, nod_fld,        &
     &    pvr_param%field_def, pvr_param%field, ierr)
      if(ierr .gt. 0) call calypso_mpi_abort(ierr,                      &
     &                   'Set vector field for LIC')
!
      if(iflag_debug .gt. 0) write(*,*) 'set_default_pvr_data_params'
      call set_default_pvr_data_params                                  &
     &   (pvr_param%outline, pvr_data%color)
!
      if(pvr_data%view%iflag_rotate_snap .gt. 0) then
        if(pvr_data%view%iflag_stereo_pvr .gt. 0) then
          call streo_lic_rendering_with_rot                             &
     &       (istep_pvr, mesh%node, mesh%ele, ele_mesh%surf, group,     &
     &        pvr_param, pvr_data)
        else
          call lic_rendering_with_rotation                              &
     &       (istep_pvr, mesh%node, mesh%ele, ele_mesh%surf, group,     &
     &        pvr_param, pvr_data)
        end if
      else
        if(pvr_data%view%iflag_stereo_pvr .gt. 0) then
          call streo_lic_rendering_fix_view                             &
     &       (istep_pvr, mesh%node, mesh%ele, ele_mesh%surf,            &
     &        group, pvr_param, pvr_data)
        else
          call lic_rendering_with_fixed_view                            &
     &       (istep_pvr, mesh%node, mesh%ele, ele_mesh%surf,            &
     &        pvr_param, pvr_data)
        end if
      end if
!
      end subroutine s_each_LIC_rendering
!
!  ---------------------------------------------------------------------
!
      subroutine dealloc_each_lic_data(pvr_param, pvr_data)
!
      use set_pvr_control
      use field_data_4_pvr
!
      type(PVR_control_params), intent(inout) :: pvr_param
      type(PVR_image_generator), intent(inout) :: pvr_data
!
!
      if(pvr_data%view%iflag_rotate_snap .eq. 0                         &
     &    .and. pvr_data%view%iflag_stereo_pvr .eq. 0) then
          call flush_rendering_4_fixed_view(pvr_data)
      end if
      call flush_pixel_on_pvr_screen                                    &
     &   (pvr_param%pixel, pvr_data%rgb)
!
      call dealloc_projected_position(pvr_data%screen)
!
      call dealloc_pvr_surf_domain_item(pvr_data%bound)
      call dealloc_nod_data_4_lic(pvr_param%field)
      call flush_each_pvr_control(pvr_data, pvr_param)
!
      end subroutine dealloc_each_lic_data
!
!  ---------------------------------------------------------------------
!
      end module each_LIC_rendering
