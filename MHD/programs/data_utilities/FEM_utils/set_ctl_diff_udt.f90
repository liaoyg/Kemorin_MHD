!
!      module set_ctl_diff_udt
!
!     Written by H. Matsui on July, 2006
!     Modified by H. Matsui on JUne, 2007
!
!      subroutine set_ctl_params_correlate_udt                          &
!     &         (mesh_file, udt_org_param, nod_fld, ucd)
!      subroutine set_ctl_params_diff_udt                               &
!     &         (mesh_file, udt_org_param, ucd)
!      subroutine s_set_ctl_4_diff_udt_steps
!
      module set_ctl_diff_udt
!
      use m_precision
!
      use calypso_mpi
      use m_machine_parameter
      use m_ctl_params_4_diff_udt
      use t_phys_data
      use t_file_IO_parameter
!
      implicit none
!
!   --------------------------------------------------------------------
!
      contains
!
!   --------------------------------------------------------------------
!
      subroutine set_ctl_params_correlate_udt                           &
     &         (mesh_file, udt_org_param, nod_fld, ucd)
!
      use t_ucd_data
      use m_ctl_data_diff_udt
      use m_fem_gauss_int_coefs
      use set_control_nodal_data
      use set_control_ele_layering
!
      type(field_IO_params), intent(inout) ::  mesh_file
      type(field_IO_params), intent(inout) :: udt_org_param
      type(phys_data), intent(inout) :: nod_fld
      type(ucd_data), intent(inout) :: ucd
      integer(kind = kint) :: ierr
!
!
      if (iflag_debug.eq.1) write(*,*) 'set_ctl_params_diff_udt'
      call set_ctl_params_diff_udt(mesh_file, udt_org_param, ucd)
!
      if (iflag_debug.eq.1) write(*,*) 's_set_control_ele_layering'
      call s_set_control_ele_layering(elayer_d_ctl)
!
      if (iflag_debug.eq.1) write(*,*) 's_set_control_nodal_data'
      call s_set_control_nodal_data(fld_d_ctl%field_ctl, nod_fld, ierr)
      if (ierr .ne. 0) call calypso_MPI_abort(ierr, e_message)
!
      if (iflag_debug.eq.1) write(*,*) 's_set_ctl_4_diff_udt_steps'
      call s_set_ctl_4_diff_udt_steps(t_d_ctl)
!
      if(fint_d_ctl%integration_points_ctl%iflag .gt. 0) then
        call maximum_integration_points                                 &
     &     (fint_d_ctl%integration_points_ctl%intvalue)
      end if
!
      end subroutine set_ctl_params_correlate_udt
!
!   --------------------------------------------------------------------
!
      subroutine set_ctl_params_diff_udt                                &
     &         (mesh_file, udt_org_param, ucd)
!
      use t_ucd_data
      use t_field_data_IO
      use m_ctl_data_4_platforms
      use m_ctl_data_4_org_data
      use m_ctl_data_diff_udt
      use m_geometry_constants
      use m_file_format_switch
      use set_ctl_parallel_platform
      use set_control_platform_data
      use set_ctl_params_2nd_files
      use ucd_IO_select
!
      type(field_IO_params), intent(inout) ::  mesh_file
      type(field_IO_params), intent(inout) :: udt_org_param
      type(ucd_data), intent(inout) :: ucd
!
!
      call turn_off_debug_flag_by_ctl(my_rank, plt1)
      call check_control_num_domains(plt1)
      call set_control_smp_def(my_rank, plt1)
      call set_control_mesh_def(plt1, mesh_file)
      call set_control_org_udt_file_def(udt_org_param)
!
!
      call set_ucd_file_define(plt1, ucd)
!
!   set field data name
!
      if (i_ref_udt_head_ctl .ne. 0) then
        ref_udt_file_head = ref_udt_head_ctl
        if (iflag_debug.gt.0)                                           &
     &   write(*,*) 'ref_udt_file_head: ', trim(ref_udt_file_head)
      end if
!
      if (i_tgt_udt_head_ctl .ne. 0) then
        tgt_udt_file_head = tgt_udt_head_ctl
        if (iflag_debug.gt.0)                                           &
     &   write(*,*) 'tgt_udt_file_head: ', trim(tgt_udt_file_head)
      end if
!
      if (i_group_mesh_head .ne. 0) then
        grouping_mesh_head = group_mesh_head_ctl
        if (iflag_debug.gt.0)                                           &
     &   write(*,*) 'grouping_mesh_head: ', trim(grouping_mesh_head)
      end if
!
!   field setting
!
      if (plt1%field_file_prefix%iflag .ne. 0) then
        diff_udt_file_head = plt1%field_file_prefix%charavalue
        ave_udt_file_head =  plt1%field_file_prefix%charavalue
        prod_udt_file_head = plt1%field_file_prefix%charavalue
      else
        diff_udt_file_head = "field_diff/out"
        ave_udt_file_head =  "out_average"
        prod_udt_file_head = "field_new/out"
      end if
!
      call choose_ucd_file_format(plt1%field_file_fmt_ctl%charavalue,   &
     &    plt1%field_file_fmt_ctl%iflag, ifmt_diff_udt_file)
!
      if (i_prod_name .ne. 0) then
        product_field_name = product_field_ctl
        if (iflag_debug.gt.0)                                           &
     &   write(*,*) 'product_field_name ', trim(product_field_name)
      end if
!
      if (i_corr_fld_name .ne. 0) then
        correlate_field_name = correlate_fld_ctl
        if (iflag_debug.gt.0)                                           &
     &   write(*,*) 'correlate_field_name ', trim(correlate_field_name)
      end if
!
      if (i_corr_cmp_name .ne. 0) then
        correlate_comp_name = correlate_cmp_ctl
        if (iflag_debug.gt.0)                                           &
     &   write(*,*) 'correlate_comp_name ', trim(correlate_comp_name)
      end if
!
      if (i_correlate_coord .ne. 0) then
        if     (correlate_coord_ctl .eq. 'cartesian'                    &
     &     .or. correlate_coord_ctl .eq. 'Cartesian'                    &
     &     .or. correlate_coord_ctl .eq. 'CARTESIAN'                    &
     &     .or. correlate_coord_ctl .eq. 'xyz'                          &
     &     .or. correlate_coord_ctl .eq. 'XYZ') then
          iflag_correlate_coord = iflag_certecian
        else if(correlate_coord_ctl .eq. 'spherical'                    &
     &     .or. correlate_coord_ctl .eq. 'Spherical'                    &
     &     .or. correlate_coord_ctl .eq. 'SPHERICAL'                    &
     &     .or. correlate_coord_ctl .eq. 'rtp'                          &
     &     .or. correlate_coord_ctl .eq. 'RTP') then
          iflag_correlate_coord = iflag_spherical
        else if(correlate_coord_ctl .eq. 'cyrindrical'                  &
     &     .or. correlate_coord_ctl .eq. 'Cyrindrical'                  &
     &     .or. correlate_coord_ctl .eq. 'CYRINDRICAL'                  &
     &     .or. correlate_coord_ctl .eq. 'spz'                          &
     &     .or. correlate_coord_ctl .eq. 'SPZ') then
          iflag_correlate_coord = iflag_cylindrical
        end if
      else
        iflag_correlate_coord = 0
      end if
      if (iflag_debug.gt.0)                                             &
     &   write(*,*) 'iflag_correlate_coord ', iflag_correlate_coord
!
!
      end subroutine set_ctl_params_diff_udt
!
!  ---------------------------------------------------------------------
!
      subroutine s_set_ctl_4_diff_udt_steps(tctl)
!
      use m_t_step_parameter
      use m_error_IDs
      use m_t_int_parameter
      use t_ctl_data_4_time_steps
      use cal_num_digits
!
      type(time_data_control), intent(in) :: tctl
!
!   parameters for time evolution
!
        i_step_init   = 0
        if (tctl%i_step_init_ctl%iflag .gt. 0) then
          i_step_init   = tctl%i_step_init_ctl%intvalue
        end if
!
        if (tctl%i_step_number_ctl%iflag .eq. 0) then
          e_message = 'Set step number to finish'
            call calypso_MPI_abort(ierr_evo, e_message)
        else
          i_step_number = tctl%i_step_number_ctl%intvalue
        end if
!
        i_step_output_ucd = 1
        if (tctl%i_step_ucd_ctl%iflag .gt. 0) then
          i_step_output_ucd = tctl%i_step_ucd_ctl%intvalue
        end if
!
        i_diff_steps = 1
        if (tctl%i_diff_steps_ctl%iflag .gt. 0) then
          i_diff_steps = tctl%i_diff_steps_ctl%intvalue
        end if
!
        dt = 1.0d0
        if (tctl%dt_ctl%iflag .gt. 0) then
          dt = tctl%dt_ctl%realvalue
        end if
!
        if (dt .eq. 0.0d0) then
          ddt = 1.0d30
        else
          ddt = 1.0d0 / dt
        end if
!
        if (iflag_debug.eq.1) then
          write(*,*) 'i_step_init ',       i_step_init
          write(*,*) 'i_step_number ',     i_step_number
          write(*,*) 'i_step_output_ucd ', i_step_output_ucd
          write(*,*) 'i_diff_steps ',      i_diff_steps
        end if
!
      end subroutine s_set_ctl_4_diff_udt_steps
!
! -----------------------------------------------------------------------
!
      end module set_ctl_diff_udt
