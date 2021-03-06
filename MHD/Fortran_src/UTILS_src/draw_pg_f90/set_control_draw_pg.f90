!
!      module set_control_draw_pg
!
!      Written by H. Matsui
!
!      subroutine s_set_control_draw_pg
!      subroutine set_control_draw_zplane
!      subroutine set_control_draw_map
!      subroutine set_psffield_id_4_plot_pg(psf_phys)
!
      module set_control_draw_pg
!
      use m_precision
!
      use m_machine_parameter
      use m_ctl_param_plot_pg
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine s_set_control_draw_pg
!
      use m_ctl_data_plot_pg
      use m_isoline_dat_pg
      use m_file_format_switch
      use m_field_file_format
      use set_components_flags
      use skip_comment_f
!
      integer(kind= kint) :: i, j
      character(len=kchara) :: tmpchara
!
!
      start_time_pg = 0.0d0
      if(t_pg_ctl%time_init_ctl%iflag .gt. 0) then
        start_time_pg = t_pg_ctl%time_init_ctl%realvalue
      end if
!
      delta_time_pg = 0.0d0
      if(t_pg_ctl%dt_ctl%iflag .gt. 0) then
        delta_time_pg = t_pg_ctl%dt_ctl%realvalue
      end if
!
      ist_pg = 1
      if(t_pg_ctl%i_step_init_ctl%iflag .gt. 0) then
        ist_pg = t_pg_ctl%i_step_init_ctl%intvalue
      end if
!
      ied_pg = 1
      if(t_pg_ctl%i_step_number_ctl%iflag .gt. 0) then
        ied_pg = t_pg_ctl%i_step_number_ctl%intvalue
      end if
!
      inc_pg = 1
      if(t_pg_ctl%i_step_psf_ctl%iflag .gt. 0) then
        inc_pg = t_pg_ctl%i_step_psf_ctl%intvalue
      end if
!
!
      npanel_window = 1
      if(num_panels_ctl%iflag .gt. 0) then
        npanel_window = num_panels_ctl%intvalue
      end if
!
!
      idisp_mode = 1
      if(contour_type_ctl%iflag .gt. 0) then
        tmpchara = contour_type_ctl%charavalue
        if     (cmp_no_case(tmpchara, 'Both')                           &
     &     .or. cmp_no_case(tmpchara, 'Line_and_Fill')) then
          idisp_mode = 3
        else if(cmp_no_case(tmpchara, 'Fill')                           &
     &     .or. cmp_no_case(tmpchara, 'Filled')       ) then
          idisp_mode = 2
        end if
      end if
!
!
      icolor_mode = 1
      if(color_mode_ctl%iflag .gt. 0) then
        tmpchara = color_mode_ctl%charavalue
        if     (cmp_no_case(tmpchara, 'Rainbow')                        &
     &     .or. cmp_no_case(tmpchara, 'Color')        ) then
          icolor_mode = 1
        else if(cmp_no_case(tmpchara, 'Yellow_Green') ) then
          icolor_mode = -1
        else if(cmp_no_case(tmpchara, 'Grayscale')                      &
     &     .or. cmp_no_case(tmpchara, 'Gray')         ) then
          icolor_mode = 0
        end if
      end if
!
!
      if(psf_file_head_ctl%iflag .gt. 0) then
        pg_psf_file_prefix = psf_file_head_ctl%charavalue
      else
        write(*,*) 'set file header for psf data'
        stop
      end if
      call choose_ucd_file_format(psf_data_fmt_ctl%charavalue,          &
     &    psf_data_fmt_ctl%iflag, iflag_pg_psf_fmt)
!
      if(map_grid_file_ctl%iflag .gt. 0) then
        fhead_map_grid =  map_grid_file_ctl%charavalue
      end if
!
      if(plot_field_ctl%icou .gt. 0) then
        ntot_plot_pg = plot_field_ctl%num
      else
        write(*,*) 'set number of component to plot'
        stop
      end if
!
      call allocate_plot_param_pg
!
        field_name_4_plot(1:ntot_plot_pg)                               &
     &       = plot_field_ctl%c1_tbl(1:ntot_plot_pg)
        comp_name_4_plot(1:ntot_plot_pg)                                &
     &       = plot_field_ctl%c2_tbl(1:ntot_plot_pg)
        field_label_4_plot(1:ntot_plot_pg)                              &
     &       = plot_field_ctl%c3_tbl(1:ntot_plot_pg)
!
      do i = 1, ntot_plot_pg
        call s_set_components_flags( comp_name_4_plot(i),               &
     &      field_name_4_plot(i), id_comp_4_plot(i),                    &
     &      num_comp_4_plot(i), ncomp_org_4_plot(i),                    &
     &      viz_name_4_plot(i) )
      end do
!
!
      if(contour_range_ctl%num .gt. 0) then
        do i = 1, contour_range_ctl%num
          j = contour_range_ctl%int1(i)
          num_line_pg(j) =  contour_range_ctl%int2(i)
          range_pg(1,j) =   contour_range_ctl%vec1(i)
          range_pg(2,j) =   contour_range_ctl%vec2(i)
        end do
!
        do i = 1, contour_range_ctl%num
          nmax_line = max(contour_range_ctl%int2(i),nmax_line)
        end do
        call dealloc_control_array_i2_r2(contour_range_ctl)
      end if
!
      call allocate_data_4_isoline
!
      if(vector_scale_ctl%num .gt. 0) then
        do i = 1, vector_scale_ctl%num
          j = vector_scale_ctl%int1(i)
          nskip_vect_pg(j) =  vector_scale_ctl%int2(i)
          scale_pg(j) = vector_scale_ctl%vect(i)
        end do
!
        call dealloc_control_array_i2_r(vector_scale_ctl)
      end if
!
      if (iflag_debug .gt. 0) then
        write(*,*) 'npanel_window', npanel_window
        write(*,*) 'idisp_mode',    idisp_mode
        write(*,*) 'icolor_mode',   icolor_mode
        write(*,*) 'ist_pg',        ist_pg
        write(*,*) 'ied_pg',        ied_pg
        write(*,*) 'inc_pg',        inc_pg
!
        write(*,*) 'nmax_line',            nmax_line
        write(*,*) 'iflag_pg_psf_fmt',   iflag_pg_psf_fmt
        write(*,*) 'fhead_map_grid',  trim(fhead_map_grid)
        write(*,*) 'pg_psf_file_prefix', trim(pg_psf_file_prefix)
        do i = 1, ntot_plot_pg
          write(*,*) 'field_name', i, trim(field_name_4_plot(i))
          write(*,*) 'comp_name', trim(comp_name_4_plot(i))
          write(*,*) 'field_label', trim(field_label_4_plot(i))
          write(*,*) 'viz_name_4_plot', trim(viz_name_4_plot(i))
          write(*,*) 'viz_ids', id_comp_4_plot(i), num_comp_4_plot(i),  &
     &               ncomp_org_4_plot(i)
          write(*,*)
        end do
        do i = 1, vector_scale_ctl%num
          write(*,*) 'interval, scale', i, nskip_vect_pg(i), scale_pg(i)
        end do
      end if
!
      end subroutine s_set_control_draw_pg
!
!-----------------------------------------------------------------------
!
      subroutine set_control_draw_zplane
!
      use m_ctl_data_plot_pg
!
!
      if(outer_radius_ctl%iflag .gt. 0) then
        shell_size = outer_radius_ctl%realvalue
      else
        shell_size = 20.0d0 / 13.0d0
      end if
!
      if(ro_ri_ratio_ctl%iflag .gt. 0) then
        shell_ratio = ro_ri_ratio_ctl%realvalue
      else
        shell_ratio = 0.35d0
      end if
!
      if(pg_plane_size_ctl%iflag .gt. 0) then
        plane_size(1:2) = pg_plane_size_ctl%realvalue(1:2)
      else
        plane_size(1:2) = 1.0d0
      end if
!
      if (iflag_debug .gt. 0) then
        write(*,*) 'shell_size', shell_size
        write(*,*) 'shell_ratio', shell_ratio
      end if
      if (iflag_debug .gt. 0) write(*,*) 'plane_size', plane_size(1:2)
!
      end subroutine set_control_draw_zplane
!
!-----------------------------------------------------------------------
!
      subroutine set_control_draw_map
!
      use m_spheric_constants
      use m_ctl_data_plot_pg
      use skip_comment_f
!
!
      id_shell_mode_pg = iflag_MESH_same
      if(pg_grid_type_ctl%iflag .gt. 0) then
        if(cmp_no_case(pg_grid_type_ctl%charavalue, 'No_pole'   ))      &
     &          id_shell_mode_pg = iflag_MESH_same
        if(cmp_no_case(pg_grid_type_ctl%charavalue, 'With_pole' ))      &
     &          id_shell_mode_pg = iflag_MESH_w_pole
        if(cmp_no_case(pg_grid_type_ctl%charavalue,'With_center'))      &
     &          id_shell_mode_pg = iflag_MESH_w_center
      end if
!
!
      if(radial_ID_ctl%iflag .gt. 0) then
        id_radial = radial_ID_ctl%intvalue
      else
        id_radial = 1
      end if
!
       r_sph = 20.0d0 / 13.0d0
!
      if (iflag_debug .gt. 0) then
        write(*,*) 'id_shell_mode_pg', id_shell_mode_pg
        write(*,*) 'id_radial', id_radial
      end if
!
      end subroutine set_control_draw_map
!
!-----------------------------------------------------------------------
!
      subroutine set_psffield_id_4_plot_pg(psf_phys)
!
      use t_phys_data
!
      type(phys_data), intent(in) :: psf_phys
      integer(kind = kint) :: iw, id
!
!
      do iw = 1, ntot_plot_pg
        do id = 1, psf_phys%num_phys
          if (field_name_4_plot(iw) .eq. psf_phys%phys_name(id)) then
            id_field_4_plot(iw) = id
            exit
          end if
        end do
      end do
!
      if(iflag_debug .gt. 0) then
        write(*,*) 'iw, id_field_4_plot(iw), field_name_4_plot'
        do iw = 1, ntot_plot_pg
          write(*,*) iw, id_field_4_plot(iw), id_comp_4_plot(iw),       &
     &               trim(field_name_4_plot(iw))
        end do
      end if
!
      end subroutine set_psffield_id_4_plot_pg
!
!  ---------------------------------------------------------------------
!
      end module set_control_draw_pg
