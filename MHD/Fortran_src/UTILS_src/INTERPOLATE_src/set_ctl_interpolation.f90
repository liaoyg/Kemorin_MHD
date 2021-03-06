!
!      module set_ctl_interpolation
!
!     Written by H. Matsui on July, 2006
!
!!      subroutine set_ctl_params_interpolation
!
      module set_ctl_interpolation
!
      use m_precision
!
      use calypso_mpi
      use m_error_IDs
!
      use t_file_IO_parameter
      use t_IO_step_parameter
!
      implicit none
!
!   --------------------------------------------------------------------
!
      contains
!
!   --------------------------------------------------------------------
!
      subroutine set_ctl_params_interpolation
!
      use m_machine_parameter
      use m_2nd_pallalel_vector
      use m_ctl_params_4_gen_table
      use m_ctl_data_gen_table
      use m_file_format_switch
      use m_field_file_format
      use itp_table_IO_select_4_zlib
      use set_control_platform_data
      use parallel_ucd_IO_select
!
!
      call turn_off_debug_flag_by_ctl(my_rank, src_plt)
      call set_control_smp_def(my_rank, src_plt)
!
      call set_control_mesh_def(src_plt, itp_org_mesh_file)
!
      if (dst_plt%mesh_file_prefix%iflag .ne. 0) then
        itp_dest_mesh_file%file_prefix                                  &
     &     = dst_plt%mesh_file_prefix%charavalue
      end if
!
      if (table_head_ctl%iflag .ne. 0) then
        table_file_head = table_head_ctl%charavalue
      end if
!
      if (single_itp_tbl_head_ctl%iflag .ne. 0) then
        sgl_table_file_head = single_itp_tbl_head_ctl%charavalue
      end if
!
      if (iflag_debug.eq.1) then
        write(*,*) 'np_smp', np_smp, np_smp
        write(*,*) 'org_mesh_head: ',                                   &
     &            trim(itp_org_mesh_file%file_prefix)
        write(*,*) 'dest_mesh_head: ',                                  &
     &            trim(itp_dest_mesh_file%file_prefix)
        write(*,*) 'table_file_head: ', trim(table_file_head)
      end if
!
      if (itp_node_head_ctl%iflag .ne. 0) then
        itp_node_file_head = itp_node_head_ctl%charavalue
        if (iflag_debug.eq.1)                                           &
     &   write(*,*) 'itp_node_file_head: ', trim(itp_node_file_head)
      end if
!
!
      call set_parallel_file_ctl_params(def_org_rst_prefix,             &
     &    src_plt%restart_file_prefix, src_plt%restart_file_fmt_ctl,    &
     &    org_fst_IO)
      call set_parallel_file_ctl_params(def_itp_rst_prefix,             &
     &    dst_plt%restart_file_prefix, dst_plt%restart_file_fmt_ctl,    &
     &    itp_fst_IO)
!
      call set_merged_ucd_file_define(src_plt, org_ucd_IO)
      call set_merged_ucd_file_ctl(itp_udt_file_head,                   &
     &    dst_plt%field_file_prefix, dst_plt%field_file_fmt_ctl,        &
     &    itp_ucd_IO)
!
      ndomain_org = 1
      if (src_plt%ndomain_ctl%iflag .gt. 0) then
        ndomain_org = src_plt%ndomain_ctl%intvalue
      end if
!
      nprocs_2nd = ndomain_org
      if (iflag_debug.eq.1)   write(*,*) 'ndomain_org', nprocs_2nd
!
      if (dst_plt%ndomain_ctl%iflag .gt. 0) then
        ndomain_dest = dst_plt%ndomain_ctl%intvalue
      else
        ndomain_dest = 1
      end if
!
      call choose_file_format                                           &
     &   (dst_plt%mesh_file_fmt_ctl, itp_dest_mesh_file%iflag_format)
!
      call choose_file_format                                           &
     &   (fmt_itp_table_file_ctl, ifmt_itp_table_file)
!
      if (nprocs .ne. max(ndomain_org,ndomain_dest) ) then
        write(e_message,*)                                              &
     &     'Num. of rank is larger num. of orgin or destinate dom.'
        call  calypso_MPI_abort(ierr_P_MPI, e_message)
      end if
!
      end subroutine set_ctl_params_interpolation
!
!  ---------------------------------------------------------------------
!
      end module set_ctl_interpolation
