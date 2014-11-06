!>@file   mhd_restart_file_IO_control.f90
!!@brief  module mhd_restart_file_IO_control
!!
!!@author H. Matsui
!!@date   programmed by H.Matsui and H.Okuda
!!@n                           on July 2000 (ver 1.1)
!!@n      modified by H. Matsui on Sep., 2006
!!@n      modified by H. Matsui on Dec., 2007
!
!> @brief Call restart data IO routines
!!
!!@verbatim
!!      subroutine init_MHD_restart_output
!!      subroutine init_restart_4_snapshot
!!
!!      subroutine output_MHD_restart_file_ctl
!!      subroutine elspased_MHD_restart_ctl
!!
!!      subroutine input_MHD_restart_file_ctl
!!      subroutine input_restart_4_snapshot
!!@endverbatim
!
      module mhd_restart_file_IO_control
!
      use m_precision
!
      use calypso_mpi
      use m_t_step_parameter
!
      implicit  none
!
      private :: output_restart_files, input_restart_files
      private :: input_model_coef_file, output_model_coef_file
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine init_MHD_restart_output
!
      use m_field_data_IO
      use set_field_to_restart
!
!
      call count_field_num_to_restart
      call allocate_phys_data_name_IO
!
      call copy_field_name_to_restart
      call allocate_phys_data_IO
!
      end subroutine init_MHD_restart_output
!
! -----------------------------------------------------------------------
!
      subroutine init_restart_4_snapshot
!
      use m_geometry_parameter
!
      use field_IO_select
      use set_field_to_restart
!
      integer(kind = kint) :: index_rst
!
!
      index_rst = i_step_init / i_step_output_rst
      call sel_read_alloc_FEM_fld_head(my_rank, index_rst)
!
      numgrid_phys_IO = numnod
      call allocate_phys_data_IO
!
      end subroutine init_restart_4_snapshot
!
! ----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine output_MHD_restart_file_ctl
!
      integer(kind = kint) :: index_rst
!
!
      if ( mod(istep_max_dt,i_step_output_rst) .ne. 0 ) return
!
      index_rst = istep_max_dt / i_step_output_rst
      call output_restart_files(index_rst)
      call output_model_coef_file(index_rst)
!
      end subroutine output_MHD_restart_file_ctl
!
! -----------------------------------------------------------------------
!
      subroutine elspased_MHD_restart_ctl
!
      integer(kind = kint), parameter :: index_rst = -1
!
!
      call output_restart_files(index_rst)
      call output_model_coef_file(index_rst)
!
      end subroutine elspased_MHD_restart_ctl
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine input_MHD_restart_file_ctl
!
!
      call input_restart_files
      call input_model_coef_file
!
      end subroutine input_MHD_restart_file_ctl
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine output_restart_files(index_rst)
!
      use m_node_phys_address
      use field_IO_select
      use copy_time_steps_4_restart
      use set_field_to_restart
      use nod_phys_send_recv
!
      integer(kind = kint), intent(in) :: index_rst
!
!
      if(iphys%i_pre_mom .gt. 0) call vector_send_recv(iphys%i_pre_mom)
      if(iphys%i_pre_uxb .gt. 0) call vector_send_recv(iphys%i_pre_uxb)
      if(iphys%i_pre_heat .gt. 0)                                       &
     &                      call scalar_send_recv(iphys%i_pre_heat)
      if(iphys%i_pre_composit .gt. 0)                                   &
     &                      call scalar_send_recv(iphys%i_pre_composit)
!
      call copy_time_steps_to_restart
      call copy_field_data_to_restart
!
      call sel_write_step_FEM_field_file(my_rank, index_rst)
!
      end subroutine output_restart_files
!
! -----------------------------------------------------------------------
!
      subroutine input_restart_files
!
      use m_control_parameter
      use m_t_int_parameter
      use m_file_format_switch
!
      use field_IO_select
      use set_field_to_restart
      use copy_time_steps_4_restart
!
      integer(kind = kint) :: ierr
!
!
      call check_step_FEM_field_file(my_rank, istep_rst_start, ierr)
      if(ierr .gt. 0) call calypso_MPI_abort(ierr,'No restart file.')
!
      call sel_read_alloc_step_FEM_file(my_rank, istep_rst_start)
!
      call copy_field_data_from_restart
      call deallocate_phys_data_IO
      call deallocate_phys_data_name_IO
!
      if(iflag_flexible_step .eq. iflag_flex_step) then
        call copy_time_steps_from_restart
      else
        call copy_init_time_from_restart
      end if
!
      if(my_rank .eq. 0)  write(*,*) 'delta t ', dt, dt_fact, idt_digit
!
      end subroutine input_restart_files
!
! -----------------------------------------------------------------------
!
      subroutine input_restart_4_snapshot
!
      use m_node_phys_data
!
      use set_field_to_restart
      use field_IO_select
!
      integer(kind = kint) :: index_rst
!
!
      if ( mod(istep_max_dt,i_step_output_rst) .ne. 0) return
!
      index_rst = istep_max_dt / i_step_output_rst
      call sel_read_step_FEM_field_file(my_rank, index_rst)
!
      call copy_field_data_from_restart
      time =       time_init
      i_step_MHD = istep_max_dt
!
      end subroutine input_restart_4_snapshot
!
! ----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine output_model_coef_file(index_rst)
!
      use sgs_ini_model_coefs_IO
      use set_parallel_file_name
!
      integer(kind = kint), intent(in) :: index_rst
!
      character(len=kchara) :: fn_tmp
!
!
      if(iflag_dynamic_SGS .eq. id_SGS_DYNAMIC_OFF) return
!
      if(index_rst .lt. 0) then
        call add_elaps_postfix(rst_sgs_coef_head, fn_tmp)
      else
        call add_int_suffix(index_rst, rst_sgs_coef_head, fn_tmp)
      end if
      call add_dat_extension(fn_tmp, rst_sgs_coef_name)
!
      call output_ini_model_coefs
!
      end subroutine output_model_coef_file
!
! -----------------------------------------------------------------------
!
      subroutine input_model_coef_file
!
      use m_control_parameter
!
      use set_parallel_file_name
      use sgs_ini_model_coefs_IO
!
      character(len=kchara) :: fn_tmp
!
!
      if(iflag_dynamic_SGS .eq. id_SGS_DYNAMIC_OFF) return
      if(iflag_rst_sgs_coef_code .eq. 0) return
!
      if (i_step_init .eq. -1) then
        call add_elaps_postfix(rst_sgs_coef_head, fn_tmp)
      else
        call add_int_suffix(istep_rst_start, rst_sgs_coef_head, fn_tmp)
      end if
!
      call add_dat_extension(fn_tmp, rst_sgs_coef_name)
      call input_ini_model_coefs
!
      end subroutine input_model_coef_file
!
! -----------------------------------------------------------------------
!
      end module mhd_restart_file_IO_control
