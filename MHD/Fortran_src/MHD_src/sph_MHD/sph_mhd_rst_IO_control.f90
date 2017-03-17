!>@file   sph_mhd_rst_IO_control.f90
!!@brief  module sph_mhd_rst_IO_control
!!
!!@author H. Matsui
!!@date Programmed in 2009
!!@n    Modified in June, 2015
!
!>@brief  I/O routines for restart data
!!
!!@verbatim
!!      subroutine set_rst_file_by_orignal_mesh                         &
!!     &         (rj_file_param, rst_file_param)
!!        type(field_IO_params), intent(in) :: rj_file_param
!!        type(field_IO_params), intent(in) :: rst_file_param
!!
!!      subroutine init_output_sph_restart_file(rj_fld)
!!      subroutine output_sph_restart_control(rj_fld, rst_step)
!!      subroutine output_sph_rst_by_elaps(rj_fld)
!!        type(phys_data), intent(in) :: rj_fld
!!        type(IO_step_param), intent(inout) :: rst_step
!!
!!      subroutine read_alloc_sph_restart_data(rj_fld, rst_step)
!!        type(phys_data), intent(inout) :: rj_fld
!!        type(IO_step_param), intent(inout) :: rst_step
!!
!!      subroutine init_radial_sph_interpolation                        &
!!     &         (rj_file_param, sph_params, sph_rj)
!!      subroutine read_alloc_sph_rst_4_snap                            &
!!     &         (i_step, rj_file_param, sph_rj, ipol, rj_fld, rst_step)
!!        type(field_IO_params), intent(in) :: rj_file_param
!!        type(phys_data), intent(inout) :: rj_fld
!!        type(IO_step_param), intent(inout) :: rst_step
!!      subroutine output_spectr_4_snap                                 &
!!     &         (i_step, sph_file_param, rj_fld, ucd_step)
!!      subroutine read_alloc_sph_spectr                                &
!!     &         (i_step, rj_file_param, sph_file_param,                &
!!     &          sph_rj, ipol, rj_fld, ucd_step)
!!        type(field_IO_params), intent(in) :: rj_file_param
!!        type(field_IO_params), intent(in) :: sph_file_param
!!        type(phys_data), intent(in) :: rj_fld
!!      subroutine read_alloc_sph_rst_2_modify                          &
!!     &         (i_step, rj_file_param, rst_file_param,                &
!!     &          sph_rj, ipol, rj_fld, rst_step)
!!     &          rj_file_param, rst_file_param, sph_rj, ipol, rj_fld)
!!        type(field_IO_params), intent(in) :: rj_file_param
!!        type(field_IO_params), intent(in) :: rst_file_param
!!        type(phys_data), intent(inout) :: rj_fld
!!        type(IO_step_param), intent(inout) :: rst_step
!!@endverbatim
!!
!!@n @param i_step  time step
!
      module sph_mhd_rst_IO_control
!
      use m_precision
!
      use m_machine_parameter
      use calypso_mpi
      use m_t_step_parameter
      use m_file_format_switch
!
      use t_IO_step_parameter
      use t_phys_address
      use t_phys_data
      use t_file_IO_parameter
      use t_time_data_IO
!
      use field_IO_select
!
      implicit  none
!
!
      type(time_params_IO), save, private :: sph_time_IO
      type(field_IO), save :: sph_fst_IO
!
      type(field_IO), save, private :: sph_out_IO
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine set_rst_file_by_orignal_mesh                           &
     &         (rj_file_param, rst_file_param)
!
      type(field_IO_params), intent(in) :: rj_file_param
      type(field_IO_params), intent(in) :: rst_file_param
!
!
      if( (rj_file_param%iflag_IO*rst_file_param%iflag_IO) .gt. 0)      &
     &            sph_fst_IO%file_prefix = rst_file_param%file_prefix
!
      end subroutine set_rst_file_by_orignal_mesh
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine init_output_sph_restart_file(rj_fld)
!
      use m_initial_field_control
      use set_sph_restart_IO
!
      type(phys_data), intent(in) :: rj_fld
!
!
      time_d1%i_time_step =   i_step_init
      time_d1%time = time_init
!
      call set_sph_restart_num_to_IO(rj_fld, sph_fst_IO)
!
      end subroutine init_output_sph_restart_file
!
! -----------------------------------------------------------------------
!
      subroutine output_sph_restart_control(rj_fld, rst_step)
!
      use set_sph_restart_IO
      use copy_time_steps_4_restart
!
      type(phys_data), intent(in) :: rj_fld
      type(IO_step_param), intent(inout) :: rst_step
!
!
      if (output_IO_flag(time_d1%i_time_step,rst_step) .ne. 0) return
!
      rst_step%istep_file = time_d1%i_time_step/rst_step%increment
!
      call copy_time_steps_to_restart(sph_time_IO)
      call set_sph_restart_data_to_IO(rj_fld, sph_fst_IO)
!
      call sel_write_step_SPH_field_file(nprocs, my_rank,               &
     &    rst_step%istep_file, sph_time_IO, sph_fst_IO)
!
      end subroutine output_sph_restart_control
!
! -----------------------------------------------------------------------
!
      subroutine output_sph_rst_by_elaps(rj_fld)
!
      use set_sph_restart_IO
      use copy_time_steps_4_restart
!
      type(phys_data), intent(in) :: rj_fld
!
      integer(kind = kint), parameter :: negaone = -1
!
!
      call copy_time_steps_to_restart(sph_time_IO)
      call set_sph_restart_data_to_IO(rj_fld, sph_fst_IO)

      call sel_write_step_SPH_field_file                                &
     &   (nprocs, my_rank, negaone, sph_time_IO, sph_fst_IO)
!
      end subroutine output_sph_rst_by_elaps
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine read_alloc_sph_restart_data(rj_fld, rst_step)
!
      use set_sph_restart_IO
      use copy_time_steps_4_restart
!
      type(phys_data), intent(inout) :: rj_fld
      type(IO_step_param), intent(inout) :: rst_step
!
!
      if (i_step_init .eq. -1) then
        call sel_read_alloc_step_SPH_file                               &
     &     (nprocs, my_rank, i_step_init, sph_time_IO, sph_fst_IO)
      else
        rst_step%istep_file = i_step_init / rst_step%increment
        call sel_read_alloc_step_SPH_file(nprocs, my_rank,              &
     &      rst_step%istep_file, sph_time_IO, sph_fst_IO)
      end if
!
      call set_sph_restart_from_IO(sph_fst_IO, rj_fld)
      call copy_init_time_from_restart(sph_time_IO)
!
      call dealloc_phys_data_IO(sph_fst_IO)
      call dealloc_phys_name_IO(sph_fst_IO)
!
      end subroutine read_alloc_sph_restart_data
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine init_radial_sph_interpolation                          &
     &         (rj_file_param, sph_params, sph_rj)
!
      use t_spheric_parameter
      use r_interpolate_sph_data
!
      type(field_IO_params), intent(in) :: rj_file_param
      type(sph_shell_parameters), intent(inout) :: sph_params
      type(sph_rj_grid), intent(inout) ::  sph_rj
!
!
      if(rj_file_param%iflag_IO .gt. 0) then
        if(iflag_debug .gt. 0) write(*,*) 'input_old_rj_sph_trans'
        call input_old_rj_sph_trans                                     &
     &     (rj_file_param, sph_params%l_truncation, sph_rj)
      end if
!
      call copy_cmb_icb_radial_point                                    &
     &   (sph_params%nlayer_ICB, sph_params%nlayer_CMB)
!
      end subroutine init_radial_sph_interpolation
!
! -----------------------------------------------------------------------
!
      subroutine read_alloc_sph_rst_4_snap                              &
     &         (i_step, rj_file_param, sph_rj, ipol, rj_fld, rst_step)
!
      use t_spheric_rj_data
      use set_sph_restart_IO
      use r_interpolate_sph_data
      use copy_time_steps_4_restart
!
      integer(kind = kint), intent(in) :: i_step
      type(field_IO_params), intent(in) :: rj_file_param
      type(sph_rj_grid), intent(in) ::  sph_rj
      type(phys_address), intent(in) :: ipol
      type(phys_data), intent(inout) :: rj_fld
      type(IO_step_param), intent(inout) :: rst_step
!
!
      rst_step%istep_file = i_step / rst_step%increment
      call sel_read_alloc_step_SPH_file(nprocs, my_rank,                &
     &    rst_step%istep_file, sph_time_IO, sph_fst_IO)
!
      call copy_init_time_from_restart(sph_time_IO)
      time_d1%time = time_init
!
      if(rj_file_param%iflag_IO .eq. 0) then
        if (iflag_debug.gt.0) write(*,*) 'set_sph_restart_from_IO'
        call set_sph_restart_from_IO(sph_fst_IO, rj_fld)
      else
        if (iflag_debug.gt.0)                                           &
     &            write(*,*) 'r_interpolate_sph_rst_from_IO'
        call r_interpolate_sph_rst_from_IO                              &
     &     (sph_fst_IO, sph_rj, ipol, rj_fld)
      end if
!
      call dealloc_phys_data_IO(sph_fst_IO)
      call dealloc_phys_name_IO(sph_fst_IO)
!
      end subroutine read_alloc_sph_rst_4_snap
!
! -----------------------------------------------------------------------
!
      subroutine output_spectr_4_snap                                   &
     &         (i_step, sph_file_param, rj_fld, ucd_step)
!
      use copy_rj_phys_data_4_IO
      use copy_time_steps_4_restart
      use const_global_element_ids
!
      type(phys_data), intent(in) :: rj_fld
      type(field_IO_params), intent(in) :: sph_file_param
      integer(kind = kint), intent(in) :: i_step
      type(IO_step_param), intent(inout) :: ucd_step
!
!
      if(sph_file_param%iflag_IO .eq. 0) return
      if(output_IO_flag(i_step,ucd_step) .ne. 0) return
!
      call set_field_file_fmt_prefix                                    &
     &   (sph_file_param%iflag_format, sph_file_param%file_prefix,      &
     &    sph_out_IO)
!
      ucd_step%istep_file = i_step / ucd_step%increment
      call copy_time_steps_to_restart(sph_time_IO)
      call copy_rj_phys_name_to_IO                                      &
     &   (rj_fld%num_phys_viz, rj_fld, sph_out_IO)
      call alloc_phys_data_IO(sph_out_IO)
      call copy_rj_phys_data_to_IO                                      &
     &   (rj_fld%num_phys_viz, rj_fld, sph_out_IO)
!
      call alloc_merged_field_stack(nprocs, sph_out_IO)
      call count_number_of_node_stack                                   &
     &   (sph_out_IO%nnod_IO, sph_out_IO%istack_numnod_IO)
!
      call sel_write_step_SPH_field_file(nprocs, my_rank,               &
     &    ucd_step%istep_file, sph_time_IO, sph_out_IO)
!
      call dealloc_merged_field_stack(sph_out_IO)
      call dealloc_phys_data_IO(sph_out_IO)
      call dealloc_phys_name_IO(sph_out_IO)
!
      end subroutine output_spectr_4_snap
!
! -----------------------------------------------------------------------
!
      subroutine read_alloc_sph_spectr                                  &
     &         (i_step, rj_file_param, sph_file_param,                  &
     &          sph_rj, ipol, rj_fld, ucd_step)
!
      use t_spheric_rj_data
      use copy_rj_phys_data_4_IO
      use r_interpolate_sph_data
      use copy_time_steps_4_restart
!
      integer(kind = kint), intent(in) :: i_step
      type(field_IO_params), intent(in) :: sph_file_param
      type(field_IO_params), intent(in) :: rj_file_param
      type(sph_rj_grid), intent(in) ::  sph_rj
      type(phys_address), intent(in) :: ipol
      type(phys_data), intent(inout) :: rj_fld
      type(IO_step_param), intent(inout) :: ucd_step
!
!
      ucd_step%istep_file = i_step / ucd_step%increment
      call set_field_file_fmt_prefix                                    &
     &   (sph_file_param%iflag_format, sph_file_param%file_prefix,      &
     &    sph_out_IO)
      call sel_read_alloc_step_SPH_file(nprocs, my_rank,                &
     &    ucd_step%istep_file, sph_time_IO, sph_out_IO)
!
      call copy_init_time_from_restart(sph_time_IO)
      time_d1%time = time_init
!
      if(rj_file_param%iflag_IO .eq. 0) then
        if (iflag_debug.gt.0) write(*,*) 'set_rj_phys_data_from_IO'
        call set_rj_phys_data_from_IO(sph_out_IO, rj_fld)
      else
        if (iflag_debug.gt.0) write(*,*)                                &
     &                        'r_interpolate_sph_fld_from_IO'
        call r_interpolate_sph_fld_from_IO                              &
     &    (sph_out_IO, sph_rj, ipol, rj_fld)
      end if
!
      call dealloc_merged_field_stack(sph_out_IO)
      call dealloc_phys_data_IO(sph_out_IO)
      call dealloc_phys_name_IO(sph_out_IO)
!
      end subroutine read_alloc_sph_spectr
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine read_alloc_sph_rst_2_modify                            &
     &         (i_step, rj_file_param, rst_file_param,                  &
     &          sph_rj, ipol, rj_fld, rst_step)
!
      use t_spheric_rj_data
!
      integer(kind = kint), intent(in) :: i_step
      type(field_IO_params), intent(in) :: rj_file_param
      type(field_IO_params), intent(in) :: rst_file_param
      type(sph_rj_grid), intent(in) ::  sph_rj
      type(phys_address), intent(in) :: ipol
      type(phys_data), intent(inout) :: rj_fld
      type(IO_step_param), intent(inout) :: rst_step
!
      integer(kind = kint) :: ifmt_rst_file_tmp
      character(len=kchara) :: restart_tmp_prefix
!
!
      ifmt_rst_file_tmp =  sph_fst_IO%iflag_file_fmt
      restart_tmp_prefix = sph_fst_IO%file_prefix
      sph_fst_IO%file_prefix = rst_file_param%file_prefix
      call set_field_file_fmt_prefix                                    &
     &   (ifmt_rst_file_tmp, rst_file_param%file_prefix, sph_fst_IO)
      call read_alloc_sph_rst_4_snap                                    &
     &   (i_step, rj_file_param, sph_rj, ipol, rj_fld, rst_step)
!
      call set_field_file_fmt_prefix                                    &
     &   (ifmt_rst_file_tmp, restart_tmp_prefix, sph_fst_IO)
!
      end subroutine read_alloc_sph_rst_2_modify
!
! -----------------------------------------------------------------------
!
      end module sph_mhd_rst_IO_control
