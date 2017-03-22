!analyzer_ene_sph_shell.f90
!      module analyzer_ene_sph_shell
!..................................................
!
!      modified by H. Matsui on Jan., 2008
!
!      subroutine initialize_ene_sph_shell
!      subroutine analyze_ene_sph_shell
!
      module analyzer_ene_sph_shell
!
      use m_precision
      use m_constants
      use m_machine_parameter
      use m_spheric_data_sph_spetr
      use m_schmidt_poly_on_rtm
      use calypso_mpi
!
      use t_time_data
      use cal_rms_fields_by_sph
      use field_IO_select
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine initialize_ene_sph_shell
!
      use m_ctl_data_4_sph_utils
      use m_ctl_params_sph_utils
      use parallel_load_data_4_sph
      use set_sph_phys_address
      use copy_rj_phys_data_4_IO
      use count_num_sph_smp
      use schmidt_poly_on_rtm_grid
      use cal_rms_fields_by_sph
!
!     --------------------- 
!
!     read controls
!
      if (iflag_debug.gt.0) write(*,*) 'read_control_data_sph_utils'
      call read_control_data_sph_utils
!
      if (iflag_debug.gt.0) write(*,*) 'set_ctl_data_4_sph_utils'
      call set_ctl_data_4_sph_utils(t_SHR, rj_fld_spec, pwr_spec)
!
!       set spectr grids
!
      if (iflag_debug.gt.0) write(*,*) 'load_para_sph_mesh'
      call load_para_sph_mesh(sph_mesh_spec%sph,                        &
     &    sph_mesh_spec%sph_comms, sph_mesh_spec%sph_grps)
!
!  ------  initialize spectr data
!
      if (iflag_debug.gt.0) write(*,*) 'sel_read_alloc_step_SPH_file'
      call set_field_file_fmt_prefix                                    &
     &  (sph_file_spec_p%iflag_format, sph_file_spec_p%file_prefix,     &
     &   sph_spec_IO)
      call sel_read_alloc_step_SPH_file(nprocs, my_rank,                &
     &    t_SHR%init_d%i_time_step, spec_time_IO, sph_spec_IO)
!
!  -------------------------------
!
      call set_sph_sprctr_data_address(sph_mesh_spec%sph%sph_rj,        &
     &    ipol_spec, idpdr_spec, itor_spec, rj_fld_spec)
!
      call init_rms_4_sph_spectr                                        &
     &   (sph_mesh_spec%sph%sph_params,                                 &
     &    sph_mesh_spec%sph%sph_rj, rj_fld_spec, pwr_spec, WK_pwr_spec)
!
      call alloc_schmidt_normalize                                      &
     &   (sph_mesh_spec%sph%sph_rlm%nidx_rlm(2),                        &
     &    sph_mesh_spec%sph%sph_rj%nidx_rj(2), leg_s)
      call copy_sph_normalization_2_rj                                  &
     &   (sph_mesh_spec%sph%sph_rj, leg_s%g_sph_rj)
!
      end subroutine initialize_ene_sph_shell
!
! ----------------------------------------------------------------------
!
      subroutine analyze_ene_sph_shell
!
      use m_ctl_params_sph_utils
      use m_schmidt_poly_on_rtm
      use copy_rj_phys_data_4_IO
      use output_sph_m_square_file
      use volume_average_4_sph
!
!
      integer(kind = kint) :: i_step
!
!
      call set_field_file_fmt_prefix                                    &
     &  (sph_file_spec_p%iflag_format, sph_file_spec_p%file_prefix,     &
     &   sph_spec_IO)
!
      do i_step = t_SHR%init_d%i_time_step, t_SHR%finish_d%i_end_step,  &
     &           t_SHR%ucd_step%increment
        t_SHR%time_d%i_time_step = i_step
!
!   Input spectr data
!
      call sel_read_step_SPH_field_file (nprocs, my_rank,               &
     &    t_SHR%time_d%i_time_step, spec_time_IO, sph_spec_IO)
!
        call set_rj_phys_data_from_IO(sph_spec_IO, rj_fld_spec)
        call copy_time_step_data(spec_time_IO, t_SHR%time_d)
!
!  evaluate energies
!
        if (iflag_debug.gt.0) write(*,*) 'cal_mean_squre_in_shell'
        call cal_mean_squre_in_shell                                    &
     &     (sph_mesh_spec%sph%sph_params%l_truncation,                  &
     &      sph_mesh_spec%sph%sph_rj, ipol_spec, rj_fld_spec,           &
     &      leg_s%g_sph_rj, pwr_spec, WK_pwr_spec)
!
        call write_sph_vol_ave_file                                     &
     &     (t_SHR%time_d, sph_mesh_spec%sph%sph_params,                 &
     &      sph_mesh_spec%sph%sph_rj, pwr_spec)
        call write_sph_vol_ms_file(my_rank, t_SHR%time_d,               &
     &     sph_mesh_spec%sph%sph_params, sph_mesh_spec%sph%sph_rj,      &
     &     pwr_spec)
        call write_sph_vol_ms_spectr_file(my_rank, t_SHR%time_d,        &
     &      sph_mesh_spec%sph%sph_params, sph_mesh_spec%sph%sph_rj,     &
     &      pwr_spec)
        call write_sph_layer_ms_file(my_rank, t_SHR%time_d,             &
     &      sph_mesh_spec%sph%sph_params, pwr_spec)
        call write_sph_layer_spectr_file(my_rank, t_SHR%time_d,         &
     &      sph_mesh_spec%sph%sph_params, pwr_spec)
      end do
!
      end subroutine analyze_ene_sph_shell
!
! ----------------------------------------------------------------------
!
      end module analyzer_ene_sph_shell
