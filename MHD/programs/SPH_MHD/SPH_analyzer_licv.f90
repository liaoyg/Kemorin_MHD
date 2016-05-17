!
!     module SPH_analyzer_licv
!
!      subroutine SPH_initialize_linear_conv
!      subroutine SPH_analyze_linear_conv(i_step, iflag_finish)
!
!      Written by H. Matsui
!
      module SPH_analyzer_licv
!
      use m_precision
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine SPH_initialize_linear_conv
!
      use calypso_mpi
      use m_constants
      use m_array_for_send_recv
      use m_machine_parameter
      use m_control_parameter
!
      use m_spheric_parameter
      use m_sph_trans_comm_table
      use m_group_data_sph_specr
      use m_sph_spectr_data
      use m_sph_phys_address
      use m_rms_4_sph_spectr
!
      use set_control_sph_mhd
      use set_initial_sph_dynamo
      use adjust_reference_fields
      use set_bc_sph_mhd
      use adjust_reference_fields
      use material_property
      use sph_transforms_4_MHD
      use set_radius_func
      use const_radial_mat_4_sph
      use sph_mhd_rms_IO
      use cal_sol_sph_MHD_crank
      use cal_nonlinear
!
      use m_work_time
!
!   Allocate spectr field data
!
      call set_sph_sprctr_data_address(sph_rj1, rj_fld1)
!
!
      if (iflag_debug.gt.0 ) write(*,*) 'allocate_vector_for_solver'
      call allocate_vector_for_solver(isix, sph_rtp1%nnod_rtp)
!
! ---------------------------------
!
      if (iflag_debug.gt.0) write(*,*) 'set_radius_rot_reft_dat_4_sph'
      call set_radius_rot_reft_dat_4_sph(depth_high_t, depth_low_t,     &
     &    high_temp, low_temp, angular, sph_rlm1, sph_rj1,              &
     &    radial_rj_grp1, sph_param1, rj_fld1)
!
      if (iflag_debug.gt.0) write(*,*) 'const_2nd_fdm_matrices'
      call const_2nd_fdm_matrices(sph_param1, sph_rj1)
!
! ---------------------------------
!
      if(iflag_debug.gt.0) write(*,*)' set_material_property'
      call set_material_property
!
!  -------------------------------
!
      if(iflag_debug.gt.0) write(*,*) 's_set_bc_sph_mhd'
      call s_set_bc_sph_mhd(sph_param1, sph_rj1, radial_rj_grp1,        &
     &    CTR_nod_grp_name, CTR_sf_grp_name)
      call init_reference_fields(sph_param1, sph_rj1)
!
!  -------------------------------
!
      if (iflag_debug.gt.0) write(*,*) 'init_sph_transform_MHD'
      call init_sph_transform_MHD                                       &
     &   (sph_param1, sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1,            &
     &    comm_rtp1, comm_rtm1, comm_rlm1, comm_rj1, rj_fld1)
!
! ---------------------------------
!
      if(iflag_debug.gt.0) write(*,*)' sph_initial_data_control'
      call sph_initial_data_control                                     &
     &   (sph_param1, sph_rj1, reftemp_rj, rj_fld1)
!
      if(iflag_debug.gt.0) write(*,*)' sync_temp_by_per_temp_sph'
      call sync_temp_by_per_temp_sph(reftemp_rj, sph_rj1, rj_fld1)
!
!  -------------------------------
!
      if(iflag_debug.gt.0) write(*,*)' const_radial_mat_sph_mhd'
      call const_radial_mat_sph_mhd(sph_rj1)
!*
!* obtain linear terms for starting
!*
      if(iflag_debug .gt. 0) write(*,*) 'set_sph_field_to_start'
      call set_sph_field_to_start(sph_rj1, rj_fld1)
!
!*  ----------------lead nonlinear term ... ----------
!*
      if(iflag_debug .gt. 0) write(*,*) 'first licv_exp'
      call licv_exp(reftemp_rj,                                         &
     &    sph_rlm1, sph_rj1, comm_rlm1, comm_rj1, rj_fld1)
!
!* -----  Open Volume integration data files -----------------
!*
      if(iflag_debug .gt. 0) write(*,*) 'open_sph_vol_rms_file_mhd'
      call open_sph_vol_rms_file_mhd(sph_param1, sph_rj1, rj_fld1)
!
      end subroutine SPH_initialize_linear_conv
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine SPH_analyze_linear_conv(i_step, iflag_finish)
!
      use m_work_time
      use m_t_step_parameter
      use m_spheric_parameter
      use m_sph_trans_comm_table
      use m_sph_spectr_data
!
      use cal_momentum_eq_explicit
      use cal_sol_sph_MHD_crank
      use adjust_reference_fields
      use lead_fields_4_sph_mhd
      use sph_mhd_rst_IO_control
      use sph_mhd_rms_IO
      use cal_nonlinear
!
      integer(kind = kint), intent(in) :: i_step
      integer(kind = kint), intent(inout) :: iflag_finish
!
!*  ----------  add time evolution -----------------
!*
      if(i_step .eq. 1) then
        if(iflag_debug.gt.0) write(*,*) 'cal_expricit_sph_euler'
        call cal_expricit_sph_euler(i_step, sph_rj1, rj_fld1)
      else
        if(iflag_debug.gt.0) write(*,*) 'cal_expricit_sph_adams'
        call cal_expricit_sph_adams(sph_rj1, rj_fld1)
      end if
!*
!*  ----------  time evolution by inplicit method ----------
!*
      call s_cal_sol_sph_MHD_crank(sph_rj1, rj_fld1)
!*
!* ----  Update fields after time evolution ------------------------=
!*
!
      call start_eleps_time(9)
      if(iflag_debug.gt.0) write(*,*) 'trans_per_temp_to_temp_sph'
      call trans_per_temp_to_temp_sph(reftemp_rj, sph_rj1, rj_fld1)
!*
      if(iflag_debug.gt.0) write(*,*) 's_lead_fields_4_sph_mhd'
      call s_lead_fields_4_sph_mhd                                      &
     &   (sph_param1, sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1,            &
     &    comm_rtp1, comm_rtm1, comm_rlm1, comm_rj1, rj_fld1)
      call end_eleps_time(9)
!
!*  ----------------lead nonlinear term ... ----------
!*
        call licv_exp(reftemp_rj,                                       &
     &      sph_rlm1, sph_rj1, comm_rlm1, comm_rj1, rj_fld1)
!
!*  -----------  output restart data --------------
!*
      call start_eleps_time(4)
      call start_eleps_time(10)
      if(iflag_debug.gt.0) write(*,*) 'output_sph_restart_control'
      call output_sph_restart_control(rj_fld1)
!
      total_time = MPI_WTIME() - total_start
      if      (istep_rst_end .eq. -1                                    &
    &   .and. total_time.gt.elapsed_time) then
        call output_sph_rst_by_elaps(rj_fld1)
        iflag_finish = 1
      end if
      call end_eleps_time(10)
!
!*  -----------  lead energy data --------------
!*
      call start_eleps_time(11)
      if(iflag_debug.gt.0)  write(*,*) 'output_rms_sph_mhd_control'
      call output_rms_sph_mhd_control(sph_param1, sph_rj1, rj_fld1)
      call end_eleps_time(11)
!
      if(iflag_debug.gt.0) write(*,*) 'sync_temp_by_per_temp_sph'
      call sync_temp_by_per_temp_sph(reftemp_rj, sph_rj1, rj_fld1)
      call end_eleps_time(4)
!
      if(i_step .ge. i_step_number .and. i_step_number.gt.0) then
        iflag_finish = 1
      end if
!
      end subroutine SPH_analyze_linear_conv
!
! ----------------------------------------------------------------------
!
!      subroutine SPH_finalize_licv
!
!
!      end subroutine SPH_finalize_licv
!
! ----------------------------------------------------------------------
!
      end module SPH_analyzer_licv
