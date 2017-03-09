!>@file   set_control_4_time_steps.f90
!!@brief  module set_control_4_time_steps
!!
!!@author H. Matsui and H. Okuda
!!@date Programmed by H. Okuda in 2000
!!@n    modified by H. Matsui in 2001
!!@n    modified by H. Matsui in Sep., 2006
!
!> @brief set parameters for time stepping
!!
!!@verbatim
!!      subroutine s_set_control_4_time_steps                           &
!!     &         (flex_p, SGS_par, MHD_step, mr_ctl, tctl)
!!        type(mhd_restart_control), intent(in) :: mr_ctl
!!        type(SGS_paremeters), intent(inout) :: SGS_par
!!        type(MHD_IO_step_param), intent(inout) :: MHD_step
!!        type(flexible_stepping_parameter), intent(inout) :: flex_p
!!        type(time_data_control), intent(inout) :: tctl
!!@endverbatim
!
      module set_control_4_time_steps
!
      use m_precision
!
      use calypso_mpi
      use m_error_IDs
      use m_machine_parameter
      use m_t_step_parameter
      use m_MHD_step_parameter
      use t_SGS_control_parameter
      use t_ctl_data_4_time_steps
      use t_VIZ_step_parameter
      use t_MHD_step_parameter
      use t_flex_delta_t_data
!
      implicit  none
!
      private :: set_flex_time_step_controls
      private :: set_fixed_time_step_controls
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine s_set_control_4_time_steps                             &
     &         (flex_p, SGS_par, MHD_step, mr_ctl, tctl)
!
      use t_ctl_data_mhd_evo_scheme
      use m_initial_field_control
      use cal_num_digits
      use skip_comment_f
!
      type(mhd_restart_control), intent(in) :: mr_ctl
      type(SGS_paremeters), intent(inout) :: SGS_par
      type(MHD_IO_step_param), intent(inout) :: MHD_step
      type(flexible_stepping_parameter), intent(inout) :: flex_p
      type(time_data_control), intent(inout) :: tctl
!
!
!  control for restert
!
      call set_initial_field_id(mr_ctl%restart_flag_ctl, tctl)
!
        flex_p%iflag_flexible_step = iflag_fixed_step
        if(tctl%flexible_step_ctl%iflag .gt. 0                          &
     &     .and. yes_flag(tctl%flexible_step_ctl%charavalue)) then
          flex_p%iflag_flexible_step = iflag_flex_step
        end if
!
        if (tctl%dt_ctl%iflag .eq. 0) then
          e_message = 'Set delta t'
          call calypso_MPI_abort(ierr_evo, e_message)
        else
          dt = tctl%dt_ctl%realvalue
          call cal_num_digit_real(dt, flex_p%dt_fact, flex_p%idt_digit)
        end if
!
        if(flex_p%iflag_flexible_step .eq. iflag_flex_step) then
          if (tctl%min_delta_t_ctl%iflag .eq. 0) then
            e_message = 'Set maximum delta t'
            call calypso_MPI_abort(ierr_evo, e_message)
          else
            flex_p%dt_min = tctl%min_delta_t_ctl%realvalue
          end if
!
          if (tctl%max_delta_t_ctl%iflag .eq. 0) then
            e_message = 'Set maximum delta t'
            call calypso_MPI_abort(ierr_evo, e_message)
          else
            flex_p%dt_max = tctl%max_delta_t_ctl%realvalue
          end if
!
          if (tctl%max_eps_to_shrink_ctl%iflag .eq. 0) then
            e_message = 'Set maximum error to shrink delta t'
            call calypso_MPI_abort(ierr_evo, e_message)
          else
            flex_p%max_eps_to_shrink                                    &
     &                = tctl%max_eps_to_shrink_ctl%realvalue
          end if
!
          if (tctl%min_eps_to_expand_ctl%iflag .eq. 0) then
            e_message = 'Set minimum error to expand delta t'
            call calypso_MPI_abort(ierr_evo, e_message)
          else
            flex_p%min_eps_to_expand                                    &
     &                = tctl%min_eps_to_expand_ctl%realvalue
          end if
!
          flex_p%istep_flex_to_max = izero
!
          if(dt .gt. zero) then
            flex_p%interval_flex_2_max = nint(flex_p%dt_max / dt)
          end if
        else
          flex_p%dt_max = dt
          flex_p%dt_min = dt
          flex_p%interval_flex_2_max = ione
          flex_p%istep_flex_to_max = izero
        end if
!
!   parameters for time evolution
!
      SGS_par%i_step_sgs_coefs = 1
!
      if(flex_p%iflag_flexible_step .eq. iflag_flex_step) then
        if (iflag_debug .ge. iflag_routine_msg)                         &
     &    write(*,*) 'set_flex_time_step_controls'
        call set_flex_time_step_controls                                &
     &     (flex_p, SGS_par, tctl, MHD_step)
      else
        if (iflag_debug .ge. iflag_routine_msg)                         &
     &    write(*,*) 'set_fixed_time_step_controls'
        call set_fixed_time_step_controls(SGS_par, tctl, MHD_step)
      end if
!
      if (i_step_number.eq.-1) then
        if (tctl%elapsed_time_ctl%iflag .eq. 0) then
          e_message                                                     &
     &      = 'Set elapsed time to finish (second)'
          call calypso_MPI_abort(ierr_evo, e_message)
        else
          elapsed_time  = tctl%elapsed_time_ctl%realvalue
        end if
      end if
!
      if (iflag_debug .ge. iflag_routine_msg) then
        write(*,*) 'dt', dt, flex_p%dt_fact, flex_p%idt_digit
        write(*,*) 'i_step_init ',i_step_init
        write(*,*) 'i_step_number ',i_step_number
        write(*,*) 'istep_rst_start ', istep_rst_start
        write(*,*) 'istep_rst_end ',  istep_rst_end
        write(*,*) 'elapsed_time ', elapsed_time
        write(*,*) 'i_step_check ', rms_step1%increment
        write(*,*) 'i_step_output_rst ', MHD_step%rst_step%increment
        write(*,*) 'i_step_output_ucd ', MHD_step%ucd_step%increment
      end if
!
      end subroutine s_set_control_4_time_steps
!
! -----------------------------------------------------------------------
!
      subroutine set_fixed_time_step_controls(SGS_par, tctl, MHD_step)
!
      use set_fixed_time_step_params
!
      type(SGS_paremeters), intent(inout) :: SGS_par
      type(time_data_control), intent(inout) :: tctl
      type(MHD_IO_step_param), intent(inout) :: MHD_step
!
      integer(kind = kint) :: ierr
!
!
      call s_set_fixed_time_step_params(tctl,                           &
     &    MHD_step%rst_step, MHD_step%ucd_step, MHD_step%viz_step,      &
     &    ierr, e_message)
      if(ierr .gt. 0) call calypso_MPI_abort(ierr, e_message)
!
      call set_output_step_4_fixed_step(ione, dt,                       &
     &    tctl%i_step_check_ctl, tctl%delta_t_check_ctl, rms_step1)
!
      if(SGS_par%model_p%iflag_dynamic .ne. id_SGS_DYNAMIC_OFF) then
        call set_output_step_4_fixed_step(ione, dt,                     &
     &      tctl%i_step_sgs_coefs_ctl, tctl%delta_t_sgs_coefs_ctl,      &
     &      SGS_par%sgs_step)
      end if
!
      call set_output_step_4_fixed_step                                 &
     &   (izero, dt, tctl%i_step_monitor_ctl, tctl%delta_t_monitor_ctl, &
     &    point_step1)
!
      call set_output_step_4_fixed_step(izero, dt,                      &
     &    tctl%i_step_boundary_ctl, tctl%delta_t_boundary_ctl,          &
     &    boundary_step1)
!
      end subroutine set_fixed_time_step_controls
!
! -----------------------------------------------------------------------
!
      subroutine set_flex_time_step_controls                            &
     &         (flex_p, SGS_par, tctl, MHD_step)
!
      type(SGS_paremeters), intent(inout) :: SGS_par
      type(flexible_stepping_parameter), intent(inout) :: flex_p
      type(time_data_control), intent(inout) :: tctl
      type(MHD_IO_step_param), intent(inout) :: MHD_step
!
!
      call set_flex_time_step_params                                    &
     &   (flex_p, SGS_par, tctl, MHD_step%rst_step, MHD_step%ucd_step,  &
     &    MHD_step%viz_step)
!
      call set_output_step_4_flex_step(ione, flex_p%dt_max,             &
     &    tctl%i_step_check_ctl, tctl%delta_t_check_ctl, rms_step1)
!
      if(SGS_par%model_p%iflag_dynamic .ne. id_SGS_DYNAMIC_OFF) then
        call set_output_step_4_flex_step(ione, flex_p%dt_max,           &
     &      tctl%i_step_sgs_coefs_ctl, tctl%delta_t_sgs_coefs_ctl,      &
     &      SGS_par%sgs_step)
      end if
!
      call set_output_step_4_flex_step(izero, flex_p%dt_max,            &
     &    tctl%i_step_monitor_ctl, tctl%delta_t_monitor_ctl,            &
     &    point_step1)
!
      call set_output_step_4_flex_step(izero, flex_p%dt_max,            &
     &    tctl%i_step_boundary_ctl, tctl%delta_t_boundary_ctl,          &
     &    boundary_step1)
!
      end subroutine set_flex_time_step_controls
!
! -----------------------------------------------------------------------
!
      subroutine set_flex_time_step_params                              &
     &         (flex_p, SGS_par, tctl, rst_step, ucd_step, viz_step)
!
      type(SGS_paremeters), intent(inout) :: SGS_par
      type(flexible_stepping_parameter), intent(inout) :: flex_p
      type(time_data_control), intent(inout) :: tctl
      type(IO_step_param), intent(inout) :: rst_step, ucd_step
      type(VIZ_step_params), intent(inout) :: viz_step
!
!
      istep_rst_start   = 0
      if (tctl%start_rst_step_ctl%iflag .gt. 0) then
        istep_rst_start   = tctl%start_rst_step_ctl%intvalue
      end if
!
      if (tctl%end_rst_step_ctl%iflag .eq. 0) then
        e_message = 'Set time to finish'
          call calypso_MPI_abort(ierr_evo, e_message)
      else
        istep_rst_end = tctl%end_rst_step_ctl%intvalue
      end if
!
!
      call set_output_step_4_flex_step(ione, flex_p%dt_max,             &
     &    tctl%i_step_rst_ctl, tctl%delta_t_rst_ctl, rst_step)
!
      call set_output_step_4_flex_step(ione, flex_p%dt_max,             &
     &   tctl%i_step_ucd_ctl, tctl%delta_t_field_ctl, ucd_step)
!
      call set_start_stop_by_restart(rst_step)
!
      call viz_flex_time_step_controls(tctl, dt, viz_step)
!
      if (istep_rst_end .eq. -1) then
        if (tctl%elapsed_time_ctl%iflag .eq. 0) then
          e_message                                                     &
     &      = 'Set elapsed time to finish (second)'
          call calypso_MPI_abort(ierr_evo, e_message)
        else
          elapsed_time  = tctl%elapsed_time_ctl%realvalue
        end if
      end if
!
      end subroutine set_flex_time_step_params
!
! -----------------------------------------------------------------------
!
      end module set_control_4_time_steps
