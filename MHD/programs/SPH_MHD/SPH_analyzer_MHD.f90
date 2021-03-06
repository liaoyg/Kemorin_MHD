!>@file   SPH_analyzer_MHD
!!@brief  module SPH_analyzer_MHD
!!
!!@author H. Matsui
!!@date    programmed by H.Matsui in Oct., 2009
!
!>@brief Evolution loop for spherical MHD
!!
!!@verbatim
!!      subroutine SPH_initialize_MHD(MHD_files, SPH_model,             &
!!     &          iphys, MHD_step, sph_fst_IO, SPH_MHD)
!!        type(MHD_file_IO_params), intent(in) :: MHD_files
!!        type(phys_address), intent(in) :: iphys
!!        type(SPH_MHD_model_data), intent(inout) :: SPH_model
!!        type(MHD_step_param), intent(inout) :: MHD_step
!!        type(SPH_mesh_field_data), intent(inout) :: SPH_MHD
!!        type(work_SPH_MHD), intent(inout) :: SPH_WK
!!        type(field_IO), intent(inout) :: sph_fst_IO
!!      subroutine SPH_analyze_MHD(i_step, MHD_files, SPH_model,        &
!!     &          iflag_finish, MHD_step, sph_fst_IO, SPH_MHD, SPH_WK)
!!        type(MHD_file_IO_params), intent(in) :: MHD_files
!!        type(MHD_step_param), intent(inout) :: MHD_step
!!        type(SPH_mesh_field_data), intent(inout) :: SPH_MHD
!!        type(work_SPH_MHD), intent(inout) :: SPH_WK
!!        type(field_IO), intent(inout) :: sph_fst_IO
!!@endverbatim
!
      module SPH_analyzer_MHD
!
      use m_precision
      use m_constants
      use m_MHD_step_parameter
      use t_phys_address
      use t_control_parameter
      use t_MHD_step_parameter
      use t_MHD_file_parameter
      use t_SPH_mesh_field_data
      use t_boundary_data_sph_MHD
      use t_work_SPH_MHD
      use t_field_data_IO
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine SPH_initialize_MHD(MHD_files, SPH_model,               &
     &          iphys, MHD_step, sph_fst_IO, SPH_MHD, SPH_WK)
!
      use calypso_mpi
      use m_machine_parameter
!
      use t_sph_boundary_input_data
!
      use set_control_sph_mhd
      use set_sph_phys_address
      use const_fdm_coefs
      use set_initial_sph_dynamo
      use adjust_reference_fields
      use set_bc_sph_mhd
      use adjust_reference_fields
      use material_property
      use init_radial_infos_sph_mhd
      use const_radial_mat_4_sph
      use cal_sol_sph_MHD_crank
      use cal_nonlinear
      use init_sphrical_transform_MHD
      use check_dependency_for_MHD
      use input_control_sph_MHD
!
      use m_work_time
!
      type(MHD_file_IO_params), intent(in) :: MHD_files
      type(phys_address), intent(in) :: iphys
!
      type(SPH_MHD_model_data), intent(inout) :: SPH_model
      type(MHD_step_param), intent(inout) :: MHD_step
      type(SPH_mesh_field_data), intent(inout) :: SPH_MHD
      type(work_SPH_MHD), intent(inout) :: SPH_WK
      type(field_IO), intent(inout) :: sph_fst_IO
!
!
!   Allocate spectr field data
!
      call set_sph_MHD_sprctr_data                                      &
     &   (SPH_MHD%sph%sph_rj, SPH_model%MHD_prop,                       &
     &    SPH_MHD%ipol, SPH_MHD%idpdr, SPH_MHD%itor, SPH_MHD%fld)
!
! ---------------------------------
!
      if (iflag_debug.gt.0) write(*,*) 'init_r_infos_sph_mhd_evo'
      call init_r_infos_sph_mhd_evo(SPH_model, SPH_WK%r_2nd, SPH_MHD)
!
! ---------------------------------
!
      if (iflag_debug.gt.0) write(*,*) 'init_sph_transform_MHD'
      call init_sph_transform_MHD                                       &
     &   (SPH_model, iphys, SPH_WK%trans_p, SPH_WK%trns_WK, SPH_MHD)
!
!  -------------------------------
!
      if(iflag_debug.gt.0) write(*,*)' sph_initial_data_control'
      call sph_initial_data_control                                     &
     &   (MHD_files, SPH_model, SPH_MHD, MHD_step, sph_fst_IO)
      MHD_step%iflag_initial_step = 0
!
      if(iflag_debug.gt.0) write(*,*)' sync_temp_by_per_temp_sph'
      call sync_temp_by_per_temp_sph(SPH_model,                         &
     &    SPH_MHD%sph%sph_rj, SPH_MHD%ipol, SPH_MHD%idpdr, SPH_MHD%fld)
!
!  -------------------------------
!
      if(iflag_debug.gt.0) write(*,*)' const_radial_mat_sph_mhd'
      call const_radial_mat_sph_mhd                                     &
     &   (MHD_step%time_d%dt, SPH_model%MHD_prop, SPH_model%sph_MHD_bc, &
     &    SPH_MHD%sph%sph_rj, SPH_WK%r_2nd, SPH_WK%trans_p%leg,         &
     &    SPH_WK%MHD_mats)
!*
!* obtain linear terms for starting
!*
      if(iflag_debug .gt. 0) write(*,*) 'set_sph_field_to_start'
      call set_sph_field_to_start                                       &
     &   (SPH_MHD%sph%sph_rj, SPH_WK%r_2nd, SPH_model%MHD_prop,         &
     &    SPH_model%sph_MHD_bc, SPH_WK%trans_p%leg,                     &
     &    SPH_MHD%ipol, SPH_MHD%itor, SPH_MHD%fld)
!
!* obtain nonlinear terms for starting
!*
      if(iflag_debug .gt. 0) write(*,*) 'first nonlinear'
      call nonlinear(SPH_WK%r_2nd, SPH_model,                           &
     &    SPH_WK%trans_p, SPH_WK%trns_WK, SPH_MHD)
!
!* -----  Open Volume integration data files -----------------
!*
      if(iflag_debug .gt. 0) write(*,*) 'open_sph_vol_rms_file_mhd'
      call start_elapsed_time(4)
      call open_sph_vol_rms_file_mhd                                    &
     &   (SPH_MHD%sph, SPH_MHD%ipol, SPH_MHD%fld, SPH_WK%monitor)
      call end_elapsed_time(4)
!
      end subroutine SPH_initialize_MHD
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine SPH_analyze_MHD(i_step, MHD_files, SPH_model,          &
     &          iflag_finish, MHD_step, sph_fst_IO, SPH_MHD, SPH_WK)
!
      use m_work_time
!
      use cal_momentum_eq_explicit
      use cal_sol_sph_MHD_crank
      use cal_nonlinear
      use adjust_reference_fields
      use lead_fields_4_sph_mhd
      use sph_mhd_rst_IO_control
      use output_viz_file_control
!
      integer(kind = kint), intent(in) :: i_step
      type(MHD_file_IO_params), intent(in) :: MHD_files
      type(SPH_MHD_model_data), intent(in) :: SPH_model
!
      integer(kind = kint), intent(inout) :: iflag_finish
      type(MHD_step_param), intent(inout) :: MHD_step
      type(SPH_mesh_field_data), intent(inout) :: SPH_MHD
      type(work_SPH_MHD), intent(inout) :: SPH_WK
      type(field_IO), intent(inout) :: sph_fst_IO
!
      integer(kind = kint) :: iflag
      real(kind = kreal) :: total_max
!
!*  ----------  add time evolution -----------------
!*
!
      call start_elapsed_time(5)
      call start_elapsed_time(6)
      if(iflag_debug.gt.0) write(*,*) 'sel_explicit_sph'
      call sel_explicit_sph(i_step, MHD_step%time_d%dt,                 &
     &    SPH_model%MHD_prop, SPH_model%sph_MHD_bc, SPH_MHD%sph%sph_rj, &
     &    SPH_MHD%ipol, SPH_MHD%itor, SPH_MHD%fld)
!*
!*  ----------  time evolution by inplicit method ----------
!*
      call start_elapsed_time(7)
      call s_cal_sol_sph_MHD_crank                                      &
     &   (MHD_step%time_d%dt, SPH_MHD%sph%sph_rj, SPH_WK%r_2nd,         &
     &    SPH_model%MHD_prop, SPH_model%sph_MHD_bc, SPH_WK%trans_p%leg, &
     &    SPH_MHD%ipol, SPH_MHD%idpdr, SPH_MHD%itor, SPH_WK%MHD_mats,   &
     &    SPH_MHD%fld)
      call end_elapsed_time(7)
      call end_elapsed_time(6)
!*
!*  ----------------lead nonlinear term ... ----------
!*
      call start_elapsed_time(8)
      call nonlinear(SPH_WK%r_2nd, SPH_model,                           &
     &    SPH_WK%trans_p, SPH_WK%trns_WK, SPH_MHD)
      call end_elapsed_time(8)
      call end_elapsed_time(5)
!
!* ----  Update fields after time evolution ------------------------=
!*
      call start_elapsed_time(9)
      if(iflag_debug.gt.0) write(*,*) 'trans_per_temp_to_temp_sph'
      call trans_per_temp_to_temp_sph(SPH_model,                        &
     &    SPH_MHD%sph%sph_rj, SPH_MHD%ipol, SPH_MHD%idpdr, SPH_MHD%fld)
!*
      iflag = lead_field_data_flag(i_step, MHD_step)
      if(iflag .eq. 0) then
        if(iflag_debug.gt.0) write(*,*) 's_lead_fields_4_sph_mhd'
        call s_lead_fields_4_sph_mhd                                    &
     &     (SPH_MHD%sph, SPH_MHD%comms, SPH_WK%r_2nd,                   &
     &      SPH_model%MHD_prop, SPH_model%sph_MHD_bc, SPH_WK%trans_p,   &
     &      SPH_MHD%ipol, SPH_WK%MHD_mats, SPH_WK%trns_WK, SPH_MHD%fld)
      end if
      call end_elapsed_time(9)
!
!*  -----------  output restart data --------------
!*
      call start_elapsed_time(4)
      call start_elapsed_time(10)
      iflag = set_IO_step_flag(MHD_step%time_d%i_time_step,             &
     &                         MHD_step%rst_step)
      if(iflag .eq. 0) then
        if(iflag_debug.gt.0) write(*,*) 'output_sph_restart_control'
        call output_sph_restart_control(MHD_files%fst_file_IO,          &
     &     MHD_step%time_d, SPH_MHD%fld, MHD_step%rst_step, sph_fst_IO)
      end if
!
      total_time = MPI_WTIME() - total_start
      call MPI_allREDUCE (total_time, total_max, ione, CALYPSO_REAL,    &
     &    MPI_MAX, CALYPSO_COMM, ierr_MPI)
      if      (MHD_step%finish_d%i_end_step .eq. -1                     &
     &   .and. total_max .gt. MHD_step%finish_d%elapsed_time) then
        MHD_step%rst_step%istep_file = MHD_step%finish_d%i_end_step
        iflag_finish = 1
        call output_sph_restart_control(MHD_files%fst_file_IO,          &
     &     MHD_step%time_d, SPH_MHD%fld, MHD_step%rst_step, sph_fst_IO)
      end if
      call end_elapsed_time(10)
!
!*  -----------  lead energy data --------------
!*
      call start_elapsed_time(11)
      iflag = output_IO_flag(i_step, MHD_step%rms_step)
      if(iflag .eq. 0) then
        if(iflag_debug.gt.0)  write(*,*) 'output_rms_sph_mhd_control'
        call output_rms_sph_mhd_control(MHD_step%time_d, SPH_MHD,       &
     &      SPH_model%sph_MHD_bc, SPH_WK%trans_p%leg, SPH_WK%monitor)
      end if
      call end_elapsed_time(11)
!
      if(iflag_debug.gt.0) write(*,*) 'sync_temp_by_per_temp_sph'
      call sync_temp_by_per_temp_sph(SPH_model,                         &
     &    SPH_MHD%sph%sph_rj, SPH_MHD%ipol, SPH_MHD%idpdr, SPH_MHD%fld)
!
      if(i_step .ge. MHD_step%finish_d%i_end_step                       &
     &    .and. MHD_step%finish_d%i_end_step .gt. 0) then
        iflag_finish = 1
      end if
      call end_elapsed_time(4)
!
      end subroutine SPH_analyze_MHD
!
! ----------------------------------------------------------------------
!
!      subroutine SPH_finalize_MHD
!
!      end subroutine SPH_finalize_MHD
!
! ----------------------------------------------------------------------
!
      end module SPH_analyzer_MHD
