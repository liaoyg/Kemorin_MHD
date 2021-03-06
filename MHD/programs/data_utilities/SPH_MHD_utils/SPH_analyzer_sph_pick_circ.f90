!SPH_analyzer_sph_pick_circ.f90
!     module SPH_analyzer_sph_pick_circ
!
!      Written by H. Matsui
!
!>@file   SPH_analyzer_d_bench.f90
!!        module SPH_analyzer_d_bench
!!
!!@author H. Matsui
!!@date   Programmed in 2012
!!@n      modified in 2013
!
!>@brief spherical harmonics part of 
!!       Initialzation and evolution loop to pick up data on circle
!!
!!@verbatim
!!      subroutine SPH_init_sph_pick_circle                             &
!!     &         (MHD_files, femmesh, iphys, SPH_model,                 &
!!     &          SPH_SGS, SPH_MHD, SPH_WK, cdat)
!!        type(MHD_file_IO_params), intent(in) :: MHD_files
!!        type(mesh_data), intent(in) :: femmesh
!!        type(phys_address), intent(in) :: iphys
!!        type(SPH_MHD_model_data), intent(inout) :: SPH_model
!!        type(SPH_SGS_structure), intent(inout) :: SPH_SGS
!!        type(SPH_mesh_field_data), intent(inout) :: SPH_MHD
!!        type(work_SPH_MHD), intent(inout) :: SPH_WK
!!        type(circle_fld_maker), intent(inout) :: cdat
!!      subroutine SPH_analyze_pick_circle(i_step, MHD_files, SPH_model,&
!!     &          SPH_SGS, SPH_MHD, SPH_WK, cdat)
!!        type(MHD_file_IO_params), intent(in) :: MHD_files
!!        type(SPH_MHD_model_data), intent(in) :: SPH_model
!!        type(SPH_SGS_structure), intent(inout) :: SPH_SGS
!!        type(SPH_mesh_field_data), intent(inout) :: SPH_MHD
!!        type(work_SPH_MHD), intent(inout) :: SPH_WK
!!      subroutine SPH_finalize_pick_circle
!!@endverbatim
!
      module SPH_analyzer_sph_pick_circ
!
      use m_precision
      use m_MHD_step_parameter
      use t_SPH_MHD_model_data
      use t_SPH_mesh_field_data
      use t_mesh_data
      use t_phys_address
      use t_MHD_file_parameter
      use t_SPH_SGS_structure
      use t_boundary_data_sph_MHD
      use t_field_on_circle
      use t_work_SPH_MHD
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine SPH_init_sph_pick_circle                               &
     &         (MHD_files, femmesh, iphys, SPH_model,                   &
     &          SPH_SGS, SPH_MHD, SPH_WK, cdat)
!
      use m_constants
      use m_array_for_send_recv
      use calypso_mpi
      use m_machine_parameter
!
      use t_sph_boundary_input_data
!
      use set_control_sph_mhd
      use set_sph_phys_address
      use const_fdm_coefs
      use adjust_reference_fields
      use set_bc_sph_mhd
      use adjust_reference_fields
      use material_property
      use sph_transforms_4_MHD
      use init_radial_infos_sph_mhd
      use const_radial_mat_4_sph
      use cal_rms_fields_by_sph
      use r_interpolate_sph_data
      use sph_mhd_rst_IO_control
      use cal_SGS_nonlinear
      use init_sph_trans_SGS_MHD
      use nod_phys_send_recv
      use sph_filtering
      use check_dependency_SGS_MHD
      use input_control_sph_MHD
!
      type(MHD_file_IO_params), intent(in) :: MHD_files
      type(mesh_data), intent(in) :: femmesh
      type(phys_address), intent(in) :: iphys
!
      type(SPH_MHD_model_data), intent(inout) :: SPH_model
      type(SPH_SGS_structure), intent(inout) :: SPH_SGS
      type(SPH_mesh_field_data), intent(inout) :: SPH_MHD
      type(circle_fld_maker), intent(inout) :: cdat
      type(work_SPH_MHD), intent(inout) :: SPH_WK
!
!   Allocate spectr field data
!
      call set_sph_SGS_MHD_spectr_data                                  &
     &   (SPH_SGS%SGS_par, SPH_model%MHD_prop, SPH_MHD)
!
      if (iflag_debug.gt.0 ) write(*,*) 'allocate_vector_for_solver'
      call allocate_vector_for_solver                                   &
     &   (isix, SPH_MHD%sph%sph_rtp%nnod_rtp)
!
      if(iflag_debug.gt.0) write(*,*)' init_nod_send_recv'
      call init_nod_send_recv(femmesh%mesh)
!
      if ( iflag_debug.gt.0 ) write(*,*) 'init_rms_4_sph_spectr_4_mhd'
      call init_rms_4_sph_spectr_4_mhd                                  &
     &   (SPH_MHD%sph, SPH_MHD%fld, SPH_WK%monitor)
!
! ---------------------------------
!
      if (iflag_debug.gt.0) write(*,*) 'init_r_infos_sph_mhd_evo'
      call init_r_infos_sph_mhd_evo(SPH_model, SPH_WK%r_2nd, SPH_MHD)
!
!  -------------------------------
!
      if (iflag_debug.gt.0) write(*,*) 'init_sph_transform_SGS_MHD'
      call init_sph_transform_SGS_MHD(SPH_SGS%SGS_par%model_p,          &
     &    SPH_model, iphys, SPH_WK%trans_p, SPH_WK%trns_WK, SPH_MHD)
!
! ---------------------------------
!
      call init_SGS_model_sph_mhd                                       &
     &   (SPH_SGS%SGS_par, SPH_MHD%sph, SPH_MHD%groups,                 &
     &    SPH_model%MHD_prop, SPH_SGS%dynamic)
!
!  -------------------------------
!
      if (iflag_debug.eq.1) write(*,*) 'const_radial_mat_sph_snap'
      call const_radial_mat_sph_snap                                    &
     &   (SPH_model%MHD_prop, SPH_model%sph_MHD_bc,                     &
     &    SPH_MHD%sph%sph_rj, SPH_WK%r_2nd, SPH_WK%trans_p%leg,         &
     &    SPH_WK%MHD_mats)
!
!     --------------------- 
!  set original spectr mesh data for extension of B
!
      call init_radial_sph_interpolation(MHD_files%org_rj_file_IO,      &
     &    SPH_MHD%sph%sph_params, SPH_MHD%sph%sph_rj)
!
!* -----  find mid-equator point -----------------
!
      call const_circle_point_global                                    &
     &   (SPH_MHD%sph%sph_params%l_truncation, SPH_MHD%sph%sph_rtp,     &
     &     SPH_MHD%sph%sph_rj, cdat)
!
      end subroutine SPH_init_sph_pick_circle
!
! ----------------------------------------------------------------------
!
      subroutine SPH_analyze_pick_circle(i_step, MHD_files, SPH_model,  &
     &          SPH_SGS, SPH_MHD, SPH_WK, cdat)
!
      use m_work_time
!
      use cal_SGS_nonlinear
      use cal_sol_sph_MHD_crank
      use adjust_reference_fields
      use lead_fields_SPH_SGS_MHD
      use sph_mhd_rst_IO_control
      use input_control_sph_MHD
      use output_viz_file_control
      use field_on_circle_IO
!
      integer(kind = kint), intent(in) :: i_step
      type(MHD_file_IO_params), intent(in) :: MHD_files
      type(SPH_MHD_model_data), intent(in) :: SPH_model
!
      type(SPH_SGS_structure), intent(inout) :: SPH_SGS
      type(SPH_mesh_field_data), intent(inout) :: SPH_MHD
      type(work_SPH_MHD), intent(inout) :: SPH_WK
      type(circle_fld_maker), intent(inout) :: cdat
!
      integer(kind = kint) :: iflag
!
!
      call read_alloc_sph_rst_4_snap(i_step,                            &
     &    MHD_files%org_rj_file_IO, MHD_files%fst_file_IO,              &
     &    SPH_MHD%sph, SPH_MHD%ipol, SPH_MHD%fld,                       &
     &    MHD_step1%rst_step, MHD_step1%init_d)
      call copy_time_data(MHD_step1%init_d, MHD_step1%time_d)
!
      call sync_temp_by_per_temp_sph(SPH_model,                         &
     &    SPH_MHD%sph%sph_rj, SPH_MHD%ipol, SPH_MHD%idpdr, SPH_MHD%fld)
!
!* obtain linear terms for starting
!*
      if(iflag_debug .gt. 0) write(*,*) 'set_sph_field_to_start'
      call set_sph_field_to_start(SPH_MHD%sph%sph_rj, SPH_WK%r_2nd,     &
     &    SPH_model%MHD_prop, SPH_model%sph_MHD_bc, SPH_WK%trans_p%leg, &
     &    SPH_MHD%ipol, SPH_MHD%itor, SPH_MHD%fld)
!
!*  ----------------lead nonlinear term ... ----------
!*
      call start_elapsed_time(8)
      call nonlinear_with_SGS                                           &
     &   (i_step, SPH_SGS%SGS_par, SPH_WK%r_2nd, SPH_model,             &
     &    SPH_WK%trans_p, SPH_WK%trns_WK, SPH_SGS%dynamic, SPH_MHD)
      call end_elapsed_time(8)
!
!* ----  Update fields after time evolution ------------------------=
!*
      call start_elapsed_time(9)
      if(iflag_debug.gt.0) write(*,*) 'trans_per_temp_to_temp_sph'
      call trans_per_temp_to_temp_sph(SPH_model,                        &
     &    SPH_MHD%sph%sph_rj, SPH_MHD%ipol, SPH_MHD%idpdr, SPH_MHD%fld)
!*
      iflag = lead_field_data_flag(i_step, MHD_step1)
      if(iflag .eq. 0) then
        if(iflag_debug.gt.0) write(*,*) 'lead_fields_4_SPH_SGS_MHD'
        call lead_fields_4_SPH_SGS_MHD                                  &
     &     (SPH_SGS%SGS_par, SPH_WK%r_2nd, SPH_model%MHD_prop,          &
     &      SPH_model%sph_MHD_bc, SPH_WK%trans_p, SPH_WK%MHD_mats,      &
     &      SPH_WK%trns_WK, SPH_SGS%dynamic, SPH_MHD)
      end if
      call end_elapsed_time(9)
!
!*  -----------  lead mid-equator field --------------
!*
      call start_elapsed_time(4)
      if(iflag_debug.gt.0)  write(*,*) 'sph_transfer_on_circle'
      call sph_transfer_on_circle                                       &
     &   (SPH_MHD%sph%sph_rj, SPH_MHD%fld, cdat)
      call write_field_data_on_circle                                   &
     &   (i_step, MHD_step1%time_d%time, cdat%circle, cdat%d_circle)
      call end_elapsed_time(4)
!
      end subroutine SPH_analyze_pick_circle
!
! ----------------------------------------------------------------------
!
!      subroutine SPH_finalize_pick_circle
!
!      end subroutine SPH_finalize_pick_circle
!
! ----------------------------------------------------------------------
!
      end module SPH_analyzer_sph_pick_circ
