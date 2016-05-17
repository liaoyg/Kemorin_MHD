!>@file   SPH_analyzer_const_initial.f90
!!@brief  module SPH_analyzer_const_initial
!!
!!@author H. Matsui
!!@date   Programmed  H. Matsui in Apr., 2010
!
!>@brief  Main loop to generate initial field
!!@n       Define initial field at const_sph_initial_spectr.f90
!!
!!@verbatim
!!      subroutine initialize_const_sph_initial
!!@endverbatim
!
!
      module SPH_analyzer_const_initial
!
      use m_precision
      use m_constants
!
      use calypso_mpi
      use m_machine_parameter
      use m_work_time
      use m_control_parameter
      use m_t_int_parameter
      use m_t_step_parameter
!
      implicit none
!
      private :: SPH_const_initial_field
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine initialize_const_sph_initial
!
      use m_ctl_data_sph_MHD_noviz
      use m_spheric_parameter
      use m_sph_trans_comm_table
      use m_group_data_sph_specr
      use m_sph_spectr_data
      use set_control_sph_mhd
      use init_sph_MHD_elapsed_label
      use input_control_sph_MHD
!
!
      write(*,*) 'Simulation start: PE. ', my_rank
      total_start = MPI_WTIME()
      call set_sph_MHD_elapsed_label
!
!   Load parameter file
!
      call start_eleps_time(1)
      call start_eleps_time(4)
      if (iflag_debug.eq.1) write(*,*) 'read_control_4_sph_MHD_noviz'
      call read_control_4_sph_MHD_noviz
!
      if (iflag_debug.eq.1) write(*,*) 'input_control_4_SPH_make_init'
      call input_control_4_SPH_make_init                                &
     &   (sph_param1, sph1%sph_rtp, sph1%sph_rtm, sph1%sph_rlm, sph1%sph_rj,        &
     &    comm_rtp1, comm_rtm1, comm_rlm1, comm_rj1, bc_rtp_grp1,       &
     &    radial_rtp_grp1, theta_rtp_grp1, zonal_rtp_grp,               &
     &    radial_rj_grp1, sphere_rj_grp1, rj_fld1)
      call end_eleps_time(4)
!
!        Initialize spherical transform dynamo
!
      call start_eleps_time(2)
      if(iflag_debug .gt. 0) write(*,*) 'SPH_initialize_MHD'
      call SPH_const_initial_field
!
      call end_eleps_time(2)
      call reset_elapse_4_init_sph_mhd
!
      end subroutine initialize_const_sph_initial
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine SPH_const_initial_field
!
      use m_spheric_parameter
      use m_group_data_sph_specr
      use m_sph_spectr_data
      use m_sph_phys_address
      use m_physical_property
!
      use set_control_sph_mhd
      use parallel_load_data_4_sph
      use const_sph_initial_spectr
      use set_reference_sph_mhd
      use set_bc_sph_mhd
      use adjust_reference_fields
      use material_property
      use sph_transforms_4_MHD
      use set_radius_func
      use const_radial_mat_4_sph
      use set_initial_sph_dynamo
!
!
!   Allocate spectr field data
!
      call set_sph_sprctr_data_address(sph1%sph_rj, rj_fld1)
!
! ---------------------------------
!
      if (iflag_debug.gt.0) write(*,*) 'set_radius_rot_reft_dat_4_sph'
      call set_radius_rot_reft_dat_4_sph(depth_high_t, depth_low_t,     &
     &    high_temp, low_temp, angular, sph1%sph_rlm, sph1%sph_rj,      &
     &    radial_rj_grp1, sph_param1, rj_fld1)
!
      if(iflag_debug.gt.0) write(*,*) 's_set_bc_sph_mhd'
      call s_set_bc_sph_mhd(sph_param1, sph1%sph_rj, radial_rj_grp1,    &
     &    CTR_nod_grp_name, CTR_sf_grp_name)
      call init_reference_fields(sph_param1, sph1%sph_rj)
!
! ---------------------------------
!
      if(iflag_debug.gt.0) write(*,*)' sph_initial_spectrum'
      call sph_initial_spectrum(rj_fld1)
!
      end subroutine SPH_const_initial_field
!
! ----------------------------------------------------------------------
!
      end module SPH_analyzer_const_initial
