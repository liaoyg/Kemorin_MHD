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
!!      subroutine SPH_const_initial_field
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
      use m_ctl_data_sph_MHD
      use m_SPH_MHD_model_data
      use m_MHD_step_parameter
      use t_MHD_file_parameter
      use t_SPH_mesh_field_data
      use t_field_data_IO
!
      implicit none
!
      type(field_IO), save, private :: rst_IO1
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
      use t_ctl_data_sph_MHD_psf
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
      call start_elapsed_time(1)
      call start_elapsed_time(4)
      if (iflag_debug.eq.1) write(*,*) 'read_control_4_sph_MHD_noviz'
      call read_control_4_sph_MHD_noviz(MHD_ctl_name, DNS_MHD_ctl1)
!
      if (iflag_debug.eq.1) write(*,*) 'input_control_4_SPH_make_init'
      call input_control_4_SPH_make_init                                &
     &   (MHD_files1, DNS_MHD_ctl1, MHD_step1, SPH_model1,              &
     &    SPH_WK1%trns_WK, SPH_WK1%monitor, SPH_MHD1, FEM_d1)
      call copy_delta_t(MHD_step1%init_d, MHD_step1%time_d)
      call end_elapsed_time(4)
!
!        Initialize spherical transform dynamo
!
      call start_elapsed_time(2)
      if(iflag_debug .gt. 0) write(*,*) 'SPH_const_initial_field'
      call SPH_const_initial_field(SPH_model1, SPH_MHD1)
!
      call end_elapsed_time(2)
      call reset_elapse_4_init_sph_mhd
!
      end subroutine initialize_const_sph_initial
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine SPH_const_initial_field(SPH_model, SPH_MHD)
!
      use set_control_sph_mhd
      use set_sph_phys_address
      use parallel_load_data_4_sph
      use const_sph_initial_spectr
      use set_reference_sph_mhd
      use set_bc_sph_mhd
      use adjust_reference_fields
      use material_property
      use sph_transforms_4_MHD
      use init_radial_infos_sph_mhd
      use const_radial_mat_4_sph
      use set_initial_sph_dynamo
      use check_dependency_for_MHD
      use input_control_sph_MHD
!
      type(SPH_MHD_model_data), intent(inout) :: SPH_model
      type(SPH_mesh_field_data), intent(inout) :: SPH_MHD
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
      if (iflag_debug.gt.0) write(*,*) 'init_r_infos_make_sph_initial'
      call init_r_infos_make_sph_initial(SPH_model, SPH_MHD)
!
! ---------------------------------
!
      if(iflag_debug.gt.0) write(*,*)' sph_initial_spectrum'
      call sph_initial_spectrum(MHD_files1%fst_file_IO,                 &
     &    SPH_model%sph_MHD_bc, SPH_MHD, MHD_step1%rst_step, rst_IO1)
!
      end subroutine SPH_const_initial_field
!
! ----------------------------------------------------------------------
!
      end module SPH_analyzer_const_initial
