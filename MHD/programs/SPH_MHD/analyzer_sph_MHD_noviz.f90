!>@file   analyzer_sph_MHD_noviz.f90
!!@brief  module analyzer_sph_MHD_noviz
!!
!!@author H. Matsui
!!@date   Programmed  H. Matsui in Apr., 2010
!
!>@brief  Main loop for MHD dynamo simulation
!!        without cross sectioning routines
!!
!!@verbatim
!!      subroutine initialize_sph_MHD_noviz
!!      subroutine evolution_sph_MHD_noviz
!!@endverbatim
!
      module analyzer_sph_MHD_noviz
!
      use m_precision
      use calypso_mpi
!
      use m_machine_parameter
      use m_work_time
      use m_SPH_MHD_model_data
      use m_MHD_step_parameter
!
      use FEM_analyzer_sph_MHD
      use SPH_analyzer_MHD
      use init_sph_MHD_elapsed_label
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine initialize_sph_MHD_noviz
!
      use t_ctl_data_sph_MHD_psf
      use m_ctl_data_sph_MHD
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
      if (iflag_debug.eq.1) write(*,*) 'input_control_SPH_MHD_psf'
      call input_control_SPH_MHD_psf                                    &
     &   (MHD_files1, DNS_MHD_ctl1, MHD_step1, SPH_model1,              &
     &    SPH_WK1%trns_WK, SPH_WK1%monitor, SPH_MHD1, FEM_d1)
      call copy_delta_t(MHD_step1%init_d, MHD_step1%time_d)
      call end_elapsed_time(4)
!
!        Initialize FEM mesh data for field data IO
      call start_elapsed_time(2)
      if(iflag_debug .gt. 0) write(*,*) 'FEM_initialize_sph_MHD'
      call FEM_initialize_sph_MHD                                       &
     &   (MHD_files1, MHD_step1, FEM_d1%geofem, FEM_d1%ele_mesh,        &
     &    FEM_d1%iphys, FEM_d1%field, MHD_IO1)
!
!        Initialize spherical transform dynamo
      if(iflag_debug .gt. 0) write(*,*) 'SPH_initialize_MHD'
      call SPH_initialize_MHD(MHD_files1, SPH_model1,                   &
     &    FEM_d1%iphys, MHD_step1, MHD_IO1%rst_IO, SPH_MHD1, SPH_WK1)
!
      call calypso_MPI_barrier
!
      call end_elapsed_time(2)
      call reset_elapse_4_init_sph_mhd
!
      end subroutine initialize_sph_MHD_noviz
!
! ----------------------------------------------------------------------
!
      subroutine evolution_sph_MHD_noviz
!
      use output_viz_file_control
!
      integer(kind = kint) :: visval, iflag_finish
      integer(kind = kint) :: iflag
!
!     ---------------------
!
      call start_elapsed_time(3)
!
!*  -----------  set initial step data --------------
!*
      call copy_time_step_data(MHD_step1%init_d, MHD_step1%time_d)
      iflag_finish = 0
!*
!*  -------  time evelution loop start -----------
!*
      do
        call evolve_time_data(MHD_step1%time_d)
!
!*  ----------  time evolution by spectral methood -----------------
!*
        if (iflag_debug.eq.1) write(*,*) 'SPH_analyze_MHD'
        call SPH_analyze_MHD(MHD_step1%time_d%i_time_step,              &
     &      MHD_files1, SPH_model1, iflag_finish,                       &
     &      MHD_step1, MHD_IO1%rst_IO, SPH_MHD1, SPH_WK1)
!*
!*  -----------  output field data --------------
!*
        call start_elapsed_time(4)
        iflag = lead_field_data_flag(MHD_step1%time_d%i_time_step,      &
     &                               MHD_step1)
        if(iflag .eq. 0) then
          if (iflag_debug.eq.1) write(*,*) 'SPH_to_FEM_bridge_MHD'
          call SPH_to_FEM_bridge_MHD                                    &
     &       (SPH_MHD1%sph%sph_params, SPH_MHD1%sph%sph_rtp,            &
     &        SPH_WK1%trns_WK, FEM_d1%geofem%mesh, FEM_d1%iphys,        &
     &        FEM_d1%field)
        end if
!
        if (iflag_debug.eq.1) write(*,*) 'FEM_analyze_sph_MHD'
        call FEM_analyze_sph_MHD(MHD_files1,                            &
     &      FEM_d1%geofem, FEM_d1%field, MHD_step1, visval, MHD_IO1)
!
        call end_elapsed_time(4)
!
!*  -----------  exit loop --------------
!*
        if(iflag_finish .gt. 0) exit
      end do
!
!  time evolution end
!
      call end_elapsed_time(3)
!
      if (iflag_debug.eq.1) write(*,*) 'FEM_finalize'
      call FEM_finalize(MHD_files1, MHD_step1, MHD_IO1)
!
!      if (iflag_debug.eq.1) write(*,*) 'SPH_finalize_MHD'
!      call SPH_finalize_MHD
!
      call copy_COMM_TIME_to_elaps(num_elapsed)
      call end_elapsed_time(1)
!
      if (iflag_debug.eq.1) write(*,*) 'write_resolution_data'
      call write_resolution_data(SPH_MHD1%sph)
      call output_elapsed_times
!
      call calypso_MPI_barrier
      if (iflag_debug.eq.1) write(*,*) 'exit evolution'
!
      end subroutine evolution_sph_MHD_noviz
!
! ----------------------------------------------------------------------
!
      end module analyzer_sph_MHD_noviz
