!>@file   analyzer_sph_zm_snap_w_psf.f90
!!@brief  module analyzer_sph_zm_snap_w_psf
!!
!!@author H. Matsui
!!@date   Programmed  H. Matsui in Apr., 2010
!
!>@brief  Main loop to evaluate zonal mean field
!!        including visualization module
!!
!!@verbatim
!!      subroutine initialize_sph_zm_snap_w_psf
!!      subroutine evolution_sph_zm_snap_w_psf
!!@endverbatim
!
      module analyzer_sph_zm_snap_w_psf
!
      use m_precision
      use calypso_mpi
!
      use m_machine_parameter
      use m_MHD_step_parameter
      use m_SGS_control_parameter
      use m_work_time
      use m_spheric_parameter
      use m_mesh_data
      use m_node_phys_data
      use m_sph_trans_arrays_MHD
      use t_step_parameter
!
      use FEM_analyzer_sph_MHD
      use sections_for_1st
!
      implicit none
!
      character(len=kchara), parameter, private                         &
     &                      :: snap_ctl_name = 'control_snapshot'
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine initialize_sph_zm_snap_w_psf
!
      use t_ctl_data_sph_MHD_psf
      use m_ctl_data_sph_MHD
      use m_sph_spectr_data
      use m_rms_4_sph_spectr
      use m_cal_max_indices
      use m_bc_data_list
      use init_sph_MHD_elapsed_label
      use SPH_analyzer_snap
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
      if (iflag_debug.eq.1) write(*,*) 'read_control_4_sph_MHD_w_psf'
      call read_control_4_sph_MHD_w_psf(snap_ctl_name, MHD_ctl1)
!
      if (iflag_debug.eq.1) write(*,*) 'input_control_SPH_dynamo'
      call input_control_SPH_dynamo                                     &
     &   (MHD_files1, MHD_ctl1, sph1, comms_sph1, sph_grps1, rj_fld1,   &
     &    nod_fld1, pwr1, SGS_par1, dynamic_SPH1, MHD_step1,            &
     &    MHD_prop1, MHD_BC1, trns_WK1, mesh1, group1, ele_mesh1)
      call copy_delta_t(MHD_step1%init_d, MHD_step1%time_d)
      call end_eleps_time(4)
!
!     --------------------- 
!
      call start_eleps_time(2)
      if(iflag_debug .gt. 0) write(*,*) 'FEM_initialize_sph_MHD'
      call FEM_initialize_sph_MHD(MHD_files1%ucd_file_IO, MHD_step1,    &
     &    mesh1, group1, ele_mesh1, iphys, nod_fld1, range)
!
!        Initialize spherical transform dynamo
      if(iflag_debug .gt. 0) write(*,*) 'SPH_init_sph_snap'
      call SPH_init_sph_snap(MHD_files1, iphys)
!
!        Initialize visualization
!
      if(iflag_debug .gt. 0) write(*,*) 'init_visualize_surface'
      call init_visualize_surface(mesh1, group1, ele_mesh1, nod_fld1)
!
      call calypso_MPI_barrier
      call end_eleps_time(2)
      call reset_elapse_4_init_sph_mhd
!
      end subroutine initialize_sph_zm_snap_w_psf
!
! ----------------------------------------------------------------------
!
      subroutine evolution_sph_zm_snap_w_psf
!
      use SPH_analyzer_zm_snap
      use output_viz_file_control
!
      integer(kind = kint) :: visval
      integer(kind = kint) :: iflag
!
!*  -----------  set initial step data --------------
!*
      call start_eleps_time(3)
      call s_initialize_time_step(MHD_step1%init_d, MHD_step1%time_d)
!*
!*  -------  time evelution loop start -----------
!*
      do
        call add_one_step(MHD_step1%time_d)
!
        iflag = output_IO_flag(MHD_step1%time_d%i_time_step,            &
     &                         MHD_step1%rst_step)
        if(iflag .ne. 0) cycle
!
!*  ----------  time evolution by spectral methood -----------------
!*
        if (iflag_debug.eq.1) write(*,*) 'SPH_analyze_zm_snap'
        call SPH_analyze_zm_snap                                        &
     &     (MHD_step1%time_d%i_time_step, MHD_files1, MHD_step1)
!*
!*  -----------  output field data --------------
!*
        call start_eleps_time(1)
        call start_eleps_time(4)
        iflag = lead_field_data_flag(MHD_step1%time_d%i_time_step,      &
     &                               MHD_step1)
        if(iflag .eq. 0) then
          if(iflag_debug.eq.1) write(*,*) 'SPH_to_FEM_bridge_zm_snap'
          call SPH_to_FEM_bridge_zm_snap                                &
     &       (sph1%sph_params, sph1%sph_rtp, trns_WK1,                  &
     &        mesh1, iphys, nod_fld1)
        end if
!
        if (iflag_debug.eq.1) write(*,*) 'FEM_analyze_sph_MHD'
        call FEM_analyze_sph_MHD(MHD_files1%ucd_file_IO,                &
     &      mesh1,  nod_fld1, MHD_step1, visval)
!
        call end_eleps_time(4)
!
!*  ----------- Visualization --------------
!*
        if(visval .eq. 0) then
          if (iflag_debug.eq.1) write(*,*) 'visualize_surface'
          call start_eleps_time(8)
          call visualize_surface(MHD_step1%viz_step, MHD_step1%time_d,  &
     &        mesh1, ele_mesh1, nod_fld1)
          call end_eleps_time(8)
        end if
        call end_eleps_time(1)
!
!*  -----------  exit loop --------------
!*
        if(MHD_step1%time_d%i_time_step                                 &
     &        .ge. MHD_step1%finish_d%i_end_step) exit
      end do
!
!  time evolution end
!
      call end_eleps_time(3)
!
      if (iflag_debug.eq.1) write(*,*) 'FEM_finalize'
      call FEM_finalize(MHD_files1%ucd_file_IO, MHD_step1)
!
!      if (iflag_debug.eq.1) write(*,*) 'SPH_finalize_snap'
!      call SPH_finalize_snap
!
      call copy_COMM_TIME_to_elaps(num_elapsed)
      call end_eleps_time(1)
!
      call output_elapsed_times
!
      call calypso_MPI_barrier
      if (iflag_debug.eq.1) write(*,*) 'exit evolution'
!
      end subroutine evolution_sph_zm_snap_w_psf
!
! ----------------------------------------------------------------------
!
      end module analyzer_sph_zm_snap_w_psf