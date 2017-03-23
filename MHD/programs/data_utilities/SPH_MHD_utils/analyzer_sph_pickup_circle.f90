!>@file   analyzer_sph_pickup_circle.f90
!!        module analyzer_sph_pickup_circle
!!
!!@author H. Matsui
!!@date   Programmed in 2012
!!@n      modified in 2013
!
!>@brief Initialzation and evolution loop to pick up data on circle
!!
!!@verbatim
!!      subroutine initialize_sph_pick_circle
!!      subroutine evolution_sph_pick_circle
!!@endverbatim
!
      module analyzer_sph_pickup_circle
!
      use m_precision
      use calypso_mpi
!
      use m_machine_parameter
      use m_work_time
      use m_t_step_parameter
      use m_sph_trans_arrays_MHD
      use t_spheric_parameter
      use t_file_IO_parameter
      use t_step_parameter
!
      use SPH_analyzer_sph_pick_circ
!
      implicit none
!
      character(len=kchara), parameter, private                         &
     &                      :: snap_ctl_name = 'control_snapshot'
!
      type(sph_grids), private :: sph_gen
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine initialize_sph_pick_circle
!
      use t_ctl_data_sph_MHD_psf
      use m_ctl_data_sph_MHD
      use m_SGS_control_parameter
      use m_node_phys_data
      use m_spheric_parameter
      use m_sph_spectr_data
      use m_rms_4_sph_spectr
      use sph_mhd_rst_IO_control
      use set_control_sph_mhd
      use set_control_sph_data_MHD
      use init_sph_MHD_elapsed_label
      use parallel_load_data_4_sph
      use input_control_sph_MHD
!
      type(field_IO_params), save ::  mesh_file_circ
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
      call read_control_4_sph_MHD_noviz(snap_ctl_name, MHD_ctl1)
      if (iflag_debug.eq.1) write(*,*) 'set_control_SGS_SPH_MHD'
      call set_control_SGS_SPH_MHD(MHD_ctl1%plt, MHD_ctl1%org_plt,      &
     &    MHD_ctl1%model_ctl, MHD_ctl1%ctl_ctl, MHD_ctl1%smonitor_ctl,  &
     &    MHD_ctl1%nmtr_ctl, MHD_ctl1%psph_ctl,                         &
     &    sph_gen, rj_fld1, mesh_file_circ, sph_file_param1,            &
     &    MHD1_org_files, sph_fst_IO, pwr1, SGS_par1,                   &
     &    trns_WK1%dynamic_SPH%sph_filters, MHD_step1)
      call copy_delta_t(MHD_step1%init_d, MHD_step1%time_d)
!
      call set_ctl_params_pick_circle                                   &
     &   (MHD_ctl1%model_ctl%fld_ctl%field_ctl,                         &
     &    MHD_ctl1%smonitor_ctl%meq_ctl)
!
!   Load spherical harmonics data
!
      if (iflag_debug.eq.1) write(*,*) 'load_para_sph_mesh'
      call load_para_sph_mesh(sph1, comms_sph1, sph_grps1)
!
      call end_eleps_time(4)
!
!        Initialize spherical transform dynamo
!
      call start_eleps_time(2)
      if(iflag_debug .gt. 0) write(*,*) 'SPH_init_sph_pick_circle'
      call SPH_init_sph_pick_circle(iphys)
      call calypso_MPI_barrier
!
      call end_eleps_time(2)
      call reset_elapse_4_init_sph_mhd
!
      end subroutine initialize_sph_pick_circle
!
! ----------------------------------------------------------------------
!
      subroutine evolution_sph_pick_circle
!
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
        if (iflag_debug.eq.1) write(*,*) 'SPH_analyze_pick_circle'
        call SPH_analyze_pick_circle(MHD_step1%time_d%i_time_step)
!*
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
!      if (iflag_debug.eq.1) write(*,*) 'SPH_finalize_pick_circle'
!      call SPH_finalize_pick_circle
!
      call copy_COMM_TIME_to_eleps(num_elapsed)
      call end_eleps_time(1)
!
      call output_elapsed_times
!
      call calypso_MPI_barrier
      if (iflag_debug.eq.1) write(*,*) 'exit evolution'
!
      end subroutine evolution_sph_pick_circle
!
! ----------------------------------------------------------------------
!
      end module analyzer_sph_pickup_circle
