!>@file   analyzer_sph_forward_trans.f90
!!@brief  module analyzer_sph_forward_trans
!!
!!@author H. Matsui
!!@date   Programmed  H. Matsui in Apr., 2010
!
!>@brief  Main loop to evaluate snapshots from spectr data
!!        without visualization routines
!!
!!@verbatim
!!      subroutine initialize_sph_forward_trans
!!      subroutine evolution_sph_forward_trans
!!@endverbatim
!
      module analyzer_sph_forward_trans
!
      use m_precision
      use calypso_mpi
!
      use m_machine_parameter
      use calypso_mpi
      use m_work_time
      use m_control_parameter
      use m_t_int_parameter
      use m_t_step_parameter
      use m_mesh_data
      use m_sph_trans_arrays_MHD
!
      use FEM_analyzer_sph_MHD
      use SPH_analyzer_snap
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
      subroutine initialize_sph_forward_trans
!
      use m_spheric_parameter
      use m_mesh_data
      use m_node_phys_data
      use m_sph_spectr_data
      use m_rms_4_sph_spectr
      use m_cal_max_indices
      use m_ctl_data_sph_MHD_noviz
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
      if (iflag_debug.eq.1) write(*,*) 'read_control_4_sph_snap_noviz'
      call read_control_4_sph_snap_noviz(snap_ctl_name)
!
      if (iflag_debug.eq.1) write(*,*) 'input_control_SPH_mesh'
      call input_control_SPH_mesh                                       &
     &   (sph1, comms_sph1, sph_grps1, rj_fld1, nod_fld1, pwr1,         &
     &    trns_WK1%dynamic_SPH, mesh1, group1, ele_mesh1)
      call end_eleps_time(4)
!
!     --------------------- 
!
      call start_eleps_time(2)
      if(iflag_debug .gt. 0) write(*,*) 'FEM_initialize_sph_MHD'
      call FEM_initialize_sph_MHD(mesh1, group1, ele_mesh1,             &
     &    iphys, nod_fld1, range)
!
!        Initialize spherical transform dynamo
      if(iflag_debug .gt. 0) write(*,*) 'SPH_init_sph_snap'
      call SPH_init_sph_snap(iphys)
!
      call calypso_MPI_barrier
!
      call end_eleps_time(2)
      call reset_elapse_4_init_sph_mhd
!
      end subroutine initialize_sph_forward_trans
!
! ----------------------------------------------------------------------
!
      subroutine evolution_sph_forward_trans
!
      use m_spheric_parameter
      use m_node_phys_data
!
      integer(kind = kint) :: visval
      integer(kind = kint) :: istep_psf, istep_iso
      integer(kind = kint) :: istep_pvr, istep_fline
!
!     ---------------------
!
      call start_eleps_time(3)
!
!*  -----------  set initial step data --------------
!*
      i_step_MHD = i_step_init - 1
!*
!*  -------  time evelution loop start -----------
!*
      do
        i_step_MHD = i_step_MHD + 1
        istep_max_dt = i_step_MHD
!
        if( mod(i_step_MHD,i_step_output_rst) .ne. 0) cycle
!
!*  ----------  time evolution by spectral methood -----------------
!*
        if (iflag_debug.eq.1) write(*,*) 'SPH_analyze_snap'
        call SPH_analyze_snap(i_step_MHD)
!*
!*  -----------  exit loop --------------
!*
        if(i_step_MHD .ge. i_step_number) exit
      end do
!
!  time evolution end
!
      call end_eleps_time(3)
!
      if (iflag_debug.eq.1) write(*,*) 'FEM_finalize'
      call FEM_finalize
!
!      if (iflag_debug.eq.1) write(*,*) 'SPH_finalize_snap'
!      call SPH_finalize_snap
!
      call copy_COMM_TIME_to_eleps(num_elapsed)
      call end_eleps_time(1)
!
      call output_elapsed_times
!
      call calypso_MPI_barrier
      if (iflag_debug.eq.1) write(*,*) 'exit evolution'
!
      end subroutine evolution_sph_forward_trans
!
! ----------------------------------------------------------------------
!
      end module analyzer_sph_forward_trans