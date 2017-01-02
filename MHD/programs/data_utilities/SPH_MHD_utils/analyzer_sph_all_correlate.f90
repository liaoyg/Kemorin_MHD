!>@file   analyzer_sph_all_correlate.f90
!!@brief  module analyzer_sph_all_correlate
!!
!!@author H. Matsui
!!@date   Programmed  H. Matsui in Apr., 2010
!
!>@brief  Main loop to evaluate snapshots from spectr data
!!
!!@verbatim
!!      subroutine initialize_sph_all_correlate
!!      subroutine evolution_sph_all_correlate
!!@endverbatim
!
      module analyzer_sph_all_correlate
!
      use m_precision
      use calypso_mpi
!
      use m_machine_parameter
      use m_work_time
      use m_control_parameter
      use m_t_int_parameter
      use m_t_step_parameter
      use m_mesh_data
      use m_node_phys_data
      use m_element_id_4_node
      use m_jacobians
      use m_sph_trans_arrays_MHD
!
      use SPH_analyzer_back_trans
      use visualizer_all
!
      implicit none
!
      character(len=kchara), parameter, private                         &
     &                      :: corr_ctl_name = 'control_sph_correlate'
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine initialize_sph_all_correlate
!
      use m_ctl_data_sph_MHD
      use m_spheric_parameter
      use m_sph_spectr_data
      use m_rms_4_sph_spectr
!
      use init_sph_MHD_elapsed_label
      use FEM_analyzer_sph_MHD_w_viz
      use input_control_sph_MHD
      use SPH_analyzer_correle_all
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
      if (iflag_debug.eq.1) write(*,*) 'read_control_4_sph_snap'
      call read_control_4_sph_snap(corr_ctl_name)
!
      if (iflag_debug.eq.1) write(*,*) 'input_control_SPH_mesh'
      call input_control_SPH_mesh                                       &
     &   (sph1, comms_sph1, sph_grps1, rj_fld1, nod_fld1, pwr1,         &
     &    trns_WK1%dynamic_SPH, mesh1, group1, ele_mesh1)
      call set_ctl_4_second_spectr_data(sph_file_param2)
      call end_eleps_time(4)
!
!     --------------------- 
!
      call start_eleps_time(2)
      if(iflag_debug .gt. 0) write(*,*) 'FEM_initialize_w_viz'
      call FEM_initialize_w_viz(mesh1, group1, ele_mesh1,               &
     &    iphys, nod_fld1, next_tbl1, jac1_3d_q, jac1_3d_l)
!
!        Initialize spherical transform dynamo
      if(iflag_debug .gt. 0) write(*,*) 'SPH_init_sph_back_trans'
      call SPH_init_sph_back_trans(iphys)
!        Initialize visualization
      if(iflag_debug .gt. 0) write(*,*) 'init_visualize'
      call init_visualize(mesh1, group1, ele_mesh1, nod_fld1)
!
      call calypso_MPI_barrier
      call end_eleps_time(2)
      call reset_elapse_4_init_sph_mhd
!
      end subroutine initialize_sph_all_correlate
!
! ----------------------------------------------------------------------
!
      subroutine evolution_sph_all_correlate
!
      use m_spheric_parameter
      use m_node_phys_data
      use copy_all_fields_4_sph_trans
!
      use FEM_analyzer_sph_MHD
      use SPH_analyzer_correle_all
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
        if (iflag_debug.eq.1) write(*,*) 'SPH_analyze_correlate_all'
        call SPH_analyze_correlate_all(i_step_MHD)
!*
!*  -----------  output field data --------------
!*
        call start_eleps_time(1)
        call start_eleps_time(4)
!
        if (iflag_debug.gt.0) write(*,*) 'copy_all_field_from_trans'
        call copy_all_field_from_trans                                  &
     &     (sph1%sph_params%m_folding, sph1%sph_rtp, trns_WK1%trns_MHD, &
     &      mesh1%node, nod_fld1)
!
        if (iflag_debug.eq.1) write(*,*) 'FEM_analyze_sph_MHD'
        call FEM_analyze_sph_MHD(i_step_MHD, mesh1, nod_fld1,           &
     &      istep_psf, istep_iso, istep_pvr, istep_fline, visval)
!
        call end_eleps_time(4)
!
!*  ----------- Visualization --------------
!*
        if(visval .eq. 0) then
          if (iflag_debug.eq.1) write(*,*) 'visualize_all'
          call start_eleps_time(12)
          call visualize_all                                            &
     &       (istep_psf, istep_iso, istep_pvr, istep_fline,             &
     &        mesh1, group1, ele_mesh1, nod_fld1,                       &
     &        next_tbl1%neib_ele, jac1_3d_q)
          call end_eleps_time(12)
        end if
        call end_eleps_time(1)
!
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
      end subroutine evolution_sph_all_correlate
!
! ----------------------------------------------------------------------
!
      subroutine set_ctl_4_second_spectr_data(sph_file_param2)
!
      use t_file_IO_parameter
      use m_ctl_data_4_2nd_data
      use m_file_format_switch
!
      type(field_IO_params), intent(inout) :: sph_file_param2
!
!
      call choose_para_file_format                                      &
     &   (new_plt%spectr_field_fmt_ctl, sph_file_param2%iflag_format)
!
      sph_file_param2%iflag_IO = new_plt%spectr_field_file_prefix%iflag
      if(sph_file_param2%iflag_IO .gt. 0) then
        sph_file_param2%file_prefix                                     &
     &         = new_plt%spectr_field_file_prefix%charavalue
      end if
!
      end subroutine set_ctl_4_second_spectr_data
!
!  --------------------------------------------------------------------
!
      end module analyzer_sph_all_correlate