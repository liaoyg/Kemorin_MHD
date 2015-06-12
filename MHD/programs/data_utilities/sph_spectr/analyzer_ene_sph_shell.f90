!analyzer_ene_sph_shell.f90
!      module analyzer_ene_sph_shell
!..................................................
!
!      modified by H. Matsui on Jan., 2008
!
!      subroutine initialize_ene_sph_shell
!      subroutine analyze_ene_sph_shell
!
      module analyzer_ene_sph_shell
!
      use m_precision
      use m_constants
      use m_machine_parameter
      use calypso_mpi
!
      use m_schmidt_poly_on_rtm
      use t_field_data_IO
      use cal_rms_fields_by_sph
      use field_IO_select
!
      implicit none
!
      type(field_IO), save, private :: sph_fld_IO
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine initialize_ene_sph_shell
!
      use m_t_step_parameter
      use m_ctl_data_4_sph_utils
      use m_ctl_params_sph_utils
      use m_control_params_sph_data
      use m_node_phys_address
      use m_sph_spectr_data
      use m_sph_phys_address
      use load_mesh_data
      use const_mesh_info
      use parallel_load_data_4_sph
      use copy_rj_phys_data_4_IO
      use count_num_sph_smp
      use init_sph_trans
      use legendre_transform_select
      use cal_rms_fields_by_sph
!
!     --------------------- 
!
!     read controls
!
      if (iflag_debug.gt.0) write(*,*) 'read_control_data_sph_utils'
      call read_control_data_sph_utils
!
      if (iflag_debug.gt.0) write(*,*) 'set_ctl_data_4_sph_utils'
      call set_ctl_data_4_sph_utils
!
!       set mesh informations
!
      if (iflag_debug.gt.0) write(*,*) 'input_mesh'
      call input_mesh(my_rank)
!
      if (iflag_debug.gt.0) write(*,*) 'const_mesh_informations'
      call const_mesh_informations(my_rank)
!
!       set spectr grids
!
      if (iflag_debug.gt.0) write(*,*) 'load_para_sph_mesh'
      call load_para_sph_mesh
!
!  ------  initialize spectr data
!
      if (iflag_debug.gt.0) write(*,*) 'sel_read_alloc_step_SPH_file'
      call set_field_file_fmt_prefix                                    &
     &   (iflag_sph_spectr_fmt, spectr_file_head, sph_fld_IO)
      call sel_read_alloc_step_SPH_file                                 &
     &   (nprocs, my_rank, i_step_init, sph_fld_IO)
!
!  -------------------------------
!
      call allocate_phys_rj_data
      call set_sph_sprctr_data_address
!
      call init_rms_4_sph_spectr
!
!  -------------------------------
!
      if (iflag_debug.gt.0) write(*,*) 'initialize_sph_trans'
      if(id_legendre_transfer.eq.iflag_leg_undefined)                   &
     &            id_legendre_transfer = iflag_leg_orginal_loop
      call copy_sph_trans_nums_from_rtp
      call initialize_sph_trans
!
!      call check_schmidt_poly_rtm(my_rank+40)
!
      end subroutine initialize_ene_sph_shell
!
! ----------------------------------------------------------------------
!
      subroutine analyze_ene_sph_shell
!
      use m_t_step_parameter
      use m_spheric_parameter
      use m_spheric_param_smp
      use m_ctl_params_sph_utils
      use m_control_params_sph_data
      use m_sph_spectr_data
      use m_rms_4_sph_spectr
      use copy_rj_phys_data_4_IO
      use output_sph_m_square_file
      use volume_average_4_sph
!
!
      integer(kind = kint) :: i_step
!
!
      call set_field_file_fmt_prefix                                    &
     &   (iflag_sph_spectr_fmt, spectr_file_head, sph_fld_IO)
!
      do i_step = i_step_init, i_step_number, i_step_output_ucd
!
!   Input spectr data
!
      call sel_read_step_SPH_field_file                                 &
     &     (nprocs, my_rank, i_step, sph_fld_IO)
!
        call set_rj_phys_data_from_IO(sph_fld_IO)
!
!  evaluate energies
!
        if (iflag_debug.gt.0) write(*,*) 'cal_rms_sph_spec_rms_whole'
        call cal_rms_sph_spec_rms_whole
!
        call write_sph_vol_ave_file(i_step, time)
        call write_sph_vol_ms_file(my_rank, i_step, time)
        call write_sph_vol_ms_spectr_file(my_rank, i_step, time)
        call write_sph_layer_ms_file(my_rank, i_step, time)
      end do
!
      end subroutine analyze_ene_sph_shell
!
! ----------------------------------------------------------------------
!
      end module analyzer_ene_sph_shell
