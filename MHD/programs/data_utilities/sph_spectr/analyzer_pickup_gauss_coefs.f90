!analyzer_pickup_gauss_coefs.f90
!      module analyzer_pickup_gauss_coefs
!..................................................
!
!      modified by H. Matsui on Jan., 2008
!
!      subroutine initialize_pick_gauss_coef
!      subroutine analyze_pick_gauss_coef
!
      module analyzer_pickup_gauss_coefs
!
      use m_precision
      use m_machine_parameter
      use calypso_mpi
!
      use m_schmidt_poly_on_rtm
      use t_field_data_IO
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
      subroutine initialize_pick_gauss_coef
!
      use m_t_step_parameter
      use m_ctl_data_4_sph_utils
      use m_ctl_params_sph_utils
      use m_spheric_parameter
      use m_sph_trans_comm_table
      use m_group_data_sph_specr
      use m_sph_spectr_data
      use parallel_load_data_4_sph
      use copy_rj_phys_data_4_IO
      use count_num_sph_smp
      use init_sph_trans
      use legendre_transform_select
      use cal_rms_fields_by_sph
!
!     --------------------- 
!
!
!     read controls
!
      if (iflag_debug.gt.0) write(*,*) 'read_control_data_sph_utils'
      call read_control_data_sph_utils
!
      if (iflag_debug.gt.0) write(*,*) 'set_ctl_data_4_sph_utils'
      call set_ctl_data_4_sph_utils(rj_fld1)
!
!       set spectr grids
!
      if (iflag_debug.gt.0) write(*,*) 'load_para_sph_mesh'
      call load_para_sph_mesh(sph1, comms_sph1, sph_grps1)
!
!  ------  initialize spectr data
!
      if (iflag_debug.gt.0) write(*,*) 'sel_read_alloc_step_SPH_file'
      call set_field_file_fmt_prefix                                    &
     &     (iflag_org_sph_file_fmt, org_sph_file_head, sph_fld_IO)
      call sel_read_alloc_step_SPH_file                                 &
     &   (nprocs, my_rank, i_step_init, sph_fld_IO)
!
!  -------------------------------
!
      call alloc_phys_data_type(sph1%sph_rj%nnod_rj, rj_fld1)
!
      call init_rms_4_sph_spectr                                        &
     &   (sph1%sph_params%l_truncation, sph1%sph_rj, rj_fld1)
!
      end subroutine initialize_pick_gauss_coef
!
! ----------------------------------------------------------------------
!
      subroutine analyze_pick_gauss_coef
!
      use m_t_step_parameter
      use m_spheric_parameter
      use m_ctl_params_sph_utils
      use m_sph_spectr_data
      use m_rms_4_sph_spectr
      use m_gauss_coefs_monitor_data
      use copy_rj_phys_data_4_IO
!
      use pickup_gauss_coefficients
!
      integer(kind = kint) :: i_step
!
!
      call init_gauss_coefs_4_monitor                                   &
     &   (sph1%sph_params%l_truncation, sph1%sph_rj)
      do i_step = i_step_init, i_step_number, i_step_output_ucd
!
!   Input spectr data
!
        call set_field_file_fmt_prefix                                  &
     &     (iflag_org_sph_file_fmt, org_sph_file_head, sph_fld_IO)
        call sel_read_step_SPH_field_file                               &
     &     (nprocs, my_rank, i_step, sph_fld_IO)
!
        call set_rj_phys_data_from_IO                                   &
     &     (sph1%sph_rj%nnod_rj, sph_fld_IO, rj_fld1)
!
!  pickup components
!
        call cal_gauss_coefficients                                     &
     &     (sph1%sph_params%nlayer_ICB, sph1%sph_params%nlayer_CMB,     &
     &      sph1%sph_rj%nidx_rj, sph1%sph_rj%radius_1d_rj_r,            &
     &      rj_fld1%n_point, rj_fld1%ntot_phys, rj_fld1%d_fld)
        call write_gauss_coefs_4_monitor(my_rank, i_step, time)
      end do
!
      end subroutine analyze_pick_gauss_coef
!
! ----------------------------------------------------------------------
!
      end module analyzer_pickup_gauss_coefs
