!analyzer_test_sph.f90
!      module analyzer_test_sph
!..................................................
!
!      modified by H. Matsui on Aug., 2007
!
!      subroutine init_test_sph
!      subroutine analyze_test_sph
!
      module analyzer_test_sph
!
      use m_precision
      use m_constants
      use m_machine_parameter
      use calypso_mpi
!
      implicit none
!
      integer(kind = kint), parameter :: id_check = 44
      character(len=kchara), parameter :: check_header = 'comm_errors'
      private :: check_header
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine init_test_sph
!
      use m_spheric_parameter
      use m_sph_trans_comm_table
      use m_group_data_sph_specr
      use m_read_ctl_gen_sph_shell
      use set_control_platform_data
      use parallel_load_data_4_sph
      use cmp_trans_sph_tests
!
!
      if (my_rank.eq.0) then
        write(*,*) 'Construct commutation filter'
        write(*,*) 'Input file: mesh data'
      end if
!
      call set_tesh_sph_elapsed_label
!
!     --------------------- 
!
      call turn_off_debug_flag_by_ctl(my_rank)
      call read_control_4_gen_shell_grids
      call set_control_sph_mesh
!
      if (iflag_debug.gt.0) write(*,*) 'load_para_sph_mesh'
      call load_para_sph_mesh                                           &
     &   (sph_param1, sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1,            &
     &    comm_rtp1, comm_rtm1, comm_rlm1, comm_rj1, bc_rtp_grp1,       &
     &    radial_rtp_grp1, theta_rtp_grp1, zonal_rtp_grp,               &
     &    radial_rj_grp1, sphere_rj_grp1)
!
       end subroutine init_test_sph
!
! ----------------------------------------------------------------------
!
      subroutine analyze_test_sph
!
      use m_solver_SR
      use m_spheric_parameter
      use m_sph_trans_comm_table
      use cmp_trans_sph_indices
      use cmp_trans_sph_tests
      use set_parallel_file_name
      use select_calypso_SR
      use select_copy_from_recv
!
      character(len=kchara) :: fname_tmp, file_name
      integer(kind = kint) :: itype
      integer(kind = kint), parameter :: NB = 8
!
!
      call allocate_idx_sph_recieve                                     &
     &     (sph_rtp1%nnod_rtp, sph_rtm1%nnod_rtm,                       &
     &      sph_rlm1%nnod_rlm, sph_rj1%nnod_rj)
      call allocate_real_sph_test                                       &
     &   (NB, sph_rtp1%nnod_rtp, sph_rtm1%nnod_rtm,                     &
     &        sph_rlm1%nnod_rlm, sph_rj1%nnod_rj)
!
      call add_int_suffix(my_rank, check_header, fname_tmp)
      call add_dat_extension(fname_tmp, file_name)
      write(*,*) 'error check result: ', trim(file_name)
      open(id_check, file=file_name)
!
      do itype = iflag_import_item, iflag_import_rev
!
       if(itype .eq. iflag_import_item) then
         write(id_check,*)  'USING IMPORT_ITEM'
       else if(itype .eq. iflag_import_rev) then
         write(id_check,*)  'USING IMPORT_REVERSE'
       end if
!
        call sph_indices_transfer                                       &
     &     (itype, sph_rtp1%nnod_rtp, sph_rtm1%nnod_rtm,                &
     &      sph_rlm1%nnod_rlm, sph_rj1%nnod_rj,                         &
     &      sph_rtp1%idx_global_rtp, sph_rtm1%idx_global_rtm,           &
     &      sph_rlm1%idx_global_rlm, sph_rj1%idx_global_rj)
        call check_missing_sph_indices(id_check,                        &
     &      sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1)
        call compare_transfer_sph_indices(id_check,                     &
     &      sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1)
!
        call sph_transfer_test_N(NB, itype,                             &
     &      sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1,                      &
     &      comm_rtp1, comm_rtm1, comm_rlm1, comm_rj1)
        call compare_transfer_sph_reals(NB, id_check,                   &
     &      sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1)
!
        call sph_transfer_test_6(itype,                                 &
     &      sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1,                      &
     &      comm_rtp1, comm_rtm1, comm_rlm1, comm_rj1)
        call compare_transfer_sph_reals(isix, id_check,                 &
     &      sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1)
!
        call sph_transfer_test_3(itype,                                 &
     &      sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1,                      &
     &      comm_rtp1, comm_rtm1, comm_rlm1, comm_rj1)
        call compare_transfer_sph_reals(ithree, id_check,               &
     &      sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1)
!
        call sph_transfer_test_2(itype,                                 &
     &      sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1,                      &
     &      comm_rtp1, comm_rtm1, comm_rlm1, comm_rj1)
        call compare_transfer_sph_reals(itwo, id_check,                 &
     &      sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1)
!
        call sph_transfer_test_1(itype,                                 &
     &      sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1,                      &
     &      comm_rtp1, comm_rtm1, comm_rlm1, comm_rj1)
        call compare_transfer_sph_reals(ione, id_check,                 &
     &      sph_rtp1, sph_rtm1, sph_rlm1, sph_rj1)
      end do
!
      close(id_check)
!
      call deallocate_real_sph_test
      call deallocate_idx_sph_recieve
!
      if (iflag_debug.eq.1) write(*,*) 'exit analyze_test_sph'
!
      end subroutine analyze_test_sph
!
! ----------------------------------------------------------------------
!
      end module analyzer_test_sph
