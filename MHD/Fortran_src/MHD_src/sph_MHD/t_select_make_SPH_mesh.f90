!>@file   t_select_make_SPH_mesh.f90
!!@brief  module t_select_make_SPH_mesh
!!
!!@author H.Matsui
!!@date     Programmed by H.Matsui in March, 2015
!
!>@brief  Load mesh and filtering data for MHD simulation
!!
!!@verbatim
!!      subroutine select_make_SPH_mesh                                 &
!!     &         (iflag_make_SPH, sph, comms_sph, sph_grps, sph_maker,  &
!!     &          geofem, ele_mesh, MHD_files)
!!        type(sph_grids), intent(inout) :: sph
!!        type(sph_comm_tables), intent(inout) :: comms_sph
!!        type(sph_group_data), intent(inout) ::  sph_grps
!!        type(mesh_data), intent(inout) ::   geofem
!!        type(element_geometry), intent(inout) :: ele_mesh
!!        type(MHD_file_IO_params), intent(inout) ::  MHD_files
!!        type(sph_grid_maker_in_sim), intent(inout) :: sph_maker
!!@endverbatim
!
!
      module t_select_make_SPH_mesh
!
      use m_precision
!
      use m_machine_parameter
      use calypso_mpi
!
      use t_control_parameter
      use t_const_spherical_grid
      use t_MHD_file_parameter
      use t_MHD_step_parameter
      use t_spheric_parameter
      use t_spheric_group
      use t_mesh_data
      use t_phys_data
      use t_rms_4_sph_spectr
      use t_file_IO_parameter
      use t_sph_trans_arrays_MHD
!
      implicit none
!
      type sph_grid_maker_in_sim
!>        Structure to construct grid
        type(construct_spherical_grid) :: gen_sph
!         Structure for temporal spherical grid
        type(sph_grids) :: sph_tmp
      end type sph_grid_maker_in_sim
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine select_make_SPH_mesh                                   &
     &         (iflag_make_SPH, sph, comms_sph, sph_grps, sph_maker,    &
     &          geofem, ele_mesh, MHD_files)
!
      use m_error_IDs
      use parallel_load_data_4_sph
      use parallel_gen_sph_grids
!
      integer(kind = kint), intent(in) :: iflag_make_SPH
      type(sph_grids), intent(inout) :: sph
      type(sph_comm_tables), intent(inout) :: comms_sph
      type(sph_group_data), intent(inout) ::  sph_grps
!
      type(mesh_data), intent(inout) ::   geofem
      type(element_geometry), intent(inout) :: ele_mesh
      type(MHD_file_IO_params), intent(inout) ::  MHD_files
      type(sph_grid_maker_in_sim), intent(inout) :: sph_maker
!
      integer(kind = kint) :: iflag_lc
!
!
      call calypso_mpi_barrier
!
      if(my_rank .eq. izero) then
        iflag_lc = 0
        if     (check_exsist_rtp_file(my_rank) .ne. 0                   &
     &     .or. check_exsist_rtm_file(my_rank) .ne. 0                   &
     &     .or. check_exsist_rlm_file(my_rank) .ne. 0                   &
     &     .or. check_exsist_rj_file(my_rank) .ne.  0) iflag_lc = 1
      end if
      call MPI_BCAST(iflag_lc, ione, CALYPSO_INTEGER, izero,            &
     &    CALYPSO_COMM, ierr_MPI)
!
      if(iflag_lc .eq. 0) then
        if(my_rank.eq.0) write(*,*) 'spherical harmonics table exists'
      else if(iflag_make_SPH .eq. 0) then
        call calypso_mpi_abort(ierr_file,                               &
     &     'Set parameters for spherical shell')
      else
        if (my_rank.eq.0) write(*,*) 'Make spherical harmonics table'
        call para_gen_sph_grids(sph_maker%sph_tmp, sph_maker%gen_sph)
        call deallocate_gen_mesh_params(sph_maker%gen_sph)
      end if
      call calypso_mpi_barrier
!
      if (iflag_debug.eq.1) write(*,*) 'load_para_SPH_and_FEM_mesh'
      call load_para_SPH_and_FEM_mesh(MHD_files%iflag_access_FEM,       &
     &    sph, comms_sph, sph_grps, geofem%mesh, geofem%group,          &
     &    ele_mesh, MHD_files%mesh_file_IO, sph_maker%gen_sph)
!
      end subroutine select_make_SPH_mesh
!
! ----------------------------------------------------------------------
!
      end module t_select_make_SPH_mesh
