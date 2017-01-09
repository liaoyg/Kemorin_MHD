!>@file   set_control_sph_mhd.f90
!!@brief  module set_control_sph_mhd
!!
!!@author H. Matsui
!!@date    programmed by H.Matsui in Sep., 2009
!
!>@brief Set control data for spherical transform MHD dynamo simulation
!!
!!@verbatim
!!      subroutine set_control_SGS_SPH_MHD                              &
!!     &         (sph_gen, rj_fld, mesh_file, sph_file_param,           &
!!     &          MHD_org_files, sph_fst_IO, pwr, sph_filters)
!!      subroutine set_control_4_SPH_MHD(sph_gen, rj_fld,               &
!!     &          mesh_file, sph_file_param, MHD_org_files,             &
!!     &          sph_fst_IO, pwr)
!!        type(sph_grids), intent(inout) :: sph_gen
!!        type(phys_data), intent(inout) :: rj_fld
!!        type(field_IO_params), intent(inout) :: mesh_file
!!        type(field_IO_params), intent(inout) :: sph_file_param
!!        type(file_params_4_sph_mhd), intent(inout) :: MHD_org_files
!!        type(field_IO), intent(inout) :: sph_fst_IO
!!        type(sph_mean_squares), intent(inout) :: pwr
!!        type(sph_filters_type), intent(inout) :: sph_filters(1)
!!@endverbatim
!
      module set_control_sph_mhd
!
      use m_precision
!
      use m_machine_parameter
      use calypso_mpi
!
      use t_file_IO_parameter
      use t_field_data_IO
      use t_SPH_MHD_file_parameters
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine set_control_SGS_SPH_MHD                                &
     &         (sph_gen, rj_fld, mesh_file, sph_file_param,             &
     &          MHD_org_files, sph_fst_IO, pwr, sph_filters)
!
      use m_spheric_global_ranks
      use m_ucd_data
      use m_read_ctl_gen_sph_shell
      use m_ctl_data_SGS_model
      use sph_mhd_rms_IO
!
      use t_spheric_parameter
      use t_phys_data
      use t_rms_4_sph_spectr
      use t_sph_filtering_data
!
      use set_control_4_SGS
!
      type(sph_grids), intent(inout) :: sph_gen
      type(phys_data), intent(inout) :: rj_fld
      type(field_IO_params), intent(inout) :: mesh_file
      type(field_IO_params), intent(inout) :: sph_file_param
      type(file_params_4_sph_mhd), intent(inout) :: MHD_org_files
      type(field_IO), intent(inout) :: sph_fst_IO
      type(sph_mean_squares), intent(inout) :: pwr
      type(sph_filters_type), intent(inout) :: sph_filters(1)
!
!
!   set parameters for SGS model
!
      if (iflag_debug.gt.0) write(*,*) 'set_control_SGS_model'
      call set_control_SGS_model(sgs_ctl1)
      call set_control_SPH_SGS                                          &
     &   (sgs_ctl1%num_sph_filter_ctl, sgs_ctl1%sph_filter_ctl(1),      &
     &    sph_filters(1))
      if(sgs_ctl1%num_sph_filter_ctl .gt. 0) then
        call dealloc_sph_filter_ctl(sgs_ctl1)
      end if
!
      call set_control_4_SPH_MHD(sph_gen, rj_fld,                       &
     &    mesh_file, sph_file_param, MHD_org_files, sph_fst_IO, pwr)
!
      end subroutine set_control_SGS_SPH_MHD
!
! ----------------------------------------------------------------------
!
      subroutine set_control_4_SPH_MHD(sph_gen, rj_fld,                 &
     &          mesh_file, sph_file_param, MHD_org_files,               &
     &          sph_fst_IO, pwr)
!
      use m_spheric_global_ranks
      use m_ucd_data
      use m_ctl_data_4_platforms
      use m_ctl_data_4_fields
      use m_ctl_data_4_time_steps
      use m_ctl_data_4_sphere_model
      use m_ctl_data_4_pickup_sph
      use m_ctl_data_node_boundary
      use m_ctl_data_surf_boundary
      use m_ctl_data_mhd_forces
      use m_ctl_data_temp_model
      use m_ctl_data_mhd_evolution
      use m_ctl_data_mhd_evo_scheme
      use m_read_ctl_gen_sph_shell
      use sph_mhd_rms_IO
!
      use t_spheric_parameter
      use t_phys_data
      use t_rms_4_sph_spectr
!
      use gen_sph_grids_modes
      use set_control_platform_data
      use set_ctl_parallel_platform
      use set_control_4_model
      use set_control_sph_data_MHD
      use set_control_4_force
      use set_control_4_normalize
      use set_control_4_time_steps
!
      use set_control_4_velo
      use set_control_4_press
      use set_control_4_temp
      use set_control_4_magne
      use set_control_4_composition
      use set_control_4_pickup_sph
      use set_ctl_params_2nd_files
      use set_ctl_gen_shell_grids

      use check_read_bc_file
!
      type(sph_grids), intent(inout) :: sph_gen
      type(phys_data), intent(inout) :: rj_fld
      type(field_IO_params), intent(inout) :: mesh_file
      type(field_IO_params), intent(inout) :: sph_file_param
      type(file_params_4_sph_mhd), intent(inout) :: MHD_org_files
      type(field_IO), intent(inout) :: sph_fst_IO
      type(sph_mean_squares), intent(inout) :: pwr
!
      integer(kind = kint) :: ierr
!
!
!   set parameters for data files
!
      call turn_off_debug_flag_by_ctl(my_rank, plt1)
      call check_control_num_domains(plt1)
      call set_control_smp_def(my_rank, plt1)
      call set_control_mesh_def(plt1, mesh_file)
      call set_FEM_mesh_switch_4_SPH(plt1, iflag_output_mesh)
      call set_control_sph_mesh(plt1, mesh_file, sph_file_param)
      call set_control_restart_file_def(plt1, sph_fst_IO)
      call set_control_MHD_field_file
      call set_control_org_sph_files(MHD_org_files)
!
      call s_set_control_4_model(reft_ctl1, mevo_ctl1, evo_ctl1)
!
!   set spherical shell parameters
!
      iflag_make_SPH = i_sph_shell
      if(iflag_make_SPH .gt. 0) then
        if (iflag_debug.gt.0) write(*,*) 'set_control_4_shell_grids'
        call set_control_4_shell_grids                                  &
     &     (nprocs, spctl1, sdctl1, sph_gen, ierr)
      end if
!
!   set forces
!
      if (iflag_debug.gt.0) write(*,*) 's_set_control_4_force'
      call s_set_control_4_force(frc_ctl1, g_ctl1, cor_ctl1, mcv_ctl1)
!
!   set parameters for general information
!
      if (iflag_debug.gt.0) write(*,*) 's_set_control_sph_data_MHD'
      call s_set_control_sph_data_MHD                                   &
     &   (plt1, fld_ctl1%field_ctl, mevo_ctl1,                          &
     &    MHD_org_files%rj_file_param, MHD_org_files%rst_file_param,    &
     &    rj_fld)
!
!   set control parameters
!
      if (iflag_debug.gt.0) write(*,*) 's_set_control_4_normalize'
      call s_set_control_4_normalize
!
!   set boundary conditions for temperature
!
      if (iflag_debug.gt.0) write(*,*) 's_set_control_4_temp'
      call s_set_control_4_temp                                         &
     &   (nbc_ctl1%node_bc_T_ctl, sbc_ctl1%surf_bc_HF_ctl)
!
!   set boundary conditions for velocity
!
      if (iflag_debug.gt.0) write(*,*) 's_set_control_4_velo'
      call s_set_control_4_velo                                         &
     &   (nbc_ctl1%node_bc_U_ctl, sbc_ctl1%surf_bc_ST_ctl)
!
!  set boundary conditions for pressure
!
      if (iflag_debug.gt.0) write(*,*) 's_set_control_4_press'
      call s_set_control_4_press                                        &
     &   (nbc_ctl1%node_bc_P_ctl, sbc_ctl1%surf_bc_PN_ctl)
!
!   set boundary conditions for composition variation
!
      if (iflag_debug.gt.0) write(*,*) 's_set_control_4_composition'
      call s_set_control_4_composition                                  &
     &   (nbc_ctl1%node_bc_C_ctl, sbc_ctl1%surf_bc_CF_ctl)
!
!   set boundary_conditons for magnetic field
!
      if (iflag_debug.gt.0) write(*,*) 's_set_control_4_magne'
      call s_set_control_4_magne                                        &
     &   (nbc_ctl1%node_bc_B_ctl, sbc_ctl1%surf_bc_BN_ctl)
!
!   set flag to read boundary condition file
!
      if (iflag_debug.gt.0) write(*,*) 'check_read_boundary_files'
      call check_read_boundary_files
!
!   set control parameters
!
      if (iflag_debug.gt.0) write(*,*) 's_set_control_4_time_steps'
      call s_set_control_4_time_steps(tctl1)
      call s_set_control_4_crank(mevo_ctl1)
!
!   set_pickup modes
!
      call set_ctl_params_layered_spectr(smonitor_ctl1%lp_ctl, pwr)
      call set_ctl_params_sph_spectr(smonitor_ctl1, pwr)
!
      call set_ctl_params_pick_sph(smonitor_ctl1%pspec_ctl,             &
     &    pickup_sph_head, pick_list1, pick1)
!
      call set_ctl_params_pick_gauss(smonitor_ctl1%g_pwr,               &
     &    gauss_coefs_file_head, gauss_list1, gauss1)
!
      call set_ctl_params_no_heat_Nu(smonitor_ctl1%Nusselt_file_prefix, &
     &    rj_fld, Nu_type1)
!
      end subroutine set_control_4_SPH_MHD
!
! ----------------------------------------------------------------------
!
      end module set_control_sph_mhd
