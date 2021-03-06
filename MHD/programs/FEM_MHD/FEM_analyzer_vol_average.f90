!FEM_analyzer_vol_average.f90
!      module FEM_analyzer_vol_average
!
!      modified by H. Matsui on June, 2005 
!
!!      subroutine FEM_initialize_vol_average                           &
!!     &         (MHD_files, MHD_step, femmesh, ele_mesh,               &
!!     &          iphys_nod, nod_fld, FEM_model, ak_MHD, FEM_SGS,       &
!!     &          SGS_MHD_wk, MHD_IO, fem_sq, label_sim)
!!        type(mesh_data), intent(inout) :: femmesh
!!        type(element_geometry), intent(inout) :: ele_mesh
!!        type(phys_address), intent(inout) :: iphys_nod
!!        type(phys_data), intent(inout) :: nod_fld
!!        type(FEM_MHD_model_data), intent(inout) :: FEM_model
!!        type(coefs_4_MHD_type), intent(inout) :: ak_MHD
!!        type(FEM_SGS_structure), intent(inout) :: FEM_SGS
!!        type(work_FEM_SGS_MHD), intent(inout) :: SGS_MHD_wk
!!        type(MHD_file_IO_params), intent(inout) :: MHD_files
!!        type(MHD_step_param), intent(inout) :: MHD_step
!!        type(MHD_IO_data), intent(inout) :: MHD_IO
!!      subroutine FEM_analyze_vol_average                              &
!!     &         (i_step, MHD_files, femmesh, iphys_nod, FEM_model,     &
!!     &          MHD_step, SGS_MHD_wk, nod_fld, fem_sq)
!!        type(MHD_file_IO_params), intent(in) :: MHD_files
!!        type(MHD_step_param), intent(inout) :: MHD_step
!!        type(FEM_MHD_mean_square), intent(inout) :: fem_sq
!
      module FEM_analyzer_vol_average
!
      use m_precision
      use m_machine_parameter
      use m_work_time
      use t_time_data
      use t_mesh_data
      use t_phys_data
      use t_phys_address
      use t_material_property
      use t_FEM_MHD_model_data
      use t_MHD_file_parameter
      use t_MHD_step_parameter
      use t_ucd_file
      use t_FEM_MHD_mean_square
      use t_FEM_SGS_structure
      use t_work_FEM_SGS_MHD
      use t_MHD_IO_data
!
      use calypso_mpi
!
      implicit none
!
      type(time_data), save, private :: SNAP_time_IO
!
! ----------------------------------------------------------------------
!
       contains
!
! ----------------------------------------------------------------------
!
      subroutine FEM_initialize_vol_average                             &
     &         (MHD_files, MHD_step, femmesh, ele_mesh,                 &
     &          iphys_nod, nod_fld, FEM_model, ak_MHD, FEM_SGS,         &
     &          SGS_MHD_wk, MHD_IO, fem_sq, label_sim)
!
      use t_boundary_field_IO
!
      use initialize_4_snapshot
!
      use node_monitor_IO
      use open_sgs_model_coefs
!
      type(MHD_file_IO_params), intent(inout) :: MHD_files
!
      type(mesh_data), intent(inout) :: femmesh
      type(element_geometry), intent(inout) :: ele_mesh
      type(phys_address), intent(inout) :: iphys_nod
      type(phys_data), intent(inout) :: nod_fld
      type(FEM_MHD_model_data), intent(inout) :: FEM_model
      type(coefs_4_MHD_type), intent(inout) :: ak_MHD
      type(FEM_SGS_structure), intent(inout) :: FEM_SGS
      type(work_FEM_SGS_MHD), intent(inout) :: SGS_MHD_wk
      type(FEM_MHD_mean_square), intent(inout) :: fem_sq
      type(MHD_step_param), intent(inout) :: MHD_step
      type(MHD_IO_data), intent(inout) :: MHD_IO
      character(len=kchara), intent(inout)   :: label_sim
!
!   matrix assembling
!
      if (iflag_debug.eq.1)  write(*,*) 'init_analyzer_snap'
      call init_analyzer_snap(MHD_files, FEM_model%FEM_prm,             &
     &   FEM_SGS%SGS_par, FEM_model%bc_FEM_IO, MHD_step,                &
     &   femmesh%mesh, femmesh%group, ele_mesh, FEM_model%MHD_mesh,     &
     &   FEM_SGS%FEM_filters, FEM_model%MHD_prop, ak_MHD,               &
     &   FEM_model%MHD_BC, FEM_model%FEM_MHD_BCs, FEM_SGS%Csims,        &
     &   iphys_nod, nod_fld, SNAP_time_IO, MHD_step%rst_step,           &
     &   SGS_MHD_wk, fem_sq, MHD_IO%rst_IO, label_sim)
!
      end subroutine FEM_initialize_vol_average
!
! ----------------------------------------------------------------------
!
      subroutine FEM_analyze_vol_average                                &
     &         (i_step, MHD_files, femmesh, iphys_nod, FEM_model,       &
     &          MHD_step, SGS_MHD_wk, nod_fld, fem_sq)
!
      use nod_phys_send_recv
      use lead_physical_values
      use copy_nodal_fields
      use input_control
!
      use time_step_data_IO_control
      use output_parallel_ucd_file
      use FEM_MHD_evolution
!
      integer(kind = kint), intent(in) :: i_step
      type(MHD_file_IO_params), intent(in) :: MHD_files
      type(mesh_data), intent(in) :: femmesh
      type(phys_address), intent(in) :: iphys_nod
      type(FEM_MHD_model_data), intent(in) :: FEM_model
!
      type(phys_data), intent(inout) :: nod_fld
      type(work_FEM_SGS_MHD), intent(inout) :: SGS_MHD_wk
      type(FEM_MHD_mean_square), intent(inout) :: fem_sq
      type(MHD_step_param), intent(inout) :: MHD_step
!
      integer(kind = kint) :: iflag
!
!     ---- Load field data --- 
!
      if (my_rank.eq.0) write(*,*) 'step: ', i_step
!
      if (iflag_debug.eq.1)  write(*,*) 'read_udt_4_snap'
      call read_udt_4_snap(i_step, MHD_files%org_ucd_file_IO,           &
     &    nod_fld, SNAP_time_IO, MHD_step%ucd_step)
      MHD_step%time_d%time = MHD_step%init_d%time                       &
     &                      + MHD_step%time_d%dt * dble(i_step)
!
!     ---------------------
!
      call set_perturbation_to_scalar                                   &
     &   (FEM_model%MHD_prop, iphys_nod, nod_fld)
!
!     ---------------------
!
      if (iflag_debug.eq.1)  write(*,*) 'phys_send_recv_all'
      call nod_fields_send_recv(femmesh%mesh, nod_fld)
!
!     -----Output monitor date
!
      call output_time_step_control                                     &
     &   (MHD_step%flex_p%istep_max_dt, MHD_step%rms_step,              &
     &    FEM_model%FEM_prm, MHD_step%time_d, femmesh%mesh,             &
     &    FEM_model%MHD_mesh, FEM_model%MHD_prop,                       &
     &    iphys_nod, nod_fld, SGS_MHD_wk%iphys_ele,                     &
     &    SGS_MHD_wk%ele_fld, SGS_MHD_wk%fem_int%jcs,                   &
     &    fem_sq%i_rms, fem_sq%j_ave, fem_sq%i_msq,                     &
     &    SGS_MHD_wk%rhs_mat, SGS_MHD_wk%mhd_fem_wk, fem_sq%msq)
!
      end subroutine FEM_analyze_vol_average
!
! ----------------------------------------------------------------------
!
!      subroutine FEM_finalize_vol_average
!
!      end subroutine FEM_finalize_vol_average
!
!-----------------------------------------------------------------------
!
      end module FEM_analyzer_vol_average
