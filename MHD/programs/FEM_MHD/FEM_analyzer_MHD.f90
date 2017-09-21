!
!      module FEM_analyzer_MHD
!
!      modified by H. Matsui on June, 2005 
!
!!      subroutine FEM_initialize_MHD                                   &
!!     &        (MHD_files, bc_FEM_IO, flex_p, flex_data, MHD_step,     &
!!     &         femmesh, ele_mesh, iphys_nod, nod_fld, MHD_CG,         &
!!     &         FEM_SGS, SGS_MHD_wk, range, fem_ucd, fem_sq, label_sim)
!!        type(MHD_file_IO_params), intent(in) :: MHD_files
!!        type(IO_boundary), intent(in) :: bc_FEM_IO
!!        type(mesh_data), intent(inout) :: femmesh
!!        type(element_geometry), intent(inout) :: ele_mesh
!!        type(phys_address), intent(inout) :: iphys_nod
!!        type(phys_data), intent(inout) :: nod_fld
!!        type(FEM_MHD_solvers), intent(inout) :: MHD_CG
!!        type(flexible_stepping_parameter), intent(inout) :: flex_p
!!        type(flexible_stepping_data), intent(inout) :: flex_data
!!        type(MHD_step_param), intent(inout) :: MHD_step
!!        type(FEM_SGS_structure), intent(inout) :: FEM_SGS
!!        type(work_FEM_SGS_MHD), intent(inout) :: SGS_MHD_wk
!!        type(maximum_informations), intent(inout) :: range
!!        type(ucd_file_data), intent(inout) :: fem_ucd
!!        type(FEM_MHD_mean_square), intent(inout) :: fem_sq
!!      subroutine FEM_analyze_MHD                                      &
!!     &         (MHD_files, femmesh, ele_mesh, iphys_nod,              &
!!     &          MHD_step, visval, retval, MHD_CG, FEM_SGS,            &
!!     &          SGS_MHD_wk, nod_fld, fem_ucd, fem_sq)
!!        type(MHD_step_param), intent(inout) :: MHD_step
!!        type(flexible_stepping_parameter), intent(inout) :: flex_p
!!        type(flexible_stepping_data), intent(inout) :: flex_data
!!        type(mesh_data), intent(in) :: femmesh
!!        type(element_geometry), intent(in) :: ele_mesh
!!        type(phys_address), intent(in) :: iphys_nod
!!        type(IO_boundary), intent(in) :: bc_FEM_IO
!!        type(phys_data), intent(inout) :: nod_fld
!!        type(FEM_MHD_solvers), intent(inout) :: MHD_CG
!!        type(FEM_SGS_structure), intent(inout) :: FEM_SGS
!!        type(work_FEM_SGS_MHD), intent(inout) :: SGS_MHD_wk
!!        type(maximum_informations), intent(inout) :: range
!!        type(ucd_file_data), intent(inout) :: fem_ucd
!!        type(FEM_MHD_mean_square), intent(inout) :: fem_sq
!!      subroutine FEM_finalize_MHD(MHD_files, MHD_step, range, fem_ucd)
!!        type(MHD_file_IO_params), intent(in) :: MHD_files
!!        type(MHD_step_param), intent(in) :: MHD_step
!
      module FEM_analyzer_MHD
!
      use m_precision
      use m_work_time
      use m_machine_parameter
!
      use m_control_parameter
      use m_finite_element_matrix
      use m_physical_property
      use t_mesh_data
      use t_phys_data
      use t_phys_address
      use t_material_property
      use t_ucd_file
      use t_MHD_step_parameter
      use t_MHD_file_parameter
      use t_flex_delta_t_data
      use t_cal_max_indices
      use t_FEM_MHD_solvers
      use t_FEM_SGS_structure
      use t_FEM_MHD_mean_square
      use t_work_FEM_SGS_MHD
!
      use calypso_mpi
!
      implicit none
!
! ----------------------------------------------------------------------
!
       contains
!
! ----------------------------------------------------------------------
!
      subroutine FEM_initialize_MHD                                     &
     &        (MHD_files, bc_FEM_IO, flex_p, flex_data, MHD_step,       &
     &         femmesh, ele_mesh, iphys_nod, nod_fld, MHD_CG,           &
     &         FEM_SGS, SGS_MHD_wk, range, fem_ucd, fem_sq, label_sim)
!
      use m_geometry_data_MHD
      use m_bc_data_velo
      use m_flexible_time_step
      use t_boundary_field_IO
!
      use initialization_4_MHD
      use lead_physical_values
      use update_after_evolution
      use FEM_MHD_evolution
!
      use nod_phys_send_recv
      use check_deltat_by_prev_rms
      use construct_matrices
!
      use chenge_step_4_dynamic
      use output_viz_file_control
      use FEM_MHD_ucd_data
!
      type(MHD_file_IO_params), intent(in) :: MHD_files
      type(IO_boundary), intent(in) :: bc_FEM_IO
!
      type(MHD_step_param), intent(inout) :: MHD_step
      type(flexible_stepping_parameter), intent(inout) :: flex_p
      type(flexible_stepping_data), intent(inout) :: flex_data
!
      type(mesh_data), intent(inout) :: femmesh
      type(element_geometry), intent(inout) :: ele_mesh
      type(phys_address), intent(inout) :: iphys_nod
      type(phys_data), intent(inout) :: nod_fld
      type(FEM_MHD_solvers), intent(inout) :: MHD_CG
!
      type(FEM_SGS_structure), intent(inout) :: FEM_SGS
      type(work_FEM_SGS_MHD), intent(inout) :: SGS_MHD_wk
      type(maximum_informations), intent(inout) :: range
      type(ucd_file_data), intent(inout) :: fem_ucd
      type(FEM_MHD_mean_square), intent(inout) :: fem_sq
      character(len=kchara), intent(inout)   :: label_sim
!
      integer(kind = kint) :: iflag
!
!   matrix assembling
!
      call init_analyzer_fl(MHD_files, bc_FEM_IO, FEM_prm1,             &
     &    FEM_SGS%SGS_par, flex_p, flex_data, MHD_step,                 &
     &    femmesh%mesh, femmesh%group, ele_mesh, MHD_mesh1,             &
     &    FEM_SGS%FEM_filters, MHD_prop1, FEM_MHD1_BCs, FEM_SGS%Csims,  &
     &    iphys_nod, nod_fld, fem_int1, mk_MHD1, MHD_CG, SGS_MHD_wk,    &
     &    fem_sq, label_sim)
!
      call nod_fields_send_recv(femmesh%mesh, nod_fld)
!
!   obtain elemental averages
!
      call reset_update_flag(nod_fld,                                   &
     &    FEM_SGS%Csims%sgs_coefs, FEM_SGS%Csims%diff_coefs)
      if (iflag_debug.eq.1) write(*,*) 'update_FEM_fields'
      call update_FEM_fields(MHD_step%time_d,                           &
     &    FEM_prm1, FEM_SGS%SGS_par, femmesh, ele_mesh, MHD_mesh1,      &
     &    FEM_MHD1_BCs%nod_bcs, FEM_MHD1_BCs%surf_bcs, iphys_nod,       &
     &    fem_int1, FEM_SGS%FEM_filters, mk_MHD1, SGS_MHD_wk,           &
     &    nod_fld, FEM_SGS%Csims)
!
      call copy_model_coef_2_previous                                   &
     &   (FEM_SGS%SGS_par%model_p, FEM_SGS%SGS_par%commute_p,           &
     &    SGS_MHD_wk%FEM_SGS_wk)
!
!   construct matrix for Poisson and diffusion terms
!
      if (iflag_debug.eq.1) write(*,*) 'set_data_4_const_matrices'
      call set_data_4_const_matrices                                    &
     &   (femmesh, MHD_mesh1, MHD_prop1, fem_int1, MHD_CG%MGCG_WK,      &
     &    MHD1_mat_tbls, MHD_CG%MHD_mat, MHD_CG%solver_pack)
      if (iflag_debug.eq.1) write(*,*) 'set_aiccg_matrices'
      call set_aiccg_matrices(MHD_step%time_d%dt,                       &
     &    FEM_prm1, FEM_SGS%SGS_par, femmesh, ele_mesh,                 &
     &    MHD_mesh1, FEM_MHD1_BCs, MHD_prop1, fem_int1,                 &
     &    FEM_SGS%FEM_filters%FEM_elens, FEM_SGS%Csims, MHD1_mat_tbls,  &
     &    mk_MHD1, SGS_MHD_wk%rhs_mat, MHD_CG)
!
!   time evolution loop start!
!  
      call cal_FEM_model_coefficients                                   &
     &   (MHD_step%time_d, FEM_prm1, FEM_SGS%SGS_par,                   &
     &    femmesh, ele_mesh, MHD_mesh1, MHD_prop1,                      &
     &    FEM_MHD1_BCs%nod_bcs, FEM_MHD1_BCs%surf_bcs,                  &
     &    iphys_nod, fem_int1, FEM_SGS%FEM_filters,                     &
     &    mk_MHD1, SGS_MHD_wk, nod_fld, FEM_SGS%Csims)
!
      iflag = lead_field_data_flag(flex_p1%istep_max_dt, MHD_step)
      if(iflag .eq. 0) then
        if (iflag_debug.eq.1) write(*,*) 'lead_fields_by_FEM'
        call lead_fields_by_FEM(MHD_step%time_d,                        &
     &      FEM_prm1, FEM_SGS%SGS_par, femmesh, ele_mesh, MHD_mesh1,    &
     &      MHD_prop1, FEM_MHD1_BCs, iphys_nod, MHD_CG%ak_MHD,          &
     &      fem_int1, FEM_SGS%FEM_filters, mk_MHD1, SGS_MHD_wk,         &
     &      nod_fld, FEM_SGS%Csims)
      end if
!
!     ---------------------
!
      FEM_SGS%SGS_par%iflag_SGS_initial = 0
!
      call s_check_deltat_by_prev_rms                                   &
     &   (flex_p1, MHD_step%time_d, femmesh%mesh,                       &
     &    MHD_mesh1, MHD_prop1%cd_prop, iphys_nod, nod_fld,             &
     &    fem_int1%jcs, SGS_MHD_wk%rhs_mat, flex_data1)
!
!
!    Open monitor files
      call end_elapsed_time(2)
      call start_elapsed_time(4)
!
      call output_grd_file_w_org_connect                               &
     &   (MHD_step%ucd_step, femmesh%mesh, MHD_mesh1, nod_fld,         &
     &    MHD_files%ucd_file_IO, fem_ucd)
!
      call alloc_phys_range(nod_fld%ntot_phys_viz, range)
!       call s_open_boundary_monitor(my_rank, femmesh%group%sf_grp)
      call end_elapsed_time(4)
!
      end subroutine FEM_initialize_MHD
!
! ----------------------------------------------------------------------
!
      subroutine FEM_analyze_MHD                                        &
     &         (MHD_files, femmesh, ele_mesh, iphys_nod,                &
     &          MHD_step, visval, retval, MHD_CG, FEM_SGS,              &
     &          SGS_MHD_wk, nod_fld, fem_ucd, fem_sq)
!
      use m_geometry_data_MHD
      use m_bc_data_velo
      use m_flexible_time_step
!
      use construct_matrices
      use lead_physical_values
      use update_after_evolution
      use chenge_step_4_dynamic
      use copy_nodal_fields
!
      use time_step_data_IO_control
      use node_monitor_IO
      use sgs_model_coefs_IO
      use fem_mhd_rst_IO_control
      use output_viz_file_control
!
      use init_iccg_matrices
      use check_deltat_by_prev_rms
      use output_viz_file_control
      use FEM_flexible_time_step
!
      type(MHD_file_IO_params), intent(in) :: MHD_files
      type(mesh_data), intent(in) :: femmesh
      type(element_geometry), intent(in) :: ele_mesh
      type(phys_address), intent(in) :: iphys_nod
!
      integer(kind=kint ), intent(inout) :: visval
      type(MHD_step_param), intent(inout) :: MHD_step
!
      type(phys_data), intent(inout) :: nod_fld
      type(FEM_MHD_solvers), intent(inout) :: MHD_CG
      type(FEM_SGS_structure), intent(inout) :: FEM_SGS
      type(work_FEM_SGS_MHD), intent(inout) :: SGS_MHD_wk
!
      integer(kind=kint ), intent(inout) :: retval
      type(ucd_file_data), intent(inout) :: fem_ucd
      type(FEM_MHD_mean_square), intent(inout) :: fem_sq
!
      integer(kind = kint) :: iflag
      real(kind = kreal) :: total_max
!
!     ---- step to next time!! --- 
!
      if (iflag_debug.eq.1) write(*,*) 'set_new_time_and_step'
      call set_new_time_and_step                                        &
     &   (MHD_prop1%cd_prop, iphys_nod, nod_fld,                        &
     &    flex_p1, MHD_step%time_d)
!
!     ----- Time integration
!
      if (iflag_debug.eq.1) write(*,*) 'FEM_fields_evolution'
      call FEM_fields_evolution(MHD_step%time_d, FEM_prm1,              &
     &   FEM_SGS%SGS_par, femmesh, ele_mesh, MHD_mesh1, MHD_prop1,      &
     &   FEM_MHD1_BCs%nod_bcs, FEM_MHD1_BCs%surf_bcs,                   &
     &   iphys_nod, MHD_CG%ak_MHD, fem_int1, FEM_SGS%FEM_filters,       &
     &   mk_MHD1, MHD_CG%solver_pack, MHD_CG%MGCG_WK, SGS_MHD_wk,       &
     &   nod_fld, FEM_SGS%Csims, fem_sq)
!
!     ----- Evaluate model coefficients
!
      call cal_FEM_model_coefficients                                   &
     &   (MHD_step%time_d, FEM_prm1, FEM_SGS%SGS_par,                   &
     &    femmesh, ele_mesh, MHD_mesh1, MHD_prop1,                      &
     &    FEM_MHD1_BCs%nod_bcs, FEM_MHD1_BCs%surf_bcs,                  &
     &    iphys_nod, fem_int1, FEM_SGS%FEM_filters,                     &
     &    mk_MHD1, SGS_MHD_wk, nod_fld, FEM_SGS%Csims)
!
!     ---------------------
!
      if (flex_p1%iflag_flexible_step .eq. iflag_flex_step) then
        if (iflag_debug.eq.1) write(*,*) 's_check_flexible_time_step'
        call s_check_flexible_time_step(femmesh%mesh, MHD_mesh1,        &
     &      MHD_prop1%cd_prop, iphys_nod, nod_fld, fem_int1%jcs,        &
     &      SGS_MHD_wk%rhs_mat, flex_data1, flex_p1, MHD_step%time_d)
      end if
!
!     ========  Data output
!
      if(flex_p1%istep_flex_to_max .eq. 0) then
        iflag = lead_field_data_flag(flex_p1%istep_max_dt, MHD_step)
        if(iflag .eq. 0) then
          call lead_fields_by_FEM(MHD_step%time_d, FEM_prm1,            &
     &        FEM_SGS%SGS_par, femmesh, ele_mesh, MHD_mesh1,            &
     &        MHD_prop1, FEM_MHD1_BCs, iphys_nod, MHD_CG%ak_MHD,        &
     &        fem_int1, FEM_SGS%FEM_filters, mk_MHD1, SGS_MHD_wk,       &
     &        nod_fld, FEM_SGS%Csims)
        end if
!
!     -----Output monitor date
!
        call end_elapsed_time(3)
        call start_elapsed_time(4)
!
        iflag = output_IO_flag(flex_p1%istep_max_dt, MHD_step%rms_step)
        if(iflag .eq. 0) then
          if (iflag_debug.eq.1) write(*,*) 'output_time_step_control'
          call output_time_step_control                                 &
     &       (FEM_prm1, MHD_step%time_d, femmesh%mesh, MHD_mesh1,       &
     &        MHD_prop1%fl_prop, MHD_prop1%cd_prop, iphys_nod, nod_fld, &
     &        SGS_MHD_wk%iphys_ele, SGS_MHD_wk%ele_fld, fem_int1%jcs,   &
     &        fem_sq%i_rms, fem_sq%j_ave, fem_sq%i_msq,                 &
     &        SGS_MHD_wk%rhs_mat, SGS_MHD_wk%mhd_fem_wk, fem_sq%msq)
        end if
!
        iflag= output_IO_flag(flex_p1%istep_max_dt,MHD_step%point_step)
        if(iflag .eq. 0) then
          if (iflag_debug.eq.1) write(*,*) 'output_monitor_control'
          call output_monitor_control                                   &
     &       (MHD_step%time_d, femmesh%mesh%node, nod_fld)
        end if
!
        if (iflag_debug.eq.1) write(*,*) 's_output_sgs_model_coefs'
        call s_output_sgs_model_coefs(flex_p1%istep_max_dt,             &
     &      MHD_step, FEM_SGS%SGS_par, MHD_prop1%cd_prop,               &
     &      SGS_MHD_wk%FEM_SGS_wk)
!
!     ---- Output voulme field data
!
        if (iflag_debug.eq.1) write(*,*) 's_output_ucd_file_control'
        call s_output_ucd_file_control                                  &
     &     (MHD_files%ucd_file_IO, flex_p1%istep_max_dt,                &
     &      MHD_step%time_d, MHD_step%ucd_step, fem_ucd)
!
        call end_elapsed_time(4)
        call start_elapsed_time(3)
      end if
!
!
!     ---- Output restart field data
!
      iflag = set_IO_step_flag(flex_p1%istep_max_dt,MHD_step%rst_step)
      if(iflag .eq. 0) then
        if (iflag_debug.eq.1) write(*,*) 'output_MHD_restart_file_ctl'
        call output_MHD_restart_file_ctl(FEM_SGS%SGS_par, MHD_files,    &
     &      MHD_step%time_d, MHD_step%rst_step, femmesh%mesh,           &
     &      iphys_nod, SGS_MHD_wk%FEM_SGS_wk, nod_fld)
      end if
!
!     ----
!
      total_time = MPI_WTIME() - total_start
      if(iflag_debug.gt.0) write(*,*) 'total_time',                     &
     &                       total_time, MHD_step%finish_d%elapsed_time
!
      call MPI_allREDUCE (total_time, total_max, ione, CALYPSO_REAL,    &
     &    MPI_MAX, CALYPSO_COMM, ierr_MPI)
!
!
!   Finish by elapsed time
      if(MHD_step%finish_d%i_end_step .eq. -1) then
        if(total_max .gt. MHD_step%finish_d%elapsed_time) then
          MHD_step%rst_step%istep_file = MHD_step%finish_d%i_end_step
          retval = 0
          call start_elapsed_time(4)
          call output_MHD_restart_file_ctl(FEM_SGS%SGS_par, MHD_files,  &
     &        MHD_step%time_d, MHD_step%rst_step, femmesh%mesh,         &
     &        iphys_nod, SGS_MHD_wk%FEM_SGS_wk, nod_fld)
          call end_elapsed_time(4)
        end if
!
!   Finish by specific time
      else
        if(flex_p1%iflag_flexible_step .eq. iflag_flex_step) then
          if(MHD_step%time_d%time .gt. flex_p1%time_to_finish)          &
     &        retval = 0
        else
          if(flex_p1%istep_max_dt                                       &
     &        .ge. MHD_step%finish_d%i_end_step) retval = 0
        end if
      end if
!
!   Set visualization flag
      if(flex_p1%iflag_flexible_step .eq. iflag_flex_step) then
        visval = viz_file_step_4_flex(MHD_step%time_d,                  &
     &                                MHD_step%viz_step)
      else
        visval = viz_file_step_4_fix(flex_p1%istep_max_dt,              &
     &                               MHD_step%viz_step)
      end if
!
!     --------------------- 
!
      call s_chenge_step_4_dynamic                                      &
     &   (my_rank, MHD_step%time_d%i_time_step,                         &
     &    FEM_SGS%SGS_par, SGS_MHD_wk)
!
      if ( retval .ne. 0 ) then
        if (iflag_debug.eq.1) write(*,*) 'update_matrices'
        call update_matrices                                            &
     &     (MHD_step%time_d, FEM_prm1, FEM_SGS%SGS_par,                 &
     &      femmesh, ele_mesh, MHD_mesh1, FEM_MHD1_BCs, MHD_prop1,      &
     &      fem_int1, FEM_SGS%FEM_filters%FEM_elens, FEM_SGS%Csims,     &
     &      MHD1_mat_tbls, flex_p1, mk_MHD1, SGS_MHD_wk%rhs_mat,        &
     &      MHD_CG)
      end if
!
      end subroutine FEM_analyze_MHD
!
! ----------------------------------------------------------------------
!
      subroutine FEM_finalize_MHD(MHD_files, MHD_step, range, fem_ucd)
!
      type(MHD_file_IO_params), intent(in) :: MHD_files
      type(MHD_step_param), intent(in) :: MHD_step
      type(maximum_informations), intent(inout) :: range
      type(ucd_file_data), intent(inout) :: fem_ucd
!
!
      if(MHD_step%ucd_step%increment .gt. 0) then
        call finalize_output_ucd(MHD_files%ucd_file_IO, fem_ucd)
        call dealloc_phys_range(range)
      end if
!      call close_boundary_monitor(my_rank)
!
      end subroutine FEM_finalize_MHD
!
!-----------------------------------------------------------------------
!
      end module FEM_analyzer_MHD
