!>@file   check_dependency_for_MHD.f90
!!@brief  module check_dependency_for_MHD
!!
!!@author H. Matsui
!!@date Programmed in Sep., 2007
!
!>@brief  Check dependecy of field list fro MHD dynamo
!!
!!@verbatim
!!      subroutine set_FEM_MHD_field_data                               &
!!     &         (SGS_param, cmt_param, node, iphys, nod_fld)
!!        type(commutation_control_params), intent(in) :: cmt_param
!!        type(node_data), intent(in) :: node
!!        type(phys_address), intent(inout) :: iphys
!!        type(phys_data), intent(inout) :: nod_fld
!!      subroutine set_sph_MHD_sprctr_data                              &
!!     &         (SGS_param, sph_rj, ipol, idpdr, itor, rj_fld)
!!@endverbatim
!
      module check_dependency_for_MHD
!
      use m_precision
      use m_error_IDs
      use m_machine_parameter
!
      use calypso_mpi
      use m_phys_labels
!
      use t_time_stepping_parameter
      use t_SGS_control_parameter
      use t_phys_data
      use t_phys_address
      use t_physical_property
!
      implicit none
!
      private :: check_missing_field_w_msg
      private :: check_field_dependencies
      private :: check_dependence_FEM_evo, check_dependence_SPH_evo
      private :: check_dependence_4_FEM_SGS, check_dependence_4_SPH_SGS
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine set_FEM_MHD_field_data                                 &
     &         (SGS_param, cmt_param, node, iphys, nod_fld)
!
      use m_control_parameter
      use m_physical_property
      use t_geometry_data
      use t_FEM_phys_data
      use check_MHD_dependency_by_id
!
      type(SGS_model_control_params), intent(in) :: SGS_param
      type(commutation_control_params), intent(in) :: cmt_param
      type(node_data), intent(in) :: node
      type(phys_address), intent(inout) :: iphys
      type(phys_data), intent(inout) :: nod_fld
!
!
      if (iflag_debug.ge.1)  write(*,*) 'init_field_address'
      call init_field_address(node%numnod, nod_fld, iphys)
!
      call check_field_dependencies                                     &
     &   (evo_magne, evo_vect_p, evo_temp, evo_comp,                    &
     &    fl_prop1, iphys, nod_fld)
      call check_dependencies_by_id(evo_magne, iphys, nod_fld)
      call check_dependence_FEM_MHD_by_id(iphys, nod_fld)
      call check_dependence_FEM_evo(fl_prop1, iphys, nod_fld)
      call check_dependence_4_FEM_SGS                                   &
     &   (evo_magne, evo_vect_p, evo_temp, evo_comp,                    &
     &    SGS_param, cmt_param, fl_prop1, iphys, nod_fld)
!
      end subroutine set_FEM_MHD_field_data
!
! -----------------------------------------------------------------------
!
      subroutine set_sph_MHD_sprctr_data                                &
     &         (SGS_param, sph_rj, ipol, idpdr, itor, rj_fld)
!
      use m_control_parameter
      use m_physical_property
      use t_spheric_rj_data
!
      use set_sph_phys_address
      use check_MHD_dependency_by_id
!
      type(SGS_model_control_params), intent(in) :: SGS_param
      type(sph_rj_grid), intent(in) :: sph_rj
      type(phys_address), intent(inout) :: ipol, idpdr, itor
      type(phys_data), intent(inout) :: rj_fld
!
!
      call set_sph_sprctr_data_address                                  &
     &   (sph_rj, ipol, idpdr, itor, rj_fld)
!
      call check_field_dependencies                                     &
     &   (evo_magne, evo_vect_p, evo_temp, evo_comp,                    &
     &    fl_prop1, ipol, rj_fld)
      call check_dependencies_by_id(evo_magne, ipol, rj_fld)
      call check_dependence_SPH_MHD_by_id(ipol, rj_fld)
      call check_dependence_SPH_evo(fl_prop1, ipol, rj_fld)
      call check_dependence_4_SPH_SGS                                   &
     &   (evo_magne, evo_temp, evo_comp,                                &
     &    SGS_param, fl_prop1, ipol, rj_fld)
!
      end subroutine set_sph_MHD_sprctr_data
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine check_field_dependencies                               &
     &         (evo_B, evo_A, evo_T, evo_C, fl_prop, iphys, fld)
!
      type(time_evolution_params), intent(in) :: evo_B, evo_A
      type(time_evolution_params), intent(in) :: evo_T, evo_C
      type(fluid_property), intent(in) :: fl_prop
      type(phys_address), intent(in) :: iphys
      type(phys_data), intent(in) :: fld
!
      character(len=kchara) :: msg
!
!
!   check dependencies for time evolution
!
      if (cd_prop%iflag_Aevo_scheme .gt. id_no_evolution                &
     &     .and. evo_B%iflag_scheme .gt. id_no_evolution) then
         call calypso_MPI_abort(ierr_fld,                               &
     &        'You should choose vector potential OR magnetic field for &
     & time evolution')
      end if
!
      if (fl_prop%iflag_scheme .gt. id_no_evolution) then
        msg = 'time integration for velocity needs'
        call check_missing_field_w_msg(fld, msg, iphys%i_velo)
        call check_missing_field_w_msg(fld, msg, iphys%i_press)
      end if
!
      if (fl_prop%iflag_scheme .gt. id_no_evolution) then
        msg = 'time integration for velocity needs'
        call check_missing_field_w_msg(fld, msg, iphys%i_vort)
      end if
!
      if (evo_T%iflag_scheme .gt. id_no_evolution) then
        msg = 'Time integration for temperature needs'
        call check_missing_field_w_msg(fld, msg, iphys%i_velo)
        call check_missing_field_w_msg(fld, msg, iphys%i_temp)
      end if
!
      if (evo_C%iflag_scheme .ne. id_no_evolution) then
        msg =  'Time integration for composition needs'
        call check_missing_field_w_msg(fld, msg, iphys%i_velo)
        call check_missing_field_w_msg(fld, msg, iphys%i_light)
      end if
!
      if (evo_B%iflag_scheme .ne. id_no_evolution) then
        msg = 'Time integration for magnetic field needs'
        call check_missing_field_w_msg(fld, msg, iphys%i_magne)
        call check_missing_field_w_msg(fld, msg, iphys%i_velo)
      end if
!
      if (cd_prop%iflag_Aevo_scheme .gt. id_no_evolution) then
        msg = 'Time integration for vector potential needs'
        call check_missing_field_w_msg(fld, msg, iphys%i_vecp)
        call check_missing_field_w_msg(fld, msg, iphys%i_velo)
      end if
!
!
      if ( fl_prop%iflag_scheme .gt. id_no_evolution) then
        if (fl_prop%iflag_4_gravity .gt. id_turn_OFF) then
          msg = 'Buoyancy needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_temp)
        end if
!
        if (fl_prop%iflag_4_composit_buo .gt. id_turn_OFF) then
          msg = 'Compositional buoyancy needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_light)
        end if
!
        if (fl_prop%iflag_4_filter_gravity .gt. id_turn_OFF) then
          msg = 'Filtered buoyancy needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_filter_temp)
        end if
!
        if (fl_prop%iflag_4_lorentz .gt. id_turn_OFF) then
          msg = 'Lorentz force needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_magne)
        end if
      end if
!
      end subroutine check_field_dependencies
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine check_dependence_FEM_evo(fl_prop, iphys, fld)
!
      type(fluid_property), intent(in) :: fl_prop
      type(phys_address), intent(in) :: iphys
      type(phys_data), intent(in) :: fld
!
      character(len=kchara) :: msg
!
!
      if (fl_prop%iflag_scheme .gt. id_no_evolution) then
        msg = 'time integration for velocity needs'
        call check_missing_field_w_msg(fld, msg, iphys%i_velo)
        call check_missing_field_w_msg(fld, msg, iphys%i_press)
      end if
!
      end subroutine check_dependence_FEM_evo
!
! -----------------------------------------------------------------------
!
      subroutine check_dependence_SPH_evo(fl_prop, iphys, fld)
!
      type(fluid_property), intent(in) :: fl_prop
      type(phys_address), intent(in) :: iphys
      type(phys_data), intent(in) :: fld
!
      character(len=kchara) :: msg
!
!
      if (fl_prop%iflag_scheme .gt. id_no_evolution) then
        msg = 'time integration for velocity needs'
        call check_missing_field_w_msg(fld, msg, iphys%i_vort)
      end if
!
      end subroutine check_dependence_SPH_evo
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine check_dependence_4_FEM_SGS                             &
     &         (evo_B, evo_A, evo_T, evo_C,                             &
     &          SGS_param, cmt_param, fl_prop, iphys, fld)
!
      type(time_evolution_params), intent(in) :: evo_B, evo_A
      type(time_evolution_params), intent(in) :: evo_T, evo_C
      type(SGS_model_control_params), intent(in) :: SGS_param
      type(commutation_control_params), intent(in) :: cmt_param
      type(fluid_property), intent(in) :: fl_prop
      type(phys_address), intent(in) :: iphys
      type(phys_data), intent(in) :: fld
!
!
      character(len=kchara) :: msg
!
!
      if (fl_prop%iflag_scheme .gt. id_no_evolution) then
        if ( SGS_param%iflag_SGS_m_flux .ne. id_SGS_none) then
          msg = 'solving SGS momentum flux needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_SGS_m_flux)
        end if
!
        if (SGS_param%iflag_SGS_lorentz .ne. id_SGS_none) then
          msg = 'solving SGS lorentz term needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_SGS_maxwell)
        end if
      end if
!
!
      if ( evo_T%iflag_scheme .gt. id_no_evolution) then
        if ( SGS_param%iflag_SGS_h_flux .ne. id_SGS_none) then
          msg = 'solving SGS heat flux needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_SGS_h_flux)
        end if
      end if
!
!
      if ( evo_B%iflag_scheme .gt. id_no_evolution) then
        if ( SGS_param%iflag_SGS_uxb .ne. id_SGS_none) then
          msg = 'solving SGS magnetic induction needs'
          call check_missing_field_w_msg                                &
     &       (fld, msg, iphys%i_SGS_induct_t)
          call check_missing_field_w_msg                                &
     &       (fld, msg, iphys%i_SGS_vp_induct)
        end if
      end if
!
!
      if ( cd_prop%iflag_Aevo_scheme .gt. id_no_evolution) then
        if ( SGS_param%iflag_SGS_uxb .ne. id_SGS_none) then
          msg = 'solving SGS induction needs'
          call check_missing_field_w_msg                                &
     &       (fld, msg, iphys%i_SGS_vp_induct)
        end if
      end if
!
!
      if ( evo_C%iflag_scheme .gt. id_no_evolution) then
        if (SGS_param%iflag_SGS_c_flux .ne. id_SGS_none) then
          msg = 'solving SGS compsition flux needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_SGS_c_flux)
        end if
      end if
!
      if ( fl_prop%iflag_scheme .gt. id_no_evolution) then
        if ( SGS_param%iflag_SGS_m_flux .eq. id_SGS_similarity          &
     &     .and. SGS_param%iflag_dynamic .ne. id_SGS_DYNAMIC_OFF) then
          msg = 'SGS momentum flux needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_filter_velo)
        end if
!
        if (     SGS_param%iflag_SGS_lorentz .eq. id_SGS_similarity     &
     &     .and. SGS_param%iflag_dynamic .ne. id_SGS_DYNAMIC_OFF) then
          msg = 'SGS Lorentz term needs'
          call check_missing_field_w_msg                                &
     &       (fld, msg, iphys%i_filter_magne)
        end if
      end if
!
      if(SGS_param%iflag_SGS_gravity .gt. id_SGS_none) then
        if(fl_prop%iflag_4_gravity .eq. id_turn_OFF                     &
     &     .and. fl_prop%iflag_4_composit_buo .eq. id_turn_OFF) then
          call calypso_MPI_abort(ierr_fld,                              &
     &       'set one of buoyancy sources')
        end if
        if(fl_prop%iflag_4_gravity .gt. id_turn_OFF) then
          if(SGS_param%iflag_SGS_m_flux.eq.id_SGS_none                  &
     &       .or. SGS_param%iflag_SGS_h_flux.eq.id_SGS_none) then
            call calypso_MPI_abort(ierr_fld,                            &
     &          'Turn on SGS momentum flux and heat flux')
          end if
        end if
        if(fl_prop%iflag_4_composit_buo .gt. id_turn_OFF) then
          if(SGS_param%iflag_SGS_m_flux .eq. id_SGS_none                &
     &       .or. SGS_param%iflag_SGS_c_flux .eq. id_SGS_none) then
              call calypso_MPI_abort(ierr_fld,                          &
     &          'Turn on SGS momentum flux and composition flux')
          end if
        end if
      end if
!
      if ( evo_T%iflag_scheme .gt. id_no_evolution) then
        if    (SGS_param%iflag_SGS_h_flux .eq. id_SGS_similarity        &
     &   .and. SGS_param%iflag_dynamic .ne. id_SGS_DYNAMIC_OFF) then
          msg = 'SGS heat flux needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_filter_temp)
        end if
      end if
!
!
      if (    evo_B%iflag_scheme .gt. id_no_evolution                   &
     &   .or. cd_prop%iflag_Aevo_scheme .gt. id_no_evolution) then
        if    (SGS_param%iflag_SGS_uxb .eq. id_SGS_similarity           &
     &   .and. SGS_param%iflag_dynamic .ne. id_SGS_DYNAMIC_OFF) then
          msg = 'SGS induction needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_filter_velo)
          call check_missing_field_w_msg                                &
     &       (fld, msg, iphys%i_filter_magne)
        end if
      end if
!
!
      if ( cd_prop%iflag_Aevo_scheme .gt. id_no_evolution) then
        if (   cmt_param%iflag_commute .gt. id_SGS_commute_OFF          &
     &   .and. SGS_param%iflag_dynamic .ne. id_SGS_DYNAMIC_OFF) then
          msg = 'filterd A is required for dynamic model'
          call check_missing_field_w_msg(fld, msg, iphys%i_filter_vecp)
        end if
      end if
!
!
      end subroutine check_dependence_4_FEM_SGS
!
! -----------------------------------------------------------------------
!
      subroutine check_dependence_4_SPH_SGS(evo_B, evo_T, evo_C,        &
                SGS_param, fl_prop, iphys, fld)
!
      type(time_evolution_params), intent(in) :: evo_B
      type(time_evolution_params), intent(in) :: evo_T, evo_C
      type(SGS_model_control_params), intent(in) :: SGS_param
      type(fluid_property), intent(in) :: fl_prop
      type(phys_address), intent(in) :: iphys
      type(phys_data), intent(in) :: fld
!
      character(len=kchara) :: msg
!
!
      if (fl_prop%iflag_scheme .gt. id_no_evolution) then
        if ( SGS_param%iflag_SGS_m_flux .ne. id_SGS_none) then
          msg = 'solving SGS momentum flux needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_SGS_inertia)
        end if
!
        if(SGS_param%iflag_SGS_lorentz .ne. id_SGS_none) then
          msg = 'solving SGS lorentz term needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_SGS_Lorentz)
        end if
      end if
!
!
      if ( evo_T%iflag_scheme .gt. id_no_evolution) then
        if ( SGS_param%iflag_SGS_h_flux .ne. id_SGS_none) then
          msg = 'solving SGS heat flux needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_SGS_h_flux)
        end if
      end if
!
!
      if ( evo_B%iflag_scheme .gt. id_no_evolution) then
        if ( SGS_param%iflag_SGS_uxb .ne. id_SGS_none) then
          msg = 'solving SGS magnetic induction needs'
          call check_missing_field_w_msg                                &
     &       (fld, msg, iphys%i_SGS_induction)
          call check_missing_field_w_msg                                &
     &       (fld, msg, iphys%i_SGS_vp_induct)
        end if
      end if
!
!
      if ( evo_C%iflag_scheme .gt. id_no_evolution) then
        if (SGS_param%iflag_SGS_c_flux .ne. id_SGS_none) then
          msg = 'solving SGS compsition flux needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_SGS_c_flux)
        end if
      end if
!
!
      if ( fl_prop%iflag_scheme .gt. id_no_evolution) then
        if    (SGS_param%iflag_SGS_m_flux .eq. id_SGS_similarity        &
     &   .and. SGS_param%iflag_dynamic .ne. id_SGS_DYNAMIC_OFF) then
          msg = 'SGS momentum flux needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_filter_velo)
        end if
!
        if    (SGS_param%iflag_SGS_lorentz .eq. id_SGS_similarity       &
     &   .and. SGS_param%iflag_dynamic .ne. id_SGS_DYNAMIC_OFF) then
          msg = 'SGS Lorentz term needs'
          call check_missing_field_w_msg                                &
     &       (fld, msg, iphys%i_filter_magne)
        end if
      end if
!
      if(SGS_param%iflag_SGS_gravity .gt. id_SGS_none) then
        if(fl_prop%iflag_4_gravity .eq. id_turn_OFF                     &
     &     .and. fl_prop%iflag_4_composit_buo .eq. id_turn_OFF) then
          call calypso_MPI_abort(ierr_fld,                              &
     &       'set one of buoyancy sources')
        end if
        if(fl_prop%iflag_4_gravity .gt. id_turn_OFF) then
          if(SGS_param%iflag_SGS_m_flux.eq.id_SGS_none                  &
     &       .or. SGS_param%iflag_SGS_h_flux.eq.id_SGS_none) then
            call calypso_MPI_abort(ierr_fld,                            &
     &          'Turn on SGS momentum flux and heat flux')
          end if
        end if
        if(fl_prop%iflag_4_composit_buo .gt. id_turn_OFF) then
          if(SGS_param%iflag_SGS_m_flux.eq.id_SGS_none                  &
     &       .or. SGS_param%iflag_SGS_c_flux.eq.id_SGS_none) then
              call calypso_MPI_abort(ierr_fld,                          &
     &          'Turn on SGS momentum flux and composition flux')
          end if
        end if
      end if
!
      if ( evo_T%iflag_scheme .gt. id_no_evolution) then
        if    (SGS_param%iflag_SGS_h_flux .eq. id_SGS_similarity        &
     &   .and. SGS_param%iflag_dynamic .ne. id_SGS_DYNAMIC_OFF) then
          msg = 'SGS heat flux needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_filter_temp)
        end if
      end if
!
!
      if ( evo_B%iflag_scheme .gt. id_no_evolution) then
        if    (SGS_param%iflag_SGS_uxb .eq. id_SGS_similarity           &
     &   .and. SGS_param%iflag_dynamic .ne. id_SGS_DYNAMIC_OFF) then
          msg = 'SGS induction needs'
          call check_missing_field_w_msg(fld, msg, iphys%i_filter_velo)
          call check_missing_field_w_msg                                &
     &       (fld, msg, iphys%i_filter_magne)
        end if
      end if
!
      end subroutine check_dependence_4_SPH_SGS
!
! -----------------------------------------------------------------------
!
      subroutine check_missing_field_w_msg(fld, target_msg, iphys_ref)
!
      type(phys_data), intent(in) :: fld
      character(len=kchara), intent(in) :: target_msg
      integer(kind = kint), intent(in) :: iphys_ref
!
      if(iphys_ref .gt. 0) return
      write(*,*) target_msg, ': ',                                      &
     &    trim(field_name_by_address(fld, iphys_ref))
      call calypso_MPI_abort(ierr_fld,'Stop program.')
!
      end subroutine check_missing_field_w_msg
!
! -----------------------------------------------------------------------
!
      end module check_dependency_for_MHD
