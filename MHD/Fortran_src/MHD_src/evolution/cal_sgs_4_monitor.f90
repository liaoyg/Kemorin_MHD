!cal_sgs_4_monitor.f90
!     module cal_sgs_4_monitor
!
!     Written by H. Matsui
!
!      subroutine cal_sgs_terms_4_monitor
!      subroutine cal_diff_of_sgs_terms
!      subroutine cal_work_4_sgs_terms
!
      module cal_sgs_4_monitor
!
      use m_precision
!
      use m_machine_parameter
      use m_control_parameter
      use m_phys_labels
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine cal_sgs_terms_4_monitor
!
      use m_node_phys_address
      use cal_sgs_fluxes
      use monitor_sgs_terms
!
!
      if (iphys%i_SGS_h_flux .gt. 0) then
        if(iflag_debug.gt.0) write(*,*) 'lead ', trim(fhd_SGS_h_flux)
        call cal_sgs_heat_flux
      end if
!
      if (iphys%i_SGS_m_flux .gt. 0) then
        if(iflag_debug.gt.0) write(*,*) 'lead ', trim(fhd_SGS_m_flux)
        call cal_sgs_momentum_flux
      end if
!
      if (iphys%i_SGS_maxwell .gt. 0) then
        if(iflag_debug.gt.0) write(*,*)                                 &
     &        'lead ', trim(fhd_SGS_maxwell_t)
        call cal_sgs_maxwell
      end if
!
      if (iphys%i_SGS_induct_t .gt. 0) then
        if(iflag_debug.gt.0) write(*,*) 'lead ', trim(fhd_induct_t)
        call cal_sgs_magne_induction
      end if
!
      if (iphys%i_SGS_vp_induct .gt. 0) then
        if(iflag_debug.gt.0) write(*,*)                                 &
     &        'lead ', trim(fhd_SGS_vp_induct)
        call cal_sgs_uxb_2_monitor
      end if
!
      end subroutine cal_sgs_terms_4_monitor
!
!-----------------------------------------------------------------------
!
      subroutine cal_diff_of_sgs_terms
!
      use m_node_phys_address
      use cal_terms_for_heat
      use cal_momentum_terms
      use cal_magnetic_terms
!
!
      if (iphys%i_SGS_div_h_flux .gt. 0) then
        if(iflag_debug.gt.0) write(*,*)                                 &
     &        'lead ', trim(fhd_div_SGS_h_flux)
        call cal_terms_4_heat(iphys%i_SGS_div_h_flux)
      end if
!
      if (iphys%i_SGS_div_m_flux .gt. 0) then
        if(iflag_debug.gt.0) write(*,*)                                 &
     &        'lead ', trim(fhd_div_SGS_m_flux)
        call cal_terms_4_momentum(iphys%i_SGS_div_m_flux)
      end if
!
      if (iphys%i_SGS_Lorentz .gt. 0) then
        if(iflag_debug.gt.0) write(*,*)                                 &
     &        'lead ', trim(fhd_SGS_Lorentz)
        call cal_terms_4_momentum(iphys%i_SGS_Lorentz)
      end if
!
      if (      iphys%i_SGS_induction .gt. 0                            &
     &   .and. iflag_t_evo_4_magne .gt. id_no_evolution) then
        if(iflag_debug.gt.0) write(*,*)                                 &
     &        'lead ', trim(fhd_SGS_induction)
        call cal_terms_4_magnetic(iphys%i_SGS_induction)
      end if
!
!
!      if (iphys%i_SGS_buoyancy .gt. 0) then
!        if(iflag_debug.gt.0) write(*,*)                                &
!     &        'lead ', trim(fhd_SGS_buoyancy)
!         call cal_terms_4_magnetic(iphys%i_SGS_induction)
!      end if
!
!      if (iphys%i_SGS_comp_buo .gt. 0) then
!        if(iflag_debug.gt.0) write(*,*)                                &
!     &        'lead ', trim(fhd_SGS_comp_buo)
!         call cal_terms_4_magnetic(iphys%i_SGS_induction)
!      end if
!
      end subroutine cal_diff_of_sgs_terms
!
!-----------------------------------------------------------------------
!
      subroutine cal_work_4_sgs_terms
!
      use m_geometry_data
      use m_node_phys_address
      use m_node_phys_data
      use m_physical_property
!
      use products_nodal_fields_smp
      use int_sgs_induction
      use sgs_buoyancy_flux
!
!
      if (     iphys%i_SGS_induction .gt. 0                             &
     &   .and. iflag_t_evo_4_vect_p .gt. id_no_evolution) then
        if(iflag_debug.gt.0) write(*,*)                                 &
     &        'lead ', trim(fhd_SGS_induction)
        call int_vol_sgs_induction
      end if
!
!
!$omp parallel
      if (iphys%i_SGS_temp_gen .gt. 0) then
        call cal_phys_product_4_scalar(node1, nod_fld1,                 &
     &      iphys%i_temp, iphys%i_SGS_div_h_flux, iphys%i_SGS_temp_gen)
      end if
!
!
      if (iphys%i_reynolds_wk .gt. 0) then
        call cal_phys_dot_product(node1, nod_fld1,                      &
     &      iphys%i_velo, iphys%i_SGS_div_m_flux, iphys%i_reynolds_wk)
      end if
!
      if (iphys%i_SGS_Lor_wk .gt. 0) then
        call cal_phys_dot_product(node1, nod_fld1,                      &
     &      iphys%i_velo, iphys%i_SGS_Lorentz, iphys%i_SGS_Lor_wk)
      end if
!
      if (iphys%i_SGS_me_gen .gt. 0) then
        call cal_phys_dot_product(node1, nod_fld1,                      &
     &      iphys%i_magne, iphys%i_SGS_induction, iphys%i_SGS_me_gen)
      end if
!$omp end parallel
!
      if (iphys%i_SGS_buo_wk .gt. 0) then
        if(iflag_debug.gt.0) write(*,*)                                 &
     &        'lead ', trim(fhd_SGS_buo_flux)
        call cal_SGS_gravity_flux(coef_buo, iphys%i_SGS_h_flux,         &
     &      iphys%i_SGS_buo_wk)
      end if
!
      if (iphys%i_SGS_comp_buo_wk .gt. 0) then
        if(iflag_debug.gt.0) write(*,*)                                 &
     &        'lead ', trim(fhd_SGS_comp_buo_flux)
        call cal_SGS_gravity_flux(coef_comp_buo,                        &
     &      iphys%i_SGS_h_flux, iphys%i_SGS_comp_buo_wk)
      end if
!
!
!$omp parallel
      if (iphys%i_SGS_Lor_wk_tr .gt. 0) then
          call cal_phys_dot_product(node1, nod_fld1,                    &
     &        iphys%i_filter_velo, iphys%i_SGS_Lor_true,                &
     &        iphys%i_SGS_Lor_wk_tr)
      end if
!
      if (iphys%i_reynolds_wk_tr .gt. 0) then
          call cal_phys_dot_product(node1, nod_fld1,                    &
     &        iphys%i_filter_velo, iphys%i_SGS_div_mf_true,             &
     &        iphys%i_reynolds_wk_tr)
      end if
!
      if (iphys%i_SGS_t_gen_tr .gt. 0) then
          call cal_phys_product_4_scalar(node1, nod_fld1,               &
     &        iphys%i_filter_temp, iphys%i_SGS_div_hf_true,             &
     &        iphys%i_SGS_t_gen_tr)
      end if
!
      if (iphys%i_SGS_me_gen_tr .gt. 0) then
          call cal_phys_dot_product(node1, nod_fld1,                    &
     &        iphys%i_filter_magne, iphys%i_SGS_idct_true,              &
     &        iphys%i_SGS_me_gen_tr)
      end if
!$omp end parallel
!
      end subroutine cal_work_4_sgs_terms
!
!-----------------------------------------------------------------------
!
      end module cal_sgs_4_monitor
