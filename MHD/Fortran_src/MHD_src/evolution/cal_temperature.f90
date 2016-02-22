!
!      module cal_temperature
!
!      subroutine cal_temperature_field
!
!        programmed by H.Matsui and H.Okuda
!                                    on July 2000 (ver 1.1)
!        modieied by H. Matsui on Sep., 2005
!
      module cal_temperature
!
      use m_precision
!
      use m_machine_parameter
      use m_phys_constants
      use m_control_parameter
      use m_t_int_parameter
      use m_nod_comm_table
      use m_geometry_data
      use m_geometry_data_MHD
      use m_node_phys_data
      use m_element_phys_data
      use m_jacobians
      use m_jacobian_sf_grp
      use m_element_id_4_node
      use m_finite_element_matrix
      use m_int_vol_data
      use m_filter_elength
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine cal_temperature_field
!
      use m_group_data
      use m_bc_data_ene
      use m_SGS_address
!
      use nod_phys_send_recv
      use cal_sgs_fluxes
      use set_boundary_scalars
      use int_vol_diffusion_ele
      use int_surf_temp
      use int_vol_thermal_ele
      use cal_stratification_by_temp
      use copy_nodal_fields
      use evolve_by_1st_euler
      use evolve_by_adams_bashforth
      use evolve_by_lumped_crank
      use evolve_by_consist_crank
!
!      call check_surface_param_smp('cal_velocity_pre start',           &
!     &    my_rank, sf_grp1, sf_grp_nod1)
!      call check_jacobians_triquad(ele1, jac1_3d_q)
!
      if (iflag_SGS_heat .ne. id_SGS_none) then
        call cal_sgs_heat_flux(icomp_sgs_hf, ie_dvx,                    &
     &      nod_comm, node1, ele1, fluid1, iphys, iphys_ele, fld_ele1,  &
     &      jac1_3d_q, rhs_tbl1, FEM1_elen, mhd_fem1_wk, fem1_wk,       &
     &      f1_l, f1_nl, nod_fld1)
      end if
!
!      call check_nodal_data(my_rank, nod_fld1, 3, iphys%i_SGS_h_flux)
!
!  ----------  clear the vector and lumped mass matrix
!
      call reset_ff_smps(node1%max_nod_smp, f1_l, f1_nl)
!
!  ----------  lead diffusion term
!
      if (coef_temp.gt.zero .and. coef_exp_t.gt.zero) then
        call int_vol_scalar_diffuse_ele(fluid1%istack_ele_fld_smp,      &
     &      node1, ele1, nod_fld1, jac1_3d_q, rhs_tbl1, FEM1_elen,      &
     &      iak_diff_t, coef_exp_t, ak_d_temp, iphys%i_temp,            &
     &      fem1_wk, f1_l)
      end if
!
!  ----------  lead advection term
!
      if (iflag_temp_supg .gt. id_turn_OFF) then
        call int_vol_temp_ele_upw                                       &
     &    (node1, ele1, fluid1, iphys, nod_fld1,                        &
     &     fld_ele1%ntot_phys, iphys_ele%i_velo, fld_ele1%d_fld,        &
     &     iak_diff_hf, jac1_3d_q, rhs_tbl1, FEM1_elen, mhd_fem1_wk,    &
     &     fem1_wk, f1_nl)
      else
        call int_vol_temp_ele                                           &
     &    (node1, ele1, fluid1, iphys, nod_fld1,                        &
     &     fld_ele1%ntot_phys, iphys_ele%i_velo, fld_ele1%d_fld,        &
     &     iak_diff_hf, jac1_3d_q, rhs_tbl1, FEM1_elen, mhd_fem1_wk,    &
     &     fem1_wk, f1_nl)
      end if
!
!      call check_ff_smp(my_rank, n_scalar, node1%max_nod_smp, f1_l)
!      call check_ff_smp(my_rank, n_scalar, node1%max_nod_smp, f1_nl)
!
      call int_surf_temp_ele(iak_diff_hf, node1, ele1, surf1, sf_grp1,  &
     &    iphys, nod_fld1, jac1_sf_grp_2d_q, rhs_tbl1, FEM1_elen,       &
     &    fem1_wk, f1_l, f1_nl)
!
!      call check_nodal_data(my_rank, nod_fld1, n_scalar, iphys%i_temp)
!      call check_nodal_data(my_rank, fld_ele1,                         &
!     &    n_vector, iphys_ele%i_velo)
!      call check_ff_smp(my_rank, n_scalar, node1%max_nod_smp, f1_l)
!      call check_ff_smp(my_rank, n_scalar, node1%max_nod_smp, f1_nl)
!
      if (iflag_t_strat .gt. id_turn_OFF) then
        if (iflag_temp_supg .gt. id_turn_OFF) then
          call cal_stratified_layer_upw                                 &
     &       (node1, ele1, fluid1, iphys, nod_fld1,                     &
     &        fld_ele1%ntot_phys, iphys_ele%i_velo, fld_ele1%d_fld,     &
     &        jac1_3d_q, rhs_tbl1, mhd_fem1_wk, fem1_wk, f1_nl)
        else
          call cal_stratified_layer                                     &
     &       (node1, ele1, fluid1, iphys, nod_fld1,                     &
     &        fld_ele1%ntot_phys, iphys_ele%i_velo, fld_ele1%d_fld,     &
     &        jac1_3d_q, rhs_tbl1, mhd_fem1_wk, fem1_wk, f1_nl)
        end if
      end if
!
!
      if (iflag_t_evo_4_temp .eq. id_explicit_euler) then
        call cal_scalar_pre_euler(iflag_temp_supg, iphys%i_temp,        &
     &      nod_comm, node1, ele1, fluid1, iphys_ele, fld_ele1,         &
     &      jac1_3d_q, rhs_tbl1, mhd_fem1_wk, fem1_wk,  f1_l, f1_nl,    &
     &      nod_fld1)
      else if (iflag_t_evo_4_temp .eq. id_explicit_adams2) then
        call cal_scalar_pre_adams                                       &
     &     (iflag_temp_supg, iphys%i_temp, iphys%i_pre_heat,            &
     &      nod_comm, node1, ele1, fluid1, iphys_ele, fld_ele1,         &
     &      jac1_3d_q, rhs_tbl1, mhd_fem1_wk, fem1_wk,  f1_l, f1_nl,    &
     &      nod_fld1)
      else if (iflag_t_evo_4_temp .eq. id_Crank_nicolson) then
        call cal_temp_pre_lumped_crank                                  &
     &     (iphys%i_temp, iphys%i_pre_heat, iak_diff_t, nod_bc1_t,      &
     &      nod_comm, node1, ele1, fluid1, iphys_ele, fld_ele1,         &
     &      jac1_3d_q, rhs_tbl1,  FEM1_elen, mhd_fem1_wk, fem1_wk,      &
     &      f1_l, f1_nl, nod_fld1)
      else if (iflag_t_evo_4_temp .eq. id_Crank_nicolson_cmass) then 
        call cal_temp_pre_consist_crank                                 &
     &     (iphys%i_temp, iphys%i_pre_heat, iak_diff_t, nod_bc1_t,      &
     &      node1, ele1, fluid1, jac1_3d_q, rhs_tbl1, FEM1_elen,        &
     &      mhd_fem1_wk, fem1_wk, f1_l, f1_nl, nod_fld1)
      end if
!
      call set_boundary_scalar(nod_bc1_t, iphys%i_temp, nod_fld1)
!
      call scalar_send_recv(iphys%i_temp, node1, nod_comm, nod_fld1)
!
      if (iphys%i_par_temp .gt. 0) then
        call subtract_2_nod_scalars(node1, nod_fld1,                    &
     &      iphys%i_temp, iphys%i_ref_t, iphys%i_par_temp)
      end if
!
      end subroutine cal_temperature_field
!
! ----------------------------------------------------------------------
!
      end module cal_temperature
