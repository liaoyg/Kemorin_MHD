!
!      module cal_electric_potential
!
!        programmed by H.Matsui and H.Okuda
!                                    on July 2000 (ver 1.1)
!        modified by H.Matsui on July, 2006
!
!      subroutine cal_scalar_potential
!
      module cal_electric_potential
!
      use m_precision
!
      use m_geometry_data
      use m_node_phys_data
      use m_jacobians
      use m_element_id_4_node
      use m_finite_element_matrix
      use m_filter_elength
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine cal_scalar_potential
!
      use m_machine_parameter
      use m_group_data
      use m_phys_constants
      use m_SGS_model_coefs
      use m_SGS_address
      use m_int_surface_data
      use m_jacobian_sf_grp
      use m_surf_data_vector_p
      use m_bc_data_magne
      use m_geometry_data
      use m_node_phys_data
      use m_finite_element_matrix
      use m_solver_djds_MHD
      use m_array_for_send_recv
      use m_type_AMG_data
      use m_type_AMG_data_4_MHD
!
      use int_vol_fractional_div
      use int_sk_4_fixed_boundary
      use int_surf_div_sgs
      use set_boundary_scalars
      use int_surf_normal_fields
      use cal_solver_MHD
!
!
      call reset_ff(node1%numnod, f1_l)
      call reset_ff_smps(node1%max_nod_smp, f1_l, f1_nl)
!
      if (iflag_debug .gt. 0)  write(*,*) 'int_vol_divergence_vect_p'
      call int_vol_fractional_div_ele                                   &
     &   (ele1%istack_ele_smp, iphys%i_vecp, iak_diff_b,                &
     &    node1, ele1, nod_fld1, jac1_3d_q, jac1_3d_l,                  &
     &    rhs_tbl1, FEM1_elen, fem1_wk, f1_l)
!
!      if (iflag_commute_magne .eq. id_SGS_commute_ON) then
!        call int_surf_sgs_div_velo_ele                                 &
!     &     (node1, ele1, surf1, sf_grp1, nod_fld1,                     &
!     &      jac1_sf_grp_2d_q, jac1_sf_grp_2d_l, rhs_tbl1, FEM1_elen,   &
!     &      intg_point_poisson, sf_sgs1_grad_a%nmax_sf_dat,            &
!     &      sf_sgs1_grad_a%ngrp_sf_dat, sf_sgs1_grad_a%id_grp_sf_dat,  &
!     &      ifilter_final, ak_diff(1,iak_diff_b),                      &
!     &      iphys%i_vecp, fem1_wk, f1_l)
!      end if
!
      call int_surf_normal_vector_p(iphys%i_vecp, node1, ele1, surf1,   &
     &    sf_grp1, nod_fld1, jac1_sf_grp_2d_l, rhs_tbl1, fem1_wk, f1_l)
!
      call int_vol_sk_mp_bc(iphys%i_m_phi, node1, ele1, nod_fld1,       &
     &    jac1_3d_l, rhs_tbl1, FEM1_elen, fem1_wk, f1_l)
!
      call set_boundary_ff(node1, nod_bc1_f, f1_l)
!
      if (iflag_debug .gt. 0)  write(*,*) 'cal_sol_mag_po'
      call cal_sol_mag_po                                               &
     &         (node1, DJDS_comm_etr, DJDS_linear, Fmat_DJDS,           &
     &          num_MG_level, MG_itp, MG_comm, MG_djds_tbl_l,           &
     &          MG_mat_magp, MG_vector, iphys%i_m_phi, f1_l, b_vec,     &
     &          x_vec, nod_fld1)
!
      if (iflag_debug .gt. 0)  write(*,*) 'set_boundary_m_phi'
      call set_boundary_scalar(nod_bc1_f, iphys%i_m_phi, nod_fld1)
!
      end subroutine cal_scalar_potential
!
!-----------------------------------------------------------------------
!
      end module cal_electric_potential
