!
!      module m_count_num_surf_vector
!
!      Written by H. Matsui on Sep. 2005
!      Modified by H. Matsui on Feb., 2009
!
!      subroutine count_num_surf_velo
!      subroutine count_num_surf_vect_p
!      subroutine count_num_surf_magne
!      subroutine count_num_surf_current
!
!      subroutine count_num_surf_torque
!      subroutine count_num_surf_grad_vecp
!      subroutine count_num_surf_grad_b
!      subroutine count_num_surf_grad_j
!
      module m_count_num_surf_vector
!
      use m_precision
!
      use m_surface_group
      use m_surf_data_list
      use m_header_4_surface_bc
!
      implicit  none
!
!-----------------------------------------------------------------------
!
      contains 
!
!-----------------------------------------------------------------------
!
      subroutine count_num_surf_velo
!
      use m_surface_group_connect
      use m_surf_data_torque
      use m_boundary_condition_IDs
      use set_surf_vector_id
      use set_stress_free_surf_id
!
!
      call s_count_num_surf_vector                                      &
     &   (num_surf, inod_stack_sf_grp, surf_name,                       &
     &    torque_surf%num_bc, torque_surf%bc_name,                      &
     &    torque_surf%ibc_type, name_svn,                               &
     &    nmax_sf_sgs_velo, ngrp_sf_sgs_velo,                           &
     &    ngrp_sf_fix_vn, nnod_sf_fix_vn)
!
      call count_num_stress_free_surf(num_surf, surf_name,              &
     &    torque_surf%num_bc, torque_surf%bc_name,                      &
     &    torque_surf%ibc_type,                                         &
     &    iflag_surf_free_sph_in, iflag_surf_free_sph_out,              &
     &    ngrp_sf_fr_in, ngrp_sf_fr_out)
!
      end subroutine count_num_surf_velo
!
!-----------------------------------------------------------------------
!
      subroutine count_num_surf_vect_p
!
      use m_surface_group_connect
      use m_surf_data_vector_p
      use m_boundary_condition_IDs
      use set_surf_vector_id
      use set_stress_free_surf_id
!
!
      call s_count_num_surf_vector                                      &
     &    (num_surf, inod_stack_sf_grp, surf_name,                      &
     &     a_potential_surf%num_bc, a_potential_surf%bc_name,           &
     &     a_potential_surf%ibc_type, name_san,                         &
     &     nmax_sf_sgs_vect_p, ngrp_sf_sgs_vect_p,                      &
     &     ngrp_sf_fix_vpn, nnod_sf_fix_vpn)
!
      call count_num_stress_free_surf(num_surf, surf_name,              &
     &    torque_surf%num_bc, torque_surf%bc_name,                      &
     &    torque_surf%ibc_type,                                         &
     &    iflag_surf_qvc_sph_in, iflag_surf_qvc_sph_out,                &
     &    ngrp_sf_qvc_in, ngrp_sf_qvc_out)
!
      end subroutine count_num_surf_vect_p
!
!-----------------------------------------------------------------------
!
      subroutine count_num_surf_magne
!
      use m_surface_group_connect
      use m_surf_data_magne
      use set_surf_vector_id
!
!
      call s_count_num_surf_vector                                      &
     &    (num_surf, inod_stack_sf_grp, surf_name,                      &
     &     magne_surf%num_bc, magne_surf%bc_name, magne_surf%ibc_type,  &
     &     name_sbn, nmax_sf_sgs_magne, ngrp_sf_sgs_magne,              &
     &     ngrp_sf_fix_bn, nnod_sf_fix_bn)
!
      end subroutine count_num_surf_magne
!
!-----------------------------------------------------------------------
!
      subroutine count_num_surf_current
!
      use m_surface_group_connect
      use m_surf_data_current
      use set_surf_vector_id
!
!
      call s_count_num_surf_vector                                      &
     &    (num_surf, inod_stack_sf_grp, surf_name,                      &
     &     current_surf%num_bc, current_surf%bc_name,                   &
     &     current_surf%ibc_type, name_sjn, nmax_sf_sgs_current,        &
     &     ngrp_sf_sgs_current, ngrp_sf_fix_jn, nnod_sf_fix_jn)
!
      end subroutine count_num_surf_current
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine count_num_surf_torque
!
      use m_surf_data_torque
      use set_sf_grad_vector_id
!
!
      call count_num_sf_grad_vector                                     &
     &   (num_surf, surf_istack, surf_name,                             &
     &    torque_surf%num_bc, torque_surf%bc_name,                      &
     &    torque_surf%ibc_type, name_vxg, name_vyg, name_vzg,           &
     &    nmax_sf_fix_tq, nmax_ele_sf_fix_tq, nmax_sf_lead_tq,          &
     &    ngrp_sf_fix_tq, nele_sf_fix_tq, ngrp_sf_lead_tq)
!
      end subroutine count_num_surf_torque
!
!-----------------------------------------------------------------------
!
      subroutine count_num_surf_grad_vecp
!
      use m_surf_data_vector_p
      use set_sf_grad_vector_id
!
!
      call count_num_sf_grad_vector                                     &
     &   (num_surf, surf_istack, surf_name, a_potential_surf%num_bc,    &
     &    a_potential_surf%bc_name, a_potential_surf%ibc_type,          &
     &    name_axg, name_ayg, name_azg, nmax_sf_fix_grad_a,             &
     &    nmax_ele_sf_fix_grad_a, nmax_sf_lead_vect_p,                  &
     &    ngrp_sf_fix_grad_a, nele_sf_fix_grad_a, ngrp_sf_lead_vect_p)
!
      end subroutine count_num_surf_grad_vecp
!
!-----------------------------------------------------------------------
!
      subroutine count_num_surf_grad_b
!
      use m_surf_data_magne
      use set_sf_grad_vector_id
!
!
      call count_num_sf_grad_vector                                     &
     &    (num_surf, surf_istack, surf_name,                            &
     &     magne_surf%num_bc, magne_surf%bc_name, magne_surf%ibc_type,  &
     &     name_bxg, name_byg, name_bzg, nmax_sf_fix_grad_b,            &
     &     nmax_ele_sf_fix_grad_b, nmax_sf_lead_b, ngrp_sf_fix_grad_b,  &
     &     nele_sf_fix_grad_b, ngrp_sf_lead_b)
!
      end subroutine count_num_surf_grad_b
!
!-----------------------------------------------------------------------
!
      subroutine count_num_surf_grad_j
!
      use m_surf_data_current
      use set_sf_grad_vector_id
!
!
      call count_num_sf_grad_vector                                     &
     &    (num_surf, surf_istack, surf_name, current_surf%num_bc,       &
     &     current_surf%bc_name, current_surf%ibc_type,                 &
     &     name_jxg, name_jyg, name_jzg, nmax_sf_fix_grad_j,            &
     &     nmax_ele_sf_fix_grad_j, nmax_sf_lead_j, ngrp_sf_fix_grad_j,  &
     &     nele_sf_fix_grad_j, ngrp_sf_lead_j)
!
      end subroutine count_num_surf_grad_j
!
!-----------------------------------------------------------------------
!
      end module m_count_num_surf_vector
