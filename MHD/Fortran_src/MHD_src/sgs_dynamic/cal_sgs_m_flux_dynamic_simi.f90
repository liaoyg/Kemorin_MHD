!cal_sgs_m_flux_dynamic_simi.f90
!      module cal_sgs_m_flux_dynamic_simi
!
!     Written by H. Matsui on Oct. 2005
!     Modified by H. Matsui on Aug., 2007
!
!!      subroutine s_cal_sgs_m_flux_dynamic_simi                        &
!!     &         (iak_sgs_mf, icomp_sgs_mf, nod_comm, node, ele, iphys, &
!!     &          layer_tbl, jac_3d_q, jac_3d_l, rhs_tbl, m_lump,       &
!!     &          fem_wk, f_l, nod_fld)
!!      subroutine cal_sgs_maxwell_dynamic_simi                         &
!!     &        (iak_sgs_lor, icomp_sgs_lor, nod_comm, node, ele, iphys,&
!!     &         layer_tbl, jac_3d_q, jac_3d_l, rhs_tbl, m_lump,        &
!!     &         fem_wk, f_l, nod_fld)
!!        type(communication_table), intent(in) :: nod_comm
!!        type(node_data), intent(in) :: node
!!        type(element_data), intent(in) :: ele
!!        type(phys_address), intent(in) :: iphys
!!        type(layering_tbl), intent(in) :: layer_tbl
!!        type(jacobians_3d), intent(in) :: jac_3d_q, jac_3d_l
!!        type(tables_4_FEM_assembles), intent(in) :: rhs_tbl
!!        type(lumped_mass_matrices), intent(in) :: m_lump
!!        type(work_finite_element_mat), intent(inout) :: fem_wk
!!        type(finite_ele_mat_node), intent(inout) :: f_l
!!        type(phys_data), intent(inout) :: nod_fld
!
      module cal_sgs_m_flux_dynamic_simi
!
      use m_precision
!
      use m_machine_parameter
      use m_control_parameter
      use m_phys_constants
!
      use t_comm_table
      use t_geometry_data
      use t_phys_data
      use t_phys_address
      use t_jacobian_3d
      use t_table_FEM_const
      use t_layering_ele_list
!
      implicit none
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine s_cal_sgs_m_flux_dynamic_simi                          &
     &         (iak_sgs_mf, icomp_sgs_mf, nod_comm, node, ele, iphys,   &
     &          layer_tbl, jac_3d_q, jac_3d_l, rhs_tbl, m_lump,         &
     &          fem_wk, f_l, nod_fld)
!
      use m_SGS_model_coefs
      use reset_dynamic_model_coefs
      use copy_nodal_fields
      use cal_sgs_fluxes_simi
      use cal_filtering_tensors
      use cal_model_diff_coefs
      use int_element_field_2_node
      use cal_similarity_terms
!
      use clear_work_4_dynamic_model
      use cvt_dynamic_scheme_coord
!
      integer(kind = kint), intent(in) :: iak_sgs_mf, icomp_sgs_mf
!
      type(communication_table), intent(in) :: nod_comm
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(phys_address), intent(in) :: iphys
      type(layering_tbl), intent(in) :: layer_tbl
      type(jacobians_3d), intent(in) :: jac_3d_q, jac_3d_l
      type(tables_4_FEM_assembles), intent(in) :: rhs_tbl
      type(lumped_mass_matrices), intent(in) :: m_lump
!
      type(work_finite_element_mat), intent(inout) :: fem_wk
      type(finite_ele_mat_node), intent(inout) :: f_l
      type(phys_data), intent(inout) :: nod_fld
!
!
!
!    reset model coefficients
!
      call reset_tensor_sgs_model_coefs                                 &
     &   (layer_tbl, icomp_sgs_mf, ele%istack_ele_smp)
      call reset_tensor_sgs_nod_m_coefs                                 &
     &   (icomp_sgs_mf, node%istack_nod_smp)
      call s_clear_work_4_dynamic_model(node, iphys, nod_fld)
!
!   similarity model with wider filter
!
      if (iflag_debug.gt.0)                                             &
     &     write(*,*) 'cal_sgs_mf_simi_wide i_wide_fil_velo'
      call cal_sgs_mf_simi_wide(iphys%i_sgs_grad_f,                     &
     &    iphys%i_filter_velo, iphys%i_wide_fil_velo, icomp_sgs_mf,     &
     &    nod_comm, node, nod_fld)
!
!    SGS term by similarity model
!
      if (iflag_debug.gt.0)                                             &
     &     write(*,*) 'cal_sgs_mf_simi iphys%i_SGS_m_flux'
      call cal_sgs_mf_simi(iphys%i_SGS_m_flux, iphys%i_velo,            &
     &    iphys%i_filter_velo, icomp_sgs_mf, nod_comm, node, nod_fld)
!
!    copy to work array
!
       call copy_tensor_component(node, nod_fld,                        &
     &     iphys%i_SGS_m_flux, iphys%i_sgs_simi)
!      call check_nodal_data                                            &
!     &   (my_rank, nod_fld, n_sym_tensor, iphys%i_sgs_simi)
!
!      filtering
!
      call cal_filtered_sym_tensor(nod_comm, node,                      &
     &    iphys%i_sgs_grad, iphys%i_SGS_m_flux, nod_fld)
!
!      call check_nodal_data                                            &
!     &   (my_rank, nod_fld, n_sym_tensor, iphys%i_sgs_grad)
!
!   Change coordinate
!
      call cvt_tensor_dynamic_scheme_coord(node, iphys, nod_fld)
!
!     obtain model coefficient
!
      if (iflag_debug.gt.0)  write(*,*)                                 &
     &    'cal_model_coefs', n_sym_tensor, iak_sgs_mf, icomp_sgs_mf
      call cal_model_coefs(layer_tbl,                                   &
     &    node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,                &
     &    itype_SGS_m_flux_coef, n_sym_tensor,                          &
     &    iak_sgs_mf, icomp_sgs_mf, intg_point_t_evo)
!
      call cal_ele_sym_tensor_2_node                                    &
     &   (node, ele, jac_3d_q, rhs_tbl, m_lump,                         &
     &    ak_sgs(1,icomp_sgs_mf), ak_sgs_nod(1,icomp_sgs_mf),           &
     &    fem_wk, f_l)
!
      end subroutine s_cal_sgs_m_flux_dynamic_simi
!
!  ---------------------------------------------------------------------
!
      subroutine cal_sgs_maxwell_dynamic_simi                           &
     &        (iak_sgs_lor, icomp_sgs_lor, nod_comm, node, ele, iphys,  &
     &         layer_tbl, jac_3d_q, jac_3d_l, rhs_tbl, m_lump,          &
     &         fem_wk, f_l, nod_fld)
!
      use m_SGS_model_coefs
      use reset_dynamic_model_coefs
      use copy_nodal_fields
      use cal_sgs_fluxes_simi
      use cal_filtering_tensors
      use cal_model_diff_coefs
      use int_element_field_2_node
      use cal_similarity_terms
!
      use clear_work_4_dynamic_model
      use cvt_dynamic_scheme_coord
!
      integer(kind = kint), intent(in) :: iak_sgs_lor, icomp_sgs_lor
!
      type(communication_table), intent(in) :: nod_comm
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(phys_address), intent(in) :: iphys
      type(layering_tbl), intent(in) :: layer_tbl
      type(jacobians_3d), intent(in) :: jac_3d_q, jac_3d_l
      type(tables_4_FEM_assembles), intent(in) :: rhs_tbl
      type(lumped_mass_matrices), intent(in) :: m_lump
!
      type(work_finite_element_mat), intent(inout) :: fem_wk
      type(finite_ele_mat_node), intent(inout) :: f_l
      type(phys_data), intent(inout) :: nod_fld
!
!
!
!    reset model coefficients
!
      call reset_tensor_sgs_model_coefs                                 &
     &   (layer_tbl, icomp_sgs_lor, ele%istack_ele_smp)
      call reset_tensor_sgs_nod_m_coefs                                 &
     &   (icomp_sgs_lor, node%istack_nod_smp)
      call s_clear_work_4_dynamic_model(node, iphys, nod_fld)
!
!   similarity model with wider filter
!
      if (iflag_debug.gt.0)                                             &
     &     write(*,*) 'cal_sgs_mf_simi_wide i_wide_fil_magne'
      call cal_sgs_mf_simi_wide(iphys%i_sgs_grad_f,                     &
     &    iphys%i_filter_magne, iphys%i_wide_fil_magne, icomp_sgs_lor,  &
     &    nod_comm, node, nod_fld)
!
!      call check_nodal_data                                            &
!     &   (my_rank, nod_fld, n_sym_tensor, iphys%i_sgs_grad_f)
!
!    SGS term by similarity model
!
      if (iflag_debug.gt.0)                                             &
     &     write(*,*) 'cal_sgs_mf_simi iphys%i_SGS_maxwell'
      call cal_sgs_mf_simi(iphys%i_SGS_maxwell, iphys%i_magne,          &
     &    iphys%i_filter_magne, icomp_sgs_lor,                          &
     &    nod_comm, node, nod_fld)
!
!    copy to work array
!
       call copy_tensor_component(node, nod_fld,                        &
     &     iphys%i_SGS_maxwell, iphys%i_sgs_simi)
!
!    filtering
!
      call cal_filtered_sym_tensor(nod_comm, node,                      &
     &    iphys%i_sgs_grad, iphys%i_SGS_maxwell, nod_fld)
!
!   Change coordinate
!
      call cvt_tensor_dynamic_scheme_coord(node, iphys, nod_fld)
!
!     obtain model coefficient
!
      if (iflag_debug.gt.0)  write(*,*)                                 &
     &   'cal_model_coefs', n_sym_tensor, iak_sgs_lor, icomp_sgs_lor
      call cal_model_coefs(layer_tbl,                                   &
     &    node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,                &
     &    itype_SGS_maxwell_coef, n_sym_tensor,                         &
     &    iak_sgs_lor, icomp_sgs_lor, intg_point_t_evo)
!
      call cal_ele_sym_tensor_2_node                                    &
     &   (node, ele, jac_3d_q, rhs_tbl, m_lump,                         &
     &    ak_sgs(1,icomp_sgs_lor), ak_sgs_nod(1,icomp_sgs_lor),         &
     &    fem_wk, f_l)
!
      end subroutine cal_sgs_maxwell_dynamic_simi
!
!  ---------------------------------------------------------------------
!
      end module cal_sgs_m_flux_dynamic_simi
