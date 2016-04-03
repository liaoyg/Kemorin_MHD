!
!      module cal_model_diff_coefs
!
!     Written by H. Matsui on Aug., 2007
!
!!      subroutine cal_model_coefs(layer_tbl,                           &
!!     &          node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,        &
!!     &          itype_csim, n_tensor, ifield_d, icomp_f, n_int,       &
!!     &          sgs_coefs)
!!
!!      subroutine cal_diff_coef(layer_tbl,                             &
!!     &          node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,        &
!!     &          numdir, ifield_d, icomp_f, n_int, diff_coefs)
!!      subroutine cal_diff_coef_fluid(layer_tbl, node, ele, fluid,     &
!!     &          iphys, nod_fld, jac_3d_q, jac_3d_l,                   &
!!     &          numdir, ifield_d, icomp_f, n_int, diff_coefs)
!!      subroutine cal_diff_coef_conduct(layer_tbl, node, ele, conduct, &
!!     &          iphys, nod_fld, jac_3d_q, jac_3d_l,                   &
!!     &          numdir, ifield_d, icomp_f, n_int, diff_coefs)
!!        type(node_data), intent(in) :: node
!!        type(element_data), intent(in) :: ele
!!        type(field_geometry_data), intent(in) :: fluid
!!        type(field_geometry_data), intent(in) :: conduct
!!        type(layering_tbl), intent(in) :: layer_tbl
!!        type(phys_address), intent(in) :: iphys
!!        type(phys_data), intent(in) :: nod_fld
!!        type(jacobians_3d), intent(in) :: jac_3d_q, jac_3d_l
!!        type(MHD_coefficients_type), intent(inout) :: sgs_coefs
!!        type(MHD_coefficients_type), intent(inout) :: diff_coefs
!
!
      module cal_model_diff_coefs
!
      use m_precision
!
      use t_geometry_data_MHD
      use t_geometry_data
      use t_phys_address
      use t_phys_data
      use t_layering_ele_list
      use t_jacobians
      use t_material_property
!
      implicit none
!
      private :: cal_layerd_diff_coef, cal_whole_diff_coef
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine cal_model_coefs(layer_tbl,                             &
     &          node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,          &
     &          itype_csim, n_tensor, ifield_d, icomp_f, n_int,         &
     &          sgs_coefs)
!
      use cal_lsq_model_coefs
      use cal_ave_rms_4_dynamic
      use cal_correlate_4_dynamic
!
      integer (kind = kint), intent(in) :: itype_csim, n_tensor
      integer (kind = kint), intent(in) :: n_int, ifield_d, icomp_f
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(layering_tbl), intent(in) :: layer_tbl
      type(phys_address), intent(in) :: iphys
      type(phys_data), intent(in) :: nod_fld
      type(jacobians_3d), intent(in) :: jac_3d_q, jac_3d_l
!
      type(MHD_coefficients_type), intent(inout) :: sgs_coefs
!
!
      call cal_ave_rms_sgs_dynamic(layer_tbl,                           &
     &    node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,                &
     &    n_tensor, icomp_f, n_int)
!
      call cal_correlate_sgs_dynamic(layer_tbl,                         &
     &    node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,                &
     &    n_tensor, icomp_f, n_int)
!
      call cal_model_coef_4_flux(layer_tbl,                             &
     &    node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,                &
     &    itype_csim, n_tensor, ifield_d, icomp_f, n_int, sgs_coefs)
!
      end subroutine cal_model_coefs
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine cal_diff_coef(layer_tbl,                               &
     &          node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,          &
     &          numdir, ifield_d, icomp_f, n_int, diff_coefs)
!
      use m_control_parameter
!
      integer (kind = kint), intent(in) :: numdir
      integer (kind = kint), intent(in) :: ifield_d, icomp_f, n_int
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(layering_tbl), intent(in) :: layer_tbl
      type(phys_address), intent(in) :: iphys
      type(phys_data), intent(in) :: nod_fld
      type(jacobians_3d), intent(in) :: jac_3d_q, jac_3d_l
!
      type(MHD_coefficients_type), intent(inout) :: diff_coefs
!
!
      if (iset_DIFF_model_coefs .eq. 1) then
        call cal_layerd_diff_coef(layer_tbl,                            &
     &      node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,              &
     &      numdir, ifield_d, icomp_f, n_int, diff_coefs)
      else
        call cal_whole_diff_coef(layer_tbl, ele%istack_ele_smp,         &
     &      node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,              &
     &      numdir, ifield_d, icomp_f, n_int, ele%volume, diff_coefs)
      end if
!
      end subroutine cal_diff_coef
!
!-----------------------------------------------------------------------
!
      subroutine cal_diff_coef_fluid(layer_tbl, node, ele, fluid,       &
     &          iphys, nod_fld, jac_3d_q, jac_3d_l,                     &
     &          numdir, ifield_d, icomp_f, n_int, diff_coefs)
!
      use m_control_parameter
!
      integer (kind = kint), intent(in) :: numdir
      integer (kind = kint), intent(in) :: ifield_d, icomp_f, n_int
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(field_geometry_data), intent(in) :: fluid
      type(layering_tbl), intent(in) :: layer_tbl
      type(phys_address), intent(in) :: iphys
      type(phys_data), intent(in) :: nod_fld
      type(jacobians_3d), intent(in) :: jac_3d_q, jac_3d_l
!
      type(MHD_coefficients_type), intent(inout) :: diff_coefs
!
!
      if (iset_DIFF_model_coefs .eq. 1) then
        call cal_layerd_diff_coef(layer_tbl,                            &
     &      node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,              &
     &      numdir, ifield_d, icomp_f, n_int, diff_coefs)
      else
        call cal_whole_diff_coef(layer_tbl, fluid%istack_ele_fld_smp,   &
     &      node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,              &
     &      numdir, ifield_d, icomp_f, n_int, fluid%volume, diff_coefs)
      end if
!
      end subroutine cal_diff_coef_fluid
!
!-----------------------------------------------------------------------
!
      subroutine cal_diff_coef_conduct(layer_tbl, node, ele, conduct,   &
     &          iphys, nod_fld, jac_3d_q, jac_3d_l,                     &
     &          numdir, ifield_d, icomp_f, n_int, diff_coefs)
!
      use m_control_parameter
!
      integer (kind = kint), intent(in) :: numdir
      integer (kind = kint), intent(in) :: ifield_d, icomp_f, n_int
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(field_geometry_data), intent(in) :: conduct
      type(layering_tbl), intent(in) :: layer_tbl
      type(phys_address), intent(in) :: iphys
      type(phys_data), intent(in) :: nod_fld
      type(jacobians_3d), intent(in) :: jac_3d_q, jac_3d_l
!
      type(MHD_coefficients_type), intent(inout) :: diff_coefs
!
!
      if (iset_DIFF_model_coefs .eq. 1) then
        call cal_layerd_diff_coef(layer_tbl,                            &
     &      node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,              &
     &      numdir, ifield_d, icomp_f, n_int, diff_coefs)
      else
        call cal_whole_diff_coef                                        &
     &     (layer_tbl, conduct%istack_ele_fld_smp,                      &
     &      node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,              &
     &      numdir, ifield_d, icomp_f, n_int, conduct%volume,           &
     &      diff_coefs)
      end if
!
      end subroutine cal_diff_coef_conduct
!
!-----------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine cal_layerd_diff_coef(layer_tbl,                        &
     &          node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,          &
     &          numdir, ifield_d, icomp_f, n_int, diff_coefs)
!
      use cal_lsq_model_coefs
      use cal_ave_rms_4_dynamic
      use cal_correlate_4_dynamic
!
      integer (kind = kint), intent(in) :: numdir
      integer (kind = kint), intent(in) :: n_int, ifield_d, icomp_f
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(layering_tbl), intent(in) :: layer_tbl
      type(phys_address), intent(in) :: iphys
      type(phys_data), intent(in) :: nod_fld
      type(jacobians_3d), intent(in) :: jac_3d_q, jac_3d_l
!
      type(MHD_coefficients_type), intent(inout) :: diff_coefs
!
!
      call cal_ave_rms_diff_layerd(layer_tbl,                           &
     &    node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,                &
     &    numdir, icomp_f, n_int)
      call cal_correlate_diff_dynamic (layer_tbl,                       &
     &    node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,                &
     &    numdir, icomp_f, n_int)
      call cal_lsq_layerd_diff_coef(layer_tbl,                          &
     &    node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,                &
     &    numdir, ifield_d, icomp_f, n_int, diff_coefs)
!
      end subroutine cal_layerd_diff_coef
!
!  ---------------------------------------------------------------------
!
      subroutine cal_whole_diff_coef(layer_tbl, iele_fsmp_stack,        &
     &          node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,          &
     &          numdir, ifield_d, icomp_f, n_int, volume_d, diff_coefs)
!
      use cal_lsq_model_coefs
      use cal_ave_rms_4_dynamic
      use cal_correlate_4_dynamic
!
!
      integer (kind = kint), intent(in) :: numdir
      integer (kind = kint), intent(in) :: n_int, ifield_d, icomp_f
      integer (kind = kint), intent(in) :: iele_fsmp_stack(0:np_smp)
      real(kind = kreal), intent(in) :: volume_d
!
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(layering_tbl), intent(in) :: layer_tbl
      type(phys_address), intent(in) :: iphys
      type(phys_data), intent(in) :: nod_fld
      type(jacobians_3d), intent(in) :: jac_3d_q, jac_3d_l
!
      type(MHD_coefficients_type), intent(inout) :: diff_coefs
!
!
      call cal_ave_rms_diff_area(iele_fsmp_stack,                       &
     &    node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,                &
     &    numdir, icomp_f, n_int, volume_d)
!
      call cal_correlate_diff_area(layer_tbl, iele_fsmp_stack,          &
     &    node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,                &
     &    numdir, icomp_f, n_int)
!
      call cal_lsq_diff_coef(iele_fsmp_stack,                           &
     &    node, ele, iphys, nod_fld, jac_3d_q, jac_3d_l,                &
     &    numdir, ifield_d, icomp_f, n_int, diff_coefs)
!
      end subroutine cal_whole_diff_coef
!
!-----------------------------------------------------------------------
!
      end module cal_model_diff_coefs
