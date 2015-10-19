!
!      module cal_3d_filter_phys_smp
!
!      Written by H. Matsui on Sep., 2007
!
!!      subroutine cal_3d_filter_scalar_phys_smp(num_filter_grp,        &
!!     &          id_filter_grp, i_field, numnod, internal_node,        &
!!     &          ntot_comp, i_filter, d_nod)
!!      subroutine cal_3d_filter_vector_phys_smp(num_filter_grp,        &
!!     &          id_filter_grp, i_field, numnod, internal_node,        &
!!     &          ntot_comp, i_filter, d_nod)
!!      subroutine cal_3d_filter_tensor_phys_smp(num_filter_grp,        &
!!     &          id_filter_grp, i_field, numnod, internal_node,        &
!!     &          ntot_comp, i_filter, d_nod)
!!         i_filter: field ID for filtered field
!!         i_field:  field ID to be filtered
!!         num_filter_grp:  num. of filtereing area
!!         id_filter_grp:   table id for filtering
!!
!!      subroutine cal_3d_ez_filter_scalar_smp(num_filter_grp,          &
!!     &          id_filter_grp, i_field, numnod, internal_node,        &
!!     &          ntot_comp, i_filter, d_nod)
!!      subroutine cal_3d_ez_filter_vector_smp(num_filter_grp,          &
!!     &          id_filter_grp, i_field, numnod, internal_node,        &
!!     &          ntot_comp, i_filter, d_nod)
!!      subroutine cal_3d_ez_filter_tensor_smp(num_filter_grp,          &
!!     &          id_filter_grp, i_field, numnod, internal_node,        &
!!     &          ntot_comp, i_filter, d_nod)
!
      module cal_3d_filter_phys_smp
!
      use m_precision
!
      use m_array_for_send_recv
      use m_3d_filter_coef_smp
      use m_nod_filter_comm_table
      use prepare_field_2_filter
      use send_recv_3d_filtering
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine cal_3d_filter_scalar_phys_smp(num_filter_grp,          &
     &          id_filter_grp, i_field, numnod, internal_node,          &
     &          ntot_comp, i_filter, d_nod)
!
      use sum_3d_filter_phys_smp
!
      integer(kind = kint), intent(in) :: num_filter_grp
      integer(kind = kint), intent(in) :: id_filter_grp(num_filter_grp)
!
      integer(kind = kint), intent(in) :: numnod, internal_node
      integer(kind = kint), intent(in) :: ntot_comp, i_field, i_filter
      real(kind = kreal), intent(inout) :: d_nod(numnod,ntot_comp)
!
!
      call prepare_scalar_2_filter                                      &
     &   (numnod, internal_node, ntot_comp, i_field, d_nod)
!
      call sum_3d_filter_scalar_phys_smp(num_filter_grp, id_filter_grp, &
     &    ngrp_nod_3d_filter_smp, istack_nod_3d_fil_smp,                &
     &    ntot_nod_3d_filter_smp, inod_3d_filter_smp,                   &
     &    istack_near_nod_3d_f_smp, ntot_near_nod_3d_filter_smp,        &
     &    inod_near_nod_3d_smp, filter_weight_3d_smp,                   &
     &    max_nsum_3d_fil_smp, istack_nsum_3d_fil_smp,                  &
     &    ntot_nsum_3d_fil_smp, ist_nsum_3d_fil_smp,                    &
     &    ied_nsum_3d_fil_smp, nnod_filtering, x_vec_filtering(1),      &
     &    numnod, x_vec(1))
!
      call scalar_send_recv_3d_filter(numnod, x_vec(1),                 &
     &    ntot_comp, i_filter, d_nod)
!
      end subroutine cal_3d_filter_scalar_phys_smp
!
! ----------------------------------------------------------------------
!
      subroutine cal_3d_filter_vector_phys_smp(num_filter_grp,          &
     &          id_filter_grp, i_field, numnod, internal_node,          &
     &          ntot_comp, i_filter, d_nod)
!
      use sum_3d_filter_phys_smp
!
      integer(kind = kint), intent(in) :: num_filter_grp
      integer(kind = kint), intent(in) :: id_filter_grp(num_filter_grp)
!
      integer(kind = kint), intent(in) :: numnod, internal_node
      integer(kind = kint), intent(in) :: ntot_comp, i_field, i_filter
      real(kind = kreal), intent(inout) :: d_nod(numnod,ntot_comp)
!
!
      call prepare_vector_2_filter                                      &
     &   (numnod, internal_node, ntot_comp, i_field, d_nod)
!
      call sum_3d_filter_vector_phys_smp(num_filter_grp, id_filter_grp, &
     &    ngrp_nod_3d_filter_smp, istack_nod_3d_fil_smp,                &
     &    ntot_nod_3d_filter_smp, inod_3d_filter_smp,                   &
     &    istack_near_nod_3d_f_smp, ntot_near_nod_3d_filter_smp,        &
     &    inod_near_nod_3d_smp, filter_weight_3d_smp,                   &
     &    max_nsum_3d_fil_smp, istack_nsum_3d_fil_smp,                  &
     &    ntot_nsum_3d_fil_smp, ist_nsum_3d_fil_smp,                    &
     &    ied_nsum_3d_fil_smp, nnod_filtering, x_vec_filtering(1),      &
     &    numnod, x_vec(1))
!
      call vector_send_recv_3d_filter(numnod, x_vec(1),                 &
     &    ntot_comp, i_filter, d_nod)
!
      end subroutine cal_3d_filter_vector_phys_smp
!
! ----------------------------------------------------------------------
!
      subroutine cal_3d_filter_tensor_phys_smp(num_filter_grp,          &
     &          id_filter_grp, i_field, numnod, internal_node,          &
     &          ntot_comp, i_filter, d_nod)
!
      use sum_3d_filter_phys_smp
!
      integer(kind = kint), intent(in) :: num_filter_grp
      integer(kind = kint), intent(in) :: id_filter_grp(num_filter_grp)
!
      integer(kind = kint), intent(in) :: numnod, internal_node
      integer(kind = kint), intent(in) :: ntot_comp, i_field, i_filter
      real(kind = kreal), intent(inout) :: d_nod(numnod,ntot_comp)
!
!
      call prepare_sym_tensor_2_filter                                  &
     &   (numnod, internal_node, ntot_comp, i_field, d_nod)
!
      call sum_3d_filter_tensor_phys_smp(num_filter_grp, id_filter_grp, &
     &    ngrp_nod_3d_filter_smp, istack_nod_3d_fil_smp,                &
     &    ntot_nod_3d_filter_smp, inod_3d_filter_smp,                   &
     &    istack_near_nod_3d_f_smp, ntot_near_nod_3d_filter_smp,        &
     &    inod_near_nod_3d_smp, filter_weight_3d_smp,                   &
     &    max_nsum_3d_fil_smp, istack_nsum_3d_fil_smp,                  &
     &    ntot_nsum_3d_fil_smp, ist_nsum_3d_fil_smp,                    &
     &    ied_nsum_3d_fil_smp, nnod_filtering, x_vec_filtering(1),      &
     &    numnod, x_vec(1))
!
      call tensor_send_recv_3d_filter(numnod, x_vec(1),                 &
     &    ntot_comp, i_filter, d_nod)
!
      end subroutine cal_3d_filter_tensor_phys_smp
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine cal_3d_ez_filter_scalar_smp(num_filter_grp,            &
     &          id_filter_grp, i_field, numnod, internal_node,          &
     &          ntot_comp, i_filter, d_nod)
!
      use sum_3d_ez_filter_phys_smp
!
      integer(kind = kint), intent(in) :: num_filter_grp
      integer(kind = kint), intent(in) :: id_filter_grp(num_filter_grp)
!
      integer(kind = kint), intent(in) :: numnod, internal_node
      integer(kind = kint), intent(in) :: ntot_comp, i_field, i_filter
      real(kind = kreal), intent(inout) :: d_nod(numnod,ntot_comp)
!
!
      call prepare_scalar_2_filter                                      &
     &   (numnod, internal_node, ntot_comp, i_field, d_nod)
!
      call sum_3d_ez_filter_scalar_smp(num_filter_grp, id_filter_grp,   &
     &    ngrp_nod_3d_filter_smp, istack_nod_3d_fil_smp,                &
     &    ntot_nod_3d_filter_smp, inod_3d_filter_smp,                   &
     &    istack_near_nod_3d_f_smp, ntot_near_nod_3d_filter_smp,        &
     &    inod_near_nod_3d_smp, filter_weight_3d_smp,                   &
     &    nnod_filtering, x_vec_filtering(1), numnod, x_vec(1))

!
      call scalar_send_recv_3d_filter(numnod, x_vec(1),                 &
     &    ntot_comp, i_filter, d_nod)
!
      end subroutine cal_3d_ez_filter_scalar_smp
!
! ----------------------------------------------------------------------
!
      subroutine cal_3d_ez_filter_vector_smp(num_filter_grp,            &
     &          id_filter_grp, i_field, numnod, internal_node,          &
     &          ntot_comp, i_filter, d_nod)
!
      use sum_3d_ez_filter_phys_smp
!
      integer(kind = kint), intent(in) :: num_filter_grp
      integer(kind = kint), intent(in) :: id_filter_grp(num_filter_grp)
!
      integer(kind = kint), intent(in) :: numnod, internal_node
      integer(kind = kint), intent(in) :: ntot_comp, i_field, i_filter
      real(kind = kreal), intent(inout) :: d_nod(numnod,ntot_comp)
!
!
      call prepare_vector_2_filter                                      &
     &   (numnod, internal_node, ntot_comp, i_field, d_nod)
!
      call sum_3d_ez_filter_vector_smp(num_filter_grp, id_filter_grp,   &
     &    ngrp_nod_3d_filter_smp, istack_nod_3d_fil_smp,                &
     &    ntot_nod_3d_filter_smp, inod_3d_filter_smp,                   &
     &    istack_near_nod_3d_f_smp, ntot_near_nod_3d_filter_smp,        &
     &    inod_near_nod_3d_smp, filter_weight_3d_smp,                   &
     &    nnod_filtering, x_vec_filtering(1), numnod, x_vec(1))
!
      call vector_send_recv_3d_filter(numnod, x_vec(1),                 &
     &    ntot_comp, i_filter, d_nod)
!
      end subroutine cal_3d_ez_filter_vector_smp
!
! ----------------------------------------------------------------------
!
      subroutine cal_3d_ez_filter_tensor_smp(num_filter_grp,            &
     &          id_filter_grp, i_field, numnod, internal_node,          &
     &          ntot_comp, i_filter, d_nod)
!
      use sum_3d_ez_filter_phys_smp
!
      integer(kind = kint), intent(in) :: num_filter_grp
      integer(kind = kint), intent(in) :: id_filter_grp(num_filter_grp)
!
      integer(kind = kint), intent(in) :: numnod, internal_node
      integer(kind = kint), intent(in) :: ntot_comp, i_field, i_filter
      real(kind = kreal), intent(inout) :: d_nod(numnod,ntot_comp)
!
!
      call prepare_sym_tensor_2_filter                                  &
     &   (numnod, internal_node, ntot_comp, i_field, d_nod)
!
      call sum_3d_ez_filter_tensor_smp(num_filter_grp, id_filter_grp,   &
     &    ngrp_nod_3d_filter_smp, istack_nod_3d_fil_smp,                &
     &    ntot_nod_3d_filter_smp, inod_3d_filter_smp,                   &
     &    istack_near_nod_3d_f_smp, ntot_near_nod_3d_filter_smp,        &
     &    inod_near_nod_3d_smp, filter_weight_3d_smp,                   &
     &    nnod_filtering, x_vec_filtering(1), numnod, x_vec(1))
!
      call tensor_send_recv_3d_filter(numnod, x_vec(1),                 &
     &    ntot_comp, i_filter, d_nod)
!
      end subroutine cal_3d_ez_filter_tensor_smp
!
! ----------------------------------------------------------------------
!
      end module cal_3d_filter_phys_smp
