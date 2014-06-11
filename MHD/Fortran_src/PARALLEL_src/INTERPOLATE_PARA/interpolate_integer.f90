!>@file   interpolate_integer.f90
!!@brief  module interpolate_integer
!!
!!@author H. Matsui
!!@date Programmed in Sep., 2006
!
!>@brief  interpolation for integers
!!
!!@verbatim
!!      subroutine s_interpolate_integer(i_dest, i_origin)
!!      subroutine interpolate_integer8(i8_vector_dest, i8_vector_org)
!!@endverbatim
!!
!!@n @param  i_dest      Output integer
!!@n @param  i_origin    Orinal integer data
!
      module interpolate_integer
!
      use m_precision
      use m_constants

      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine s_interpolate_integer(i_vector_dest, i_vector_org)
!
      use calypso_mpi
      use m_geometry_parameter
      use m_machine_parameter
      use m_geometry_data
      use m_2nd_pallalel_vector
      use m_2nd_nod_comm_table
      use m_2nd_geometry_param
      use m_interpolate_table_orgin
      use m_interpolate_table_dest
!
      use m_array_for_send_recv
      use m_work_4_interpolation
!
      use interpolate_imark_1pe
      use select_calypso_SR
      use solver_SR_int
!
!
      integer(kind = kint), intent(in) :: i_vector_org(numnod)
      integer(kind = kint), intent(inout) :: i_vector_dest(nnod_2nd)
!
!     initialize
!
      call verify_2nd_iccg_int_mat(nnod_2nd)
!
      call verifty_work_4_itp_int(ntot_table_org)
!
!
      ix_vec(1:numnod) = i_vector_org(1:numnod)
!
!    interpolation
!
      if (num_dest_domain.gt.0) then
        call s_interporate_imark_para(np_smp, numnod, numele,           &
     &    nnod_4_ele, ie, ix_vec(1), istack_tbl_type_org_smp,           &
     &    ntot_table_org, iele_org_4_org,                               &
     &    itype_inter_org, i_inter_org(1) )
      end if
!
!
!     communication
!
      call sel_calypso_send_recv_int                                    &
     &          (iflag_import_item, ntot_table_org, nnod_2nd,           &
     &           num_dest_domain, iflag_self_itp_send,                  &
     &           id_dest_domain, istack_nod_tbl_org, inod_itp_send,     &
     &           num_org_domain, iflag_self_itp_recv,                   &
     &           id_org_domain, istack_nod_tbl_dest,                    &
     &           inod_dest_4_dest, irev_dest_4_dest,                    &
     &           i_inter_org(1), ivec_2nd(1) )
!
!
      if (num_neib_2.gt.0) then
        call solver_send_recv_i                                         &
     &                (nnod_2nd, num_neib_2, id_neib_2,                 &
     &                 istack_import_2, item_import_2,                  &
     &                 istack_export_2, item_export_2,                  &
     &                 ivec_2nd(1) )
      end if
!
      i_vector_dest(1:nnod_2nd) = ivec_2nd(1:nnod_2nd)
!
      end subroutine s_interpolate_integer
!
! ----------------------------------------------------------------------
!
      subroutine interpolate_integer8(i8_vector_dest, i8_vector_org)
!
      use calypso_mpi
      use m_geometry_parameter
      use m_machine_parameter
      use m_geometry_data
      use m_2nd_pallalel_vector
      use m_2nd_nod_comm_table
      use m_2nd_geometry_param
      use m_interpolate_table_orgin
      use m_interpolate_table_dest
!
      use m_array_for_send_recv
      use m_work_4_interpolation
!
      use interpolate_imark_1pe
      use select_calypso_SR
      use solver_SR_int
!
!
      integer(kind = kint_d), intent(in) :: i8_vector_org(numnod)
      integer(kind = kint_d), intent(inout) :: i8_vector_dest(nnod_2nd)
!
!     initialize
!
      call verify_2nd_iccg_int8_mat(nnod_2nd)
!
      call verifty_work_4_itp_int8(ntot_table_org)
!
!
      i8x_vec(1:numnod) = i8_vector_org(1:numnod)
!
!    interpolation
!
      if (num_dest_domain.gt.0) then
        call s_interporate_i8mark_para(np_smp, numnod, numele,          &
     &    nnod_4_ele, ie, i8x_vec(1), istack_tbl_type_org_smp,          &
     &    ntot_table_org, iele_org_4_org,                               &
     &    itype_inter_org, i8_inter_org(1) )
      end if
!
!
!     communication
!
      call sel_calypso_send_recv_int8                                   &
     &          (iflag_import_item, ntot_table_org, nnod_2nd,           &
     &           num_dest_domain, iflag_self_itp_send,                  &
     &           id_dest_domain, istack_nod_tbl_org, inod_itp_send,     &
     &           num_org_domain, iflag_self_itp_recv,                   &
     &           id_org_domain, istack_nod_tbl_dest,                    &
     &           inod_dest_4_dest, irev_dest_4_dest,                    &
     &           i8_inter_org(1), i8vec_2nd(1) )
!
!
      if (num_neib_2.gt.0) then
        call solver_send_recv_i8                                        &
     &                (nnod_2nd, num_neib_2, id_neib_2,                 &
     &                 istack_import_2, item_import_2,                  &
     &                 istack_export_2, item_export_2,                  &
     &                 i8vec_2nd(1) )
      end if
!
      i8_vector_dest(1:nnod_2nd) = i8vec_2nd(1:nnod_2nd)
!
      end subroutine interpolate_integer8
!
! ----------------------------------------------------------------------
!
      end module interpolate_integer
