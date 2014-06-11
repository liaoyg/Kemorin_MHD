!>@file   interpolate_i8mark_1pe.f90
!!@brief  module interpolate_i8mark_1pe
!!
!!@author H. Matsui
!!@date Programmed in Sep., 2006
!
!>@brief  interpolation for each domain
!!
!!@verbatim
!!      subroutine s_interporate_i8mark_para(np_smp, numnod, numele,    &
!!     &          nnod_4_ele, ie, i8mark_org, istack_tbl_wtype_smp,     &
!!     &          num_points, iele_gauss, itype_gauss, i8mark)
!!@endverbatim
!
      module interpolate_i8mark_1pe
!
      use m_precision
!
      use interpolate_on_node
      use interpolate_imark_edge
      use interpolate_imark_surf
      use interpolate_imark_ele
!
      implicit none
!
      private :: interpolate_i8mark_8, interpolate_i8mark_20
      private :: interpolate_i8mark_27
!
! ----------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine s_interporate_i8mark_para(np_smp, numnod, numele,      &
     &          nnod_4_ele, ie, i8mark_org, istack_tbl_wtype_smp,       &
     &          num_points, iele_gauss, itype_gauss, i8mark)
!
      integer (kind = kint), intent(in) :: np_smp
      integer (kind = kint), intent(in) :: nnod_4_ele
      integer (kind = kint), intent(in) :: numnod, numele
      integer (kind = kint), intent(in) :: ie(numele,nnod_4_ele)
      integer (kind = kint), intent(in)                                 &
     &       :: istack_tbl_wtype_smp(0:4*np_smp)
      integer (kind = kint), intent(in) :: num_points
      integer (kind = kint), intent(in) :: iele_gauss(num_points)
      integer (kind = kint), intent(in) :: itype_gauss(num_points)
      integer (kind = kint_d), intent(in) :: i8mark_org(numnod)
!
      integer (kind = kint_d), intent(inout) :: i8mark(num_points)
!
!
      if (nnod_4_ele .eq. 8)then
        call interpolate_i8mark_8(np_smp, numnod, numele,               &
     &      nnod_4_ele, ie, i8mark_org, istack_tbl_wtype_smp,           &
     &      num_points, iele_gauss, itype_gauss, i8mark)
      else if (nnod_4_ele .eq. 20)then
        call interpolate_i8mark_20(np_smp, numnod, numele,              &
     &      nnod_4_ele, ie, i8mark_org, istack_tbl_wtype_smp,           &
     &      num_points, iele_gauss, itype_gauss, i8mark)
      else if (nnod_4_ele .eq. 27)then
        call interpolate_i8mark_27(np_smp, numnod, numele,              &
     &      nnod_4_ele, ie, i8mark_org, istack_tbl_wtype_smp,           &
     &      num_points, iele_gauss, itype_gauss, i8mark)
      end if
!
      end subroutine s_interporate_i8mark_para
!
!  ---------------------------------------------------------------------
!
      subroutine interpolate_i8mark_8(np_smp, numnod, numele,           &
     &          nnod_4_ele, ie, i8mark_org, istack_wtype_smp,           &
     &          num_points, iele_gauss, itype_gauss, i8mark)
!
      integer (kind = kint), intent(in) :: np_smp
      integer (kind = kint), intent(in) :: numnod, numele, nnod_4_ele
      integer (kind = kint), intent(in) :: ie(numele,nnod_4_ele)
      integer (kind = kint), intent(in) :: istack_wtype_smp(0:4*np_smp)
      integer (kind = kint), intent(in) :: num_points
      integer (kind = kint), intent(in) :: iele_gauss(num_points)
      integer (kind = kint), intent(in) :: itype_gauss(num_points)
      integer (kind = kint_d), intent(in) :: i8mark_org(numnod)
!
      integer (kind = kint_d), intent(inout) :: i8mark(num_points)
!
      integer(kind = kint) :: ist
!
!
      ist = 0
      call s_interpolate_i8mark_node(np_smp, numnod, numele,            &
     &    nnod_4_ele, ie, i8mark_org, istack_wtype_smp(ist),            &
     &    num_points, iele_gauss, itype_gauss, i8mark)
!
      ist = np_smp
      call s_interpolate_i8mark_edge2(np_smp, numnod, numele, ie,       &
     &    i8mark_org, istack_wtype_smp(ist), num_points, iele_gauss,    &
     &    itype_gauss, i8mark)
!
      ist = 2*np_smp
      call s_interpolate_i8mark_surf4(np_smp, numnod, numele, ie,       &
     &    i8mark_org, istack_wtype_smp(ist), num_points, iele_gauss,    &
     &    itype_gauss, i8mark)
!
      ist = 3*np_smp
      call s_interpolate_i8mark_ele8(np_smp, numnod, numele, ie,        &
     &    i8mark_org, istack_wtype_smp(ist), num_points, iele_gauss,    &
     &    i8mark)
!
      end subroutine interpolate_i8mark_8
!
! ----------------------------------------------------------------------
!
      subroutine interpolate_i8mark_20(np_smp, numnod, numele,          &
     &          nnod_4_ele, ie, i8mark_org, istack_wtype_smp,           &
     &          num_points, iele_gauss, itype_gauss, i8mark)
!
      integer (kind = kint), intent(in) :: np_smp
      integer (kind = kint), intent(in) :: numnod, numele, nnod_4_ele
      integer (kind = kint), intent(in) :: ie(numele,nnod_4_ele)
      integer (kind = kint), intent(in) :: istack_wtype_smp(0:4*np_smp)
      integer (kind = kint), intent(in) :: num_points
      integer (kind = kint), intent(in) :: iele_gauss(num_points)
      integer (kind = kint), intent(in) :: itype_gauss(num_points)
      integer (kind = kint_d), intent(in) :: i8mark_org(numnod)
!
      integer (kind = kint_d), intent(inout) :: i8mark(num_points)
!
      integer(kind = kint) :: ist
!
!
      ist = 0
      call s_interpolate_i8mark_node(np_smp, numnod, numele,            &
     &    nnod_4_ele, ie, i8mark_org, istack_wtype_smp(ist),            &
     &    num_points, iele_gauss, itype_gauss, i8mark)
!
      ist = np_smp
      call s_interpolate_i8mark_edge3(np_smp, numnod, numele, ie,       &
     &    i8mark_org, istack_wtype_smp(ist), num_points, iele_gauss,    &
     &    itype_gauss, i8mark)
!
      ist = 2*np_smp
      call s_interpolate_i8mark_surf8(np_smp, numnod, numele, ie,       &
     &    i8mark_org, istack_wtype_smp(ist), num_points, iele_gauss,    &
     &    itype_gauss, i8mark)
!
      ist = 3*np_smp
      call s_interpolate_i8mark_ele20(np_smp, numnod, numele, ie,       &
     &    i8mark_org, istack_wtype_smp(ist), num_points, iele_gauss,    &
     &    i8mark)
!
      end subroutine interpolate_i8mark_20
!
!  ---------------------------------------------------------------------
!
      subroutine interpolate_i8mark_27(np_smp, numnod, numele,          &
     &          nnod_4_ele, ie, i8mark_org, istack_wtype_smp,           &
     &          num_points, iele_gauss, itype_gauss, i8mark)
!
      integer (kind = kint), intent(in) :: np_smp
      integer (kind = kint), intent(in) :: numnod, numele, nnod_4_ele
      integer (kind = kint), intent(in) :: ie(numele,nnod_4_ele)
      integer (kind = kint), intent(in) :: istack_wtype_smp(0:4*np_smp)
      integer (kind = kint), intent(in) :: num_points
      integer (kind = kint), intent(in) :: iele_gauss(num_points)
      integer (kind = kint), intent(in) :: itype_gauss(num_points)
      integer (kind = kint_d), intent(in) :: i8mark_org(numnod)
!
      integer (kind = kint_d), intent(inout) :: i8mark(num_points)
!
      integer(kind = kint) :: ist
!
!
      ist = 0
      call s_interpolate_i8mark_node(np_smp, numnod, numele,            &
     &    nnod_4_ele, ie, i8mark_org, istack_wtype_smp(ist),            &
     &    num_points, iele_gauss, itype_gauss, i8mark)
!
      ist = np_smp
      call s_interpolate_i8mark_edge3(np_smp, numnod, numele, ie,       &
     &    i8mark_org, istack_wtype_smp(ist), num_points, iele_gauss,    &
     &    itype_gauss, i8mark)
!
      ist = 2*np_smp
      call s_interpolate_i8mark_surf9(np_smp, numnod, numele, ie,       &
     &    i8mark_org, istack_wtype_smp(ist), num_points, iele_gauss,    &
     &    itype_gauss, i8mark)
!
      ist = 3*np_smp
      call s_interpolate_i8mark_ele27(np_smp, numnod, numele, ie,       &
     &    i8mark_org, istack_wtype_smp(ist), num_points, iele_gauss,    &
     &    i8mark)
!
      end subroutine interpolate_i8mark_27
!
!  ---------------------------------------------------------------------
!
      end module interpolate_i8mark_1pe
