!>@file   interpolate_imark_edge.f90
!!@brief  module interpolate_imark_edge
!!
!!@author H. Matsui
!!@date Programmed in Sep., 2006
!
!>@brief  interpolation on edge
!!
!!@verbatim
!!      subroutine s_interpolate_imark_edge2(np_smp, numnod, numele,    &
!!     &          ie, imark_org, istack_smp, num_points, iele_gauss,    &
!!     &          iedge_gauss, imark)
!!      subroutine s_interpolate_imark_edge3(np_smp, numnod, numele,    &
!!     &          ie, imark_org, istack_smp, num_points, iele_gauss,    &
!!     &          iedge_gauss, imark)
!!
!!      subroutine s_interpolate_i8mark_edge2(np_smp, numnod, numele,   &
!!     &          ie, i8mark_org, istack_smp, num_points, iele_gauss,   &
!!     &          iedge_gauss, i8mark)
!!      subroutine s_interpolate_i8mark_edge3(np_smp, numnod, numele,   &
!!     &          ie, i8mark_org, istack_smp, num_points, iele_gauss,   &
!!     &          iedge_gauss, i8mark)
!!@endverbatim
!
      module interpolate_imark_edge
!
      use m_precision
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine s_interpolate_imark_edge2(np_smp, numnod, numele,      &
     &          ie, imark_org, istack_smp, num_points, iele_gauss,      &
     &          iedge_gauss, imark)
!
      use m_geometry_constants
!
      integer (kind = kint), intent(in) :: np_smp
      integer (kind = kint), intent(in) :: numnod, numele
      integer (kind = kint), intent(in) :: ie(numele,8)
      integer (kind = kint), intent(in) :: istack_smp(0:np_smp)
      integer (kind = kint), intent(in) :: num_points
      integer (kind = kint), intent(in) :: iele_gauss(num_points)
      integer (kind = kint), intent(in) :: iedge_gauss(num_points)
      integer (kind = kint), intent(in) :: imark_org(numnod)
!
      integer (kind = kint), intent(inout) :: imark(num_points)
!
      integer (kind = kint) :: ip, ist, ied
      integer (kind = kint) :: iele, i1, iedge, k1, ig
!
!
!$omp parallel do private(ist,ied,ig,iele,iedge,k1,i1)
      do ip = 1, np_smp
        ist = istack_smp(ip-1) + 1
        ied = istack_smp(ip)
        do ig = ist, ied
!
          iele =  iele_gauss(ig)
          iedge = mod(iedge_gauss(ig),100)
!
          k1 = node_on_edge_l(1,iedge)
          i1 = ie(iele,k1)
!
          imark(ig) =  imark_org(i1)
        end do
      end do
!$omp end parallel do
!
      end subroutine s_interpolate_imark_edge2
!
! ----------------------------------------------------------------------
!
      subroutine s_interpolate_imark_edge3(np_smp, numnod, numele,      &
     &          ie, imark_org, istack_smp, num_points, iele_gauss,      &
     &          iedge_gauss, imark)
!
      use m_geometry_constants
!
      integer (kind = kint), intent(in) :: np_smp
      integer (kind = kint), intent(in) :: numnod, numele
      integer (kind = kint), intent(in) :: ie(numele,20)
      integer (kind = kint), intent(in) :: istack_smp(0:np_smp)
      integer (kind = kint), intent(in) :: num_points
      integer (kind = kint), intent(in) :: iele_gauss(num_points)
      integer (kind = kint), intent(in) :: iedge_gauss(num_points)
      integer (kind = kint), intent(in) :: imark_org(numnod)
!
      integer (kind = kint), intent(inout) :: imark(num_points)
!
!
      integer (kind = kint) :: ip, ist, ied
      integer (kind = kint) :: iele, i2, iedge, k2, ig
!
!
!$omp parallel do private(ist,ied,ig,iele,iedge, k2, i2)
      do ip = 1, np_smp
        ist = istack_smp(ip-1) + 1
        ied = istack_smp(ip)
        do ig = ist, ied
!
          iele =  iele_gauss(ig)
          iedge = mod(iedge_gauss(ig),100)
!
          k2 = node_on_edge_q(2,iedge)
          i2 = ie(iele,k2)
!
          imark(ig) = imark_org(i2 )
        end do
      end do
!$omp end parallel do
!
      end subroutine s_interpolate_imark_edge3
!
! ----------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine s_interpolate_i8mark_edge2(np_smp, numnod, numele,     &
     &          ie, i8mark_org, istack_smp, num_points, iele_gauss,     &
     &          iedge_gauss, i8mark)
!
      use m_geometry_constants
!
      integer (kind = kint), intent(in) :: np_smp
      integer (kind = kint), intent(in) :: numnod, numele
      integer (kind = kint), intent(in) :: ie(numele,8)
      integer (kind = kint), intent(in) :: istack_smp(0:np_smp)
      integer (kind = kint), intent(in) :: num_points
      integer (kind = kint), intent(in) :: iele_gauss(num_points)
      integer (kind = kint), intent(in) :: iedge_gauss(num_points)
      integer (kind = kint_d), intent(in) :: i8mark_org(numnod)
!
      integer (kind = kint_d), intent(inout) :: i8mark(num_points)
!
      integer (kind = kint) :: ip, ist, ied
      integer (kind = kint) :: iele, i1, iedge, k1, ig
!
!
!$omp parallel do private(ist,ied,ig,iele,iedge,k1,i1)
      do ip = 1, np_smp
        ist = istack_smp(ip-1) + 1
        ied = istack_smp(ip)
        do ig = ist, ied
!
          iele =  iele_gauss(ig)
          iedge = mod(iedge_gauss(ig),100)
!
          k1 = node_on_edge_l(1,iedge)
          i1 = ie(iele,k1)
!
          i8mark(ig) =  i8mark_org(i1)
        end do
      end do
!$omp end parallel do
!
      end subroutine s_interpolate_i8mark_edge2
!
! ----------------------------------------------------------------------
!
      subroutine s_interpolate_i8mark_edge3(np_smp, numnod, numele,     &
     &          ie, i8mark_org, istack_smp, num_points, iele_gauss,     &
     &          iedge_gauss, i8mark)
!
      use m_geometry_constants
!
      integer (kind = kint), intent(in) :: np_smp
      integer (kind = kint), intent(in) :: numnod, numele
      integer (kind = kint), intent(in) :: ie(numele,20)
      integer (kind = kint), intent(in) :: istack_smp(0:np_smp)
      integer (kind = kint), intent(in) :: num_points
      integer (kind = kint), intent(in) :: iele_gauss(num_points)
      integer (kind = kint), intent(in) :: iedge_gauss(num_points)
      integer (kind = kint_d), intent(in) :: i8mark_org(numnod)
!
      integer (kind = kint_d), intent(inout) :: i8mark(num_points)
!
!
      integer (kind = kint) :: ip, ist, ied
      integer (kind = kint) :: iele, i2, iedge, k2, ig
!
!
!$omp parallel do private(ist,ied,ig,iele,iedge, k2, i2)
      do ip = 1, np_smp
        ist = istack_smp(ip-1) + 1
        ied = istack_smp(ip)
        do ig = ist, ied
!
          iele =  iele_gauss(ig)
          iedge = mod(iedge_gauss(ig),100)
!
          k2 = node_on_edge_q(2,iedge)
          i2 = ie(iele,k2)
!
          i8mark(ig) = i8mark_org(i2 )
        end do
      end do
!$omp end parallel do
!
      end subroutine s_interpolate_i8mark_edge3
!
! ----------------------------------------------------------------------
!
      end module interpolate_imark_edge
