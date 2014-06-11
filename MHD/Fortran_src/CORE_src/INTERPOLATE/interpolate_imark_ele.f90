!>@file   interpolate_imark_ele.f90
!!@brief  module interpolate_imark_ele
!!
!!@author H. Matsui
!!@date Programmed in Sep., 2006
!
!>@brief  interpolation in element
!!
!!@verbatim
!!      subroutine s_interpolate_imark_ele8(np_smp, numnod, numele,     &
!!     &          ie, imark_org, istack_smp, num_points, iele_gauss,    &
!!     &          imark)
!!      subroutine s_interpolate_imark_ele20(np_smp, numnod, numele,    &
!!     &          ie, imark_org, istack_smp, num_points, iele_gauss,    &
!!     &          imark)
!!      subroutine s_interpolate_imark_ele27(np_smp, numnod, numele,    &
!!     &          ie, imark_org, istack_smp, num_points, iele_gauss,    &
!!     &          imark)
!!
!!      subroutine s_interpolate_i8mark_ele8(np_smp, numnod, numele,    &
!!     &          ie, i8mark_org, istack_smp, num_points, iele_gauss,   &
!!     &          i8mark)
!!      subroutine s_interpolate_i8mark_ele20(np_smp, numnod, numele,   &
!!     &          ie, i8mark_org, istack_smp, num_points, iele_gauss,   &
!!     &          i8mark)
!!      subroutine s_interpolate_i8mark_ele27(np_smp, numnod, numele,   &
!!     &          ie, i8mark_org, istack_smp, num_points, iele_gauss,   &
!!     &          i8mark)
!!@endverbatim
!
      module interpolate_imark_ele
!
      use m_precision
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine s_interpolate_imark_ele8(np_smp, numnod, numele,       &
     &          ie, imark_org, istack_smp, num_points, iele_gauss,      &
     &          imark)
!
      integer (kind = kint), intent(in) :: np_smp
      integer (kind = kint), intent(in) :: numnod, numele
      integer (kind = kint), intent(in) :: ie(numele,8)
      integer (kind = kint), intent(in) :: istack_smp(0:np_smp)
      integer (kind = kint), intent(in) :: num_points
      integer (kind = kint), intent(in) :: iele_gauss(num_points)
      integer (kind = kint), intent(in) :: imark_org(numnod)
!
      integer (kind = kint), intent(inout) :: imark(num_points)
!
      integer (kind = kint) :: ip, ist, ied
      integer (kind = kint) :: iele, i1, ig
!
!
!$omp parallel do private(ist,ied,ig,iele,i1)
      do ip = 1, np_smp
        ist = istack_smp(ip-1) + 1
        ied = istack_smp(ip)
        do ig = ist, ied
!
          iele = iele_gauss(ig)
          i1 = ie(iele,1)
!
          imark(ig) = imark_org(i1 )
        end do
      end do
!$omp end parallel do
!
      end subroutine s_interpolate_imark_ele8
!
! ----------------------------------------------------------------------
!
      subroutine s_interpolate_imark_ele20(np_smp, numnod, numele,      &
     &          ie, imark_org, istack_smp, num_points, iele_gauss,      &
     &          imark)
!
      integer (kind = kint), intent(in) :: np_smp
      integer (kind = kint), intent(in) :: numnod, numele
      integer (kind = kint), intent(in) :: ie(numele,20)
      integer (kind = kint), intent(in) :: istack_smp(0:np_smp)
      integer (kind = kint), intent(in) :: num_points
      integer (kind = kint), intent(in) :: iele_gauss(num_points)
      integer (kind = kint), intent(in) :: imark_org(numnod)
!
      integer (kind = kint), intent(inout) :: imark(num_points)
!
      integer (kind = kint) :: ip, ist, ied
      integer (kind = kint) :: iele, i1
      integer (kind = kint) :: ig
!
!
!$omp parallel do private(ist,ied,ig,iele,i1)
      do ip = 1, np_smp
        ist = istack_smp(ip-1) + 1
        ied = istack_smp(ip)
        do ig = ist, ied
!
          iele = iele_gauss(ig)
          i1 = ie(iele,1)
!
          imark(ig)   =  imark_org(i1)
        end do
      end do
!$omp end parallel do
!
      end subroutine s_interpolate_imark_ele20
!
! ----------------------------------------------------------------------
!
      subroutine s_interpolate_imark_ele27(np_smp, numnod, numele,      &
     &          ie, imark_org, istack_smp, num_points, iele_gauss,      &
     &          imark)
!
      integer (kind = kint), intent(in) :: np_smp
      integer (kind = kint), intent(in) :: numnod, numele
      integer (kind = kint), intent(in) :: ie(numele,27)
      integer (kind = kint), intent(in) :: istack_smp(0:np_smp)
      integer (kind = kint), intent(in) :: num_points
      integer (kind = kint), intent(in) :: iele_gauss(num_points)
      integer (kind = kint), intent(in) :: imark_org(numnod)
!
      integer (kind = kint), intent(inout) :: imark(num_points)
!
      integer (kind = kint) :: ip, ist, ied
      integer (kind = kint) :: iele, i27
      integer (kind = kint) :: ig
!
!
!$omp parallel do private(ist,ied,ig,iele,i27)
      do ip = 1, np_smp
        ist = istack_smp(ip-1) + 1
        ied = istack_smp(ip)
        do ig = ist, ied
!
          iele = iele_gauss(ig)
          i27 = ie(iele,27)
!
          imark(ig) = imark_org(i27)
        end do
      end do
!$omp end parallel do
!
      end subroutine s_interpolate_imark_ele27
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine s_interpolate_i8mark_ele8(np_smp, numnod, numele,      &
     &          ie, i8mark_org, istack_smp, num_points, iele_gauss,     &
     &          i8mark)
!
      integer (kind = kint), intent(in) :: np_smp
      integer (kind = kint), intent(in) :: numnod, numele
      integer (kind = kint), intent(in) :: ie(numele,8)
      integer (kind = kint), intent(in) :: istack_smp(0:np_smp)
      integer (kind = kint), intent(in) :: num_points
      integer (kind = kint), intent(in) :: iele_gauss(num_points)
      integer (kind = kint_d), intent(in) :: i8mark_org(numnod)
!
      integer (kind = kint_d), intent(inout) :: i8mark(num_points)
!
      integer (kind = kint) :: ip, ist, ied
      integer (kind = kint) :: iele, i1, ig
!
!
!$omp parallel do private(ist,ied,ig,iele,i1)
      do ip = 1, np_smp
        ist = istack_smp(ip-1) + 1
        ied = istack_smp(ip)
        do ig = ist, ied
!
          iele = iele_gauss(ig)
          i1 = ie(iele,1)
!
          i8mark(ig) = i8mark_org(i1 )
        end do
      end do
!$omp end parallel do
!
      end subroutine s_interpolate_i8mark_ele8
!
! ----------------------------------------------------------------------
!
      subroutine s_interpolate_i8mark_ele20(np_smp, numnod, numele,     &
     &          ie, i8mark_org, istack_smp, num_points, iele_gauss,     &
     &          i8mark)
!
      integer (kind = kint), intent(in) :: np_smp
      integer (kind = kint), intent(in) :: numnod, numele
      integer (kind = kint), intent(in) :: ie(numele,20)
      integer (kind = kint), intent(in) :: istack_smp(0:np_smp)
      integer (kind = kint), intent(in) :: num_points
      integer (kind = kint), intent(in) :: iele_gauss(num_points)
      integer (kind = kint_d), intent(in) :: i8mark_org(numnod)
!
      integer (kind = kint_d), intent(inout) :: i8mark(num_points)
!
      integer (kind = kint) :: ip, ist, ied
      integer (kind = kint) :: iele, i1
      integer (kind = kint) :: ig
!
!
!$omp parallel do private(ist,ied,ig,iele,i1)
      do ip = 1, np_smp
        ist = istack_smp(ip-1) + 1
        ied = istack_smp(ip)
        do ig = ist, ied
!
          iele = iele_gauss(ig)
          i1 = ie(iele,1)
!
          i8mark(ig)   =  i8mark_org(i1)
        end do
      end do
!$omp end parallel do
!
      end subroutine s_interpolate_i8mark_ele20
!
! ----------------------------------------------------------------------
!
      subroutine s_interpolate_i8mark_ele27(np_smp, numnod, numele,     &
     &          ie, i8mark_org, istack_smp, num_points, iele_gauss,     &
     &          i8mark)
!
      integer (kind = kint), intent(in) :: np_smp
      integer (kind = kint), intent(in) :: numnod, numele
      integer (kind = kint), intent(in) :: ie(numele,27)
      integer (kind = kint), intent(in) :: istack_smp(0:np_smp)
      integer (kind = kint), intent(in) :: num_points
      integer (kind = kint), intent(in) :: iele_gauss(num_points)
      integer (kind = kint_d), intent(in) :: i8mark_org(numnod)
!
      integer (kind = kint_d), intent(inout) :: i8mark(num_points)
!
      integer (kind = kint) :: ip, ist, ied
      integer (kind = kint) :: iele, i27
      integer (kind = kint) :: ig
!
!
!$omp parallel do private(ist,ied,ig,iele,i27)
      do ip = 1, np_smp
        ist = istack_smp(ip-1) + 1
        ied = istack_smp(ip)
        do ig = ist, ied
!
          iele = iele_gauss(ig)
          i27 = ie(iele,27)
!
          i8mark(ig) = i8mark_org(i27)
        end do
      end do
!$omp end parallel do
!
      end subroutine s_interpolate_i8mark_ele27
!
! ----------------------------------------------------------------------
!
      end module interpolate_imark_ele
