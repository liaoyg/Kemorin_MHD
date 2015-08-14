!
!     module set_initial_rotation
!
!     Written by H. Matsui on June, 2005
!
!      subroutine set_initial_velo_1(numnod, xx)
!      subroutine set_initial_velo_2(numnod, xx)
!      subroutine set_initial_velo_3(numnod, xx)
!
      module set_initial_rotation
!
      use m_precision
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine set_initial_velo_1(numnod, xx)
!
      use m_node_phys_address
      use m_node_phys_data
!
      integer(kind = kint), intent(in) :: numnod
      real(kind = kreal), intent(in) :: xx(numnod,3)
!
      integer (kind = kint) :: inod
!
!$omp parallel do
      do inod = 1, numnod
       d_nod(inod,iphys%i_press)  = 0.0d0
       d_nod(inod,iphys%i_velo  ) = 0.0d0
       d_nod(inod,iphys%i_velo+1) = -xx(inod,3)
       d_nod(inod,iphys%i_velo+2) =  xx(inod,2)
      end do
!$omp end parallel do
!
      end subroutine set_initial_velo_1
!
!-----------------------------------------------------------------------
!
      subroutine set_initial_velo_2(numnod, xx)
!
      use m_node_phys_address
      use m_node_phys_data
!
      integer(kind = kint), intent(in) :: numnod
      real(kind = kreal), intent(in) :: xx(numnod,3)
!
      integer (kind = kint) :: inod
!
!$omp parallel do
      do inod = 1, numnod
       d_nod(inod,iphys%i_press)  = 0.0d0
       d_nod(inod,iphys%i_velo  ) =  xx(inod,3)
       d_nod(inod,iphys%i_velo+1) =  0.0d0
       d_nod(inod,iphys%i_velo+2) = -xx(inod,1)
      end do
!$omp end parallel do
!
      end subroutine set_initial_velo_2
!
!-----------------------------------------------------------------------
!
      subroutine set_initial_velo_3(numnod, xx)
!
      use m_node_phys_address
      use m_node_phys_data
!
      integer(kind = kint), intent(in) :: numnod
      real(kind = kreal), intent(in) :: xx(numnod,3)
!
      integer (kind = kint) :: inod
!
!$omp parallel do
      do inod = 1, numnod
       d_nod(inod,iphys%i_press)  = 0.0d0
       d_nod(inod,iphys%i_velo  ) = -xx(inod,2)
       d_nod(inod,iphys%i_velo+1) =  xx(inod,1)
       d_nod(inod,iphys%i_velo+2) = 0.0d0
      end do
!$omp end parallel do
!
      end subroutine set_initial_velo_3
!
!-----------------------------------------------------------------------
!
      end module set_initial_rotation
