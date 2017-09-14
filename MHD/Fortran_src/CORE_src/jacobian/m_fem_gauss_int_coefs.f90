!m_fem_gauss_int_coefs
!   module   m_fem_gauss_int_coefs
!.......................................................................
!
!
!      Written by H. Matsui on Dec. 2003
!      Modified by H. Matsui on Oct., 2006
!
!      subroutine max_int_point_by_etype(nnod_4_ele)
!      subroutine maximum_integration_points(num_int)
!      subroutine allocate_gauss_coef_4_fem
!
!      subroutine deallocate_gauss_coef_4_fem
!      subroutine set_num_of_int_points
!
      module   m_fem_gauss_int_coefs
!
      use m_precision
      use m_constants
      use t_fem_gauss_int_coefs
!
      implicit  none
!
      type(FEM_gauss_int_coefs), save :: g_FEM1
!
      integer (kind=kint), save :: max_int_point = 4
      integer (kind=kint), save :: maxtot_int_3d= 100
      integer (kind=kint), save :: maxtot_int_2d= 30
      integer (kind=kint), save :: maxtot_int_1d= 10
!
!
      integer (kind=kint), allocatable :: int_start1(:)
      integer (kind=kint), allocatable :: int_start2(:)
      integer (kind=kint), allocatable :: int_start3(:)
!
      real (kind=kreal), allocatable :: owe(:)
      real (kind=kreal), allocatable :: owe2d(:)
      real (kind=kreal), allocatable :: owe3d(:)
! 
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!> Set maximum number for integration points of FEM
!
      subroutine max_int_point_by_etype(nnod_4_ele)
!
      use m_geometry_constants
!
      integer(kind = kint), intent(in) :: nnod_4_ele
!
!
      if(nnod_4_ele.eq.num_t_quad  .or. nnod_4_ele.eq.num_t_lag) then
        call maximum_integration_points(ithree)
      else if(nnod_4_ele .eq. ione) then
        call maximum_integration_points(ione)
      else
        call maximum_integration_points(itwo)
      end if
!
      end subroutine max_int_point_by_etype
!
! ----------------------------------------------------------------------
!
      subroutine maximum_integration_points(num_int)
!
      integer(kind = kint), intent(in) :: num_int
!
      max_int_point = num_int
!
      end subroutine maximum_integration_points
!
!-----------------------------------------------------------------------
!
      subroutine allocate_gauss_coef_4_fem
!
      integer(kind = kint) :: n
!
      allocate( owe(maxtot_int_1d)   )
      allocate( owe2d(maxtot_int_2d) )
      allocate( owe3d(maxtot_int_3d) )
!
      allocate( int_start1(max_int_point) )
      allocate( int_start2(max_int_point) )
      allocate( int_start3(max_int_point) )
!
      owe =   0.0d0
      owe2d = 0.0d0
      owe3d = 0.0d0
!
      int_start3(1) = 0
      int_start2(1) = 0
      int_start1(1) = 0
      do n = 2, max_int_point
        int_start3(n) = int_start3(n-1) + (n-1)*(n-1)*(n-1)
        int_start2(n) = int_start2(n-1) + (n-1)*(n-1)
        int_start1(n) = int_start1(n-1) + (n-1)
      end do

!
      end subroutine allocate_gauss_coef_4_fem
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine deallocate_gauss_coef_4_fem
!
      deallocate( owe, owe2d, owe3d )
      deallocate( int_start1, int_start2, int_start3 )
!
      end subroutine deallocate_gauss_coef_4_fem
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine set_num_of_int_points
!
      integer(kind = kint) :: n
!
!
      maxtot_int_3d = 0
      maxtot_int_2d = 0
      maxtot_int_1d = 0
      do n = 1, max_int_point
        maxtot_int_3d = maxtot_int_3d + n*n*n
        maxtot_int_2d = maxtot_int_2d + n*n
        maxtot_int_1d = maxtot_int_1d + n
      end do
!
      end subroutine set_num_of_int_points
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine copy_fem_gauss_int_coefs(g_FEM)
!
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
!
!
      max_int_point = g_FEM%max_int_point
!
      maxtot_int_3d = g_FEM%maxtot_int_3d
      maxtot_int_2d = g_FEM%maxtot_int_2d
      maxtot_int_1d = g_FEM%maxtot_int_1d
!
      allocate( owe(maxtot_int_1d)   )
      allocate( owe2d(maxtot_int_2d) )
      allocate( owe3d(maxtot_int_3d) )
!
      allocate( int_start1(max_int_point) )
      allocate( int_start2(max_int_point) )
      allocate( int_start3(max_int_point) )
!
      int_start1(1:max_int_point) = g_FEM%int_start1(1:max_int_point)
      int_start2(1:max_int_point) = g_FEM%int_start2(1:max_int_point)
      int_start3(1:max_int_point) = g_FEM%int_start3(1:max_int_point)
!
      owe(1:maxtot_int_1d) =   g_FEM%owe(1:maxtot_int_1d)
      owe2d(1:maxtot_int_2d) = g_FEM%owe2d(1:maxtot_int_2d)
      owe3d(1:maxtot_int_3d) = g_FEM%owe3d(1:maxtot_int_3d)
!
      end subroutine copy_fem_gauss_int_coefs
!
!-----------------------------------------------------------------------
!
!
      subroutine copy_fem_gauss_int_coef_type(g_FEM)
!
      type(FEM_gauss_int_coefs), intent(inout) :: g_FEM
!
!
      g_FEM%max_int_point = max_int_point
!
      g_FEM%maxtot_int_3d = maxtot_int_3d
      g_FEM%maxtot_int_2d = maxtot_int_2d
      g_FEM%maxtot_int_1d = maxtot_int_1d
!
      call alloc_gauss_coef_4_fem(g_FEM)
!
      g_FEM%owe(1:maxtot_int_1d) =   owe(1:maxtot_int_1d)
      g_FEM%owe2d(1:maxtot_int_2d) = owe2d(1:maxtot_int_2d)
      g_FEM%owe3d(1:maxtot_int_3d) = owe3d(1:maxtot_int_3d)
!
      end subroutine copy_fem_gauss_int_coef_type
!
!-----------------------------------------------------------------------
!
      end module   m_fem_gauss_int_coefs
