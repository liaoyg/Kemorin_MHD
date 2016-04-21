!>@file   t_spheric_parameter.f90
!!@brief  module t_spheric_parameter
!!
!!@author H. Matsui
!!@date Programmed in July, 2007
!
!>@brief  Structure for indexing table of speherical harmonics transform
!!
!!@verbatim
!!      iflag_radial_grid:: radial grid type
!!        igrid_Chebyshev =    2 :: Chebyshev collocation points
!!        igrid_non_euqidist = 1 :: non-equi-distance
!!        igrid_euqidistance = 0 :: equi-distance
!!
!!      subroutine alloc_type_spheric_parameter(sph)
!!        type(sph_grids), intent(inout) :: sph
!!      subroutine dealloc_type_spheric_parameter(sph)
!!        type(sph_grids), intent(inout) :: sph
!!
!!      subroutine alloc_type_spheric_param_rtp(rtp)
!!        type(sph_rtp_grid), intent(inout) :: rtp
!!      subroutine alloc_type_spheric_param_rtm(rtm)
!!        type(sph_rtm_grid), intent(inout) :: rtm
!!      subroutine alloc_type_spheric_param_rlm(rlm)
!!        type(sph_rlm_grid), intent(inout) :: rlm
!!      subroutine alloc_type_spheric_param_rj(rj)
!!        type(sph_rj_grid), intent(inout) :: rj
!!
!!      subroutine alloc_type_sph_1d_index_rtp(rtp)
!!        type(sph_rtp_grid), intent(inout) :: rtp
!!      subroutine alloc_type_sph_1d_index_rtm(rtm)
!!        type(sph_rtm_grid), intent(inout) :: rtm
!!      subroutine alloc_type_sph_1d_index_rlm(rlm)
!!        type(sph_rlm_grid), intent(inout) :: rlm
!!      subroutine alloc_type_sph_1d_index_rj(rj)
!!        type(sph_rj_grid), intent(inout) :: rj
!!
!!      subroutine dealloc_type_spheric_param_rtp(rtp)
!!        type(sph_rtp_grid), intent(inout) :: rtp
!!      subroutine dealloc_type_spheric_param_rtm(rtm)
!!        type(sph_rtm_grid), intent(inout) :: rtm
!!      subroutine dealloc_type_spheric_param_rlm(rlm)
!!        type(sph_rlm_grid), intent(inout) :: rlm
!!      subroutine dealloc_spheric_param_rj(rj)
!!        type(sph_rj_grid), intent(inout) :: rj
!!
!!      subroutine dealloc_type_sph_1d_index_rtp(rtp)
!!        type(sph_rtp_grid), intent(inout) :: rtp
!!      subroutine dealloc_type_sph_1d_index_rtm(rtm)
!!        type(sph_rtm_grid), intent(inout) :: rtm
!!      subroutine dealloc_type_sph_1d_index_rlm(rlm)
!!        type(sph_rlm_grid), intent(inout) :: rlm
!!      subroutine dealloc_type_sph_1d_index_rj(rj)
!!        type(sph_rj_grid), intent(inout) :: rj
!!
!!      subroutine check_type_spheric_para_gl_part(sph)
!!        type(sph_grids), intent(in) :: sph
!!      subroutine check_type_spheric_parameter(my_rank, sph)
!!        integer(kind = kint), intent(in) :: my_rank
!!        type(sph_grids), intent(in) :: sph
!!      subroutine check_type_spheric_param_rtp(my_rank, rtp)
!!        integer(kind = kint), intent(in) :: my_rank
!!        type(sph_rtp_grid), intent(in) :: rtp
!!      subroutine check_type_spheric_param_rtm(my_rank, rtm)
!!        integer(kind = kint), intent(in) :: my_rank
!!        type(sph_rtm_grid), intent(in) :: rtm
!!      subroutine check_type_spheric_param_rlm(my_rank, rlm)
!!        integer(kind = kint), intent(in) :: my_rank
!!        type(sph_rlm_grid), intent(in) :: rlm
!!      subroutine check_type_spheric_param_rj(my_rank, rj)
!!        integer(kind = kint), intent(in) :: my_rank
!!        type(sph_rj_grid), intent(in) :: rj
!!
!!      integer(kind = kint) function find_local_sph_address(rj, l, m)
!!        type(sph_rj_grid), intent(in) :: rj
!!        integer(kind = 4), intent(in) :: l, m
!!      integer(kind = kint) function local_sph_node_address            &
!!     &                            (rj, kr, j_lc)
!!        type(sph_rj_grid), intent(in) :: rj
!!@endverbatim
!!
!!@n @param  my_rank     Running rank ID
!!
      module t_spheric_parameter
!
      use m_precision
      use m_spheric_constants
!
      implicit none
!
!    Group names for spherical shell dynamos
!
!>      Group name for ICB
      character(len=kchara), parameter :: ICB_nod_grp_name = 'ICB'
!>      Group name for CMB
      character(len=kchara), parameter :: CMB_nod_grp_name = 'CMB'
!>      Group name for innermost radius
      character(len=kchara), parameter                                  &
     &                      :: CTR_nod_grp_name = 'to_Center'
!
!>      Element Group name for inner core
      character(len=kchara), parameter                                  &
     &                      :: IC_ele_grp_name = 'inner_core'
!>      Element Group name for outer core
      character(len=kchara), parameter                                  &
     &                      :: OC_ele_grp_name = 'outer_core'
!>      Element Group name for outer core
      character(len=kchara), parameter                                  &
     &                      :: MT_ele_grp_name = 'external'
!
!>      Surface Group name for ICB
      character(len=kchara), parameter :: ICB_sf_grp_name = 'ICB_surf'
!>      Surface Group name for CMB
      character(len=kchara), parameter :: CMB_sf_grp_name = 'CMB_surf'
!>      Group name for innermost radius
      character(len=kchara), parameter                                  &
     &                      :: CTR_sf_grp_name = 'to_Center_surf'
!
!
!>        structure of index table for @f$ f(r,\theta,\phi) @f$
      type sph_rtp_grid
!>        number of global 1d data points for @f$ f(r,\theta,\phi) @f$
        integer(kind = kint) :: nidx_global_rtp(3)
!>        1d subdomain ID for @f$ f(r,\theta,\phi) @f$ (start from 0)
        integer(kind = kint) :: sph_rank_rtp(3)
!
!>        number of data points for @f$ f(r,\theta,\phi) @f$
        integer(kind = kint) :: nnod_rtp
!
!>        number of 1d data points for @f$ f(r,\theta,\phi) @f$
        integer(kind = kint) :: nidx_rtp(3)
!>        number of increments for @f$ f(r,\theta,\phi) @f$
        integer(kind = kint) :: istep_rtp(3)
!>        1d start address of global data for @f$ f(r,\theta,\phi) @f$
        integer(kind = kint) :: ist_rtp(3)
!>        1d end address of global data for @f$ f(r,\theta,\phi) @f$
        integer(kind = kint) :: ied_rtp(3)
!
!>        SMP stack for spectr data @f$ f(r,t,p) @f$
      integer(kind = kint), pointer :: istack_inod_rtp_smp(:)
!
!>        SMP stacks for indexing @f$ r@f$
      integer(kind = kint), pointer :: istack_rtp_kr_smp(:)
!>        SMP stacks for indexing @f$ t @f$
      integer(kind = kint), pointer :: istack_rtp_lt_smp(:)
!>        SMP stacks for indexing @f$ p @f$
      integer(kind = kint), pointer :: istack_rtp_mp_smp(:)
!
!>        SMP stacks for indexing @f$ r, t@f$
      integer(kind = kint), pointer :: istack_rtp_rt_smp(:)
!
!>        global address for each direction @f$ f(r,\theta,\phi) @f$
        integer(kind = kint), pointer :: idx_global_rtp(:,:)
!
!>        radial global address for @f$ f(r,\theta,\phi) @f$
        integer(kind = kint), pointer :: idx_gl_1d_rtp_r(:)
!>        meridional global address for @f$ f(r,\theta,\phi) @f$
        integer(kind = kint), pointer :: idx_gl_1d_rtp_t(:)
!>        zonal global address for @f$ f(r,\theta,\phi) @f$
        integer(kind = kint), pointer :: idx_gl_1d_rtp_p(:,:)
!
!>        1d radius data for @f$ f(r,\theta,\phi) @f$
        real(kind = kreal), pointer :: radius_1d_rtp_r(:)
!>        1 / radius_1d_rtp_r
        real(kind = kreal), pointer :: a_r_1d_rtp_r(:)
      end type sph_rtp_grid
!
!
!>        structure of index table for @f$ f(r,\theta,m) @f$
      type sph_rtm_grid
!>        number of global 1d data points for @f$ f(r,\theta,m) @f$
        integer(kind = kint) :: nidx_global_rtm(3)
!>        1d subdomain ID for @f$ f(r,\theta,m) @f$ (start from 0)
        integer(kind = kint) :: sph_rank_rtm(3)
!
!>        number of data points for @f$ f(r,\theta,m) @f$
        integer(kind = kint) :: nnod_rtm
!>        number of 1d data points for @f$ f(r,\theta,m) @f$
        integer(kind = kint) :: nidx_rtm(3)
!>        number of increments for @f$ f(r,\theta,m) @f$
        integer(kind = kint) :: istep_rtm(3)
!>        1d start address of global data for @f$ f(r,\theta,m) @f$
        integer(kind = kint) :: ist_rtm(3)
!>        1d end address of global data for @f$ f(r,\theta,m) @f$
        integer(kind = kint) :: ied_rtm(3)
!
!>        SMP stack for spectr data @f$ f(r,t,m) @f$
      integer(kind = kint), pointer :: istack_inod_rtm_smp(:)
!
!>        SMP stacks for indexing @f$ r@f$
      integer(kind = kint), pointer :: istack_rtm_kr_smp(:)
!>        SMP stacks for indexing @f$ t @f$
      integer(kind = kint), pointer :: istack_rtm_lt_smp(:)
!>        SMP stacks for indexing @f$ m @f$
      integer(kind = kint), pointer :: istack_rtm_m_smp(:)
!
!>        SMP stacks for indexing @f$ r, t@f$
      integer(kind = kint), pointer :: istack_rtm_rt_smp(:)
!
!>        Start address for @f$ m = 0 @f$ for @f$ f(r,\theta,m) @f$
        integer (kind=kint) :: ist_rtm_order_zero =   0
!
!>        global address for each direction @f$ f(r,\theta,m) @f$
        integer(kind = kint), pointer :: idx_global_rtm(:,:)
!
!>        radial global address for @f$ f(r,\theta,m) @f$
        integer(kind = kint), pointer :: idx_gl_1d_rtm_r(:)
!>        meridional global address for @f$ f(r,\theta,m) @f$
        integer(kind = kint), pointer :: idx_gl_1d_rtm_t(:)
!>        Zonal wave number for @f$ f(r,\theta,m) @f$
        integer(kind = kint), pointer :: idx_gl_1d_rtm_m(:,:)
!
!>        1d radius data for @f$ f(r,\theta,m) @f$
        real(kind = kreal), pointer :: radius_1d_rtm_r(:)
!>        1 / radius_1d_rtm_r
        real(kind = kreal), pointer :: a_r_1d_rtm_r(:)
      end type sph_rtm_grid
!
!
!>        structure of index table for @f$ f(r,l,m) @f$
      type sph_rlm_grid
!>        number of global 1d data points for @f$ f(r,l,m) @f$
        integer(kind = kint) :: nidx_global_rlm(2)
!>        1d subdomain ID for @f$ f(r,l,m) @f$ (start from 0)
        integer(kind = kint) :: sph_rank_rlm(2)
!
!>        number of data points for @f$ f(r,l,m) @f$
        integer(kind = kint) :: nnod_rlm
!>        number of 1d data points for @f$ f(r,l,m) @f$
        integer(kind = kint) :: nidx_rlm(2)
!>        number of increments for @f$ f(r,l,m) @f$
        integer(kind = kint) :: istep_rlm(2)
!>        1d start address of global data for @f$ f(r,l,m) @f$
        integer(kind = kint) :: ist_rlm(2)
!>        1d end address of global data for @f$ f(r,l,m) @f$
        integer(kind = kint) :: ied_rlm(2)
!
!>        SMP stack for spectr data @f$ f(r,l,m) @f$
      integer(kind = kint), pointer :: istack_inod_rlm_smp(:)
!
!>        SMP stacks for indexing @f$ r@f$
      integer(kind = kint), pointer :: istack_rlm_kr_smp(:)
!>        SMP stacks for indexing @f$ j @f$
      integer(kind = kint), pointer :: istack_rlm_j_smp(:)
!
!>        global address for each direction @f$ f(r,l,m) @f$
        integer(kind = kint), pointer :: idx_global_rlm(:,:)
!
!>        radial global address for @f$ f(r,l,m) @f$
        integer(kind = kint), pointer :: idx_gl_1d_rlm_r(:)
!>        spherical harmonics mode for  @f$ f(r,l,m) @f$
!!@n        idx_gl_1d_rj_j(j,1): global ID for spherical harmonics
!!@n        idx_gl_1d_rj_j(j,2): spherical hermonincs degree
!!@n        idx_gl_1d_rj_j(j,3): spherical hermonincs order
        integer(kind = kint), pointer :: idx_gl_1d_rlm_j(:,:)
!
!>        1d radius data for @f$ f(r,l,m) @f$
        real(kind = kreal), pointer :: radius_1d_rlm_r(:)
!>        1 / radius_1d_rlm_r
        real(kind = kreal), pointer :: a_r_1d_rlm_r(:)
      end type sph_rlm_grid
!
!
!>        structure of index table for @f$ f(r,j) @f$
      type sph_rj_grid
!>        number of global 1d data points for @f$ f(r,j) @f$
        integer(kind = kint) :: nidx_global_rj(2)
!>        1d subdomain ID for @f$ f(r,j) @f$ (start from 0)
        integer(kind = kint) :: irank_sph_rj(2)
!
!>        local spectr index for @f$ l = m = 0 @f$
!!@n      If @f$ l = m = 0 @f$ mode does not exist in subdomain, 
!!@n      idx_rj_degree_zero = 0.
        integer (kind=kint) :: idx_rj_degree_zero
!>        local spectr index for @f$ l = 1@f$ and  @f$ m = -1, 0, 1@f$.
!!        for @f$ f(r,j) @f$
!!@n        If spectr data do not exist in subdomain,
!!@n        idx_rj_degree_one(m) = 0.
        integer (kind=kint) :: idx_rj_degree_one(-1:1)
!
!>        number of data points for @f$ f(r,j) @f$
        integer(kind = kint) :: nnod_rj
!>        number of 1d data points for @f$ f(r,j) @f$
        integer(kind = kint) :: nidx_rj(2)
!>        number of increments for @f$ f(r,j) @f$
        integer(kind = kint) :: istep_rj(2)
!>        1d start address of global data for @f$ f(r,j) @f$
        integer(kind = kint) :: ist_rj(2)
!>        1d end address of global data for @f$ f(r,j) @f$
        integer(kind = kint) :: ied_rj(2)
!
!>        SMP stack for spectr data @f$ f(r,j) @f$
      integer(kind = kint), pointer :: istack_inod_rj_smp(:)
!
!>        SMP stacks for indexing @f$ (r,j) @f$
      integer(kind = kint), pointer :: istack_rj_kr_smp(:)
!>        SMP stacks for indexing @f$ (r,j) @f$
      integer(kind = kint), pointer :: istack_rj_j_smp(:)
!
!>        global address for each direction @f$ f(r,j) @f$
        integer(kind = kint), pointer :: idx_global_rj(:,:)
!
!>        radial global address @f$ f(r,j) @f$
        integer(kind = kint), pointer :: idx_gl_1d_rj_r(:)
!>        spherical harmonics mode for  @f$ f(r,j) @f$
!!@n        idx_gl_1d_rj_j(j,1): global ID for spherical harmonics
!!@n        idx_gl_1d_rj_j(j,2): spherical hermonincs degree
!!@n        idx_gl_1d_rj_j(j,3): spherical hermonincs order
        integer(kind = kint), pointer :: idx_gl_1d_rj_j(:,:)
!
!>        1d radius data for @f$ f(r,j) @f$
        real(kind = kreal), pointer :: radius_1d_rj_r(:)
!>        1d @f$1 / r @f$ for @f$ f(r,j) @f$
        real(kind = kreal), pointer :: a_r_1d_rj_r(:)
!
!>        1d @f$1 / r @f$ for @f$ f(r,j) @f$
!!@n@see  set_radius_func_cheby or set_radius_func_cheby
        real(kind = kreal), pointer :: ar_1d_rj(:,:)
!
!>        1d radius between grids for @f$ f(r,j) @f$
        real(kind = kreal), pointer :: r_ele_rj(:)
!>        1d @f$1 / r @f$ between grids for @f$ f(r,j) @f$
        real(kind = kreal), pointer :: ar_ele_rj(:,:)
      end type sph_rj_grid
!
!
!>  Structure of grid and spectr data for spherical spectr method
      type sph_grids
!>        integer flag for FEM mesh type
!!@n      iflag_MESH_same:     same grid point as Gauss-Legendre points
!!@n      iflag_MESH_w_pole:   Gauss-Legendre points with poles
!!@n      iflag_MESH_w_center: Gauss-Legendre points with center and poles
        integer (kind=kint) :: iflag_shell_mode =   iflag_MESH_same
!>        integer flag for center point in spectr data
!!@n      This flag should have same value for all processes
!!@n      0: No center point
!!@n      1: Center point is there
        integer (kind=kint) :: iflag_rj_center =  0
!>        radial grid type flag
!!@n      igrid_Chebyshev =    2 :: Chebyshev collocation points
!!@n      igrid_non_euqidist = 1 :: non-equi-distance
!!@n      igrid_euqidistance = 0 :: equi-distance
        integer (kind=kint) :: iflag_radial_grid =  igrid_non_euqidist
!
!>        local spectr index for @f$ l = m = 0 @f$ at center
!!@n      if center does not exist in subdomain, inod_rj_center = 0.
        integer (kind=kint) :: inod_rj_center =   0
!
!>        Truncation for spherical harmonics
        integer(kind = kint) :: l_truncation
!>        m-folding symmetry for longitudinal direction
        integer(kind = kint) :: m_folding = 1
!
!>        structure of index table for @f$ f(r,\theta,\phi) @f$
        type(sph_rtp_grid) :: sph_rtp
!>        structure of index table for @f$ f(r,\theta,m) @f$
        type(sph_rtm_grid) :: sph_rtm
!>        structure of index table for @f$ f(r,l,m) @f$
        type(sph_rlm_grid) :: sph_rlm
!>        structure of index table for @f$ f(r,j) @f$
        type(sph_rj_grid) ::  sph_rj
      end type sph_grids
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine alloc_type_spheric_parameter(sph)
!
      type(sph_grids), intent(inout) :: sph
!
!
      call alloc_type_spheric_param_rtp(sph%sph_rtp)
      call alloc_type_spheric_param_rtm(sph%sph_rtm)
      call alloc_type_spheric_param_rlm(sph%sph_rlm)
      call alloc_type_spheric_param_rj(sph%sph_rj)
!
      end subroutine alloc_type_spheric_parameter
!
! ----------------------------------------------------------------------
!
      subroutine dealloc_type_spheric_parameter(sph)
!
      type(sph_grids), intent(inout) :: sph
!
!
      call dealloc_type_spheric_param_rtp(sph%sph_rtp)
      call dealloc_type_spheric_param_rtm(sph%sph_rtm)
      call dealloc_type_spheric_param_rlm(sph%sph_rlm)
      call dealloc_spheric_param_rj(sph%sph_rj)
!
      end subroutine dealloc_type_spheric_parameter
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine alloc_type_spheric_param_rtp(rtp)
!
      type(sph_rtp_grid), intent(inout) :: rtp
!
      allocate(rtp%idx_global_rtp(rtp%nnod_rtp,3))
      if(rtp%nnod_rtp .gt. 0) rtp%idx_global_rtp = 0
!
      end subroutine alloc_type_spheric_param_rtp
!
! ----------------------------------------------------------------------
!
      subroutine alloc_type_spheric_param_rtm(rtm)
!
      type(sph_rtm_grid), intent(inout) :: rtm
!
      allocate(rtm%idx_global_rtm(rtm%nnod_rtm,3))
!
      if(rtm%nnod_rtm .gt. 0) rtm%idx_global_rtm = 0
!
      end subroutine alloc_type_spheric_param_rtm
!
! ----------------------------------------------------------------------
!
      subroutine alloc_type_spheric_param_rlm(rlm)
!
      type(sph_rlm_grid), intent(inout) :: rlm
!
      allocate(rlm%idx_global_rlm(rlm%nnod_rlm,2))
!
      if(rlm%nnod_rlm .gt. 0) then
        rlm%idx_global_rlm = 0
      end if
!
      end subroutine alloc_type_spheric_param_rlm
!
! ----------------------------------------------------------------------
!
      subroutine alloc_type_spheric_param_rj(rj)
!
      type(sph_rj_grid), intent(inout) :: rj
!
      allocate(rj%idx_global_rj(rj%nnod_rj,2))
!
      if(rj%nnod_rj .gt. 0) rj%idx_global_rj =  0
!
      end subroutine alloc_type_spheric_param_rj
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine alloc_type_sph_1d_index_rtp(rtp)
!
      type(sph_rtp_grid), intent(inout) :: rtp
      integer(kind = kint) :: num
!
      num = rtp%nidx_rtp(1)
      allocate(rtp%idx_gl_1d_rtp_r(num))
      allocate(rtp%radius_1d_rtp_r(num))
      allocate(rtp%a_r_1d_rtp_r(num))
      num = rtp%nidx_rtp(2)
      allocate(rtp%idx_gl_1d_rtp_t(num))
      num = rtp%nidx_rtp(3)
      allocate(rtp%idx_gl_1d_rtp_p(num,2))
!
      if(rtp%nidx_rtp(3) .gt. 0) rtp%idx_gl_1d_rtp_p = 0
      if(rtp%nidx_rtp(2) .gt. 0) rtp%idx_gl_1d_rtp_t = 0
      if(rtp%nidx_rtp(1) .gt. 0) then
        rtp%idx_gl_1d_rtp_r = 0
        rtp%radius_1d_rtp_r = 0.0d0
        rtp%a_r_1d_rtp_r = 0.0d0
      end if
!
      end subroutine alloc_type_sph_1d_index_rtp
!
! ----------------------------------------------------------------------
!
      subroutine alloc_type_sph_1d_index_rtm(rtm)
!
      type(sph_rtm_grid), intent(inout) :: rtm
      integer(kind = kint) :: num
!
      num = rtm%nidx_rtm(1)
      allocate(rtm%idx_gl_1d_rtm_r(num))
      allocate(rtm%radius_1d_rtm_r(num))
      allocate(rtm%a_r_1d_rtm_r(num))
      num = rtm%nidx_rtm(2)
      allocate(rtm%idx_gl_1d_rtm_t(num))
      num = rtm%nidx_rtm(3)
      allocate(rtm%idx_gl_1d_rtm_m(num,2))
!
      if(rtm%nidx_rtm(3) .gt. 0) rtm%idx_gl_1d_rtm_m = 0
      if(rtm%nidx_rtm(2) .gt. 0) rtm%idx_gl_1d_rtm_t = 0
      if(rtm%nidx_rtm(1) .gt. 0) then
        rtm%idx_gl_1d_rtm_r = 0
        rtm%radius_1d_rtm_r = 0.0d0
        rtm%a_r_1d_rtm_r = 0.0d0
      end if
!
      end subroutine alloc_type_sph_1d_index_rtm
!
! ----------------------------------------------------------------------
!
      subroutine alloc_type_sph_1d_index_rlm(rlm)
!
      type(sph_rlm_grid), intent(inout) :: rlm
      integer(kind = kint) :: num
!
      num = rlm%nidx_rlm(1)
      allocate(rlm%idx_gl_1d_rlm_r(num))
      allocate(rlm%radius_1d_rlm_r(num))
      allocate(rlm%a_r_1d_rlm_r(num))
      num = rlm%nidx_rlm(2)
      allocate(rlm%idx_gl_1d_rlm_j(num,3))
!
      if(rlm%nidx_rlm(2) .gt. 0) rlm%idx_gl_1d_rlm_j = 0
      if(rlm%nidx_rlm(1) .gt. 0) then
        rlm%idx_gl_1d_rlm_r = 0
        rlm%radius_1d_rlm_r = 0.0d0
        rlm%a_r_1d_rlm_r    = 0.0d0
      end if
!
      end subroutine alloc_type_sph_1d_index_rlm
!
! ----------------------------------------------------------------------
!
      subroutine alloc_type_sph_1d_index_rj(rj)
!
      type(sph_rj_grid), intent(inout) :: rj
      integer(kind = kint) :: num
!
      num = rj%nidx_rj(1)
      allocate(rj%idx_gl_1d_rj_r(num))
      allocate(rj%radius_1d_rj_r(num))
      allocate(rj%a_r_1d_rj_r(num))
!
      allocate(rj%ar_1d_rj(num,3))
      allocate(rj%r_ele_rj(num))
      allocate(rj%ar_ele_rj(num,3))
!
      num = rj%nidx_rj(2)
      allocate(rj%idx_gl_1d_rj_j(num,3))
!
      rj%idx_rj_degree_zero = 0
      rj%idx_rj_degree_one(-1:1) = 0
!
      if(rj%nidx_rj(2) .gt. 0) rj%idx_gl_1d_rj_j = 0
      if(rj%nidx_rj(1) .gt. 0) then
        rj%idx_gl_1d_rj_r = 0
        rj%radius_1d_rj_r = 0.0d0
        rj%a_r_1d_rj_r = 0.0d0
!
        rj%ar_1d_rj = 0.0d0
        rj%r_ele_rj = 0.0d0
        rj%ar_ele_rj = 0.0d0
      end if
!
      end subroutine alloc_type_sph_1d_index_rj
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine dealloc_type_spheric_param_rtp(rtp)
!
      type(sph_rtp_grid), intent(inout) :: rtp
!
      deallocate(rtp%idx_global_rtp)
!
      end subroutine dealloc_type_spheric_param_rtp
!
! ----------------------------------------------------------------------
!
      subroutine dealloc_type_spheric_param_rtm(rtm)
!
      type(sph_rtm_grid), intent(inout) :: rtm
!
      deallocate(rtm%idx_global_rtm)
!
      end subroutine dealloc_type_spheric_param_rtm
!
! ----------------------------------------------------------------------
!
      subroutine dealloc_type_spheric_param_rlm(rlm)
!
      type(sph_rlm_grid), intent(inout) :: rlm
!
      deallocate(rlm%idx_global_rlm)
!
      end subroutine dealloc_type_spheric_param_rlm
!
! ----------------------------------------------------------------------
!
      subroutine dealloc_spheric_param_rj(rj)
!
      type(sph_rj_grid), intent(inout) :: rj
!
      deallocate(rj%idx_global_rj)
!
      end subroutine dealloc_spheric_param_rj
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine dealloc_type_sph_1d_index_rtp(rtp)
!
      type(sph_rtp_grid), intent(inout) :: rtp
!
      deallocate(rtp%radius_1d_rtp_r)
      deallocate(rtp%a_r_1d_rtp_r)
      deallocate(rtp%idx_gl_1d_rtp_r)
      deallocate(rtp%idx_gl_1d_rtp_t)
      deallocate(rtp%idx_gl_1d_rtp_p)
!
      end subroutine dealloc_type_sph_1d_index_rtp
!
! ----------------------------------------------------------------------
!
      subroutine dealloc_type_sph_1d_index_rtm(rtm)
!
      type(sph_rtm_grid), intent(inout) :: rtm
!
!
      deallocate(rtm%radius_1d_rtm_r)
      deallocate(rtm%a_r_1d_rtm_r)
      deallocate(rtm%idx_gl_1d_rtm_r)
      deallocate(rtm%idx_gl_1d_rtm_t)
      deallocate(rtm%idx_gl_1d_rtm_m)
!
      end subroutine dealloc_type_sph_1d_index_rtm
!
! ----------------------------------------------------------------------
!
      subroutine dealloc_type_sph_1d_index_rlm(rlm)
!
      type(sph_rlm_grid), intent(inout) :: rlm
!
!
      deallocate(rlm%radius_1d_rlm_r)
      deallocate(rlm%a_r_1d_rlm_r   )
      deallocate(rlm%idx_gl_1d_rlm_r)
      deallocate(rlm%idx_gl_1d_rlm_j)
!
      end subroutine dealloc_type_sph_1d_index_rlm
!
! ----------------------------------------------------------------------
!
      subroutine dealloc_type_sph_1d_index_rj(rj)
!
      type(sph_rj_grid), intent(inout) :: rj
!
!
      deallocate(rj%radius_1d_rj_r, rj%a_r_1d_rj_r)
      deallocate(rj%ar_1d_rj, rj%r_ele_rj, rj%ar_ele_rj)
      deallocate(rj%idx_gl_1d_rj_r, rj%idx_gl_1d_rj_j)
!
      end subroutine dealloc_type_sph_1d_index_rj
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine check_type_spheric_para_gl_part(sph)
!
      type(sph_grids), intent(in) :: sph
!
      write(50,*) 'l_truncation ', sph%l_truncation
      write(50,*) 'm_folding ',    sph%m_folding
      write(50,*) 'nidx_global_rtm ', sph%sph_rtm%nidx_global_rtm(1:3)
      write(50,*) 'nidx_global_rlm ', sph%sph_rlm%nidx_global_rlm(1:2)
      write(50,*) 'nidx_global_rj ',  sph%sph_rj%nidx_global_rj(1:2)
!
!
      end subroutine check_type_spheric_para_gl_part
!
! ----------------------------------------------------------------------
!
      subroutine check_type_spheric_parameter(my_rank, sph)
!
      integer(kind = kint), intent(in) :: my_rank
      type(sph_grids), intent(in) :: sph
!
      call check_type_spheric_param_rtp(my_rank, sph%sph_rtp)
      call check_type_spheric_param_rtm(my_rank, sph%sph_rtm)
      call check_type_spheric_param_rlm(my_rank, sph%sph_rlm)
      call check_type_spheric_param_rj(my_rank, sph%sph_rj)
!
      end subroutine check_type_spheric_parameter
!
! -----------------------------------------------------------------------
!
      subroutine check_type_spheric_param_rtp(my_rank, rtp)
!
      integer(kind = kint), intent(in) :: my_rank
      type(sph_rtp_grid), intent(in) :: rtp
      integer(kind = kint) :: i
!
!
      write(my_rank+50,*) 'sph_rank_rtp ', rtp%sph_rank_rtp(1:3)
      write(my_rank+50,*) 'nidx_rtp ', rtp%nidx_rtp(1:3)
      write(my_rank+50,*) 'nnod_rtp ', rtp%nnod_rtp
!
      write(my_rank+50,*)  'i, idx_global_rtp(r,t,p)'
      do i = 1, rtp%nnod_rtp
        write(my_rank+50,*) i, rtp%idx_global_rtp(i,1:3)
      end do
!
      end subroutine check_type_spheric_param_rtp
!
! -----------------------------------------------------------------------
!
      subroutine check_type_spheric_param_rtm(my_rank, rtm)
!
      integer(kind = kint), intent(in) :: my_rank
      type(sph_rtm_grid), intent(in) :: rtm
      integer(kind = kint) :: i
!
!
      write(my_rank+50,*) 'sph_rank_rtm ', rtm%sph_rank_rtm(1:3)
      write(my_rank+50,*) 'nidx_rtm ', rtm%nidx_rtm(1:3)
      write(my_rank+50,*) 'nnod_rtm ', rtm%nnod_rtm
!
      write(my_rank+50,*) 'i, idx_global_rtm(r,t,p)'
      do i = 1, rtm%nnod_rtm
        write(my_rank+50,*) i, rtm%idx_global_rtm(i,1:3)
      end do
!
      end subroutine check_type_spheric_param_rtm
!
! -----------------------------------------------------------------------
!
      subroutine check_type_spheric_param_rlm(my_rank, rlm)
!
      integer(kind = kint), intent(in) :: my_rank
      type(sph_rlm_grid), intent(in) :: rlm
      integer(kind = kint) :: i
!
!
      write(my_rank+50,*) 'sph_rank_rlm ', rlm%sph_rank_rlm(1:2)
      write(my_rank+50,*) 'nidx_rlm ', rlm%nidx_rlm(1:2)
      write(my_rank+50,*) 'nnod_rlm ', rlm%nnod_rlm
!
      write(my_rank+50,*) 'i, idx_global_rlm(r,j)'
      do i = 1, rlm%nnod_rlm
        write(my_rank+50,*) i, rlm%idx_global_rlm(i,1:2)
      end do
!
      end subroutine check_type_spheric_param_rlm
!
! -----------------------------------------------------------------------
!
      subroutine check_type_spheric_param_rj(my_rank, rj)
!
      integer(kind = kint), intent(in) :: my_rank
      type(sph_rj_grid), intent(in) :: rj
      integer(kind = kint) :: i
!
!
      write(my_rank+50,*) 'irank_sph_rj ',  rj%irank_sph_rj(1:2)
      write(my_rank+50,*) 'nidx_rj  ',  rj%nidx_rj(1:2)
      write(my_rank+50,*) 'nnod_rj ',  rj%nnod_rj
!
      write(my_rank+50,*) 'i, idx_global_rj(r,j)'
      do i = 1, rj%nnod_rj
        write(my_rank+50,*) i, rj%idx_global_rj(i,1:2)
      end do
!
      end subroutine check_type_spheric_param_rj
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      integer(kind = kint) function find_local_sph_address(rj, l, m)
!
      type(sph_rj_grid), intent(in) :: rj
      integer(kind = 4), intent(in) :: l, m
!
      integer(kind = kint) :: j
!
!
      find_local_sph_address = 0
      do j = 1, rj%nidx_rj(2)
        if (   rj%idx_gl_1d_rj_j(j,2) .eq. l                            &
     &   .and. rj%idx_gl_1d_rj_j(j,3) .eq. m) then
          find_local_sph_address = j
          return
        end if
      end do
!
      end function find_local_sph_address
!
!-----------------------------------------------------------------------
!
      integer(kind = kint) function local_sph_node_address              &
     &                            (rj, kr, j_lc)
!
      type(sph_rj_grid), intent(in) :: rj
      integer(kind = kint), intent(in) :: kr, j_lc
!
!
      local_sph_node_address = j_lc + (kr-1)*rj%nidx_rj(2)
!
      end function local_sph_node_address
!
!-----------------------------------------------------------------------
!
      end module t_spheric_parameter
