!>@file   t_legendre_work_matmul.f90
!!@brief  module t_legendre_work_matmul
!!
!!@author H. Matsui
!!@date Programmed in Aug., 2013
!
!>@brief  Work array for forward Legendre transform useing mat multi 
!>@n      data are strored communication buffer
!!
!!@verbatim
!!      subroutine alloc_leg_trns_matmul(nth_rtm, maxidx_rtm_r_smp,     &
!!     &          ntensor, nvector, nscalar, idx_trns, WK_l_mtl)
!!      subroutine dealloc_leg_vec_matmul(WK_l_mtl)
!!
!!     field data for Legendre transform
!!       original layout: vr_rtm(l_rtm,m_rtm,k_rtm,icomp)
!!       size: vr_rtm(nth_rtm,nidx_rtm(1)*ncomp,nidx_rtm(3))
!!      real(kind = kreal), allocatable :: vr_rtm(:,:,:)
!!
!!     spectr data for Legendre transform
!!       original layout: sp_rlm(j_rlm,k_rtm,icomp)
!!        size: sp_rlm(nidx_rlm(2),nidx_rtm(1)*ncomp)
!!      real(kind = kreal), allocatable :: sp_rlm(:,:)
!!@endverbatim
!!
!!@param   ncomp    Total number of components for spherical transform
!!@param   nvector  Number of vector for spherical transform
!!@param   nscalar  Number of scalar (including tensor components)
!!                  for spherical transform
!
      module t_legendre_work_matmul
!
      use m_precision
      use m_constants
      use m_machine_parameter
      use calypso_mpi
!
      use t_work_4_sph_trans
!
      use matmul_for_legendre_trans
!
      implicit none
!
!
!>      Work structure for Legendre trasform by matmul
      type leg_trns_matmul_work
!>        Maximum matrix size for spectr tensor
        integer(kind = kint) :: ntsr_jk
!!@n       @f$e_w@f$ component with 
!!          @f$ (-l(l+1) + (2 m^2 / \sin \theta)) P_{l}^{m} 
!!             - 2 (\cos \theta / \sin^2 \theta) dP_{l}^{m}/d\theta @f$
!!          and  evem (l-m)
        real(kind = kreal), allocatable :: tsrw_e(:,:)
!!@n       @f$e_\chi@f$ component with 
!!          @f$ (-l(l+1) + (2 m^2 / \sin \theta)) P_{l}^{m} 
!!             - 2 (\cos \theta / \sin^2 \theta) dP_{l}^{m}/d\theta @f$
!!          and  evem (l-m)
        real(kind = kreal), allocatable :: tsrchi_e(:,:)
!!@n       @f$e_\tau@f$ component with evem (l-m)
        real(kind = kreal), allocatable :: tsrtau_e(:,:)
!!@n       @f$e_rr@f$ component with evem (l-m)  
        real(kind = kreal), allocatable :: tsrr_e(:,:)
!!@n       Tensor toroidal component with phi derivative and evem (l-m)
        real(kind = kreal), allocatable :: dtsrtdp_e(:,:)
!!@n       Tensor poloidal component with phi derivative and evem (l-m)
        real(kind = kreal), allocatable :: dtsrhdp_e(:,:)
!!@n
!!@n       @f$e_w@f$ component with 
!>          @f$ 2m (-(\cos \theta / \sin^2 \theta)) P_{l}^{m}
!!             + (1 / \sin \theta) dP_{l}^{m}/d\theta @f$
!!          and  evem (l-m)
        real(kind = kreal), allocatable :: dtsrwdt_e(:,:)
!!@n       @f$e_\chi@f$ component with 
!>          @f$ 2m (-(\cos \theta / \sin^2 \theta)) P_{l}^{m}
!!             + (1 / \sin \theta) dP_{l}^{m}/d\theta @f$
!!          and  evem (l-m)
        real(kind = kreal), allocatable :: dtsrchidt_e(:,:)
!!@n       Tensor toroidal component with theta derivative and evem (l-m)
        real(kind = kreal), allocatable :: dtsrtdt_e(:,:)
!!@n       Tensor poloidal component with theta derivative and evem (l-m)
        real(kind = kreal), allocatable :: dtsrhdt_e(:,:)
!!
!!
!>       Maximum matrix size for symmeytric tensor field data
        integer(kind = kint) :: ntsr_lk
!>         Symmetric tensor P  component
        real(kind = kreal), allocatable :: symp_dp(:,:)
!>         Symmetric tensor theta-phi component
        real(kind = kreal), allocatable :: symp_tp(:,:)
!>         Symmetric tensor tau component
        real(kind = kreal), allocatable :: symp_tau(:,:)
!>         Symmetric tensor r-r component
        real(kind = kreal), allocatable :: symp_rr(:,:)
!>         Symmetric tensor r-theta-component with condugate order
        real(kind = kreal), allocatable :: symn_rt(:,:)
!>         Symmetric tensor r-phi-component with condugate order
        real(kind = kreal), allocatable :: symn_rp(:,:)
!
!>         Anti-symmetric tensor theta-phi component with condugate order
        real(kind = kreal), allocatable :: asmn_tp(:,:)
!>         Anti-symmetric tensor P  component with condugate order
        real(kind = kreal), allocatable :: asmn_dp(:,:)
!>         Anti-symmetric tensor r-phi-component
        real(kind = kreal), allocatable :: asmp_rp(:,:)
!>         Anti-symmetric tensor r-theta-component
        real(kind = kreal), allocatable :: asmp_rt(:,:)
!!
!>        Maximum matrix size for spectr data
        integer(kind = kint) :: nvec_jk
!>        Poloidal component with evem (l-m)
        real(kind = kreal), allocatable :: pol_e(:,:)
!>        radial difference of Poloidal component with evem (l-m)
        real(kind = kreal), allocatable :: dpoldt_e(:,:)
!>        radial difference of Poloidal component with evem (l-m)
        real(kind = kreal), allocatable :: dpoldp_e(:,:)
!>        Toroidal component with evem (l-m)
        real(kind = kreal), allocatable :: dtordt_e(:,:)
!>        Toroidal component with evem (l-m)
        real(kind = kreal), allocatable :: dtordp_e(:,:)
!
!>        Maximum matrix size for field data
        integer(kind = kint) :: nvec_lk
!>        Symmetric radial component
        real(kind = kreal), allocatable :: symp_r(:,:)
!>        Anti-symmetric theta-component
        real(kind = kreal), allocatable :: asmp_t(:,:)
!>        Anti-symmetric phi-component
        real(kind = kreal), allocatable :: asmp_p(:,:)
!>        Symmetric theta-component with condugate order
        real(kind = kreal), allocatable :: symn_t(:,:)
!>        Symmetric phi-component with condugate order
        real(kind = kreal), allocatable :: symn_p(:,:)
!
!>        Maximum matrix size for spectr data
        integer(kind = kint) :: nscl_jk
!>        Scalar with evem (l-m)
        real(kind = kreal), allocatable :: scl_e(:,:)
!
!>        Maximum matrix size for field data
        integer(kind = kint) :: nscl_lk
!>        Symmetric scalar component
        real(kind = kreal), allocatable :: symp(:,:)
      end type leg_trns_matmul_work
!
      private :: alloc_leg_tsr_matmul, alloc_leg_vec_matmul
      private :: alloc_leg_scl_matmul
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine alloc_leg_trns_matmul(nth_rtm, maxidx_rtm_r_smp,       &
     &          ntensor, nvector, nscalar, idx_trns, WK_l_mtl)
!
      integer(kind = kint), intent(in) :: nth_rtm
      integer(kind = kint), intent(in) :: maxidx_rtm_r_smp
      integer(kind = kint), intent(in) :: ntensor, nvector, nscalar
      type(index_4_sph_trans), intent(in) :: idx_trns
!
      type(leg_trns_matmul_work), intent(inout) :: WK_l_mtl
!
!
      WK_l_mtl%ntsr_jk = idx_trns%maxdegree_rlm                         &
     &                  * maxidx_rtm_r_smp * ntensor
      WK_l_mtl%nvec_jk = idx_trns%maxdegree_rlm                         &
     &                  * maxidx_rtm_r_smp * nvector
      WK_l_mtl%nscl_jk = idx_trns%maxdegree_rlm                         &
     &                  * maxidx_rtm_r_smp * nscalar
      WK_l_mtl%ntsr_lk = nth_rtm * maxidx_rtm_r_smp * ntensor
      WK_l_mtl%nvec_lk = nth_rtm * maxidx_rtm_r_smp * nvector
      WK_l_mtl%nscl_lk = nth_rtm * maxidx_rtm_r_smp * nscalar
!
      call alloc_leg_vec_matmul(WK_l_mtl)
      call alloc_leg_scl_matmul(WK_l_mtl)
!
      if(ntensor .le. 0) return
      call alloc_leg_tsr_matmul(WK_l_mtl)
!
      end subroutine alloc_leg_trns_matmul
!
! -----------------------------------------------------------------------
!
      subroutine dealloc_leg_vec_matmul(WK_l_mtl)
!
      type(leg_trns_matmul_work), intent(inout) :: WK_l_mtl
!
!
      deallocate(WK_l_mtl%pol_e, WK_l_mtl%dpoldt_e, WK_l_mtl%dpoldp_e)
      deallocate(WK_l_mtl%dtordt_e, WK_l_mtl%dtordp_e)
      deallocate(WK_l_mtl%symp_r, WK_l_mtl%symn_t, WK_l_mtl%symn_p)
      deallocate(WK_l_mtl%asmp_t, WK_l_mtl%asmp_p)
!
      deallocate(WK_l_mtl%scl_e, WK_l_mtl%symp)
!
      if(allocated(WK_l_mtl%tsrw_e)) then
        deallocate(WK_l_mtl%tsrw_e, WK_l_mtl%tsrchi_e)
        deallocate(WK_l_mtl%tsrtau_e, WK_l_mtl%tsrr_e)
        deallocate(WK_l_mtl%dtsrtdp_e, WK_l_mtl%dtsrhdp_e)
        deallocate(WK_l_mtl%dtsrwdt_e, WK_l_mtl%dtsrchidt_e)
        deallocate(WK_l_mtl%dtsrtdt_e, WK_l_mtl%dtsrhdt_e)
!
        deallocate(WK_l_mtl%symp_dp, WK_l_mtl%symp_tp)
        deallocate(WK_l_mtl%symp_tau, WK_l_mtl%symp_rr)
        deallocate(WK_l_mtl%symn_rt, WK_l_mtl%symn_rp)
        deallocate(WK_l_mtl%asmn_tp, WK_l_mtl%asmn_dp)
        deallocate(WK_l_mtl%asmp_rp, WK_l_mtl%asmp_rt)
      end if
!
      end subroutine dealloc_leg_vec_matmul
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine alloc_leg_tsr_matmul(WK_l_mtl)
!
      type(leg_trns_matmul_work), intent(inout) :: WK_l_mtl
!
!
      allocate(WK_l_mtl%tsrw_e(WK_l_mtl%ntsr_jk,np_smp))
      allocate(WK_l_mtl%tsrchi_e(WK_l_mtl%ntsr_jk,np_smp))
      allocate(WK_l_mtl%tsrtau_e(WK_l_mtl%ntsr_jk,np_smp))
      allocate(WK_l_mtl%tsrr_e(WK_l_mtl%ntsr_jk,np_smp))
      allocate(WK_l_mtl%dtsrtdp_e(WK_l_mtl%ntsr_jk,np_smp))
      allocate(WK_l_mtl%dtsrhdp_e(WK_l_mtl%ntsr_jk,np_smp))
!
      allocate(WK_l_mtl%dtsrwdt_e(WK_l_mtl%ntsr_jk,np_smp))
      allocate(WK_l_mtl%dtsrchidt_e(WK_l_mtl%ntsr_jk,np_smp))
      allocate(WK_l_mtl%dtsrtdt_e(WK_l_mtl%ntsr_jk,np_smp))
      allocate(WK_l_mtl%dtsrhdt_e(WK_l_mtl%ntsr_jk,np_smp))
!
      allocate(WK_l_mtl%symp_dp(WK_l_mtl%ntsr_lk,np_smp))
      allocate(WK_l_mtl%symp_tp(WK_l_mtl%ntsr_lk,np_smp))
      allocate(WK_l_mtl%symp_tau(WK_l_mtl%ntsr_lk,np_smp))
      allocate(WK_l_mtl%symp_rr(WK_l_mtl%ntsr_lk,np_smp))
      allocate(WK_l_mtl%symn_rt(WK_l_mtl%ntsr_lk,np_smp))
      allocate(WK_l_mtl%symn_rp(WK_l_mtl%ntsr_lk,np_smp))
!
      allocate(WK_l_mtl%asmn_tp(WK_l_mtl%ntsr_lk,np_smp))
      allocate(WK_l_mtl%asmn_dp(WK_l_mtl%ntsr_lk,np_smp))
      allocate(WK_l_mtl%asmp_rp(WK_l_mtl%ntsr_lk,np_smp))
      allocate(WK_l_mtl%asmp_rt(WK_l_mtl%ntsr_lk,np_smp))
!
      end subroutine alloc_leg_tsr_matmul
!
! -----------------------------------------------------------------------
!
      subroutine alloc_leg_vec_matmul(WK_l_mtl)
!
      type(leg_trns_matmul_work), intent(inout) :: WK_l_mtl
!
!
      allocate(WK_l_mtl%pol_e(WK_l_mtl%nvec_jk,np_smp))
      allocate(WK_l_mtl%dpoldt_e(WK_l_mtl%nvec_jk,np_smp))
      allocate(WK_l_mtl%dpoldp_e(WK_l_mtl%nvec_jk,np_smp))
      allocate(WK_l_mtl%dtordt_e(WK_l_mtl%nvec_jk,np_smp))
      allocate(WK_l_mtl%dtordp_e(WK_l_mtl%nvec_jk,np_smp))
!
      allocate(WK_l_mtl%symp_r(WK_l_mtl%nvec_lk,np_smp))
      allocate(WK_l_mtl%symn_t(WK_l_mtl%nvec_lk,np_smp))
      allocate(WK_l_mtl%symn_p(WK_l_mtl%nvec_lk,np_smp))
      allocate(WK_l_mtl%asmp_t(WK_l_mtl%nvec_lk,np_smp))
      allocate(WK_l_mtl%asmp_p(WK_l_mtl%nvec_lk,np_smp))
!
      end subroutine alloc_leg_vec_matmul
!
! -----------------------------------------------------------------------
!
      subroutine alloc_leg_scl_matmul(WK_l_mtl)
!
      type(leg_trns_matmul_work), intent(inout) :: WK_l_mtl
!
!
      allocate(WK_l_mtl%scl_e(WK_l_mtl%nscl_jk,np_smp))
      allocate(WK_l_mtl%symp(WK_l_mtl%nscl_lk,np_smp))
!
      end subroutine alloc_leg_scl_matmul
!
! -----------------------------------------------------------------------
!
      end module t_legendre_work_matmul
