!>@file   t_tsr_leg_trans_sym_matmul.f90
!!@brief  module t_tsr_leg_trans_sym_matmul
!!
!!@author H. Matsui
!!@date Programmed in Aug., 2013
!
!>@brief  Work array for forward Legendre transform useing mat multi 
!>@n      data are strored communication buffer
!!
!!@verbatim
!!      subroutine init_legendre_w_tsr_sym_matmul                       &
!!     &         (sph_rtm, sph_rlm, leg, idx_trns,                      &
!!     &          ntensor, nvector, nscalar, WK_l_sml, WK_tl_sml)
!!      subroutine init_legendre_w_tsr_symmetry                         &
!!     &         (ntensor, sph_rtm, sph_rlm, leg, idx_trns,             &
!!     &          WK_l_sml, WK_tl_sml)
!!        type(sph_rtm_grid), intent(in) :: sph_rtm
!!        type(sph_rlm_grid), intent(in) :: sph_rlm
!!        type(legendre_4_sph_trans), intent(in) :: leg
!!        type(index_4_sph_trans), intent(in) :: idx_trns
!!        type(leg_trns_sym_mul_work), intent(inout) :: WK_l_sml
!!        type(tsr_leg_trns_sym_mul_work), intent(inout) :: WK_tl_sml
!!      subroutine finalize_leg_w_tsr_sym_matmul(WK_l_sml, WK_tl_sml)
!!        type(leg_trns_sym_mul_work), intent(inout) :: WK_l_sml
!!        type(tsr_leg_trns_sym_mul_work), intent(inout) :: WK_tl_sml
!!
!!      subroutine dealloc_tsr_hemi_schmidt_rtm(WK_tl_sml)
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
      module t_tsr_leg_trans_sym_matmul
!
      use m_precision
      use m_constants
      use m_machine_parameter
!
      use calypso_mpi
!
      use t_spheric_rtm_data
      use t_spheric_rlm_data
      use t_schmidt_poly_on_rtm
      use t_work_4_sph_trans
      use t_legendre_work_sym_matmul
!
      use matmul_for_legendre_trans
!
      implicit none
!
!>      Work structure for tensor Legendre trasform
!!        by matmul with symmetry
      type tsr_leg_trns_sym_mul_work
!>        Number of meridional grid points in northern hemisphere
        integer(kind = kint) :: nth_hemi_rtm
!>         @$f P_{l}{m} @$f
!!         at gouss points in northen hemisphere
        real(kind = kreal), allocatable :: P2s_lj(:,:)
!>         @$f dP_{l}{m}/d\theta @$f  with even (l-m) 
!!         at gouss points in northen hemisphere
        real(kind = kreal), allocatable :: dP2sdt_lj(:,:)
!
!>         @$f P_{l}{m} @$f
!!         at gouss points in northen hemisphere
        real(kind = kreal), allocatable :: P2s_jl(:,:)
!>         @$f dP_{l}{m}/d\theta @$f  with even (l-m) 
!!         at gouss points in northen hemisphere
         real(kind = kreal), allocatable :: dP2sdt_jl(:,:)
!
!
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
!!@n
!!@n       @f$e_w@f$ component with 
!!          @f$ (-l(l+1) + (2 m^2 / \sin \theta)) P_{l}^{m} 
!!             - 2 (\cos \theta / \sin^2 \theta) dP_{l}^{m}/d\theta @f$
!!          and  odd (l-m)
        real(kind = kreal), allocatable :: tsrw_o(:,:)
!!@n       @f$e_\chi@f$ component with 
!!          @f$ (-l(l+1) + (2 m^2 / \sin \theta)) P_{l}^{m} 
!!             - 2 (\cos \theta / \sin^2 \theta) dP_{l}^{m}/d\theta @f$
!!          and  odd (l-m)
        real(kind = kreal), allocatable :: tsrchi_o(:,:)
!!@n       @f$e_\tau@f$ component with odd (l-m)
        real(kind = kreal), allocatable :: tsrtau_o(:,:)
!!@n       @f$e_rr@f$ component with odd (l-m)  
        real(kind = kreal), allocatable :: tsrr_o(:,:)
!!@n       Tensor toroidal component with phi derivative and odd (l-m)
        real(kind = kreal), allocatable :: dtsrtdp_o(:,:)
!!@n       Tensor poloidal component with phi derivative and odd (l-m)
        real(kind = kreal), allocatable :: dtsrhdp_o(:,:)
!!@n
!!@n       @f$e_w@f$ component with 
!>          @f$ 2m (-(\cos \theta / \sin^2 \theta)) P_{l}^{m}
!!             + (1 / \sin \theta) dP_{l}^{m}/d\theta @f$
!!          and  odd (l-m)
        real(kind = kreal), allocatable :: dtsrwdt_o(:,:)
!!@n       @f$e_\chi@f$ component with 
!>          @f$ 2m (-(\cos \theta / \sin^2 \theta)) P_{l}^{m}
!!             + (1 / \sin \theta) dP_{l}^{m}/d\theta @f$
!!          and  odd (l-m)
        real(kind = kreal), allocatable :: dtsrchidt_o(:,:)
!!@n       Tensor toroidal component with theta derivative and odd (l-m)
        real(kind = kreal), allocatable :: dtsrtdt_o(:,:)
!!@n       Tensor poloidal component with theta derivative and odd (l-m)
        real(kind = kreal), allocatable :: dtsrhdt_o(:,:)
!
!
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
!
!>         Anti-symmetric tensor P  component
        real(kind = kreal), allocatable :: asmp_dp(:,:)
!>         Anti-symmetric tensor theta-phi component
        real(kind = kreal), allocatable :: asmp_tp(:,:)
!>         Anti-symmetric tensor tau component
        real(kind = kreal), allocatable :: asmp_tau(:,:)
!>         Anti-symmetric tensor r-r component
        real(kind = kreal), allocatable :: asmp_rr(:,:)
!>         Anti-symmetric tensor r-theta-component with condugate order
        real(kind = kreal), allocatable :: asmn_rt(:,:)
!>         Anti-symmetric tensor r-phi-component with condugate order
        real(kind = kreal), allocatable :: asmn_rp(:,:)
!
!>         Symmetric tensor theta-phi component with condugate order
        real(kind = kreal), allocatable :: symn_tp(:,:)
!>         Symmetric tensor P  component with condugate order
        real(kind = kreal), allocatable :: symn_dp(:,:)
!>         Symmetric tensor r-phi-component
        real(kind = kreal), allocatable :: symp_rp(:,:)
!>         Symmetric tensor r-theta-component
        real(kind = kreal), allocatable :: symp_rt(:,:)
      end type tsr_leg_trns_sym_mul_work
!
      private :: const_tsr_symmetric_legendres
      private :: alloc_tsr_hemi_schmidt_rtm
      private :: alloc_leg_tsr_sym_matmul, dealloc_leg_tsr_sym_matmul
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine init_legendre_w_tsr_sym_matmul                         &
     &         (sph_rtm, sph_rlm, leg, idx_trns,                        &
     &          ntensor, nvector, nscalar, WK_l_sml, WK_tl_sml)
!
      type(sph_rtm_grid), intent(in) :: sph_rtm
      type(sph_rlm_grid), intent(in) :: sph_rlm
      type(legendre_4_sph_trans), intent(in) :: leg
      type(index_4_sph_trans), intent(in) :: idx_trns
      integer(kind = kint), intent(in) :: ntensor, nvector, nscalar
!
      type(leg_trns_sym_mul_work), intent(inout) :: WK_l_sml
      type(tsr_leg_trns_sym_mul_work), intent(inout) :: WK_tl_sml
!
!
      call init_legendre_sym_matmul(sph_rtm, sph_rlm, leg,              &
     &    idx_trns, nvector, nscalar, WK_l_sml)
      if(ntensor .le. 0) return
!
      call const_tsr_symmetric_legendres                                &
     &   (sph_rlm%nidx_rlm(2), sph_rlm%idx_gl_1d_rlm_j,                 &
     &    sph_rtm%nidx_rtm(2), sph_rtm%nidx_rtm(3),                     &
     &    leg, idx_trns, WK_l_sml, WK_tl_sml)
!
      WK_tl_sml%ntsr_jk = ((idx_trns%maxdegree_rlm+1)/2)                &
     &                  * sph_rtm%maxidx_rtm_smp(1) * ntensor
      WK_tl_sml%ntsr_lk = ((sph_rtm%nidx_rtm(2) + 1)/2)                 &
     &                  * sph_rtm%maxidx_rtm_smp(1) * ntensor
      call alloc_leg_tsr_sym_matmul(WK_tl_sml)
!
      end subroutine init_legendre_w_tsr_sym_matmul
!
! -----------------------------------------------------------------------
!
      subroutine init_legendre_w_tsr_symmetry                           &
     &         (ntensor, sph_rtm, sph_rlm, leg, idx_trns,               &
     &          WK_l_sml, WK_tl_sml)
!
      integer(kind = kint), intent(in) :: ntensor
      type(sph_rtm_grid), intent(in) :: sph_rtm
      type(sph_rlm_grid), intent(in) :: sph_rlm
      type(legendre_4_sph_trans), intent(in) :: leg
      type(index_4_sph_trans), intent(in) :: idx_trns
!
      type(leg_trns_sym_mul_work), intent(inout) :: WK_l_sml
      type(tsr_leg_trns_sym_mul_work), intent(inout) :: WK_tl_sml
!
!
      call init_legendre_symmetry                                       &
     &         (sph_rtm, sph_rlm, leg, idx_trns, WK_l_sml)
      if(ntensor .le. 0) return
!
      call const_tsr_symmetric_legendres                                &
     &   (sph_rlm%nidx_rlm(2), sph_rlm%idx_gl_1d_rlm_j,                 &
     &    sph_rtm%nidx_rtm(2), sph_rtm%nidx_rtm(3),                     &
     &    leg, idx_trns, WK_l_sml, WK_tl_sml)
!
      WK_tl_sml%ntsr_jk = (idx_trns%maxdegree_rlm+1)/2
      WK_tl_sml%ntsr_lk = (sph_rtm%nidx_rtm(2) + 1)/2
      call alloc_leg_tsr_sym_matmul(WK_tl_sml)
!
      end subroutine init_legendre_w_tsr_symmetry
!
! -----------------------------------------------------------------------
!
      subroutine finalize_leg_w_tsr_sym_matmul(WK_l_sml, WK_tl_sml)
!
      type(leg_trns_sym_mul_work), intent(inout) :: WK_l_sml
      type(tsr_leg_trns_sym_mul_work), intent(inout) :: WK_tl_sml
!
!
      call finalize_legendre_sym_matmul(WK_l_sml)
!
      if(allocated(WK_tl_sml%P2s_lj)) then
        call dealloc_leg_tsr_sym_matmul(WK_tl_sml)
        call dealloc_tsr_hemi_schmidt_rtm(WK_tl_sml)
      end if
!
      end subroutine finalize_leg_w_tsr_sym_matmul
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine const_tsr_symmetric_legendres                          &
     &         (jmax_rlm, idx_gl_1d_rlm_j, nth_rtm, mphi_rtm,           &
     &          leg, idx_trns, WK_l_sml, WK_tl_sml)
!
      use set_legendre_matrices
!
      integer(kind = kint), intent(in) :: nth_rtm, mphi_rtm, jmax_rlm
      integer(kind = kint), intent(in) :: idx_gl_1d_rlm_j(jmax_rlm,3)
      type(legendre_4_sph_trans), intent(in) :: leg
      type(index_4_sph_trans), intent(in) :: idx_trns
      type(leg_trns_sym_mul_work), intent(in) :: WK_l_sml
!
      type(tsr_leg_trns_sym_mul_work), intent(inout) :: WK_tl_sml
!
!
      call alloc_tsr_hemi_schmidt_rtm(jmax_rlm, WK_l_sml, WK_tl_sml)
!
      call symmetric_legendre_4_tensor_lj                               &
     &  (nth_rtm, mphi_rtm, jmax_rlm, WK_tl_sml%nth_hemi_rtm,           &
     &   idx_trns%lstack_rlm, idx_trns%lstack_even_rlm,                 &
     &   idx_gl_1d_rlm_j, leg%g_sph_rlm, leg%asin_t_rtm, leg%cos_t_rtm, &
     &   WK_l_sml%Ps_rtm, WK_l_sml%dPsdt_rtm,                           &
     &   WK_tl_sml%P2s_lj, WK_tl_sml%dP2sdt_lj)
      call symmetric_legendre_4_tensor_jl                               &
     &  (nth_rtm, mphi_rtm, jmax_rlm, WK_tl_sml%nth_hemi_rtm,           &
     &   idx_trns%lstack_rlm, idx_trns%lstack_even_rlm,                 &
     &   idx_gl_1d_rlm_j, leg%g_sph_rlm, leg%asin_t_rtm, leg%cos_t_rtm, &
     &   WK_l_sml%Ps_jl, WK_l_sml%dPsdt_jl,                             &
     &   WK_tl_sml%P2s_jl, WK_tl_sml%dP2sdt_jl)
!
      end subroutine const_tsr_symmetric_legendres
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine alloc_tsr_hemi_schmidt_rtm                             &
     &         (jmax_rlm, WK_l_sml, WK_tl_sml)
!
      integer(kind = kint), intent(in) :: jmax_rlm
      type(leg_trns_sym_mul_work), intent(in) :: WK_l_sml
!
      type(tsr_leg_trns_sym_mul_work), intent(inout) :: WK_tl_sml
!
!
      WK_tl_sml%nth_hemi_rtm = WK_l_sml%nth_hemi_rtm
      allocate( WK_tl_sml%P2s_lj(WK_tl_sml%nth_hemi_rtm,jmax_rlm) )
      allocate( WK_tl_sml%dP2sdt_lj(WK_tl_sml%nth_hemi_rtm,jmax_rlm) )
!
      allocate( WK_tl_sml%P2s_jl(jmax_rlm,WK_tl_sml%nth_hemi_rtm) )
      allocate( WK_tl_sml%dP2sdt_jl(jmax_rlm,WK_tl_sml%nth_hemi_rtm) )
!
      WK_tl_sml%P2s_lj =    0.0d0
      WK_tl_sml%dP2sdt_lj = 0.0d0
!
      WK_tl_sml%P2s_jl =    0.0d0
      WK_tl_sml%dP2sdt_jl = 0.0d0
!
      end subroutine alloc_tsr_hemi_schmidt_rtm
!
! -----------------------------------------------------------------------
!
      subroutine dealloc_tsr_hemi_schmidt_rtm(WK_tl_sml)
!
      type(tsr_leg_trns_sym_mul_work), intent(inout) :: WK_tl_sml
!
!
      deallocate(WK_tl_sml%P2s_lj, WK_tl_sml%dP2sdt_lj)
      deallocate(WK_tl_sml%P2s_jl,  WK_tl_sml%dP2sdt_jl)
!
      end subroutine dealloc_tsr_hemi_schmidt_rtm
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine alloc_leg_tsr_sym_matmul(WK_tl_sml)
!
      type(tsr_leg_trns_sym_mul_work), intent(inout) :: WK_tl_sml
!
!
      allocate(WK_tl_sml%tsrw_e(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%tsrchi_e(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%tsrtau_e(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%tsrr_e(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%dtsrtdp_e(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%dtsrhdp_e(WK_tl_sml%ntsr_jk,np_smp))
!
      allocate(WK_tl_sml%dtsrwdt_e(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%dtsrchidt_e(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%dtsrtdt_e(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%dtsrhdt_e(WK_tl_sml%ntsr_jk,np_smp))
!
      allocate(WK_tl_sml%tsrw_o(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%tsrchi_o(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%tsrtau_o(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%tsrr_o(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%dtsrtdp_o(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%dtsrhdp_o(WK_tl_sml%ntsr_jk,np_smp))
!
      allocate(WK_tl_sml%dtsrwdt_o(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%dtsrchidt_o(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%dtsrtdt_o(WK_tl_sml%ntsr_jk,np_smp))
      allocate(WK_tl_sml%dtsrhdt_o(WK_tl_sml%ntsr_jk,np_smp))
!
!
      allocate(WK_tl_sml%symp_dp(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%symp_tp(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%symp_tau(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%symp_rr(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%symn_rt(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%symn_rp(WK_tl_sml%ntsr_lk,np_smp))
!
      allocate(WK_tl_sml%asmn_tp(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%asmn_dp(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%asmp_rp(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%asmp_rt(WK_tl_sml%ntsr_lk,np_smp))
!
      allocate(WK_tl_sml%asmp_dp(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%asmp_tp(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%asmp_tau(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%asmp_rr(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%asmp_rt(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%asmp_rp(WK_tl_sml%ntsr_lk,np_smp))
!
      allocate(WK_tl_sml%symn_tp(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%symn_dp(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%symp_rp(WK_tl_sml%ntsr_lk,np_smp))
      allocate(WK_tl_sml%symp_rt(WK_tl_sml%ntsr_lk,np_smp))
!
      end subroutine alloc_leg_tsr_sym_matmul
!
! -----------------------------------------------------------------------
!
      subroutine dealloc_leg_tsr_sym_matmul(WK_tl_sml)
!
      type(tsr_leg_trns_sym_mul_work), intent(inout) :: WK_tl_sml
!
!
      deallocate(WK_tl_sml%tsrw_e,    WK_tl_sml%tsrchi_e)
      deallocate(WK_tl_sml%tsrtau_e,  WK_tl_sml%tsrr_e)
      deallocate(WK_tl_sml%dtsrtdp_e, WK_tl_sml%dtsrhdp_e)
      deallocate(WK_tl_sml%dtsrwdt_e, WK_tl_sml%dtsrchidt_e)
      deallocate(WK_tl_sml%dtsrtdt_e, WK_tl_sml%dtsrhdt_e)
!
      deallocate(WK_tl_sml%tsrw_o,    WK_tl_sml%tsrchi_o)
      deallocate(WK_tl_sml%tsrtau_o,  WK_tl_sml%tsrr_o)
      deallocate(WK_tl_sml%dtsrtdp_o, WK_tl_sml%dtsrhdp_o)
      deallocate(WK_tl_sml%dtsrwdt_o, WK_tl_sml%dtsrchidt_o)
      deallocate(WK_tl_sml%dtsrtdt_o, WK_tl_sml%dtsrhdt_o)
!
      deallocate(WK_tl_sml%symp_dp,  WK_tl_sml%symp_tp)
      deallocate(WK_tl_sml%symp_tau, WK_tl_sml%symp_rr)
      deallocate(WK_tl_sml%symn_rt,  WK_tl_sml%symn_rp)
      deallocate(WK_tl_sml%asmn_tp,  WK_tl_sml%asmn_dp)
      deallocate(WK_tl_sml%asmp_rp,  WK_tl_sml%asmp_rt)
!
      deallocate(WK_tl_sml%asmp_dp,  WK_tl_sml%asmp_tp)
      deallocate(WK_tl_sml%asmp_tau, WK_tl_sml%asmp_rr)
      deallocate(WK_tl_sml%asmp_rt,  WK_tl_sml%asmp_rp)
      deallocate(WK_tl_sml%symn_tp,  WK_tl_sml%symn_dp)
      deallocate(WK_tl_sml%symp_rp,  WK_tl_sml%symp_rt)
!
      end subroutine dealloc_leg_tsr_sym_matmul
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      end module t_tsr_leg_trans_sym_matmul
