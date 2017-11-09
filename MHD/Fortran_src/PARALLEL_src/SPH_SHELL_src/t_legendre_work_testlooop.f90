!>@file   t_legendre_work_testlooop.f90
!!@brief  module t_legendre_work_testlooop
!!
!!@author H. Matsui
!!@date Programmed in Aug., 2013
!
!>@brief  Work array for forward Legendre transform useing mat multi 
!>@n      data are strored communication buffer
!!
!!@verbatim
!!      subroutine init_legendre_testloop(sph_rtm, sph_rlm, leg,        &
!!     &          idx_trns, ntensor, nvector, nscalar, WK_l_tst)
!!        type(sph_rtm_grid), intent(in) :: sph_rtm
!!        type(sph_rlm_grid), intent(in) :: sph_rlm
!!        type(legendre_4_sph_trans), intent(in) :: leg
!!        type(index_4_sph_trans), intent(in) :: idx_trns
!!        type(leg_trns_testloop_work), intent(inout) :: WK_l_tst
!!
!!      subroutine alloc_leg_vec_test                                   &
!!     &         (nri_rtm, maxidx_rtm_t_smp, nvector, nscalar, WK_l_tst)
!!      subroutine dealloc_leg_vec_test(WK_l_tst)
!!      subroutine dealloc_leg_scl_test(WK_l_tst)
!!        type(leg_trns_testloop_work), intent(inout) :: WK_l_tst
!!
!!     field data for Legendre transform
!!       original layout: vr_rtm(l_rtm,m_rtm,k_rtm,icomp)
!!       size: vr_rtm(nidx_rtm(2),nidx_rtm(1)*ncomp,nidx_rtm(3))
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
      module t_legendre_work_testlooop
!
      use m_precision
      use m_constants
      use calypso_mpi
!
      use m_machine_parameter
      use matmul_for_legendre_trans
!
      use t_spheric_rtm_data
      use t_spheric_rlm_data
      use t_schmidt_poly_on_rtm
      use t_work_4_sph_trans
!
      implicit none
!
!>      Structure for Legendre trasdorm for test
      type leg_trns_testloop_work
!>        Number of meridional grid points in northern hemisphere
        integer(kind = kint) :: nth_sym
!>          @$f P_{l}{m} @$f
!!          at gouss points in northen hemisphere
        real(kind = kreal), allocatable :: Ps_tj(:,:)
!>          @$f dP_{l}{m}/d\theta @$f  with even (l-m) 
!!          at gouss points in northen hemisphere
        real(kind = kreal), allocatable :: dPsdt_tj(:,:)
!>          @f$ (-l(l+1) + (2 m^2 / \sin \theta)) P_{l}^{m} 
!!             - 2 (\cos \theta / \sin^2 \theta) dP_{l}^{m}/d\theta @f$
!!          at gouss points in northen hemisphere
        real(kind = kreal), allocatable :: P2s_tj(:,:)
!>          @f$ 2m (-(\cos \theta / \sin^2 \theta)) P_{l}^{m}
!!             + (1 / \sin \theta) dP_{l}^{m}/d\theta @f$
!!          at gouss points in northen hemisphere
        real(kind = kreal), allocatable :: dP2sdt_tj(:,:)
!
!
!>        Maximum matrix size for spectr tensor
        integer(kind = kint) :: ntsr_jk
!>       Maximum matrix size for spectr data
        integer(kind = kint) :: nvec_jk
!>       Maximum matrix size for spectr data
        integer(kind = kint) :: nscl_jk
!
!>       size for work area of pol_e and pol_o
        integer(kind = kint) :: n_pol_e
!>       size for work area of tor_e and tor_o
        integer(kind = kint) :: n_tor_e
!
!!@n       @f$e_w@f$ component with 
!!          @f$ (-l(l+1) + (2 m^2 / \sin \theta)) P_{l}^{m} 
!!             - 2 (\cos \theta / \sin^2 \theta) dP_{l}^{m}/d\theta @f$
!!          and  evem (l-m)
!!@n        real(kind = kreal), allocatable :: tsrw_e(:,:)
!!@n       @f$e_\chi@f$ component with 
!!          @f$ (-l(l+1) + (2 m^2 / \sin \theta)) P_{l}^{m} 
!!             - 2 (\cos \theta / \sin^2 \theta) dP_{l}^{m}/d\theta @f$
!!          and  evem (l-m)
!!@n        real(kind = kreal), allocatable :: tsrchi_e(:,:)
!!@n       @f$e_\tau@f$ component with evem (l-m)
!!@n        real(kind = kreal), allocatable :: dtsrtau_e(:,:)
!!@n       @f$e_rr@f$ component with evem (l-m)  
!!@n        real(kind = kreal), allocatable :: tsrr_e(:,:)
!!@n       Tensor toroidal component with phi derivative and evem (l-m)
!!@n        real(kind = kreal), allocatable :: dtsrtdp_e(:,:)
!!@n       Tensor poloidal component with phi derivative and evem (l-m)
!!@n        real(kind = kreal), allocatable :: dtsrhdp_e(:,:)
!!@n
!!@n       Poloidal component with evem (l-m)  
!!@n        real(kind = kreal), allocatable :: pol_e(:,:)
!!@n       Toroidal component with phi derivative and evem (l-m)
!!@n        real(kind = kreal), allocatable :: dtordp_e(:,:)
!!@n       Poloidal component with phi derivative and evem (l-m)
!!@n        real(kind = kreal), allocatable :: dpoldp_e(:,:)
!!@n       Scalar with evem (l-m)
!!@n        real(kind = kreal), allocatable :: scl_e(:,:)
!!@n
!!@n       tsrw_e =    Pol_e(          1:  ntsr_jk,ip)
!!@n       tsrchi_e =  Pol_e(  ntsr_jk+1:2*ntsr_jk,ip)
!!@n       dtsrtau_e = Pol_e(2*ntsr_jk+1:3*ntsr_jk,ip)
!!@n       tsrr_e =    Pol_e(3*ntsr_jk+1:4*ntsr_jk,ip)
!!@n       dtsrtdp_e = Pol_e(4*ntsr_jk+1:5*ntsr_jk,ip)
!!@n       dtsrhdp_e = Pol_e(5*ntsr_jk+1:6*ntsr_jk,ip)
!!@n
!!@n       pol_e =    Pol_e(          6*ntsr_jk+1:  nvec_jk+6*ntsr_jk,ip)
!!@n       dtordp_e = Pol_e(  nvec_jk+6*ntsr_jk+1:2*nvec_jk+6*ntsr_jk,ip)
!!@n       dpoldp_e = Pol_e(2*nvec_jk+6*ntsr_jk+1:3*nvec_jk+6*ntsr_jk,ip)
!!@n
!!@n       scl_e 
!!@n        = Pol_e(3*nvec_jk+6*ntsr_jk+1:nscl_jk+3*nvec_jk+6*ntsr_jk,ip)
        real(kind = kreal), allocatable :: pol_e(:,:)
!
!!@n       @f$e_w@f$ component with 
!>          @f$ 2m (-(\cos \theta / \sin^2 \theta)) P_{l}^{m}
!!             + (1 / \sin \theta) dP_{l}^{m}/d\theta @f$
!!          and  evem (l-m)
!!@n        real(kind = kreal), allocatable :: dtsrwdt_e(:,:)
!!@n       @f$e_\chi@f$ component with 
!>          @f$ 2m (-(\cos \theta / \sin^2 \theta)) P_{l}^{m}
!!             + (1 / \sin \theta) dP_{l}^{m}/d\theta @f$
!!          and  evem (l-m)
!!@n        real(kind = kreal), allocatable :: dtsrchidt_e(:,:)
!!@n       Tensor toroidal component with theta derivative and evem (l-m)
!!@n        real(kind = kreal), allocatable :: dtsrtdt_e(:,:)
!!@n       Tensor poloidal component with theta derivative and evem (l-m)
!!@n        real(kind = kreal), allocatable :: dtsrhdt_e(:,:)
!!@n
!!@n       Toroidal component with theta derivative and evem (l-m)
!!@n        real(kind = kreal), allocatable :: dtordt_e(:,:)
!!@n       Poloidal component with theta derivative and evem (l-m)
!!@n        real(kind = kreal), allocatable :: dpoldt_e(:,:)
!!
!!@n       dtsrwdt_e = tor_e(          1:  ntsr_jk,ip)
!!@n       dtsrwdt_e = tor_e(  ntsr_jk+1:2*ntsr_jk,ip)
!!@n       dtsrtdt_e = tor_e(2*ntsr_jk+1:3*ntsr_jk,ip)
!!@n       dtsrhdt_e = tor_e(3*ntsr_jk+1:4*ntsr_jk,ip)
!!@n
!!@n       dtordt_e = tor_e(          4*ntsr_jk+1:  nvec_jk+4*ntsr_jk,ip)
!!@n       dpoldt_e = tor_e(  nvec_jk+4*ntsr_jk+1:2*nvec_jk+4*ntsr_jk,ip)
        real(kind = kreal), allocatable :: tor_e(:,:)
!
!!@n       @f$e_w@f$ component with 
!!          @f$ (-l(l+1) + (2 m^2 / \sin \theta)) P_{l}^{m} 
!!             - 2 (\cos \theta / \sin^2 \theta) dP_{l}^{m}/d\theta @f$
!!          and  odd (l-m)
!!@n        real(kind = kreal), allocatable :: tsrw_o(:,:)
!!@n       @f$e_\chi@f$ component with 
!!          @f$ (-l(l+1) + (2 m^2 / \sin \theta)) P_{l}^{m} 
!!             - 2 (\cos \theta / \sin^2 \theta) dP_{l}^{m}/d\theta @f$
!!          and  odd (l-m)
!!@n        real(kind = kreal), allocatable :: tsrchi_o(:,:)
!!@n       @f$e_\tau@f$ component with odd (l-m)
!!@n        real(kind = kreal), allocatable :: dtsrtau_o(:,:)
!!@n       @f$e_rr@f$ component with odd (l-m)  
!!@n        real(kind = kreal), allocatable :: tsrr_o(:,:)
!!@n       Tensor toroidal component with phi derivative and odd (l-m)
!!@n        real(kind = kreal), allocatable :: dtsrtdp_o(:,:)
!!@n       Tensor poloidal component with phi derivative and odd (l-m)
!!@n        real(kind = kreal), allocatable :: dtsrhdp_o(:,:)
!!@n
!!@n       Poloidal component with odd (l-m)  
!!@n        real(kind = kreal), allocatable :: pol_o(:,:)
!!@n       Toroidal component with phi derivative and odd (l-m)
!!@n        real(kind = kreal), allocatable :: dtordp_o(:,:)
!!@n       Poloidal component with phi derivative and odd (l-m)
!!@n        real(kind = kreal), allocatable :: dpoldp_o(:,:)
!!@n       Scalar with odd (l-m)
!!@n        real(kind = kreal), allocatable :: scl_o(:,:)
!!@n
!!@n       tsrw_o =    Pol_o(          1:  ntsr_jk,ip)
!!@n       tsrchi_o =  Pol_o(  ntsr_jk+1:2*ntsr_jk,ip)
!!@n       dtsrtau_o = Pol_o(2*ntsr_jk+1:3*ntsr_jk,ip)
!!@n       tsrr_o =    Pol_o(3*ntsr_jk+1:4*ntsr_jk,ip)
!!@n       dtsrtdp_o = Pol_o(4*ntsr_jk+1:5*ntsr_jk,ip)
!!@n       dtsrhdp_o = Pol_o(5*ntsr_jk+1:6*ntsr_jk,ip)
!!@n
!!@n       pol_o =    Pol_o(          6*ntsr_jk+1:  nvec_jk+6*ntsr_jk,ip)
!!@n       dtordp_o = Pol_o(  nvec_jk+6*ntsr_jk+1:2*nvec_jk+6*ntsr_jk,ip)
!!@n       dpoldp_o = Pol_o(2*nvec_jk+6*ntsr_jk+1:3*nvec_jk+6*ntsr_jk,ip)
!!@n
!!@n       scl_o
!!@n        = Pol_o(3*nvec_jk+6*ntsr_jk+1:nscl_jk+3*nvec_jk+6*ntsr_jk,ip)
        real(kind = kreal), allocatable :: pol_o(:,:)
!
!!@n       @f$e_w@f$ component with 
!>          @f$ 2m (-(\cos \theta / \sin^2 \theta)) P_{l}^{m}
!!             + (1 / \sin \theta) dP_{l}^{m}/d\theta @f$
!!          and  odd (l-m)
!!@n        real(kind = kreal), allocatable :: dtsrwdt_o(:,:)
!!@n       @f$e_\chi@f$ component with 
!>          @f$ 2m (-(\cos \theta / \sin^2 \theta)) P_{l}^{m}
!!             + (1 / \sin \theta) dP_{l}^{m}/d\theta @f$
!!          and  odd (l-m)
!!@n        real(kind = kreal), allocatable :: dtsrchidt_o(:,:)
!!@n       Tensor toroidal component with theta derivative and odd (l-m)
!!@n        real(kind = kreal), allocatable :: dtsrtdt_o(:,:)
!!@n       Tensor poloidal component with theta derivative and odd (l-m)
!!@n        real(kind = kreal), allocatable :: dtsrhdt_o(:,:)
!!@n
!!@n       Toroidal component with theta derivative and odd (l-m)
!!@n        real(kind = kreal), allocatable :: dtordt_o(:,:)
!!@n       Poloidal component with theta derivative and odd (l-m)
!!@n        real(kind = kreal), allocatable :: dpoldt_o(:,:)
!!
!!@n       dtsrwdt_o = tor_o(          1:  ntsr_jk,ip)
!!@n       dtsrwdt_o = tor_o(  ntsr_jk+1:2*ntsr_jk,ip)
!!@n       dtsrtdt_o = tor_o(2*ntsr_jk+1:3*ntsr_jk,ip)
!!@n       dtsrhdt_o = tor_o(3*ntsr_jk+1:4*ntsr_jk,ip)
!!@n
!!@n       dtordt_o = tor_o(          4*ntsr_jk+1:  nvec_jk+4*ntsr_jk,ip)
!!@n       dpoldt_o = tor_o(  nvec_jk+4*ntsr_jk+1:2*nvec_jk+4*ntsr_jk,ip)
        real(kind = kreal), allocatable :: tor_o(:,:)
!
!
!>       Maximum matrix size for symmeytric tensor field data
        integer(kind = kint) :: ntsr_lk
!>       Maximum matrix size for vector field data
        integer(kind = kint) :: nvec_lk
!>       Maximum matrix size for scalar field data
        integer(kind = kint) :: nscl_lk
!
!>       size for work area of symp_r and asmp_r
        integer(kind = kint) :: n_sym_r
!>       size for work area of symp_p and asmp_p
        integer(kind = kint) :: n_sym_p
!
!>         Symmetric tensor P  component
!!@n        real(kind = kreal), allocatable :: symp_dp(:,:)
!!@n       Symmetric tensor theta-phi component
!!@n        real(kind = kreal), allocatable :: symp_tp(:,:)
!!@n       Symmetric tensor tau component
!!@n        real(kind = kreal), allocatable :: symp_tau(:,:)
!!@n       Symmetric tensor r-r component
!!@n        real(kind = kreal), allocatable :: symp_rr(:,:)
!!@n       Symmetric tensor r-theta-component with condugate order
!!@n        real(kind = kreal), allocatable :: symn_rt(:,:)
!!@n       Symmetric tensor r-phi-component with condugate order
!!@n        real(kind = kreal), allocatable :: symn_rp(:,:)
!!@n
!!@n       Symmetric radial component
!!@n        real(kind = kreal), allocatable :: symp_r(:,:)
!!@n       Symmetric theta-component with condugate order
!!@n        real(kind = kreal), allocatable :: symn_t(:,:)
!!@n       Symmetric phi-component with condugate order
!!@n        real(kind = kreal), allocatable :: symn_p(:,:)
!!@n
!!@n       Symmetric scalar component
!!@n        real(kind = kreal), allocatable :: symp(:,:)
!!@n
!!@n       symp_dp =  symp_r(          1:  ntsr_lk,ip)
!!@n       symp_tp =  symp_r(  ntsr_lk+1:2*ntsr_lk,ip)
!!@n       symp_tau = symp_r(2*ntsr_lk+1:3*ntsr_lk,ip)
!!@n       symp_rr =  symp_r(3*ntsr_lk+1:4*ntsr_lk,ip)
!!@n       symn_rt =  symp_r(4*ntsr_lk+1:5*ntsr_lk,ip)
!!@n       symn_rp =  symp_r(5*ntsr_lk+1:6*ntsr_lk,ip)
!!@n
!!@n       symp_r = symp_r(          6*ntsr_lk+1:  nvec_lk+6*ntsr_lk,ip)
!!@n       symn_t = symp_r(  nvec_lk+6*ntsr_lk+1:2*nvec_lk+6*ntsr_lk,ip)
!!@n       symn_p = symp_r(2*nvec_lk+6*ntsr_lk+1:3*nvec_lk+6*ntsr_lk,ip)
!!@n
!!@n       symp
!!@n       = symp_r(3*nvec_lk+6*ntsr_lk+1:3*nvec_lk+6*ntsr_lk+nscl_lk,ip)
        real(kind = kreal), allocatable :: symp_r(:,:)
!
!>         Anti-symmetric tensor theta-phi component with condugate order
!!@n        real(kind = kreal), allocatable :: asmn_tp(:,:)
!!@n       Anti-symmetric tensor P  component with condugate order
!!@n        real(kind = kreal), allocatable :: asmn_dp(:,:)
!!@n       Anti-symmetric tensor r-phi-component
!!@n        real(kind = kreal), allocatable :: asmp_rp(:,:)
!!@n       Anti-symmetric tensor r-theta-component
!!@n        real(kind = kreal), allocatable :: asmp_rt(:,:)
!!@n
!!@n       Anti-symmetric phi-component
!!@n        real(kind = kreal), allocatable :: asmp_p(:,:)
!!@n       Anti-symmetric theta-component
!!@n        real(kind = kreal), allocatable :: asmp_t(:,:)
!!@n
!!@n       asmn_tp = asmp_p(          1:  ntsr_lk,ip)
!!@n       asmn_dp = asmp_p(  ntsr_lk+1:2*ntsr_lk,ip)
!!@n       asmp_rp = asmp_p(2*ntsr_lk+1:3*ntsr_lk,ip)
!!@n       asmp_rt = asmp_p(3*ntsr_lk+1:4*ntsr_lk,ip)
!!@n
!!@n       asmp_p = asmp_p(          1:  nvec_lk,ip)
!!@n       asmp_t = asmp_p(  nvec_lk+1:2*nvec_lk,ip)
        real(kind = kreal), allocatable :: asmp_p(:,:)
!
!>         Anti-symmetric tensor P  component
!!@n        real(kind = kreal), allocatable :: asmp_dp(:,:)
!!@n       Anti-symmetric tensor theta-phi component
!!@n        real(kind = kreal), allocatable :: asmp_tau(:,:)
!!@n       Anti-symmetric tensor tau component
!!@n        real(kind = kreal), allocatable :: asmp_tau(:,:)
!!@n       Anti-symmetric tensor r-r component
!!@n        real(kind = kreal), allocatable :: asmp_rr(:,:)
!!@n       Anti-symmetric tensor r-theta-component with condugate order
!!@n        real(kind = kreal), allocatable :: asmp_rt(:,:)
!!@n       Anti-symmetric tensor r-phi-component with condugate order
!!@n        real(kind = kreal), allocatable :: asmp_rp(:,:)
!!@n
!!@n       Anti-symmetric radial component
!!@n        real(kind = kreal), allocatable :: asmp_r(:,:)
!!@n       Anti-symmetric theta-component with condugate order
!!@n        real(kind = kreal), allocatable :: asmn_t(:,:)
!!@n       Anti-symmetric phi-component with condugate order
!!@n        real(kind = kreal), allocatable :: asmn_p(:,:)
!!@n
!!@n       Anti-symmetric scalar component
!!@n        real(kind = kreal), allocatable :: asmp(:,:)
!!@n
!!@n       asmp_dp =  asmp_r(          1:  ntsr_lk,ip)
!!@n       asmp_tp =  asmp_r(  ntsr_lk+1:2*ntsr_lk,ip)
!!@n       asmp_tau = asmp_r(2*ntsr_lk+1:3*ntsr_lk,ip)
!!@n       asmp_rr =  asmp_r(3*ntsr_lk+1:4*ntsr_lk,ip)
!!@n       asmp_rp =  asmp_r(4*ntsr_lk+1:5*ntsr_lk,ip)
!!@n       asmp_rp =  asmp_r(5*ntsr_lk+1:6*ntsr_lk,ip)
!!@n
!!@n       asmp_r = asmp_r(          6*ntsr_lk+1:  nvec_lk+6*ntsr_lk,ip)
!!@n       asmp_t = asmp_r(  nvec_lk+6*ntsr_lk+1:2*nvec_lk+6*ntsr_lk,ip)
!!@n       asmp_p = asmp_r(2*nvec_lk+6*ntsr_lk+1:3*nvec_lk+6*ntsr_lk,ip)
!!@n
!!@n       asmp
!!@n       = asmp_r(3*nvec_lk+6*ntsr_lk+1:3*nvec_lk+6*ntsr_lk+nscl_lk,ip)
        real(kind = kreal), allocatable :: asmp_r(:,:)
!
!>         Symmetric tensor theta-phi component with condugate order
!!@n        real(kind = kreal), allocatable :: symn_tp(:,:)
!!@n       Symmetric tensor P  component with condugate order
!!@n        real(kind = kreal), allocatable :: symn_dp(:,:)
!!@n       Symmetric tensor r-phi-component
!!@n        real(kind = kreal), allocatable :: symp_rp(:,:)
!!@n       Symmetric tensor r-theta-component
!!@n        real(kind = kreal), allocatable :: symp_rt(:,:)
!!@n
!!@n       Symmetric phi-component
!!@n        real(kind = kreal), allocatable :: symp_p(:,:)
!!@n       Symmetric theta-component
!!@n        real(kind = kreal), allocatable :: symp_t(:,:)
!!@n
!!@n       symn_tp = symp_p(          1:  ntsr_lk,ip)
!!@n       symn_dt = symp_p(  ntsr_lk+1:2*ntsr_lk,ip)
!!@n       symp_rp = symp_p(2*ntsr_lk+1:3*ntsr_lk,ip)
!!@n       symp_rt = symp_p(3*ntsr_lk+1:4*ntsr_lk,ip)
!!@n
!!@n       symp_p = symp_p(          4*ntsr_lk+1:  nvec_lk+4*ntsr_lk,ip)
!!@n       symp_t = symp_p(  nvec_lk+4*ntsr_lk*1:2*nvec_lk+4*ntsr_lk,ip)
        real(kind = kreal), allocatable :: symp_p(:,:)
      end type leg_trns_testloop_work
!
      private :: const_legendre_testloop, alloc_leg_vec_test
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine init_legendre_testloop(sph_rtm, sph_rlm, leg,          &
     &          idx_trns, ntensor, nvector, nscalar, WK_l_tst)
!
      type(sph_rtm_grid), intent(in) :: sph_rtm
      type(sph_rlm_grid), intent(in) :: sph_rlm
      type(legendre_4_sph_trans), intent(in) :: leg
      type(index_4_sph_trans), intent(in) :: idx_trns
      integer(kind = kint), intent(in) :: ntensor, nvector, nscalar
!
      type(leg_trns_testloop_work), intent(inout) :: WK_l_tst
!
!
      call const_legendre_testloop                                      &
     &   (ntensor, sph_rlm%nidx_rlm(2), sph_rlm%idx_gl_1d_rlm_j,        &
     &    sph_rtm%nidx_rtm(2), sph_rtm%nidx_rtm(3),                     &
     &    leg, idx_trns, WK_l_tst)
      call alloc_leg_vec_test                                           &
     &   (sph_rtm%nidx_rtm(2), sph_rtm%maxidx_rtm_smp(1),               &
     &    ntensor, nvector, nscalar, idx_trns, WK_l_tst)
!
      end subroutine init_legendre_testloop
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine const_legendre_testloop                                &
     &         (ntensor, jmax_rlm, idx_gl_1d_rlm_j, nth_rtm, mphi_rtm,  &
     &          leg, idx_trns, WK_l_tst)
!
      use set_legendre_matrices
!
      integer(kind = kint), intent(in) :: ntensor
      integer(kind = kint), intent(in) :: nth_rtm, mphi_rtm, jmax_rlm
      integer(kind = kint), intent(in) :: idx_gl_1d_rlm_j(jmax_rlm,3)
      type(legendre_4_sph_trans), intent(in) :: leg
      type(index_4_sph_trans), intent(in) :: idx_trns
!
      type(leg_trns_testloop_work), intent(inout) :: WK_l_tst
!
!
      WK_l_tst%nth_sym = (nth_rtm+1) / 2
      allocate( WK_l_tst%Ps_tj(WK_l_tst%nth_sym,jmax_rlm) )
      allocate( WK_l_tst%dPsdt_tj(WK_l_tst%nth_sym,jmax_rlm) )
!
      call set_symmetric_legendre_lj                                    &
     &   (nth_rtm, mphi_rtm, jmax_rlm, WK_l_tst%nth_sym,                &
     &    idx_trns%lstack_rlm, idx_trns%lstack_even_rlm,                &
     &    leg%P_rtm, leg%dPdt_rtm, WK_l_tst%Ps_tj, WK_l_tst%dPsdt_tj)
!
      if(ntensor .le. 0) return
!
      allocate( WK_l_tst%P2s_tj(WK_l_tst%nth_sym,jmax_rlm) )
      allocate( WK_l_tst%dP2sdt_tj(WK_l_tst%nth_sym,jmax_rlm) )
      call symmetric_legendre_4_tensor_lj                               &
     &  (nth_rtm, mphi_rtm, jmax_rlm, WK_l_tst%nth_sym,                 &
     &   idx_trns%lstack_rlm, idx_trns%lstack_even_rlm,                 &
     &   idx_gl_1d_rlm_j, leg%g_sph_rlm, leg%asin_t_rtm, leg%cos_t_rtm, &
     &   WK_l_tst%Ps_tj, WK_l_tst%dPsdt_tj,                             &
     &   WK_l_tst%P2s_tj, WK_l_tst%dP2sdt_tj)
!
      end subroutine const_legendre_testloop
!
! -----------------------------------------------------------------------
!
      subroutine alloc_leg_vec_test (nri_rtm, maxidx_rtm_t_smp,         &
     &          ntensor, nvector, nscalar, idx_trns, WK_l_tst)
!
      integer(kind = kint), intent(in) :: nri_rtm, maxidx_rtm_t_smp
      integer(kind = kint), intent(in) :: ntensor, nvector, nscalar
      type(index_4_sph_trans), intent(in) :: idx_trns
!
      type(leg_trns_testloop_work), intent(inout) :: WK_l_tst
!
      integer(kind = kint) :: n_hemi
!
!
      n_hemi = (idx_trns%maxdegree_rlm+1)/2
      WK_l_tst%ntsr_jk = n_hemi * nri_rtm * ntensor
      WK_l_tst%nvec_jk = n_hemi * nri_rtm * nvector
      WK_l_tst%nscl_jk = n_hemi * nri_rtm * nscalar
!
      WK_l_tst%n_pol_e =    WK_l_tst%nscl_jk                            &
     &                  + 3*WK_l_tst%nvec_jk + 6*WK_l_tst%ntsr_jk
      WK_l_tst%n_tor_e =  2*WK_l_tst%nvec_jk + 4*WK_l_tst%ntsr_jk
      allocate(WK_l_tst%pol_e(WK_l_tst%n_pol_e,np_smp))
      allocate(WK_l_tst%tor_e(WK_l_tst%n_tor_e,np_smp))
      allocate(WK_l_tst%pol_o(WK_l_tst%n_pol_e,np_smp))
      allocate(WK_l_tst%tor_o(WK_l_tst%n_tor_e,np_smp))
!
      n_hemi = (maxidx_rtm_t_smp+1)/2
      WK_l_tst%ntsr_lk = n_hemi * nri_rtm * ntensor
      WK_l_tst%nvec_lk = n_hemi * nri_rtm * nvector
      WK_l_tst%nscl_lk = n_hemi * nri_rtm * nscalar
!
      WK_l_tst%n_sym_r =    WK_l_tst%nscl_lk                            &
     &                  + 3*WK_l_tst%nvec_lk + 6*WK_l_tst%ntsr_lk 
      WK_l_tst%n_sym_p =  2*WK_l_tst%nvec_lk + 4*WK_l_tst%ntsr_lk
      allocate(WK_l_tst%symp_r(WK_l_tst%n_sym_r,np_smp))
      allocate(WK_l_tst%symp_p(WK_l_tst%n_sym_p,np_smp))
      allocate(WK_l_tst%asmp_r(WK_l_tst%n_sym_r,np_smp))
      allocate(WK_l_tst%asmp_p(WK_l_tst%n_sym_p,np_smp))
!
      end subroutine alloc_leg_vec_test
!
! -----------------------------------------------------------------------
!
      subroutine dealloc_leg_vec_test(WK_l_tst)
!
      type(leg_trns_testloop_work), intent(inout) :: WK_l_tst
!
!
      deallocate(WK_l_tst%pol_e, WK_l_tst%tor_e)
      deallocate(WK_l_tst%pol_o, WK_l_tst%tor_o)
      deallocate(WK_l_tst%symp_r, WK_l_tst%symp_p)
      deallocate(WK_l_tst%asmp_r, WK_l_tst%asmp_p)
!
      deallocate(WK_l_tst%Ps_tj, WK_l_tst%dPsdt_tj)
      if(allocated(WK_l_tst%P2s_tj)) then
        deallocate(WK_l_tst%P2s_tj, WK_l_tst%dP2sdt_tj)
      end if
!
      end subroutine dealloc_leg_vec_test
!
! -----------------------------------------------------------------------
!
      end module t_legendre_work_testlooop
