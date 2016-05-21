!>@file   t_sph_trans_arrays_MHD.f90
!!@brief  module t_sph_trans_arrays_MHD
!!
!!@author H. Matsui
!!@date Programmed in Jan., 2010
!
!!@verbatim
!!      subroutine alloc_sph_trans_address(sph_rtp, WK)
!!      subroutine dealloc_sph_trans_address(WK)
!!        type(sph_rtp_grid), intent(in) :: sph_rtp
!!        type(works_4_sph_trans_MHD), intent(inout) :: WK
!!@endverbatim
!
!>@brief Field addresses for spherical harmonics transform
!!       in MHD dynamo simulation
!!
      module t_sph_trans_arrays_MHD
!
      use m_precision
!
      use t_phys_address
      use t_addresses_sph_transform
!
      implicit none
!
!>        strucutres for spherical transform dor MHD dynamo
      type works_4_sph_trans_MHD
!>        strucutre for spherical transform data addresses
        type(address_4_sph_trans) :: trns_MHD
!>      strucutre for spherical transform data addresses
        type(address_4_sph_trans) :: trns_snap
!>      strucutre for spherical transform data addresses
        type(address_4_sph_trans) :: trns_tmp
!
!>      Nonoliear terms data using simulation
        real(kind = kreal), pointer :: frm_rtp(:,:)
!
!>      field data to evaluate nonliear terms at pole
        real(kind = kreal), pointer :: fls_pl(:,:)
!>        local field data to evaluate nonliear terms at pole
        real(kind = kreal), pointer :: flc_pl(:,:)
!
!>        field data to evaluate nonliear terms at pole
        real(kind = kreal), pointer :: frs_pl(:,:)
!>        field data to evaluate nonliear terms at pole
        real(kind = kreal), pointer :: frm_pl(:,:)
      end type works_4_sph_trans_MHD
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine alloc_sph_trans_address(sph_rtp, WK)
!
      use t_spheric_rtp_data
      use m_work_pole_sph_trans
!
      type(sph_rtp_grid), intent(in) :: sph_rtp
      type(works_4_sph_trans_MHD), intent(inout) :: WK
!
!
      call alloc_nonlinear_data(sph_rtp%nnod_rtp, wk%trns_MHD)
      call alloc_nonlinear_data(sph_rtp%nnod_rtp, WK%trns_snap)
      call alloc_nonlinear_data(sph_rtp%nnod_rtp, wk%trns_tmp)
!
      allocate(WK%frm_rtp(sph_rtp%nnod_rtp,WK%trns_MHD%ncomp_rtp_2_rj))
      if(WK%trns_MHD%ncomp_rtp_2_rj .gt. 0) WK%frm_rtp = 0.0d0
!
!
      allocate(WK%fls_pl(nnod_pole,WK%trns_snap%ncomp_rj_2_rtp))
      allocate(WK%flc_pl(nnod_pole,WK%trns_snap%ncomp_rj_2_rtp))
      allocate(WK%frs_pl(nnod_pole,WK%trns_snap%ncomp_rtp_2_rj))
      allocate(WK%frm_pl(nnod_pole,WK%trns_MHD%ncomp_rtp_2_rj))
!
      if(WK%trns_snap%ncomp_rj_2_rtp .gt. 0) WK%fls_pl = 0.0d0
      if(WK%trns_snap%ncomp_rj_2_rtp .gt. 0) WK%flc_pl = 0.0d0
      if(WK%trns_snap%ncomp_rtp_2_rj .gt. 0) WK%frs_pl = 0.0d0
      if(WK%trns_MHD%ncomp_rtp_2_rj .gt. 0)  WK%frm_pl = 0.0d0
!
      end subroutine alloc_sph_trans_address
!
!-----------------------------------------------------------------------
!
      subroutine dealloc_sph_trans_address(WK)
!
      type(works_4_sph_trans_MHD), intent(inout) :: WK
!
      deallocate(WK%frm_rtp, WK%fls_pl, WK%flc_pl)
      deallocate(WK%frs_pl, WK%frm_pl)
!
      call dealloc_nonlinear_data(WK%trns_tmp)
      call dealloc_nonlinear_data(WK%trns_snap)
      call dealloc_nonlinear_data(WK%trns_MHD)
!
      end subroutine dealloc_sph_trans_address
!
!-----------------------------------------------------------------------
!
      end module t_sph_trans_arrays_MHD
