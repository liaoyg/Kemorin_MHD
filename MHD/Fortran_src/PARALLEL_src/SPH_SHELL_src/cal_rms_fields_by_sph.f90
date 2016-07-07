!>@file   cal_rms_fields_by_sph.f90
!!@brief      module cal_rms_fields_by_sph
!!
!!@author H. Matsui and H. Okuda
!!@date Programmed in  Dec., 2012
!
!> @brief evaluate mean square data from spectr data
!!
!!@verbatim
!!      subroutine init_rms_4_sph_spectr(l_truncation, sph_rj, rj_fld)
!!
!!      subroutine cal_mean_squre_in_shell(kr_st, kr_ed, l_truncation,  &
!!     &          sph_rj, ipol, rj_fld, g_sph_rj)
!!        type(sph_rj_grid), intent(in) :: sph_rj
!!        type(phys_data), intent(in) :: rj_fld
!!@endverbatim
!
      module cal_rms_fields_by_sph
!
      use m_precision
      use m_constants
      use m_machine_parameter
!
      use t_spheric_rj_data
      use t_phys_data
      use t_phys_address
      use t_sum_sph_rms_data
!
      implicit none
!
      type(sph_mean_square_work), save, private :: WK_pwr
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine init_rms_4_sph_spectr(l_truncation, sph_rj, rj_fld)
!
      use calypso_mpi
      use m_rms_4_sph_spectr
!
      use sum_sph_rms_data
      use volume_average_4_sph
      use quicksort
!
      integer(kind = kint), intent(in) :: l_truncation
      type(sph_rj_grid), intent(in) :: sph_rj
      type(phys_data), intent(in) :: rj_fld
!
      integer(kind = kint) :: i_fld, j_fld
      integer(kind = kint) :: k, knum, num_field
!
!
      num_field = 0
      do i_fld = 1, rj_fld%num_phys
        num_field = num_field + rj_fld%iflag_monitor(i_fld)
      end do
!
      call alloc_rms_name_sph_spec(num_field, pwr1)
!
      j_fld = 0
      do i_fld = 1, rj_fld%num_phys
        if(rj_fld%iflag_monitor(i_fld) .gt. 0) then
          j_fld = j_fld + 1
          pwr1%id_field(j_fld) =   i_fld
          pwr1%num_comp_sq(j_fld) =    rj_fld%num_component(i_fld)
          pwr1%istack_comp_sq(j_fld) = pwr1%istack_comp_sq(j_fld-1)     &
     &                              + rj_fld%num_component(i_fld)
          pwr1%pwr_name(j_fld) =   rj_fld%phys_name(i_fld)
        end if
      end do
!
      call quicksort_int                                                &
     &   (pwr1%nri_rms, pwr1%kr_4_rms, ione, pwr1%nri_rms)
!
      call alloc_rms_4_sph_spectr(my_rank, l_truncation, pwr1)
      call alloc_ave_4_sph_spectr                                       &
     &   (sph_rj%idx_rj_degree_zero, sph_rj%nidx_rj(1), pwr1)
      call allocate_rms_sph_local_data                                  &
     &   (l_truncation, sph_rj%nidx_rj, WK_pwr)
!
      call set_sum_table_4_sph_spectr                                   &
     &   (l_truncation, sph_rj%nidx_rj, sph_rj%idx_gl_1d_rj_j,          &
     &    WK_pwr%num_mode_sum_l,  WK_pwr%num_mode_sum_m,                &
     &    WK_pwr%num_mode_sum_lm, WK_pwr%istack_mode_sum_l,             &
     &    WK_pwr%istack_mode_sum_m, WK_pwr%istack_mode_sum_lm,          &
     &    WK_pwr%item_mode_sum_l, WK_pwr%item_mode_sum_m,               &
     &    WK_pwr%item_mode_sum_lm)
!
!
      do knum = 1, pwr1%nri_rms
        k = pwr1%kr_4_rms(knum)
        if(k .le. 0) then
          pwr1%r_4_rms(knum) = 0.0d0
        else
          pwr1%r_4_rms(knum) = sph_rj%radius_1d_rj_r(k)
        end if
      end do
!
      end subroutine init_rms_4_sph_spectr
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine cal_mean_squre_in_shell(kr_st, kr_ed, l_truncation,    &
     &          sph_rj, ipol, rj_fld, g_sph_rj)
!
      use calypso_mpi
      use m_rms_4_sph_spectr
!
      use volume_average_4_sph
      use cal_ave_4_rms_vector_sph
      use sum_sph_rms_data
!
      type(sph_rj_grid), intent(in) ::  sph_rj
      type(phys_data), intent(in) :: rj_fld
      type(phys_address), intent(in) :: ipol
      integer(kind = kint), intent(in) :: l_truncation
      integer(kind = kint), intent(in) :: kr_st, kr_ed
      real(kind = kreal), intent(in) :: g_sph_rj(sph_rj%nidx_rj(2),13)
!
      real(kind = kreal) :: avol
!
!
      if(pwr1%num_comp_rms .eq. 0) return

      if(iflag_debug .gt. 0) write(*,*) 'cal_one_over_volume'
      call cal_one_over_volume(kr_st, kr_ed,                            &
     &   sph_rj%nidx_rj(1), sph_rj%radius_1d_rj_r, avol)
      if(iflag_debug .gt. 0) write(*,*) 'sum_sph_layerd_rms'
      call sum_sph_layerd_rms(kr_st, kr_ed, l_truncation,               &
     &    sph_rj, ipol, g_sph_rj, rj_fld,                               &
     &    pwr1%nri_rms, pwr1%num_fld_rms, pwr1%num_comp_rms,            &
     &    pwr1%istack_comp_sq, pwr1%id_field,                           &
     &    WK_pwr%istack_mode_sum_l,  WK_pwr%istack_mode_sum_m,          &
     &    WK_pwr%istack_mode_sum_lm, WK_pwr%item_mode_sum_l,            &
     &    WK_pwr%item_mode_sum_m,    WK_pwr%item_mode_sum_lm,           &
     &    pwr1%kr_4_rms, WK_pwr%shl_rj, WK_pwr%volume_j,                &
     &    WK_pwr%shl_l_local, WK_pwr%shl_m_local, WK_pwr%shl_lm_local,  &
     &    WK_pwr%vol_l_local, WK_pwr%vol_m_local, WK_pwr%vol_lm_local)
!
      call global_sum_sph_layerd_rms                                    &
     &    (l_truncation, pwr1%nri_rms, pwr1%num_comp_rms,               &
     &     WK_pwr%shl_l_local, WK_pwr%shl_m_local, WK_pwr%shl_lm_local, &
     &     WK_pwr%vol_l_local, WK_pwr%vol_m_local, WK_pwr%vol_lm_local, &
     &     pwr1%shl_l, pwr1%shl_m, pwr1%shl_lm,                         &
     &     pwr1%vol_l, pwr1%vol_m, pwr1%vol_lm,                         &
     &     pwr1%shl_sq, pwr1%shl_m0, pwr1%ratio_shl_m0,                 &
     &     pwr1%vol_sq, pwr1%vol_m0, pwr1%ratio_vol_m0)
!
      if(my_rank .eq. 0) then
        if(iflag_debug .gt. 0) write(*,*) 'surf_ave_4_sph_rms_int'
        call surf_ave_4_sph_rms_int                                     &
     &     (l_truncation, sph_rj%nidx_rj(1), sph_rj%a_r_1d_rj_r,        &
     &      pwr1%nri_rms, pwr1%num_comp_rms, pwr1%kr_4_rms,             &
     &      pwr1%shl_l, pwr1%shl_m, pwr1%shl_lm, pwr1%shl_sq, pwr1%shl_m0)
        call vol_ave_4_rms_sph(l_truncation, pwr1%num_comp_rms, avol,   &
     &      pwr1%vol_l, pwr1%vol_m, pwr1%vol_lm,                        &
     &      pwr1%vol_sq, pwr1%vol_m0)
      end if
!
      if(iflag_debug .gt. 0) write(*,*) 'cal_volume_average_sph'
      call cal_volume_average_sph(kr_st, kr_ed, avol, sph_rj, rj_fld)
!
      end subroutine cal_mean_squre_in_shell
!
! ----------------------------------------------------------------------
!
      subroutine global_sum_sph_layerd_rms                              &
     &         (l_truncation, nri_rms, ntot_rms_rj,                     &
     &          rms_sph_l_local, rms_sph_m_local, rms_sph_lm_local,     &
     &          rms_sph_vl_local, rms_sph_vm_local, rms_sph_vlm_local,  &
     &          rms_sph_l, rms_sph_m, rms_sph_lm,                       &
     &          rms_sph_vol_l, rms_sph_vol_m, rms_sph_vol_lm,           &
     &          rms_sph, rms_sph_m0, ratio_sph_m0,                      &
     &          rms_sph_vol, rms_sph_vol_m0, ratio_sph_vol_m0)
!
      use calypso_mpi
!
      use cal_rms_by_sph_spectr
      use cal_ave_4_rms_vector_sph
      use radial_int_for_sph_spec
!
      integer(kind = kint), intent(in) :: l_truncation
      integer(kind = kint), intent(in) :: nri_rms
      integer(kind = kint), intent(in) :: ntot_rms_rj
!
      real(kind = kreal), intent(in)                                    &
     &          :: rms_sph_l_local(nri_rms,0:l_truncation,ntot_rms_rj)
      real(kind = kreal), intent(in)                                    &
     &          :: rms_sph_m_local(nri_rms,0:l_truncation,ntot_rms_rj)
      real(kind = kreal), intent(in)                                    &
     &          :: rms_sph_lm_local(nri_rms,0:l_truncation,ntot_rms_rj)
!
      real(kind = kreal), intent(in)                                    &
     &          :: rms_sph_vl_local(0:l_truncation,ntot_rms_rj)
      real(kind = kreal), intent(in)                                    &
     &          :: rms_sph_vm_local(0:l_truncation,ntot_rms_rj)
      real(kind = kreal), intent(in)                                    &
     &          :: rms_sph_vlm_local(0:l_truncation,ntot_rms_rj)
!
      real(kind = kreal), intent(inout)                                 &
     &          :: rms_sph_l(nri_rms,0:l_truncation,ntot_rms_rj)
      real(kind = kreal), intent(inout)                                 &
     &          :: rms_sph_m(nri_rms,0:l_truncation,ntot_rms_rj)
      real(kind = kreal), intent(inout)                                 &
     &          :: rms_sph_lm(nri_rms,0:l_truncation,ntot_rms_rj)
!
      real(kind = kreal), intent(inout)                                 &
     &          :: rms_sph_vol_l(0:l_truncation,ntot_rms_rj)
      real(kind = kreal), intent(inout)                                 &
     &          :: rms_sph_vol_m(0:l_truncation,ntot_rms_rj)
      real(kind = kreal), intent(inout)                                 &
     &          :: rms_sph_vol_lm(0:l_truncation,ntot_rms_rj)
!
      real(kind = kreal), intent(inout) :: rms_sph(nri_rms,ntot_rms_rj)
      real(kind = kreal), intent(inout)                                 &
     &          :: rms_sph_m0(nri_rms,ntot_rms_rj)
      real(kind = kreal), intent(inout)                                 &
     &          :: ratio_sph_m0(nri_rms,ntot_rms_rj)
!
      real(kind = kreal), intent(inout) :: rms_sph_vol(ntot_rms_rj)
      real(kind = kreal), intent(inout) :: rms_sph_vol_m0(ntot_rms_rj)
      real(kind = kreal), intent(inout) :: ratio_sph_vol_m0(ntot_rms_rj)
!
      integer(kind = kint) :: num
!
!
      num = ntot_rms_rj * (l_truncation + 1)
      call MPI_REDUCE (rms_sph_vl_local, rms_sph_vol_l,                 &
     &    num, CALYPSO_REAL, MPI_SUM, izero, CALYPSO_COMM, ierr_MPI)
      call MPI_REDUCE (rms_sph_vm_local, rms_sph_vol_m,                 &
     &    num, CALYPSO_REAL, MPI_SUM, izero, CALYPSO_COMM, ierr_MPI)
      call MPI_REDUCE (rms_sph_vlm_local, rms_sph_vol_lm,               &
     &    num, CALYPSO_REAL, MPI_SUM, izero, CALYPSO_COMM, ierr_MPI)
!
      if(my_rank .eq. 0) then
        call sum_sph_vol_rms_all_modes(l_truncation, ntot_rms_rj,       &
     &      rms_sph_vol_l, rms_sph_vol)
        call pick_axis_sph_vol_pwr                                      &
     &     (l_truncation, ntot_rms_rj, rms_sph_vol_m, rms_sph_vol,      &
     &      rms_sph_vol_m0, ratio_sph_vol_m0)
      end if
!
      if(nri_rms .le. 0) return
      num = ntot_rms_rj * nri_rms * (l_truncation + 1)
      call MPI_REDUCE (rms_sph_l_local, rms_sph_l, num, CALYPSO_REAL,   &
     &    MPI_SUM, izero, CALYPSO_COMM, ierr_MPI)
      call MPI_REDUCE (rms_sph_m_local, rms_sph_m, num, CALYPSO_REAL,   &
     &    MPI_SUM, izero, CALYPSO_COMM, ierr_MPI)
      call MPI_REDUCE (rms_sph_lm_local, rms_sph_lm, num, CALYPSO_REAL, &
     &    MPI_SUM, izero, CALYPSO_COMM, ierr_MPI)
!
      if(my_rank .gt. 0) return
      call sum_sph_rms_all_modes(l_truncation, nri_rms, ntot_rms_rj,    &
     &    rms_sph_l, rms_sph)
      call pick_axis_sph_power(l_truncation, nri_rms, ntot_rms_rj,      &
     &    rms_sph_m, rms_sph, rms_sph_m0, ratio_sph_m0)
!
      end subroutine global_sum_sph_layerd_rms
!
! -----------------------------------------------------------------------
!
      end module cal_rms_fields_by_sph
