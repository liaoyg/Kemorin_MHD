!>@file   global_field_4_dynamobench.f90
!!@brief  module global_field_4_dynamobench
!!
!!@author H. Matsui
!!@date    programmed by H.Matsui in June., 2012
!
!>@brief Evaluate global data for dynamo benchmark test
!!
!!@verbatim
!!      subroutine copy_energy_4_dynamobench
!!      subroutine copy_icore_energy_4_dbench
!!
!!      subroutine pick_inner_core_rotation(idx_rj_degree_one, nidx_rj, &
!!     &          nlayer_ICB, ar_1d_rj, it_velo,                        &
!!     &          nnod_rj, ntot_phys_rj, d_rj)
!!      subroutine pick_mag_torque_inner_core(idx_rj_degree_one,        &
!!     &          nidx_rj, nlayer_ICB, radius_1d_rj_r, it_lorentz,      &
!!     &          nnod_rj, ntot_phys_rj, d_rj)
!!@endverbatim
!
      module global_field_4_dynamobench
!
      use m_precision
!
      use m_constants
      use m_field_4_dynamobench
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine copy_energy_4_dynamobench
!
      use m_phys_labels
      use m_rms_4_sph_spectr
!
      integer(kind = kint) :: i_fld, i_comp
!
!
      do i_fld = 1, num_rms_rj
        if(rms_name_rj(i_fld) .eq. fhd_velo) then
          i_comp = istack_rms_comp_rj(i_fld-1) + 1
          KE_bench(1) = pwr1%vol_sq(i_comp  )
          KE_bench(2) = pwr1%vol_sq(i_comp+1)
          KE_bench(3) = pwr1%vol_sq(i_comp+2)
          exit
        end if
      end do
!
      do i_fld = 1, num_rms_rj
        if(rms_name_rj(i_fld) .eq. fhd_magne) then
          i_comp = istack_rms_comp_rj(i_fld-1) + 1
          ME_bench(1) = pwr1%vol_sq(i_comp  )
          ME_bench(2) = pwr1%vol_sq(i_comp+1)
          ME_bench(3) = pwr1%vol_sq(i_comp+2)
          exit
        end if
      end do
!
      end subroutine copy_energy_4_dynamobench
!
! ----------------------------------------------------------------------
!
      subroutine copy_icore_energy_4_dbench
!
      use m_phys_labels
      use m_rms_4_sph_spectr
!
      integer(kind = kint) :: i_fld, i_comp
!
!
      do i_fld = 1, num_rms_rj
        if(rms_name_rj(i_fld) .eq. fhd_magne) then
          i_comp = istack_rms_comp_rj(i_fld-1) + 1
          mene_icore(1) = pwr1%vol_sq(i_comp  )
          mene_icore(2) = pwr1%vol_sq(i_comp+1)
          mene_icore(3) = pwr1%vol_sq(i_comp+2)
!
          exit
        end if
      end do
!
      end subroutine copy_icore_energy_4_dbench
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine pick_inner_core_rotation(idx_rj_degree_one, nidx_rj,   &
     &          nlayer_ICB, ar_1d_rj, it_velo,                          &
     &          nnod_rj, ntot_phys_rj, d_rj)
!
      use calypso_mpi
!
      integer(kind = kint), intent(in) :: idx_rj_degree_one(-1:1)
      integer(kind = kint), intent(in) :: nidx_rj(2)
      integer(kind = kint), intent(in) :: nlayer_ICB
      real(kind = kreal), intent(in) :: ar_1d_rj(nidx_rj(1),3)
!
      integer(kind = kint), intent(in) :: it_velo
      integer(kind = kint), intent(in) :: nnod_rj, ntot_phys_rj
      real (kind=kreal), intent(in) :: d_rj(nnod_rj,ntot_phys_rj)
!
      integer(kind = kint) :: i, i10c_o
      real(kind = kreal) :: rotate_ic_local(-1:1)
!
!
      do i = -1, 1
        if(idx_rj_degree_one(i) .gt. 0) then
          i10c_o = idx_rj_degree_one(i) + (nlayer_ICB-1)*nidx_rj(2)
          rotate_ic_local(i) = d_rj(i10c_o,it_velo)                     &
     &                       * ar_1d_rj(nlayer_ICB,2)
        else
          rotate_ic_local(i) = zero
        end if
      end do
!
      call MPI_allREDUCE (rotate_ic_local, rotate_icore, ithree,        &
     &    CALYPSO_REAL, MPI_SUM, CALYPSO_COMM, ierr_MPI)
!
      end subroutine pick_inner_core_rotation
!
! ----------------------------------------------------------------------
!
      subroutine pick_mag_torque_inner_core(idx_rj_degree_one,          &
     &          nidx_rj, nlayer_ICB, radius_1d_rj_r, it_lorentz,        &
     &          nnod_rj, ntot_phys_rj, d_rj)
!
      use calypso_mpi
!
      integer(kind = kint), intent(in) :: idx_rj_degree_one(-1:1)
      integer(kind = kint), intent(in) :: nidx_rj(2)
      integer(kind = kint), intent(in) :: nlayer_ICB
      integer(kind = kint), intent(in) :: it_lorentz
      integer(kind = kint), intent(in) :: nnod_rj, ntot_phys_rj
      real(kind = kreal), intent(in) :: radius_1d_rj_r(nidx_rj(1))
      real (kind=kreal), intent(in) :: d_rj(nnod_rj,ntot_phys_rj)
!
      integer(kind = kint) :: i, i10c_o
      real(kind = kreal) :: m_torque_local(-1:1)
!
!
      do i = -1, 1
        if(idx_rj_degree_one(i) .gt. 0) then
          i10c_o = idx_rj_degree_one(i) + (nlayer_ICB-1)*nidx_rj(2)
          m_torque_local(i) = d_rj(i10c_o,it_lorentz)                   &
     &                       * (radius_1d_rj_r(nlayer_ICB)**3)          &
     &                       * eight * four*atan(one) / (five*three)
        else
          m_torque_local(i) = zero
        end if
      end do
!
      call MPI_allREDUCE (m_torque_local, m_torque_icore, ithree,       &
     &    CALYPSO_REAL, MPI_SUM, CALYPSO_COMM, ierr_MPI)
!
      end subroutine pick_mag_torque_inner_core
!
! ----------------------------------------------------------------------
!
      end module global_field_4_dynamobench
