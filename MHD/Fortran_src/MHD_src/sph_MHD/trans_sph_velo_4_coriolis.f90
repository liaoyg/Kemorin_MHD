!>@file   trans_sph_velo_4_coriolis.f90
!!@brief  module trans_sph_velo_4_coriolis
!!
!!@author H. Matsui
!!@date Programmed in Jan., 2010
!
!>@brief  Data transfer to evaluate Coriolis term
!!
!!@verbatim
!!      subroutine s_trans_sph_velo_4_coriolis(kr_in, kr_out)
!!@endverbatim
!!
!!@n @param kr_in       Radial ID for inner boundary
!!@n @param kr_out      Radial ID for outer boundary
!
      module trans_sph_velo_4_coriolis
!
      use m_precision
!
      use calypso_mpi
      use m_spheric_parameter
      use m_sph_spectr_data
      use m_sph_phys_address
      use m_integrals_4_sph_coriolis
!
      implicit none
!
      private :: copy_sph_diff_velocties
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine s_trans_sph_velo_4_coriolis(kr_in, kr_out)
!
      use solver_sph_coriolis_sr
!
      integer(kind = kint), intent(in) :: kr_in, kr_out
!
!
      call copy_sph_diff_velocties(kr_in, kr_out)
!
      call solver_sph_coriolis_sr_5(nshift_j_cor,                       &
     &     nidx_rj(2), nidx_rj(1), d_rj_cor(1,1),                       &
     &     nidx_j_cor, d_sph_cor(1,1))
!
      end subroutine s_trans_sph_velo_4_coriolis
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine copy_sph_diff_velocties(kr_in, kr_out)
!
      use m_schmidt_poly_on_rtm
!
      integer(kind = kint), intent(in) :: kr_in, kr_out
      integer(kind = kint) :: ist, ied, i, k, j
!
!
      ist = (kr_in-1)*nidx_rj(2) + 1
      ied = kr_out * nidx_rj(2)
!$omp parallel do private(k,j)
      do i = ist, ied
        j = mod((i-1),nidx_rj(2)) + 1
        k =  1 + (i- j) / nidx_rj(2)
!
        d_rj_cor(i,ic_vp) =    d_rj(i,ipol%i_velo )
        d_rj_cor(i,ic_dvp) =   d_rj(i,idpdr%i_velo)
        d_rj_cor(i,ic_d2vp) = -d_rj(i,itor%i_vort )                     &
     &                       + g_sph_rj(j,3)*ar_1d_rj(k,2)              &
     &                        *d_rj(i,ipol%i_velo)
        d_rj_cor(i,ic_vt) =    d_rj(i,ipol%i_vort )
        d_rj_cor(i,ic_dvt) =   d_rj(i,idpdr%i_vort)
      end do
!$omp end parallel do
!
!
      end subroutine copy_sph_diff_velocties
!
! -----------------------------------------------------------------------
!
      end module trans_sph_velo_4_coriolis
