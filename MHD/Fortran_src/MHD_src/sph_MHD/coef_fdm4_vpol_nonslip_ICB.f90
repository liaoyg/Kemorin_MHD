!>@file   coef_fdm4_vpol_nonslip_ICB.f90
!!@brief  module coef_fdm4_vpol_nonslip_ICB
!!
!!@author H. Matsui
!!@date Programmed in May., 2013
!
!>@brief Matrix to evaluate radial derivative for non-slip at ICB
!!
!!@verbatim
!!      subroutine cal_fdm4_ICB0_nonslip_vp(r_from_ICB, fdm4_noslip_ICB)
!!      subroutine cal_fdm4_ICB1_nonslip_vp(r_from_ICB, fdm4_noslip_ICB)
!!        type(fdm4_ICB_vpol), intent(inout) :: fdm4_noslip_ICB
!!
!!      subroutine check_4th_ICB_nonslip_vp_fdm(fdm4_noslip_ICB)
!!        type(fdm4_ICB_vpol), intent(in) :: fdm4_noslip_ICB
!!
!!   Matrix for poloidal velocity with non-slip boundary at ICB
!!      d2fdr2 =    fdm4_noslip_ICB%dmat_vp0( 2,3) * d_rj(ICB-2)
!!                + fdm4_noslip_ICB%dmat_vp0( 1,3) * d_rj(ICB-1)
!!                + fdm4_noslip_ICB%dmat_vp0( 0,3) * d_rj(ICB  )
!!      d3fdr3 =    fdm4_noslip_ICB%dmat_vp0( 2,4) * d_rj(ICB-2)
!!                + fdm4_noslip_ICB%dmat_vp0( 1,4) * d_rj(ICB-1)
!!                + fdm4_noslip_ICB%dmat_vp0( 0,4) * d_rj(ICB  )
!!
!!   Matrix for poloidal velocity with non-slip boundary at next of ICB
!!      dfdr =      fdm4_noslip_ICB%dmat_vp1( 2,2) * d_rj(ICB-3)
!!                + fdm4_noslip_ICB%dmat_vp1( 1,2) * d_rj(ICB-2)
!!                + fdm4_noslip_ICB%dmat_vp1( 0,2) * d_rj(ICB-1)
!!                + fdm4_noslip_ICB%dmat_vp1(-1,2) * d_rj(ICB  )
!!      d2fdr2 =    fdm4_noslip_ICB%dmat_vp1( 2,3) * d_rj(ICB-3)
!!                + fdm4_noslip_ICB%dmat_vp1( 1,3) * d_rj(ICB-2)
!!                + fdm4_noslip_ICB%dmat_vp1( 0,3) * d_rj(ICB-1)
!!                + fdm4_noslip_ICB%dmat_vp1(-1,3) * d_rj(ICB  )
!!      d3fdr3 =    fdm4_noslip_ICB%dmat_vp1( 2,4) * d_rj(ICB-3)
!!                + fdm4_noslip_ICB%dmat_vp1( 1,4) * d_rj(ICB-2)
!!                + fdm4_noslip_ICB%dmat_vp1( 0,4) * d_rj(ICB-1)
!!                + fdm4_noslip_ICB%dmat_vp1(-1,4) * d_rj(ICB  )
!!      d4fdr4 =    fdm4_noslip_ICB%dmat_vp1( 2,5) * d_rj(ICB-3)
!!                + fdm4_noslip_ICB%dmat_vp1( 1,5) * d_rj(ICB-2)
!!                + fdm4_noslip_ICB%dmat_vp1( 0,5) * d_rj(ICB-1)
!!                + fdm4_noslip_ICB%dmat_vp1(-1,5) * d_rj(ICB  )
!!@endverbatim
!!
!!@n @param r_from_ICB(0:3) radius to three next points of ICB
!!
      module coef_fdm4_vpol_nonslip_ICB
!
      use m_precision
      use m_constants
!
      use t_coef_fdm4_MHD_boundaries
      use cal_inverse_small_matrix
!
      implicit none
!
!
!>      Work matrix to evaluate fdm4_noslip_ICB%dmat_vp0
!!@verbatim
!!      d2fdr2 =    mat_fdm_noslip_ICB_4(3,4) * d_rj(ICB+2)
!!                + mat_fdm_noslip_ICB_4(3,3) * d_rj(ICB+1)
!!                + mat_fdm_noslip_ICB_4(3,1) * d_rj(ICB  )
!!                + mat_fdm_noslip_ICB_4(3,2) * dfdr(ICB)
!!      d3fdr3 =    mat_fdm_noslip_ICB_4(4,4) * d_rj(ICB+2)
!!                + mat_fdm_noslip_ICB_4(4,3) * d_rj(ICB+1)
!!                + mat_fdm_noslip_ICB_4(4,1) * d_rj(ICB  )
!!                + mat_fdm_noslip_ICB_4(4,2) * dfdr(ICB)
!!@endverbatim
      real(kind = kreal) :: mat_fdm_noslip_ICB_4(4,4)
!
!>      Work matrix to evaluate fdm4_noslip_ICB%dmat_vp1
!!@verbatim
!!      dfdr =      mat_fdm_noslip_ICB1_4(2,5) * d_rj(ICB+3)
!!                + mat_fdm_noslip_ICB1_4(2,4) * d_rj(ICB+2)
!!                + mat_fdm_noslip_ICB1_4(2,1) * d_rj(ICB+1)
!!                + mat_fdm_noslip_ICB1_4(2,3) * d_rj(ICB  )
!!                + mat_fdm_noslip_ICB1_4(2,2) * dfdr(ICB)
!!      d2fdr2 =    mat_fdm_noslip_ICB1_4(3,5) * d_rj(ICB+3)
!!                + mat_fdm_noslip_ICB1_4(3,4) * d_rj(ICB+2)
!!                + mat_fdm_noslip_ICB1_4(3,1) * d_rj(ICB+1)
!!                + mat_fdm_noslip_ICB1_4(3,3) * d_rj(ICB  )
!!                + mat_fdm_noslip_ICB1_4(3,2) * dfdr(ICB)
!!      d3fdr3 =    mat_fdm_noslip_ICB1_4(4,5) * d_rj(ICB+3)
!!                + mat_fdm_noslip_ICB1_4(4,4) * d_rj(ICB+2)
!!                + mat_fdm_noslip_ICB1_4(4,1) * d_rj(ICB+1)
!!                + mat_fdm_noslip_ICB1_4(4,3) * d_rj(ICB  )
!!                + mat_fdm_noslip_ICB1_4(4,2) * dfdr(ICB)
!!      d4fdr4 =    mat_fdm_noslip_ICB1_4(5,5) * d_rj(ICB+3)
!!                + mat_fdm_noslip_ICB1_4(5,4) * d_rj(ICB+2)
!!                + mat_fdm_noslip_ICB1_4(5,1) * d_rj(ICB+1)
!!                + mat_fdm_noslip_ICB1_4(5,3) * d_rj(ICB  )
!!                + mat_fdm_noslip_ICB1_4(5,2) * dfdr(ICB)
!!@endverbatim
      real(kind = kreal) :: mat_fdm_noslip_ICB1_4(5,5)
!
      private :: mat_fdm_noslip_ICB_4, mat_fdm_noslip_ICB1_4
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine cal_fdm4_ICB0_nonslip_vp(r_from_ICB, fdm4_noslip_ICB)
!
      real(kind = kreal), intent(in) :: r_from_ICB(0:3)
      type(fdm4_ICB_vpol), intent(inout) :: fdm4_noslip_ICB
!
      integer(kind = kint) :: ierr
      real(kind = kreal) :: mat_taylor_4(4,4)
      real(kind = kreal) :: dr_p1, dr_p2, dr_p3
!
!
      dr_p1 = r_from_ICB(1) - r_from_ICB(0)
      dr_p2 = r_from_ICB(2) - r_from_ICB(0)
      dr_p3 = r_from_ICB(3) - r_from_ICB(0)
!
      mat_taylor_4(1,1) = one
      mat_taylor_4(1,2) = zero
      mat_taylor_4(1,3) = zero
      mat_taylor_4(1,4) = zero
!
      mat_taylor_4(2,1) = zero
      mat_taylor_4(2,2) =  one
      mat_taylor_4(2,3) = zero
      mat_taylor_4(2,4) = zero
!
      mat_taylor_4(3,1) = one
      mat_taylor_4(3,2) = dr_p1
      mat_taylor_4(3,3) = dr_p1*dr_p1 / two
      mat_taylor_4(3,4) = dr_p1**3 / six
!
      mat_taylor_4(4,1) = one
      mat_taylor_4(4,2) = dr_p2
      mat_taylor_4(4,3) = dr_p2*dr_p2 / two
      mat_taylor_4(4,4) = dr_p2**3 / six
!
      call cal_inverse_44_matrix(mat_taylor_4,                          &
     &    mat_fdm_noslip_ICB_4, ierr)
!
      if(ierr .eq. 1) then
        write(*,*) 'singular matrix cal_fdm4_ICB0_nonslip_vp ',         &
     &            r_from_ICB(0)
      end if
!
      fdm4_noslip_ICB%dmat_vp0(0:2,2) = zero
!
      fdm4_noslip_ICB%dmat_vp0( 2,3:4) = mat_fdm_noslip_ICB_4(3:4,4)
      fdm4_noslip_ICB%dmat_vp0( 1,3:4) = mat_fdm_noslip_ICB_4(3:4,3)
      fdm4_noslip_ICB%dmat_vp0( 0,3:4) = mat_fdm_noslip_ICB_4(3:4,1)
!
      end subroutine cal_fdm4_ICB0_nonslip_vp
!
! -----------------------------------------------------------------------
!
      subroutine cal_fdm4_ICB1_nonslip_vp(r_from_ICB, fdm4_noslip_ICB)
!
      real(kind = kreal), intent(in) :: r_from_ICB(0:3)
      type(fdm4_ICB_vpol), intent(inout) :: fdm4_noslip_ICB
!
      integer(kind = kint) :: ierr
      real(kind = kreal) :: mat_taylor_5(5,5)
      real(kind = kreal) :: dr_n1, dr_p1, dr_p2
!
!
      dr_n1 = r_from_ICB(1) - r_from_ICB(0)
      dr_p1 = r_from_ICB(2) - r_from_ICB(1)
      dr_p2 = r_from_ICB(3) - r_from_ICB(1)
!
      mat_taylor_5(1,1) = one
      mat_taylor_5(1,2) = zero
      mat_taylor_5(1,3) = zero
      mat_taylor_5(1,4) = zero
      mat_taylor_5(1,5) = zero
!
      mat_taylor_5(2,1) = zero
      mat_taylor_5(2,2) =  one
      mat_taylor_5(2,3) =-dr_n1
      mat_taylor_5(2,4) = dr_n1*dr_n1 / two
      mat_taylor_5(2,5) =-dr_n1**3 / six
!
      mat_taylor_5(3,1) = one
      mat_taylor_5(3,2) =-dr_n1
      mat_taylor_5(3,3) = dr_n1*dr_n1 / two
      mat_taylor_5(3,4) =-dr_n1**3 / six
      mat_taylor_5(3,5) = dr_n1**4 / (six*four)
!
      mat_taylor_5(4,1) = one
      mat_taylor_5(4,2) = dr_p1
      mat_taylor_5(4,3) = dr_p1*dr_p1 / two
      mat_taylor_5(4,4) = dr_p1**3 / six
      mat_taylor_5(4,5) = dr_p1**4 / (six*four)
!
      mat_taylor_5(5,1) = one
      mat_taylor_5(5,2) = dr_p2
      mat_taylor_5(5,3) = dr_p2*dr_p2 / two
      mat_taylor_5(5,4) = dr_p2**3 / six
      mat_taylor_5(5,5) = dr_p2**4 / (six*four)
!
      call cal_inverse_nn_matrix(ifive, mat_taylor_5,                   &
     &    mat_fdm_noslip_ICB1_4, ierr)
!
      if(ierr .eq. 1) then
        write(*,*) 'singular matrix mat_fdm_noslip_ICB1_4 ',            &
     &            r_from_ICB(0)
      end if
!
      fdm4_noslip_ICB%dmat_vp1(-1,2:5) = mat_fdm_noslip_ICB1_4(2:5,3)
      fdm4_noslip_ICB%dmat_vp1( 0,2:5) = mat_fdm_noslip_ICB1_4(2:5,1)
      fdm4_noslip_ICB%dmat_vp1( 1,2:5) = mat_fdm_noslip_ICB1_4(2:5,4)
      fdm4_noslip_ICB%dmat_vp1( 2,2:5) = mat_fdm_noslip_ICB1_4(2:5,5)
!
      end subroutine cal_fdm4_ICB1_nonslip_vp
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine check_4th_ICB_nonslip_vp_fdm(fdm4_noslip_ICB)
!
      type(fdm4_ICB_vpol), intent(in) :: fdm4_noslip_ICB
!
!
      write(50,*) ' non slip boundary'
      call check_4th_ICB_vpol_fdm(fdm4_noslip_ICB)
!
      end subroutine check_4th_ICB_nonslip_vp_fdm
!
! -----------------------------------------------------------------------
!
      end module coef_fdm4_vpol_nonslip_ICB
