!
!      module set_aiccg_boundaries
!
!     Written by H. Matsui on June, 2005
!
!      subroutine set_bc_4_velo_mat(nst, ned, iele, k1, k2)
!      subroutine set_bc_4_press_mat(iele, k1, k2)
!      subroutine set_bc_4_temp_mat(iele, k1, k2)
!      subroutine set_bc_4_magne_mat(nst, ned, iele, k1, k2)
!      subroutine set_bc_4_m_potential_mat(iele, k1, k2)
!      subroutine set_bc_4_m_potential_mat_cd(iele, k1, k2)
!      subroutine set_bc_4_m_potential_mat_ins(iele, k1, k2)
!      subroutine set_bc_4_d_scalar_mat(iele, k1, k2)
!
      module set_aiccg_boundaries
!
      use m_precision
!
      use m_geometry_parameter
      use m_geometry_data
!
      use correct_matrix_4_boundary
      use set_off_diagonals
!
      implicit none
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine set_bc_4_velo_mat(nst, ned, iele, k1, k2)
!
      use m_solver_djds_fluid
      use m_velo_matrix
!
      integer (kind = kint), intent(in) :: nst, ned, iele, k1, k2
      integer (kind = kint) :: nd2, nod1, nod2, mat_num
!
      nod1 = ie(iele,k1)
      nod2 = ie(iele,k2)
      call set_off_diag_fluid( nod1, nod2, mat_num )
!
      do nd2 = nst, ned
        call correct_matrix33_4_boundary(numnod, itotal_fl_u,           &
     &      itotal_fl_l, nd2, mat_num, Vmat_DJDS%aiccg)
      end do
!
      end subroutine set_bc_4_velo_mat
!
! -----------------------------------------------------------------------
!
      subroutine set_bc_4_temp_mat(iele, k1, k2)
!
      use m_solver_djds_fluid
      use m_temp_matrix
!
      integer (kind = kint), intent(in) :: iele, k1, k2
      integer (kind = kint) :: nod1, nod2, mat_num
!
      nod1 = ie(iele,k1)
      nod2 = ie(iele,k2)
!
      call set_off_diag_fluid ( nod1, nod2, mat_num )
      call correct_matrix11_4_boundary                                  &
     &    (numnod, itotal_fl_u, itotal_fl_l, mat_num, Tmat_DJDS%aiccg)
!
      end subroutine set_bc_4_temp_mat
!
! -----------------------------------------------------------------------
!
      subroutine set_bc_4_press_mat(iele, k1, k2)
!
      use m_solver_djds_linear_fl
      use m_velo_matrix
!
      integer (kind = kint), intent(in) :: iele, k1, k2
      integer (kind = kint) :: nod1, nod2, mat_num
!
!
      nod1 = ie(iele,k1)
      nod2 = ie(iele,k2)
!
      call set_off_diag_linear_fl( nod1, nod2, mat_num )
      call correct_matrix11_4_boundary                                  &
     &   (numnod, itotal1_fl_u, itotal1_fl_l, mat_num, Pmat_DJDS%aiccg)
!
      end subroutine set_bc_4_press_mat
!
! -----------------------------------------------------------------------
!
      subroutine set_bc_4_magne_mat(nst, ned, iele, k1, k2)
!
      use m_solver_djds
      use m_magne_matrix
!
      integer (kind = kint), intent(in) :: nst, ned, iele, k1, k2
      integer (kind = kint) :: nd2, nod1, nod2, mat_num
!
!
      nod1 = ie(iele,k1)
      nod2 = ie(iele,k2)
      call set_off_diag ( nod1, nod2, mat_num )
!
      do nd2 = nst, ned
        call correct_matrix33_4_boundary                                &
     &       (numnod, itotal_u, itotal_l, nd2, mat_num, aiccg_magne )
       end do
!
      end subroutine set_bc_4_magne_mat
!
! -----------------------------------------------------------------------
!
      subroutine set_bc_4_m_potential_mat(iele, k1, k2)
!
      use m_solver_djds_linear
!      use m_solver_djds_linear_cd
!      use m_solver_djds_linear_ins
      use m_mag_potential_matrix
!
      integer (kind = kint), intent(in) :: iele, k1, k2
      integer (kind = kint) :: nod1, nod2, mat_num
!
!
      nod1 = ie(iele,k1)
      nod2 = ie(iele,k2)
!
      call set_off_diag_linear (nod1, nod2, mat_num)
      call correct_matrix11_4_boundary                                  &
     &    (numnod, itotal1_u, itotal1_l, mat_num, aiccg_mag_p)
!
!      call set_off_diag_linear_ins ( nod1, nod2, mat_num )
!      call correct_matrix11_4_boundary                                 &
!     &    (numnod, itotal1_ins_u, itotal1_ins_l, mat_num, aiccg_mag_pi)
!
!      call set_off_diag_linear_cd ( nod1, nod2, mat_num )
!      call correct_matrix11_4_boundary                                 &
!     &    (numnod, itotal1_cd_u, itotal1_cd_l, mat_num, aiccg_mag_pc)
!
      end subroutine set_bc_4_m_potential_mat
!
! -----------------------------------------------------------------------
!
!      subroutine set_bc_4_m_potential_mat_cd(iele, k1, k2)
!
!      use m_solver_djds_linear_cd
!      use m_mag_potential_matrix
!
!      integer (kind = kint), intent(in) :: iele, k1, k2
!      integer (kind = kint) :: nod1, nod2, mat_num
!
!
!      nod1 = ie(iele,k1)
!      nod2 = ie(iele,k2)
!
!      call set_off_diag_linear_cd ( nod1, nod2, mat_num )
!      call correct_matrix11_4_boundary                                 &
!     &    (numnod, itotal1_cd_u, itotal1_cd_l, mat_num, aiccg_mag_pc)
!
!      end subroutine set_bc_4_m_potential_mat_cd
!
! -----------------------------------------------------------------------
!
!      subroutine set_bc_4_m_potential_mat_ins(iele, k1, k2)
!
!      use m_solver_djds_linear_ins
!      use m_mag_potential_matrix
!
!      integer (kind = kint), intent(in) :: iele, k1, k2
!      integer (kind = kint) :: nod1, nod2, mat_num
!
!
!      nod1 = ie(iele,k1)
!      nod2 = ie(iele,k2)
!
!      call set_off_diag_linear_ins ( nod1, nod2, mat_num )
!      call correct_matrix11_4_boundary                                 &
!     &    (numnod, itotal1_ins_u, itotal1_ins_l, mat_num, aiccg_mag_pi)
!
!      end subroutine set_bc_4_m_potential_mat_ins
!
! -----------------------------------------------------------------------
!
      subroutine set_bc_4_d_scalar_mat(iele, k1, k2)
!
      use m_solver_djds_fluid
      use m_light_element_matrix
!
      integer (kind = kint), intent(in) :: iele, k1, k2
      integer (kind = kint) :: nod1, nod2, mat_num
!
!
      nod1 = ie(iele,k1)
      nod2 = ie(iele,k2)
!
      call set_off_diag_fluid ( nod1, nod2, mat_num )
      call correct_matrix11_4_boundary                                  &
     &    (numnod, itotal_fl_u, itotal_fl_l, mat_num, Cmat_DJDS%aiccg)
!
      end subroutine set_bc_4_d_scalar_mat
!
! -----------------------------------------------------------------------
!
      end module set_aiccg_boundaries
