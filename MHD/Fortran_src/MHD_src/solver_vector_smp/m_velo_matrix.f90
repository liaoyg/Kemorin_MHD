!>@file   m_velo_matrix.f90
!!@brief  module m_velo_matrix
!!
!!@author H. Matsui
!!@date Programmed in Dec., 2002
!!@date Modified in Nov., 2013
!
!>     DJDS matrix for velocity
!!
!!@verbatim
!!      subroutine allocate_aiccg_velo
!!      subroutine deallocate_aiccg_velo
!!      subroutine reset_aiccg_velo
!!
!!       subroutine allocate_aiccg_press
!!       subroutine reset_aiccg_press
!!       subroutine deallocate_aiccg_press
!!@endverbatim
!
      module m_velo_matrix
!
      use m_precision
      use t_solver_djds
!
      implicit  none
!
!
!>      Structure of matrix for time evolution of velocity
      type(DJDS_MATRIX), save :: Vmat_DJDS
!>      Structure of matrix for poission equation of pressure
      type(DJDS_MATRIX), save :: Pmat_DJDS
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine allocate_aiccg_velo
!
       use m_geometry_parameter
       use m_solver_djds_fluid
!
!
       Vmat_DJDS%num_diag =      numnod
       Vmat_DJDS%internal_diag = internal_node
!
       Vmat_DJDS%istart_diag = 1
       Vmat_DJDS%istart_l =    9*numnod + 1
       Vmat_DJDS%istart_u =    9*(numnod+itotal_fl_l) + 1
!
       Vmat_DJDS%num_non0 =    9 * (numnod+itotal_fl_u+itotal_fl_l)
!
       allocate(Vmat_DJDS%aiccg(-8:Vmat_DJDS%num_non0) )
       allocate(Vmat_DJDS%ALUG_U(9*internal_node) )
       allocate(Vmat_DJDS%ALUG_L(9*internal_node) )
!
       Vmat_DJDS%D =>  Vmat_DJDS%aiccg(Vmat_DJDS%istart_diag:Vmat_DJDS%istart_l-1)
       Vmat_DJDS%AL => Vmat_DJDS%aiccg(Vmat_DJDS%istart_l:Vmat_DJDS%istart_u-1)
       Vmat_DJDS%AU => Vmat_DJDS%aiccg(Vmat_DJDS%istart_u:Vmat_DJDS%num_non0)
!
       call reset_aiccg_velo
!
       end subroutine allocate_aiccg_velo
!
! ----------------------------------------------------------------------
!
      subroutine allocate_aiccg_press
!
      use m_geometry_parameter
      use m_solver_djds_linear_fl
!
!
      Pmat_DJDS%num_diag =      numnod
      Pmat_DJDS%internal_diag = internal_node
!
      Pmat_DJDS%istart_diag = 1
      Pmat_DJDS%istart_l = numnod + 1
      Pmat_DJDS%istart_u = numnod + itotal1_fl_l + 1
!
      Pmat_DJDS%num_non0 = numnod + itotal1_fl_u + itotal1_fl_l
!
      allocate(Pmat_DJDS%aiccg(0:Pmat_DJDS%num_non0))
!
      allocate (Pmat_DJDS%ALUG_U(internal_node) )
      allocate (Pmat_DJDS%ALUG_L(internal_node) )
!
       Pmat_DJDS%D =>  Pmat_DJDS%aiccg(Pmat_DJDS%istart_diag:Pmat_DJDS%istart_l-1)
       Pmat_DJDS%AL => Pmat_DJDS%aiccg(Pmat_DJDS%istart_l:Pmat_DJDS%istart_u-1)
       Pmat_DJDS%AU => Pmat_DJDS%aiccg(Pmat_DJDS%istart_u:Pmat_DJDS%num_non0)
!
      call reset_aiccg_press
!
      end subroutine allocate_aiccg_press
!
! ----------------------------------------------------------------------
!
      subroutine deallocate_aiccg_velo
!
!
       call dealloc_type_djds_mat(Vmat_DJDS)
!
       end subroutine deallocate_aiccg_velo
!
! ----------------------------------------------------------------------
!
       subroutine deallocate_aiccg_press
!
!
       call dealloc_type_djds_mat(Pmat_DJDS)
!
       end subroutine deallocate_aiccg_press
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine reset_aiccg_velo
!
      use m_geometry_parameter
      use m_geometry_data
      use m_geometry_data_MHD
      use m_solver_djds_fluid
!
      integer (kind = kint) :: inod, iele, k1, in
!
!
      Vmat_DJDS%aiccg = 0.0d00
!
      do inod = 1, numnod
        Vmat_DJDS%aiccg(inod*9-8) = 1.0d0
        Vmat_DJDS%aiccg(inod*9-4) = 1.0d0
        Vmat_DJDS%aiccg(inod*9  ) = 1.0d0
      end do
!
      do k1 = 1, nnod_4_ele
        do iele = iele_fl_start, iele_fl_end
          inod = ie(iele,k1)
          in = OLDtoNEW(inod)
          Vmat_DJDS%aiccg(in*9-8) = 0.0d0
          Vmat_DJDS%aiccg(in*9-4) = 0.0d0
          Vmat_DJDS%aiccg(in*9  ) = 0.0d0
        end do
      end do
!
      Vmat_DJDS%ALUG_U= 1.d0
      Vmat_DJDS%ALUG_L= 1.d0
!
      end subroutine reset_aiccg_velo
!
! ----------------------------------------------------------------------
!
      subroutine reset_aiccg_press
!
      use m_geometry_constants
      use m_geometry_parameter
      use m_geometry_data
      use m_geometry_data_MHD
      use m_solver_djds_linear_fl
!
      integer(kind = kint) :: inod, iele, in, k1
!
!
      Pmat_DJDS%aiccg = 0.0d0
!
      do inod = 1, numnod
        Pmat_DJDS%aiccg(inod) = 1.0d0
      end do
!
      do k1 = 1, num_t_linear
        do iele = iele_fl_start, iele_fl_end
          inod = ie(iele,k1)
          in = OLDtoNEW1(inod)
          Pmat_DJDS%aiccg(in) = 0.0d0
        end do
      end do
!
      Pmat_DJDS%ALUG_U= 1.d0
      Pmat_DJDS%ALUG_L= 1.d0
!
      end subroutine reset_aiccg_press
!
! ----------------------------------------------------------------------
!
      end module m_velo_matrix
