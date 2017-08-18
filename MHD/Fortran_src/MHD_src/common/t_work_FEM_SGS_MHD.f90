!>@file  t_work_FEM_MHD.f90
!!       module t_work_FEM_MHD
!!
!!@author H. Matsui
!!@date   Programmed in Aug., 2011
!
!> @brief Structures for Dynamic SGS model in FEM dynamo
!!
!!@verbatim
!!@endverbatim
      module t_work_FEM_SGS_MHD
!
      use m_precision
      use m_constants
!
      use t_MHD_finite_element_mat
      use t_work_FEM_integration
      use t_work_FEM_dynamic_SGS
!
      implicit  none
!
!
      type work_FEM_SGS_MHD
!>        Stracture for FEM assembling
        type(arrays_finite_element_mat) :: rhs_mat
!>        Work array for FEM assemble in MHD model
        type(work_MHD_fe_mat) :: mhd_fem_wk
!>        Work area for SGS model
        type(work_FEM_dynamic_SGS) :: FEM_SGS_wk
      end type work_FEM_SGS_MHD
!
      end module t_work_FEM_SGS_MHD