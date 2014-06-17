!
!      module mesh_interpolation
!
      module mesh_interpolation
!
!     Written by H. Matsui on Sep., 2006
!
      use m_precision
!
        implicit none
!
!      subroutine interpolation_4_mesh_test
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine interpolation_4_mesh_test
!
      use calypso_mpi
      use m_machine_parameter
      use m_ctl_params_4_gen_table
      use m_interpolated_geometry
      use interpolate_position
      use m_read_mesh_data
!
      use check_ineterppolated_mesh
!
!     return global node from table
!
      if (iflag_debug.eq.1)   write(*,*) 's_interpolate_global_node'
      call s_interpolate_global_node(nnod_2nd)
!
!     interpolate 2nd mesh from 1st mesh
!
      if (iflag_debug.eq.1)   write(*,*) 's_interpolate_position'
      call s_interpolate_position
!      if (iflag_debug.eq.1)   write(*,*) 's_interpolate_position_by_N'
!      call s_interpolate_position_by_N
!      if (iflag_debug.eq.1)   write(*,*) 's_interpolate_position_by_s'
!      call s_interpolate_position_by_s
!
!
      if (iflag_debug.gt.0)  write(*,*) 's_check_ineterppolated_mesh'
      call s_check_ineterppolated_mesh
!
      end subroutine interpolation_4_mesh_test
!
! ----------------------------------------------------------------------
!
      end module mesh_interpolation
