!set_export_cube.f90
!     module set_export_cube
!
!     Written by H. Matsui
!     modified by H. Matsui on Aug., 2007
!
!      subroutine set_export_data(ipe, jpe, kpe)
!      subroutine set_export_data_quad(ipe, jpe, kpe)
!
      module set_export_cube
!
      use m_precision
!
      use m_size_of_cube
      use m_comm_data_cube_kemo
      use count_export_inside_cube
      use set_export_inside_cube
      use count_export_peri
      use set_export_peri_cube
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
       subroutine set_export_data(ipe, jpe, kpe)
!
      integer (kind = kint) :: ipe, jpe, kpe
      integer (kind = kint) :: inod
!
!
            neibpetot = 0
            inod   = 0

            call count_export_inside(inod)
!
            call count_export_peri_linear(ipe, jpe, inod)

            num_export = stack_export(neibpetot)
!
!
            inod = 0
            neibpetot = 0

            call set_export_inside(inod)
!
            call set_export_peri(ipe, jpe, inod)
!
          end subroutine set_export_data
!
! ----------------------------------------------------------------------
!
      subroutine set_export_data_quad(ipe, jpe, kpe)
!
      integer (kind = kint) :: ipe, jpe, kpe
      integer (kind = kint) :: inod
!
! ***** set and write export nodes
!                                     .... count nodes 
            inod = 0
            neibpetot = 0
!
            call count_export_inside_quad(kpe, inod)

            call count_export_peri_quad(ipe, jpe, kpe, inod)

            write(*,*) stack_export
            num_export = stack_export(neibpetot)
!
            inod = 0
            neibpetot = 0

            call set_export_inside_quad(kpe, inod)

!
            call set_export_peri_quad(ipe, jpe, kpe, inod)
!
      end subroutine set_export_data_quad
!
! ----------------------------------------------------------------------
!
      end module set_export_cube
