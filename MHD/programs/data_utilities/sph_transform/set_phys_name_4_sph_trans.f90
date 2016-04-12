!>@file   set_phys_name_4_sph_trans.f90
!!@brief  module set_phys_name_4_sph_trans
!!
!!@author H. Matsui
!!@date Programmed in Jan., 2008
!
!!@verbatim
!!      subroutine copy_sph_name_rj_to_rtp
!!@endverbatim
!
      module set_phys_name_4_sph_trans
!
      use m_precision
!
      implicit  none
!
!>      Number of fields on spherical grid @f$ f(r,\theta,\phi) @f$
      integer (kind=kint) :: num_phys_rtp
!>      Field name for @f$ f(r,\theta,\phi) @f$
      character (len=kchara), allocatable :: phys_name_rtp(:)
!
!>      Start field address of scalar fields @f$ f(r,\theta,\phi) @f$
      integer (kind=kint) :: istart_scalar_rtp
!>      Start field address of vector fields @f$ f(r,\theta,\phi) @f$
      integer (kind=kint) :: istart_vector_rtp
!>      Start field address of tensor fields @f$ f(r,\theta,\phi) @f$
      integer (kind=kint) :: istart_tensor_rtp
!
      private :: allocate_phys_rtp_name
!
! -------------------------------------------------------------------
!
      contains
!
! -------------------------------------------------------------------
!
      subroutine allocate_phys_rtp_name
!
      allocate( phys_name_rtp(num_phys_rtp) )
!
      end subroutine allocate_phys_rtp_name
!
!  --------------------------------------------------------------------
!
      subroutine deallocate_phys_rtp_name
!
      deallocate( phys_name_rtp )
!
      end subroutine deallocate_phys_rtp_name
!
!  --------------------------------------------------------------------
!  --------------------------------------------------------------------
!
      subroutine copy_sph_name_rj_to_rtp
!
      use m_constants
      use m_machine_parameter
      use m_phys_constants
      use m_sph_spectr_data
      use m_work_4_sph_trans
!
      integer(kind = kint) :: i, i0
!
      num_phys_rtp =  num_phys_rj
      call allocate_phys_rtp_name
!
      i0 = 0
      num_scalar_rtp = 0
      do i = 1, num_phys_rj
        if (num_phys_comp_rj(i) .eq. n_scalar) then
          i0 = i0 + 1
          num_scalar_rtp = num_scalar_rtp + 1
          phys_name_rtp(i0) =        rj_fld1%phys_name(i)
        end if
      end do
      istart_scalar_rtp = 1
!
      num_vector_rtp = 0
      do i = 1, num_phys_rj
        if (num_phys_comp_rj(i) .eq. n_vector) then
          i0 = i0 + 1
          num_vector_rtp = num_vector_rtp + 1
          phys_name_rtp(i0) =        rj_fld1%phys_name(i)
        end if
      end do
      istart_vector_rtp = istart_scalar_rtp + num_scalar_rtp
!
      num_tensor_rtp = 0
      do i = 1, num_phys_rj
        if (num_phys_comp_rj(i) .eq. n_sym_tensor) then
          i0 = i0 + 1
          num_tensor_rtp = num_tensor_rtp + 1
          phys_name_rtp(i0) =        rj_fld1%phys_name(i)
        end if
      end do
      istart_tensor_rtp = istart_vector_rtp + num_vector_rtp
!
      if (iflag_debug .gt. 0) then
!        write(*,*) 'num_phys_rj', num_phys_rj
!        write(*,*) 'id, components, stack, phys_name_rj'
!        do i = 1, num_phys_rj
!          write(*,*) i, num_phys_comp_rj(i), istack_phys_comp_rj(i),   &
!     &              trim(rj_fld1%phys_name(i))
!        end do
        write(*,*)
        write(*,*) 'num_phys_rtp', num_phys_rtp
        write(*,*) 'phys_name_rtp', size(phys_name_rtp)
        write(*,*) 'id, components, stack, phys_name_rtp'
        do i = 1, num_phys_rtp
          write(*,*) i, trim(phys_name_rtp(i))
        end do
        write(*,*) 'istart_scalar_rtp',                                 &
     &              istart_scalar_rtp, num_scalar_rtp
        write(*,*) 'istart_vector_rtp',                                 &
     &              istart_vector_rtp, num_vector_rtp
        write(*,*) 'istart_tensor_rtp',                                 &
     &              istart_tensor_rtp, num_tensor_rtp
      end if
!
!
      end subroutine copy_sph_name_rj_to_rtp
!
! -------------------------------------------------------------------
!
      end module set_phys_name_4_sph_trans
 