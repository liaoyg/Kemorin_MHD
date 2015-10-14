!set_nodal_bc_4_velo.f90
!      module set_nodal_bc_4_velo
!
!      Written by H. Matsui on july, 2005
!
!!      subroutine set_rotation_boundary(numnod, xx,                    &
!!     &          num_bc_v10_nod, ibc_v10_id, bc_v10_id_apt,            &
!!     &          ncomp_nod, i_velo, d_nod)
!!      subroutine set_specific_boundary_velo(numnod, xx,               &
!!     &          num_bc_vsp_nod, ibc_vsp_id, bc_vsp_id_apt,            &
!!     &          ncomp_nod, i_velo, d_nod)
!!      subroutine set_fixed_bc_zero_ff_rot(num_phys_bc, ibc_id)
!!      subroutine set_specific_boundary_velo_rhs(num_bc_vsp_nod,       &
!!     &    ibc_vsp_id)
!
      module set_nodal_bc_4_velo
!
      use m_precision
      use m_constants
!
      implicit none
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine set_rotation_boundary(numnod, xx,                      &
     &          num_bc_v10_nod, ibc_v10_id, bc_v10_id_apt,              &
     &          ncomp_nod, i_velo, d_nod)
!
      integer(kind = kint), intent(in) :: numnod, ncomp_nod, i_velo
      real(kind = kreal), intent(in) :: xx(numnod,3)
!
       integer(kind=kint), intent(in) :: num_bc_v10_nod
       integer(kind=kint), intent(in) :: ibc_v10_id(num_bc_v10_nod)
       real(kind=kreal), intent(in) :: bc_v10_id_apt(num_bc_v10_nod,3)
!
      real(kind = kreal), intent(inout) :: d_nod(numnod,ncomp_nod)
!
      integer (kind = kint) :: inum, inod
!
!
       do inum=1, num_bc_v10_nod
         inod = ibc_v10_id(inum)
         d_nod(inod,i_velo  ) =  bc_v10_id_apt(inum,2)*xx(inod,3)       &
     &                         - bc_v10_id_apt(inum,3)*xx(inod,2)
         d_nod(inod,i_velo+1) =  bc_v10_id_apt(inum,3)*xx(inod,1)       &
     &                         - bc_v10_id_apt(inum,1)*xx(inod,3)
         d_nod(inod,i_velo+2) =  bc_v10_id_apt(inum,1)*xx(inod,2)       &
     &                         - bc_v10_id_apt(inum,2)*xx(inod,1)
       end do
!
      end subroutine set_rotation_boundary
!
!  ---------------------------------------------------------------------
!
      subroutine set_specific_boundary_velo(numnod, xx,                 &
     &          num_bc_vsp_nod, ibc_vsp_id, bc_vsp_id_apt,              &
     &          ncomp_nod, i_velo, d_nod)
!
       use m_constants
       use m_t_step_parameter
!
      integer(kind = kint), intent(in) :: numnod, ncomp_nod, i_velo
      real(kind = kreal), intent(in) :: xx(numnod,3)
!
      integer(kind=kint), intent(in) :: num_bc_vsp_nod
      integer(kind=kint), intent(in) :: ibc_vsp_id(num_bc_vsp_nod)
      real(kind=kreal), intent(in) :: bc_vsp_id_apt(num_bc_vsp_nod)
!
      real(kind = kreal), intent(inout) :: d_nod(numnod,ncomp_nod)
!
      integer (kind = kint) :: inum, inod
!
!
       do inum=1, num_bc_vsp_nod
         inod = ibc_vsp_id(inum)
         d_nod(inod,i_velo  ) = 0.0d0
         d_nod(inod,i_velo+1) = 0.0d0
         d_nod(inod,i_velo+2) = 0.01 * cos( four*xx(inod,1)             &
     &                               - bc_vsp_id_apt(inum) * time )
       end do
!
      end subroutine set_specific_boundary_velo
!
!  ---------------------------------------------------------------------
!
      subroutine set_fixed_bc_zero_ff_rot(numnod,                       &
     &          num_phys_bc, ibc_id, ff, ff_nl)
!
      integer(kind = kint), intent(in) :: numnod
      integer(kind = kint), intent(in) :: num_phys_bc
      integer(kind = kint), intent(in) :: ibc_id(num_phys_bc)
!
      real(kind = kreal), intent(inout) :: ff(numnod,3)
      real(kind = kreal), intent(inout) :: ff_nl(numnod,3)
!
      integer(kind = kint) :: inum, inod
!
!
       do inum=1, num_phys_bc
         inod = ibc_id(inum)
         ff(inod,1)=    0.0d0
         ff(inod,2)=    0.0d0
         ff(inod,3)=    0.0d0
         ff_nl(inod,1)= 0.0d0
         ff_nl(inod,2)= 0.0d0
         ff_nl(inod,3)= 0.0d0
       end do
!
      end subroutine set_fixed_bc_zero_ff_rot
!
!  ---------------------------------------------------------------------
!
      subroutine set_specific_boundary_velo_rhs(numnod,                 &
     &          num_bc_vsp_nod, ibc_vsp_id, ff)
!
      integer(kind = kint), intent(in) :: numnod
      integer(kind = kint), intent(in) :: num_bc_vsp_nod
      integer(kind = kint), intent(in) :: ibc_vsp_id(num_bc_vsp_nod)
      real(kind = kreal), intent(inout) :: ff(numnod,3)
!
      integer (kind = kint) :: inum, inod
!
!
       do inum=1, num_bc_vsp_nod
         inod = ibc_vsp_id(inum)
         ff(inod,1)= two * ff(inod,1)
         ff(inod,2)= two * ff(inod,2)
         ff(inod,3)= two * ff(inod,3)
       end do
!
      end subroutine set_specific_boundary_velo_rhs
!
!  ---------------------------------------------------------------------
!
      end module set_nodal_bc_4_velo
