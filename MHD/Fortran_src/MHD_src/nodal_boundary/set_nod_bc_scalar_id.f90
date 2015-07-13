!set_nod_bc_scalar_id.f90
!      module set_nod_bc_scalar_id
!
!        programmed by H.Matsui and H.Okuda
!                                    on July 2000 (ver 1.1)
!        modified by Kemorin on Jan. 2004
!        modified by Kemorin on Oct. 2005
!
!      subroutine set_fixed_bc_scalar_id(num_bc_field, bc_field_name,   &
!     &          ibc_field_type, bc_field_mag, ibc, ibc2, num_bc_nod,   &
!     &          ibc_id, bc_apt, field_name, ii)
!      subroutine set_bc_scalar_id(num_bc_field, bc_field_name,         &
!     &          ibc_field_type, bc_field_mag, ibc, ibc2, num_bc_nod,   &
!     &          ibc_id, bc_apt, iref, ii)
!      subroutine set_bc_sph_magne_p_id(num_bc_field, bc_field_name,    &
!     &          ibc_field_type, ii)
!
      module set_nod_bc_scalar_id
!
      use m_precision
!
      use m_geometry_parameter
      use m_node_group
      use set_nodal_boundary
!
      implicit none
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine set_fixed_bc_scalar_id(num_bc_field, bc_field_name,    &
     &          ibc_field_type, bc_field_mag, ibc, ibc2, num_bc_nod,    &
     &          ibc_id, bc_apt, field_name, ii)
!
      integer (kind=kint), intent(in) :: num_bc_field
      real (kind=kreal), intent(in) :: bc_field_mag(num_bc_field)
      integer (kind=kint), intent(in) :: ibc_field_type(num_bc_field)
      character(len=kchara), intent(in) :: bc_field_name(num_bc_field)
      character(len=kchara), intent(in) :: field_name
!
      integer (kind=kint), intent(inout) :: ibc(numnod)
      integer (kind=kint), intent(inout) :: ibc2(numnod)
!
      integer (kind=kint), intent(in) :: num_bc_nod
      integer (kind=kint), intent(inout) :: ibc_id(num_bc_nod)
      real (kind=kreal),   intent(inout) :: bc_apt(num_bc_nod)
! 
!
      integer (kind = kint), intent(inout) :: ii
      integer (kind = kint) :: i, j
!
!
      do i=1, nod_grp1%num_grp 
!
! ----------- loop for boundary conditions
!
        do j=1, num_bc_field
!
! ----------- check node group
!
          if (bc_name(i)==bc_field_name(j)) then
!
! -----------set boundary from control file
!
            if ( ibc_field_type(j) .eq. 1 ) then
!
              call set_nod_bc_from_ctl (num_bc_nod, ii, i, ibc_id,      &
     &              ibc, ibc2, bc_apt, bc_field_mag(j) )
!
! -----------set boundary from data file
!
            else if ( ibc_field_type(j).eq. -1 ) then
!
              call set_nod_bc_from_data (num_bc_nod, ii, i, ibc_id,     &
     &              ibc, ibc2, bc_apt, field_name)
            end if
!
          end if
        end do
      end do
!
      end subroutine set_fixed_bc_scalar_id
!
!  ---------------------------------------------------------------------
!
      subroutine set_bc_scalar_id(num_bc_field, bc_field_name,          &
     &          ibc_field_type, bc_field_mag, ibc, ibc2, num_bc_nod,    &
     &          ibc_id, bc_apt, iref, ii)
!
      integer (kind=kint), intent(in) :: iref
!
      integer (kind=kint), intent(in) :: num_bc_field
      real (kind=kreal), intent(in) :: bc_field_mag(num_bc_field)
      integer (kind=kint), intent(in) :: ibc_field_type(num_bc_field)
      character (len=kchara), intent(in) :: bc_field_name(num_bc_field)
!
      integer (kind=kint), intent(inout) :: ibc(numnod)
      integer (kind=kint), intent(inout) :: ibc2(numnod)
!
      integer (kind=kint), intent(in) :: num_bc_nod
      integer (kind=kint), intent(inout) :: ibc_id(num_bc_nod)
      real (kind=kreal),   intent(inout) :: bc_apt(num_bc_nod)
!
      integer (kind = kint), intent(inout) :: ii
      integer (kind = kint) :: i, j
!
      ii = 0
!
      do i=1, nod_grp1%num_grp 
        do j=1, num_bc_field
          if (bc_name(i)==bc_field_name(j)) then
!
            if ( ibc_field_type(j) .eq. iref) then
              call set_nod_bc_from_ctl (num_bc_nod, ii, i, ibc_id,      &
     &              ibc, ibc2, bc_apt, bc_field_mag(j) )
            end if
!
          end if
        end do
      end do
!
      end subroutine set_bc_scalar_id
!
!  ---------------------------------------------------------------------
!
      subroutine set_bc_sph_magne_p_id(num_bc_field, bc_field_name,     &
     &          ibc_field_type, ii)
!
      use set_mag_p_sph
!
      integer (kind=kint), intent(in) :: num_bc_field
      integer (kind=kint), intent(in) :: ibc_field_type(num_bc_field)
      character (len=kchara), intent(in) :: bc_field_name(num_bc_field)
!
      integer (kind = kint), intent(inout) :: ii
      integer (kind = kint) :: i, j
!
      do i=1, nod_grp1%num_grp 
        do j=1, num_bc_field
          if (bc_name(i)==bc_field_name(j)) then
!
            if ( ibc_field_type(j) .eq. 999 ) then
              call s_set_mag_p_sph(ii, i, j)
            end if
!
          end if
        end do
      end do
!
      end subroutine set_bc_sph_magne_p_id
!
!  ---------------------------------------------------------------------
!
      end module set_nod_bc_scalar_id
