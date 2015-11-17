!
!      module set_nodal_boundary
!
!      Written by H. Matsui on july, 2005
!      Modified by H. Matsui on Jan., 2009
!
!!      subroutine set_nod_bc_from_ctl(nod_grp, numnod, num_phys_bc,    &
!!     &          ii, i, ibc_id, ibc, ibc2, bc_id_apt, bc_magnitude )
!!      subroutine set_nod_bc_from_data(nod_grp, numnod, num_phys_bc,   &
!!     &          ii, i, ibc_id, ibc, ibc2, bc_id_apt, field_name)
!!      subroutine set_fixed_bc_4_par_temp                              &
!!     &         (numnod, ncomp_nod, i_ref_t, d_nod)
!!      subroutine set_potential_4_fixed_press                          &
!!     &         (numnod, ncomp_nod, i_press, i_p_phi, d_nod)
!!      subroutine set_potential_4_sgs_press                            &
!!     &         (numnod, ncomp_nod, i_press, i_p_phi, d_nod)
!
      module set_nodal_boundary
!
      use m_precision
!
      implicit none
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine set_nod_bc_from_ctl(nod_grp, numnod, num_phys_bc,      &
     &          ii, i, ibc_id, ibc, ibc2, bc_id_apt, bc_magnitude )
!
      use t_group_data
!
      integer(kind = kint), intent(in) :: numnod
      type(group_data), intent(in) :: nod_grp
!
      integer(kind = kint), intent(in) :: i
      integer(kind = kint), intent(in) :: num_phys_bc
      real ( kind = kreal), intent(in) :: bc_magnitude
!
      integer(kind = kint), intent(inout) :: ii
      integer(kind = kint), intent(inout) :: ibc_id(num_phys_bc)
      integer(kind = kint), intent(inout) :: ibc(numnod)
      integer(kind = kint), intent(inout) :: ibc2(numnod)
      real ( kind = kreal), intent(inout) :: bc_id_apt(num_phys_bc)
!
      integer(kind = kint) :: k
!
!
      do k=1, nod_grp%istack_grp(i)-nod_grp%istack_grp(i-1)
        ii=ii+1
!
        ibc_id(ii)=nod_grp%item_grp(k+nod_grp%istack_grp(i-1))
        bc_id_apt(ii)=bc_magnitude
!
      end do
!
      if ( bc_magnitude .ne. 0.0d0 ) then
        do k=1, nod_grp%istack_grp(i)-nod_grp%istack_grp(i-1)
         ibc(nod_grp%item_grp(k+nod_grp%istack_grp(i-1)) ) = 1
        end do
      end if
!
      do k=1, nod_grp%istack_grp(i)-nod_grp%istack_grp(i-1)
        ibc2(nod_grp%item_grp(k+nod_grp%istack_grp(i-1)) ) = 1
      end do
!
      end subroutine set_nod_bc_from_ctl
!
!  ---------------------------------------------------------------------
!
      subroutine set_nod_bc_from_data(nod_grp, numnod, num_phys_bc,     &
     &          ii, i, ibc_id, ibc, ibc2, bc_id_apt, field_name)
!
      use t_group_data
      use m_boundary_field_IO
!
      integer(kind = kint), intent(in) :: numnod
      type(group_data), intent(in) :: nod_grp
!
      character(len=kchara), intent(in) :: field_name
      integer(kind = kint), intent(in) :: i
      integer(kind = kint), intent(in) :: num_phys_bc
!
      integer(kind = kint), intent(inout) :: ii
      integer(kind = kint), intent(inout) :: ibc_id(num_phys_bc)
      integer(kind = kint), intent(inout) :: ibc(numnod), ibc2(numnod)
      real ( kind = kreal), intent(inout) :: bc_id_apt(num_phys_bc)
!
      integer(kind = kint) :: k, ja, ia
!
!
      do ia = 1, num_bc_group_IO
        if(bc_group_type_IO(ia) .eq. flag_nod_grp) then
          if ( bc_data_group_IO(ia) .eq. nod_grp%grp_name(i)            &
     &       .and. bc_field_type_IO(ia) .eq. field_name ) then
!
            do k=1, nod_grp%istack_grp(i)-nod_grp%istack_grp(i-1)
              ja = istack_bc_data_IO(ia-1) + k
              ii=ii+1
!
              ibc_id(ii)=nod_grp%item_grp(k+nod_grp%istack_grp(i-1))
              bc_id_apt(ii)=boundary_field_IO(ja)
!
              ibc( nod_grp%item_grp(k+nod_grp%istack_grp(i-1)) ) = 1
              ibc2(nod_grp%item_grp(k+nod_grp%istack_grp(i-1)) ) = 1
            end do
!
          end if
        end if
      end do
!
      end subroutine set_nod_bc_from_data
!
!  ---------------------------------------------------------------------
!
      subroutine set_fixed_bc_4_par_temp                                &
     &         (numnod, ncomp_nod, i_ref_t, d_nod)
!
      use m_bc_data_ene
!
      integer (kind = kint), intent(in) :: numnod, ncomp_nod, i_ref_t
      real(kind = kreal), intent(inout) :: d_nod(numnod,ncomp_nod)
!
      integer (kind = kint) :: inum, inod
!
      do inum = 1, nod_bc1_t%num_bc_nod
        inod = nod_bc1_t%ibc_id(inum)
        nod_bc1_t%bc_apt(inum) = nod_bc1_t%bc_apt(inum)                 &
     &                          - d_nod(inod,i_ref_t)
      end do
!
      end subroutine set_fixed_bc_4_par_temp
!
!  ---------------------------------------------------------------------
!
      subroutine set_potential_4_fixed_press                            &
     &         (numnod, ncomp_nod, i_press, i_p_phi, d_nod)
!
      use m_bc_data_velo
      use m_t_int_parameter
      use m_physical_property
!
      integer (kind = kint), intent(in) :: numnod, ncomp_nod
      integer (kind = kint), intent(in) :: i_press, i_p_phi
      real(kind = kreal), intent(inout) :: d_nod(numnod,ncomp_nod)
!
      integer (kind = kint) :: inum, inod
!
      do inum = 1, nod_bc1_p%num_bc_nod
        inod = nod_bc1_p%ibc_id(inum)
        nod_bc1_p%bc_apt(inum)                                          &
     &        =   -dt * nod_bc1_p%bc_apt(inum) * coef_press
        d_nod(inod,i_p_phi) = -dt * coef_press * d_nod(inod,i_press)
      end do
!
      end subroutine set_potential_4_fixed_press
!
!  ---------------------------------------------------------------------
!
      subroutine set_potential_4_sgs_press                              &
     &         (numnod, ncomp_nod, i_press, i_p_phi, d_nod)
!
      use m_bc_data_velo
      use m_t_int_parameter
      use m_physical_property
!
      integer (kind = kint), intent(in) :: numnod, ncomp_nod
      integer (kind = kint), intent(in) :: i_press, i_p_phi
      real(kind = kreal), intent(inout) :: d_nod(numnod,ncomp_nod)
!
       integer (kind = kint) :: inum, inod
!
      do inum = 1, sgs_bc1_p%num_bc_nod
        inod = sgs_bc1_p%ibc_id(inum)
        sgs_bc1_p%bc_apt(inum)                                          &
     &       = -dt * sgs_bc1_p%bc_apt(inum) * coef_press
        d_nod(inod,i_p_phi) =  -dt * coef_press * d_nod(inod,i_press)
      end do
!
      end subroutine set_potential_4_sgs_press
!
!  ---------------------------------------------------------------------
!
      end module set_nodal_boundary
