!copy_nodal_type_4_sph_trans.f90
!      module copy_nodal_type_4_sph_trans
!
!      Written by H. Matsui on Nov., 2012
!
!      subroutine copy_scalar_t_from_trans_wpole(nscalar_trans, i_trns, &
!     &          i_field, node, nod_fld)
!      subroutine copy_xyz_vec_t_from_trans_wpole(nvector_trans, i_trns,&
!     &          i_field, node, nod_fld)
!      subroutine copy_xyz_tsr_t_from_trans_wpole(ncomp_trans, i_trns,  &
!     &         i_field, node, nod_fld)
!
!      subroutine copy_scalar_t_from_sph_trans(nscalar_trans, i_trns,   &
!     &          i_field, node, nod_fld)
!      subroutine copy_xyz_vec_t_from_sph_trans(nvector_trans, i_trns,  &
!     &          i_field, node, nod_fld)
!      subroutine copy_xyz_tsr_t_from_sph_trans(ncomp_trans, i_trns,    &
!     &          i_field, node, nod_fld)
!
!      subroutine copy_scalar_t_to_sph_trans(nscalar_trans, i_trns,     &
!     &          i_field, node, nod_fld)
!      subroutine copy_xyz_vec_t_to_sph_trans(ncomp_trans, i_trns,      &
!     &          i_field, node, nod_fld)
!      subroutine copy_xyz_tsr_t_to_sph_trans(ncomp_trans, i_trns,      &
!     &           i_field, node, nod_fld)
!
      module copy_nodal_type_4_sph_trans
!
      use m_precision
!
      use m_constants
      use m_machine_parameter
      use m_phys_constants
!
      use t_geometry_data
      use t_phys_data
!
      implicit  none
! 
! -------------------------------------------------------------------
!
      contains
!
! -------------------------------------------------------------------
!
      subroutine copy_scalar_t_from_trans_wpole(nscalar_trans, i_trns,  &
     &          i_field, node, nod_fld)
!
      use copy_xyz_field_4_sph_trans
!
      type(node_data), intent(in) :: node
      type(phys_data), intent(inout) :: nod_fld
!
      integer(kind = kint), intent(in) :: nscalar_trans
      integer(kind = kint), intent(in) :: i_field, i_trns
!
!
      call copy_scalar_from_trans_w_pole(node%numnod,                   &
     &    node%internal_node, node%xx, nscalar_trans, i_trns, i_field,  &
     &    nod_fld%ntot_phys, nod_fld%d_fld)
!
      end subroutine copy_scalar_t_from_trans_wpole
!
! -------------------------------------------------------------------
!
      subroutine copy_xyz_vec_t_from_trans_wpole(nvector_trans, i_trns, &
     &          i_field, node, nod_fld)
!
      use copy_xyz_field_4_sph_trans
!
      type(node_data), intent(in) :: node
      type(phys_data), intent(inout) :: nod_fld
!
      integer(kind = kint), intent(in) :: nvector_trans
      integer(kind = kint), intent(in) :: i_field, i_trns
!
!
      call copy_xyz_vec_from_trans_w_pole(node%numnod,                  &
     &    node%internal_node, node%istack_nod_smp, node%xx,             &
     &    node%theta, node%phi, nvector_trans, i_trns, i_field,         &
     &    nod_fld%ntot_phys, nod_fld%d_fld)
!
      end subroutine copy_xyz_vec_t_from_trans_wpole
!
! -------------------------------------------------------------------
!
      subroutine copy_xyz_tsr_t_from_trans_wpole(ncomp_trans, i_trns,   &
     &         i_field, node, nod_fld)
!
      use copy_xyz_field_4_sph_trans
!
      type(node_data), intent(in) :: node
      type(phys_data), intent(inout) :: nod_fld
!
      integer(kind = kint), intent(in) :: ncomp_trans
      integer(kind = kint), intent(in) :: i_field, i_trns
!
!
      call copy_xyz_tsr_from_trans_w_pole(node%numnod,                  &
     &    node%internal_node, node%istack_nod_smp, node%xx, node%rr,    &
     &    node%ss, node%a_r, node%a_s, ncomp_trans, i_trns, i_field,    &
     &    nod_fld%ntot_phys, nod_fld%d_fld)
!
      end subroutine copy_xyz_tsr_t_from_trans_wpole
!
! -------------------------------------------------------------------
! -------------------------------------------------------------------
!
      subroutine copy_scalar_t_from_sph_trans(nscalar_trans, i_trns,    &
     &          i_field, node, nod_fld)
!
      use copy_xyz_field_4_sph_trans
!
      type(node_data), intent(in) :: node
      type(phys_data), intent(inout) :: nod_fld
!
      integer(kind = kint), intent(in) :: nscalar_trans
      integer(kind = kint), intent(in) :: i_field, i_trns
!
!
      call copy_scalar_from_sph_trans(node%numnod, nscalar_trans,       &
     &    i_trns, i_field, nod_fld%ntot_phys, nod_fld%d_fld)
!
      end subroutine copy_scalar_t_from_sph_trans
!
! -------------------------------------------------------------------
!
      subroutine copy_xyz_vec_t_from_sph_trans(nvector_trans, i_trns,   &
     &          i_field, node, nod_fld)
!
      use copy_xyz_field_4_sph_trans
!
      type(node_data), intent(in) :: node
      type(phys_data), intent(inout) :: nod_fld
!
      integer(kind = kint), intent(in) :: nvector_trans
      integer(kind = kint), intent(in) :: i_field, i_trns
!
!
      call copy_xyz_vec_from_sph_trans(node%numnod,                     &
     &    node%istack_nod_smp, node%theta, node%phi, nvector_trans,     &
     &    i_trns, i_field, nod_fld%ntot_phys, nod_fld%d_fld)
!
      end subroutine copy_xyz_vec_t_from_sph_trans
!
! -------------------------------------------------------------------
!
      subroutine copy_xyz_tsr_t_from_sph_trans(ncomp_trans, i_trns,     &
     &          i_field, node, nod_fld)
!
      use copy_xyz_field_4_sph_trans
!
      type(node_data), intent(in) :: node
      type(phys_data), intent(inout) :: nod_fld
!
      integer(kind = kint), intent(in) :: ncomp_trans
      integer(kind = kint), intent(in) :: i_field, i_trns
!
!
      call copy_xyz_tsr_from_sph_trans                                  &
     &   (node%numnod, node%istack_nod_smp, node%xx, node%rr, node%ss,  &
     &    node%a_r, node%a_s, ncomp_trans, i_trns, i_field,             &
     &    nod_fld%ntot_phys, nod_fld%d_fld)
!
      end subroutine copy_xyz_tsr_t_from_sph_trans
!
! -------------------------------------------------------------------
! -------------------------------------------------------------------
!
      subroutine copy_scalar_t_to_sph_trans(nscalar_trans, i_trns,      &
     &          i_field, node, nod_fld)
!
      use copy_xyz_field_4_sph_trans
!
      type(node_data), intent(in) :: node
      type(phys_data), intent(in) :: nod_fld
!
      integer(kind = kint), intent(in) :: nscalar_trans
      integer(kind = kint), intent(in) :: i_field, i_trns
!
!
      call copy_scalar_to_sph_trans(node%numnod, nscalar_trans, i_trns, &
     &    i_field, nod_fld%ntot_phys, nod_fld%d_fld)
!
      end subroutine copy_scalar_t_to_sph_trans
!
! -------------------------------------------------------------------
!
      subroutine copy_xyz_vec_t_to_sph_trans(ncomp_trans, i_trns,       &
     &          i_field, node, nod_fld)
!
      use copy_xyz_field_4_sph_trans
!
      type(node_data), intent(in) :: node
      type(phys_data), intent(in) :: nod_fld
!
      integer(kind = kint), intent(in) :: ncomp_trans
      integer(kind = kint), intent(in) :: i_field, i_trns
!
!
      call copy_xyz_vec_to_sph_trans(node%numnod, node%istack_nod_smp,  &
     &    node%xx, node%rr, node%ss, node%a_r, node%a_s, ncomp_trans,   &
     &    i_trns, i_field, nod_fld%ntot_phys, nod_fld%d_fld)
!
      end subroutine copy_xyz_vec_t_to_sph_trans
!
! -------------------------------------------------------------------
!
      subroutine copy_xyz_tsr_t_to_sph_trans(ncomp_trans, i_trns,       &
     &           i_field, node, nod_fld)
!
      use copy_xyz_field_4_sph_trans
!
      type(node_data), intent(in) :: node
      type(phys_data), intent(in) :: nod_fld
!
      integer(kind = kint), intent(in) :: ncomp_trans
      integer(kind = kint), intent(in) :: i_field, i_trns
!
!
      call copy_xyz_tsr_to_sph_trans(node%numnod, node%istack_nod_smp,  &
     &    node%xx, node%rr, node%ss, node%a_r, node%a_s, ncomp_trans,   &
     &    i_trns, i_field, nod_fld%ntot_phys, nod_fld%d_fld)
!
      end subroutine copy_xyz_tsr_t_to_sph_trans
!
! -------------------------------------------------------------------
!
      end module copy_nodal_type_4_sph_trans
