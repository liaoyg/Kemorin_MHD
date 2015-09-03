!add_nodal_field_type.f90
!     module add_nodal_field_type
!
!      Written by H. Matsui
!
!      subroutine add_scalar_type_smp(node, fld, i_v1, i_v2, i_r)
!      subroutine add_vector_type_smp(node, fld, i_v1, i_v2, i_r)
!      subroutine add_tensor_type_smp(node, fld, i_v1, i_v2, i_r)
!
!         d_nod(inod,i_r) =  d_nod(inod,i_v1) + d_nod(inod,i_v2)
!        i_r: result field ID
!        i_v1, i_v2: source field IDs
!
      module add_nodal_field_type
!
      use m_precision
!
      implicit none
!
!      private :: add_scalar_type_smp, add_vector_array_smp
!      private :: add_tensor_array_smp
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine add_scalar_type_smp(node, fld, i_v1, i_v2, i_r)
!
      use m_machine_parameter
      use t_geometry_data
      use t_phys_data
!
      integer(kind = kint), intent(in) :: i_r, i_v1, i_v2
      type(node_data), intent(in) :: node
      type(phys_data), intent(inout) :: fld
!
!
      call add_scalar_array_smp(np_smp,                                 &
     &    node%numnod, node%istack_nod_smp, fld%ntot_phys,              &
     &    fld%d_fld, i_v1, i_v2, i_r)
!
       end subroutine add_scalar_type_smp
!
!-----------------------------------------------------------------------
!
      subroutine add_vector_type_smp(node, fld, i_v1, i_v2, i_r)
!
      use m_machine_parameter
      use t_geometry_data
      use t_phys_data
!
      integer(kind = kint), intent(in) :: i_r, i_v1, i_v2
      type(node_data), intent(in) :: node
      type(phys_data), intent(inout) :: fld
!
!
      call add_vector_array_smp(np_smp,                                 &
     &    node%numnod, node%istack_nod_smp, fld%ntot_phys, fld%d_fld,   &
     &    i_v1, i_v2, i_r)
!
       end subroutine add_vector_type_smp
!
!-----------------------------------------------------------------------
!
      subroutine add_tensor_type_smp(node, fld, i_v1, i_v2, i_r)
!
      use m_machine_parameter
      use t_geometry_data
      use t_phys_data
!
      integer(kind = kint), intent(in) :: i_r, i_v1, i_v2
      type(node_data), intent(in) :: node
      type(phys_data), intent(inout) :: fld
!
!
      call add_tensor_array_smp(np_smp,                                 &
     &    node%numnod, node%istack_nod_smp, fld%ntot_phys, fld%d_fld,   &
     &    i_v1, i_v2, i_r)
!
      end subroutine add_tensor_type_smp
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine add_scalar_array_smp(np_smp, nnod, inod_smp_stack,     &
     &          ntot_comp, d_phys, i_v1, i_v2, i_r)
!
      use cal_add_smp
!
      integer (kind=kint), intent(in) :: np_smp, nnod, ntot_comp
      integer (kind=kint), intent(in) :: inod_smp_stack(0:np_smp)
      integer(kind = kint), intent(in) :: i_r, i_v1, i_v2
!
      real (kind=kreal), intent(inout) :: d_phys(nnod,ntot_comp)
!
!
!$omp parallel
      call add_scalars_smp(np_smp, nnod, inod_smp_stack,               &
     &    d_phys(1,i_v1), d_phys(1,i_v2), d_phys(1,i_r))
!$omp end parallel
!
       end subroutine add_scalar_array_smp
!
!-----------------------------------------------------------------------
!
      subroutine add_vector_array_smp(np_smp, nnod, inod_smp_stack,     &
     &          ntot_comp, d_phys, i_v1, i_v2, i_r)
!
      use cal_add_smp
!
      integer (kind=kint), intent(in) :: np_smp, nnod, ntot_comp
      integer (kind=kint), intent(in) :: inod_smp_stack(0:np_smp)
      integer(kind = kint), intent(in) :: i_r, i_v1, i_v2
!
      real (kind=kreal), intent(inout) :: d_phys(nnod,ntot_comp)
!
!
!$omp parallel
      call add_vectors_smp(np_smp, nnod, inod_smp_stack,                &
     &    d_phys(1,i_v1), d_phys(1,i_v2), d_phys(1,i_r))
!$omp end parallel
!
       end subroutine add_vector_array_smp
!
!-----------------------------------------------------------------------
!
      subroutine add_tensor_array_smp(np_smp, nnod, inod_smp_stack,     &
     &          ntot_comp, d_phys, i_v1, i_v2, i_r)
!
      use cal_add_smp
!
      integer (kind=kint), intent(in) :: np_smp, nnod, ntot_comp
      integer (kind=kint), intent(in) :: inod_smp_stack(0:np_smp)
      integer(kind = kint), intent(in) :: i_r, i_v1, i_v2
!
      real (kind=kreal), intent(inout) :: d_phys(nnod,ntot_comp)
!
!
!$omp parallel
      call add_tensors_smp(np_smp, nnod, inod_smp_stack,                &
     &    d_phys(1,i_v1), d_phys(1,i_v2), d_phys(1,i_r))
!$omp end parallel
!
      end subroutine add_tensor_array_smp
!
!-----------------------------------------------------------------------
!
      end module add_nodal_field_type
