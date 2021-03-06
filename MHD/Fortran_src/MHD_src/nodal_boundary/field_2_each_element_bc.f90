!field_2_each_element_bc.f90
!     module field_2_each_element_bc
!
!      Written by H. Matsui and H.Okuda  on July 2001
!      Modified by H. Matsui on Sep., 2005
!
!!      subroutine scalar_2_element_4_boundary(numele, nnod_4_ele, ie,  &
!!     &          num_idx_ibc, ele_bc_id, ibc_stack_smp,                &
!!     &          k2, i_comp, scalar_e )
!!      subroutine vector_2_element_4_boundary(numnod, numele,          &
!!     &          nnod_4_ele, ie, num_idx_ibc, ele_bc_id, ibc_stack_smp,&
!!     &          k2, ncomp_nod, i_vect, d_nod, vector_e)
!
      module field_2_each_element_bc
!
      use m_precision
      use m_machine_parameter
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine scalar_2_element_4_boundary(numnod, numele,            &
     &          nnod_4_ele, ie, num_idx_ibc, ele_bc_id, ibc_stack_smp,  &
     &          k2, ncomp_nod, i_comp, d_nod, scalar_e)
!
      integer(kind = kint), intent(in) :: numele, nnod_4_ele
      integer(kind = kint), intent(in) :: ie(numele,nnod_4_ele)
!
      integer (kind=kint), intent(in) :: num_idx_ibc
      integer (kind=kint), intent(in) :: ele_bc_id(num_idx_ibc)
      integer (kind=kint), intent(in) :: ibc_stack_smp(0:np_smp)
      integer(kind = kint), intent(in) :: numnod, ncomp_nod, i_comp
      real(kind = kreal), intent(in) :: d_nod(numnod,ncomp_nod)
!
      real (kind=kreal), intent(inout) :: scalar_e(numele)
!
      integer (kind=kint) :: iproc, inod, iele, inum, k2
      integer (kind=kint) :: istart, iend
!
!
!$omp parallel do private(inod,inum,iele,istart,iend) 
        do iproc = 1, np_smp
         istart = ibc_stack_smp( iproc-1 ) + 1
         iend   = ibc_stack_smp( iproc )
!
!cdir nodep
         do inum = istart, iend
          iele = ele_bc_id(inum)
          inod = ie(iele,k2)
          scalar_e(iele) = d_nod(inod,i_comp)
         end do
        end do
!$omp end parallel do
!
!
      end subroutine scalar_2_element_4_boundary
!
!-----------------------------------------------------------------------
!
      subroutine vector_2_element_4_boundary(numnod, numele,            &
     &          nnod_4_ele, ie, num_idx_ibc, ele_bc_id, ibc_stack_smp,  &
     &          k2, ncomp_nod, i_vect, d_nod, vector_e)
!
      integer(kind = kint), intent(in) :: numele, nnod_4_ele
      integer(kind = kint), intent(in) :: ie(numele,nnod_4_ele)
!
      integer(kind = kint), intent(in) :: num_idx_ibc
      integer(kind = kint), intent(in) :: ele_bc_id(num_idx_ibc)
      integer(kind = kint), intent(in) :: ibc_stack_smp(0:np_smp)
      integer(kind = kint), intent(in) :: numnod, ncomp_nod, i_vect
      real(kind = kreal), intent(in) :: d_nod(numnod,ncomp_nod)
!
!
      real (kind=kreal), intent(inout) :: vector_e(numele,3)
!
      integer (kind=kint) :: iproc, inod, iele, inum, k2
      integer (kind=kint) :: istart, iend
!
!
!$omp parallel do private(inod,inum,iele,istart,iend) 
        do iproc = 1, np_smp
         istart = ibc_stack_smp( iproc-1 ) + 1
         iend   = ibc_stack_smp( iproc )
!
!cdir nodep
!VOPTION INDEP, VEC
         do inum = istart, iend
          iele = ele_bc_id(inum)
          inod = ie(iele,k2)
          vector_e(iele,1) = d_nod(inod,i_vect  )
          vector_e(iele,2) = d_nod(inod,i_vect+1)
          vector_e(iele,3) = d_nod(inod,i_vect+2)
         end do
        end do
!$omp end parallel do
!
!
      end subroutine vector_2_element_4_boundary
!
!-----------------------------------------------------------------------
!
      end module field_2_each_element_bc
