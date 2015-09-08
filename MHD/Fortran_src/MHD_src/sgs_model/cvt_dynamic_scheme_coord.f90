!cvt_dynamic_scheme_coord.f90
!      module cvt_dynamic_scheme_coord
!
!     Written by H. Matsui on Oct. 2005
!     Modified by H. Matsui on Aug., 2007
!
!> @brief Change coordinate system for dynamic model
!
!      subroutine cvt_vector_dynamic_scheme_coord
!      subroutine cvt_tensor_dynamic_scheme_coord
!
      module cvt_dynamic_scheme_coord
!
      use m_precision
      use m_machine_parameter
!
      implicit none
!
      private :: convert_dynamic_vectors_2_sph
      private :: convert_dynamic_vectors_2_cyl
      private :: convert_dynamic_tensors_2_sph
      private :: convert_dynamic_tensors_2_cyl
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine cvt_vector_dynamic_scheme_coord
!
!
      use m_geometry_constants
      use m_control_parameter
      use m_geometry_data
      use m_node_phys_address
      use m_node_phys_data
!
!
      if(icoord_SGS_model_coef .eq. iflag_spherical) then
        call convert_dynamic_vectors_2_sph                              &
     &     (node1%numnod, node1%istack_nod_smp, node1%xx,               &
     &      node1%rr, node1%ss, node1%a_r, node1%a_s, num_tot_nod_phys, &
     &      iphys%i_sgs_simi, iphys%i_sgs_grad, iphys%i_sgs_grad_f,     &
     &      d_nod)
      else if(icoord_SGS_model_coef .eq. iflag_cylindrical) then
        call convert_dynamic_vectors_2_cyl                              &
     &     (node1%numnod, node1%istack_nod_smp, node1%xx,               &
     &      node1%ss, node1%a_s, num_tot_nod_phys,                      &
     &      iphys%i_sgs_simi, iphys%i_sgs_grad, iphys%i_sgs_grad_f,     &
     &      d_nod)
      end if
!
      end subroutine cvt_vector_dynamic_scheme_coord
!
!  ---------------------------------------------------------------------
!
      subroutine cvt_tensor_dynamic_scheme_coord
!
      use m_machine_parameter
      use m_geometry_constants
      use m_control_parameter
      use m_geometry_data
      use m_node_phys_address
      use m_node_phys_data
!
!
      if(icoord_SGS_model_coef .eq. iflag_spherical) then
        call convert_dynamic_tensors_2_sph                              &
     &     (node1%numnod, node1%istack_nod_smp, node1%xx,               &
     &      node1%rr, node1%ss, node1%a_r, node1%a_s, num_tot_nod_phys, &
     &      iphys%i_sgs_simi, iphys%i_sgs_grad, iphys%i_sgs_grad_f,     &
     &      d_nod)
      else if(icoord_SGS_model_coef .eq. iflag_cylindrical) then
      call convert_dynamic_tensors_2_cyl                                &
     &     (node1%numnod, node1%istack_nod_smp, node1%xx,               &
     &      node1%ss, node1%a_s, num_tot_nod_phys,                      &
     &      iphys%i_sgs_simi, iphys%i_sgs_grad, iphys%i_sgs_grad_f,     &
     &      d_nod)
      end if
!
      end subroutine cvt_tensor_dynamic_scheme_coord
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine convert_dynamic_vectors_2_cyl(numnod, inod_smp_stack,  &
     &          xx, s, a_s, ncomp_nod, i_sgs_simi, i_sgs_grad,          &
     &          i_sgs_grad_f, d_nod)
!
      use cvt_xyz_vector_2_cyl_smp
!
      integer (kind = kint), intent(in) :: numnod
      integer (kind = kint), intent(in) :: inod_smp_stack(0:np_smp)
      real(kind=kreal), intent(in) :: xx(numnod,3)
      real(kind=kreal), intent(in) :: s(numnod)
      real(kind=kreal), intent(in) :: a_s(numnod)
!
      integer (kind = kint), intent(in) :: ncomp_nod, i_sgs_simi 
      integer (kind = kint), intent(in) :: i_sgs_grad, i_sgs_grad_f
      real(kind=kreal), intent(inout) :: d_nod(numnod,ncomp_nod)
!
!
      if(iflag_debug.gt.0) write(*,*) 'convert cylindrical corrdinate'
!$omp parallel
      call overwrite_vector_2_cyl_smp(np_smp, numnod, inod_smp_stack,   &
     &    d_nod(1,i_sgs_simi), xx(1,1), xx(1,2), s, a_s)

      call overwrite_vector_2_cyl_smp(np_smp, numnod, inod_smp_stack,   &
     &    d_nod(1,i_sgs_grad), xx(1,1), xx(1,2), s, a_s)

      call overwrite_vector_2_cyl_smp(np_smp, numnod, inod_smp_stack,   &
     &    d_nod(1,i_sgs_grad_f), xx(1,1), xx(1,2), s, a_s)
!$omp end parallel
!
      end subroutine convert_dynamic_vectors_2_cyl
!
!  ---------------------------------------------------------------------
!
      subroutine convert_dynamic_vectors_2_sph(numnod,                  &
     &          inod_smp_stack, xx, r, s, a_r, a_s, ncomp_nod,          &
     &          i_sgs_simi, i_sgs_grad, i_sgs_grad_f, d_nod)
!
      use cvt_xyz_vector_2_sph_smp
!
      integer (kind = kint), intent(in) :: numnod
      integer (kind = kint), intent(in) :: inod_smp_stack(0:np_smp)
      real(kind=kreal), intent(in) :: xx(numnod,3)
      real(kind=kreal), intent(in) :: r(numnod)
      real(kind=kreal), intent(in) :: s(numnod)
      real(kind=kreal), intent(in) :: a_r(numnod)
      real(kind=kreal), intent(in) :: a_s(numnod)
!
      integer (kind = kint), intent(in) :: ncomp_nod, i_sgs_simi 
      integer (kind = kint), intent(in) :: i_sgs_grad, i_sgs_grad_f
      real(kind=kreal), intent(inout) :: d_nod(numnod,ncomp_nod)
!
!
      if(iflag_debug .gt. 0) write(*,*) 'convert spherical corrdinate'
!$omp parallel
      call overwrite_vector_2_sph_smp(np_smp, numnod, inod_smp_stack,   &
     &    d_nod(1,i_sgs_simi), xx(1,1), xx(1,2), xx(1,3),               &
     &    r, s, a_r, a_s)

      call overwrite_vector_2_sph_smp(np_smp, numnod, inod_smp_stack,   &
     &    d_nod(1,i_sgs_grad), xx(1,1), xx(1,2), xx(1,3),               &
     &    r, s, a_r, a_s)

      call overwrite_vector_2_sph_smp(np_smp, numnod, inod_smp_stack,   &
     &    d_nod(1,i_sgs_grad_f), xx(1,1), xx(1,2), xx(1,3),             &
     &    r, s, a_r, a_s)
!$omp end parallel
!
      end subroutine convert_dynamic_vectors_2_sph
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine convert_dynamic_tensors_2_cyl(numnod,                  &
     &          inod_smp_stack, xx, s, a_s, ncomp_nod,                  &
     &          i_sgs_simi, i_sgs_grad, i_sgs_grad_f, d_nod)
!
      use cvt_xyz_tensor_2_cyl_smp
!
      integer (kind = kint), intent(in) :: numnod
      integer (kind = kint), intent(in) :: inod_smp_stack(0:np_smp)
      real(kind=kreal), intent(in) :: xx(numnod,3)
      real(kind=kreal), intent(in) :: s(numnod)
      real(kind=kreal), intent(in) :: a_s(numnod)
!
      integer (kind = kint), intent(in) :: ncomp_nod, i_sgs_simi 
      integer (kind = kint), intent(in) :: i_sgs_grad, i_sgs_grad_f
      real(kind=kreal), intent(inout) :: d_nod(numnod,ncomp_nod)
!
!
      if(iflag_debug.gt.0) write(*,*) 'convert cylindrical corrdinate'
!$omp parallel
      call overwrite_cyl_tensor_smp(np_smp, numnod, inod_smp_stack,     &
     &    d_nod(1,i_sgs_simi), xx(1,1), xx(1,2), s, a_s)

      call overwrite_cyl_tensor_smp(np_smp, numnod, inod_smp_stack,     &
     &    d_nod(1,i_sgs_grad), xx(1,1), xx(1,2), s, a_s)

      call overwrite_cyl_tensor_smp(np_smp, numnod, inod_smp_stack,     &
     &    d_nod(1,i_sgs_grad_f), xx(1,1), xx(1,2), s, a_s)
!$omp end parallel
!
      end subroutine convert_dynamic_tensors_2_cyl
!
!  ---------------------------------------------------------------------
!
      subroutine convert_dynamic_tensors_2_sph(numnod,                  &
     &          inod_smp_stack, xx, r, s, a_r, a_s, ncomp_nod,          &
     &          i_sgs_simi, i_sgs_grad, i_sgs_grad_f, d_nod)
!
      use cvt_xyz_tensor_2_sph_smp
!
      integer (kind = kint), intent(in) :: numnod
      integer (kind = kint), intent(in) :: inod_smp_stack(0:np_smp)
      real(kind=kreal), intent(in) :: xx(numnod,3)
      real(kind=kreal), intent(in) :: r(numnod)
      real(kind=kreal), intent(in) :: s(numnod)
      real(kind=kreal), intent(in) :: a_r(numnod)
      real(kind=kreal), intent(in) :: a_s(numnod)
!
      integer (kind = kint), intent(in) :: ncomp_nod, i_sgs_simi 
      integer (kind = kint), intent(in) :: i_sgs_grad, i_sgs_grad_f
      real(kind=kreal), intent(inout) :: d_nod(numnod,ncomp_nod)
!
      if(iflag_debug .gt. 0) write(*,*) 'convert spherical corrdinate'
!$omp parallel
      call overwrite_sph_tensor_smp(np_smp, numnod, inod_smp_stack,     &
     &    d_nod(1,i_sgs_simi), xx(1,1), xx(1,2), xx(1,3),               &
     &    r, s, a_r, a_s)
!
      call overwrite_sph_tensor_smp(np_smp, numnod, inod_smp_stack,     &
     &    d_nod(1,i_sgs_grad), xx(1,1), xx(1,2), xx(1,3),               &
     &    r, s, a_r, a_s)
!
      call overwrite_sph_tensor_smp(np_smp, numnod, inod_smp_stack,     &
     &    d_nod(1,i_sgs_grad_f), xx(1,1), xx(1,2), xx(1,3),             &
     &    r, s, a_r, a_s)
!$omp end parallel
!
      end subroutine convert_dynamic_tensors_2_sph
!
!  ---------------------------------------------------------------------
!
      end module cvt_dynamic_scheme_coord
