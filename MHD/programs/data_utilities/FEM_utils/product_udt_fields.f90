!
!      module product_udt_fields
!
!      Written by H. Matsui on Dec., 2007
!
!!      subroutine allocate_product_data(numnod)
!!      subroutine allocate_product_result(nod_fld)
!!        type(phys_data), intent(inout) :: nod_fld
!!      subroutine deallocate_product_data
!!
!!      subroutine set_field_id_4_product                               &
!!     &         (init_d, numnod, t_IO, ucd_step)
!!        type(time_data), intent(in) :: init_d
!!        type(time_data), intent(inout) :: t_IO
!!        type(IO_step_param), intent(inout) :: ucd_step
!!
!!      subroutine set_data_for_product(numnod, istep_ucd, t_IO)
!!      subroutine cal_rev_of_2nd_field(numnod)
!!      subroutine cal_products_of_fields                               &
!!     &         (nod_comm, node, ncomp_nod, d_nod)
!
      module product_udt_fields
!
      use m_precision
      use m_constants
      use m_machine_parameter
!
      use t_IO_step_parameter
      use t_time_data
      use t_comm_table
      use t_geometry_data
      use t_phys_data
!
      implicit none
!
      real(kind = kreal), allocatable :: d_prod1(:,:)
      real(kind = kreal), allocatable :: d_prod2(:,:)
!
      private :: d_prod1, d_prod2
      private :: find_field_id_in_read_ucd
      private :: set_one_field_by_read_ucd_once
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine allocate_product_data(numnod)
!
      use m_ctl_params_4_prod_udt
!
      integer(kind = kint), intent(in) :: numnod
!
!
      allocate(d_prod1(numnod,ncomp_4_product1))
      allocate(d_prod2(numnod,ncomp_4_product2))
      d_prod1 = 0.0d0
      d_prod2 = 0.0d0
!
      end subroutine allocate_product_data
!
!-----------------------------------------------------------------------
!
      subroutine deallocate_product_data
!
!
      deallocate(d_prod1, d_prod2)
!
      end subroutine deallocate_product_data
!
!-----------------------------------------------------------------------
!
      subroutine allocate_product_result(nod_fld)
!
      use m_ctl_params_4_prod_udt
!
      type(phys_data), intent(inout) :: nod_fld
!
!
      nod_fld%num_phys =     ione
      nod_fld%num_phys_viz = ione
      call alloc_phys_name_type(nod_fld)
!
      nod_fld%num_component(1) =    ncomp_4_result
      nod_fld%istack_component(1) = ncomp_4_result
      nod_fld%iflag_monitor(1) =    1
      nod_fld%phys_name(1) =        result_field_name
!
      nod_fld%ntot_phys =     ncomp_4_result
      nod_fld%ntot_phys_viz = ncomp_4_result
!
      call alloc_phys_data_type(nod_fld%n_point, nod_fld)
!
      end subroutine allocate_product_result
!
!-----------------------------------------------------------------------
!
      subroutine set_field_id_4_product                                 &
     &         (init_d, numnod, t_IO, ucd_step)
!
      use calypso_mpi
      use t_time_data
      use m_error_IDs
      use m_ctl_params_4_prod_udt
      use ucd_IO_select
!
      integer(kind = kint), intent(in) :: numnod
      type(time_data), intent(in) :: init_d
      type(time_data), intent(inout) :: t_IO
      type(IO_step_param), intent(inout) :: ucd_step
!
!
      ucd_step%istep_file = init_d%i_time_step / ucd_step%increment
      call find_field_id_in_read_ucd                                    &
     &   (my_rank, ucd_step%istep_file, prod1_ucd_param,                &
     &    numnod, product_field_1_name, i_field_product1,               &
     &    ncomp_4_product1, t_IO)
!
      call find_field_id_in_read_ucd                                    &
     &   (my_rank, ucd_step%istep_file, prod2_ucd_param,                &
     &    numnod, product_field_2_name, i_field_product2,               &
     &    ncomp_4_product2, t_IO)
!
      if( (i_field_product1*i_field_product2) .eq. 0) then
        call calypso_MPI_abort(ierr_fld,'Field does not excist')
      end if
!
      if(ncomp_4_product1.eq.1) then
        if(ncomp_4_product2 .eq. 1) then
          ncomp_4_result = 1
        else if(ncomp_4_product2 .eq.3) then
          ncomp_4_result = 3
        else if(ncomp_4_product2 .eq.6) then
          ncomp_4_result = 6
        end if
      else if(ncomp_4_product1.eq.3) then
        if(ncomp_4_product2 .eq. 1) then
          ncomp_4_result = 3
        else if(ncomp_4_product2 .eq.3) then
          if(iflag_product_type .eq. 2) then
            ncomp_4_result = 3
          else
            ncomp_4_result = 1
          end if
        else if(ncomp_4_product2 .eq.6) then
            ncomp_4_result = 3
        end if
      else if(ncomp_4_product1.eq.6) then
        if(ncomp_4_product2 .eq. 1) then
          ncomp_4_result = 6
        else if(ncomp_4_product2 .eq.3) then
          ncomp_4_result = 3
        else if(ncomp_4_product2 .eq.6) then
          ncomp_4_result = 6
        end if
      end if
!
      end subroutine set_field_id_4_product
!
!-----------------------------------------------------------------------
!
      subroutine set_data_for_product(numnod, istep_ucd, t_IO)
!
      use calypso_mpi
      use m_ctl_params_4_prod_udt
!
      integer(kind = kint), intent(in) :: numnod
      integer(kind = kint), intent(in) :: istep_ucd
      type(time_data), intent(inout) :: t_IO
!
!
      call set_one_field_by_read_ucd_once                               &
     &   (my_rank, istep_ucd, prod1_ucd_param,                          &
     &    i_field_product1, ncomp_4_product1, numnod, d_prod1, t_IO)
!
      call set_one_field_by_read_ucd_once                               &
     &   (my_rank, istep_ucd, prod2_ucd_param,                          &
     &    i_field_product2, ncomp_4_product2, numnod, d_prod2, t_IO)
!
      end subroutine set_data_for_product
!
!-----------------------------------------------------------------------
!
      subroutine cal_rev_of_2nd_field(numnod)
!
      use m_ctl_params_4_prod_udt
!
      integer(kind = kint), intent(in) :: numnod
      integer(kind = kint) :: nd, inod
!
!
       do nd = 1, ncomp_4_product2
         do inod = 1, numnod
           if(d_prod2(inod,nd) .eq. zero) then
             d_prod2(inod,nd) = zero
           else
             d_prod2(inod,nd) = one / d_prod2(inod,nd)
           end if
         end do
       end do
!
!
      end subroutine cal_rev_of_2nd_field
!
!-----------------------------------------------------------------------
!
      subroutine cal_products_of_fields                                 &
     &         (nod_comm, node, ncomp_nod, d_nod)
!
      use m_ctl_params_4_prod_udt
      use cal_products_smp
      use nod_phys_send_recv
!
      type(node_data), intent(in) :: node
      type(communication_table), intent(in) :: nod_comm
      integer(kind = kint), intent(in) :: ncomp_nod
      real(kind = kreal), intent(inout) :: d_nod(node%numnod,ncomp_nod)
!
!
      if(ncomp_4_product1.eq.1) then
        if(ncomp_4_product2 .eq. 1) then
!$omp parallel
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,1), d_prod2(1,1),  d_nod(1,1))
!$omp end parallel
        else if(ncomp_4_product2 .eq.3) then
!$omp parallel
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,1), d_prod2(1,1),  d_nod(1,1))
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,1), d_prod2(1,2),  d_nod(1,2))
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,1), d_prod2(1,3),  d_nod(1,3))
!$omp end parallel
        else if(ncomp_4_product2 .eq.6) then
!$omp parallel
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,1), d_prod2(1,1),  d_nod(1,1))
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,1), d_prod2(1,2),  d_nod(1,2))
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,1), d_prod2(1,3),  d_nod(1,3))
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,1), d_prod2(1,4),  d_nod(1,4))
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,1), d_prod2(1,5),  d_nod(1,5))
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,1), d_prod2(1,6),  d_nod(1,6))
!$omp end parallel
        end if
      else if(ncomp_4_product1.eq.3) then
        if(ncomp_4_product2 .eq. 1) then
!$omp parallel
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,1), d_prod2(1,1),  d_nod(1,1))
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,2), d_prod2(1,1),  d_nod(1,2))
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,3), d_prod2(1,1),  d_nod(1,3))
!$omp end parallel
        else if(ncomp_4_product2 .eq.3) then
          if(iflag_product_type .eq. 2) then
!$omp parallel
           call cal_cross_prod_no_coef_smp(node%numnod,                 &
     &         d_prod1(1,1), d_prod2(1,1), d_nod(1,1))
!$omp end parallel
          else
!$omp parallel
           call cal_dot_prod_no_coef_smp(node%numnod,                   &
     &         d_prod1(1,1), d_prod2(1,1),  d_nod(1,1))
!$omp end parallel
          end if
        else if(ncomp_4_product2 .eq.6) then
!$omp parallel
           call cal_tensor_vec_prod_no_coef_smp(node%numnod,            &
     &         d_prod2(1,1), d_prod1(1,1),  d_nod(1,1))
!$omp end parallel
        end if
      else if(ncomp_4_product1.eq.6) then
        if(ncomp_4_product2 .eq. 1) then
!$omp parallel
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,1), d_prod2(1,1),  d_nod(1,1))
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,2), d_prod2(1,1),  d_nod(1,2))
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,3), d_prod2(1,1),  d_nod(1,3))
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,4), d_prod2(1,1),  d_nod(1,4))
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,5), d_prod2(1,1),  d_nod(1,5))
           call cal_scalar_prod_no_coef_smp(node%numnod,                &
     &         d_prod1(1,6), d_prod2(1,1),  d_nod(1,6))
!$omp end parallel
        else if(ncomp_4_product2 .eq.3) then
!$omp parallel
           call cal_tensor_vec_prod_no_coef_smp(node%numnod,            &
     &         d_prod1(1,1), d_prod2(1,1),  d_nod(1,1))
!$omp end parallel
        else if(ncomp_4_product2 .eq.6) then
!$omp parallel
           call cal_tensor_vec_prod_no_coef_smp(node%numnod,            &
     &         d_prod1(1,1), d_prod2(1,1),  d_nod(1,1))
           call cal_tensor_vec_prod_no_coef_smp(node%numnod,            &
     &         d_prod1(1,1), d_prod2(1,4),  d_nod(1,4))
!$omp end parallel
        end if
      end if
!
!
      if(ncomp_4_result .eq. ione) then
        call nod_scalar_send_recv(node%numnod, nod_comm, d_nod)
      else if(ncomp_4_result .eq. ithree) then
        call nod_vector_send_recv(node%numnod, nod_comm, d_nod)
      else if(ncomp_4_result .eq. isix) then
        call nod_tensor_send_recv(node%numnod, nod_comm, d_nod)
      end if
!
      end subroutine cal_products_of_fields
!
!-----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine find_field_id_in_read_ucd(my_rank, istep_ucd,          &
     &          ucd_param, numnod, field_name,                          &
     &          i_field, ncomp_field, t_IO)
!
      use set_and_cal_udt_data
      use ucd_IO_select
!
      integer(kind = kint),  intent(in) :: my_rank, istep_ucd, numnod
      type(field_IO_params), intent(in) :: ucd_param
!
      character(len = kchara), intent(in) :: field_name
      integer(kind = kint),  intent(inout) :: i_field, ncomp_field
      type(time_data), intent(inout) :: t_IO
!
      type(ucd_data) :: local_ucd
!
!
      local_ucd%nnod = numnod
      call sel_read_udt_param                                           &
     &   (my_rank, istep_ucd, ucd_param, t_IO, local_ucd)
      call find_field_id_in_ucd(local_ucd, field_name,                  &
     &    i_field, ncomp_field)
      call deallocate_ucd_data(local_ucd)
!
      end subroutine find_field_id_in_read_ucd
!
! -----------------------------------------------------------------------
!
      subroutine set_one_field_by_read_ucd_once(my_rank, istep_ucd,     &
     &          ucd_param, i_field, ncomp_field, numnod, d_fld, t_IO)
!
      use set_and_cal_udt_data
      use ucd_IO_select
!
      type(field_IO_params), intent(in) :: ucd_param
      integer(kind = kint),  intent(in) :: my_rank, istep_ucd
      integer(kind = kint),  intent(in) :: numnod, i_field, ncomp_field
!
      real(kind = kreal), intent(inout) :: d_fld(numnod,ncomp_field)
      type(time_data), intent(inout) :: t_IO
!
      type(ucd_data) :: local_ucd
!
!
!
      local_ucd%nnod =      numnod
      call sel_read_alloc_udt_file                                      &
     &   (my_rank, istep_ucd, ucd_param, t_IO, local_ucd)
      call set_one_field_by_udt_data(numnod, ncomp_field,               &
     &    i_field, d_fld, local_ucd)
      call deallocate_ucd_data(local_ucd)
!
      end subroutine set_one_field_by_read_ucd_once
!
! -----------------------------------------------------------------------
!
      end module product_udt_fields
