!analyzer_MHD_lengthscale.f90
!      module analyzer_MHD_lengthscale
!
!      Written by H. Matsui on Dec., 2007
!
!
!      subroutine initialize_MHD_lscale
!      subroutine analyze_MHD_lscale
!
!..................................................
!
      module analyzer_MHD_lengthscale
!
      use m_precision
      use calypso_mpi
!
      use m_constants
      use m_machine_parameter
!
      use m_FEM_utils
      use t_VIZ_step_parameter
      use ucd_IO_select
      use FEM_MHD_length_scale
!
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine initialize_MHD_lscale
!
      use m_array_for_send_recv
      use m_ctl_params_4_prod_udt
      use m_ctl_data_product_udt
      use t_step_parameter
      use product_udt_fields
!
      integer(kind = kint) :: ierr
!
!
      if (my_rank.eq.0) then
        write(*,*) 'Get length scale'
        write(*,*) 'Input file: mesh data, udt data'
      end if
!
!     --------------------- 
!
      if (iflag_debug.eq.1) write(*,*) 'read_control_4_prod_udt'
      call read_control_4_prod_udt
!
      if (iflag_debug.eq.1) write(*,*) 'set_ctl_params_prod_udt'
      call set_ctl_params_prod_udt                                      &
     &   (mesh_file_FUTIL, udt_param_FUTIL)
      call set_fixed_time_step_params                                   &
     &   (t_pu_ctl, time_U, ierr, e_message)
!
!     ---------------------
!
      call mesh_setup_4_FEM_UTIL(mesh_file_FUTIL)
!
!     --------------------- 
!
      call allocate_work_4_lscale(femmesh_FUTIL%mesh%node%numnod)
      write(*,*) 'find_field_address_4_lscale'
      call find_field_address_4_lscale(field_FUTIL, iphys_FUTIL)
!
      end subroutine initialize_MHD_lscale
!
! ----------------------------------------------------------------------
!
      subroutine analyze_MHD_lscale
!
      use m_ctl_params_4_diff_udt
      use set_ucd_data_to_type
      use FEM_MHD_length_scale
!
      integer(kind=kint ) :: istep
!
!
      do istep = time_U%init_d%i_time_step, time_U%finish_d%i_end_step
        if(set_IO_step_flag(istep,time_U%ucd_step) .ne. izero) cycle
!
        call set_data_by_read_ucd_once                                  &
     &     (my_rank, time_U%ucd_step%istep_file,                        &
     &      first_ucd_param, field_FUTIL, time_IO_FUTIL)
!
        call const_MHD_length_scales                                    &
     &     (femmesh_FUTIL%mesh%node, iphys_FUTIL, field_FUTIL,          &
     &      time_U%ucd_step%istep_file, time_IO_FUTIL, ucd_FUTIL)
      end do
!
      call deallocate_work_4_lscale
      call calypso_MPI_barrier
!
      end subroutine analyze_MHD_lscale
!
! ----------------------------------------------------------------------
!
      end module analyzer_MHD_lengthscale

