!collect_psf_data
!      module collect_psf_data
!
!      Written by H. Matsui on July, 2006
!
!!      subroutine collect_numbers_4_psf(num_psf, psf_prefix, ifmt_psf, &
!!     &     istack_nod_psf_smp, istack_patch_psf_smp,                  &
!!     &     psf_fld, collect, psf_ucd)
!
!      subroutine collect_mesh_4_psf
!      subroutine collect_field_4_psf
!
!      subroutine collect_numbers_4_iso
!      subroutine collect_data_4_iso
!
      module collect_psf_data
!
      use m_precision
!
      use m_constants
      use calypso_mpi
      use m_machine_parameter
!
      implicit  none
!
!  ---------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine collect_numbers_4_psf(num_psf, psf_prefix, ifmt_psf,   &
     &     istack_nod_psf_smp, istack_patch_psf_smp,                    &
     &     psf_fld, collect, psf_ucd)
!
      use m_geometry_constants
      use t_psf_outputs
      use t_phys_data
      use t_ucd_data
!
      use count_numbers_collected_psf
!
      integer(kind = kint), intent(in) :: num_psf
      integer(kind = kint), intent(in) :: ifmt_psf(num_psf)
      character(len = kchara), intent(in) :: psf_prefix(num_psf)
!
      integer(kind = kint), intent(in)                                  &
     &      :: istack_nod_psf_smp(0:np_smp*num_psf)
      integer(kind = kint), intent(in)                                  &
     &      :: istack_patch_psf_smp(0:np_smp*num_psf)
!
      type(phys_data), intent(in) :: psf_fld(num_psf)
      type(psf_collect_type), intent(inout) :: collect
      type(ucd_data), intent(inout) :: psf_ucd(num_psf)
!
      integer(kind = kint) :: i_psf, i_fld, nmax
!
!
      call count_numbers_4_psf_out(num_psf, istack_nod_psf_smp,         &
     &    collect%ntot_nod_output_psf, nmax,                            &
     &    collect%istack_nod_para_psf, collect%istack_nod_recv_psf,     &
     &    collect%istack_nod_output_psf)
!
      call count_numbers_4_psf_out(num_psf, istack_patch_psf_smp,       &
     &    collect%ntot_ele_output_psf, nmax,                            &
     &    collect%istack_ele_para_psf, collect%istack_ele_recv_psf,     &
     &    collect%istack_ele_output_psf)
!
      if(my_rank .gt. 0) return
!
      psf_ucd(1:num_psf)%file_prefix = psf_prefix(1:num_psf)
      psf_ucd(1:num_psf)%ifmt_file = ifmt_psf(1:num_psf)
      psf_ucd(1:num_psf)%nnod_4_ele = num_triangle
!
      do i_psf = 1, num_psf
        psf_ucd(i_psf)%nnod = collect%istack_nod_output_psf(i_psf)      &
     &                       - collect%istack_nod_output_psf(i_psf-1)
        psf_ucd(i_psf)%nele = collect%istack_ele_output_psf(i_psf)      &
     &                       - collect%istack_ele_output_psf(i_psf-1)
        psf_ucd(i_psf)%num_field = psf_fld(i_psf)%num_phys
        psf_ucd(i_psf)%ntot_comp = psf_fld(i_psf)%ntot_phys
!
        call allocate_ucd_phys_name( psf_ucd(i_psf) )
!
        do i_fld = 1, psf_ucd(i_psf)%num_field
          psf_ucd(i_psf)%phys_name(i_fld)                               &
     &         =  psf_fld(i_psf)%phys_name(i_fld)
          psf_ucd(i_psf)%num_comp(i_fld)                                &
     &         =  psf_fld(i_psf)%num_component(i_fld)
        end do
!
        call allocate_ucd_node(psf_ucd(i_psf))
        call allocate_ucd_ele(psf_ucd(i_psf))
        call allocate_ucd_phys_data(psf_ucd(i_psf))
      end do
!
!
      end subroutine collect_numbers_4_psf
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine collect_mesh_4_psf(num_psf, patch, collect, psf_ucd)
!
      use t_psf_geometry_list
      use t_psf_outputs
      use t_ucd_data
!
      use psf_send_recv
      use reconnect_psf_overlap_nod
!
      integer(kind = kint), intent(in) :: num_psf
      type(psf_patch_data), intent(in) ::   patch
      type(psf_collect_type), intent(in) :: collect
      type(ucd_data), intent(inout) :: psf_ucd(num_psf)
!
!
      call psf_grids_send_recv(num_psf, patch%nnod_psf_tot,             &
     &    collect%ntot_nod_output_psf, collect%istack_nod_para_psf,     &
     &    collect%istack_nod_recv_psf, patch%xyz_psf, collect%send_psf, &
     &    collect%recv_psf, psf_ucd)
!
      call psf_hash_send_recv(num_psf, patch%nnod_psf_tot,              &
     &    collect%ntot_nod_output_psf, collect%istack_nod_para_psf,     &
     &    collect%istack_nod_recv_psf, patch%inod_hash_psf,             &
     &    collect%isend_psf, collect%irecv_psf,                         &
     &    collect%ihash_output_psf)
!
      call set_global_psf_node_id(num_psf, psf_ucd)
!
      call psf_connect_send_recv(num_psf, patch%npatch_tot,             &
     &    collect%ntot_ele_output_psf, collect%istack_nod_para_psf,     &
     &    collect%istack_ele_para_psf, collect%istack_ele_recv_psf,     &
     &    patch%ie_tri, collect%isend_psf, collect%irecv_psf, psf_ucd)
!
      call s_reconnect_psf_overlap_nod(num_psf,                         &
     &    collect%ntot_nod_output_psf, collect%istack_nod_output_psf,   &
     &    collect%ihash_output_psf, psf_ucd)
!
      end subroutine collect_mesh_4_psf
!
! ----------------------------------------------------------------------
!
      subroutine collect_field_4_psf(patch, collect, psf_ucd)
!
      use m_control_params_4_psf
      use t_psf_geometry_list
      use t_psf_outputs
      use t_ucd_data
!
      use psf_send_recv
!
      type(psf_patch_data), intent(in) ::   patch
      type(psf_collect_type), intent(in) :: collect
      type(ucd_data), intent(inout) :: psf_ucd(num_psf)
!
!
      call psf_results_send_recv(num_psf, patch%nnod_psf_tot,           &
     &    collect%ntot_nod_output_psf, collect%istack_nod_para_psf,     &
     &    collect%istack_nod_recv_psf, max_ncomp_psf_out,               &
     &    patch%dat_psf, collect%send_psf, collect%recv_psf, psf_ucd)
!
      end subroutine collect_field_4_psf
!
! ----------------------------------------------------------------------
!
      subroutine deallocate_psf_outputs_data(my_rank, num_psf)
!
      use m_psf_data
!
      integer(kind = kint), intent(in) :: my_rank, num_psf
      integer(kind = kint) :: i_psf
!
!
      if(my_rank .ne. 0) return
      do i_psf = 1, num_psf
        call deallocate_ucd_mesh(psf_out(i_psf))
      end do
!
      end subroutine deallocate_psf_outputs_data
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine collect_data_4_iso(patch, collect, psf_ucd)
!
      use m_control_params_4_iso
      use t_psf_geometry_list
      use t_psf_outputs
      use t_ucd_data
!
      use psf_send_recv
      use reconnect_psf_overlap_nod
!
      type(psf_patch_data), intent(in) ::   patch
      type(psf_collect_type), intent(in) :: collect
      type(ucd_data), intent(inout) :: psf_ucd(num_iso)
!
!
      call psf_grids_send_recv(num_iso, patch%nnod_psf_tot,           &
     &    collect%ntot_nod_output_psf, collect%istack_nod_para_psf,             &
     &    collect%istack_nod_recv_psf, patch%xyz_psf,                         &
     &    collect%send_psf, collect%recv_psf, psf_ucd)
!
      call psf_hash_send_recv(num_iso, patch%nnod_psf_tot,            &
     &    collect%ntot_nod_output_psf, collect%istack_nod_para_psf,             &
     &    collect%istack_nod_recv_psf, patch%inod_hash_psf,                   &
     &    collect%isend_psf, collect%irecv_psf,                         &
     &    collect%ihash_output_psf)
!
      call set_global_psf_node_id(num_iso, psf_ucd)
!
      call psf_connect_send_recv(num_iso, patch%npatch_tot,           &
     &    collect%ntot_ele_output_psf, collect%istack_nod_para_psf,             &
     &    collect%istack_ele_para_psf, collect%istack_ele_recv_psf, patch%ie_tri,     &
     &    collect%isend_psf, collect%irecv_psf, psf_ucd)
!
      call s_reconnect_psf_overlap_nod(num_iso, collect%ntot_nod_output_psf,    &
     &    collect%istack_nod_output_psf, collect%ihash_output_psf, psf_ucd)
!
      call psf_results_send_recv(num_iso, patch%nnod_psf_tot,         &
     &    collect%ntot_nod_output_psf, collect%istack_nod_para_psf,             &
     &    collect%istack_nod_recv_psf, max_ncomp_iso_out, patch%dat_psf,      &
     &    collect%send_psf, collect%recv_psf, psf_ucd)
!
      end subroutine collect_data_4_iso
!
! ----------------------------------------------------------------------
!
      subroutine deallocate_iso_outputs_data(my_rank, num_iso)
!
      use m_iso_data
!
      integer(kind = kint), intent(in) :: my_rank, num_iso
      integer(kind = kint) :: i_iso
!
!
      if(my_rank .ne. 0) return
      do i_iso = 1, num_iso
        call deallocate_ucd_mesh(iso_out(i_iso))
      end do
!
      end subroutine deallocate_iso_outputs_data
!
! ----------------------------------------------------------------------
!
      end module collect_psf_data
