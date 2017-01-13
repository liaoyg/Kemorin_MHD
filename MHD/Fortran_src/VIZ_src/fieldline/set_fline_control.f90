!set_fline_control.f90
!      module set_fline_control
!
!     Written by H. Matsui on Aug., 2011
!
!!      subroutine s_set_fline_control(numele, interior_ele,            &
!!     &       num_mat, num_mat_bc, mat_name, mat_istack, mat_item,     &
!!     &       num_surf, num_surf_bc, surf_name, surf_istack, surf_item,&
!!     &       num_nod_phys, phys_nod_name)
!
      module set_fline_control
!
      use m_precision
!
      use m_machine_parameter
      use m_control_data_flines
      use m_control_params_4_fline
      use m_source_4_filed_line
!
      implicit none
!
      character(len=kchara) :: hd_fline_ctl = 'field_line'
      private :: hd_fline_ctl
!
      private :: read_control_4_fline
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine s_set_fline_control(numele, interior_ele,              &
     &       num_mat, num_mat_bc, mat_name, mat_istack, mat_item,       &
     &       num_surf, num_surf_bc, surf_name, surf_istack, surf_item,  &
     &       num_nod_phys, phys_nod_name)
!
      use set_control_each_fline
!
      integer(kind=kint), intent(in) :: numele
      integer(kind=kint), intent(in) :: interior_ele(numele)
!
      integer(kind=kint), intent(in) :: num_mat, num_mat_bc
      integer(kind=kint), intent(in) :: mat_istack(0:num_mat)
      integer(kind=kint), intent(in) :: mat_item(num_mat_bc)
      character(len=kchara), intent(in) :: mat_name(num_mat)
!
      integer(kind=kint), intent(in) :: num_surf, num_surf_bc
      integer(kind=kint), intent(in) :: surf_istack(0:num_surf)
      integer(kind=kint), intent(in) :: surf_item(2,num_surf_bc)
      character(len=kchara), intent(in) :: surf_name(num_surf)
!
      integer(kind = kint), intent(in) :: num_nod_phys
      character(len=kchara), intent(in) :: phys_nod_name(num_nod_phys)
!
      integer(kind = kint) :: i_fln
!
!
      ctl_file_code = fline_ctl_file_code
!
      call allocate_control_params_fline
      call allocate_local_start_grp_num
!
      do i_fln = 1, num_fline
        call read_control_4_fline(i_fln)
      end do
!
      do i_fln = 1, num_fline
        call count_control_4_fline(i_fln, fline_ctl_struct(i_fln),      &
     &      numele, interior_ele, num_mat, mat_name,                    &
     &      num_surf, num_surf_bc, surf_name, surf_istack, surf_item)
      end do
!
      call allocate_iflag_fline_used_ele(numele)
      call allocate_fline_starts_ctl
      call allocate_local_start_grp_item
!
      do i_fln = 1, num_fline
        call set_control_4_fline(i_fln, fline_ctl_struct(i_fln),        &
     &      numele, interior_ele, num_mat, mat_name,                    &
     &      num_surf, num_surf_bc, surf_istack, surf_item,              &
     &      num_nod_phys, phys_nod_name)
        call set_iflag_fline_used_ele(i_fln, numele, interior_ele,      &
     &      num_mat, num_mat_bc, mat_istack, mat_item)
        call deallocate_cont_dat_fline(fline_ctl_struct(i_fln))
!
        if(iflag_debug .gt. 0) call check_control_params_fline(i_fln)
      end do
!
      call deallocate_fline_fhead_ctl
!
      end subroutine s_set_fline_control
!
!   --------------------------------------------------------------------
!   --------------------------------------------------------------------
!
      subroutine read_control_4_fline(i_fline)
!
      use calypso_mpi
      use t_control_data_4_fline
!
      integer(kind = kint), intent(in) :: i_fline
!
      if(fname_fline_ctl(i_fline) .eq. 'NO_FILE') return
      call reset_fline_control_flags(fline_ctl_struct(i_fline))
!
      if(my_rank .eq. 0) then
        open(fline_ctl_file_code, file=fname_fline_ctl(i_fline),        &
     &         status='old')
        call load_ctl_label_and_line
        call read_field_line_ctl                                        &
     &     (hd_fline_ctl, fline_ctl_struct(i_fline))
        close(fline_ctl_file_code)
      end if
!
      call bcast_field_line_ctl(fline_ctl_struct(i_fline))
!
      end subroutine read_control_4_fline
!
!  ---------------------------------------------------------------------
!
      end module set_fline_control
