!m_ctl_data_test_bc_temp.f90
!      module m_ctl_data_test_bc_temp
!
!      Written by H. Matsui on July, 2006
!      Mmodified by H. Matsui on June, 2007
!
!      subroutine read_control_4_bc_temp
!
!
!!   --------------------------------------------------------------------
!!    Example of control block
!!
!!  begin mesh_test
!!    begin data_files_def
!!      num_subdomain_ctl           2
!!      num_smp_ctl                 1
!!  
!!      mesh_file_prefix         'mesh/in'
!!      mesh_file_fmt_ctl        'gzip'
!!    end data_files_def
!!  
!!    begin boundary_ctl
!!      node_grp_name_ctl      'CMB'
!!      harmonics_degree_ctl      2
!!      harmonics_order_ctl      -2
!!    end boundary_ctl
!!  end mesh_test
!!
!!    -------------------------------------------------------------------
!
      module m_ctl_data_test_bc_temp
!
      use m_precision
      use t_ctl_data_4_platforms
      use t_control_elements
!
      implicit  none
!
!
      integer(kind = kint), parameter :: test_mest_ctl_file_code = 11
      character(len = kchara), parameter                                &
     &                        :: fname_test_mesh_ctl = "ctl_bc_temp"
!
!>      Structure for file settings
      type(platform_data_control), save :: bc_test_plt
!
      type(read_character_item), save :: temp_nod_grp_name
      type(read_integer_item), save :: hermonic_degree_ctl
      type(read_integer_item), save :: hermonic_order_ctl
!
!     Top level
!
      character(len=kchara), parameter                                  &
     &         :: hd_mesh_test_ctl = 'mesh_test'
      integer (kind=kint) :: i_mesh_test_ctl = 0
!
!     1st level
!
      character(len=kchara), parameter                                  &
     &                    :: hd_platform = 'data_files_def'
      character(len=kchara), parameter :: hd_bc_def =   'boundary_ctl'
!
      integer (kind=kint) :: i_platform =   0
      integer (kind=kint) :: i_bc_def =    0
!
!     2nd level for boundary defeine
!
      character(len=kchara), parameter                                  &
     &                      :: hd_nod_grp_t =  'node_grp_name_ctl'
      character(len=kchara), parameter                                  &
     &                      :: hd_sph_degree = 'harmonics_degree_ctl'
      character(len=kchara), parameter                                  &
     &                      :: hd_sph_order =  'harmonics_order_ctl'
!
      private :: test_mest_ctl_file_code, fname_test_mesh_ctl
!
      private :: hd_platform, i_platform
      private :: hd_bc_def, i_bc_def
      private :: read_test_mesh_ctl_data
      private :: hd_nod_grp_t, hd_sph_degree, hd_sph_order
!
      private :: read_ctl_data_4_temp_nod_bc
!
!   --------------------------------------------------------------------
!
      contains
!
!   --------------------------------------------------------------------
!
      subroutine read_control_4_bc_temp
!
      use m_machine_parameter
      use m_read_control_elements
      use skip_comment_f
!
!
      ctl_file_code = test_mest_ctl_file_code
!
      open(ctl_file_code, file = fname_test_mesh_ctl, status='old')
!
      call load_ctl_label_and_line
      call read_test_mesh_ctl_data
!
      close(ctl_file_code)
!
      end subroutine read_control_4_bc_temp
!
!  ---------------------------------------------------------------------
!   --------------------------------------------------------------------
!
      subroutine read_test_mesh_ctl_data
!
      use m_machine_parameter
      use m_read_control_elements
      use skip_comment_f
!
!
      if(right_begin_flag(hd_mesh_test_ctl) .eq. 0) return
      if (i_mesh_test_ctl .gt. 0) return
      do
        call load_ctl_label_and_line
!
        call find_control_end_flag(hd_mesh_test_ctl, i_mesh_test_ctl)
        if(i_mesh_test_ctl .gt. 0) exit
!
!
        call read_control_platforms                                     &
     &     (hd_platform, i_platform, bc_test_plt)
        call read_ctl_data_4_temp_nod_bc
      end do
!
      end subroutine read_test_mesh_ctl_data
!
!   --------------------------------------------------------------------
!
      subroutine read_ctl_data_4_temp_nod_bc
!
      use m_machine_parameter
      use m_read_control_elements
      use skip_comment_f
!
!
      if(right_begin_flag(hd_bc_def) .eq. 0) return
      if (i_bc_def .gt. 0) return
      do
        call load_ctl_label_and_line
!
        call find_control_end_flag(hd_bc_def, i_bc_def)
        if(i_bc_def .gt. 0) exit
!
!
        call read_chara_ctl_type(hd_nod_grp_t, temp_nod_grp_name)
!
        call read_integer_ctl_type(hd_sph_degree, hermonic_degree_ctl)
        call read_integer_ctl_type(hd_sph_order, hermonic_order_ctl)
      end do
!
      end subroutine read_ctl_data_4_temp_nod_bc
!
!   --------------------------------------------------------------------
!
      end module m_ctl_data_test_bc_temp
