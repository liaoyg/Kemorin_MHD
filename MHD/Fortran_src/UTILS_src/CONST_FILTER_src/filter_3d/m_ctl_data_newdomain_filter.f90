!m_ctl_data_newdomain_filter.f90
!      module m_ctl_data_newdomain_filter
!
!      subroutine read_control_filter_newdomain
!
!      begin org_filter_filtes_ctl
!        org_filter_file_header       'org/filter_node'
!        org_filter_elength_header    'org/filter_elength'
!        org_filter_moment_header     'org/filter_moms'
!        org_filter_coefs_header      'org/filter_coef'
!      end
!
!      Written by H. Matsui on Apr., 2008
!
      module m_ctl_data_newdomain_filter
!
      use m_precision
      use t_ctl_data_4_platforms
      use t_ctl_data_filter_files
!
      implicit  none
!
      integer(kind = kint), parameter :: id_filter_ctl_file = 11
      character(len = kchara), parameter                                &
     &             :: fname_trans_flt_ctl = "ctl_new_domain_filter"
!
      type(platform_data_control), save :: org_filter_plt
      type(platform_data_control), save :: new_filter_plt
!
!>      Structure for filtering files
      type(filter_file_control), save :: ffile_ndom_ctl
!
!     Top level
!
      character(len=kchara), parameter                                  &
     &         :: hd_filter_newdomain_ctl = 'change_filter_domain_ctl'
      integer (kind=kint) :: i_filter_newdomain_ctl = 0
!
!
      private :: id_filter_ctl_file, fname_trans_flt_ctl
      private :: hd_filter_newdomain_ctl, i_filter_newdomain_ctl
      private :: read_ctl_filter_newdomain_data
!
      character(len=kchara), parameter                                  &
     &                    :: hd_platform = 'data_files_def'
!
      character(len=kchara), parameter                                  &
     &                    :: hd_new_data = 'new_data_files_def'
      character(len=kchara), parameter :: hd_filter_fnames              &
     &                        = 'filter_files_def'
      integer (kind=kint) :: i_platform =   0
      integer (kind=kint) :: i_new_data =      0
      integer (kind=kint) :: i_filter_fnames = 0
!
      private :: hd_platform, i_platform
      private :: hd_new_data, i_new_data
      private :: hd_filter_fnames, i_filter_fnames
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine read_control_filter_newdomain
!
      use m_machine_parameter
      use m_read_control_elements
      use skip_comment_f
!
      integer(kind = kint) :: iflag
!
!
      ctl_file_code = id_filter_ctl_file
!
      open(ctl_file_code, file=fname_trans_flt_ctl, status='old')
!
      call load_ctl_label_and_line
      call read_ctl_filter_newdomain_data
!
      close(ctl_file_code)
!
      end subroutine read_control_filter_newdomain
!
!   --------------------------------------------------------------------
!
      subroutine read_ctl_filter_newdomain_data
!
      use m_machine_parameter
      use m_read_control_elements
      use skip_comment_f
!
      use m_ctl_data_org_filter_name
!
!
      if(right_begin_flag(hd_filter_newdomain_ctl) .eq. 0) return
      if (i_filter_newdomain_ctl .gt. 0) return
      do
        call load_ctl_label_and_line
!
        call find_control_end_flag(hd_filter_newdomain_ctl,             &
     &      i_filter_newdomain_ctl)
        if(i_filter_newdomain_ctl .gt. 0) exit
!
!
        call read_control_platforms                                     &
     &     (hd_platform, i_platform, org_filter_plt)
        call read_control_platforms                                     &
     &     (hd_new_data, i_new_data, new_filter_plt)
        call read_filter_fnames_control                                 &
     &     (hd_filter_fnames, i_filter_fnames, ffile_ndom_ctl)
        call read_org_filter_fnames_ctl
      end do
!
      end subroutine read_ctl_filter_newdomain_data
!
!   --------------------------------------------------------------------
!
      end module m_ctl_data_newdomain_filter
