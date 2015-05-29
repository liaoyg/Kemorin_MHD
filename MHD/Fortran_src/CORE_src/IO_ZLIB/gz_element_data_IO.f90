!
!      module gz_element_data_IO
!
!     Written by H. Matsui on Oct., 2006
!
!      subroutine write_element_comm_table_gz
!      subroutine read_element_comm_table_gz
!
!      subroutine write_element_geometry_gz
!      subroutine write_element_geometry_sph_gz
!      subroutine write_element_geometry_cyl_gz
!      subroutine read_element_geometries_gz
!
      module gz_element_data_IO
!
      use m_precision
!
      use skip_gz_comment
!
      implicit none
!
!------------------------------------------------------------------
!
       contains
!
!------------------------------------------------------------------
!
      subroutine write_element_comm_table_gz
!
      use gz_domain_data_IO
!
!
      write(textbuf,'(a,a1)') '!' , char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '!  element position ', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '!  and communication table ', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '! 1.parallel information', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
!
      call write_domain_info_gz
!
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '! 2.import / export information ',       &
     &      char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '! 2.1 element ID for import ', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
!
      call write_import_data_gz
!
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '! 2.2 element ID for export ', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '! ', char(0)
      call gz_write_textbuf_w_lf
!
      call write_export_data_gz
!
      end subroutine write_element_comm_table_gz
!
!------------------------------------------------------------------
!
      subroutine read_element_comm_table_gz
!
      use gz_domain_data_IO
!
!
!      write(id_file,*) '! 1.parallel information'
      call read_domain_info_gz
!
!      write(id_file,*) '! 2.import / export information '
!      write(id_file,*) '! 2.1 element ID for import '
      call read_import_data_gz
!      write(id_file,*) '! 2.2 element ID for export '
      call read_export_data_gz
!
      end subroutine read_element_comm_table_gz
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
      subroutine write_element_geometry_gz
!
      use gz_node_geometry_IO
!
!
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '! 3.element information', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '! 3.1 center of element (position) ',    &
     &      char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
!
      call write_geometry_info_gz
!
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '! 3.2 Volume of element ', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
!
      call write_scalar_in_element_gz
!
      end subroutine write_element_geometry_gz
!
!------------------------------------------------------------------
!
      subroutine write_element_geometry_sph_gz
!
      use gz_node_geometry_IO
!
!
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '! 3.element information', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '! 3.1 center of element (r,theta,phi)'   &
     &      , char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
!
      call write_geometry_info_gz
!
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '! 3.2 Volume of element ', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
!
      call write_scalar_in_element_gz
!
      end subroutine write_element_geometry_sph_gz
!
!------------------------------------------------------------------
!
      subroutine write_element_geometry_cyl_gz
!
      use gz_node_geometry_IO
!
!
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '! 3.element information', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '! 3.1 center of element (r,theta,phi)'   &
     &      , char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
!
      call write_geometry_info_gz
!
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '! 3.2 Volume of element ', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)') '!', char(0)
      call gz_write_textbuf_w_lf
!
      call write_scalar_in_element_gz
!
      end subroutine write_element_geometry_cyl_gz
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
      subroutine read_element_geometries_gz
!
      use gz_node_geometry_IO
!
!
!      write(id_file,*) '! 3.element information'
!      write(id_file,*) '! 3.1 center of element'
      call read_number_of_node_gz
      call read_geometry_info_gz
!
!      write(id_file,*) '! 3.2 Volume of element '
      call read_scalar_in_element_gz
!
      end subroutine read_element_geometries_gz
!
!------------------------------------------------------------------
!
      end module gz_element_data_IO
