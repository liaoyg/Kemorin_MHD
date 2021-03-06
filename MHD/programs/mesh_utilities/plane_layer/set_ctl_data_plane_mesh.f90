!
!      module set_ctl_data_plane_mesh
!
!        programmed by H.Matsui on Aug., 2007
!
!      subroutine s_set_ctl_data_plane_mesh
!
      module set_ctl_data_plane_mesh
!
      use m_precision
      use m_constants
!
      implicit  none
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine s_set_ctl_data_plane_mesh
!
      use m_size_4_plane
      use m_size_of_cube
      use m_cube_position
      use m_grp_data_cub_kemo
      use m_filtering_nod_4_cubmesh
      use m_cube_files_data
      use m_ctl_data_4_plane_model
      use m_ctl_data_4_cub_kemo
      use m_spheric_constants
      use skip_comment_f
!
!
      if (cubmesh_plt%mesh_file_prefix%iflag .gt. 0) then
        mesh_file_header = cubmesh_plt%mesh_file_prefix%charavalue
      else
        mesh_file_header = 'mesh/in'
      end if
!
      if (ffile_cub_ctl%filter_head_ctl%iflag .eq. 1) then
        filter_file_header = ffile_cub_ctl%filter_head_ctl%charavalue
      else
        filter_file_header = 'mesh/filter_node_l'
      end if
!
      if (z_filter_head_ctl%iflag .eq. 1) then
        z_filter_header = z_filter_head_ctl%charavalue
      else
        z_filter_header = 'filter_info'
      end if
!
!
!
      if (plane_size_ctl%iflag .eq. 0) then
        plane_size_ctl%realvalue(1:3) = 1.0d0
      end if
!
      pi = 4.0d0 * atan(1.0d0)
      if (cmp_no_case(unit_len_plane_ctl%charavalue(1),'pi')) then
        xsize = pi*plane_size_ctl%realvalue(1)
      else
        xsize = plane_size_ctl%realvalue(1)
      end if
!
      if (cmp_no_case(unit_len_plane_ctl%charavalue(2),'pi')) then
        ysize = pi*plane_size_ctl%realvalue(2)
      else
        ysize = plane_size_ctl%realvalue(2)
      end if
!
      if (cmp_no_case(unit_len_plane_ctl%charavalue(3),'pi')) then
        zsize = pi*plane_size_ctl%realvalue(3)
      else
        zsize = plane_size_ctl%realvalue(3)
      end if
!
      xmax = xsize / two
      ymax = ysize / two
      zmax = zsize / two
      xmin = -xmax
      ymin = -ymax
      zmin = -zmax
!
      if      (cmp_no_case(horizontal_grid_ctl%charavalue,              &
     &                     label_Chebyshev)) then
        iradi = igrid_Chebyshev
      else if (cmp_no_case(horizontal_grid_ctl%charavalue,              &
     &                     label_half_Cbyv)) then
        iradi = igrid_half_Chebyshev
      else
        iradi = igrid_equidistance
      end if
!
!
      if (nnod_plane_ctl%iflag .gt. 0) then
        nx_all = nnod_plane_ctl%intvalue(1)
        ny_all = nnod_plane_ctl%intvalue(2)
        nz_all = nnod_plane_ctl%intvalue(3)
      else
        nx_all = 2
        ny_all = 2
        nz_all = 2
      end if
!
!
      if (ndomain_plane_ctl%iflag .gt. 0) then
        ndx = ndomain_plane_ctl%intvalue(1)
        ndy = ndomain_plane_ctl%intvalue(2)
        ndz = ndomain_plane_ctl%intvalue(3)
      else
        ndx = 1
        ndy = 1
        ndz = 1
      end if
!
      if ( ( mod(nx_all,ndx) /= 0 ) .or.                                &
     &     ( mod(ny_all,ndy) /= 0 ) .or.                                &
     &     ( mod(nz_all,ndz) /= 0 )      ) then
        stop ' ***** error : illegal input '
      endif
      if (neib .ge. (nz_all/ndz) ) then
        write(*,*) 'neighbouring should be less than', (nz_all/ndz)
        stop ' ***** error : illegal input '
      end if
!
      if (num_of_sleeve_ctl%iflag .gt. 0) then
        ndepth = num_of_sleeve_ctl%intvalue
      else
        ndepth = 1
      end if
!
!
      if (num_z_filter_ctl%iflag .gt. 0) then
        iflag_filter = num_z_filter_ctl%intvalue
      else
        iflag_filter = 0
      end if
!
      iflag_z_filter = 0
!      write(*,*) ' Setting of vertical filter'
!      call skip_comment(character_4_read, l_in)
!      read (character_4_read, * )   iflag_z_filter
!      write(*,*) ' iflag_z_filter    = ',  iflag_z_filter
!
!
!
      if (omitting_value_ctl%iflag .gt. 0) then
        eps_filter = omitting_value_ctl%realvalue
      else
        eps_filter = 0.0d0
      end if
!
!
      write(*,*) 'mesh_file_header:   ', trim(mesh_file_header)
      write(*,*) 'filter_file_header: ', trim(filter_file_header)
      write(*,*) 'z_filter_header:    ', trim(z_filter_header)
!
      write(*,'(a,1p3E25.15e3)') 'size of domain: ', xsize,ysize,zsize
      write(*,*) 'grid type: ', iradi
!
      write(*,*) 'num. of grid: ',   nx_all, ny_all, nz_all
      write(*,*) 'num. of domain: ', ndx, ndy, ndz
      write(*,*) 'num. of sleeve: ', ndepth
!
      write(*,*) 'num. of filter: ', iflag_filter
      write(*,*) 'omitting parameter: ', eps_filter
!
      end subroutine s_set_ctl_data_plane_mesh
!
!  ---------------------------------------------------------------------
!
      end module set_ctl_data_plane_mesh
