!
!     program  pick_surface_single
!
!      program for pick up surface connectivities form subdomain mesh
!         programmed  by  H.Matsui (U. Chicago)  on Oct. 2003 (ver 1.0)
!         Modified  by  H.Matsui (U. Chicago)  on Jan. 2007 (ver 2.0)
!
      program   pick_surface_single
!
      use m_precision
!
      use t_file_IO_parameter
      use t_mesh_data_4_merge
      use single_const_surface_mesh
      use getarg_kemo
!
      implicit    none
!
      character(len=kchara) :: file_head
      integer(kind = kint) :: icount
!
      type(field_IO_params), save ::  pick_mesh_file
!
      icount = iargc_kemo()
      if(icount .eq. 0) then
        write(*,*) ' Please input header of mesh  !!'
        read (*,*) file_head
      else
        call getarg_k(1, file_head)
      end if
!
      pick_mesh_file%file_prefix = file_head
      call choose_surface_mesh_sgl(pick_mesh_file)
!
      stop ' //// program normally finished //// '
!
      end program pick_surface_single 