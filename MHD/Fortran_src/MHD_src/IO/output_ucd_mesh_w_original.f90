!output_ucd_mesh_w_original.f90
!      module output_ucd_mesh_w_original
!
!      Programmed by H.Matsui on July 2000 (ver 1.1)
!      Modified by H. Matsui on Aug, 2007
!
!     subroutine output_grd_file_w_org_connect(my_rank)
!
      module output_ucd_mesh_w_original
!
      use m_precision
      use m_constants
!
      implicit none
!
      private :: link_local_org_mesh_4_ucd
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine output_grd_file_w_org_connect
!
      use m_ucd_data
      use m_geometry_data
      use m_field_file_format
      use m_t_step_parameter
      use set_ucd_data
!
      use merged_udt_vtk_file_IO
      use parallel_ucd_IO_select
!
!
      if(i_step_output_ucd .eq. 0) return
!
      call link_fem_num_field_2_ucd_out
      call link_local_org_mesh_4_ucd
      call link_fem_field_data_2_ucd_out
!
      if (fem_ucd%ifmt_file/icent .eq. iflag_single/icent) then
        call init_merged_ucd(fem_ucd, merged_ucd)
      end if
!
      call sel_write_parallel_ucd_mesh(fem_ucd, merged_ucd)
!
      if(   mod(fem_ucd%ifmt_file,icent)/iten .eq. iflag_udt/iten       &
     & .or. mod(fem_ucd%ifmt_file,icent)/iten .eq. iflag_vtd/iten) then
        call deallocate_ucd_ele(fem_ucd)
      end if
!
      if(mod(fem_ucd%ifmt_file,icent)/iten .eq. iflag_vtd/iten) then
        call deallocate_ucd_node(fem_ucd)
      end if
!
      end subroutine output_grd_file_w_org_connect
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine link_global_org_mesh_4_ucd
!
      use m_geometry_data
      use m_geometry_data_MHD
      use m_ucd_data
      use set_and_cal_udt_data
!
!
      call link_fem_node_data_2_ucd_out
      call const_udt_global_connect                                     &
     &   (node1%internal_node, ele1%numele, ele1%nnod_4_ele,            &
     &    iele_global_org, ie_org, fem_ucd)
!
      end subroutine link_global_org_mesh_4_ucd
!
!-----------------------------------------------------------------------
!
      subroutine link_local_org_mesh_4_ucd
!
      use m_geometry_data
      use m_geometry_data_MHD
      use m_ucd_data
      use set_and_cal_udt_data
!
!
      call const_udt_local_nodes(node1%numnod, node1%xx, fem_ucd)
      call const_udt_local_connect                                      &
     &   (node1%internal_node, ele1%numele, ele1%nnod_4_ele,            &
     &    ie_org, fem_ucd)
!
      end subroutine link_local_org_mesh_4_ucd
!
!-----------------------------------------------------------------------
!
      end module output_ucd_mesh_w_original