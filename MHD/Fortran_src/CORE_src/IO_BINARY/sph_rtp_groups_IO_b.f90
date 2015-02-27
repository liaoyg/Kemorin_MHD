!sph_rtp_groups_IO_b.f90
!      module sph_rtp_groups_IO_b
!
!     Written by H. Matsui on July, 2007
!
!      subroutine read_geom_rtp_groups_b(mesh_file_id)
!      subroutine write_geom_rtp_groups_b(mesh_file_id)
!
      module sph_rtp_groups_IO_b
!
      use m_precision
!
      use m_group_data_sph_specr_IO
      use group_data_IO_b
!
      implicit none
!
      private :: read_rtp_node_grp_data_b, write_rtp_node_grp_data_b
      private :: read_rtp_radial_grp_data_b
      private :: write_rtp_radial_grp_data_b
      private :: read_rtp_theta_grp_data_b, write_rtp_theta_grp_data_b
      private :: read_rtp_zonal_grp_data_b, write_rtp_zonal_grp_data_b
!
!------------------------------------------------------------------
!
      contains
!
!------------------------------------------------------------------
!
      subroutine read_geom_rtp_groups_b(mesh_file_id)
!
      integer(kind = kint), intent(in) :: mesh_file_id
!
      call read_rtp_node_grp_data_b(mesh_file_id)
      call read_rtp_radial_grp_data_b(mesh_file_id)
      call read_rtp_theta_grp_data_b(mesh_file_id)
      call read_rtp_zonal_grp_data_b(mesh_file_id)
!
      end subroutine read_geom_rtp_groups_b
!
!------------------------------------------------------------------
!
      subroutine write_geom_rtp_groups_b(mesh_file_id)
!
      integer(kind = kint), intent(in) :: mesh_file_id
!
      call write_rtp_node_grp_data_b(mesh_file_id)
      call write_rtp_radial_grp_data_b(mesh_file_id)
      call write_rtp_theta_grp_data_b(mesh_file_id)
      call write_rtp_zonal_grp_data_b(mesh_file_id)
!
      end subroutine write_geom_rtp_groups_b
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
      subroutine read_rtp_node_grp_data_b(mesh_file_id)
!
      integer(kind = kint), intent(in) :: mesh_file_id
!
!
      read(mesh_file_id) num_bc_grp_rtp_IO
!
      call allocate_rtp_nod_grp_IO_stack
!
      if (num_bc_grp_rtp_IO .gt. 0) then
        call read_group_stack_b(mesh_file_id, num_bc_grp_rtp_IO,        &
     &      ntot_bc_grp_rtp_IO, istack_bc_grp_rtp_IO)
!
        call allocate_rtp_nod_grp_IO_item
        call read_group_item_b(mesh_file_id, num_bc_grp_rtp_IO,         &
     &      ntot_bc_grp_rtp_IO, istack_bc_grp_rtp_IO,                   &
     &      name_bc_grp_rtp_IO,item_bc_grp_rtp_IO)
!
      else
        ntot_bc_grp_rtp_IO = 0
        call allocate_rtp_nod_grp_IO_item
      end if
!
      end subroutine read_rtp_node_grp_data_b
!
!------------------------------------------------------------------
!
      subroutine read_rtp_radial_grp_data_b(mesh_file_id)
!
      integer(kind = kint), intent(in) :: mesh_file_id
!
!
      read(mesh_file_id) num_radial_grp_rtp_IO
!
      call allocate_rtp_r_grp_IO_stack
!
      if (num_radial_grp_rtp_IO .gt. 0) then
        call read_group_stack_b(mesh_file_id, num_radial_grp_rtp_IO,    &
     &      ntot_radial_grp_rtp_IO, istack_radial_grp_rtp_IO)
!
        call allocate_rtp_r_grp_IO_item
        call read_group_item_b(mesh_file_id, num_radial_grp_rtp_IO,     &
     &      ntot_radial_grp_rtp_IO, istack_radial_grp_rtp_IO,           &
     &      name_radial_grp_rtp_IO, item_radial_grp_rtp_IO)
!
      else
        ntot_radial_grp_rtp_IO = 0
        call allocate_rtp_r_grp_IO_item
      end if
!
      end subroutine read_rtp_radial_grp_data_b
!
!------------------------------------------------------------------
!
      subroutine read_rtp_theta_grp_data_b(mesh_file_id)
!
      integer(kind = kint), intent(in) :: mesh_file_id
!
!
      read(mesh_file_id) num_theta_grp_rtp_IO
!
      call allocate_rtp_t_grp_IO_stack
!
      if (num_theta_grp_rtp_IO .gt. 0) then
        call read_group_stack_b(mesh_file_id, num_theta_grp_rtp_IO,     &
     &      ntot_theta_grp_rtp_IO, istack_theta_grp_rtp_IO)
!
        call allocate_rtp_t_grp_IO_item
        call read_group_item_b(mesh_file_id, num_theta_grp_rtp_IO,      &
     &      ntot_theta_grp_rtp_IO, istack_theta_grp_rtp_IO,             &
     &      name_theta_grp_rtp_IO, item_theta_grp_rtp_IO)
!
      else
        ntot_theta_grp_rtp_IO = 0
        call allocate_rtp_t_grp_IO_item
      end if
!
      end subroutine read_rtp_theta_grp_data_b
!
!------------------------------------------------------------------
!
      subroutine read_rtp_zonal_grp_data_b(mesh_file_id)
!
      integer(kind = kint), intent(in) :: mesh_file_id
!
!
      read(mesh_file_id) num_zonal_grp_rtp_IO
!
      call allocate_rtp_p_grp_IO_stack
!
      if (num_zonal_grp_rtp_IO .gt. 0) then
        call read_group_stack_b(mesh_file_id, num_zonal_grp_rtp_IO,     &
     &      ntot_zonal_grp_rtp_IO, istack_zonal_grp_rtp_IO)
!
        call allocate_rtp_p_grp_IO_item
        call read_group_item_b(mesh_file_id, num_zonal_grp_rtp_IO,      &
     &      ntot_zonal_grp_rtp_IO, istack_zonal_grp_rtp_IO,             &
     &      name_zonal_grp_rtp_IO, item_zonal_grp_rtp_IO)
!
      else
        ntot_zonal_grp_rtp_IO = 0
        call allocate_rtp_p_grp_IO_item
      end if
!
      end subroutine read_rtp_zonal_grp_data_b
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
      subroutine write_rtp_node_grp_data_b(mesh_file_id)
!
      integer(kind = kint), intent(in) :: mesh_file_id
!
!
      call write_group_data_b(mesh_file_id, num_bc_grp_rtp_IO,          &
     &    ntot_bc_grp_rtp_IO, istack_bc_grp_rtp_IO, name_bc_grp_rtp_IO, &
     &    item_bc_grp_rtp_IO)
!
      call deallocate_rtp_nod_grp_IO_item
!
      end subroutine write_rtp_node_grp_data_b
!
!------------------------------------------------------------------
!
      subroutine write_rtp_radial_grp_data_b(mesh_file_id)
!
      integer(kind = kint), intent(in) :: mesh_file_id
!
!
      call write_group_data_b(mesh_file_id, num_radial_grp_rtp_IO,      &
     &    ntot_radial_grp_rtp_IO, istack_radial_grp_rtp_IO,             &
     &    name_radial_grp_rtp_IO, item_radial_grp_rtp_IO)
!
      call deallocate_rtp_r_grp_IO_item
!
      end subroutine write_rtp_radial_grp_data_b
!
!------------------------------------------------------------------
!
      subroutine write_rtp_theta_grp_data_b(mesh_file_id)
!
      integer(kind = kint), intent(in) :: mesh_file_id
!
!
      call write_group_data_b(mesh_file_id, num_theta_grp_rtp_IO,       &
     &    ntot_theta_grp_rtp_IO, istack_theta_grp_rtp_IO,               &
     &    name_theta_grp_rtp_IO, item_theta_grp_rtp_IO)
!
      call deallocate_rtp_t_grp_IO_item
!
      end subroutine write_rtp_theta_grp_data_b
!
!------------------------------------------------------------------
!
      subroutine write_rtp_zonal_grp_data_b(mesh_file_id)
!
      integer(kind = kint), intent(in) :: mesh_file_id
!
!
      call write_group_data_b(mesh_file_id, num_zonal_grp_rtp_IO,       &
     &    ntot_zonal_grp_rtp_IO, istack_zonal_grp_rtp_IO,               &
     &    name_zonal_grp_rtp_IO, item_zonal_grp_rtp_IO)
!
      call deallocate_rtp_p_grp_IO_item
!
      end subroutine write_rtp_zonal_grp_data_b
!
!------------------------------------------------------------------
!
      end module sph_rtp_groups_IO_b
