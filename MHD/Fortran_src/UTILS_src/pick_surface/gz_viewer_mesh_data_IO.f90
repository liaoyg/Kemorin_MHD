!gz_viewer_mesh_data_IO.f90
!      module gz_viewer_mesh_data_IO
!
!      Written by Kemorin on Jan., 2007
!
!!      subroutine write_domain_data_viewer_gz
!!        type(merged_viewer_mesh), intent(in) :: mgd_view_mesh
!!      subroutine read_domain_data_viewer_gz
!!        type(merged_viewer_mesh), intent(inout) :: mgd_view_mesh
!!
!!      subroutine write_node_data_viewer_gz(view_mesh)
!!        type(viewer_mesh_data), intent(in) :: view_mesh
!!      subroutine read_node_data_viewer_gz(view_mesh)
!!        type(viewer_mesh_data), intent(inout) :: view_mesh
!!
!!      subroutine write_surf_connect_viewer_gz(nnod_4_surf, view_mesh)
!!        type(viewer_mesh_data), intent(in) :: view_mesh
!!      subroutine read_surf_connect_viewer_gz                          &
!!     &         (nnod_4_ele, nnod_4_surf, nnod_4_edge, view_mesh)
!!        type(viewer_mesh_data), intent(inout) :: view_mesh
!!
!!      subroutine write_edge_connect_viewer_gz(nnod_4_edge)
!!        type(viewer_mesh_data), intent(in) :: view_mesh
!!      subroutine read_edge_connect_viewer_gz(nnod_4_edge, view_mesh)
!!        type(viewer_mesh_data), intent(inout) :: view_mesh
!
module gz_viewer_mesh_data_IO
!
use m_precision
!
use t_viewer_mesh
use t_merged_viewer_mesh
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
subroutine write_domain_data_viewer_gz(mgd_view_mesh)
!
type(merged_viewer_mesh), intent(in) :: mgd_view_mesh
!
!
write(textbuf,'(a,a1)') '!', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)') '! number of domain ', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)')                                           &
&                 '!   stack of node for domain ', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)')                                           &
&                 '!   stack of surface for domain ', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)')                                           &
&                 '!   stack of edge for domain ', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)') '!', char(0)
call gz_write_textbuf_w_lf
!
write(textbuf,'(i15,a1)') mgd_view_mesh%num_pe_sf, char(0)
call gz_write_textbuf_w_lf
!
call write_gz_multi_int_8i10(mgd_view_mesh%num_pe_sf,             &
&    mgd_view_mesh%inod_sf_stack(1))
call write_gz_multi_int_8i10(mgd_view_mesh%num_pe_sf,             &
&    mgd_view_mesh%isurf_sf_stack(1))
call write_gz_multi_int_8i10(mgd_view_mesh%num_pe_sf,             &
&    mgd_view_mesh%iedge_sf_stack(1))
!
end subroutine write_domain_data_viewer_gz
!
!------------------------------------------------------------------
!
subroutine read_domain_data_viewer_gz(mgd_view_mesh)
!
type(merged_viewer_mesh), intent(inout) :: mgd_view_mesh
!
integer(kind = kint) :: num_pe
!
!
call skip_gz_comment_int(num_pe)
!
call alloc_num_mesh_sf(num_pe, mgd_view_mesh)
!
call read_gz_multi_int                                            &
&   (mgd_view_mesh%num_pe_sf, mgd_view_mesh%inod_sf_stack)
call read_gz_multi_int                                            &
&   (mgd_view_mesh%num_pe_sf, mgd_view_mesh%isurf_sf_stack)
call read_gz_multi_int                                            &
&   (mgd_view_mesh%num_pe_sf, mgd_view_mesh%iedge_sf_stack)
!
end subroutine read_domain_data_viewer_gz
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
subroutine write_node_data_viewer_gz(view_mesh)
!
type(viewer_mesh_data), intent(in) :: view_mesh
!
integer(kind = kint) :: i
!
!
write(textbuf,'(a,a1)') '!', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)') '! 1. node information', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)')                                           &
&                 '! number_of node, intenal_node', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)') '! Global ID, x, y, z', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)') '!', char(0)
call gz_write_textbuf_w_lf
!
write(textbuf,'(i15,a1)') view_mesh%nodpetot_viewer, char(0)
call gz_write_textbuf_w_lf
!
do i = 1, view_mesh%nodpetot_viewer
write(textbuf,1002) i, view_mesh%xx_view(i,1:3), char(0)
call gz_write_textbuf_w_lf
end do
1002 format(i15, 1p3e23.12,a1)
!
end subroutine write_node_data_viewer_gz
!
!------------------------------------------------------------------
!
subroutine read_node_data_viewer_gz(view_mesh)
!
type(viewer_mesh_data), intent(inout) :: view_mesh
!
integer(kind = kint) :: i, itmp
!
!
call skip_gz_comment_int(itmp)
!
call alloc_nod_position_viewer(view_mesh)
!
do i = 1, view_mesh%nodpetot_viewer
call get_one_line_from_gz_f
read(textbuf,*) itmp, view_mesh%xx_view(i,1:3)
end do
!
end subroutine read_node_data_viewer_gz
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
subroutine write_surf_connect_viewer_gz(nnod_4_surf, view_mesh)
!
integer(kind = kint), intent(in) :: nnod_4_surf
type(viewer_mesh_data), intent(in) :: view_mesh
!
integer(kind = kint) :: i
character(len=kchara) :: fmt_txt
!
!
write(textbuf,'(a,a1)') '!', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)') '! 2. element information', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)') '! element type', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)') '! Global ID, connectivity', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)') '!', char(0)
call gz_write_textbuf_w_lf
!
write(textbuf,'(i15,a1)') view_mesh%nsurf_viewer, char(0)
call gz_write_textbuf_w_lf
call write_gz_multi_int_10i8                                      &
&   (view_mesh%nsurf_viewer, view_mesh%surftyp_viewer)
!
write(fmt_txt,'(a5,i2,a9)')                                       &
&                '(i15,', nnod_4_surf, '(i15),a1)'
do i = 1, view_mesh%nsurf_viewer
write(textbuf,fmt_txt)                                          &
&      i, view_mesh%ie_sf_viewer(i,1:nnod_4_surf), char(0)
call gz_write_textbuf_w_lf
end do
!
end subroutine write_surf_connect_viewer_gz
!
!------------------------------------------------------------------
!
subroutine read_surf_connect_viewer_gz                            &
&         (nnod_4_ele, nnod_4_surf, nnod_4_edge, view_mesh)
!
use m_geometry_constants
use m_node_quad_2_linear_sf
use set_nnod_4_ele_by_type
!
integer(kind = kint), intent(in) :: nnod_4_ele
integer(kind = kint), intent(inout) :: nnod_4_surf
integer(kind = kint), intent(inout) :: nnod_4_edge
type(viewer_mesh_data), intent(inout) :: view_mesh
!
integer(kind = kint) :: i, itmp
!
!
call skip_gz_comment_int(itmp)
!
call alloc_surf_type_viewer(view_mesh)
!
call read_gz_multi_int                                            &
&   (view_mesh%nsurf_viewer, view_mesh%surftyp_viewer)
!
call set_3D_nnod_4_sfed_by_ele                                    &
&   (nnod_4_ele, nnod_4_surf, nnod_4_edge)
call allocate_quad4_2_linear(nnod_4_ele)
!
call alloc_surf_connect_viewer(nnod_4_surf, view_mesh)
!
do i = 1, view_mesh%nsurf_viewer
call get_one_line_from_gz_f
read(textbuf,*) itmp, view_mesh%ie_sf_viewer(i,1:nnod_4_surf)
end do
!
end subroutine read_surf_connect_viewer_gz
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
subroutine write_edge_connect_viewer_gz(nnod_4_edge, view_mesh)
!
use m_geometry_constants
!
integer(kind = kint), intent(in) :: nnod_4_edge
type(viewer_mesh_data), intent(in) :: view_mesh
!
integer(kind = kint) :: i
character(len=kchara) :: fmt_txt
!
!
write(textbuf,'(a,a1)') '!', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)') '!  edge information', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)') '!  edge type', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)') '!  global ID, connectivity', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)') '!', char(0)
call gz_write_textbuf_w_lf
!
write(textbuf,'(i15,a1)') view_mesh%edgepetot_viewer, char(0)
call gz_write_textbuf_w_lf
!
!
write(fmt_txt,'(a5,i2,a9)') '(i15,', nnod_4_edge, '(i15),a1)'
do i = 1, view_mesh%edgepetot_viewer
write(textbuf,fmt_txt)                                          &
&         i, view_mesh%ie_edge_viewer(i,1:nnod_4_edge), char(0)
call gz_write_textbuf_w_lf
end do
!
write(textbuf,'(a,a1)') '!', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)') '!  edge ID for surfaces', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(a,a1)') '!', char(0)
call gz_write_textbuf_w_lf
write(textbuf,'(i15,a1)') view_mesh%nsurf_viewer, char(0)
call gz_write_textbuf_w_lf
!
write(fmt_txt,'(a5,i2,a9)')                                       &
&                '(i15,', nedge_4_surf, '(i15),a1)'
do i = 1, view_mesh%nsurf_viewer
write(textbuf,fmt_txt)                                          &
& i, view_mesh%iedge_sf_viewer(i,1:nedge_4_surf), char(0)
call gz_write_textbuf_w_lf
end do
!
end subroutine write_edge_connect_viewer_gz
!
!------------------------------------------------------------------
!
subroutine read_edge_connect_viewer_gz(nnod_4_edge, view_mesh)
!
use m_geometry_constants
use m_node_quad_2_linear_sf
!
integer(kind = kint), intent(in) :: nnod_4_edge
type(viewer_mesh_data), intent(inout) :: view_mesh
!
integer(kind = kint) :: i, itmp
!
!
call skip_gz_comment_int(itmp)
!
call alloc_edge_data_4_sf(nnod_4_edge, view_mesh)
!
do i = 1, view_mesh%edgepetot_viewer
call get_one_line_from_gz_f
read(textbuf,*) itmp, view_mesh%ie_edge_viewer(i,1:nnod_4_edge)
end do
!
call skip_gz_comment_int(itmp)
!
do i = 1, view_mesh%nsurf_viewer
call get_one_line_from_gz_f
read(textbuf,*)                                                  &
&        itmp, view_mesh%iedge_sf_viewer(i,1:nedge_4_surf)
end do
!
end subroutine read_edge_connect_viewer_gz
!
!------------------------------------------------------------------
!
end module gz_viewer_mesh_data_IO
