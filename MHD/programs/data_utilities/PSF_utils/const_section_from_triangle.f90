!>@file   const_section_from_triangle.f90
!!@brief  module const_section_from_triangle
!!
!!@author H. Matsui
!!@date Programmed in 2013
!
!>@brief Construct line data for sections of triangles
!!
!!@verbatim
!!      subroutine allocate_edge_section_flags(nedge_tri)
!!      subroutine deallocate_edge_section_flags
!!      subroutine count_section_fld_in_triangle(nnod_tri, nedge_tri,   &
!!     &          nfield, ncomptot, iedge_4_tri, ref_tri, reference,    &
!!     &          nfield_edge, ncomp_edge, nnod_edge)
!!      subroutine set_section_list_in_triangle                         &
!!     &         (nedge_tri, nnod_edge, iedge_4_line)
!!      subroutine set_section_fld_in_triangle(nnod_tri, nedge_tri,     &
!!     &          nfield, ncomptot, xx_tri, iedge_4_tri,                &
!!     &          ncomp_tri, tri_phys_name, d_tri, nnod_edge,           &
!!     &          iedge_4_line, coef_on_edge, inod_gl,                  &
!!     &          num_edge_comp, edge_phys_name, d_edge)
!!@endverbatim
!
      module const_section_from_triangle
!
      use m_precision
      use m_constants
!
      implicit none
!
      integer(kind = kint), allocatable :: inod_on_edge(:)
      private :: inod_on_edge
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine allocate_edge_section_flags(nedge_tri)
!
      integer(kind = kint), intent(in) :: nedge_tri
!
!
      allocate(inod_on_edge(nedge_tri))
      if(nedge_tri .gt. 0) inod_on_edge(1:nedge_tri) = 0
!
      end subroutine allocate_edge_section_flags
!
!-----------------------------------------------------------------------
!
      subroutine deallocate_edge_section_flags
!
!
      deallocate(inod_on_edge)
!
      end subroutine deallocate_edge_section_flags
!
!-----------------------------------------------------------------------
!
      subroutine count_section_fld_in_triangle(nnod_tri, nedge_tri,     &
     &          nfield, ncomptot, iedge_4_tri, ref_tri, reference,      &
     &          nfield_edge, ncomp_edge, nnod_edge)
!
      integer(kind = kint), intent(in) :: nnod_tri, nedge_tri
      integer(kind = kint), intent(in) :: nfield, ncomptot
      real(kind = kreal), intent(in) :: ref_tri(nnod_tri)
      real(kind = kreal), intent(in) :: reference
      integer(kind = kint), intent(in) :: iedge_4_tri(nedge_tri,2)
!
      integer(kind = kint_gl), intent(inout) :: nnod_edge
      integer(kind = kint), intent(inout) :: nfield_edge
      integer(kind = kint), intent(inout) :: ncomp_edge
!
      integer(kind = kint) :: i1, i2, icou, iedge
      real(kind = kreal) :: c1, c2
!
!
      nfield_edge = nfield
      ncomp_edge =  ncomptot
!
      icou = 0
      inod_on_edge(1:nedge_tri) = 0
      do iedge = 1, nedge_tri
        i1 = iedge_4_tri(iedge,1)
        i2 = iedge_4_tri(iedge,2)
        c1 = ref_tri(i1)
        c2 = ref_tri(i2)
!
        if( ((c1-reference)*(c2-reference)) .le. 0.0d0) then
          icou = icou + 1
          inod_on_edge(iedge) = icou
        end if
      end do
!
      do i1 = 1, nedge_tri
        if( inod_on_edge(i1) .gt. 0) nnod_edge = nnod_edge + 1
      end do
!
      end subroutine count_section_fld_in_triangle
!
!-----------------------------------------------------------------------
!
      subroutine set_section_list_in_triangle                           &
     &         (nedge_tri, nnod_edge, iedge_4_line)
!
      integer(kind = kint), intent(in) :: nedge_tri
      integer(kind = kint), intent(in) :: nnod_edge
!
      integer(kind = kint), intent(inout) :: iedge_4_line(nnod_edge)
!
      integer(kind = kint) :: icou, iedge
!
!
      icou = 0
      do iedge = 1, nedge_tri
        if(inod_on_edge(iedge) .gt. 0) then
          icou = icou + 1
          iedge_4_line(icou) = iedge
        end if
      end do
!
      end subroutine set_section_list_in_triangle
!
!-----------------------------------------------------------------------
!
      subroutine set_section_fld_in_triangle(nnod_tri, nedge_tri,       &
     &          nfield, ncomptot, xx_tri, iedge_4_tri,                  &
     &          ncomp_tri, tri_phys_name, d_tri, nnod_edge,             &
     &          iedge_4_line, coef_on_edge, inod_gl,                    &
     &          xx_edge, num_edge_comp, edge_phys_name, d_edge)
!
      integer(kind = kint), intent(in) :: nnod_tri, nedge_tri
      integer(kind = kint), intent(in) :: nfield, ncomptot
      real(kind = kreal), intent(in) :: xx_tri(nnod_tri,3)
      real(kind = kreal), intent(in) :: d_tri(nnod_tri,ncomptot)
      integer(kind = kint), intent(in) :: ncomp_tri(nfield)
      character(len=kchara), intent(in) :: tri_phys_name(nfield)
!
      integer(kind = kint), intent(in) :: iedge_4_tri(nedge_tri,2)
!
      integer(kind = kint_gl), intent(in) :: nnod_edge
      integer(kind = kint), intent(in) :: iedge_4_line(nnod_edge)
      real(kind = kreal), intent(in) ::   coef_on_edge(nnod_edge,2)
!
      integer(kind = kint), intent(inout) :: num_edge_comp(nfield)
      character(len=kchara), intent(inout) :: edge_phys_name(nfield)
      integer(kind = kint_gl), intent(inout) :: inod_gl(nnod_edge)
      real(kind = kreal), intent(inout) :: xx_edge(nnod_edge,3)
      real(kind = kreal), intent(inout) :: d_edge(nnod_edge,ncomptot)
!
      integer(kind = kint) :: i1, i2, iedge
      integer(kind = kint_gl) :: inod, icou
!
!
      num_edge_comp(1:nfield) =  ncomp_tri(1:nfield)
      edge_phys_name(1:nfield) = tri_phys_name(1:nfield)
!
!$omp parallel do
      do inod = 1, nnod_edge
        inod_gl(inod) = inod
      end do
!$omp end parallel do
!
!$omp parallel do private(icou,iedge,i1,i2)
      do icou = 1, nnod_edge
        iedge = iedge_4_line(icou)
        i1 = iedge_4_tri(iedge,1)
        i2 = iedge_4_tri(iedge,2)
        xx_edge(icou,1:3) =  coef_on_edge(icou,1) * xx_tri(i1,1:3)      &
     &                     + coef_on_edge(icou,2) * xx_tri(i2,1:3)
        d_edge(icou,1:ncomptot)                                         &
     &                 =  coef_on_edge(icou,1) * d_tri(i1,1:ncomptot)   &
     &                  + coef_on_edge(icou,2) * d_tri(i2,1:ncomptot)
      end do
!$omp end parallel do
!
      end subroutine set_section_fld_in_triangle
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine count_sections_in_triangle(nele_tri, ie_tri,           &
     &          iedge_4_tri, nele_edge, nnod_4_edge)
!
      integer(kind = kint), intent(in) :: nele_tri
      integer(kind = kint), intent(in) :: ie_tri(nele_tri,3)
      integer(kind = kint), intent(in) :: iedge_4_tri(nele_tri,3)
!
      integer(kind = kint_gl), intent(inout) :: nele_edge
      integer(kind = kint), intent(inout) :: nnod_4_edge
!
      integer(kind = kint_gl) :: icou
      integer(kind = kint) :: isurf
      integer(kind = kint) :: inod1, inod2, inod3
      integer(kind = kint) :: iege1, iege2, iege3
!
!
      icou = 0
      do isurf = 1, nele_tri
        inod1 = ie_tri(isurf,1)
        inod2 = ie_tri(isurf,2)
        inod3 = ie_tri(isurf,3)
        iege1 = iedge_4_tri(isurf,1)
        iege2 = iedge_4_tri(isurf,2)
        iege3 = iedge_4_tri(isurf,3)
!
        if     (inod_on_edge(iege1).gt.0                                &
     &    .and. inod_on_edge(iege2).gt.0                                &
     &    .and. inod_on_edge(iege3).gt.0) then
          icou = icou + 1
        else if(inod_on_edge(iege1).gt.0                                &
     &    .and. inod_on_edge(iege2).gt.0) then
          icou = icou + 1
        else if(inod_on_edge(iege2).gt.0                                &
     &    .and. inod_on_edge(iege3).gt.0) then
          icou = icou + 1
        else if(inod_on_edge(iege3).gt.0                                &
     &    .and. inod_on_edge(iege1).gt.0) then
          icou = icou + 1
        end if
      end do
      nele_edge = icou
      nnod_4_edge = 2
!
      end subroutine count_sections_in_triangle
!
!-----------------------------------------------------------------------
!
      subroutine set_sections_in_triangle(nele_tri, ie_tri,             &
     &          iedge_4_tri, nele_edge, iedge_gl, ie_edge)
!
      integer(kind = kint), intent(in) :: nele_tri
      integer(kind = kint), intent(in) :: ie_tri(nele_tri,3)
      integer(kind = kint), intent(in) :: iedge_4_tri(nele_tri,3)
!
      integer(kind = kint_gl), intent(in) :: nele_edge
      integer(kind = kint_gl), intent(inout) :: iedge_gl(nele_edge)
      integer(kind = kint_gl), intent(inout) :: ie_edge(nele_edge,2)
!
      integer(kind = kint_gl) :: iedge
      integer(kind = kint) :: icou, isurf
      integer(kind = kint) :: inod1, inod2, inod3
      integer(kind = kint) :: iege1, iege2, iege3
!
!
      do iedge = 1, nele_edge
        iedge_gl(iedge) = iedge
      end do
!
      icou = 0
      do isurf = 1, nele_tri
        inod1 = ie_tri(isurf,1)
        inod2 = ie_tri(isurf,2)
        inod3 = ie_tri(isurf,3)
        iege1 = iedge_4_tri(isurf,1)
        iege2 = iedge_4_tri(isurf,2)
        iege3 = iedge_4_tri(isurf,3)
!
        if     (inod_on_edge(iege1).gt.0                                &
     &    .and. inod_on_edge(iege2).gt.0                                &
     &    .and. inod_on_edge(iege3).gt.0) then
          icou = icou + 1
          if     (inod_on_edge(iege1) .eq. inod_on_edge(iege2)) then
            ie_edge(icou,1) = inod_on_edge(iege3)
            ie_edge(icou,2) = inod_on_edge(iege1)
          else if(inod_on_edge(iege2) .eq. inod_on_edge(iege3)) then
            ie_edge(icou,1) = inod_on_edge(iege1)
            ie_edge(icou,2) = inod_on_edge(iege2)
          else 
            ie_edge(icou,1) = inod_on_edge(iege2)
            ie_edge(icou,2) = inod_on_edge(iege3)
          end if
        else if(inod_on_edge(iege1).gt.0                                &
     &    .and. inod_on_edge(iege2).gt.0) then
          icou = icou + 1
          ie_edge(icou,1) = inod_on_edge(iege1)
          ie_edge(icou,2) = inod_on_edge(iege2)
        else if(inod_on_edge(iege2).gt.0                                &
     &    .and. inod_on_edge(iege3).gt.0) then
          icou = icou + 1
          ie_edge(icou,1) = inod_on_edge(iege2)
          ie_edge(icou,2) = inod_on_edge(iege3)
        else if(inod_on_edge(iege3).gt.0                                &
     &    .and. inod_on_edge(iege1).gt.0) then
          icou = icou + 1
          ie_edge(icou,1) = inod_on_edge(iege3)
          ie_edge(icou,2) = inod_on_edge(iege1)
        end if
      end do
!
      end subroutine set_sections_in_triangle
!
!-----------------------------------------------------------------------
!
      end module const_section_from_triangle
