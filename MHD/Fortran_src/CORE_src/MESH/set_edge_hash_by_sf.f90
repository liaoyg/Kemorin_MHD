!>@file   set_edge_hash_by_sf.f90
!!@brief  module set_edge_hash_by_sf
!!
!!@author H. Matsui
!!@date Programmed in ????
!
!>@brief Hash table using sum of local node ID
!!
!!@verbatim
!!      subroutine const_edge_hash_4_sf                                 &
!!     &         (numnod, numsurf, nnod_4_surf, nnod_4_edge, ie_surf,   &
!!     &          inum_edge_hash, istack_edge_hash,                     &
!!     &          iend_edge_hash, iedge_hash, iedge_flag)
!!
!!      subroutine const_part_edge_hash_4_sf(numnod, numsurf,           &
!!     &           numsurf_part, nnod_4_surf, nnod_4_edge, ie_surf,     &
!!     &           isf_part, inum_edge_hash, istack_edge_hash,          &
!!     &           iend_edge_hash, iedge_hash, iedge_flag)
!!@endverbatim
!
      module set_edge_hash_by_sf
!
      use m_precision
!
      use m_geometry_constants
!
      implicit none
!
      private :: set_edge_hash_4_sf, set_part_edge_hash_4_sf
      private :: mark_all_edges_by_sf
!
!------------------------------------------------------------------
!
      contains
!
!------------------------------------------------------------------
!
      subroutine const_edge_hash_4_sf                                   &
     &         (numnod, numsurf, nnod_4_surf, nnod_4_edge, ie_surf,     &
     &          inum_edge_hash, istack_edge_hash,                       &
     &          iend_edge_hash, iedge_hash, iedge_flag)
!
      integer(kind = kint), intent(in) :: numnod, numsurf
      integer(kind = kint), intent(in) :: nnod_4_surf, nnod_4_edge
      integer(kind = kint), intent(in) :: ie_surf(numsurf,nnod_4_surf)
!
      integer(kind = kint), intent(inout) :: iend_edge_hash
      integer(kind = kint), intent(inout)                               &
     &                     :: istack_edge_hash(0:nnod_4_edge*numnod)
      integer(kind = kint), intent(inout)                               &
     &                     :: inum_edge_hash(nnod_4_edge*numnod)
      integer(kind = kint), intent(inout)                               &
     &                     :: iedge_hash(nedge_4_surf*numsurf,2)
      integer(kind = kint), intent(inout)                               &
     &                     :: iedge_flag(nedge_4_surf*numsurf)
!
!
      call set_edge_hash_4_sf                                           &
     &   (numnod, numsurf, nnod_4_surf, nnod_4_edge, ie_surf,           &
     &    inum_edge_hash, istack_edge_hash, iend_edge_hash, iedge_hash)
!
      call mark_all_edges_by_sf                                         &
     &   (numnod, numsurf, nnod_4_surf, nnod_4_edge, ie_surf,           &
     &    istack_edge_hash, iend_edge_hash, iedge_hash, iedge_flag)
!
      end subroutine const_edge_hash_4_sf
!
!------------------------------------------------------------------
!
      subroutine const_part_edge_hash_4_sf(numnod, numsurf,             &
     &         numsurf_part, nnod_4_surf, nnod_4_edge, ie_surf,         &
     &         isf_part, inum_edge_hash, istack_edge_hash,              &
     &         iend_edge_hash, iedge_hash, iedge_flag)
!
      integer(kind = kint), intent(in) :: numnod, numsurf, numsurf_part
      integer(kind = kint), intent(in) :: nnod_4_surf, nnod_4_edge
      integer(kind = kint), intent(in) :: isf_part(numsurf_part)
      integer(kind = kint), intent(in) :: ie_surf(numsurf,nnod_4_surf)
!
      integer(kind = kint), intent(inout) :: iend_edge_hash
      integer(kind = kint), intent(inout)                               &
     &                     :: istack_edge_hash(0:nnod_4_edge*numnod)
      integer(kind = kint), intent(inout)                               &
     &                     :: inum_edge_hash(nnod_4_edge*numnod)
      integer(kind = kint), intent(inout)                               &
     &                     :: iedge_hash(nedge_4_surf*numsurf,2)
      integer(kind = kint), intent(inout)                               &
     &                     :: iedge_flag(nedge_4_surf*numsurf)
!
!
      call set_part_edge_hash_4_sf(numnod, numsurf, numsurf_part,       &
     &    nnod_4_surf, nnod_4_edge, ie_surf, isf_part,                  &
     &    inum_edge_hash, istack_edge_hash, iend_edge_hash, iedge_hash)
!
      call mark_all_edges_by_sf                                         &
     &   (numnod, numsurf, nnod_4_surf, nnod_4_edge, ie_surf,           &
     &    istack_edge_hash, iend_edge_hash, iedge_hash, iedge_flag)
!
      end subroutine const_part_edge_hash_4_sf
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
      subroutine set_edge_hash_4_sf                                     &
     &         (numnod, numsurf, nnod_4_surf, nnod_4_edge, ie_surf,     &
     &          inum_edge_hash, istack_edge_hash,                       &
     &          iend_edge_hash, iedge_hash)
!
      integer(kind = kint), intent(in) :: numnod, numsurf
      integer(kind = kint), intent(in) :: nnod_4_surf, nnod_4_edge
      integer(kind = kint), intent(in) :: ie_surf(numsurf,nnod_4_surf)
!
      integer(kind = kint), intent(inout)                               &
     &                     :: istack_edge_hash(0:nnod_4_edge*numnod)
      integer(kind = kint), intent(inout)                               &
     &                     :: inum_edge_hash(nnod_4_edge*numnod)
      integer(kind = kint), intent(inout)                               &
     &                     :: iedge_hash(nedge_4_surf*numsurf,2)
      integer(kind = kint), intent(inout) :: iend_edge_hash
!
      integer(kind = kint) :: isurf, is1, is2, k1
      integer(kind = kint) :: ihash, icou
!
!
! Count numbers
      do isurf = 1, numsurf
        do k1 = 1, nedge_4_surf
          is1 = node_on_edge_sf_l(1,k1)
          is2 = node_on_edge_sf_l(2,k1)
          ihash = ie_surf(isurf,is1) + ie_surf(isurf,is2)
!
          inum_edge_hash(ihash) = inum_edge_hash(ihash) + 1
        end do
      end do
!
! Set stacks
      istack_edge_hash = 0
      do ihash = 1, nnod_4_edge*numnod
        istack_edge_hash(ihash) = istack_edge_hash(ihash-1)             &
     &                               + inum_edge_hash(ihash)
        if (istack_edge_hash(ihash) .le. (nedge_4_surf*numsurf) ) then
          iend_edge_hash = ihash
        end if
      end do
!
! Set ID
      inum_edge_hash = 0
      do isurf = 1, numsurf
        do k1 = 1, nedge_4_surf
          is1 = node_on_edge_sf_l(1,k1)
          is2 = node_on_edge_sf_l(2,k1)
          ihash = ie_surf(isurf,is1) + ie_surf(isurf,is2)
!
          inum_edge_hash(ihash) = inum_edge_hash(ihash) + 1
          icou = istack_edge_hash(ihash-1) + inum_edge_hash(ihash)
          iedge_hash(icou,1) = isurf
          iedge_hash(icou,2) = k1
        end do
      end do
!
      end subroutine set_edge_hash_4_sf
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
      subroutine set_part_edge_hash_4_sf(numnod, numsurf,               &
     &         numsurf_part, nnod_4_surf, nnod_4_edge, ie_surf,         &
     &         isf_part, inum_edge_hash, istack_edge_hash,              &
     &         iend_edge_hash, iedge_hash)
!
      integer(kind = kint), intent(in) :: numnod, numsurf, numsurf_part
      integer(kind = kint), intent(in) :: nnod_4_surf, nnod_4_edge
      integer(kind = kint), intent(in) :: isf_part(numsurf_part)
      integer(kind = kint), intent(in) :: ie_surf(numsurf,nnod_4_surf)
!
      integer(kind = kint), intent(inout)                               &
     &                     :: istack_edge_hash(0:nnod_4_edge*numnod)
      integer(kind = kint), intent(inout)                               &
     &                     :: inum_edge_hash(nnod_4_edge*numnod)
      integer(kind = kint), intent(inout)                               &
     &                     :: iedge_hash(nedge_4_surf*numsurf,2)
      integer(kind = kint), intent(inout) :: iend_edge_hash
!
      integer(kind = kint) :: inum, isurf, is1, is2, k1
      integer(kind = kint) :: ihash, icou
!
! Count numbers
      do inum = 1, numsurf_part
        isurf = abs( isf_part(inum) )
        do k1 = 1, nedge_4_surf
          is1 = node_on_edge_sf_l(1,k1)
          is2 = node_on_edge_sf_l(2,k1)
          ihash = ie_surf(isurf,is1) + ie_surf(isurf,is2)
!
          inum_edge_hash(ihash) = inum_edge_hash(ihash) + 1
        end do
      end do
!
! Set stacks
      istack_edge_hash = 0
      do ihash = 1, nnod_4_edge*numnod
        istack_edge_hash(ihash) = istack_edge_hash(ihash-1)             &
     &                               + inum_edge_hash(ihash)
        if (istack_edge_hash(ihash) .le. (nedge_4_surf*numsurf_part) )  &
     &   then
          iend_edge_hash = ihash
        end if
      end do
!
! Set ID
      inum_edge_hash = 0
      do inum = 1, numsurf_part
        isurf = abs( isf_part(inum) )
        do k1 = 1, nedge_4_surf
          is1 = node_on_edge_sf_l(1,k1)
          is2 = node_on_edge_sf_l(2,k1)
          ihash = ie_surf(isurf,is1) + ie_surf(isurf,is2)
!
          inum_edge_hash(ihash) = inum_edge_hash(ihash) + 1
          icou = istack_edge_hash(ihash-1) + inum_edge_hash(ihash)
          iedge_hash(icou,1) = isurf
          iedge_hash(icou,2) = k1
        end do
      end do
!
      end subroutine set_part_edge_hash_4_sf
!
!------------------------------------------------------------------
!
      subroutine mark_all_edges_by_sf                                   &
     &       (numnod, numsurf, nnod_4_surf, nnod_4_edge, ie_surf,       &
     &        istack_edge_hash, iend_edge_hash, iedge_hash, iedge_flag)
!
      integer(kind = kint), intent(in) :: numnod, numsurf
      integer(kind = kint), intent(in) :: nnod_4_surf, nnod_4_edge
      integer(kind = kint), intent(in) :: ie_surf(numsurf,nnod_4_surf)
!
      integer(kind = kint), intent(in)                                  &
     &                     :: istack_edge_hash(0:nnod_4_edge*numnod)
      integer(kind = kint), intent(in)                                  &
     &                     :: iedge_hash(nedge_4_surf*numsurf,2)
      integer(kind = kint), intent(in) :: iend_edge_hash
!
      integer(kind = kint), intent(inout)                               &
     &                     :: iedge_flag(nedge_4_surf*numsurf)
!
      integer(kind = kint) :: isurf, iedge, inod1, inod2, is1, is2
      integer(kind = kint) :: jsurf, jedge, jnod1, jnod2, js1, js2
      integer(kind = kint) :: ihash
      integer(kind = kint) :: ist, ied, k1, k2
!      integer(kind= kint_gl) :: i1_gl, i2_gl
!
!
      iedge_flag = 0
      do ihash = 1, iend_edge_hash
        ist = istack_edge_hash(ihash-1)+1
        ied = istack_edge_hash(ihash)
        if (ied .eq. ist) then
          iedge_flag(ist) = ist
        else if (ied .gt. ist) then
          do k1 = ist, ied
            if(iedge_flag(k1) .eq. 0) then
              iedge_flag(k1) = k1
              isurf = iedge_hash(k1,1)
              iedge = iedge_hash(k1,2)
              is1 = node_on_edge_sf_l(1,iedge)
              is2 = node_on_edge_sf_l(2,iedge)
              inod1 = ie_surf(isurf,is1)
              inod2 = ie_surf(isurf,is2)
              do k2 = k1+1, ied
                jsurf = iedge_hash(k2,1)
                jedge = iedge_hash(k2,2)
                js1 = node_on_edge_sf_l(1,jedge)
                js2 = node_on_edge_sf_l(2,jedge)
                jnod1 = ie_surf(jsurf,js1)
                jnod2 = ie_surf(jsurf,js2)
                if ( (inod2-inod1) .eq. (jnod2-jnod1) ) then
                  iedge_flag(k2) = k1
                else if ( (inod2-inod1) .eq. (jnod1-jnod2) ) then
                  iedge_flag(k2) =-k1
                end if
              end do
            end if
          end do
        end if
      end do
!
      end subroutine mark_all_edges_by_sf
!
!------------------------------------------------------------------
!
      end module set_edge_hash_by_sf
