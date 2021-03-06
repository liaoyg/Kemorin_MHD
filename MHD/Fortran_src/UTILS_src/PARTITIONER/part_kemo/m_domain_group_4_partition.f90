!m_domain_group_4_partition.f90
!     module m_domain_group_4_partition
!
!      Written by H. Matsui on Aug., 2007
!
!      subroutine allocate_domain_nod_group
!      subroutine allocate_domain_nese_group
!      subroutine allocate_local_ne_id_tbl
!      subroutine allocate_local_nese_id_tbl
!      subroutine allocate_org_gl_nod_id
!      subroutine allocate_org_gl_nese_id
!      subroutine allocate_finer_domain_group
!
!      subroutine deallocate_domain_nod_group
!      subroutine deallocate_domain_nese_group
!      subroutine deallocate_local_nese_id_tbl
!      subroutine deallocate_org_gl_ne_id
!      subroutine deallocate_org_gl_nese_id
!      subroutine deallocate_finer_domain_group
!
!      subroutine allocate_work_4_rcb(nnod)
!      subroutine deallocate_work_4_rcb
!
      module m_domain_group_4_partition
!
      use m_precision
!
      implicit none
!
      integer(kind = kint) :: nnod_s_domin
      integer(kind = kint) :: nele_s_domin
      integer(kind = kint) :: nsurf_s_domin
      integer(kind = kint) :: nedge_s_domin
      integer(kind = kint) :: intnod_s_domin
!
      integer(kind = kint), allocatable :: IGROUP_nod(:)
      integer(kind = kint), allocatable :: IGROUP_ele(:)
      integer(kind = kint), allocatable :: IGROUP_surf(:)
      integer(kind = kint), allocatable :: IGROUP_edge(:)
!
      integer(kind = kint), allocatable :: inod_local_part(:)
      integer(kind = kint), allocatable :: iele_local_part(:)
      integer(kind = kint), allocatable :: isurf_local_part(:)
      integer(kind = kint), allocatable :: iedge_local_part(:)
!
      integer(kind = kint_gl), allocatable :: id_glnode_org(:)
      integer(kind = kint_gl), allocatable :: id_glelem_org(:)
      integer(kind = kint_gl), allocatable :: id_glsurf_org(:)
      integer(kind = kint_gl), allocatable :: id_gledge_org(:)
!
!
      integer(kind = kint) :: nproc_finer
      integer(kind = kint) :: nnod_group_finer, internod_group_finer
      integer(kind = kint),  allocatable :: IGROUP_FINER(:)
!
      real(kind=kreal), allocatable :: VAL(:)
      integer(kind=kint), allocatable :: IS1(:)
!
!   --------------------------------------------------------------------
!
      contains
!
!   --------------------------------------------------------------------
!
      subroutine allocate_domain_nod_group
!
      allocate(IGROUP_nod(nnod_s_domin))
      allocate(IGROUP_ele(nele_s_domin))
      IGROUP_nod = 0
      IGROUP_ele = 0
!
      end subroutine allocate_domain_nod_group
!
!   --------------------------------------------------------------------
!
      subroutine allocate_domain_nese_group
!
      allocate(IGROUP_nod(nnod_s_domin))
      allocate(IGROUP_ele(nele_s_domin))
      allocate(IGROUP_edge(nedge_s_domin))
      allocate(IGROUP_surf(nsurf_s_domin))
      IGROUP_nod = 0
      IGROUP_ele = 0
      IGROUP_edge = 0
      IGROUP_surf = 0
!
      end subroutine allocate_domain_nese_group
!
!   --------------------------------------------------------------------
!
      subroutine allocate_local_ne_id_tbl
!
      allocate(inod_local_part(nnod_s_domin))
      allocate(iele_local_part(nele_s_domin))
      inod_local_part =  0
      iele_local_part =  0
!
      end subroutine allocate_local_ne_id_tbl
!
!   --------------------------------------------------------------------
!
      subroutine allocate_local_nese_id_tbl
!
      allocate(inod_local_part(nnod_s_domin))
      allocate(iele_local_part(nele_s_domin))
      allocate(isurf_local_part(nsurf_s_domin))
      allocate(iedge_local_part(nedge_s_domin))
      inod_local_part =  0
      iele_local_part =  0
      iedge_local_part = 0
      isurf_local_part = 0
!
      end subroutine allocate_local_nese_id_tbl
!
!   --------------------------------------------------------------------
!
      subroutine allocate_org_gl_nod_id
!
      allocate(id_glnode_org(nnod_s_domin))
      id_glnode_org = 0
!
      end subroutine allocate_org_gl_nod_id
!
!   --------------------------------------------------------------------
!
      subroutine allocate_org_gl_nese_id
!
      allocate(id_glnode_org(nnod_s_domin))
      allocate(id_glelem_org(nele_s_domin))
      allocate(id_glsurf_org(nsurf_s_domin))
      allocate(id_gledge_org(nedge_s_domin))
      id_glnode_org = 0
      id_glelem_org = 0
      id_glsurf_org = 0
      id_gledge_org = 0
!
      end subroutine allocate_org_gl_nese_id
!
!   --------------------------------------------------------------------
!
      subroutine allocate_finer_domain_group
!
      allocate(IGROUP_FINER(nnod_group_finer))
      IGROUP_FINER = 0
!
      end subroutine allocate_finer_domain_group
!
!   --------------------------------------------------------------------
!   --------------------------------------------------------------------
!
      subroutine deallocate_domain_nod_group
!
      deallocate(IGROUP_nod, IGROUP_ele)
!
      end subroutine deallocate_domain_nod_group
!
!   --------------------------------------------------------------------
!
      subroutine deallocate_domain_nese_group
!
      deallocate(IGROUP_nod, IGROUP_ele)
      deallocate(IGROUP_surf, IGROUP_edge)
!
      end subroutine deallocate_domain_nese_group
!
!   --------------------------------------------------------------------
!
      subroutine deallocate_local_ne_id_tbl
!
!
      deallocate(inod_local_part, iele_local_part)
!
      end subroutine deallocate_local_ne_id_tbl
!
!   --------------------------------------------------------------------
!
      subroutine deallocate_local_nese_id_tbl
!
      deallocate(inod_local_part, iele_local_part)
      deallocate(isurf_local_part, iedge_local_part)
!
      end subroutine deallocate_local_nese_id_tbl
!
!   --------------------------------------------------------------------
!
      subroutine deallocate_org_gl_ne_id
!
      deallocate(id_glnode_org, id_glelem_org)
!
      end subroutine deallocate_org_gl_ne_id
!
!   --------------------------------------------------------------------
!
      subroutine deallocate_org_gl_nese_id
!
      deallocate(id_glnode_org, id_glelem_org)
      deallocate(id_glsurf_org, id_gledge_org)
!
      end subroutine deallocate_org_gl_nese_id
!
!   --------------------------------------------------------------------
!
      subroutine deallocate_finer_domain_group
!
      deallocate(IGROUP_FINER)
!
      end subroutine deallocate_finer_domain_group
!
!   --------------------------------------------------------------------
!   --------------------------------------------------------------------
!
      subroutine allocate_work_4_rcb(nnod)
!
      integer(kind = kint), intent(in) :: nnod
!
      allocate (VAL(nnod))
      allocate (IS1(nnod))
!
      VAL = 0.0d0
      IS1 = 0
!
      end subroutine allocate_work_4_rcb
!
!   --------------------------------------------------------------------
!
      subroutine deallocate_work_4_rcb
!
      deallocate (VAL)
      deallocate (IS1)
!
      end subroutine deallocate_work_4_rcb
!
!   --------------------------------------------------------------------
!
      end module m_domain_group_4_partition
