!
!      module set_itp_destIO_2_org
!
!        programmed by H.Matsui on Sep. 2006
!
!!      subroutine count_num_interpolation_4_orgin                      &
!!     &         (n_org_rank, n_dest_rank, itp_org)
!!      subroutine set_interpolation_4_orgin(n_org_rank, itp_org)
!!        type(interpolate_table_org), intent(inout) :: itp_org
!
      module set_itp_destIO_2_org
!
      use m_precision
!
      use m_interpolate_table_dest_IO
      use t_interpolate_tbl_org
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine count_num_interpolation_4_orgin                        &
     &         (n_org_rank, n_dest_rank, itp_org)
!
      integer(kind = kint), intent(in) :: n_org_rank, n_dest_rank
      type(interpolate_table_org), intent(inout) :: itp_org
!
      integer(kind = kint) :: i
!
      do i = 1, IO_itp_dest%num_org_domain
!
        if (IO_itp_dest%id_org_domain(i) .eq. n_org_rank) then
          itp_org%num_dest_domain = itp_org%num_dest_domain + 1
          itp_org%id_dest_domain(itp_org%num_dest_domain)               &
     &       = n_dest_rank
          itp_org%istack_nod_tbl_org(itp_org%num_dest_domain)           &
     &       = itp_org%istack_nod_tbl_org(itp_org%num_dest_domain-1)    &
     &                        + IO_itp_dest%istack_nod_tbl_dest(i)      &
     &                        - IO_itp_dest%istack_nod_tbl_dest(i-1)
        end if
      end do
!
      call deallocate_itp_nod_dst_IO
      call deallocate_itp_num_dst_IO
!
      end subroutine count_num_interpolation_4_orgin
!
!-----------------------------------------------------------------------
!
      subroutine set_interpolation_4_orgin(n_org_rank, itp_org)
!
      use m_work_const_itp_table
!
      integer(kind = kint), intent(in) :: n_org_rank
      type(interpolate_table_org), intent(inout) :: itp_org
!
      integer(kind = kint) :: i, j, nnod, inum, iorg, idest
!
!
      do i = 1, IO_itp_dest%num_org_domain
        if (IO_itp_dest%id_org_domain(i) .eq. n_org_rank) then
          itp_org%num_dest_domain = itp_org%num_dest_domain + 1
!
          do j = 1, 4
            istack_org_para_type(4*(itp_org%num_dest_domain-1)+j)       &
     &       = istack_org_para_type(4*(itp_org%num_dest_domain-1)+j-1)  &
     &            + istack_table_wtype_dest_IO(4*(i-1)+j)               &
     &            - istack_table_wtype_dest_IO(4*(i-1)+j-1)
          end do
!
          nnod = itp_org%istack_nod_tbl_org(itp_org%num_dest_domain)    &
     &        - itp_org%istack_nod_tbl_org(itp_org%num_dest_domain-1)
          do inum = 1, nnod
            iorg                                                        &
     &       =  itp_org%istack_nod_tbl_org(itp_org%num_dest_domain-1)   &
     &        + inum
            idest = IO_itp_dest%istack_nod_tbl_dest(i-1) + inum
!
            itp_org%inod_itp_send(iorg) =      iorg
            itp_org%inod_gl_dest_4_org(iorg)                            &
     &            = inod_global_dest_IO(idest)
            itp_org%iele_org_4_org(iorg) =     iele_orgin_IO(idest)
            itp_org%itype_inter_org(iorg)                               &
     &            = itype_inter_dest_IO(idest)
            itp_org%coef_inter_org(iorg,1:3)                            &
     &            = coef_inter_dest_IO(idest,1:3)
          end do
        end if
      end do
!
      call deallocate_itp_coefs_dst_IO
      call deallocate_itp_nod_dst_IO
      call deallocate_itp_num_dst_IO
!
      end subroutine set_interpolation_4_orgin
!
!-----------------------------------------------------------------------
!
      end module set_itp_destIO_2_org
