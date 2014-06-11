!>@file   calypso_SR_rev_int.f90
!!@brief  module calypso_SR_rev_int
!!
!!@author H. Matsui
!!@date Programmed in March, 2013
!
!>@brief  Integer components data communication
!!@n      using reverse import table
!!
!!@verbatim
!!      subroutine calypso_send_recv_rev_int(nnod_org, nnod_new,        &
!!     &                           npe_send, isend_self,                &
!!     &                           id_pe_send, istack_send, inod_export,&
!!     &                           npe_recv, irecv_self,                &
!!     &                           id_pe_recv, istack_recv, irev_import,&
!!     &                           iX_org, iX_new)
!!@endverbatim
!!
!!@n @param  nnod_org    Number of data points for origin
!!@n @param  nnod_new    Number of components for destination
!!@n
!!@n @param  npe_send    Number of processses to send
!!@n @param  isend_self  Integer flag to copy within own process
!!@n @param  id_pe_send(npe_send)      Process ID to send
!!@n @param  istack_send(0:npe_send)
!!                    End points of send buffer for each process
!!@n @param  inod_export(istack_send(npe_send))
!!                    local node ID to copy in send buffer
!!@n
!!@n @param  npe_recv    Number of processses to receive
!!@n @param  irecv_self  Integer flag to copy within own process
!!@n @param  id_pe_recv(npe_send)      Process ID to receive
!!@n @param  istack_recv(0:npe_send)
!!                    End points of receive buffer for each process
!!@n @param  irev_import(nnod_new)
!!                    import buffer ID for each data point
!!@n
!!@n @param  iX_org(nnod_org)   Send data
!!@n @param  iX_new(nnod_new)   Received data
!
      module calypso_SR_rev_int
!
      use m_precision
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine calypso_send_recv_rev_int(nnod_org, nnod_new,          &
     &                            npe_send, isend_self,                 &
     &                            id_pe_send, istack_send, inod_export, &
     &                            npe_recv, irecv_self,                 &
     &                            id_pe_recv, istack_recv, irev_import, &
     &                            iX_org, iX_new)
!
      use calypso_mpi
      use m_solver_SR
      use set_to_send_buffer
      use set_from_recv_buf_rev
!
      integer(kind = kint), intent(in) :: nnod_org
      integer(kind = kint), intent(in) :: nnod_new
!
      integer(kind = kint), intent(in) :: npe_send, isend_self
      integer(kind = kint), intent(in) :: id_pe_send(npe_send)
      integer(kind = kint), intent(in) :: istack_send(0:npe_send)
      integer(kind = kint), intent(in)                                  &
     &                      :: inod_export( istack_send(npe_send) )
!
      integer(kind = kint), intent(in) :: npe_recv, irecv_self
      integer(kind = kint), intent(in) :: id_pe_recv(npe_recv)
      integer(kind = kint), intent(in) :: istack_recv(0:npe_recv)
      integer(kind = kint), intent(in) :: irev_import(nnod_new)
!
      integer (kind=kint), intent(in):: iX_org(nnod_org)
!
      integer (kind=kint), intent(inout):: iX_new(nnod_new)
!
      integer (kind = kint) :: neib, ist, num
      integer (kind = kint) :: i
      integer (kind = kint) :: ncomm_send, ncomm_recv
      integer (kind = kint) :: ist_send, ist_recv
!
!
      ncomm_send = npe_send - isend_self
      ncomm_recv = npe_recv - irecv_self
!
!C-- SEND
!
      call set_to_send_buf_int(nnod_org,                                &
     &    istack_send(npe_send), inod_export, iX_org, iWS)
!C
      do neib = 1, ncomm_send
        ist= istack_send(neib-1) + 1
        num  = istack_send(neib  ) - istack_send(neib-1)
        call MPI_ISEND(iWS(ist), num, CALYPSO_INTEGER,                  &
     &      id_pe_send(neib), 0, CALYPSO_COMM, req1(neib), ierr_MPI)
      end do
!C
!C-- RECEIVE
      if(ncomm_recv .gt. 0) then
        do neib= 1, ncomm_recv
          ist= istack_recv(neib-1) + 1
          num  = istack_recv(neib  ) - istack_recv(neib-1)
          call MPI_IRECV(iWR(ist), num, CALYPSO_INTEGER,                &
     &        id_pe_recv(neib), 0, CALYPSO_COMM, req2(neib), ierr_MPI)
        end do
!
        call MPI_WAITALL (ncomm_recv, req2, sta2, ierr_MPI)
      end if
!
      if (isend_self .eq. 1) then
        ist_send= istack_send(npe_send-1)
        ist_recv= istack_recv(npe_recv-1)
        num  =   istack_send(npe_send  ) - istack_send(npe_send-1) 
!$omp parallel do
        do i = 1, num
          iWR(ist_recv+i) = iWS(ist_send+i)
        end do
!$omp end parallel do
      end if
!
      call set_from_recv_buf_rev_int(nnod_new,                          &
     &    istack_recv(npe_recv), irev_import, iWR, iX_new)
!
      end subroutine calypso_send_recv_rev_int
!
! ----------------------------------------------------------------------
!
      subroutine calypso_send_recv_rev_int8(nnod_org, nnod_new,         &
     &                            npe_send, isend_self,                 &
     &                            id_pe_send, istack_send, inod_export, &
     &                            npe_recv, irecv_self,                 &
     &                            id_pe_recv, istack_recv, irev_import, &
     &                            i8X_org, i8X_new)
!
      use calypso_mpi
      use m_solver_SR
      use set_to_send_buffer
      use set_from_recv_buf_rev
!
      integer(kind = kint), intent(in) :: nnod_org
      integer(kind = kint), intent(in) :: nnod_new
!
      integer(kind = kint), intent(in) :: npe_send, isend_self
      integer(kind = kint), intent(in) :: id_pe_send(npe_send)
      integer(kind = kint), intent(in) :: istack_send(0:npe_send)
      integer(kind = kint), intent(in)                                  &
     &                      :: inod_export( istack_send(npe_send) )
!
      integer(kind = kint), intent(in) :: npe_recv, irecv_self
      integer(kind = kint), intent(in) :: id_pe_recv(npe_recv)
      integer(kind = kint), intent(in) :: istack_recv(0:npe_recv)
      integer(kind = kint), intent(in) :: irev_import(nnod_new)
!
      integer (kind=kint_d), intent(in):: i8X_org(nnod_org)
!
      integer (kind=kint_d), intent(inout):: i8X_new(nnod_new)
!
      integer (kind = kint) :: neib, ist, num
      integer (kind = kint) :: i
      integer (kind = kint) :: ncomm_send, ncomm_recv
      integer (kind = kint) :: ist_send, ist_recv
!
!
      ncomm_send = npe_send - isend_self
      ncomm_recv = npe_recv - irecv_self
!
!C-- SEND
!
      call set_to_send_buf_int8(nnod_org,                               &
     &    istack_send(npe_send), inod_export, i8X_org, i8WS)
!C
      do neib = 1, ncomm_send
        ist= istack_send(neib-1) + 1
        num  = istack_send(neib  ) - istack_send(neib-1)
        call MPI_ISEND(i8WS(ist), num, CALYPSO_INTEGER8,                &
     &      id_pe_send(neib), 0, CALYPSO_COMM, req1(neib), ierr_MPI)
      end do
!C
!C-- RECEIVE
      if(ncomm_recv .gt. 0) then
        do neib= 1, ncomm_recv
          ist= istack_recv(neib-1) + 1
          num  = istack_recv(neib  ) - istack_recv(neib-1)
          call MPI_IRECV(i8WR(ist), num, CALYPSO_INTEGER8,              &
     &        id_pe_recv(neib), 0, CALYPSO_COMM, req2(neib), ierr_MPI)
        end do
!
        call MPI_WAITALL (ncomm_recv, req2, sta2, ierr_MPI)
      end if
!
      if (isend_self .eq. 1) then
        ist_send= istack_send(npe_send-1)
        ist_recv= istack_recv(npe_recv-1)
        num  =   istack_send(npe_send  ) - istack_send(npe_send-1) 
!$omp parallel do
        do i = 1, num
          i8WR(ist_recv+i) = i8WS(ist_send+i)
        end do
!$omp end parallel do
      end if
!
      call set_from_recv_buf_rev_int8(nnod_new,                         &
     &    istack_recv(npe_recv), irev_import, i8WR, i8X_new)
!
      end subroutine calypso_send_recv_rev_int8
!
! ----------------------------------------------------------------------
!
      end module calypso_SR_rev_int
