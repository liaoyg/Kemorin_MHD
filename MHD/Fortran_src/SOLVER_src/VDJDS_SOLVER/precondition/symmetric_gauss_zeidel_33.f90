!symmetric_gauss_zeidel_33.f90
!      module symmetric_gauss_zeidel_33
!
!     Written by Kemorin
!
!!       subroutine weak_sym_gauss_zeidel_33                            &
!!     &           (N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,       &
!!     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,  &
!!     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,   &
!!     &            ALU_U, B, S, W3, NEIBPETOT, NEIBPE,                 &
!!     &            STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT, &
!!     &            niter)
!!
!!       subroutine full_sym_gauss_zeidel_33                            &
!!     &           (N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,       &
!!     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,  &
!!     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,   &
!!     &            ALU_U, B, S, W3, NEIBPETOT, NEIBPE,                 &
!!     &            STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT, &
!!     &            niter)
!!
!!
!!       subroutine weak_sym_gauss_zeidel_3x33                          &
!!     &          (N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,        &
!!     &           PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,   &
!!     &           OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,    &
!!     &           ALU_U, B1, B2, B3, S1, S2, S3, W9, NEIBPETOT, NEIBPE,&
!!     &           STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT,  &
!!     &           niter)
!!
!!       subroutine full_sym_gauss_zeidel_3x33                          &
!!     &           (N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,       &
!!     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,  &
!!     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,   &
!!     &           ALU_U, B1, B2, B3, S1, S2, S3, W9, NEIBPETOT, NEIBPE,&
!!     &            STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT, &
!!     &            niter)
!
      module symmetric_gauss_zeidel_33
!
      use m_precision
!
      use gauss_zeidel_33
      use m_solver_count_time
      use solver_SR_3
!
      implicit none
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
       subroutine weak_sym_gauss_zeidel_33                              &
     &           (N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,         &
     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,    &
     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,     &
     &            ALU_U, B, S, W3, NEIBPETOT, NEIBPE,                   &
     &            STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT,   &
     &            niter)
!
      use calypso_mpi
!
      integer(kind=kint ), intent(in) :: niter
!
      integer(kind=kint ), intent(in) :: NEIBPETOT
      integer(kind=kint ), intent(in) :: NEIBPE(NEIBPETOT)
      integer(kind=kint ), intent(in) :: STACK_IMPORT(0:NEIBPETOT)
      integer(kind=kint ), intent(in)                                   &
     &      :: NOD_IMPORT(STACK_IMPORT(NEIBPETOT))
      integer(kind=kint ), intent(in) :: STACK_EXPORT(0:NEIBPETOT)
      integer(kind=kint ), intent(in)                                   &
     &      :: NOD_EXPORT(STACK_EXPORT(NEIBPETOT)) 
!
       integer(kind = kint), intent(in) :: N, NP, PEsmpTOT
       integer(kind = kint), intent(in) :: NL, NU, NPL, NPU
       integer(kind = kint), intent(in) :: NVECT, npLX1, npUX1
       integer(kind = kint), intent(in) :: NLhyp(NVECT)
       integer(kind = kint), intent(in) :: NUhyp(NVECT)
       integer(kind = kint), intent(in) :: STACKmcG(0:PEsmpTOT)
       integer(kind = kint), intent(in) :: STACKmc(0:PEsmpTOT*NVECT)
       integer(kind = kint), intent(in) :: OtoN_L(NP), OtoN_U(NP)
       integer(kind = kint), intent(in) :: NtoO_U(NP), LtoU(NP)
       integer(kind = kint), intent(in) :: INL(0:NL*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: INU(0:NU*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: IAL(NPL)
       integer(kind = kint), intent(in) :: IAU(NPU)
!
       real(kind = kreal), intent(in) :: ALU_U(9*N)
       real(kind = kreal), intent(in) :: AL(9*NPL)
       real(kind = kreal), intent(in) :: AU(9*NPU)
       real(kind = kreal), intent(in) :: B(3*NP)
!
       real(kind = kreal), intent(inout) :: S(3*NP)
       real(kind = kreal), intent(inout) :: W3(3*NP,3)
!
       integer(kind = kint) :: icou
!
!
      do icou = 1, niter
        call gauss_zeidel_forward_33                                    &
     &           (N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,         &
     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,    &
     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,     &
     &            ALU_U, B, S, W3)
!
        START_TIME= MPI_WTIME()
        call SOLVER_SEND_RECV_3                                         &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, S)
          END_TIME= MPI_WTIME()
          COMMtime = COMMtime + END_TIME - START_TIME
!
        call gauss_zeidel_backward_33                                   &
     &           (N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,         &
     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,    &
     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,     &
     &            ALU_U, B, S, W3)
!
        START_TIME= MPI_WTIME()
        call SOLVER_SEND_RECV_3                                         &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, S)
          END_TIME= MPI_WTIME()
          COMMtime = COMMtime + END_TIME - START_TIME
      end do
!
      end subroutine weak_sym_gauss_zeidel_33
!
!  ---------------------------------------------------------------------
!
       subroutine full_sym_gauss_zeidel_33                              &
     &           (N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,         &
     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,    &
     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,     &
     &            ALU_U, B, S, W3, NEIBPETOT, NEIBPE,                   &
     &            STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT,   &
     &            niter)
!
      use calypso_mpi
!
      integer(kind=kint ), intent(in) :: niter
!
      integer(kind=kint ), intent(in) :: NEIBPETOT
      integer(kind=kint ), intent(in) :: NEIBPE(NEIBPETOT)
      integer(kind=kint ), intent(in) :: STACK_IMPORT(0:NEIBPETOT)
      integer(kind=kint ), intent(in)                                   &
     &      :: NOD_IMPORT(STACK_IMPORT(NEIBPETOT))
      integer(kind=kint ), intent(in) :: STACK_EXPORT(0:NEIBPETOT)
      integer(kind=kint ), intent(in)                                   &
     &      :: NOD_EXPORT(STACK_EXPORT(NEIBPETOT)) 
!
       integer(kind = kint), intent(in) :: N, NP, PEsmpTOT
       integer(kind = kint), intent(in) :: NL, NU, NPL, NPU
       integer(kind = kint), intent(in) :: NVECT, npLX1, npUX1
       integer(kind = kint), intent(in) :: NLhyp(NVECT)
       integer(kind = kint), intent(in) :: NUhyp(NVECT)
       integer(kind = kint), intent(in) :: STACKmcG(0:PEsmpTOT)
       integer(kind = kint), intent(in) :: STACKmc(0:PEsmpTOT*NVECT)
       integer(kind = kint), intent(in) :: OtoN_L(NP), OtoN_U(NP)
       integer(kind = kint), intent(in) :: NtoO_U(NP), LtoU(NP)
       integer(kind = kint), intent(in) :: INL(0:NL*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: INU(0:NU*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: IAL(NPL)
       integer(kind = kint), intent(in) :: IAU(NPU)
!
       real(kind = kreal), intent(in) :: ALU_U(9*N)
       real(kind = kreal), intent(in) :: AL(9*NPL)
       real(kind = kreal), intent(in) :: AU(9*NPU)
       real(kind = kreal), intent(in) :: B(3*NP)
!
       real(kind = kreal), intent(inout) :: S(3*NP)
       real(kind = kreal), intent(inout) :: W3(3*NP,3)
!
       integer(kind = kint) :: icou, iv
!
!
      do icou = 1, niter
        do iv= 1, NVECT
         call gauss_zeidel_each_plane_33                                &
     &           (iv, N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,     &
     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,    &
     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,     &
     &            ALU_U, B, S, W3)
!
!
          START_TIME= MPI_WTIME()
          call SOLVER_SEND_RECV_3                                       &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, S)
           END_TIME= MPI_WTIME()
           COMMtime = COMMtime + END_TIME - START_TIME
        end do
!
        do iv= NVECT, 1, -1
         call gauss_zeidel_each_plane_33                                &
     &           (iv, N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,     &
     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,    &
     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,     &
     &            ALU_U, B, S, W3)
!
!
        START_TIME= MPI_WTIME()
          call SOLVER_SEND_RECV_3                                       &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, S)
          END_TIME= MPI_WTIME()
          COMMtime = COMMtime + END_TIME - START_TIME
        end do
      end do
!
      end subroutine full_sym_gauss_zeidel_33
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
       subroutine weak_sym_gauss_zeidel_3x33                            &
     &          (N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,          &
     &           PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,     &
     &           OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,      &
     &           ALU_U, B1, B2, B3, S1, S2, S3, W9, NEIBPETOT, NEIBPE,  &
     &           STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT,    &
     &           niter)
!
      use calypso_mpi
!
      integer(kind=kint ), intent(in) :: niter
!
      integer(kind=kint ), intent(in) :: NEIBPETOT
      integer(kind=kint ), intent(in) :: NEIBPE(NEIBPETOT)
      integer(kind=kint ), intent(in) :: STACK_IMPORT(0:NEIBPETOT)
      integer(kind=kint ), intent(in)                                   &
     &      :: NOD_IMPORT(STACK_IMPORT(NEIBPETOT))
      integer(kind=kint ), intent(in) :: STACK_EXPORT(0:NEIBPETOT)
      integer(kind=kint ), intent(in)                                   &
     &      :: NOD_EXPORT(STACK_EXPORT(NEIBPETOT)) 
!
       integer(kind = kint), intent(in) :: N, NP, PEsmpTOT
       integer(kind = kint), intent(in) :: NL, NU, NPL, NPU
       integer(kind = kint), intent(in) :: NVECT, npLX1, npUX1
       integer(kind = kint), intent(in) :: NLhyp(NVECT)
       integer(kind = kint), intent(in) :: NUhyp(NVECT)
       integer(kind = kint), intent(in) :: STACKmcG(0:PEsmpTOT)
       integer(kind = kint), intent(in) :: STACKmc(0:PEsmpTOT*NVECT)
       integer(kind = kint), intent(in) :: OtoN_L(NP), OtoN_U(NP)
       integer(kind = kint), intent(in) :: NtoO_U(NP), LtoU(NP)
       integer(kind = kint), intent(in) :: INL(0:NL*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: INU(0:NU*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: IAL(NPL)
       integer(kind = kint), intent(in) :: IAU(NPU)
!
       real(kind = kreal), intent(in) :: ALU_U(9*N)
       real(kind = kreal), intent(in) :: AL(9*NPL)
       real(kind = kreal), intent(in) :: AU(9*NPU)
       real(kind = kreal), intent(in) :: B1(3*NP), B2(3*NP), B3(3*NP)
!
       real(kind = kreal), intent(inout) :: S1(3*NP), S2(3*NP)
       real(kind = kreal), intent(inout) :: S3(3*NP)
       real(kind = kreal), intent(inout) :: W9(3*NP,9)
!
       integer(kind = kint) :: icou
!
!
      do icou = 1, niter
        call gauss_zeidel_forward_3x33                                  &
     &           (N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,         &
     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,    &
     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,     &
     &            ALU_U, B1, B2, B3, S1, S2, S3, W9)
!
        START_TIME= MPI_WTIME()
        call solver_send_recv_3x3                                       &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, S1, S2, S3)
        END_TIME= MPI_WTIME()
        COMMtime = COMMtime + END_TIME - START_TIME
!
        call gauss_zeidel_backward_3x33                                 &
     &           (N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,         &
     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,    &
     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,     &
     &            ALU_U, B1, B2, B3, S1, S2, S3, W9)
!
        START_TIME= MPI_WTIME()
        call solver_send_recv_3x3                                       &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, S1, S2, S3)
        END_TIME= MPI_WTIME()
        COMMtime = COMMtime + END_TIME - START_TIME
      end do
!
       end subroutine weak_sym_gauss_zeidel_3x33
!
!  ---------------------------------------------------------------------
!
       subroutine full_sym_gauss_zeidel_3x33                            &
     &          (N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,          &
     &           PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,     &
     &           OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,      &
     &           ALU_U, B1, B2, B3, S1, S2, S3, W9, NEIBPETOT, NEIBPE,  &
     &           STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT,    &
     &           niter)
!
      use calypso_mpi
!
      integer(kind=kint ), intent(in) :: niter
!
      integer(kind=kint ), intent(in) :: NEIBPETOT
      integer(kind=kint ), intent(in) :: NEIBPE(NEIBPETOT)
      integer(kind=kint ), intent(in) :: STACK_IMPORT(0:NEIBPETOT)
      integer(kind=kint ), intent(in)                                   &
     &      :: NOD_IMPORT(STACK_IMPORT(NEIBPETOT))
      integer(kind=kint ), intent(in) :: STACK_EXPORT(0:NEIBPETOT)
      integer(kind=kint ), intent(in)                                   &
     &      :: NOD_EXPORT(STACK_EXPORT(NEIBPETOT)) 
!
       integer(kind = kint), intent(in) :: N, NP, PEsmpTOT
       integer(kind = kint), intent(in) :: NL, NU, NPL, NPU
       integer(kind = kint), intent(in) :: NVECT, npLX1, npUX1
       integer(kind = kint), intent(in) :: NLhyp(NVECT)
       integer(kind = kint), intent(in) :: NUhyp(NVECT)
       integer(kind = kint), intent(in) :: STACKmcG(0:PEsmpTOT)
       integer(kind = kint), intent(in) :: STACKmc(0:PEsmpTOT*NVECT)
       integer(kind = kint), intent(in) :: OtoN_L(NP), OtoN_U(NP)
       integer(kind = kint), intent(in) :: NtoO_U(NP), LtoU(NP)
       integer(kind = kint), intent(in) :: INL(0:NL*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: INU(0:NU*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: IAL(NPL)
       integer(kind = kint), intent(in) :: IAU(NPU)
!
       real(kind = kreal), intent(in) :: ALU_U(9*N)
       real(kind = kreal), intent(in) :: AL(9*NPL)
       real(kind = kreal), intent(in) :: AU(9*NPU)
       real(kind = kreal), intent(in) :: B1(3*NP), B2(3*NP), B3(3*NP)
!
       real(kind = kreal), intent(inout) :: S1(3*NP), S2(3*NP)
       real(kind = kreal), intent(inout) :: S3(3*NP)
       real(kind = kreal), intent(inout) :: W9(3*NP,9)
!
       integer(kind = kint) :: icou, iv
!
!
      do icou = 1, niter
        do iv= 1, NVECT
          call gauss_zeidel_each_plane_3x33                             &
     &           (iv, N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,     &
     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,    &
     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,     &
     &            ALU_U, B1, B2, B3, S1, S2, S3, W9)
!
          START_TIME= MPI_WTIME()
          call solver_send_recv_3x3                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, S1, S2, S3)
          END_TIME= MPI_WTIME()
          COMMtime = COMMtime + END_TIME - START_TIME
        end do
!
        do iv= NVECT, 1, -1
          call gauss_zeidel_each_plane_3x33                             &
     &           (iv, N, NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,     &
     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,    &
     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, AL, AU,     &
     &            ALU_U, B1, B2, B3, S1, S2, S3, W9)
!
          START_TIME= MPI_WTIME()
          call solver_send_recv_3x3                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, S1, S2, S3)
          END_TIME= MPI_WTIME()
          COMMtime = COMMtime + END_TIME - START_TIME
!
        end do
      end do
!
      end subroutine full_sym_gauss_zeidel_3x33
!
!  ---------------------------------------------------------------------
!
      end module symmetric_gauss_zeidel_33
