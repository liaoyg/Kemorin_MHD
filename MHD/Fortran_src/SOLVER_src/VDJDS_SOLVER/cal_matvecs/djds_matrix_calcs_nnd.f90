!
!      module djds_matrix_calcs_nnd
!
!     Written by Kemorin
!
!!C +-------------+
!!C | {S}= [A]{X} |
!!C +-------------+
!!
!!       subroutine cal_matvec_nnd                                      &
!!     &           (NP, NB, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,      &
!!     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,  &
!!     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,&
!!     &            S, X, W3)
!!       subroutine cal_matvec_3xnnd                                    &
!!     &           (NP, NB, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,      &
!!     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,  &
!!     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,&
!!     &            S1, S2, S3, X1, X2, X3, W9)
!!
!!C +-------------- ----+
!!C | {S}= {B} - [A]{X} |
!!C +-------------------+
!!
!!       subroutine subtruct_matvec_nnd                                 &
!!     &           (NP, NB, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,      &
!!     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,  &
!!     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,&
!!     &            S, B, X, W3)
!!       subroutine subtruct_matvec_3xnnd                               &
!!     &           (NP, NB, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,      &
!!     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,  &
!!     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,&
!!     &            S1, S2, S3, X1, X2, X3, W9)
!
      module djds_matrix_calcs_nnd
!
      use m_precision
!
      use cal_4_diagonal_nnd
      use cal_4_lower_nnd
      use cal_4_upper_nnd
      use ordering_by_o2nl_nn
      use ordering_by_l2u_o2nu_nn
      use ordering_by_new2old_U_nn
!
      implicit none
!
       integer(kind = kint), parameter :: IWK1 = 1, IWK2 = 2, IWK3 = 3
       integer(kind = kint), parameter :: IWK11= 1, IWK12= 4, IWK13= 7
       integer(kind = kint), parameter :: IWK21= 2, IWK22= 5, IWK23= 8
       integer(kind = kint), parameter :: IWK31= 3, IWK32= 6, IWK33= 9
       private :: IWK1,  IWK2,  IWK3
       private :: IWK11, IWK12, IWK13
       private :: IWK21, IWK22, IWK23
       private :: IWK31, IWK32, IWK33
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
       subroutine cal_matvec_nnd                                        &
     &           (NP, NB, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,        &
     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,    &
     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,  &
     &            S, X, W3)
!
       integer(kind = kint), intent(in) :: NP, NB, PEsmpTOT
       integer(kind = kint), intent(in) :: NL, NU, NPL, NPU
       integer(kind = kint), intent(in) :: NVECT, npLX1, npUX1
       integer(kind = kint), intent(in) :: NLhyp(NVECT), NUhyp(NVECT)
       integer(kind = kint), intent(in) :: STACKmcG(0:PEsmpTOT)
       integer(kind = kint), intent(in) :: STACKmc(0:PEsmpTOT*NVECT)
       integer(kind = kint), intent(in) :: OtoN_L(NP), OtoN_U(NP)
       integer(kind = kint), intent(in) :: NtoO_U(NP), LtoU(NP)
       integer(kind = kint), intent(in) :: INL(0:NL*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: INU(0:NU*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: IAL(NPL), IAU(NPU)
!
       real(kind = kreal), intent(in) :: X(NB*NP), D(NB*NB*NP)
       real(kind = kreal), intent(in) :: AL(NB*NB*NPL), AU(NB*NB*NPU)
!
       real(kind = kreal), intent(inout) :: S(NB*NP)
       real(kind = kreal), intent(inout) :: W3(NB*NP,3)
!
!
      if(NP .le. 0) return
!
      call set_by_diagonal_nnd(NP, NB, PEsmpTOT, STACKmcG,              &
     &    W3(1,IWK3), X, D)

      call ordering_nx2_by_old2new_L(NP, NB, PEsmpTOT, STACKmcG,        &
     &    OtoN_L, W3(1,IWK1), W3(1,IWK2), X, W3(1,IWK3) )

      call add_lower_nnd(NP, NB, NL, NPL, PEsmpTOT, NVECT,              &
     &    npLX1, STACKmc, NLhyp, INL, IAL, W3(1,IWK2), AL, W3(1,IWK1) )

      call ordering_nx1_l2u_o2n_u(NP, NB, OtoN_U, LtoU,                 &
     &    W3(1,IWK1), W3(1,IWK3), X, W3(1,IWK2) )

      call add_upper_nnd(NP, NB, NU, NPU, PEsmpTOT, NVECT, npUX1,       &
     &    STACKmc, NUhyp, INU, IAU, W3(1,IWK3), AU, W3(1,IWK1) )

      call ordering_nx1_by_new2old_U(NP, NB, PEsmpTOT, STACKmcG,        &
     &    NtoO_U, S, W3(1,IWK3) )
!
       end subroutine cal_matvec_nnd
!
!  ---------------------------------------------------------------------
!
       subroutine cal_matvec_3xnnd                                      &
     &           (NP, NB, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,        &
     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,    &
     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,  &
     &            S1, S2, S3, X1, X2, X3, W9)
!
       integer(kind = kint), intent(in) :: NP, NB, PEsmpTOT
       integer(kind = kint), intent(in) :: NL, NU, NPL, NPU
       integer(kind = kint), intent(in) :: NVECT, npLX1, npUX1
       integer(kind = kint), intent(in) :: NLhyp(NVECT), NUhyp(NVECT)
       integer(kind = kint), intent(in) :: STACKmcG(0:PEsmpTOT)
       integer(kind = kint), intent(in) :: STACKmc(0:PEsmpTOT*NVECT)
       integer(kind = kint), intent(in) :: OtoN_L(NP), OtoN_U(NP)
       integer(kind = kint), intent(in) :: NtoO_U(NP), LtoU(NP)
       integer(kind = kint), intent(in) :: INL(0:NL*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: INU(0:NU*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: IAL(NPL), IAU(NPU)
!
       real(kind = kreal), intent(in) :: D(NB*NB*NP)
       real(kind = kreal), intent(in) :: AL(NB*NB*NPL), AU(NB*NB*NPU)
       real(kind = kreal), intent(in) :: X1(NB*NP), X2(NB*NP)
       real(kind = kreal), intent(in) :: X3(NB*NP)
!
       real(kind = kreal), intent(inout) :: S1(NB*NP), S2(NB*NP)
       real(kind = kreal), intent(inout) :: S3(NB*NP)
       real(kind = kreal), intent(inout) :: W9(NB*NP,9)
!
!
      if(NP .le. 0) return
!
      call set_by_diagonal_3xnnd(NP, NB, PEsmpTOT, STACKmcG,            &
     &    W9(1,IWK13), W9(1,IWK23), W9(1,IWK33), X1, X2, X3, D)

      call ordering_nx6_by_old2new_L(NP, NB, PEsmpTOT, STACKmcG,        &
     &    OtoN_L, W9(1,IWK11), W9(1,IWK21), W9(1,IWK31),                &
     &    W9(1,IWK12), W9(1,IWK22), W9(1,IWK32), X1, X2, X3,            &
     &    W9(1,IWK13), W9(1,IWK23), W9(1,IWK33) )

      call add_lower_3xnnd(NP, NB, NL, NPL, PEsmpTOT, NVECT,            &
     &    npLX1, STACKmc, NLhyp, INL, IAL, W9(1,IWK12), W9(1,IWK22),    &
     &    W9(1,IWK32), AL, W9(1,IWK11), W9(1,IWK21), W9(1,IWK31) )

       call ordering_nx3_l2u_o2n_u(NP, NB, OtoN_U, LtoU,                &
     &    W9(1,IWK11), W9(1,IWK21), W9(1,IWK31),                        &
     &    W9(1,IWK13), W9(1,IWK23), W9(1,IWK33), X1, X2, X3,            &
     &    W9(1,IWK12), W9(1,IWK22), W9(1,IWK32) )

      call add_upper_3xnnd(NP, NB, NU, NPU, PEsmpTOT, NVECT, npUX1,     &
     &    STACKmc, NUhyp, INU, IAU, W9(1,IWK13), W9(1,IWK23),           &
     &    W9(1,IWK33), AU, W9(1,IWK11), W9(1,IWK21), W9(1,IWK31) )

      call ordering_nx3_by_new2old_U(NP, NB, PEsmpTOT, STACKmcG,        &
     &    NtoO_U, S1, S2, S3, W9(1,IWK13), W9(1,IWK23), W9(1,IWK33) )

!
       end subroutine cal_matvec_3xnnd
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
       subroutine subtruct_matvec_nnd                                   &
     &           (NP, NB, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,        &
     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,    &
     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,  &
     &            S, B, X, W3)
!
       integer(kind = kint), intent(in) :: NP, NB, PEsmpTOT
       integer(kind = kint), intent(in) :: NL, NU, NPL, NPU
       integer(kind = kint), intent(in) :: NVECT, npLX1, npUX1
       integer(kind = kint), intent(in) :: NLhyp(NVECT), NUhyp(NVECT)
       integer(kind = kint), intent(in) :: STACKmcG(0:PEsmpTOT)
       integer(kind = kint), intent(in) :: STACKmc(0:PEsmpTOT*NVECT)
       integer(kind = kint), intent(in) :: OtoN_L(NP), OtoN_U(NP)
       integer(kind = kint), intent(in) :: NtoO_U(NP), LtoU(NP)
       integer(kind = kint), intent(in) :: INL(0:NL*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: INU(0:NU*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: IAL(NPL), IAU(NPU)
!
       real(kind = kreal), intent(in) :: B(NB*NP), X(NB*NP)
       real(kind = kreal), intent(in) :: D(NB*NB*NP)
       real(kind = kreal), intent(in) :: AL(NB*NB*NPL), AU(NB*NB*NPU)
!
       real(kind = kreal), intent(inout) :: S(NB*NP)
       real(kind = kreal), intent(inout) :: W3(NB*NP,3)
!
!
      if(NP .le. 0) return
!
       call subtract_diagonal_nnd(NP, NB, PEsmpTOT, STACKmcG,           &
     &     W3(1,IWK3), B, X, D)

       call ordering_nx2_by_old2new_L(NP, NB, PEsmpTOT, STACKmcG,       &
     &     OtoN_L, W3(1,IWK1), W3(1,IWK2), X, W3(1,IWK3) )

       call subtract_lower_nnd(NP, NB, NL, NPL, PEsmpTOT, NVECT, npLX1, &
     &     STACKmc, NLhyp, INL, IAL, W3(1,IWK2), AL, W3(1,IWK1) )

       call ordering_nx1_l2u_o2n_u(NP, NB, OtoN_U, LtoU,                &
     &     W3(1,IWK1), W3(1,IWK3), X, W3(1,IWK2) )

       call subtract_upper_nnd(NP, NB, NU, NPU, PEsmpTOT, NVECT, npUX1, &
     &     STACKmc, NUhyp, INU, IAU, W3(1,IWK3), AU, W3(1,IWK1) )

       call ordering_nx1_by_new2old_U(NP, NB, PEsmpTOT, STACKmcG,       &
     &     NtoO_U, S, W3(1,IWK3) )
!
       end subroutine subtruct_matvec_nnd
!
!  ---------------------------------------------------------------------
!
       subroutine subtruct_matvec_3xnnd                                 &
     &           (NP, NB, NL, NU, NPL, NPU, npLX1, npUX1, NVECT,        &
     &            PEsmpTOT, STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L,    &
     &            OtoN_U, NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,  &
     &            S1, S2, S3, B1, B2, B3, X1, X2, X3, W9)
!
       integer(kind = kint), intent(in) :: NP, NB, PEsmpTOT
       integer(kind = kint), intent(in) :: NL, NU, NPL, NPU
       integer(kind = kint), intent(in) :: NVECT, npLX1, npUX1
       integer(kind = kint), intent(in) :: NLhyp(NVECT), NUhyp(NVECT)
       integer(kind = kint), intent(in) :: STACKmcG(0:PEsmpTOT)
       integer(kind = kint), intent(in) :: STACKmc(0:PEsmpTOT*NVECT)
       integer(kind = kint), intent(in) :: OtoN_L(NP), OtoN_U(NP)
       integer(kind = kint), intent(in) :: NtoO_U(NP), LtoU(NP)
       integer(kind = kint), intent(in) :: INL(0:NL*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: INU(0:NU*NVECT*PEsmpTOT)
       integer(kind = kint), intent(in) :: IAL(NPL), IAU(NPU)
!
       real(kind = kreal), intent(in) :: D(NB*NB*NP)
       real(kind = kreal), intent(in) :: AL(NB*NB*NPL), AU(NB*NB*NPU)
       real(kind = kreal), intent(in) :: B1(NB*NP), X1(NB*NP)
       real(kind = kreal), intent(in) :: B2(NB*NP), X2(NB*NP)
       real(kind = kreal), intent(in) :: B3(NB*NP), X3(NB*NP)
!
       real(kind = kreal), intent(inout) :: S1(NB*NP), S2(NB*NP)
       real(kind = kreal), intent(inout) :: S3(NB*NP)
       real(kind = kreal), intent(inout) :: W9(NB*NP,9)
!
!
      if(NP .le. 0) return
!
       call subtract_diagonal_3xnnd(NP, NB, PEsmpTOT, STACKmcG,         &
     &     W9(1,IWK13), W9(1,IWK23), W9(1,IWK33), B1, B2, B3,           &
     &     X1, X2, X3, D)

      call ordering_nx6_by_old2new_L(NP, NB, PEsmpTOT, STACKmcG,        &
     &    OtoN_L, W9(1,IWK11), W9(1,IWK21), W9(1,IWK31),                &
     &    W9(1,IWK12), W9(1,IWK22), W9(1,IWK32), X1, X2, X3,            &
     &    W9(1,IWK13), W9(1,IWK23), W9(1,IWK33) )

       call subtract_lower_3xnnd(NP, NB, NL, NPL, PEsmpTOT, NVECT,      &
     &     npLX1, STACKmc, NLhyp, INL, IAL, W9(1,IWK12), W9(1,IWK22),   &
     &     W9(1,IWK32), AL, W9(1,IWK11), W9(1,IWK21), W9(1,IWK31) )

       call ordering_nx3_l2u_o2n_u(NP, NB, OtoN_U, LtoU,                &
     &     W9(1,IWK11), W9(1,IWK21), W9(1,IWK31),                       &
     &     W9(1,IWK13), W9(1,IWK23), W9(1,IWK33), X1, X2, X3,           &
     &     W9(1,IWK12), W9(1,IWK22), W9(1,IWK32) )

       call subtract_upper_3xnnd(NP, NB, NU, NPU, PEsmpTOT, NVECT,      &
     &     npUX1, STACKmc, NUhyp, INU, IAU, W9(1,IWK13), W9(1,IWK23),   &
     &     W9(1,IWK33), AU, W9(1,IWK11), W9(1,IWK21), W9(1,IWK31) )

       call ordering_nx3_by_new2old_U(NP, NB, PEsmpTOT, STACKmcG,       &
     &     NtoO_U, S1, S2, S3, W9(1,IWK13), W9(1,IWK23), W9(1,IWK33) )
!
       end subroutine subtruct_matvec_3xnnd
!
!  ---------------------------------------------------------------------
!
      end module djds_matrix_calcs_nnd
