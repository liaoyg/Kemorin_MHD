!djds_matrix_calcs_11
!      module djds_matrix_calcs_11
!
!     Written by Kemorin
!
!   Ordering by Multicolor
!       subroutine change_order_2_solve_bx1(NP, PEsmpTOT, STACKmcG,     &
!     &           NtoO, B, X, W2)
!   Returning from Multicolor
!       subroutine back_2_original_order_bx1(NP, NtoO, B, X, W2)
!
!C +-------------+
!C | {S} = {X}   |
!C +-------------+
!      subroutine copy_internal_vect_1_smp(NP, PEsmpTOT, STACKmcG, S, X)
!
!C +-------------+
!C | {S}= [A]{X} |
!C +-------------+
!       subroutine cal_matvec_11                                        &
!     &           (NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT, PEsmpTOT, &
!     &            STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L, OtoN_U,     &
!     &            NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,         &
!     &            S, X, W3)
!       subroutine cal_matvec_3x11                                      &
!     &           (NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT, PEsmpTOT, &
!     &            STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L, OtoN_U,     &
!     &            NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,         &
!     &            S1, S2, S3, X1, X2, X3, W9)
!
!C +-------------- ----+
!C | {S}= {B} - [A]{X} |
!C +-------------------+
!       subroutine subtruct_matvec_11                                   &
!     &           (NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT, PEsmpTOT, &
!     &            STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L, OtoN_U,     &
!     &            NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,         &
!     &            S, B, X )
!       subroutine subtruct_matvec_3x11                                 &
!     &           (NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT, PEsmpTOT, &
!     &            STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L, OtoN_U,     &
!     &            NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,         &
!     &            S1, S2, S3, B1, B2, B3, X1, X2, X3, W9)
!
      module djds_matrix_calcs_11
!
      use m_precision
!
      use cal_4_diagonal_11
      use cal_4_lower_11
      use cal_4_upper_11
      use ordering_by_o2nl_11
      use ordering_by_l2u_o2nu_11
      use ordering_by_new2old_U_11
      use order_vect_4_solver_11
!
      implicit none
!
       integer(kind = kint), parameter :: IWK1 = 1, IWK2 = 2, IWK3 = 3
       integer(kind = kint), parameter :: IWK11= 1, IWK12= 4, IWK13= 7
       integer(kind = kint), parameter :: IWK21= 2, IWK22= 5, IWK23= 8
       integer(kind = kint), parameter :: IWK31= 3, IWK32= 6, IWK33= 9
       private :: IWK1, IWK2, IWK3
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
       subroutine change_order_2_solve_bx1(NP, PEsmpTOT, STACKmcG,      &
     &           NtoO, B, X, W2)
!
       integer(kind = kint), intent(in) :: NP, PEsmpTOT
       integer(kind = kint), intent(in) :: STACKmcG(0:PEsmpTOT)
       integer(kind = kint), intent(in) :: NtoO(NP)
       real(kind = kreal), intent(inout) :: B(NP)
       real(kind = kreal), intent(inout) :: X(NP)
       real(kind = kreal), intent(inout) :: W2(NP,2)
!
!
      if(NP .le. 0) return
!
       call s_change_order_2_solve_bx1(NP, PEsmpTOT, STACKmcG,          &
     &     NtoO, B, X, W2(1,IWK1), W2(1,IWK2) )
!
      end subroutine change_order_2_solve_bx1
!
!  ---------------------------------------------------------------------
!
      subroutine back_2_original_order_bx1(NP, NtoO, B, X, W2)
!
      integer(kind = kint), intent(in) :: NP
      integer(kind = kint), intent(in) :: NtoO(NP)
      real(kind = kreal), intent(inout) :: B(NP)
      real(kind = kreal), intent(inout) :: X(NP)
       real(kind = kreal), intent(inout) :: W2(NP,2)
!
!
      if(NP .le. 0) return
!
      call s_back_2_original_order_bx1(NP, NtoO, B, X,                  &
     &    W2(1,IWK1), W2(1,IWK2) )
!
      end subroutine back_2_original_order_bx1
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine copy_internal_vect_1_smp(NP, PEsmpTOT, STACKmcG, S, X)
!
       integer(kind = kint), intent(in) :: NP, PEsmpTOT
       integer(kind = kint), intent(in) :: STACKmcG(0:PEsmpTOT)
       real(kind = kreal), intent(in) :: X(NP)
       real(kind = kreal), intent(inout) :: S(NP)
!
       integer (kind = kint) :: ip, iS, iE
!
!
      if(NP .le. 0) return
!
!$omp parallel do private(iS,iE,ip)
       do ip= 1, PEsmpTOT
         iS= STACKmcG(ip-1) + 1
         iE= STACKmcG(ip  )
         S(iS:iE) = X(iS:iE)
       end do
!$omp end parallel do
!
      end subroutine copy_internal_vect_1_smp
!
!  ---------------------------------------------------------------------
!
       subroutine cal_matvec_11                                         &
     &           (NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT, PEsmpTOT,  &
     &            STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L, OtoN_U,      &
     &            NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,          &
     &            S, X, W3)
!
       integer(kind = kint), intent(in) :: NP, PEsmpTOT
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
       real(kind = kreal), intent(in) :: X(NP)
       real(kind = kreal), intent(in) :: D(NP)
       real(kind = kreal), intent(in) :: AL(NPL)
       real(kind = kreal), intent(in) :: AU(NPU)
!
       real(kind = kreal), intent(inout) :: S(NP)
       real(kind = kreal), intent(inout) :: W3(NP,3)
!
!
      if(NP .le. 0) return
!
      call set_by_diagonal_11(NP, PEsmpTOT, STACKmcG,                   &
     &    W3(1,IWK3), X, D)

      call ordering_1x2_by_old2new_L(NP, PEsmpTOT, STACKmcG,            &
     &    OtoN_L, W3(1,IWK1), W3(1,IWK2), X, W3(1,IWK3) )

      call add_lower_11(NP, NL, NPL, PEsmpTOT, NVECT,                   &
     &    npLX1, STACKmc, NLhyp, INL, IAL, W3(1,IWK2), AL, W3(1,IWK1) )

      call ordering_1x1_l2u_o2n_u(NP, OtoN_U, LtoU,                     &
     &    W3(1,IWK1), W3(1,IWK3), X, W3(1,IWK2) )

      call add_upper_11(NP, NU, NPU, PEsmpTOT, NVECT, npUX1,            &
     &    STACKmc, NUhyp, INU, IAU, W3(1,IWK3), AU, W3(1,IWK1) )

      call ordering_1x1_by_new2old_U(NP, PEsmpTOT, STACKmcG,            &
     &    NtoO_U, S, W3(1,IWK3) )
!
       end subroutine cal_matvec_11
!
!  ---------------------------------------------------------------------
!
       subroutine cal_matvec_3x11                                       &
     &           (NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT, PEsmpTOT,  &
     &            STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L, OtoN_U,      &
     &            NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,          &
     &            S1, S2, S3, X1, X2, X3, W9)
!
       integer(kind = kint), intent(in) :: NP, PEsmpTOT
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
       real(kind = kreal), intent(in) :: D(NP)
       real(kind = kreal), intent(in) :: AL(NPL)
       real(kind = kreal), intent(in) :: AU(NPU)
       real(kind = kreal), intent(in) :: X1(NP), X2(NP), X3(NP)
!
       real(kind = kreal), intent(inout) :: S1(NP), S2(NP), S3(NP)
       real(kind = kreal), intent(inout) :: W9(NP,9)
!
!
      if(NP .le. 0) return
!
      call set_by_diagonal_3x11(NP, PEsmpTOT, STACKmcG,                 &
     &    W9(1,IWK13), W9(1,IWK23), W9(1,IWK33), X1, X2, X3, D)

      call ordering_1x6_by_old2new_L(NP, PEsmpTOT, STACKmcG, OtoN_L,    &
     &    W9(1,IWK11), W9(1,IWK21), W9(1,IWK31),                        &
     &    W9(1,IWK12), W9(1,IWK22), W9(1,IWK32), X1, X2, X3,            &
     &    W9(1,IWK13), W9(1,IWK23), W9(1,IWK33) )

      call add_lower_3x11(NP, NL, NPL, PEsmpTOT, NVECT,                 &
     &    npLX1, STACKmc, NLhyp, INL, IAL, W9(1,IWK12), W9(1,IWK22),    &
     &    W9(1,IWK32), AL, W9(1,IWK11), W9(1,IWK21), W9(1,IWK31) )

      call ordering_1x3_l2u_o2n_u(NP, OtoN_U, LtoU,                     &
     &    W9(1,IWK11), W9(1,IWK21), W9(1,IWK31),                        &
     &    W9(1,IWK13), W9(1,IWK23), W9(1,IWK33), X1, X2, X3,            &
     &    W9(1,IWK12), W9(1,IWK22), W9(1,IWK32) )

      call add_upper_3x11(NP, NU, NPU, PEsmpTOT, NVECT, npUX1,          &
     &    STACKmc, NUhyp, INU, IAU, W9(1,IWK13), W9(1,IWK23),           &
     &    W9(1,IWK33), AU, W9(1,IWK11), W9(1,IWK21), W9(1,IWK31) )

      call ordering_1x3_by_new2old_U(NP, PEsmpTOT, STACKmcG,            &
     &    NtoO_U, S1, S2, S3, W9(1,IWK13), W9(1,IWK23), W9(1,IWK33) )
!
       end subroutine cal_matvec_3x11
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
       subroutine subtruct_matvec_11                                    &
     &           (NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT, PEsmpTOT,  &
     &            STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L, OtoN_U,      &
     &            NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,          &
     &            S, B, X, W3)
!
       integer(kind = kint), intent(in) :: NP, PEsmpTOT
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
       real(kind = kreal), intent(in) :: B(NP), X(NP)
       real(kind = kreal), intent(in) :: D(NP)
       real(kind = kreal), intent(in) :: AL(NPL)
       real(kind = kreal), intent(in) :: AU(NPU)
!
       real(kind = kreal), intent(inout) :: S(NP)
       real(kind = kreal), intent(inout) :: W3(NP,3)
!
!
      if(NP .le. 0) return
!
       call subtract_diagonal_11(NP, PEsmpTOT, STACKmcG,                &
     &     W3(1,IWK3), B, X, D)

       call ordering_1x2_by_old2new_L(NP, PEsmpTOT, STACKmcG, OtoN_L,   &
     &     W3(1,IWK1), W3(1,IWK2), X, W3(1,IWK3) )

       call subtract_lower_11(NP, NL, NPL, PEsmpTOT, NVECT, npLX1,      &
     &     STACKmc, NLhyp, INL, IAL, W3(1,IWK2), AL, W3(1,IWK1) )

       call ordering_1x1_l2u_o2n_u(NP, OtoN_U, LtoU,                    &
     &     W3(1,IWK1), W3(1,IWK3), X, W3(1,IWK2) )

       call subtract_upper_11(NP, NU, NPU, PEsmpTOT, NVECT, npUX1,      &
     &     STACKmc, NUhyp, INU, IAU, W3(1,IWK3), AU, W3(1,IWK1) )

       call ordering_1x1_by_new2old_U(NP, PEsmpTOT, STACKmcG,           &
     &     NtoO_U, S, W3(1,IWK3) )
!
       end subroutine subtruct_matvec_11
!
!  ---------------------------------------------------------------------
!
       subroutine subtruct_matvec_3x11                                  &
     &           (NP, NL, NU, NPL, NPU, npLX1, npUX1, NVECT, PEsmpTOT,  &
     &            STACKmcG, STACKmc, NLhyp, NUhyp, OtoN_L, OtoN_U,      &
     &            NtoO_U, LtoU, INL, INU, IAL, IAU, D, AL, AU,          &
     &            S1, S2, S3, B1, B2, B3, X1, X2, X3, W9)
!
       integer(kind = kint), intent(in) :: NP, PEsmpTOT
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
       real(kind = kreal), intent(in) :: D(NP)
       real(kind = kreal), intent(in) :: AL(NPL)
       real(kind = kreal), intent(in) :: AU(NPU)
       real(kind = kreal), intent(in) :: B1(NP), B2(NP), B3(NP)
       real(kind = kreal), intent(in) :: X1(NP), X2(NP), X3(NP)
!
       real(kind = kreal), intent(inout) :: S1(NP), S2(NP), S3(NP)
       real(kind = kreal), intent(inout) :: W9(NP,9)
!
!
      if(NP .le. 0) return
!
       call subtract_diagonal_3x11(NP, PEsmpTOT, STACKmcG,              &
     &     W9(1,IWK13), W9(1,IWK23), W9(1,IWK33), B1, B2, B3,           &
     &     X1, X2, X3, D)

       call ordering_1x6_by_old2new_L(NP, PEsmpTOT, STACKmcG, OtoN_L,   &
     &    W9(1,IWK11), W9(1,IWK21), W9(1,IWK31),                        &
     &    W9(1,IWK12), W9(1,IWK22), W9(1,IWK32), X1, X2, X3,            &
     &    W9(1,IWK13), W9(1,IWK23), W9(1,IWK33) )

       call subtract_lower_3x11(NP, NL, NPL, PEsmpTOT, NVECT, npLX1,    &
     &     STACKmc, NLhyp, INL, IAL, W9(1,IWK12), W9(1,IWK22),          &
     &     W9(1,IWK32), AL, W9(1,IWK11), W9(1,IWK21), W9(1,IWK31) )

       call ordering_1x3_l2u_o2n_u(NP, OtoN_U, LtoU,                    &
     &    W9(1,IWK11), W9(1,IWK21), W9(1,IWK31),                        &
     &    W9(1,IWK13), W9(1,IWK23), W9(1,IWK33), X1, X2, X3,            &
     &    W9(1,IWK12), W9(1,IWK22), W9(1,IWK32) )

       call subtract_upper_3x11(NP, NU, NPU, PEsmpTOT, NVECT, npUX1,    &
     &     STACKmc, NUhyp, INU, IAU, W9(1,IWK13), W9(1,IWK23),          &
     &     W9(1,IWK33), AU, W9(1,IWK11), W9(1,IWK21), W9(1,IWK31) )

       call ordering_1x3_by_new2old_U(NP, PEsmpTOT, STACKmcG,           &
     &    NtoO_U, S1, S2, S3, W9(1,IWK13), W9(1,IWK23), W9(1,IWK33) )
!
       end subroutine subtruct_matvec_3x11
!
!  ---------------------------------------------------------------------
!
      end module djds_matrix_calcs_11
