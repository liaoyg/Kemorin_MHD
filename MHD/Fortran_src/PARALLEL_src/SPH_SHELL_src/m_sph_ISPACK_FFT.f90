!>@file   m_sph_ISPACK_FFT.f90
!!@brief  module m_sph_ISPACK_FFT
!!
!!@author H. Matsui
!!@date Programmed in 2008
!
!
!>@brief  Fourier transform using ISPACK
!!
!!@verbatim
!!  ---------------------------------------------------------------------
!!
!!      subroutine init_sph_ISPACK(nidx_rtp, maxirt_rtp_smp, ncomp)
!!      subroutine finalize_sph_ISPACK
!!      subroutine verify_sph_ISPACK(nidx_rtp, maxirt_rtp_smp, ncomp)
!! ------------------------------------------------------------------
!! wrapper subroutine for initierize FFT for ISPACK
!! ------------------------------------------------------------------
!!
!!      subroutine sph_FTTRUF_to_send                                   &
!!     &         (nnod_rtp, nidx_rtp, irt_rtp_smp_stack,                &
!!     &          ncomp, n_WS, irev_sr_rtp, X_rtp, WS)
!! ------------------------------------------------------------------
!!
!! wrapper subroutine for forward Fourier transform by ISPACK
!!
!! a_{k} = \frac{2}{Nfft} \sum_{j=0}^{Nfft-1} x_{j} \cos (\frac{2\pi j k}{Nfft})
!! b_{k} = \frac{2}{Nfft} \sum_{j=0}^{Nfft-1} x_{j} \cos (\frac{2\pi j k}{Nfft})
!!
!! a_{0} = \frac{1}{Nfft} \sum_{j=0}^{Nfft-1} x_{j}
!! K = Nfft/2....
!! a_{k} = \frac{1}{Nfft} \sum_{j=0}^{Nfft-1} x_{j} \cos (\frac{2\pi j k}{Nfft})
!!
!! ------------------------------------------------------------------
!!
!!      subroutine sph_FTTRUB_from_recv                                 &
!!     &         (nnod_rtp, nidx_rtp, irt_rtp_smp_stack,                &
!!     &         ncomp, n_WR, irev_sr_rtp, WR, X_rtp)
!! ------------------------------------------------------------------
!!
!! wrapper subroutine for backward Fourier transform by ISPACK
!!
!! x_{k} = a_{0} + (-1)^{j} a_{Nfft/2} + sum_{k=1}^{Nfft/2-1}
!! (a_{k} \cos(2\pijk/Nfft) + b_{k} \sin(2\pijk/Nfft))
!!
!! ------------------------------------------------------------------
!!
!! i = 1:     a_{0}
!! i = 2:     a_{Nfft/2}
!! i = 3:     a_{1}
!! i = 4:     b_{1}
!! ...
!! i = 2*k+1: a_{k}
!! i = 2*k+2: b_{k}
!! ...
!! i = Nfft-1:   a_{Nfft/2-1}
!! i = Nfft:     b_{Nfft/2-1}
!!
!! ------------------------------------------------------------------
!!@endverbatim
!!
!!@n @param Nsmp  Number of SMP processors
!!@n @param Nstacksmp(0:Nsmp)   End number for each SMP process
!!@n @param M           Number of components for Fourier transforms
!!@n @param Nfft        Data length for eadh FFT
!!@n @param X(M, Nfft)  Data for Fourier transform
!
      module m_sph_ISPACK_FFT
!
      use m_precision
      use m_constants
      use m_machine_parameter
      use t_sph_ISPACK_FFT
!
      implicit none
!
!>      Structure to use ISPACK
      type(work_for_ispack), save :: sph_ispack
!
      private :: sph_ispack
!
! ------------------------------------------------------------------
!
      contains
!
! ------------------------------------------------------------------
!
      subroutine init_sph_ISPACK(nidx_rtp, maxirt_rtp_smp, ncomp)
!
      integer(kind = kint), intent(in) :: maxirt_rtp_smp
      integer(kind = kint), intent(in) :: nidx_rtp(3)
      integer(kind = kint), intent(in) :: ncomp
!
!
      call init_sph_ISPACK_t(ncomp, nidx_rtp, maxirt_rtp_smp,           &
     &    sph_ispack)
!
      end subroutine init_sph_ISPACK
!
! ------------------------------------------------------------------
!
      subroutine finalize_sph_ISPACK
!
!
      call finalize_sph_ISPACK_t(sph_ispack)
!
      end subroutine finalize_sph_ISPACK
!
! ------------------------------------------------------------------
!
      subroutine verify_sph_ISPACK(nidx_rtp, maxirt_rtp_smp, ncomp)
!
      integer(kind = kint), intent(in) :: maxirt_rtp_smp
      integer(kind = kint), intent(in) :: nidx_rtp(3)
      integer(kind = kint), intent(in) :: ncomp
!
!
      call verify_sph_ISPACK_t(ncomp, nidx_rtp, maxirt_rtp_smp,         &
     &    sph_ispack)
!
      end subroutine verify_sph_ISPACK
!
! ------------------------------------------------------------------
! ------------------------------------------------------------------
!
      subroutine sph_FTTRUF_to_send                                     &
     &         (nnod_rtp, nidx_rtp, irt_rtp_smp_stack,                  &
     &          ncomp, n_WS, irev_sr_rtp, X_rtp, WS)
!
      integer(kind = kint), intent(in) :: nnod_rtp
      integer(kind = kint), intent(in) :: nidx_rtp(3)
      integer(kind = kint), intent(in) :: irt_rtp_smp_stack(0:np_smp)
!
      integer(kind = kint), intent(in) :: ncomp
      real(kind = kreal), intent(in)                                    &
     &     :: X_rtp(irt_rtp_smp_stack(np_smp),nidx_rtp(3),ncomp)
!
      integer(kind = kint), intent(in) :: n_WS
      integer(kind = kint), intent(in) :: irev_sr_rtp(nnod_rtp)
      real (kind=kreal), intent(inout):: WS(n_WS)
!
!
      call sph_FTTRUF_to_send_t(ncomp, nnod_rtp, nidx_rtp,              &
     &    irt_rtp_smp_stack, n_WS, irev_sr_rtp, X_rtp, WS, sph_ispack)
!
      end subroutine sph_FTTRUF_to_send
!
! ------------------------------------------------------------------
!
      subroutine sph_FTTRUB_from_recv                                   &
     &         (nnod_rtp, nidx_rtp, irt_rtp_smp_stack,                  &
     &          ncomp, n_WR, irev_sr_rtp, WR, X_rtp)
!
      integer(kind = kint), intent(in) :: nnod_rtp
      integer(kind = kint), intent(in) :: nidx_rtp(3)
      integer(kind = kint), intent(in) :: irt_rtp_smp_stack(0:np_smp)
!
      integer(kind = kint), intent(in) :: ncomp
      integer(kind = kint), intent(in) :: n_WR
      integer(kind = kint), intent(in) :: irev_sr_rtp(nnod_rtp)
      real (kind=kreal), intent(inout):: WR(n_WR)
!
      real(kind = kreal), intent(inout)                                 &
     &     :: X_rtp(irt_rtp_smp_stack(np_smp),nidx_rtp(3),ncomp)
!
!
      call sph_FTTRUB_from_recv_t(ncomp, nnod_rtp, nidx_rtp,            &
     &    irt_rtp_smp_stack, n_WR, irev_sr_rtp, WR, X_rtp, sph_ispack)
!
      end subroutine sph_FTTRUB_from_recv
!
! ------------------------------------------------------------------
! ------------------------------------------------------------------
!
      end module m_sph_ISPACK_FFT
