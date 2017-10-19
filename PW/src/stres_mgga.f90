!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stres_mgga( rho, rhog, rho_core, rhog_core, kinedens, & 
                       nr1, nr2, nr3, nrxx, nl, ngm, g, omega, sigmaxc )
  !----------------------------------------------------------------------------
  !
  ! Analytic stress tensor contribution from metagga 
  !
  USE kinds,                  ONLY : DP
  USE control_flags,          ONLY : gamma_only
  USE noncollin_module,       ONLY : noncolin
  USE cell_base,              ONLY : tpiba
  USE scf,                    ONLY : v
  USE wavefunctions_module,   ONLY : evc, psic
  USE funct,                  ONLY : dft_is_meta, tau_xc
  USE klist,                  ONLY : nks, xk, ngk
  USE buffers,                ONLY : get_buffer
  USE io_files,               ONLY : iunwfc, nwordwfc, iunigk
  USE wvfct,                  ONLY : nbnd, npwx, npw, wg, igk 
  USE lsda_mod,               ONLY : lsda, nspin, current_spin, isk
  USE fft_interfaces,         ONLY : fwfft, invfft
  USE mp,                     ONLY : mp_sum
  USE mp_pools,               ONLY : inter_pool_comm
  USE mp_bands,               ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  INTEGER                   :: nr1, nr2, nr3, nrxx, ngm, nl(ngm)
  REAL(DP)                  :: kinedens(nrxx,nspin), sigmaxc(3,3)
  REAL(DP)                  :: rho (nrxx, nspin), rho_core (nrxx)
  REAL(DP)                  :: omega, g(3, ngm) 
  COMPLEX(DP)               :: rhog(ngm, nspin), rhog_core(ngm)
  !
  ! Internal variables
  !
  INTEGER                   :: ix, iy, ir, ipol, iss, incr, ibnd, ik
  INTEGER                   :: ipol2xy(3,3), nspin0 
  REAL(DP), PARAMETER       :: epsr = 1.0d-6, epsg = 1.0d-10, e2 = 2.d0
  REAL(DP) , ALLOCATABLE    :: grho (:,:,:)
  COMPLEX(DP), ALLOCATABLE  :: gradwfc (:,:), crosstaus(:,:,:)
  !
  REAL(DP)                  :: w1, w2, delta 
  REAL(DP)                  :: sx, sc, v1x, v2x,v3x,v1c,v2c,v3c 
  REAL(DP)                  :: grh2, grho2 (2), fac 
  REAL(DP)                  :: sigma_mgga(3,3), sigma_gradcorr (3, 3) 
  !
  !
  if ( .not. dft_is_meta() ) return
  !
  current_spin=1
  !
  ! Stop if something is not yet implemented
  !
  if (noncolin) call errore('stres_mgga', &
                    'noncollinear stress + meta-GGA not implemented',1)
  if (nspin>1) call errore('stres_mgga', &
                    'spin polarized mgga stress not implemented',1)
  !
  ! Initialization of a set of variables
  !
  allocate (gradwfc( nrxx, 3))    
  allocate (crosstaus( nrxx,6,nspin))    
  allocate (grho( 3, nrxx, nspin))    
  !
  nspin0=nspin
  if (nspin==4) nspin0=1
  fac = 1.d0 / DBLE (nspin0)
  !
  sigma_mgga(:,:)=0_DP
  sigma_gradcorr(:,:)=0_DP
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! First, contribution from the gradient dependent part of 
  ! the exchange correlation energy
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !    calculate the gradient of rho+rhocore in real space
  !
  DO iss = 1, nspin0
     !
     rho(:,iss)  = fac * rho_core(:)  + rho(:,iss)
     rhog(:,iss) = fac * rhog_core(:) + rhog(:,iss)
     !
     CALL gradrho( nrxx, rhog(1,iss), ngm, g, nl, grho(1,1,iss) )
     !
  END DO
  !
  ! Recompute the exchange-correlation energy 
  !
  IF (nspin.eq.1) THEN
     !
     DO ir = 1, nrxx
        !
        grho2 (1) = grho(1,ir,1)**2 + grho(2,ir,1)**2 + grho(3,ir,1)**2
        !
        IF (abs (rho (ir, 1) ) .gt.epsr.and.grho2 (1) .gt.epsg) THEN
           !
           kinedens(ir,1) = kinedens(ir,1) / e2
           call tau_xc (rho(ir,1), grho2(1),kinedens(ir,1), sx, sc, v1x, v2x,v3x,v1c,v2c,v3c)
           kinedens(ir,1) = kinedens(ir,1) * e2
           !
           ! Grad rho contributio to the stress tensor 
           !
           DO ix = 1, 3
              !
              DO iy = 1, ix
                 !
                 sigma_gradcorr (ix, iy) = sigma_gradcorr (ix, iy) + &
                                grho(ix,ir,1) * grho(iy,ir,1) * e2 * (v2x + v2c)
                 !
              ENDDO
              !
           ENDDO
           !
        ENDIF
        !
     ENDDO
     !
  ENDIF
  !
  DO ix = 1, 3
     !
     DO iy = 1, ix - 1
        !
        sigma_gradcorr (iy, ix) = sigma_gradcorr (ix, iy)
        !
     ENDDO
     !
  ENDDO
  !
  call mp_sum(  sigma_gradcorr, intra_bgrp_comm )
  !
  call dscal (9, 1.d0 / (nr1 * nr2 * nr3), sigma_gradcorr, 1)
  !
  call daxpy (9, 1.d0, sigma_gradcorr, 1, sigmaxc, 1)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 
  ! Kinetic energy density contribution to the stress tensor
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! For gamma_only efficiency
  !
  incr=1
  IF ( gamma_only ) incr=2 
  !
  crosstaus(:,:,:) = 0.d0
  gradwfc(:,:) = 0.d0
  !
  !
  if (nks.gt.1) rewind (iunigk)
  !
  ! Loop over the k points
  !
  k_loop: DO ik = 1, nks
  !
    !
    IF ( lsda ) current_spin = isk(ik)
    !
    npw = ngk(ik)
    !
    ! Read the wavefunctions 
    !
    IF ( nks > 1 ) THEN
       !
       read (iunigk) igk
       CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
       !
    END IF
    !
    do ibnd = 1, nbnd, incr 
       !
       ! MCA: weights for each k point and band
       !
       w1 = wg(ibnd,ik) / omega 
       !
       IF ( (ibnd < nbnd) .and. (gamma_only) ) THEN
          !
          ! ... two ffts at the same time
          !
          w2 = wg(ibnd+1,ik) / omega
          !
       ELSE
          !
          w2 = w1
          !
       END IF
       !
       ! Gradient of the wavefunction in real space
       ! 
       CALL wfc_gradient( ibnd, ik, npw, gradwfc )
       !
       ! Cross terms of kinetic energy density
       !
       ipol=1
       !
       do ix=1,3
          !
          do iy=1,ix
             !
             ipol2xy(ix,iy)=ipol
             ipol2xy(iy,ix)=ipol
             !
             !
             do ir=1,nrxx
                !
                crosstaus(ir,ipol,current_spin) = crosstaus(ir,ipol,current_spin) +&
                                        2.0_DP*w1*DBLE(gradwfc(ir,ix))*DBLE(gradwfc(ir,iy)) +&
                                        2.0_DP*w2*AIMAG(gradwfc(ir,ix))*AIMAG(gradwfc(ir,iy))
                !
             end do
             !
             !
             ipol=ipol+1
             !
          end do
          !
       end do
       !
    end do
    !
  END DO k_loop 
  !
  call mp_sum(  crosstaus, inter_pool_comm )
  !
  ! gradwfc not used anymore
  !
  deallocate (gradwfc)    
  !
  sigma_mgga(:,:) = 0.D0
  !
  ! metagga contribution to the stress tensor
  !
  do iss=1,nspin
     !
     do ix=1,3
        !
        do iy=1,3
           !
           delta=0.
           if (ix==iy) delta=1.
           !
           do ir=1,nrxx
              !
              sigma_mgga(ix,iy) = sigma_mgga(ix,iy) + v%kin_r(ir,iss) &
                                * ( kinedens(ir,iss) * delta &
                                + crosstaus(ir,ipol2xy(ix,iy),iss) )
              !
           end do
           !
           !
        end do
        !
     end do
     !
  end do
  !
  call mp_sum(  sigma_mgga, intra_bgrp_comm )
  !
  call dscal (9, 1.d0 / (nr1 * nr2 * nr3), sigma_mgga, 1)
  !
  call daxpy (9, 1.d0, sigma_mgga, 1, sigmaxc, 1)
  !
  DO iss = 1, nspin0
     !
     rho(:,iss)  = rho(:,iss)  - fac * rho_core(:)
     rhog(:,iss) = rhog(:,iss) - fac * rhog_core(:)
     !
  END DO
  !
  deallocate(grho)
  deallocate( crosstaus )
  return
  !  
END SUBROUTINE stres_mgga

SUBROUTINE wfc_gradient ( ibnd, ik, npw, gradpsi )
  !
  ! Returns the gradient of the wavefunction in real space
  ! Slightly adapted from sum_bands.f90 
  !
  USE kinds,                  ONLY : DP
  USE control_flags,          ONLY : gamma_only
  USE wavefunctions_module,   ONLY : psic, evc
  USE wvfct,                  ONLY : npwx, nbnd, igk
  USE cell_base,              ONLY : omega, tpiba
  USE klist,                  ONLY : xk
  USE gvect,                  ONLY : g
  USE gvecs,                  ONLY : nls, nlsm
  USE fft_base,               ONLY : dffts
  USE fft_interfaces,         ONLY : invfft
  !
  IMPLICIT NONE 
  !
  INTEGER                   :: ibnd, ik, npw 
  COMPLEX(DP)               :: gradpsi(dffts%nnr,3)  
  !
  ! Internal variables
  !
  REAL(DP)               :: kplusg(npwx)
  INTEGER                :: ipol
  !
  ! Compute the gradient of the wavefunction in reciprocal space
  !
  IF ( gamma_only ) THEN
     !
     DO ipol=1,3
        !
        psic(:) = ( 0.D0, 0.D0 )
        !
        kplusg (1:npw) = (xk(ipol,ik)+g(ipol,igk(1:npw))) * tpiba
        !
        IF ( ibnd < nbnd ) THEN
           !
           ! ... two ffts at the same time
           !
           psic(nls(1:npw)) = CMPLX(0d0, kplusg(1:npw),kind=DP)* &
                                   ( evc(1:npw,ibnd) + &
                                   ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
           !
           psic(nlsm(1:npw)) = CMPLX(0d0,-kplusg(1:npw),kind=DP) * &
                                    CONJG( evc(1:npw,ibnd) - &
                                    ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
           !
        ELSE
           !
           psic(nls(1:npw)) = CMPLX(0d0, kplusg(1:npw),kind=DP)* &
                                   evc(1:npw,ibnd)
           !
           psic(nlsm(1:npw)) = CMPLX(0d0,-kplusg(1:npw),kind=DP) * &
                                    CONJG( evc(1:npw,ibnd) )
           !
        END IF
        !
        ! Gradient of the wavefunction in real space
        !
        CALL invfft ('Wave', psic, dffts)
        !
        gradpsi(:,ipol) = psic
        !
     END DO
     !
  ELSE
     !
     DO ipol=1,3
         !
         psic(:) = ( 0.D0, 0.D0 )
         !
         kplusg (1:npw) = (xk(ipol,ik)+g(ipol,igk(1:npw))) * tpiba
         psic(nls(igk(1:npw))) = CMPLX(0d0, kplusg(1:npw),kind=DP)* &
                                 evc(1:npw,ibnd)
         !
         ! Gradient of the wavefunction in real space
         !
         CALL invfft ('Wave', psic, dffts)
         !
         gradpsi(:,ipol) = psic
         !
     END DO 
     !
  END IF
  !
END SUBROUTINE wfc_gradient
