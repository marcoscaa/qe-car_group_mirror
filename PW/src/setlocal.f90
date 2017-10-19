!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE setlocal
  !----------------------------------------------------------------------
  !
  !    This routine computes the local potential in real space vltot(ir)
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps8
  USE ions_base, ONLY : zv, ntyp => nsp
  USE cell_base, ONLY : omega
  USE extfield,  ONLY : tefield, dipfield, etotefield
  USE gvect,     ONLY : igtongl, gg
  USE scf,       ONLY : rho, v_of_0, vltot
  USE vlocal,    ONLY : strf, vloc
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces,ONLY : invfft
  USE gvect,     ONLY : nl, nlm, ngm
  USE control_flags, ONLY : gamma_only
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE martyna_tuckerman, ONLY : wg_corr_loc, do_comp_mt
  USE esm,       ONLY : esm_local, esm_bc, do_comp_esm
  USE qmmm,      ONLY : qmmm_add_mm_field
  !
  IMPLICIT NONE
  COMPLEX(DP), ALLOCATABLE :: aux (:), v_corr(:)
  ! auxiliary variable
  INTEGER :: nt, ng
  ! counter on atom types
  ! counter on g vectors
  !
  ALLOCATE (aux( dfftp%nnr))
  aux(:)=(0.d0,0.d0)
  !
  IF (do_comp_mt) THEN
      ALLOCATE(v_corr(ngm))
      CALL wg_corr_loc(omega,ntyp,ngm,zv,strf,v_corr)
      aux(nl(:)) = v_corr(:)
      DEALLOCATE(v_corr)
  END IF
  !
  DO nt = 1, ntyp
      DO ng = 1, ngm
          aux (nl(ng))=aux(nl(ng)) + vloc (igtongl (ng), nt) * strf (ng, nt)
      END DO
  END DO
  IF (gamma_only) THEN
      DO ng = 1, ngm
          aux (nlm(ng)) = CONJG(aux (nl(ng)))
      END DO
  END IF
  !
  IF ( do_comp_esm .AND. ( esm_bc .NE. 'pbc' ) ) THEN
     !
     ! ... Perform ESM correction to local potential
     !
      CALL esm_local ( aux )
  ENDIF
  !
  ! ... v_of_0 is (Vloc)(G=0)
  !
  v_of_0=0.0_DP
  IF (gg(1) < eps8) v_of_0 = DBLE ( aux (nl(1)) )
  !
  CALL mp_sum( v_of_0, intra_bgrp_comm )
  !
  ! ... aux = potential in G-space . FFT to real space
  !
  CALL invfft ('Dense', aux, dfftp)
  !
  vltot (:) =  DBLE (aux (:) )
  !
  ! ... If required add an electric field to the local potential 
  !
  IF ( tefield .AND. ( .NOT. dipfield ) )  &
      CALL add_efield(vltot,etotefield,rho%of_r,.TRUE.)
  !
  !  ... Add the electrostatic field generated by MM atoms
  !  in a QM/MM calculation to the local potential
  !
  CALL qmmm_add_mm_field()
  !
  DEALLOCATE(aux)
  !
  RETURN
END SUBROUTINE setlocal

