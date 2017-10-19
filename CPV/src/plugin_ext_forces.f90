!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_ext_forces(stress)
  !----------------------------------------------------------------------------
  !
  !
  USE mp_global,        ONLY : intra_image_comm
  USE mp,               ONLY : mp_bcast
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : outdir
  !
  USE plugin_flags
  !
  !HK/BS
  USE cell_base,        ONLY : h
  USE ions_base,        ONLY : nat, amass
  USE ions_positions,   ONLY : fion,tau0 !HK/BS
  USE cp_main_variables,ONLY : nfi
  USE energies,         ONLY : etot 
  !
  IMPLICIT NONE
  !
  REAL(DP), intent(inout) :: stress(3,3) !HK/BS: stress to be modified by plumed (we do it here because the stress is local variable in cpr.f90
  INTEGER:: i,j
  REAL(DP) :: at_plumed(3,3)
  REAL(DP) :: virial(3,3)
  REAL(DP) :: volume
  REAL(DP), ALLOCATABLE :: tau_plumed(:,:)
  !
#ifdef PLUMED2
  IF(use_plumed) then
    IF(ionode)THEN
      !at_plumed=alat*at;  ! the cell, rescaled properly
      at_plumed=transpose(h);  ! the cell, rescaled properly to bohr !HK/BS: TODO check if transpose is needed
      allocate(tau_plumed(3,nat))
      !tau_plumed=alat*tau
      tau_plumed=tau0 !HK/BS: TODO: check if the sorted version or the unsorted version should be used
      volume=+at_plumed(1,1)*at_plumed(2,2)*at_plumed(3,3) &
             +at_plumed(1,2)*at_plumed(2,3)*at_plumed(3,1) &
             +at_plumed(1,3)*at_plumed(2,1)*at_plumed(3,2) &
             -at_plumed(1,1)*at_plumed(3,2)*at_plumed(2,3) &
             -at_plumed(1,2)*at_plumed(3,3)*at_plumed(2,1) &
             -at_plumed(1,3)*at_plumed(3,1)*at_plumed(2,2) 
      !virial=-sigma*volume
      virial=-stress*volume

      CALL plumed_f_gcmd("setStep"//char(0),nfi)
      CALL plumed_f_gcmd("setMasses"//char(0),amass)
      CALL plumed_f_gcmd("setForces"//char(0),fion)
      CALL plumed_f_gcmd("setPositions"//char(0),tau_plumed)
      CALL plumed_f_gcmd("setBox"//char(0),at_plumed)
      CALL plumed_f_gcmd("setVirial"//char(0),virial)
      CALL plumed_f_gcmd("setEnergy"//char(0),etot)
      CALL plumed_f_gcmd("calc"//char(0),0)

      !sigma=-virial/volume
      stress=-virial/volume

      deallocate(tau_plumed)
    ENDIF
    CALL mp_bcast(fion, ionode_id, intra_image_comm)
    CALL mp_bcast(stress, ionode_id, intra_image_comm)
  ENDIF
#endif
  !
  !
END SUBROUTINE plugin_ext_forces
