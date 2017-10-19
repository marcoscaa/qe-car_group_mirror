!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stop_run()
  !----------------------------------------------------------------------------
  !
  ! ... Close all files and synchronize processes before stopping.
  !
  USE environment,        ONLY : environment_end
  USE control_flags,      ONLY : lconstrain,ldriver
  USE constraints_module, ONLY : deallocate_constraint
  USE mp_global,          ONLY : mp_global_end
  USE io_global,          ONLY : stdout
  ! HK/MCA : we don't know what is going on the ISO_C_BINGING ...
  !! HK/MCA : use the more up-to-date ipi socket machinary (including ISO_C_BINGING) in the fsockets.f90
  !USE F90SOCKETS, ONLY : error_socket
  !
  IMPLICIT NONE
  !
  !
  CALL environment_end( 'CP' )
  !
  !call error from socket to exit socket in driver mode
  !IF(ldriver) CALL error("@ CP DRIVER MODE: exit socket")
  if (ldriver) then
    ! HK/MCA : remove the -1 error return from error and print directly from here
    write(stdout,*) "@ CP DRIVER MODE: exit socket"
  end if ! ldriver
  !
  CALL deallocate_modules_var()
  !
  IF ( lconstrain ) CALL deallocate_constraint()
  !
  CALL mp_global_end()
  !
END SUBROUTINE stop_run

SUBROUTINE do_stop( flag )
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: flag
  !
  IF ( flag ) THEN
     STOP
  ELSE
     STOP 1
  END IF
  !
END SUBROUTINE do_stop
