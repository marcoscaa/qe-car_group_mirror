MODULE cp_bo_emin_module
!
IMPLICIT NONE
!
SAVE
!
! PUBLIC variables 
!
INTEGER, PUBLIC :: nsteps_emin 
!
! PUBLIC subroutines
!
PUBLIC :: cp_bo_emin, cp_bo_print_info
!
CONTAINS
!
!----------------------------------------------------------------------------
SUBROUTINE cp_bo_emin(tfirst,tlast,enb,enbi,ccc,stress,bigr,iter)
  !----------------------------------------------------------------------------
  !
  USE autopilot,                ONLY : event_index_emin, lpilot_electron_damping_emin
  USE kinds,                    ONLY : DP
  USE gvecw,                    ONLY : ngw,ggp,ecutwfc,g2kin_init
  USE electrons_base,           ONLY : nbspx,nbsp_bgrp,nbsp,nbspx_bgrp,f
  USE wavefunctions_module,     ONLY : c0_bgrp, cm_bgrp, phi_bgrp
  USE input_parameters,         ONLY : dt_emin,emass_emin,emass_cutoff_emin,electron_damping_emin,dforce_emin
  USE input_parameters,         ONLY : electron_maxstep,conv_thr,ref_cell
  USE cp_autopilot,             ONLY : pilot_emin
  USE cp_main_variables,        ONLY : nfi,ema0bg,eigr,lambda,descla,bephi,becp_bgrp,bec_bgrp,sfac,dbec
  USE cp_main_variables,        ONLY : eigrb,taub,irb,acc,lambda,lambdam,rhor,iprint_stdout 
  USE cell_base,                ONLY : h,init_tpiba2,bg,r_to_s,ainv,alat,at,ibrav,velh,celldm,omega,ref_tpiba2,tpiba2
  USE efcalc,                   ONLY : wf_efield
  USE control_flags,            ONLY : conv_elec,ldriver,tpre,lwf,tfor,tprnfor,nbeg,thdyn,ts_vdw,iprint,nomore
  USE ions_positions,           ONLY : fion,taus,tau0,tausm,vels,velsm
  USE ions_base,                ONLY : nsp,nat,na
  USE energies,                 ONLY : enthal,etot,ekincm,exx,eht,epseu,exc
  USE uspp_param,               ONLY : nvb
  USE uspp,                     ONLY : vkb
  USE io_global,                ONLY : stdout
  USE cp_interfaces,            ONLY : move_electrons,ortho,elec_fakekine,calbec_bgrp,caldbec_bgrp
  USE cp_interfaces,            ONLY : strucf,phfacs,nlfh,prefor,printout_new
  USE cp_electronic_mass,       ONLY : emass_precond,emass,emass_cutoff
  USE orthogonalize_base,       ONLY : updatc
  USE printout_base,            ONLY : printout_base_open,printout_base_close,printout_pos,printout_cell,printout_stress
  USE mp_global,                ONLY : me_image,inter_bgrp_comm
  USE gvect,                    ONLY : mill,eigts1,eigts2,eigts3,ecutrho,gg
  USE gvecs,                    ONLY : ngms
  USE fft_base,                 ONLY : dfftp
  USE time_step,                ONLY : dt2,delt,tps
  USE wave_base,                ONLY : frice
  USE wannier_subroutines,      ONLY : wf_closing_options
  USE electrons_nose,           ONLY : xnhe0,xnhem,vnhe
  USE ions_nose,                ONLY : nhpcl,nhpdim,vnhp,xnhp0,xnhpm
  USE cell_nose,                ONLY : xnhh0,xnhhm,vnhh
  USE ensemble_dft,             ONLY : z0t
  USE funct,                    ONLY : dft_is_hybrid, start_exx, exx_is_active
  USE exx_module,               ONLY : exxalfa
  USE constants,                ONLY : au_gpa
  USE tsvdw_module,             ONLY : EtsvdW,vdw_lscreen,vdw_lstep0
  !
  IMPLICIT NONE
  !
  ! I/O variables...
  !
  LOGICAL :: tfirst,tlast
  INTEGER :: iter
  REAL(DP) :: enb,enbi,ccc,stress(3,3),bigr!,stress_gpa(3,3)
  !
  ! Local variables...
  !
  INTEGER  :: nfi_emin,ia
  REAL(DP) :: dt2bye_emin,fccc_emin,etotm,ekinc_emin,dt2bye,fccc,electron_damping,emass_cutoff_emin_new
  REAL(DP) :: fionm(3,nat),dfion(3,nat),dfion_f(nat),max_dfion,mad_dfion,delta_etot,epot,temps(nat)
  COMPLEX(DP), ALLOCATABLE :: psim_bgrp(:,:)
  LOGICAL  :: damp_adapt, tprint, tfile, tstdout
  !
  CALL start_clock( 'cp_bo' )
  !
  ! Code inconsistencies (warnings)...
  !
  IF ( wf_efield ) CALL errore( 'cp_bo_emin', 'CP-BO and wf_efield not compatible...', 1 )
  !
  ! CP-BO electron minimization steps...
  !
  conv_elec=.FALSE.
  !
  !in tsvdw calculations atomic overlaps should be computed only at the first step 
  !
  IF (ts_vdw) THEN
    !
    vdw_lscreen=.TRUE.
    vdw_lstep0=.FALSE.
    !
  END IF
  !
  IF ( ldriver ) THEN
    !
    ! special treatment for the first step during restart ... 
    ! first bead comes in during restart is the bead from the previous step ...
    ! so we use cm, the optimized and orthogonalized wavefuntion from the previous step ... 
    !
    IF(tfirst) THEN
      IF(nbeg.GE.0) c0_bgrp=cm_bgrp
    END IF
    !
    CALL r_to_s(tau0,taus,na,nsp,ainv)
    !
    CALL phfacs(eigts1,eigts2,eigts3,eigr,mill,taus,dfftp%nr1,dfftp%nr2,dfftp%nr3,nat)
    !
    CALL prefor( eigr, vkb ) ! BS: may be needed in future
    !
    CALL strucf(sfac,eigts1,eigts2,eigts3,mill,ngms)
    !
    CALL calbec_bgrp(nvb+1,nsp,eigr,c0_bgrp,bec_bgrp)
    !
    ! Orthogonalization of extrapolated wavefunction c0 ...
    !
    CALL ortho(eigr,c0_bgrp,phi_bgrp,lambda,descla,bigr,iter,ccc,bephi,becp_bgrp)
    CALL updatc(ccc,lambda,phi_bgrp,bephi,becp_bgrp,bec_bgrp,c0_bgrp,descla)
    CALL calbec_bgrp(nvb+1,nsp,eigr,c0_bgrp,bec_bgrp)
    IF(tpre) CALL caldbec_bgrp(eigr,c0_bgrp,dbec,descla)
    !
    ! Compute Wannier functions for c0 ... 
    !
    IF(lwf) CALL wf_closing_options( nfi, c0_bgrp, cm_bgrp, bec_bgrp, eigr, eigrb,&
                 taub, irb, ibrav, bg(:,1), bg(:,2), bg(:,3), taus, tausm, vels,&
                 velsm, acc, lambda, lambdam, descla, xnhe0, xnhem,  &
                 vnhe, xnhp0, xnhpm, vnhp, nhpcl, nhpdim,    &
                 ekincm, xnhh0, xnhhm, vnhh, velh, ecutrho,  &
                 ecutwfc,delt,celldm, fion, tps, z0t, f, rhor )
    !  
  END IF
  !
  ! Set variables for CP-BO electron minimization...
  !
  dt2bye_emin=(dt_emin*dt_emin)/emass_emin
  !
  ! only for 0th and 1st step in from_scratch calculations a damping 0.1 is set here ...
  !
  !IF(nfi.EQ.1.AND.nbeg.LT.0) electron_damping=0.1_DP
  IF(nfi.LE.1.AND.nbeg.LT.0) electron_damping = 0.1_DP
  !
  etotm=0.0_DP
  nsteps_emin=0
  nfi_emin=-1 
  IF(nfi.EQ.1) nfi_emin=1 
  !
  ! Temporary storage of wavefunction cm ...
  !
  ALLOCATE(psim_bgrp(ngw,nbspx)); psim_bgrp=0.0_DP
  ! 
  psim_bgrp=cm_bgrp 
  !
  ! First step during CP-BO electron minimization = steepest descent
  ! (using extrapolated wf guess) ...
  !
  cm_bgrp=c0_bgrp
  !
  ! Set mass-preconditioning for Fourier acceleration during CP-BO electron
  ! minimization...
  !
  !compute ema0bg according to emass_cutoff_emin
  IF ( ref_cell ) THEN
    CALL g2kin_init( gg, ref_tpiba2 )
    CALL emass_precond( ema0bg, ggp, ngw, ref_tpiba2, emass_cutoff_emin ) 
    CALL g2kin_init( gg, tpiba2 )
  ELSE
    CALL g2kin_init( gg, init_tpiba2 )
    CALL emass_precond( ema0bg, ggp, ngw, init_tpiba2, emass_cutoff_emin ) 
    CALL g2kin_init( gg, tpiba2 )
  END IF
  !
  ! manage friction during minimization ...
  !
 !electron_damping=electron_damping_emin
 !emass_cutoff_emin_new=emass_cutoff_emin
  damp_adapt=.true.
  !
  ! ion force ...
  !
  fionm=0.0_DP; dfion=0.0_DP; dfion_f=0.0_DP
  !
  ! Print header for energy, forces
  !
  IF(dft_is_hybrid().AND.exx_is_active()) THEN
    !
    IF (ts_vdw) THEN
      WRITE(stdout,19472)
    ELSE
      WRITE(stdout,19470)
    END IF
    !
  ELSE
    IF (ts_vdw) THEN
      WRITE(stdout,19471)
    ELSE
      WRITE(stdout,1947)
    END IF
    !
  END IF
  !
  !
  DO WHILE (.NOT.conv_elec)
    !
    ! Increment step counter... 
    !
    nsteps_emin=nsteps_emin+1
    !
    ! Autopilot (Dynamic Rules) Implimentation
    ! this part will control only:  *_emin    
    !
    CALL pilot_emin(nsteps_emin)
    !
    ! manage friction ...
    !
    IF (lpilot_electron_damping_emin) THEN
      electron_damping=electron_damping_emin
    ELSE
      IF(nbeg.GE.0.OR.nfi.GT.1) THEN
        IF(nsteps_emin.EQ.1) electron_damping=0.0_DP 
      END IF
    END IF
    !
    fccc_emin=1.0_DP/(1.0_DP+electron_damping)
    !
    !CALL emass_precond(ema0bg,ggp,ngw,init_tpiba2,emass_cutoff_emin_new) 
    !
    ! Perform single step of second-order damped dynamics...
    !
    CALL move_electrons(nfi_emin,tfirst,tlast,bg(:,1),bg(:,2),bg(:,3),fion,c0_bgrp, & 
                        cm_bgrp,phi_bgrp,enthal,enb,enbi,fccc_emin,ccc,dt2bye_emin,stress,.false.)
    !
    ! Orthogonalization steps...
    !
    CALL ortho(eigr,cm_bgrp,phi_bgrp,lambda,descla,bigr,iter,ccc,bephi,becp_bgrp)
    CALL updatc(ccc,lambda,phi_bgrp,bephi,becp_bgrp,bec_bgrp,cm_bgrp,descla)
    CALL calbec_bgrp(nvb+1,nsp,eigr,cm_bgrp,bec_bgrp)
    !
    ! Compute fictitious kinetic energy...
    !
    CALL elec_fakekine(ekinc_emin,ema0bg,emass_emin,c0_bgrp,cm_bgrp,ngw,nbsp_bgrp,1,dt_emin)
    !
    ! Swap c0 and cm... after swap c0 is the most updated wavefunction ...
    !
    CALL dswap( 2*SIZE( c0_bgrp ), c0_bgrp, 1, cm_bgrp, 1 )
    !
    ! compute change in ionic forces ...
    !
    dfion=DABS(fion-fionm)
    ! 
    DO ia=1,nat
      !
      dfion_f(ia)=DSQRT((dfion(1,ia)*dfion(1,ia)+dfion(2,ia)*dfion(2,ia)+dfion(3,ia)*dfion(3,ia)) &
                       /(fion(1,ia)*fion(1,ia)+fion(2,ia)*fion(2,ia)+fion(3,ia)*fion(3,ia))) 
      !
    END DO 
    !
    mad_dfion=SUM(dfion_f(:))/DBLE(nat)
    max_dfion=MAXVAL(dfion_f(:))
    !
    fionm=fion
    !
    delta_etot=etot-etotm
    !
    ! Print energy, forces at each emin step 
    !
    IF(dft_is_hybrid().AND.exx_is_active()) THEN
      !
      IF (ts_vdw) THEN
        WRITE(stdout,19482)nfi,nsteps_emin,etot,-delta_etot/etot,mad_dfion,ekinc_emin,electron_damping,(-exx*exxalfa),EtsvdW
      ELSE
        WRITE(stdout,19480)nfi,nsteps_emin,etot,-delta_etot/etot,mad_dfion,ekinc_emin,electron_damping,(-exx*exxalfa)
      END IF
      !
    ELSE
      IF (ts_vdw) THEN
        WRITE(stdout,19481)nfi,nsteps_emin,etot,-delta_etot/etot,mad_dfion,ekinc_emin,electron_damping,EtsvdW
      ELSE
        WRITE(stdout,1948)nfi,nsteps_emin,etot,-delta_etot/etot,mad_dfion,ekinc_emin,electron_damping
      END IF
      !
    END IF
    !
    ! Convergence loop conditionals...
    !
    IF (DABS(delta_etot).LT.(conv_thr*1.E+1_DP).AND.(mad_dfion.LT.dforce_emin)) conv_elec=.TRUE.
    !
    ! Compute Wannier functions for c0 ... 
    !
    IF(lwf) CALL wf_closing_options( nfi, c0_bgrp, cm_bgrp, bec_bgrp, eigr, eigrb,&
                 taub, irb, ibrav, bg(:,1), bg(:,2), bg(:,3), taus, tausm, vels,&
                 velsm, acc, lambda, lambdam, descla, xnhe0, xnhem,  &
                 vnhe, xnhp0, xnhpm, vnhp, nhpcl, nhpdim,    &
                 ekincm, xnhh0, xnhhm, vnhh, velh, ecutrho,  &
                 ecutwfc,delt,celldm, fion, tps, z0t, f, rhor )
    !  
    !exx_wf related
    !The following criteria is used to turn on exact exchange calculation when
    !GGA energy is converged up to 100 times of the input etot convergence
    !thereshold  
    !
    IF( .NOT.exx_is_active().AND.dft_is_hybrid() ) THEN
      !
      IF(DABS(delta_etot).LT.conv_thr*1.E+2_DP.AND.delta_etot.LT.0.0_DP) THEN
        !
        WRITE(stdout,'(/,3X,"Exact Exchange is turned on ...")')
        ! 
        CALL start_exx()
        !
      END IF
      !
    END IF
    !
    ! manage friction ...
    !
    IF (lpilot_electron_damping_emin) THEN
      electron_damping=electron_damping_emin
    ELSE
      IF(nbeg.GE.0.OR.nfi.GT.1) THEN
        IF(damp_adapt) THEN
          IF(DABS(delta_etot).LT.(conv_thr*1.E+4_DP).AND.(etot-etotm).LT.0.0_DP) THEN
            !   
            IF(ldriver) THEN 
              electron_damping=0.25_DP !PI
            ELSE  
              electron_damping=0.30_DP !CL
            END IF  
            !   
            damp_adapt=.FALSE.
            !   
          ELSE
            !   
            electron_damping=electron_damping_emin
            !   
          END IF
        END IF
      END IF
    END IF
    !
    etotm=etot
    !
    !in tsvdw calculations atomic overlap computation turned off
    !
    IF (ts_vdw) THEN
      !
      vdw_lscreen=.FALSE.
      IF (conv_elec) vdw_lstep0 = .TRUE.
      !
    END IF
    !
    IF(tpre.AND.conv_elec) CALL caldbec_bgrp(eigr,c0_bgrp,dbec,descla)
    !
    ! Maximum step conditional...
    !
    IF (nsteps_emin.GT.electron_maxstep) CALL errore( 'cp_bo_emin', 'CP-BO ran out of minimization steps...', 1 )
    !
  END DO
  !
  ! reset autopilot event_index_emin !HK
  event_index_emin = 1
  !
  ! Replace wavefunctions... cm = minimized c(t-dt) and c0 = minimized c(t)
  ! For nfi=1, cm and c0 are the last 2 steps during cp-bo minimization (to
  ! avoid ortho issues during extrapolation) ...
  !
  IF(nfi.NE.1.OR.nbeg.GE.0) cm_bgrp=psim_bgrp 
  !
  ! Clean-up temporary arrays...
  !
  IF (ALLOCATED(psim_bgrp))   DEALLOCATE(psim_bgrp)
  !
  !compute ema0bg according to emass_cutoff
  IF ( ref_cell ) THEN
    CALL g2kin_init( gg, ref_tpiba2 )
    CALL emass_precond( ema0bg, ggp, ngw, ref_tpiba2, emass_cutoff ) 
    CALL g2kin_init( gg, tpiba2 )
  ELSE
    CALL g2kin_init( gg, init_tpiba2 )
    CALL emass_precond( ema0bg, ggp, ngw, init_tpiba2, emass_cutoff ) 
    CALL g2kin_init( gg, tpiba2 )
  END IF
  !
  IF ( ldriver ) THEN
    !
    ! computation of energy, ion forces, stress, and wavefunction extrapolation ...
    ! 
    dt2bye=dt2/emass
    fccc=1.0_DP/(1.0_DP+frice)
    !
    CALL move_electrons(nfi_emin,tfirst,tlast,bg(:,1),bg(:,2),bg(:,3),fion,c0_bgrp, &
                        cm_bgrp,phi_bgrp,enthal,enb,enbi,fccc,ccc,dt2bye,stress,.false.)
    !
    delta_etot=etot-etotm
    !
    IF(tpre) THEN
      CALL caldbec_bgrp(eigr,c0_bgrp,dbec,descla)
      CALL nlfh(stress,bec_bgrp,dbec,lambda,descla)
    END IF
    !
    ! Swap c0 and cm... after swap 
    ! c0 = non-orthogonal extrapolated wavefunction at time t+dt and 
    ! cm =     orthogonal minimized wavefunction at time t ...
    !
    CALL dswap( 2*SIZE( c0_bgrp ), c0_bgrp, 1, cm_bgrp, 1 )
    !
    ! compute change in ionic forces ...
    !
    dfion=DABS(fion-fionm)
    ! 
    DO ia=1,nat
      !
      dfion_f(ia)=DSQRT((dfion(1,ia)*dfion(1,ia)+dfion(2,ia)*dfion(2,ia)+dfion(3,ia)*dfion(3,ia)) &
                       /(fion(1,ia)*fion(1,ia)+fion(2,ia)*fion(2,ia)+fion(3,ia)*fion(3,ia))) 
      !
    END DO 
    !
    mad_dfion=SUM(dfion_f(:))/DBLE(nat)
    max_dfion=MAXVAL(dfion_f(:))
    !
    ! Print energy, forces at each emin step 
    !
    IF(dft_is_hybrid().AND.exx_is_active()) THEN
      !
      IF (ts_vdw) THEN
        WRITE(stdout,19492)nfi,nsteps_emin,etot,-delta_etot/etot,mad_dfion,ekinc_emin,electron_damping,(-exx*exxalfa),EtsvdW
      ELSE
        WRITE(stdout,19490)nfi,nsteps_emin,etot,-delta_etot/etot,mad_dfion,ekinc_emin,electron_damping,(-exx*exxalfa)
      END IF
      !
    ELSE
      IF (ts_vdw) THEN
        WRITE(stdout,19491)nfi,nsteps_emin,etot,-delta_etot/etot,mad_dfion,ekinc_emin,electron_damping,EtsvdW
      ELSE
        WRITE(stdout,1949)nfi,nsteps_emin,etot,-delta_etot/etot,mad_dfion,ekinc_emin,electron_damping
      END IF
      !
    END IF
    !
    !IF (tpre) THEN
    !  stress_gpa = stress * au_gpa
    !  CALL printout_stress( stdout, stress_gpa )
    !END IF
    !
    !print files
    !
    tlast   = ( nfi == nomore ) .OR. tlast
    tprint  = ( MOD( nfi, iprint ) == 0 ) .OR. tlast 
    tfile   = ( MOD( nfi, iprint ) == 0 )
    tstdout = ( MOD( nfi, iprint_stdout ) == 0 ) .OR. tlast
    !
    epot = eht + epseu + exc
    temps=0.0_DP
    !
    CALL printout_new( nfi, tfirst, tfile, tprint, tps, h, stress, &
                       tau0, vels, fion, 0.0_DP, 0.0_DP, 0.0_DP, temps, etot, &
                       enthal, 0.0_DP, 0.0_DP, vnhh, xnhh0, vnhp, xnhp0, 0.0_DP, &
                       0.0_DP, epot, tprnfor, tpre, tstdout )
    !
  END IF
  !
  CALL stop_clock( 'cp_bo' )
  !
1947   FORMAT(/,20X,'nfi',4X,'nfe',5X,'E',17X,'dE/E',7X,'mad_dfion',5X,'ekinc',5X,'damping')                    ! GGA
19470  FORMAT(/,20X,'nfi',4X,'nfe',5X,'E',17X,'dE/E',7X,'mad_dfion',5X,'ekinc',5X,'damping',6X,'exx')           ! Hybrid 
19471  FORMAT(/,20X,'nfi',4X,'nfe',5X,'E',17X,'dE/E',7X,'mad_dfion',5X,'ekinc',5X,'damping',6X,'evdw')          ! GGA+vdW
19472  FORMAT(/,20X,'nfi',4X,'nfe',5X,'E',17X,'dE/E',7X,'mad_dfion',5X,'ekinc',5X,'damping',6X,'exx',11X,'evdw')! Hybrid+vdW
  !
1948   FORMAT(3X,"CP-BO-STEP:",1X,2I7,F18.10,3ES13.3,F6.2,1X)        ! GGA 
19480  FORMAT(3X,"CP-BO-STEP:",1X,2I7,F18.10,3ES13.3,F6.2,3X,F14.6)  ! Hybrid 
19481  FORMAT(3X,"CP-BO-STEP:",1X,2I7,F18.10,3ES13.3,F6.2,3X,F14.6)  ! GGA+vdW 
19482  FORMAT(3X,"CP-BO-STEP:",1X,2I7,F18.10,3ES13.3,F6.2,3X,2F14.6) ! Hybrid+vdW                        
  !
1949   FORMAT(2X,"!CP-BO-STEP:",1X,2I7,F18.10,3ES13.3,F6.2,1X)        ! GGA 
19490  FORMAT(2X,"!CP-BO-STEP:",1X,2I7,F18.10,3ES13.3,F6.2,3X,F14.6)  ! Hybrid 
19491  FORMAT(2X,"!CP-BO-STEP:",1X,2I7,F18.10,3ES13.3,F6.2,3X,F14.6)  ! GGA+vdW 
19492  FORMAT(2X,"!CP-BO-STEP:",1X,2I7,F18.10,3ES13.3,F6.2,3X,2F14.6) ! Hybrid+vdW                        
  !
  RETURN
  !
END SUBROUTINE cp_bo_emin
!
SUBROUTINE cp_bo_print_info( )

  USE input_parameters,         ONLY : dt_emin, emass_emin, emass_cutoff_emin, electron_damping_emin, dforce_emin, conv_thr
  USE io_global,                ONLY : stdout

  IMPLICIT NONE

  WRITE( stdout, 545 )
  WRITE( stdout, 601 ) emass_emin, emass_cutoff_emin
  WRITE( stdout, 602 ) electron_damping_emin, dt_emin
  WRITE( stdout, 603 ) conv_thr, dforce_emin

  545     FORMAT(//,3X,'CP-BO Electron Minimizatoin Parameters (from STDIN)',/ &
    ,3X,'-------------------------------------')
  601  format( 3X, 'emass for CP-BO steps (emass_emin)                        = ',F10.3, &
    &       //,3X, 'emass_cutoff for CP-BO steps (emass_cutoff_emin)          = ',F10.3 )
  602  format( 3X, 'electron_damping for CP-BO steps (electron_damping_emin)  = ',ES13.3, &
    &       //,3X, 'dt for CP-BO steps (dt_emin)                              = ',F10.3 )
  603  format( 3X, 'Total DFT energy convergence threshold (conv_thr)         = ',ES13.3, &
    &       //,3X, 'force convergence threshold (dforce_emin)                 = ',ES13.3 )

  RETURN
END SUBROUTINE cp_bo_print_info

END MODULE cp_bo_emin_module
