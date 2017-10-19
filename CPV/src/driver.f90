  SUBROUTINE driver()
    USE io_global,              ONLY : stdout,ionode,ionode_id
    USE io_files,               ONLY : srvaddress,iunupdate
    USE mp_global,              ONLY : mp_startup,intra_image_comm,inter_bgrp_comm,nbgrp
    USE mp,                     ONLY : mp_bcast
    USE control_flags,          ONLY : conv_elec,thdyn,isave,tpre,iverbosity,nomore
    USE cp_main_variables,      ONLY : nfi,acc,lambda,lambdam,descla,rhor 
    USE uspp,                   ONLY : okvan,nlcc_any
    USE cell_base,              ONLY : alat,at,omega,h,hold,velh,ainv
    USE ions_positions,         ONLY : fion,tau0,taus,tausm,vels,velsm 
    USE energies,               ONLY : etot,ekincm,eself
    USE wavefunctions_module,   ONLY : c0_bgrp,cm_bgrp
    USE electrons_nose,         ONLY : xnhe0,xnhem,vnhe 
    USE ions_nose,              ONLY : xnhp0,xnhpm,vnhp,nhpcl,nhpdim 
    USE cell_nose,              ONLY : xnhh0,xnhhm,xnhhp,vnhh
    USE time_step,              ONLY : tps, delt
    USE ensemble_dft,           ONLY : z0t
    USE electrons_base,         ONLY : f
    USE cp_interfaces,          ONLY : writefile,newinit
    USE cp_bo_emin_module,      ONLY : cp_bo_emin 
    USE constants,              ONLY : au_gpa, au_ps
    ! HK/MCA : use the more up-to-date ipi socket machinary (including ISO_C_BINGING) in the fsockets.f90
    USE F90SOCKETS, ONLY : open_socket, writebuffer, readbuffer

    IMPLICIT NONE
    
    INTEGER, PARAMETER :: MSGLEN=12
    LOGICAL :: isinit=.true., hasdata=.false., firststep=.true.
    CHARACTER*12 :: header
    CHARACTER*1024 :: parbuffer
    INTEGER socket, nat , flag ! HK/MCA
    INTEGER inet, port, ccmd, i, exst
    CHARACTER*1024 :: host  
    !REAL*8 :: cellh(3,3), cellih(3,3), vir(3,3), pot
    REAL*8 :: pot, vir(3,3)
    REAL*8 :: cellh(9), cellih(9), bufvir(9) ! HK/MCA todo : check the shape (transpose of the communication)
    REAL*8, ALLOCATABLE :: combuf(:)
    REAL*8 :: stress(3,3)=0.0d0 
    REAL*8 :: stress_gpa(3,3)=0.0d0 

    INTEGER  :: iter,ia 
    LOGICAL  :: tfirst,tlast
    REAL*8   :: enb,enbi,ccc,bigr
    INTEGER  :: j1,j2,jtmp

    !
    ! Inconsistencies with the present code ....
    !
    IF(okvan.OR.nlcc_any) CALL errore('driver','PI calculations are not working with okvan and nlcc_any...',1) 
    !
!   IF(thdyn) CALL errore('driver','PI calculations are not working with thdyn...',1) 
!   thdyn=.true.
!   tpre=.true.
    !
    port = 0
    flag=1
    tfirst = .TRUE. ! deal ..
    tlast  = .FALSE. ! deal .. 
    !
    inet=1
    host=srvaddress(1:INDEX(srvaddress,':',back=.true.)-1)//achar(0)
    
    IF (srvaddress(1:INDEX(srvaddress,':')-1).eq.('UNIX')) THEN
      inet=0
      host=srvaddress(INDEX(srvaddress,':',back=.true.)+1:)//achar(0)    
    else
      read(srvaddress(INDEX(srvaddress,':',back=.true.)+1:),*) port
    ENDIF
    
    IF (ionode) write(stdout,'(3X,"DRIVER MODE: Connecting to ",A)') trim(srvaddress)

    IF (ionode) call open_socket(socket, inet, port, host) ! HK/MCA : trivial to use fsockets.f90:open_socket (bug: fstr2cstr missing)
    !
    driver_loop: DO
      !
      CALL start_clock( 'driver_md' )
      !
      ! do communication on master node only...
      !
!       flag = 1
      if (ionode) call readbuffer(socket, header, MSGLEN, flag)
      !
      !call mp_bcast(flag,ionode_id,intra_image_comm)
!       print *, 'flag=',flag
!       flag = 1
      !
      !IF(flag.EQ.0) GO TO 1 
      !
      call mp_bcast(header,ionode_id, intra_image_comm)
      ! 
      if (ionode) write(*,*) "  DRIVER MODE: Message from server: ", header
      !
      if (trim(header) == "STATUS") then
         !
         if (ionode) then  ! does not  need init (well, maybe it should, just to check atom numbers and the like... )
            if (hasdata) then
               call writebuffer(socket,"HAVEDATA    ",MSGLEN, flag)
            else
               call writebuffer(socket,"READY       ",MSGLEN, flag)
            endif
         endif
         !
      else if (trim(header) == "POSDATA") then              
         !
         if (ionode) then        
           call readbuffer(socket, cellh, 9, flag)     ! HK/MCA : the new readbuffer interface takes care of the double size (no need to times 8)
           call readbuffer(socket, cellih, 9, flag)    ! HK/MCA : the new readbuffer interface takes care of the double size (no need to times 8)
           call readbuffer(socket, nat, flag)          ! HK/MCA : we don't need to put the length here due to the interface rule
         endif
         call mp_bcast(cellh,ionode_id, intra_image_comm)
         !call mp_bcast(cellih,ionode_id, intra_image_comm) ! HK/MCA : cellih only read not used (save communication)
         call mp_bcast(nat,ionode_id, intra_image_comm)
         !
         if (.not.allocated(combuf)) then
           allocate(combuf(3*nat))
         end if
         if (ionode) call readbuffer(socket, combuf, nat*3, flag)  ! HK/MCA : the new readbuffer interface takes care of the double size (no need to times 8)
         call mp_bcast(combuf,ionode_id, intra_image_comm)
         !
         if (ionode) write(*,*) "  DRIVER MODE: Received positions "
         if (firststep) then
            if (ionode) write(*,*) "  DRIVER MODE: Preparing first evaluation "
            call init_run()
            firststep=.false.
         end if 
         !
        !!DEBUG
        DO ia=1,3
          WRITE(stdout,'(3x,"cell: ",I7,3F18.8)') nfi,h(:,ia) 
        END DO
        !DO ia=1,3
        !  WRITE(stdout,'(3x,"ainv: ",I7,3F18.8)') nfi,ainv(:,ia) 
        !END DO
        !h=TRANSPOSE(cellh)
        !ainv=TRANSPOSE(cellih)
        !DO ia=1,3
        !  WRITE(stdout,'(3x,"h: ",I7,3F18.8)') nfi,h(:,ia) 
        !END DO
        !DO ia=1,3
        !  WRITE(stdout,'(3x,"celli: ",I7,3F18.8)') nfi,ainv(:,ia) 
        !END DO
        !!DEBUG
         !
         ! update number of steps, first step repeated in restart so skip ...
         !
         IF(.NOT.tfirst) nfi = nfi+1
         !
        !at=cellh/alat
         tau0 = RESHAPE(combuf, (/ 3 , nat /) )             
         !
         DO ia=1,nat
           WRITE(stdout,'(3x,"tau0: ",I7,3F18.10)') nfi,tau0(:,ia) 
         END DO
         !
         IF ( thdyn ) THEN
            !
            !h=TRANSPOSE(cellh) ! HK/MCA : using new fsocket the transpose is no longer needed due to the 1-D buffer
         h = TRANSPOSE(RESHAPE(cellh, (/ 3 , 3 /) ))
            !do j1 = 1,3
            !  jtmp = (j1-1)*3
            !  do j2 = 1,3
            !    h(j1,j2) = cellh(j2+jtmp)
            !  end do ! j2
            !end do ! j1
            do j1 = 1,3
              write(stdout,*) (h(j1,j2),j2=1,3)
            end do ! j1
            !
            IF( nbgrp > 1 ) THEN
               CALL mp_bcast( h, 0, inter_bgrp_comm )
            END IF
            !
            CALL newinit( h, iverbosity )
            !
            CALL newnlinit()
            !
            CALL formf( tfirst, eself ) ! should be called after newinit ...
            !
         END IF
         !
         ! increase time step
         !
         tps = tps + delt * au_ps
         !
         ! Call CP-BO to minimize electron ...
         !
         CALL cp_bo_emin(tfirst,tlast,enb,enbi,ccc,stress,bigr,iter)
         !
         tfirst=.FALSE. 
         !
         !DO ia=1,nat
         !  WRITE(stdout,'(3x,"fion: ",I7,3ES18.8)') nfi,fion(:,ia) 
         !END DO
         ! 
         combuf=RESHAPE(fion, (/ 3 * nat /) ) ! return force in atomic units
         ! 
         pot=etot                      ! return potential in atomic units
         vir=transpose(stress)*omega   ! return virial in atomic units and without the volume scaling
         !
         ! HK/MCA : todo check the virial convention (as in cellh)
         bufvir = RESHAPE(TRANSPOSE(vir), (/ 9 /) )
         !do j1 = 1,3
         !  jtmp = (j1-1)*3
         !  do j2 = 1,3
         !    bufvir(j2+jtmp) = vir(j2,j1) 
         !  end do ! j2
         !end do ! j1
         !
         stress_gpa = stress * au_gpa
         !DO ia=1,3
         !  WRITE(stdout,'(3x,"stress2: ",I7,3F18.8)') nfi,stress_gpa(:,ia) 
         !END DO
         WRITE(stdout,'(3x,"omega: ",I7,F18.8)') nfi,omega 
         !WRITE(stdout,'(3x,"pot: ",F18.10)') pot
         !       
         hasdata=.true.
         !
         !write save file every isave step ...
         !
         IF((nfi.GT.0).AND.(MOD(nfi,isave).EQ.(isave-1)).AND.(nfi.LT.nomore)) &
           CALL writefile( h, hold, nfi, c0_bgrp, cm_bgrp, taus,  &
                           tausm, vels, velsm, acc,  lambda, lambdam, descla, xnhe0,   &
                           xnhem, vnhe, xnhp0, xnhpm, vnhp, nhpcl, nhpdim, ekincm,&
                           xnhh0, xnhhm, vnhh, velh, fion, tps, z0t, f, rhor )
         !       
      else if (trim(header)=="GETFORCE") then
         !
         if (ionode) write(*,*) "  DRIVER MODE: Returning v,forces,stress "
         if (ionode) then      
            call writebuffer(socket,"FORCEREADY  ",MSGLEN, flag)            
            call writebuffer(socket,pot, flag) !HK/MCA : size taken care of by the interface
            call writebuffer(socket,nat, flag) !HK/MCA : size taken care of by the interface
            call writebuffer(socket,combuf,3*nat, flag)   ! HK/MCA : the new writebuffer interface takes care of the double size (no need to times 8)
            call writebuffer(socket,bufvir,9, flag)       !HK/MCA : 1-D buffer (to prevent confusion)
            nat=0
            call writebuffer(socket,nat, flag) !HK/MCA : size taken care of by the interface
         endif
         hasdata=.false.
         !
      endif
      !
      CALL stop_clock( 'driver_md' )
      !
      ! HK/MCA : if there is error then save and quit
      call mp_bcast(flag,ionode_id,intra_image_comm)
      IF(flag.EQ.0) GO TO 1
      !
    END DO driver_loop    
    !
    ! writefile ...
    !
1   CALL writefile( h, hold, nfi, c0_bgrp, cm_bgrp, taus,  &
                    tausm, vels, velsm, acc,  lambda, lambdam, descla, xnhe0,   &
                    xnhem, vnhe, xnhp0, xnhpm, vnhp, nhpcl, nhpdim, ekincm,&
                    xnhh0, xnhhm, vnhh, velh, fion, tps, z0t, f, rhor )
    ! 
    !
    !CALL error("leaving driver routine")
    !
    !
  END SUBROUTINE
