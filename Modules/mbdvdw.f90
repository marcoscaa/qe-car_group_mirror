!! Copyright (C) 2015 T. Markovich, M. Forsythe
!!
!! This library is free software; you can redistribute it and/or
!! modify it under the terms of the GNU Lesser General Public
!! License as published by the Free Software Foundation; either
!! version 3.0 of the License, or (at your option) any later version.
!!
!! This library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!
!! Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public
!! License along with this library.

module mbdvdw_module

  USE kinds,              ONLY: DP                 !double-precision kind (selected_real_kind(14,200))
  USE constants,          ONLY: pi                 !pi in double-precision
  USE ions_base,          ONLY: nat                !number of total atoms (all atomic species)
  USE ions_base,          ONLY: nsp                !number of unique atomic species
  USE ions_base,          ONLY: ityp               !ityp(i):= the type of i-th atom in stdin
  USE ions_base,          ONLY: atm                !atm(j) := name of j-th atomic species (3 characters): len(atm)=nsp
  USE io_global,          ONLY: stdout,ionode, ionode_id            !print/write argument for standard output (to output file)
  USE funct,              ONLY: get_iexch          !retrieves type of exchange utilized in functional
  USE funct,              ONLY: get_icorr          !retrieves type of correlation utilized in functional
  USE funct,              ONLY: get_igcx           !retrieves type of gradient correction to exchange utilized in functional
  USE funct,              ONLY: get_igcc           !retrieves type of gradient correction to correlation utilized in functional
  use control_flags,      only: iverbosity
  use cell_base,          only: omega
  USE input_parameters,   ONLY: mbd_vdw_vacuum
  USE input_parameters,   ONLY: mbd_vdw_supercell
  use input_parameters,   ONLY: vdw_self_consistent
  USE input_parameters,   ONLY: mbd_vdw_custom_beta
  use input_parameters,   ONLY: mbd_vdw_beta
  USE input_parameters,   ONLY: mbd_vdw_econv_thr
  USE input_parameters,   ONLY: mbd_vdw_n_quad_pts
  use input_parameters,   ONLY: mbd_vdw_forces
  use input_parameters,   ONLY: vdw_debug
  use quadrature_grid_module,    only: casimir_omega
  use quadrature_grid_module,    only: casimir_omega_weight
  use quadrature_grid_module,    only: generate_grid
  USE fft_base,           ONLY: dffts              !FFT derived data type
  USE fft_base,           ONLY: dfftp              !FFT derived data type
  USE mp,                 ONLY: mp_sum             !MPI collection with sum
  use mp,                 ONLY: mp_bcast
  use mp,                 ONLY: mp_barrier
  use mp,                 ONLY: mp_set_displs
  USE parallel_include                             !MPI header
  USE mp_bands,           ONLY: nproc_bgrp         !number of processors
  USE mp_bands,           ONLY: me_bgrp            !processor number (0,1,...,nproc_bgrp-1)
  USE mp_bands,           ONLY: intra_bgrp_comm    !standard MPI communicator
  !use mpi
  USE mp_images,          ONLY: nproc_image        !number of processors
  USE mp_images,          ONLY: me_image           !processor number (0,1,...,nproc_image-1)
  USE mp_images,          ONLY: intra_image_comm   !standard MPI communicator
  use mp_images,          only: root_image
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Imports: TSVDW
  ! We're going to call TSVDW_CALCULATE, but because mbdvdw is true, tsvdw only generates the volume
  ! and derivatives of the volumes. We'll need the GetVdWParam method to grab all the free atom
  ! quantities.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use tsvdw_module, ONLY: GetVdWParam
  use tsvdw_module, ONLY: tsvdw_calculate, tsvdw_initialize
  use tsvdw_module, ONLY: tsvdw_cleanup
  use tsvdw_module, ONLY: VefftsvdW, dveffdr, dveffdh, vfree
  use tsvdw_module, ONLY: tsvdw_cleanup_post_mbd
  use tsvdw_module, ONLY: somegaA, NsomegaA, dveffadn

  implicit none
  save

integer, parameter :: check = 60

INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: nstates       !number of atoms per processor
INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: sdispls       !send displacement (offset) array
INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: rdispls       !receive displacement (offset) array
INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: sendcount     !send count array
INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: recvcount     !receive count array
INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: istatus       !MPI status array

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Variable Declarations: Integration Grid
  ! npts is an integer that tells us the number of integration points for the numerical integral
  ! Casimir_omega is an array of length npts that contains the location of points to evaluate the integral on
  ! Casimir_omega_weight is an array of length npts that contains the weights for the numerical integral
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER, private :: npts
  INTEGER, PRIVATE :: nr1,nr2,nr3                              !real space grid dimensions (global first, second, and third dimensions of the 3D grid)
  INTEGER, PRIVATE :: nr1r,nr2r,nr3r                           !reduced real space grid dimensions (global first, second, and third dimensions of the 3D grid)

  !real(dp), dimension(:), allocatable, private :: casimir_omega
  !real(dp), dimension(:), allocatable, private :: casimir_omega_weight

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Variable Declarations: VdW Radii at various levels of theory and their derivatives.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(dp), dimension(:),     allocatable, private :: R_vdw_free

  real(dp), dimension(:),     allocatable, private :: R_TS_VdW
  real(dp), dimension(:,:,:), allocatable, private :: dR_TS_VdWdR
  real(dp), dimension(:,:,:), allocatable, private :: dR_TS_VdWdh
  real(dp), dimension(:, :),  allocatable, private :: dR_TS_VdWdV

  real(dp), dimension(:),     allocatable, private :: R_MBD_VdW
  real(dp), dimension(:,:,:), allocatable, private :: dR_MBD_VdWdR
  real(dp), dimension(:,:,:), allocatable, private :: dR_MBD_VdWdh
  real(dp), dimension(:, :),  allocatable, private :: dR_MBD_VdWdV

  real(dp), dimension(:),     allocatable, private :: R_MBD_VdW_sl
  real(dp), dimension(:,:,:), allocatable, private :: dR_MBD_VdWdR_sl
  real(dp), dimension(:,:,:), allocatable, private :: dR_MBD_VdWdh_sl
  real(dp), dimension(:,:),   allocatable, private :: dR_MBD_VdWdV_sl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Variable Declarations: Screened oscillator frequency
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(dp), dimension(:),     allocatable, private :: C6_free

  real(dp), dimension(:),     allocatable, private :: omega_scs
  real(dp), dimension(:,:,:), allocatable, private :: domegadR
  real(dp), dimension(:,:,:), allocatable, private :: domegadh
  real(dp), dimension(:, :),  allocatable, private :: domegadV

  real(dp), dimension(:),     allocatable, private :: omega_scs_sl
  real(dp), dimension(:,:,:), allocatable, private :: domegadR_sl
  real(dp), dimension(:,:,:), allocatable, private :: domegadh_sl
  real(dp), dimension(:,:),   allocatable, private :: domegadV_sl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Variable Declarations: Polarizabilities at various levels of theory and their derivatives
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(dp), dimension(:),     allocatable, private :: alpha_free

  real(dp), dimension(:),     allocatable, private :: alpha_ts
  real(dp), dimension(:,:,:), allocatable, private :: dalpha_tsdR
  real(dp), dimension(:,:,:), allocatable, private :: dalpha_tsdh
  real(dp), dimension(:, :),  allocatable, private :: dalpha_tsdV

  real(dp), dimension(:),     allocatable, private :: alpha_0
  real(dp), dimension(:,:,:), allocatable, private :: dalpha_0dR
  real(dp), dimension(:,:,:), allocatable, private :: dalpha_0dh
  real(dp), dimension(:, :),  allocatable, private :: dalpha_0dV

  real(dp), dimension(:),     allocatable, private :: alpha_0_sl
  real(dp), dimension(:,:,:), allocatable, private :: dalpha_0dR_sl
  real(dp), dimension(:,:,:), allocatable, private :: dalpha_0dh_sl
  real(dp), dimension(:,:),   allocatable, private :: dalpha_0dV_sl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Variable Declarations: VdW Radii at various levels of theory and their derivatives.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(dp), dimension(:),     allocatable, private :: sigma
  real(dp), dimension(:,:,:), allocatable, private :: dsigmadR
  real(dp), dimension(:,:,:), allocatable, private :: dsigmadh
  real(dp), dimension(:, :),  allocatable, private :: dsigmadV

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Variable Declarations: Working variables
  ! molecular_polarizability: This is the 3x3 tensor that's been contracted
  !                                      to be the polarizability in each cartesian direction
  ! A_matrix: This is the 3Nx3N frequency dependent matrix that we'll
  !                populate by (A + T)^-1. In the nomenclature of the MBD forces paper
  !                we'll understand this paper as having a bar over top
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(dp), dimension(:,:),     allocatable, private :: A_matrix
  real(dp), dimension(:,:,:,:), allocatable, private :: dA_matrixdR
  real(dp), dimension(:,:,:,:), allocatable, private :: dA_matrixdh
  real(dp), dimension(:,:,:),   allocatable, private :: dA_matrixdV

  real(dp), dimension(:, :),    allocatable, private :: Cpq
  real(dp), dimension(:,:,:,:), allocatable, private :: dCpqdR
  real(dp), dimension(:,:,:,:), allocatable, private :: dCpqdh
  real(dp), dimension(:,:,:),   allocatable, private :: dCpqdV

  real(dp), private :: interacting_energy
  real(dp), private :: non_interacting_energy
  real(dp), public  :: EmbdvdW

  real(dp), dimension(:,:), allocatable, private :: int_FmbdVdW
  real(dp), dimension(:,:), allocatable, private :: nonint_FmbdVdW
  real(dp), dimension(:,:), allocatable, public  :: FmbdVdW

  real(dp), dimension(:,:), allocatable, private :: int_HmbdVdW
  real(dp), dimension(:,:), allocatable, private :: nonint_HmbdVdW
  real(dp), dimension(:,:), allocatable, public  :: HmbdVdW

  real(dp), dimension(:),   allocatable, private :: int_Uprefactor
  real(dp), dimension(:),   allocatable, private :: nonint_Uprefactor
  real(dp), dimension(:),   allocatable, private :: Uprefactor
  real(dp), dimension(:),   allocatable, public  :: UmbdvdW

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Variable Declarations: Simulation parameters
  ! beta: Damping factor associated with the functional used. Only defined for
  !         PBE, beta=0.83
  !         PBE0, beta=0.85
  !         HSE, beta=0.85
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(dp), private               :: beta
  real(dp), private               :: dip_cutoff
  real(dp), private               :: supercell_cutoff
  logical,  dimension(3), private :: mbdvdw_vacuum
  logical                         :: sc_only
  logical                         :: do_forces
  logical                         :: mbd_conv_elec
  logical                         :: mbd_first_step
  real(dp)                        :: omega_to_print

  REAL(DP), PRIVATE :: h_(3,3)                                 !HK: PW and CP use difference var for cell (at) and (h), so this is wrapper between the two
  REAL(DP), PRIVATE :: ainv_(3,3)                              !HK: PW and CP use difference var for cell (at) and (h), so this is wrapper between the two


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Super Lattice Variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, private                                 :: nat_sl
  real(dp), dimension(3, 3),  private              :: h_sl
  real(dp), dimension(3, 3),  private              :: ainv_sl
  real(dp), dimension(:, :),  allocatable, private :: tau_sl, tau_sl_s
  real(dp), dimension(:,:),   allocatable          :: tau
  real(dp), dimension(:,:),   allocatable, private :: tau_s
  integer                                          :: sl_i, sl_j, sl_k
  integer, dimension(3)                            :: sl_mult

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Parallel variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type pairs
    integer p
    integer q
    real(dp) rpq(3)
    real(dp) rvdw
    logical unique
    integer CPU
  end type pairs

  integer, private                           :: num_pairs
  type(pairs), dimension(:), allocatable     :: unique_pairs
  type(pairs), dimension(:), allocatable     :: pairs_scs
  integer, dimension(:), allocatable         :: f_cpu_id
  integer                                    :: max_proc_forces
  integer                                    :: max_proc_pairs, me

  character(len=6) :: pos="append"
  character(len=8) :: stat="unknown"
  character(len=5) :: act="write"
  character(len=50)               :: dir
  logical :: mbd_debug_dh = .false.
  logical :: mbd_debug_dr = .false.
  logical :: mbd_write_check = .false.
  integer, dimension(:), allocatable :: n_pairs, n_comps
  integer, dimension(:), allocatable :: dh_displs, dr_displs, f_displs



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Method Prototypes
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  public :: mbdvdw_initialize
  public :: mbdvdw_calculate
  public :: mbdvdw_cleanup

  private :: mbdvdw_effqts
  private :: mbdvdw_SCS
  private :: mbdvdw_construct_hamiltonian
  private :: mbdvdw_TSR
  private :: mbdvdw_TLR
  private :: mbdvdw_para_init
  private :: mbdvdw_para_init_sl
  private :: mbdvdw_noninteracting_energy
  private :: mbdvdw_interacting_energy
  private :: mbdvdw_calculate_screened_pol
  private :: mbdvdw_pbc

  CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This method will initialize the MBD calc by initalizing the associated constants
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mbdvdw_initialize()
  implicit none
  integer :: i_atom, Ndim
  CALL start_clock('mbd_init')

  dir = ""
  dip_cutoff = 180.0_DP
  mbd_first_step = .true.
  mbd_conv_elec = .false.
  supercell_cutoff = mbd_vdw_supercell/0.5291772490_DP
  mbdvdw_vacuum = mbd_vdw_vacuum
  do_forces = mbd_vdw_forces
  mbd_debug_dr = vdw_debug

  if(.not.allocated(R_vdw_free))      allocate(R_vdw_free(nat))
  if(.not.allocated(R_TS_VdW))        allocate(R_TS_VdW(nat));            R_TS_VdW    = 0.0_DP
  if(.not.allocated(R_MBD_VdW))       allocate(R_MBD_VdW(nat));            R_MBD_VdW    = 0.0_DP
  if(.not.allocated(C6_free))         allocate(C6_free(nat));              C6_free = 0.0_DP
  if(.not.allocated(tau))             allocate(tau(3, nat))
  if(.not.allocated(vfree))           allocate(vfree(nat))
  if(.not.allocated(alpha_free))      allocate(alpha_free(nat))
  if(.not.allocated(alpha_ts))        allocate(alpha_ts(nat));             alpha_ts = 0.0_DP
  if(.not.allocated(alpha_0))         allocate(alpha_0(nat));              alpha_0 = 0.0_DP
  if(.not.allocated(omega_scs))       allocate(omega_scs(nat));            omega_scs = 0.0_DP
  if(.not.allocated(sigma))           allocate(sigma(nat));                sigma = 0.0_DP
  if(.not.allocated(dR_TS_VdWdV))     allocate(dR_TS_VdWdV(nat, nat)); dR_TS_VdWdV = 0.0_DP
  if(.not.allocated(dR_MBD_VdWdV))    allocate(dR_MBD_VdWdV(nat, nat));dR_MBD_VdWdv = 0.0_DP
  if(.not.allocated(domegadV))        allocate(domegadV(nat, nat));    domegadV = 0.0_DP
  if(.not.allocated(dalpha_tsdV))     allocate(dalpha_tsdV(nat, nat)); dalpha_tsdV = 0.0_DP
  if(.not.allocated(dalpha_0dV))      allocate(dalpha_0dV(nat, nat));  dalpha_0dV = 0.0_DP
  if(.not.allocated(dsigmadV))        allocate(dsigmadV(nat, nat));    dsigmadV = 0.0_DP
  if(.not.allocated(dR_TS_VdWdR))     allocate(dR_TS_VdWdR(nat, nat, 3)); dR_TS_VdWdR = 0.0_DP
  if(.not.allocated(dR_TS_VdWdh))     allocate(dR_TS_VdWdh(nat, 3  , 3)); dR_TS_VdWdh = 0.0_DP
  if(.not.allocated(dR_MBD_VdWdR))    allocate(dR_MBD_VdWdR(nat, nat, 3)); dR_MBD_VdWdR = 0.0_DP
  if(.not.allocated(dR_MBD_VdWdh))    allocate(dR_MBD_VdWdh(nat, 3  , 3)); dR_MBD_VdWdh = 0.0_DP
  if(.not.allocated(domegadR))        allocate(domegadR(nat, nat, 3));     domegadR = 0.0_DP
  if(.not.allocated(domegadh))        allocate(domegadh(nat, 3, 3));       domegadh = 0.0_DP
  if(.not.allocated(dalpha_tsdR))     allocate(dalpha_tsdR(nat, nat, 3));  dalpha_tsdR = 0.0_DP
  if(.not.allocated(dalpha_tsdh))     allocate(dalpha_tsdh(nat, 3  , 3));  dalpha_tsdh = 0.0_DP
  if(.not.allocated(dalpha_0dR))      allocate(dalpha_0dR(nat, nat, 3));   dalpha_0dR = 0.0_DP
  if(.not.allocated(dalpha_0dh))      allocate(dalpha_0dh(nat, 3, 3));     dalpha_0dh = 0.0_DP
  if(.not.allocated(dsigmadR))        allocate(dsigmadR(nat, nat, 3));     dsigmadR = 0.0_DP
  if(.not.allocated(dsigmadh))        allocate(dsigmadh(nat, 3  , 3));     dsigmadh = 0.0_DP

  if(.not.allocated(nonint_Uprefactor))  allocate(nonint_Uprefactor(nat));       nonint_Uprefactor = 0.0_DP
  if(.not.allocated(int_Uprefactor))     allocate(int_Uprefactor(nat));          int_Uprefactor = 0.0_DP
  if(.not.allocated(Uprefactor))         allocate(Uprefactor(nat));              Uprefactor = 0.0_DP

  if(.not.allocated(nonint_FmbdVdW))  allocate(nonint_FmbdVdW(nat, 3));    nonint_FmbdVdW = 0.0_DP
  if(.not.allocated(int_FmbdVdW))     allocate(int_FmbdVdW(nat, 3));       int_FmbdVdW = 0.0_DP
  if(.not.allocated(FmbdVdW))         allocate(FmbdVdW(nat, 3));           FmbdVdW = 0.0_DP

  if(.not.allocated(nonint_HmbdVdW))  allocate(nonint_HmbdVdW(3, 3));      nonint_HmbdVdW = 0.0_DP
  if(.not.allocated(int_HmbdVdW))     allocate(int_HmbdVdW(3, 3));         int_HmbdVdW = 0.0_DP
  if(.not.allocated(HmbdVdW))         allocate(HmbdVdW(3, 3));             HmbdVdW = 0.0_DP

  if(.not.allocated(tau_s))           allocate(tau_s(3, nat));            tau_s = 0.0_DP
  if(.not.allocated(n_pairs))         allocate(n_pairs(nproc_image));     n_pairs = 0
  if(.not.allocated(n_comps))         allocate(n_comps(nproc_image));     n_comps = 0

  Ndim=dfftp%nnr
  if(.not.allocated(UmbdvdW)) ALLOCATE(UmbdvdW(Ndim)); UmbdvdW=0.0_DP
  call generate_grid(mbd_vdw_n_quad_pts, npts)

  do i_atom=1, nat, 1
    call GetVdWParam(atm(ityp(i_atom)), C6_free(i_atom), alpha_free(i_atom), R_vdw_free(i_atom))
  end do

  if(.not.mbd_vdw_custom_beta) then
    IF (get_iexch().EQ.1.AND.get_icorr().EQ.4.AND.get_igcx().EQ.3.AND.get_igcc().EQ.4) THEN
      beta=0.83_DP !PBE=sla+pw+pbx+pbc
    ELSE IF (get_iexch().EQ.6.AND.get_icorr().EQ.4.AND.get_igcx().EQ.8.AND.get_igcc().EQ.4) THEN
      beta=0.85_DP !PBE0=pb0x+pw+pb0x+pbc !RAD/BS: This line will not work in CP unless PBE0 code update funct.f90...
    ELSE
      CALL errore('mbdvdw','mbd-vdW sR parameter only available for PBE and PBE0 functionals...',1)
    END IF
    else
    beta = mbd_vdw_beta
    end if

  CALL stop_clock('mbd_init')
  end subroutine mbdvdw_initialize

  subroutine mbdvdw_zeroout()
  R_TS_VdW    = 0.0_DP
  R_MBD_VdW    = 0.0_DP
  alpha_ts = 0.0_DP
  alpha_0 = 0.0_DP
  omega_scs = 0.0_DP
  sigma = 0.0_DP

  if(vdw_self_consistent) then
    dR_TS_VdWdV = 0.0_DP
    dR_MBD_VdWdv = 0.0_DP
    domegadV = 0.0_DP
    dalpha_tsdV = 0.0_DP
    dalpha_0dV = 0.0_DP
    dsigmadV = 0.0_DP
    end if

  if(do_forces) then
    dR_TS_VdWdR = 0.0_DP
    dR_TS_VdWdh = 0.0_DP
    dR_MBD_VdWdR = 0.0_DP
    dR_MBD_VdWdh = 0.0_DP
    domegadR = 0.0_DP
    domegadh = 0.0_DP
    dalpha_tsdR = 0.0_DP
    dalpha_tsdh = 0.0_DP
    dalpha_0dR = 0.0_DP
    dalpha_0dh = 0.0_DP
    dsigmadR = 0.0_DP
    dsigmadh = 0.0_DP
  end if

  nonint_Uprefactor = 0.0_DP
  int_Uprefactor = 0.0_DP
  Uprefactor = 0.0_DP

  nonint_FmbdVdW = 0.0_DP
  int_FmbdVdW = 0.0_DP
  FmbdVdW = 0.0_DP

  nonint_HmbdVdW = 0.0_DP
  int_HmbdVdW = 0.0_DP
  HmbdVdW = 0.0_DP

  tau_s = 0.0_DP
  n_pairs = 0
  n_comps = 0

  end subroutine mbdvdw_zeroout

  subroutine mbdvdw_cleanup_postscs()

  if(allocated(A_matrix))        deallocate(A_matrix)
  if(allocated(dA_matrixdR))     deallocate(dA_matrixdR)
  if(allocated(dA_matrixdh))     deallocate(dA_matrixdh)
  if(allocated(dA_matrixdV))     deallocate(dA_matrixdV)

  end subroutine mbdvdw_cleanup_postscs

  subroutine mbdvdw_cleanup()

  if(allocated(Cpq))              deallocate(Cpq)
  if(allocated(dCpqdR))           deallocate(dCpqdR)
  if(allocated(dCpqdh))           deallocate(dCpqdh)

  if(allocated(VefftsvdW))             deallocate(VefftsvdW)
  if(allocated(dveffdr))          deallocate(dveffdr)
  if(allocated(dveffdh))          deallocate(dveffdh)
  end subroutine mbdvdw_cleanup

  subroutine mbdvdw_finalize()
  if(allocated(tau_sl))           deallocate(tau_sl)

  if(allocated(R_TS_VdW))        deallocate(R_TS_VdW)
  if(allocated(dR_TS_VdWdR))     deallocate(dR_TS_VdWdR)
  if(allocated(dR_TS_VdWdh))     deallocate(dR_TS_VdWdh)
  if(allocated(dR_TS_VdWdV))     deallocate(dR_TS_VdWdV)

  if(allocated(alpha_ts))        deallocate(alpha_ts)
  if(allocated(dalpha_tsdR))     deallocate(dalpha_tsdR)
  if(allocated(dalpha_tsdh))     deallocate(dalpha_tsdh)
  if(allocated(dalpha_tsdV))     deallocate(dalpha_tsdV)

  if(allocated(R_MBD_VdW))       deallocate(R_MBD_VdW)
  if(allocated(dR_MBD_VdWdR))    deallocate(dR_MBD_VdWdR)
  if(allocated(dR_MBD_VdWdh))    deallocate(dR_MBD_VdWdh)
  if(allocated(dR_MBD_VdWdV))    deallocate(dR_MBD_VdWdV)

  if(allocated(domegadR))        deallocate(domegadR)
  if(allocated(domegadh))        deallocate(domegadh)
  if(allocated(domegadV))        deallocate(domegadV)

  if(allocated(alpha_0))         deallocate(alpha_0)
  if(allocated(dalpha_0dR))      deallocate(dalpha_0dR)
  if(allocated(dalpha_0dh))      deallocate(dalpha_0dh)
  if(allocated(dalpha_0dV))      deallocate(dalpha_0dV)

  if(allocated(omega_scs))       deallocate(omega_scs)
  if(allocated(sigma))           deallocate(sigma)
  if(allocated(dsigmadR))        deallocate(dsigmadR)
  if(allocated(dsigmadh))        deallocate(dsigmadh)
  if(allocated(dsigmadV))        deallocate(dsigmadV)

  if(allocated(alpha_0_sl))       deallocate(alpha_0_sl)
  if(allocated(dalpha_0dR_sl))    deallocate(dalpha_0dR_sl)
  if(allocated(dalpha_0dh_sl))    deallocate(dalpha_0dh_sl)

  if(allocated(omega_scs_sl))     deallocate(omega_scs_sl)
  if(allocated(domegadR_sl))      deallocate(domegadR_sl)
  if(allocated(domegadh_sl))      deallocate(domegadh_sl)

  if(allocated(R_MBD_VdW_sl))     deallocate(R_MBD_VdW_sl)
  if(allocated(dR_MBD_VdWdR_sl))  deallocate(dR_MBD_VdWdR_sl)
  if(allocated(dR_MBD_VdWdh_sl))  deallocate(dR_MBD_VdWdh_sl)

  if(allocated(vfree))            deallocate(vfree)
  if(allocated(alpha_free))       deallocate(alpha_free)
  if(allocated(C6_free))          deallocate(C6_free)
  if(allocated(R_vdw_free))       deallocate(R_vdw_free)

  if(allocated(nonint_FmbdVdW))   deallocate(nonint_FmbdVdW)
  if(allocated(int_FmbdVdW))      deallocate(int_FmbdVdW)
  if(allocated(FmbdVdW))          deallocate(FmbdVdW)

  if(allocated(nonint_HmbdVdW))   deallocate(nonint_HmbdVdW)
  if(allocated(int_HmbdVdW))      deallocate(int_HmbdVdW)
  if(allocated(HmbdVdW))          deallocate(HmbdVdW)

  if(allocated(nonint_Uprefactor))   deallocate(nonint_Uprefactor)
  if(allocated(int_Uprefactor))      deallocate(int_Uprefactor)
  if(allocated(Uprefactor))          deallocate(Uprefactor)
  if(allocated(UmbdVdW))             deallocate(UmbdvdW)

  end subroutine mbdvdw_finalize

  subroutine mbdvdw_pbc(tauin, hin, ainvin, tau_sout)
  implicit none
  integer              :: i_atom
  REAL(DP), intent(inout) :: tauin(3, nat)
  REAL(DP), intent(in) :: hin(3, 3)
  real(dp), intent(in) :: ainvin(3, 3)
  real(dp), intent(out):: tau_sout(3, nat)

  do i_atom = 1, nat, 1
    tau_sout(1,i_atom)=ainvin(1,1)*tauin(1,i_atom)+ainvin(1,2)*tauin(2,i_atom)+ainvin(1,3)*tauin(3,i_atom)   ! s = h_^-1 r
    tau_sout(2,i_atom)=ainvin(2,1)*tauin(1,i_atom)+ainvin(2,2)*tauin(2,i_atom)+ainvin(2,3)*tauin(3,i_atom)   ! s = h_^-1 r
    tau_sout(3,i_atom)=ainvin(3,1)*tauin(1,i_atom)+ainvin(3,2)*tauin(2,i_atom)+ainvin(3,3)*tauin(3,i_atom)   ! s = h_^-1 r
    !
    ! tau_sout(1,i_atom)=tau_sout(1,i_atom)-floor(tau_sout(1,i_atom))   ! impose PBC on s in range: [0,1)
    ! tau_sout(2,i_atom)=tau_sout(2,i_atom)-floor(tau_sout(2,i_atom))   ! impose PBC on s in range: [0,1)
    ! tau_sout(3,i_atom)=tau_sout(3,i_atom)-floor(tau_sout(3,i_atom))   ! impose PBC on s in range: [0,1)
    !
    tauin(1,i_atom)=hin(1,1)*tau_sout(1,i_atom)+hin(1,2)*tau_sout(2,i_atom)+hin(1,3)*tau_sout(3,i_atom)   ! r = h_ s
    tauin(2,i_atom)=hin(2,1)*tau_sout(1,i_atom)+hin(2,2)*tau_sout(2,i_atom)+hin(2,3)*tau_sout(3,i_atom)   ! r = h_ s
    tauin(3,i_atom)=hin(3,1)*tau_sout(1,i_atom)+hin(3,2)*tau_sout(2,i_atom)+hin(3,3)*tau_sout(3,i_atom)   ! r = h_ s
  end do

  end subroutine mbdvdw_pbc

  subroutine mbdvdw_effqts(omega)
  implicit none
  real(dp), intent(in) :: omega
  real(dp) :: omega_free, gamma, xi, pade_approx, lambda
  integer :: i_atom, s, i, ias
  integer, parameter :: sigmats_check     = 100
  integer, parameter :: dsigmats_dr_check = 101
  integer, parameter :: dsigmats_dh_check = 102
  integer, parameter :: alphats_check     = 103
  integer, parameter :: dalphats_dr_check = 104
  integer, parameter :: dalphats_dh_check = 105
  integer, parameter :: rts_check         = 106
  integer, parameter :: drts_dr_check     = 107
  integer, parameter :: drts_dh_check     = 108
  integer, parameter :: veff_check        = 109
  integer, parameter :: dveff_dr_check    = 110
  integer, parameter :: dveff_dh_check    = 111

  call start_clock('mbd_effqts')

  if(vdw_self_consistent) then
    dsigmadV = 0.0_DP
    dR_TS_VdWdV = 0.0_DP
    dalpha_tsdV = 0.0_DP
    end if

  if(do_forces) then
    dsigmadR = 0.0_DP
    dR_TS_VdWdR = 0.0_DP
    dalpha_tsdR = 0.0_DP
    dsigmadH = 0.0_DP
    dR_TS_VdWdH = 0.0_DP
    dalpha_tsdH = 0.0_DP
    end if

  do i_atom=1, nat, 1
    ias = ityp(i_atom)
    omega_free = ((4.0_DP/3.0_DP)*C6_free(i_atom)/(alpha_free(i_atom)**2.0_DP))
    ! Pade Approx
    pade_approx = 1.0_DP/(1.0_DP + (omega/omega_free)**2.0_DP )
    ! Computes sigma
    gamma = (1.0_DP/3.0_DP)*dsqrt(2.0_DP/pi)*pade_approx*(alpha_free(i_atom)/vfree(ias))
    gamma = gamma**(1.0_DP/3.0_DP)
    sigma(i_atom) = gamma*VefftsvdW(i_atom)**(1.0_DP/3.0_DP)

    ! Computes R_TS
    xi = R_vdw_free(i_atom)/(vfree(ias))**(1.0_DP/3.0_DP)
    R_TS_VdW(i_atom) = xi*VefftsvdW(i_atom)**(1.0_DP/3.0_DP)

    ! Computes alpha_ts
    lambda = pade_approx*alpha_free(i_atom)/vfree(ias)
    alpha_ts(i_atom) = lambda*VefftsvdW(i_atom)


    if(do_forces.or.vdw_self_consistent) then
      do s=1, nat, 1
        if(vdw_self_consistent.and.(i_atom.eq.s)) then
          dsigmadV(i_atom, s) = gamma/(3.0_DP*VefftsvdW(i_atom)**(2.0_DP/3.0_DP))
          dR_TS_VdWdV(i_atom, s) = xi/(3.0_DP*VefftsvdW(i_atom)**(2.0_DP/3.0_DP))
          dalpha_tsdV(i_atom, s) = lambda
          end if
        if(do_forces) then
          do i=1,3,1
            dsigmadR(i_atom, s, i) = gamma/(3.0_DP*VefftsvdW(i_atom)**(2.0_DP/3.0_DP))*dveffdR(i_atom, s, i)
            dR_TS_VdWdR(i_atom, s, i) = xi/(3.0_DP*VefftsvdW(i_atom)**(2.0_DP/3.0_DP))*dveffdR(i_atom, s, i)
            dalpha_tsdR(i_atom, s, i) = lambda*dveffdR(i_atom, s, i)
            if (s.le.3) then
              dsigmadh(i_atom, s, i) = gamma/(3.0_DP*VefftsvdW(i_atom)**(2.0_DP/3.0_DP))*dveffdh(i_atom, s, i)
              dR_TS_VdWdh(i_atom, s, i) = xi/(3.0_DP*VefftsvdW(i_atom)**(2.0_DP/3.0_DP))*dveffdh(i_atom, s, i)
              dalpha_tsdh(i_atom, s, i) = lambda*dveffdh(i_atom, s, i)
              end if
            end do
          end if
        end do
      end if
    end do

  if(mbd_debug_dr.or.mbd_debug_dh) then
    open (unit=sigmats_check, file=trim(dir)//"sigma_ts.tmp",action=act,status=stat,position=pos)
    open (unit=rts_check, file=trim(dir)//"r_ts.tmp",action=act,status=stat,position=pos)
    open (unit=alphats_check, file=trim(dir)//"alpha_ts.tmp",action=act,status=stat,position=pos)
    open (unit=veff_check, file=trim(dir)//"VefftsvdW.tmp",action=act,status=stat,position=pos)
    do i_atom=1,nat,1
      write(sigmats_check, *) omega, i_atom, sigma(i_atom)
      write(rts_check, *) omega, i_atom, R_TS_VdW(i_atom)
      write(alphats_check, *) omega, i_atom, alpha_ts(i_atom)
      write(veff_check, *) omega, i_atom, VefftsvdW(i_atom)
      end do
    close(sigmats_check)
    close(rts_check)
    close(alphats_check)
    close(veff_check)
    end if

  if(mbd_debug_dr) then
    open (unit=dsigmats_dr_check, file=trim(dir)//"dsigmats_dr.tmp",action=act,status=stat,position=pos)
    open (unit=drts_dr_check,     file=trim(dir)//"drts_dr.tmp",action=act,status=stat,position=pos)
    open (unit=dalphats_dr_check, file=trim(dir)//"dalpha_tsdr.tmp",action=act,status=stat,position=pos)
    open (unit=dveff_dr_check, file=trim(dir)//"dveff_dr.tmp",action=act,status=stat,position=pos)
    do i_atom=1,nat,1
      do s=1,3,1
        do i=1,3,1
          write(dsigmats_dr_check, *) omega, i_atom, s, i, dsigmadr(i_atom, s, i)
          write(drts_dr_check, *)     omega, i_atom, s, i, dR_TS_VdWdr(i_atom, s, i)
          write(dalphats_dr_check, *) omega, i_atom, s, i, dalpha_tsdr(i_atom, s, i)
          write(dveff_dr_check, *)    omega, i_atom, s, i, dveffdr(i_atom, s, i)
          end do
        end do
      end do
    close(dsigmats_dr_check)
    close(drts_dr_check)
    close(dalphats_dr_check)
    close(dveff_dr_check)
    end if

  if(mbd_debug_dh) then
    open (unit=dsigmats_dh_check, file=trim(dir)//"dsigmats_dh.tmp",action=act,status=stat,position=pos)
    open (unit=drts_dh_check,     file=trim(dir)//"drts_dh.tmp",action=act,status=stat,position=pos)
    open (unit=dalphats_dh_check, file=trim(dir)//"dalpha_tsdh.tmp",action=act,status=stat,position=pos)
    open (unit=dveff_dh_check, file=trim(dir)//"dveff_dh.tmp",action=act,status=stat,position=pos)
    do i_atom=1,nat,1
      do s=1,3,1
        do i=1,3,1
          write(dsigmats_dh_check, *) omega, i_atom, s, i, dsigmadh(i_atom, s, i)
          write(drts_dh_check, *)     omega, i_atom, s, i, dR_TS_VdWdh(i_atom, s, i)
          write(dalphats_dh_check, *) omega, i_atom, s, i, dalpha_tsdh(i_atom, s, i)
          write(dveff_dh_check, *)    omega, i_atom, s, i, dveffdh(i_atom, s, i)
          end do
        end do
      end do
    close(dsigmats_dh_check)
    close(drts_dh_check)
    close(dalphats_dh_check)
    close(dveff_dh_check)
    end if
  call stop_clock('mbd_effqts')

  end subroutine mbdvdw_effqts

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! initializes the parallel part of the mbd computation ! TODO: rewrite
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mbdvdw_para_init()
  implicit none
  integer :: p, q, counter
  integer :: n_pairs_per_proc, cpu, i, j, k

  if(.not.allocated(nstates  )) ALLOCATE(nstates(nproc_image)); nstates=0
  if(.not.allocated(sdispls  )) ALLOCATE(sdispls(nproc_image)); sdispls=0
  if(.not.allocated(sendcount)) ALLOCATE(sendcount(nproc_image)); sendcount=0
  if(.not.allocated(rdispls  )) ALLOCATE(rdispls(nproc_image)); rdispls=0
  if(.not.allocated(recvcount)) ALLOCATE(recvcount(nproc_image)); recvcount=0
  if(.not.allocated(istatus  )) ALLOCATE(istatus(nproc_image)); istatus=0

  num_pairs = (nat*nat-nat)/2
  n_pairs = 0
  if(.not.allocated(pairs_scs))        allocate(pairs_scs(num_pairs))
  n_pairs_per_proc = floor(dble(num_pairs)/nproc_image)
  cpu = 0
  n_pairs = 0
  counter = 1

  do p = 1, nat, 1
    do q = p, nat, 1
      if(p.ne.q) then
        pairs_scs(counter)%p = p
        pairs_scs(counter)%q = q
        counter = counter + 1
      end if
    end do
  end do
  do counter = 0, num_pairs-1, 1
    n_pairs(modulo(counter, nproc_image)+1) = n_pairs(modulo(counter, nproc_image)+1) + 1
    end do

  p = 1
  do counter = 1, num_pairs, 1
    pairs_scs(counter)%cpu = cpu
    if((counter.lt.num_pairs)) then
      if(p.eq.n_pairs(cpu+1)) then
        cpu = cpu + 1
        p=0
        end if
      end if
      p=p+1
    end do

  if(num_pairs.ge.nproc_image) then
    max_proc_pairs = nproc_image
    else
    max_proc_pairs = num_pairs
    end if

  ! Copied from tsvdw so that we can have the right parallelization scheme for wf forces
  me=me_image+1
  IF (nat.LE.nproc_image) THEN
    DO i=1,nat
      nstates(i)=1
    END DO
    ELSE
    k=0
  10  DO j=1,nproc_image
      nstates(j)=nstates(j)+1
      k=k+1
      IF (k.GE.nat) GO TO 20
    END DO
    IF (k.LT.nat) GO TO 10
    END IF
    20 CONTINUE

  RETURN
  end subroutine mbdvdw_para_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! initializes the parallel part of the mbd computation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mbdvdw_para_init_sl()
  implicit none
  integer :: p, q, counter
  integer :: n_pairs_per_proc, cpu

  num_pairs = (nat_sl*nat_sl-nat_sl)/2
  n_pairs = 0
  if(.not.allocated(unique_pairs))        allocate(unique_pairs(num_pairs))
  n_pairs_per_proc = floor(dble(num_pairs)/nproc_image)
  cpu = 0
  n_pairs = 0
  counter = 1

  do p = 1, nat_sl, 1
    do q = p, nat_sl, 1
      if(p.ne.q) then
        unique_pairs(counter)%p = p
        unique_pairs(counter)%q = q
        counter = counter + 1
      end if
    end do
  end do
  do counter = 0, num_pairs-1, 1
    n_pairs(modulo(counter, nproc_image)+1) = n_pairs(modulo(counter, nproc_image)+1) + 1
    end do

  p = 1
  do counter = 1, num_pairs, 1
    unique_pairs(counter)%cpu = cpu
    if((counter.lt.num_pairs)) then
      if(p.eq.n_pairs(cpu+1)) then
        cpu = cpu + 1
        p=0
        end if
      end if
      p=p+1
    end do

  if(num_pairs.ge.nproc_image) then
    max_proc_pairs = nproc_image
    else
    max_proc_pairs = num_pairs
    end if

  RETURN

  end subroutine mbdvdw_para_init_sl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! initializes the parallel part of the mbd computation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mbdvdw_para_init_forces()
  implicit none
  integer :: p, s, counter
  integer :: cpu

  if(.not.allocated(f_cpu_id))        allocate(f_cpu_id(nat));  f_cpu_id = 0
  cpu = 0
  n_comps = 0

  counter = 1
  p = 1
  do counter = 0, nat-1, 1
    n_comps(modulo(counter, nproc_image)+1) = n_comps(modulo(counter, nproc_image)+1) + 1
    end do

  do s = 1, nat, 1
    f_cpu_id(s) = cpu
    if((s.lt.nat).and.(p.eq.n_comps(cpu+1))) then
      cpu = cpu + 1
      p = 0
      end if
    p = p + 1
    end do

  if(nat.ge.nproc_image) then
    max_proc_forces = nproc_image
    else
    max_proc_forces = nat
    end if

  end subroutine mbdvdw_para_init_forces

  subroutine mbdvdw_noninteracting_energy(energy, forcedR, forcedh, forcedV)
    implicit none
    real(dp), intent(out) :: energy
    real(dp), dimension(nat, 3), intent(out) :: forcedR
    real(dp), dimension(3,3), intent(out) :: forcedh
    real(dp), dimension(nat), intent(out) :: forcedV
    ! Local variable
    integer :: p

    call start_clock('mbd_nonint_energy')
    energy = 0.0_DP
    forcedR = 0.0_DP
    forcedh = 0.0_DP
    forcedV = 0.0_DP

    do p=1, nat, 1
      energy = energy + omega_scs(p)
      if(do_forces) then
        forcedR = forcedR + domegadR(p, :, :)
        forcedh = forcedh + domegadh(p, :, :)
      end if
      if(vdw_self_consistent) forcedV(:) = forcedV(:) + domegadV(p, :)
    end do

    call stop_clock('mbd_nonint_energy')
    return
  end subroutine mbdvdw_noninteracting_energy

  subroutine mbdvdw_interacting_energy(energy, forcedR, forcedh, forcedV)
    implicit none
    real(dp), intent(out) :: energy
    real(dp), dimension(nat, 3), intent(out) :: forcedR
    real(dp), dimension(3,3), intent(out) :: forcedh
    real(dp), dimension(nat), intent(out) :: forcedV
    real(dp), dimension(3*nat_sl, 3*nat_sl) :: temp, temp1
    integer :: num_negative, i_atom, s, i, j, counter
    integer, parameter :: eigs_check = 200

    ! lapack work variables
    integer :: LWORK, errorflag
    real(dp) :: WORK((3*nat_sl)*(3+(3*nat_sl)/2)), eigenvalues(3*nat_sl), eval_inv(3*nat_sl)

    call start_clock('mbd_int_energy')
    eigenvalues = 0.0_DP
    forcedR = 0.0_DP
    forcedh = 0.0_DP
    energy = 0.0_DP
    num_negative = 0
    forcedV = 0.0_DP

    errorflag=0
    LWORK=3*nat_sl*(3+(3*nat_sl)/2)
    call DSYEV('V', 'U', 3*nat_sl, Cpq, 3*nat_sl, eigenvalues, WORK, LWORK, errorflag)

    if(errorflag.eq.0) then
      do i_atom=1, 3*nat_sl, 1
        if(eigenvalues(i_atom).ge.0.0_DP) then
          eigenvalues(i_atom) = dsqrt(eigenvalues(i_atom))
          eval_inv(i_atom) = 1.0_DP/(2.0_DP*eigenvalues(i_atom))
          energy = energy + eigenvalues(i_atom)
          else
          num_negative = num_negative + 1
          end if
        end do

      if(num_negative.ge.1) then
        write(stdout, '(3X," WARNING: Found ", I3, " Negative Eigenvalues.")') num_negative
        end if
      end if
    energy = energy*nat/nat_sl

    ! Do self consistent calculation
    if(vdw_self_consistent) then
      counter = 1
      forcedV = 0.0_DP
      do s=1,nat,1
        if(me_image.eq.f_cpu_id(s)) then
          temp = 0.0_DP
          temp1 = 0.0_DP
          call dgemm('N', 'N', 3*nat_sl, 3*nat_sl, 3*nat_sl, 1.0_DP, &
                dCpqdV(:, :, counter), 3*nat_sl,&
                Cpq, 3*nat_sl, 0.0_DP, temp1, 3*nat_sl)

          call dgemm('T', 'N', 3*nat_sl, 3*nat_sl, 3*nat_sl, 1.0_DP, &
                Cpq, 3*nat_sl,&
                temp1, 3*nat_sl, 0.0_DP, temp, 3*nat_sl)

          do j=1,3*nat_sl,1
            if(eigenvalues(j).ge.0.0_DP) then
              forcedV(s) = forcedV(s) + eval_inv(j)*temp(j,j)
              end if
            end do
          counter = counter + 1
          else
          forcedV(s) = 0.0_DP
          end if
        end do
      call mp_sum(forcedV, intra_image_comm)
      forcedV = forcedV*nat/nat_sl
      end if


    !!!!!!!!!!!!!!!!!!!!
    ! Forces below here. There's going to be some long parallelization business.
    !!!!!!!!!!!!!!!!!!!!

    call start_clock('mbd_int_forces')

    if(do_forces) then
      counter = 1
      forcedR = 0.0_DP
      do s=1,nat,1
        if(me_image.eq.f_cpu_id(s)) then
          do i=1,3,1
            temp = 0.0_DP
            temp1 = 0.0_DP
            call dgemm('N', 'N', 3*nat_sl, 3*nat_sl, 3*nat_sl, 1.0_DP, &
                  dCpqdR(:, :, counter, i), 3*nat_sl,&
                  Cpq, 3*nat_sl, 0.0_DP, temp1, 3*nat_sl)

            call dgemm('T', 'N', 3*nat_sl, 3*nat_sl, 3*nat_sl, 1.0_DP, &
                  Cpq, 3*nat_sl,&
                  temp1, 3*nat_sl, 0.0_DP, temp, 3*nat_sl)

            do j=1,3*nat_sl,1
              if(eigenvalues(j).ge.0.0_DP) then
                forcedR(s, i) = forcedR(s, i) + eval_inv(j)*temp(j,j)
                end if
              end do
            end do
            counter = counter + 1
          end if
        end do
      call mp_sum(forcedR, intra_image_comm)
      forcedR = forcedR*nat/nat_sl

      do s=1,3,1
        if(f_cpu_id(s).eq.me_image) then
          do i=1,3,1
            temp = 0.0_DP
            temp1 = 0.0_DP
            call dgemm('N', 'N', 3*nat_sl, 3*nat_sl, 3*nat_sl, 1.0_DP, &
                  dCpqdh(:, :, s, i), 3*nat_sl,&
                  Cpq, 3*nat_sl, 0.0_DP, temp1, 3*nat_sl)
            call dgemm('T', 'N', 3*nat_sl, 3*nat_sl, 3*nat_sl, 1.0_DP, &
                  Cpq, 3*nat_sl,&
                  temp1, 3*nat_sl, 0.0_DP, temp, 3*nat_sl)
            do j=1,3*nat_sl,1
              if(eigenvalues(j).ge.0.0_DP) then
                  forcedh(s, i) = forcedh(s, i) + 1.0_DP/(2.0_DP*dsqrt(eigenvalues(j)))*temp(j,j)
                end if
              end do
            end do
          end if
        end do
      call mp_sum(forcedh, intra_image_comm)
      forcedh = forcedh*nat/nat_sl
    end if

    call stop_clock('mbd_int_forces')
    call stop_clock('mbd_int_energy')
    return
  end subroutine mbdvdw_interacting_energy

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This method simply performs the trace of the A matrix found in the SCS
  ! procedure. The equation is
  ! $\bar{\alpha}_p(i \omega) = \frac{1}{3} Tr \left[ \sum_q \bar{A}_{pq}(i \omega) \right]$
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mbdvdw_calculate_screened_pol(alpha_isotropic, dalpha_isotropicdR, dalpha_isotropicdh, dalpha_isotropicdV)
    implicit none

    ! variables to be returned
    real(dp), dimension(nat), intent(out) :: alpha_isotropic
    real(dp), dimension(nat, nat, 3), intent(out) :: dalpha_isotropicdR
    real(dp), dimension(nat, 3  , 3), intent(out) :: dalpha_isotropicdh
    real(dp), dimension(nat, nat), intent(out) :: dalpha_isotropicdV
    ! local variables
    integer :: p_atom, q_atom, i_index, j_index, s, i, i_idx, counter
    call start_clock('mbd_calculate_screened_pol')

    alpha_isotropic = 0.0_DP
    if(vdw_self_consistent) dalpha_isotropicdV = 0.0_DP
    if(do_forces) dalpha_isotropicdR = 0.0_DP
    if(do_forces) dalpha_isotropicdh = 0.0_DP

    ! This set of loops goes over all atoms, and traces along the rows. We'll create a 3x3 matrix that takes
    ! care to preserve the directional information. Then, we'll compute the eigenvalues of the 3x3
    ! matrix, sum them, and throw them into the isotropic alpha.
    do p_atom=1, nat, 1
      do q_atom=1, nat, 1
        do i=1,3,1
          counter = 1
          do s=1, nat, 1
            if(f_cpu_id(s).eq.me_image) then
              do i_idx=1,3,1
                ! Helper variables defined to help make the code more readable
                i_index = (3*p_atom - 3)
                j_index = (3*q_atom - 3)
                ! Simply takes the trace of the directionally dependent alphas
                if(do_forces) then
                  dalpha_isotropicdR(p_atom, s, i)=dalpha_isotropicdR(p_atom, s, i) &
                                    +dA_matrixdR(i_index+i_idx, j_index+i_idx, counter, i)
                  if(s.le.3) then
                    dalpha_isotropicdh(p_atom, s, i) = dalpha_isotropicdh(p_atom, s, i) &
                                  +dA_matrixdh(i_index+i_idx, j_index+i_idx, s, i)
                    end if
                  end if
                if(vdw_self_consistent.and.(i.eq.1)) then
                  dalpha_isotropicdV(p_atom, s) = dalpha_isotropicdV(p_atom, s) &
                                    +dA_matrixdV(i_index+i_idx, j_index+i_idx, counter)
                  end if
                end do
              counter = counter + 1
              end if
            end do
          end do
          do i_idx=1,3,1
            i_index = (3*p_atom - 3)
            j_index = (3*q_atom - 3)
            alpha_isotropic(p_atom)=alpha_isotropic(p_atom)+A_matrix(i_index+i_idx, j_index+i_idx)
            end do
        end do
      end do

    alpha_isotropic = alpha_isotropic/3.0_DP
    if(do_forces) dalpha_isotropicdR = dalpha_isotropicdR/3.0_DP
    if(do_forces) dalpha_isotropicdh = dalpha_isotropicdh/3.0_DP
    if(vdw_self_consistent) dalpha_isotropicdV = dalpha_isotropicdV/3.0_DP

    call stop_clock('mbd_calculate_screened_pol')
    return
  end subroutine mbdvdw_calculate_screened_pol

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This method is going to compute TSR and its relevant derivatives.
  ! The formula in latex for TSR is:
  ! T^{ij}_{SR} = \left[  1  - f^{TS}_{damp} \right] T^{ij}
  ! with
  ! T^{ij} = T_{dip}^{ij} \left[  erf\left[ \frac{R}{\sigma} \right] - \frac{2 R}{\sqrt{\pi} \sigma} e^{- \frac{R^2}{\sigma^2} } \right] + \frac{4}{\sqrt{\pi}} \frac{R^iR^j}{\sigma^3 R^2} e^{-\frac{R^2}{\sigma^2}}
  ! and
  ! T^{ij}_{dip} = \frac{-3 R^i R^j + R^2 \delta_{ij}}{ R^5 }
  !
  ! and \partial T_{SR} is given by
  ! \partial \mathbf{T}^{ij}_{\rm SR} &=& \mathbf{T}^{ij} \; \partial f\left(Z^{\rm TS}\right) + \left[1- f\left(Z^{\rm TS}\right) \right] \partial \mathbf{T}^{ij},
  ! \partial \mathbf{T}^{ij} &=& -3 \left[{\rm erf}\left[\zeta\right]   -\frac{h_(\zeta)}{2\zeta}   \right]  \partial \mathbf{T}_{\rm dip}^{ij}  \\
  ! &&  - \frac{1}{3} \zeta\, h_(\zeta)  \left[   \partial \mathbf{T}_{\rm dip}^{ij} - \frac{\delta_{ij}}{R^4} \partial R \right]    \nonumber \\
  ! && +  \left[ \mathbf{T}_{\rm dip}^{ij} + \frac{R^iR^j}{ R^5} \Big[ 3  - 2\zeta^2  \Big]\right]  h_(\zeta)\,  \partial \zeta,   \nonumber
  !
  ! Yeah... that's super ugly.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mbdvdw_compute_TSR(p, q, i_lat, j_lat, k_lat, RPQ, Spq_lat, TSR, dTSRdR, dTSRdh, dTSRdV)
    implicit none
    integer, intent(in) :: p, q, i_lat, j_lat, k_lat
    ! Variables needed for Tsr
    real(dp), dimension(3, 3), intent(out):: TSR
    real(dp) :: Rpq_norm, Spq, Sigma_pq, R_vdw_pq, Z, fermi_fn, zeta, U, W
    real(dp), dimension(3), intent(in) :: Rpq, Spq_lat
    real(dp), dimension(3,3) :: Tdip, Rmat, TGG

    ! Variables needed for dTsrdR
    real(dp), dimension(3, 3, nat, 3), intent(out) :: dTSRdR
    real(dp) :: dRpq_normdR, dsigma_pqdR, dZetadR, dZdR, dSpqdR, dFermi_fndR, dUdR, dWdR
    real(dp), dimension(3) :: dRpqdR
    real(dp), dimension(3, 3) :: dRmatdR, dTdipdR

    ! Variables needed for dTsrdh
    real(dp), dimension(3, 3, 3, 3), intent(out) :: dTSRdh
    real(dp) :: dRpq_normdH, dsigma_pqdh, dZetadh, dZdh, dSpqdh, dFermi_fndh, dUdh, dWdh
    real(dp), dimension(3) :: dRpqdh
    real(dp), dimension(3, 3) :: dRmatdh, dTdipdh

    ! Variables needed for dTsrdV
    real(dp), dimension(3, 3, nat), intent(out) :: dTSRdV
    real(dp) :: dRpq_normdV, dsigma_pqdV, dZetadV, dZdV, dSpqdV, dFermi_fndV, dUdV, dWdV
    real(dp), dimension(3) :: dRpqdV
    real(dp), dimension(3, 3) :: dRmatdV, dTdipdV

    ! Loop variables
    integer :: i, j, s, i_idx, j_idx

    ! Debug variables
    integer, parameter :: Rpq_check          = 200
    integer, parameter :: dRpq_check_dh      = 201
    integer, parameter :: dRpq_check_dr      = 202
    integer, parameter :: Rpq_norm_check     = 203
    integer, parameter :: dRpq_norm_check_dh = 204
    integer, parameter :: dRpq_norm_check_dr = 205
    integer, parameter :: zeta_check         = 206
    integer, parameter :: dzeta_check_dh     = 207
    integer, parameter :: dzeta_check_dr     = 208
    integer, parameter :: z_check            = 209
    integer, parameter :: dz_check_dh        = 210
    integer, parameter :: dz_check_dr        = 211
    integer, parameter :: f_check            = 212
    integer, parameter :: df_check_dh        = 213
    integer, parameter :: df_check_dr        = 214
    integer, parameter :: Rmat_check         = 215
    integer, parameter :: dRmat_check_dh     = 216
    integer, parameter :: dRmat_check_dr     = 217
    integer, parameter :: Tsr_check          = 218
    integer, parameter :: dTsr_check_dh      = 219
    integer, parameter :: dTsr_check_dr      = 220
    integer, parameter :: Tdip_check         = 221
    integer, parameter :: dTdip_check_dh     = 222
    integer, parameter :: dTdip_check_dr     = 223
    integer, parameter :: dW_check_dh        = 224
    integer, parameter :: dW_check_dr        = 225
    integer, parameter :: W_check            = 226
    integer, parameter :: dU_check_dh        = 227
    integer, parameter :: dU_check_dr        = 228
    integer, parameter :: U_check            = 229
    integer, parameter :: dSigma_pq_check_dh = 230
    integer, parameter :: dSigma_pq_check_dr = 231
    integer, parameter :: sigma_pq_check     = 232
    integer, parameter :: ds_pq_check_dh     = 233
    integer, parameter :: ds_pq_check_dr     = 234
    integer, parameter :: s_pq_check         = 235

    TSR = 0.0_DP
    dTSRdr = 0.0_DP
    dTSRdh = 0.0_DP
    dTSRdV = 0.0_DP
    zeta = 0.0_DP

    ! Computes the cartesian distance from the vector distance
    Rpq_norm = dsqrt(Rpq(1) **2.0_DP + Rpq(2)**2.0_DP + Rpq(3)**2.0_DP)

    ! Computes the effective correlation length of the interaction potential
    ! defined from the widths of the QHO Gaussians
    Sigma_pq = dsqrt(sigma(p)**2.0_DP + sigma(q)**2.0_DP)
    ! Computes the damping radius
    R_VdW_pq = R_TS_VdW(p) + R_TS_VdW(q)
    Spq = beta*R_VdW_pq
    Z = 6.0_DP*(Rpq_norm/Spq - 1.0_DP)
    zeta = Rpq_norm/Sigma_pq
    fermi_fn = 1.0_DP

    ! computes the fermi damping function. The latex for this is
    ! f_{damp}(R_{pq}) = \frac{1}{ 1 + exp( - Z(R_{pq}) ) }
    ! where Z = 6 \left( \frac{R_{pq} }{ S_{pq}} - 1 \right)
    ! and S_{pq} = \beta \left(  R_{p, VdW} + R_{q, VdW} \right)
    if(Z.le.35.0_DP) then
      fermi_fn = 1.0_DP/(1.0_DP+exp( -Z))
      ! Computes the factors for U
      ! U = {\rm erf}\left[\zeta\right] -  \frac{2}{\sqrt{\pi}}\zeta \exp\left[-\zeta^2\right]
      zeta = Rpq_norm/Sigma_pq
      end if

    if(zeta.ge.6.0_DP) then
      U = 1.0_DP
      W = 0.0_DP
      else
      U = derf(zeta)-(2.0_DP*zeta)/dsqrt(pi)*dexp(-zeta*zeta)
      ! Computes the first half of the factor for W before we multiply by the R^i R^j tensor
      ! \mathbf{W}^{ij} &\equiv&   \left(\frac{R^i R^j}{R^5}\right) \, \frac{4}{\sqrt{\pi}}  \zeta^3  \exp\left[-\zeta^2\right]
      W = 4.0_DP*zeta**3.0_DP/dsqrt(pi)*dexp(-zeta*zeta)
      end if

    ! Loops over the cartesian coordinates to compute the R^i*R^j
    ! matrix that goes in to constructing the dipole matrix
    do i=1,3,1
      do j=1,3,1
        Rmat(i, j) = Rpq(i)*Rpq(j)/(Rpq_norm**5.0_DP)
        Tdip(i, j) = -3.0_DP*Rmat(i, j)
        end do
      ! This just applies the kronheker delta center for the dipole tensor. Recall
      ! that T^{ij}_{dip} = \frac{-3 R^i R^j + R^2 \delta_{ij}}{ R^5 }
      Tdip(i, i) = Tdip(i, i) + 1.0_DP/(Rpq_norm**3.0_DP)
      end do

    ! Computes the short range dipole coupling quantity using the fact that
    ! T = T_{dip}\left[ U \right] + W
    TGG = (Tdip*U + Rmat*W)
    TSR = (1.0_DP - fermi_fn)*TGG

    ! loops over all the nuclear derivatives
    ! The first loop loops over all the nuclei
    if(do_forces.or.vdw_self_consistent) then
      do s = 1, nat, 1
        if(vdw_self_consistent) then
          ! Does dV
          dRmatdV = 0.0_DP
          dTdipdV = 0.0_DP
          dRpqdV = 0.0_DP
          dRpq_normdV = 0.0_DP
          dsigma_pqdV = (sigma(p)*dsigmadV(p, s) + sigma(q)*dsigmadV(q, s))/Sigma_pq
          dSpqdV = beta*(dR_TS_VdWdV(p, s) + dR_TS_VdWdV(q, s))
          dZetadV = dRpq_normdV/Sigma_pq - Rpq_norm*dsigma_pqdV/Sigma_pq**2.0_DP
          dZdV = 6.0_DP*( dRpq_normdV/Spq - Rpq_norm*dSpqdV/Spq**2.0_DP)
          dFermi_fndV = exp(-Z)/(1.0_DP + exp(-Z))**2.0_DP*dZdV

          if(zeta.ge.6.0_DP) then
            dUdV = 0.0_DP
            dWdV = 0.0_DP
            else
            dUdV = dZetadV*zeta*zeta*4.0_DP*exp(-zeta*zeta)/dsqrt(pi)
            dWdV = 4.0_DP/dsqrt(pi)*zeta*zeta*exp(-zeta*zeta)*dzetadV*(3.0_DP - 2.0_DP*zeta*zeta)
            end if
          dTSRdV(:, :, s) = -TGG*dFermi_fndV + (1.0_DP - fermi_fn)*(dTdipdV*U+Tdip*dUdV+W*dRmatdV+Rmat*dWdV)
          end if

        if(do_forces) then
          do i=1,3,1
            ! Does dR
            dRmatdR = 0.0_DP
            dTdipdR = 0.0_DP
            dRpqdR = 0.0_DP
            dRpq_normdR = 0.0_DP
            if(q.eq.s) then
              dRpqdR(i) = 1.0_DP
              dRpq_normdR = Rpq(i)/(Rpq_norm)
              else if(p.eq.s) then
              dRpqdR(i) = -1.0_DP
              dRpq_normdR = -Rpq(i)/(Rpq_norm)
              end if
            do i_idx=1,3,1
              do j_idx=1,3,1
                dRmatdR(i_idx, j_idx) = (Rpq(i_idx)*dRpqdR(j_idx) + dRpqdR(i_idx)*Rpq(j_idx))/Rpq_norm**5.0_DP&
                                - 5.0_DP*Rpq(i_idx)*Rpq(j_idx)*dRpq_normdR/Rpq_norm**6.0_DP
                dTdipdR(i_idx, j_idx) = -3.0_DP*dRmatdR(i_idx, j_idx)
                end do
              dTdipdR(i_idx, i_idx) = dTdipdR(i_idx, i_idx) - 3.0_DP*dRpq_normdR/Rpq_norm**4.0_DP
              end do
            dsigma_pqdR = (sigma(p)*dsigmadR(p, s, i) + sigma(q)*dsigmadR(q, s, i))/Sigma_pq
            dZetadR = dRpq_normdR/Sigma_pq - Rpq_norm*dsigma_pqdR/Sigma_pq**2.0_DP
            dSpqdR = beta*(dR_TS_VdWdR(p, s, i) + dR_TS_VdWdR(q, s, i))
            dZdR = 6.0_DP*( dRpq_normdR/Spq - Rpq_norm*dSpqdR/Spq**2.0_DP)
            dFermi_fndR = exp(-Z)/(1.0_DP + exp(-Z))**2.0_DP*dZdR
            if(zeta.ge.6.0_DP) then
              dUdR = 0.0_DP
              dWdR = 0.0_DP
              else
              dUdR = dZetadR*zeta*zeta*4.0_DP*exp(-zeta*zeta)/dsqrt(pi)
              dWdR = 4.0_DP/dsqrt(pi)*zeta*zeta*exp(-zeta*zeta)*dzetadR*(3.0_DP - 2.0_DP*zeta*zeta)
              end if
            dTSRdr(:, :, s, i) = -TGG*dFermi_fndR + (1.0_DP - fermi_fn)*(dTdipdR*U + Tdip*dUdR + W*dRmatdR + Rmat*dWdR)

            if(mbd_debug_dr) then
              open (unit=dRpq_check_dr     , file=trim(dir)//"dRpq_dr.tmp      ", action=act,status=stat,position=pos)
              open (unit=dRpq_norm_check_dr, file=trim(dir)//"dRpq_norm_dr.tmp ", action=act,status=stat,position=pos)
              open (unit=dzeta_check_dr    , file=trim(dir)//"dzeta_dr.tmp     ", action=act,status=stat,position=pos)
              open (unit=dz_check_dr       , file=trim(dir)//"dz_dr.tmp        ", action=act,status=stat,position=pos)
              open (unit=dU_check_dr       , file=trim(dir)//"dU_dr.tmp        ", action=act,status=stat,position=pos)
              open (unit=dW_check_dr       , file=trim(dir)//"dW_dr.tmp        ", action=act,status=stat,position=pos)
              open (unit=dSigma_pq_check_dr, file=trim(dir)//"dsigma_pq_dr.tmp ", action=act,status=stat,position=pos)
              open (unit=ds_pq_check_dr    , file=trim(dir)//"ds_pq_dr.tmp"     , action=act,status=stat,position=pos)
              open (unit=df_check_dr       , file=trim(dir)//"df_dr.tmp        ", action=act,status=stat,position=pos)
              open (unit=dRmat_check_dr    , file=trim(dir)//"dRmat_dr.tmp     ", action=act,status=stat,position=pos)
              open (unit=dTsr_check_dr     , file=trim(dir)//"dTsr_dr.tmp      ", action=act,status=stat,position=pos)
              open (unit=dTdip_check_dr    , file=trim(dir)//"dTdip_dr.tmp     ", action=act,status=stat,position=pos)

              write(dRpq_check_dr     , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dRpqdr
              write(dRpq_norm_check_dr, *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dRpq_normdr
              write(dzeta_check_dr    , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dzetadr
              write(dz_check_dr       , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dZdr
              write(dSigma_pq_check_dr, *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dsigma_pqdr
              write(ds_pq_check_dr,     *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dSpqdr
              write(dU_check_dr       , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dUdr
              write(dW_check_dr       , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dWdr
              write(df_check_dr       , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dFermi_fndr
              write(dRmat_check_dr    , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dRmatdr
              write(dTsr_check_dr     , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dTsrdr(:, :, s, i)
              write(dTdip_check_dr    , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dTdipdr

              close(dRpq_check_dr     )
              close(dRpq_norm_check_dr)
              close(dzeta_check_dr    )
              close(dSigma_pq_check_dr)
              close(ds_pq_check_dr    )
              close(dz_check_dr       )
              close(du_check_dr       )
              close(dw_check_dr       )
              close(df_check_dr       )
              close(dRmat_check_dr    )
              close(dTsr_check_dr     )
              close(dTdip_check_dr    )
              end if

            ! Does dH
            if(s.le.3) then
              dRmatdh = 0.0_DP
              dTdipdh = 0.0_DP
              dRpqdh = 0.0_DP
              dRpq_normdh = 0.0_DP
              dRpq_normdh = -Rpq(i)*Spq_lat(s)/(rpq_norm)
              dRpqdh(i) = -Spq_lat(s)
              do i_idx=1,3,1
                do j_idx=1,3,1
                  dRmatdH(i_idx, j_idx) = (Rpq(i_idx)*dRpqdH(j_idx) + dRpqdH(i_idx)*Rpq(j_idx))/Rpq_norm**5.0_DP&
                                  - 5.0_DP*Rpq(i_idx)*Rpq(j_idx)*dRpq_normdH/Rpq_norm**6.0_DP
                  dTdipdH(i_idx, j_idx) = -3.0_DP*dRmatdH(i_idx, j_idx)
                  end do
                dTdipdH(i_idx, i_idx) = dTdipdH(i_idx, i_idx) - 3.0_DP*dRpq_normdH/Rpq_norm**4.0_DP
                end do
              dsigma_pqdH = (sigma(p)*dsigmadH(p, s, i) + sigma(q)*dsigmadH(q, s, i))/Sigma_pq
              dZetadH = dRpq_normdH/Sigma_pq - Rpq_norm*dsigma_pqdH/Sigma_pq**2.0_DP
              dSpqdH = beta*(dR_ts_VdWdH(p, s, i) + dR_ts_VdWdH(q, s, i))
              dZdH = 6.0_DP*( dRpq_normdH/Spq - Rpq_norm*dSpqdH/Spq**2.0_DP)
              dFermi_fndH = exp(-Z)/(1.0_DP + exp(-Z))**2.0_DP*dZdH
              if(zeta.ge.6.0_DP) then
                dUdH = 0.0_DP
                dWdH = 0.0_DP
                else
                dUdH = dZetadH*zeta*zeta*4.0_DP*exp(-zeta*zeta)/dsqrt(pi)
                dWdH = 4.0_DP/dsqrt(pi)*zeta*zeta*exp(-zeta*zeta)*dzetadH*(3.0_DP - 2.0_DP*zeta*zeta)
                end if
              dTSRdH(:, :, s, i) = -TGG*dFermi_fndH + (1.0_DP - fermi_fn)*(dTdipdH*U + Tdip*dUdH + W*dRmatdH + Rmat*dWdH)

              if(mbd_debug_dh) then
                open (unit=dRpq_check_dh     , file=trim(dir)//"dRpq_dh.tmp      ", action=act,status=stat,position=pos)
                open (unit=dRpq_norm_check_dh, file=trim(dir)//"dRpq_norm_dh.tmp ", action=act,status=stat,position=pos)
                open (unit=dzeta_check_dh    , file=trim(dir)//"dzeta_dh.tmp     ", action=act,status=stat,position=pos)
                open (unit=dz_check_dh       , file=trim(dir)//"dz_dh.tmp        ", action=act,status=stat,position=pos)
                open (unit=dU_check_dh       , file=trim(dir)//"dU_dh.tmp        ", action=act,status=stat,position=pos)
                open (unit=dW_check_dh       , file=trim(dir)//"dW_dh.tmp        ", action=act,status=stat,position=pos)
                open (unit=dSigma_pq_check_dh, file=trim(dir)//"dsigma_pq_dh.tmp ", action=act,status=stat,position=pos)
                open (unit=ds_pq_check_dh    , file=trim(dir)//"ds_pq_dh.tmp"     , action=act,status=stat,position=pos)
                open (unit=df_check_dh       , file=trim(dir)//"df_dh.tmp        ", action=act,status=stat,position=pos)
                open (unit=dRmat_check_dh    , file=trim(dir)//"dRmat_dh.tmp     ", action=act,status=stat,position=pos)
                open (unit=dTsr_check_dh     , file=trim(dir)//"dTsr_dh.tmp      ", action=act,status=stat,position=pos)
                open (unit=dTdip_check_dh    , file=trim(dir)//"dTdip_dh.tmp     ", action=act,status=stat,position=pos)

                write(dRpq_check_dh     , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dRpqdh
                write(dRpq_norm_check_dh, *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dRpq_normdh
                write(dzeta_check_dh    , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dzetadH
                write(dz_check_dh       , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dZdh
                write(dSigma_pq_check_dh, *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dsigma_pqdH
                write(ds_pq_check_dh,     *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dSpqdh
                write(dU_check_dh       , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dUdh
                write(dW_check_dh       , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dWdh
                write(df_check_dh       , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dFermi_fndH
                write(dRmat_check_dh    , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dRmatdh
                write(dTsr_check_dh     , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dTsrdh(:, :, s, i)
                write(dTdip_check_dh    , *) omega_to_print, p, q, i_lat, j_lat, k_lat, s, i, dTdipdh

                close(dRpq_check_dh     )
                close(dRpq_norm_check_dh)
                close(dzeta_check_dh    )
                close(dSigma_pq_check_dh)
                close(ds_pq_check_dh    )
                close(dz_check_dh       )
                close(du_check_dh       )
                close(dw_check_dh       )
                close(df_check_dh       )
                close(dRmat_check_dh    )
                close(dTsr_check_dh     )
                close(dTdip_check_dh    )
                end if
              end if
            end do
          end if
        end do
      else
      end if

      if(mbd_debug_dh.or.mbd_debug_dr) then
        open (unit=Rpq_check      , file=trim(dir)//"Rpq.tmp       ", action=act,status=stat,position=pos)
        open (unit=Rpq_norm_check , file=trim(dir)//"Rpq_norm.tmp  ", action=act,status=stat,position=pos)
        open (unit=zeta_check     , file=trim(dir)//"zeta.tmp      ", action=act,status=stat,position=pos)
        open (unit=sigma_pq_check , file=trim(dir)//"sigma_pq.tmp  ", action=act,status=stat,position=pos)
        open (unit=s_pq_check     , file=trim(dir)//"s_pq.tmp      ", action=act,status=stat,position=pos)
        open (unit=z_check        , file=trim(dir)//"z.tmp         ", action=act,status=stat,position=pos)
        open (unit=w_check        , file=trim(dir)//"w.tmp         ", action=act,status=stat,position=pos)
        open (unit=u_check        , file=trim(dir)//"u.tmp         ", action=act,status=stat,position=pos)
        open (unit=f_check        , file=trim(dir)//"f.tmp         ", action=act,status=stat,position=pos)
        open (unit=Rmat_check     , file=trim(dir)//"Rmat.tmp      ", action=act,status=stat,position=pos)
        open (unit=Tsr_check      , file=trim(dir)//"Tsr.tmp       ", action=act,status=stat,position=pos)
        open (unit=Tdip_check     , file=trim(dir)//"Tdip.tmp      ", action=act,status=stat,position=pos)

        write(Rpq_check     , *) omega_to_print, p, q, i_lat, j_lat, k_lat, Rpq
        write(Rpq_norm_check, *) omega_to_print, p, q, i_lat, j_lat, k_lat, Rpq_norm
        write(zeta_check    , *) omega_to_print, p, q, i_lat, j_lat, k_lat, zeta
        write(sigma_pq_check, *) omega_to_print, p, q, i_lat, j_lat, k_lat, Sigma_pq
        write(s_pq_check    , *) omega_to_print, p, q, i_lat, j_lat, k_lat, Spq
        write(z_check       , *) omega_to_print, p, q, i_lat, j_lat, k_lat, z
        write(w_check       , *) omega_to_print, p, q, i_lat, j_lat, k_lat, w
        write(u_check       , *) omega_to_print, p, q, i_lat, j_lat, k_lat, u
        write(f_check       , *) omega_to_print, p, q, i_lat, j_lat, k_lat, fermi_fn
        write(Rmat_check    , *) omega_to_print, p, q, i_lat, j_lat, k_lat, Rmat
        write(Tsr_check     , *) omega_to_print, p, q, i_lat, j_lat, k_lat, Tsr
        write(Tdip_check    , *) omega_to_print, p, q, i_lat, j_lat, k_lat, Tdip

        close(Rpq_check     )
        close(Rpq_norm_check)
        close(sigma_pq_check)
        close(s_pq_check    )
        close(zeta_check    )
        close(z_check       )
        close(u_check       )
        close(w_check       )
        close(f_check       )
        close(Rmat_check    )
        close(Tsr_check     )
        close(Tdip_check    )
        end if

  end subroutine mbdvdw_compute_TSR

  subroutine mbdvdw_TSR(p, q, TSR, dTSRdR, dTSRdh, dTSRdV)
    implicit none

    ! IO Variables
    integer, intent(in) :: p, q
    real(dp), dimension(3, 3), intent(out):: TSR
    real(dp), dimension(3, 3, nat, 3), intent(out) :: dTSRdR
    real(dp), dimension(3, 3, nat), intent(out) :: dTSRdV
    real(dp), dimension(3, 3, 3, 3), intent(out) :: dTSRdh
    real(dp), dimension(3) :: Rpq
    real(dp), dimension(3) :: Spq_lat, Spq
    real(dp), dimension(3, 3) :: TSR_temp
    real(dp), dimension(3, 3, nat, 3) :: dTSRdR_temp
    real(dp), dimension(3, 3, nat)  :: dTSRdV_temp
    real(dp), dimension(3, 3, 3, 3) :: dTSRdh_temp
    real(dp) :: max_value
    integer :: n_period
    logical :: converged
    integer :: i_lat, j_lat, k_lat
    integer :: i_idx, j_idx, k_idx

    call start_clock('mbd_TSR')

    TSR = 0.0_DP
    dTSRdh = 0.0_DP
    dTSRdR = 0.0_DP
    dTSRdV = 0.0_DP

    Spq = tau_s(:, p) - tau_s(:, q)

    ! Computes the vector distance between two atoms
    n_period = 0
    max_value = 0.0_DP
    converged = .false.
    do while(.not.converged)
      max_value = 0.0_DP
      do i_idx = -n_period, n_period
        do j_idx = -n_period, n_period
          do k_idx = -n_period, n_period
            if((abs(i_idx).eq.n_period).or.(abs(j_idx).eq.n_period).or.(abs(k_idx).eq.n_period)) then
             if(mbdvdw_vacuum(1)) then
              if(i_idx.ne.0) cycle
              i_lat = 0
              else
              i_lat = i_idx
              end if
             if(mbdvdw_vacuum(2)) then
              if(j_idx.ne.0) cycle
              j_lat = 0
              else
              j_lat = j_idx
              end if
             if(mbdvdw_vacuum(3)) then
              if(k_idx.ne.0) cycle
              k_lat = 0
              else
              k_lat = k_idx
              end if

             Spq_lat(1) = Spq(1) + DBLE(i_lat)
             Spq_lat(2) = Spq(2) + DBLE(j_lat)
             Spq_lat(3) = Spq(3) + DBLE(k_lat)

             Rpq(1)=h_(1,1)*Spq_lat(1)+h_(1,2)*Spq_lat(2)+h_(1,3)*Spq_lat(3)   ! r_AB = h_ s_AB (MIC only if n_period == 0)
             Rpq(2)=h_(2,1)*Spq_lat(1)+h_(2,2)*Spq_lat(2)+h_(2,3)*Spq_lat(3)   ! r_AB = h_ s_AB (MIC only if n_period == 0)
             Rpq(3)=h_(3,1)*Spq_lat(1)+h_(3,2)*Spq_lat(2)+h_(3,3)*Spq_lat(3)   ! r_AB = h_ s_AB (MIC only if n_period == 0)

             TSR_temp = 0.0_DP
             dTSRdR_temp = 0.0_DP
             dTSRdh_temp = 0.0_DP

             call mbdvdw_compute_TSR(p, q, i_lat, j_lat, k_lat, Rpq, Spq_lat, TSR_temp, dTSRdR_temp, dTSRdh_temp, dTSRdV_temp)
             if(max_value.le.maxval(abs(TSR_temp))) then
              max_value = abs(maxval(TSR_temp))
             end if
             TSR = TSR + TSR_temp
             if(vdw_self_consistent) dTSRdV = dTSRdV + dTSRdV_temp
             if(do_forces) dTSRdR = dTSRdR + dTSRdR_temp
             if(do_forces) dTSRdh = dTSRdh + dTSRdh_temp
             end if
          end do
        end do
      end do
      if((max_value.le.mbd_vdw_econv_thr).or.(mbdvdw_vacuum(1).and.mbdvdw_vacuum(2).and.mbdvdw_vacuum(3))) then
        converged = .true.
      end if
      n_period = n_period + 1
    end do

    call stop_clock('mbd_TSR')
  end subroutine mbdvdw_TSR

  subroutine mbdvdw_compute_TLR(p, q, i_lat, j_lat, k_lat, Rpq, Spq_lat, TLR, dTLRdR, dTLRdh, dTLRdV)
    implicit none

    ! IO Variables
    integer, intent(in) :: p, q, i_lat, j_lat, k_lat
    real(dp), dimension(3), intent(in) :: Rpq, Spq_lat
    real(dp), dimension(3, 3), intent(out):: TLR

    ! Variables needed for TLR
    real(dp) :: Rpq_norm, Spq
    real(dp) :: R_vdw_pq
    real(dp), dimension(3,3) :: Tdip
    real(dp) :: Z, fermi_fn

    ! Variables needed for dTLRdR
    real(dp), dimension(3, 3, nat, 3), intent(out) :: dTLRdr
    real(dp) :: dRpq_normdR, dZdR, dSpqdR, dFermi_fndR
    real(dp), dimension(3) :: dRpqdR
    real(dp), dimension(3,3) :: dTdipdR, dRmatdR

    ! Variables needed for dTLRdH
    real(dp), dimension(3, 3, 3, 3), intent(out) :: dTLRdh
    real(dp) :: dRpq_normdh, dZdh, dSpqdh, dFermi_fndh
    real(dp), dimension(3) :: dRpqdh
    real(dp), dimension(3,3) :: dTdipdh, dRmatdh

    ! Variables needed for dTLRdH
    real(dp), dimension(3, 3, nat), intent(out) :: dTLRdV
    real(dp) :: dRpq_normdV, dZdV, dSpqdV, dFermi_fndV
    real(dp), dimension(3) :: dRpqdV
    real(dp), dimension(3,3) :: dTdipdV, dRmatdV

    ! Loop variables
    integer :: i, j, s, i_idx, j_idx

    integer, parameter :: r_pq_tlr_check          = 300
    integer, parameter :: dr_pq_tlr_check_dh      = 301
    integer, parameter :: dr_pq_tlr_check_dr      = 301
    integer, parameter :: r_pq_norm_tlr_check     = 302
    integer, parameter :: dr_pq_norm_tlr_check_dh = 303
    integer, parameter :: dr_pq_norm_tlr_check_dr = 303
    integer, parameter :: s_pq_tlr_check          = 304
    integer, parameter :: ds_pq_tlr_check_dh      = 305
    integer, parameter :: ds_pq_tlr_check_dr      = 305
    integer, parameter :: z_tlr_check             = 306
    integer, parameter :: dz_tlr_check_dh         = 307
    integer, parameter :: dz_tlr_check_dr         = 307
    integer, parameter :: f_tlr_check             = 308
    integer, parameter :: df_tlr_check_dh         = 309
    integer, parameter :: df_tlr_check_dr         = 309
    integer, parameter :: tdip_tlr_check          = 310
    integer, parameter :: dtdip_tlr_check_dh      = 311
    integer, parameter :: dtdip_tlr_check_dr      = 311
    integer, parameter :: tlr_check               = 312
    integer, parameter :: dtlr_check_dh           = 313
    integer, parameter :: dtlr_check_dr           = 313

    ! Computes the cartesian distance from the vector distance
    Rpq_norm = dsqrt(Rpq(1)**2.0_DP + Rpq(2)**2.0_DP + Rpq(3)**2.0_DP)

    ! Computes the damping radius: note the use of screened effective vdW radii
    R_VdW_pq = R_MBD_VdW_sl(p) + R_MBD_VdW_sl(q)
    Spq = beta*R_VdW_pq
    Z = 6.0_DP*(Rpq_norm/Spq - 1.0_DP)

    ! computes the fermi damping function. The latex for this is
    ! f_{damp}(R_{pq}) = \frac{1}{ 1 + exp( - Z(R_{pq}) ) }
    ! where Z = 6 \left( \frac{R_{pq} }{ S_{pq}} - 1 \right)
    ! and S_{pq} = \beta \left(  R_{p, VdW} + R_{q, VdW} \right)
    if(Z.le.35.0_DP) then
      fermi_fn = 1.0_DP/(1.0_DP+exp( -Z))
    else
      call start_clock('mbd_big_z')
      fermi_fn = 1.0_DP
      end if

    ! Loops over the cartesian coordinates to compute the R^i*R^j
    ! matrix that goes in to constructing the dipole matrix
    do i=1,3,1
      do j=1,3,1
        Tdip(i, j) = -3.0*Rpq(i)*Rpq(j)/(Rpq_norm**5.0_DP)
        end do
      ! This just applies the kronheker delta center for the dipole tensor. Recall
      ! that T^{ij}_{dip} = \frac{-3 R^i R^j + R^2 \delta_{ij}}{ R^5 }
      Tdip(i, i) = Tdip(i, i) + 1.0_DP/(Rpq_norm**3.0_DP)
      end do

    ! Computes the short range dipole coupling quantity using the fact that
    ! T = T_{dip}\left[ U \right] + W
    TLR = fermi_fn*Tdip

    ! loops over all the nuclear derivatives
    ! The first loop loops over all the nuclei

    if(do_forces.or.vdw_self_consistent) then
      do s=1,nat,1
        if(vdw_self_consistent) then
          dTdipdV = 0.0_DP
          dRpqdV = 0.0_DP
          dRpq_normdV = 0.0_DP
          dRmatdV = 0.0_DP
          dTdipdV = 0.0_DP
          dSpqdV = beta*(dR_MBD_VdWdV_sl(p, s) + dR_MBD_VdWdV_sl(q, s))
          dZdV = 6.0_DP*(dRpq_normdV/Spq - Rpq_norm*dSpqdV/Spq**2.0_DP)
          dFermi_fndV = exp(-Z)/(1.0_DP + exp(-Z))**2.0_DP*dZdV
          dTLRdV(:, :, s) = Tdip*dFermi_fndV + fermi_fn*dTdipdV
          end if
        if(do_forces) then
          do i=1,3,1
            dTdipdR = 0.0_DP
            dRpqdR = 0.0_DP
            dRpq_normdR = 0.0_DP
            if(q.eq.s) then
              dRpqdR(i) = 1.0_DP
              dRpq_normdR = Rpq(i)/(Rpq_norm)
            else if(p.eq.s) then
              dRpqdR(i) = -1.0_DP
              dRpq_normdR = -Rpq(i)/(Rpq_norm)
              end if
            do i_idx=1,3,1
              do j_idx=1,3,1
                dRmatdR(i_idx, j_idx) = (Rpq(i_idx)*dRpqdR(j_idx) + dRpqdR(i_idx)*Rpq(j_idx))/Rpq_norm**5.0_DP
                dRmatdR(i_idx, j_idx) = dRmatdR(i_idx, j_idx)-5.0_DP*Rpq(i_idx)*Rpq(j_idx)*dRpq_normdR/Rpq_norm**6.0_DP
                end do
              dTdipdR(i_idx, i_idx) = dRpq_normdR/Rpq_norm**4.0_DP
              end do
            dTdipdR = -3.0_DP*(dTdipdR + dRmatdR)
            dSpqdR = beta*(dR_MBD_VdWdR_sl(p, s, i) + dR_MBD_VdWdR_sl(q, s, i))
            dZdR = 6.0_DP*(dRpq_normdR/Spq - Rpq_norm*dSpqdR/Spq**2.0_DP)
            dFermi_fndR = exp(-Z)/(1.0_DP + exp(-Z))**2.0_DP*dZdR
            dTLRdr(:, :, s, i) = Tdip*dFermi_fndR + fermi_fn*dTdipdR

            if(mbd_debug_dr) then
              open(unit=dr_pq_tlr_check_dr     ,file=trim(dir)//"drpq_tlr_dr.tmp",     action=act,status=stat,position=pos)
              open(unit=dr_pq_norm_tlr_check_dr,file=trim(dir)//"drpq_norm_tlr_dr.tmp",action=act,status=stat,position=pos)
              open(unit=ds_pq_tlr_check_dr     ,file=trim(dir)//"dspq_tlr_dr.tmp",     action=act,status=stat,position=pos)
              open(unit=dz_tlr_check_dr        ,file=trim(dir)//"dz_tlr_dr.tmp",       action=act,status=stat,position=pos)
              open(unit=df_tlr_check_dr        ,file=trim(dir)//"df_tlr_dr.tmp",       action=act,status=stat,position=pos)
              open(unit=dtdip_tlr_check_dr     ,file=trim(dir)//"dtdip_tlr_dr.tmp",    action=act,status=stat,position=pos)
              open(unit=dtlr_check_dr          ,file=trim(dir)//"dtlr_dr.tmp",         action=act,status=stat,position=pos)

              write(dr_pq_tlr_check_dr     , *) 0.0_DP, p, q, i_lat, j_lat, k_lat, s, i, drpqdr
              write(dr_pq_norm_tlr_check_dr, *) 0.0_DP, p, q, i_lat, j_lat, k_lat, s, i, dRpq_normdr
              write(ds_pq_tlr_check_dr     , *) 0.0_DP, p, q, i_lat, j_lat, k_lat, s, i, dSpqdr
              write(dz_tlr_check_dr        , *) 0.0_DP, p, q, i_lat, j_lat, k_lat, s, i, dZdr
              write(df_tlr_check_dr        , *) 0.0_DP, p, q, i_lat, j_lat, k_lat, s, i, dFermi_fndr
              write(dtdip_tlr_check_dr     , *) 0.0_DP, p, q, i_lat, j_lat, k_lat, s, i, dTdipdr
              write(dtlr_check_dr          , *) 0.0_DP, p, q, i_lat, j_lat, k_lat, s, i, dTlrdr(:, :, s, i)

              close(dr_pq_tlr_check_dr     )
              close(dr_pq_norm_tlr_check_dr)
              close(ds_pq_tlr_check_dr     )
              close(dz_tlr_check_dr        )
              close(df_tlr_check_dr        )
              close(dtdip_tlr_check_dr     )
              close(dtlr_check_dr          )
              end if

            if(s.le.3) then
              dTdipdh = 0.0_DP
              dRpqdh = 0.0_DP
              dRpq_normdh = 0.0_DP
              dRpq_normdh = -Rpq(i)*Spq_lat(s)/(rpq_norm)
              dRpqdh(i) = -Spq_lat(s)
              do i_idx=1,3,1
                do j_idx=1,3,1
                  dRmatdh(i_idx, j_idx) = (Rpq(i_idx)*dRpqdh(j_idx) + dRpqdh(i_idx)*Rpq(j_idx))/Rpq_norm**5.0_DP
                  dRmatdh(i_idx, j_idx) = dRmatdh(i_idx, j_idx)&
                                    -5.0_DP*Rpq(i_idx)*Rpq(j_idx)*dRpq_normdh/Rpq_norm**6.0_DP
                  end do
                dTdipdh(i_idx, i_idx) = dRpq_normdH/Rpq_norm**4.0_DP
                end do
              dTdipdh = -3.0_DP*(dTdipdh + dRmatdh)
              dSpqdh = beta*(dR_MBD_VdWdh_sl(p, s, i) + dR_MBD_VdWdh_sl(q, s, i))
              dZdh = 6.0_DP*(dRpq_normdh/Spq - Rpq_norm*dSpqdh/Spq**2.0_DP)
              dFermi_fndh = exp(-Z)/(1.0_DP + exp(-Z))**2.0_DP*dZdh
              dTLRdh(:, :, s, i) = Tdip*dFermi_fndh + fermi_fn*dTdipdh
              if(mbd_debug_dh) then
                open(unit=dr_pq_tlr_check_dh     ,file=trim(dir)//"drpq_tlr_dh.tmp",     action=act,status=stat,position=pos)
                open(unit=dr_pq_norm_tlr_check_dh,file=trim(dir)//"drpq_norm_tlr_dh.tmp",action=act,status=stat,position=pos)
                open(unit=ds_pq_tlr_check_dh     ,file=trim(dir)//"dspq_tlr_dh.tmp",     action=act,status=stat,position=pos)
                open(unit=dz_tlr_check_dh        ,file=trim(dir)//"dz_tlr_dh.tmp",       action=act,status=stat,position=pos)
                open(unit=df_tlr_check_dh        ,file=trim(dir)//"df_tlr_dh.tmp",       action=act,status=stat,position=pos)
                open(unit=dtdip_tlr_check_dh     ,file=trim(dir)//"dtdip_tlr_dh.tmp",    action=act,status=stat,position=pos)
                open(unit=dtlr_check_dh          ,file=trim(dir)//"dtlr_dh.tmp",         action=act,status=stat,position=pos)

                write(dr_pq_tlr_check_dh     , *) 0.0_DP, p, q, i_lat, j_lat, k_lat, s, i, drpqdh
                write(dr_pq_norm_tlr_check_dh, *) 0.0_DP, p, q, i_lat, j_lat, k_lat, s, i, dRpq_normdh
                write(ds_pq_tlr_check_dh     , *) 0.0_DP, p, q, i_lat, j_lat, k_lat, s, i, dSpqdh
                write(dz_tlr_check_dh        , *) 0.0_DP, p, q, i_lat, j_lat, k_lat, s, i, dZdh
                write(df_tlr_check_dh        , *) 0.0_DP, p, q, i_lat, j_lat, k_lat, s, i, dFermi_fndh
                write(dtdip_tlr_check_dh     , *) 0.0_DP, p, q, i_lat, j_lat, k_lat, s, i, dTdipdh
                write(dtlr_check_dh          , *) 0.0_DP, p, q, i_lat, j_lat, k_lat, s, i, dTlrdh(:, :, s, i)

                close(dr_pq_tlr_check_dh     )
                close(dr_pq_norm_tlr_check_dh)
                close(ds_pq_tlr_check_dh     )
                close(dz_tlr_check_dh        )
                close(df_tlr_check_dh        )
                close(dtdip_tlr_check_dh     )
                close(dtlr_check_dh          )
                end if
              end if
            end do
          end if
        end do
      end if

    if(mbd_debug_dh.or.mbd_debug_dr) then
      open(unit=r_pq_tlr_check,      file=trim(dir)//"rpq_tlr.tmp",      action=act,status=stat,position=pos)
      open(unit=r_pq_norm_tlr_check, file=trim(dir)//"rpq_norm_tlr.tmp", action=act,status=stat,position=pos)
      open(unit=s_pq_tlr_check,      file=trim(dir)//"spq_tlr.tmp",      action=act,status=stat,position=pos)
      open(unit=z_tlr_check,         file=trim(dir)//"z_tlr.tmp",        action=act,status=stat,position=pos)
      open(unit=f_tlr_check,         file=trim(dir)//"f_tlr.tmp",        action=act,status=stat,position=pos)
      open(unit=tdip_tlr_check,      file=trim(dir)//"tdip_tlr.tmp",     action=act,status=stat,position=pos)
      open(unit=tlr_check,           file=trim(dir)//"tlr.tmp",          action=act,status=stat,position=pos)

      write(r_pq_tlr_check,      *) 0.0_DP, p, q, i_lat, j_lat, k_lat, rpq
      write(r_pq_norm_tlr_check, *) 0.0_DP, p, q, i_lat, j_lat, k_lat, Rpq_norm
      write(s_pq_tlr_check,      *) 0.0_DP, p, q, i_lat, j_lat, k_lat, Spq
      write(z_tlr_check,         *) 0.0_DP, p, q, i_lat, j_lat, k_lat, Z
      write(f_tlr_check,         *) 0.0_DP, p, q, i_lat, j_lat, k_lat, fermi_fn
      write(tdip_tlr_check,      *) 0.0_DP, p, q, i_lat, j_lat, k_lat, Tdip
      write(tlr_check,           *) 0.0_DP, p, q, i_lat, j_lat, k_lat, Tlr

      close(r_pq_tlr_check     )
      close(r_pq_norm_tlr_check)
      close(s_pq_tlr_check     )
      close(z_tlr_check        )
      close(f_tlr_check        )
      close(tdip_tlr_check     )
      close(tlr_check          )
      end if

    if(Z.ge.35.0_DP) then
      call stop_clock('mbd_big_z')
    end if

  end subroutine mbdvdw_compute_TLR

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! And the formula in latex for TLR is
  ! T^{ij}_{LR} = \bar{f}_{damp} T^{ij}_{dip}
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mbdvdw_TLR(p, q, TLR, dTLRdr, dTLRdh, dTLRdV)

    ! IO Variables
    integer, intent(in) :: p, q
    real(dp), dimension(3, 3), intent(out):: TLR
    real(dp), dimension(3, 3, nat, 3), intent(out) :: dTLRdR
    real(dp), dimension(3, 3, 3, 3), intent(out) :: dTLRdh
    real(dp), dimension(3, 3, nat), intent(out) :: dTLRdV
    real(dp), dimension(3) :: Rpq, Spq_lat, Spq
    real(dp), dimension(3,3) :: TLR_temp
    real(dp), dimension(3,3,nat,3) :: dTLRdR_temp
    real(dp), dimension(3,3,nat)   :: dTLRdV_temp
    real(dp), dimension(3,3,3 ,3)  :: dTLRdh_temp
    real(dp) :: max_value
    integer :: n_period
    logical :: converged
    integer :: i_lat, j_lat, k_lat
    integer :: i_idx, j_idx, k_idx

    call start_clock('mbd_TLR')
    Spq = tau_sl_s(:, p) - tau_sl_s(:, q)

    ! Computes the vector distance between two atoms
    n_period = 0
    max_value = 0.0_DP
    converged = .false.
    TLR = 0.0_DP
    dTLRdH = 0.0_DP
    dTLRdR = 0.0_DP
    dTLRdV = 0.0_DP
    do while(.not.converged)
      max_value = 0.0_DP
      do i_idx = -n_period, n_period
        do j_idx = -n_period, n_period
          do k_idx = -n_period, n_period
            if((abs(i_idx).eq.n_period).or.(abs(j_idx).eq.n_period).or.(abs(k_idx).eq.n_period)) then
             if(mbdvdw_vacuum(1)) then
              if(i_idx.ne.0) cycle
              i_lat = 0
              else
              i_lat = i_idx
              end if
             if(mbdvdw_vacuum(2)) then
              if(j_idx.ne.0) cycle
              j_lat = 0
              else
              j_lat = j_idx
              end if
             if(mbdvdw_vacuum(3)) then
              if(k_idx.ne.0) cycle
              k_lat = 0
             else
              k_lat = k_idx
             end if

             Spq_lat(1) = Spq(1) + DBLE(i_lat)
             Spq_lat(2) = Spq(2) + DBLE(j_lat)
             Spq_lat(3) = Spq(3) + DBLE(k_lat)

             Rpq(1)=h_(1,1)*Spq_lat(1)+h_(1,2)*Spq_lat(2)+h_(1,3)*Spq_lat(3)   ! r_AB = h_ s_AB (MIC only if n_period == 0)
             Rpq(2)=h_(2,1)*Spq_lat(1)+h_(2,2)*Spq_lat(2)+h_(2,3)*Spq_lat(3)   ! r_AB = h_ s_AB (MIC only if n_period == 0)
             Rpq(3)=h_(3,1)*Spq_lat(1)+h_(3,2)*Spq_lat(2)+h_(3,3)*Spq_lat(3)   ! r_AB = h_ s_AB (MIC only if n_period == 0)

             TLR_temp = 0.0_DP
             dTLRdR_temp = 0.0_DP
             dTLRdh_temp = 0.0_DP
             dTLRdV_temp = 0.0_DP

             call mbdvdw_compute_TLR(p, q, i_lat, j_lat, k_lat, Rpq, Spq_lat, TLR_temp, dTLRdR_temp, dTLRdh_temp, dTLRdv_temp)
             if(max_value.le.abs(maxval(TLR_temp))) then
              max_value = abs(maxval(TLR_temp))
             end if
             TLR = TLR + TLR_temp
             if(vdw_self_consistent) dTLRdV = dTLRdV + dTLRdV_temp
             if(do_forces) dTLRdR = dTLRdR + dTLRdR_temp
             if(do_forces) dTLRdh = dTLRdh + dTLRdh_temp
             end if
          end do
        end do
      end do
      if((max_value.le.mbd_vdw_econv_thr).or.(mbdvdw_vacuum(1).and.mbdvdw_vacuum(2).and.mbdvdw_vacuum(3))) then
        converged = .true.
      end if
      n_period = n_period + 1
    end do
    call stop_clock('mbd_TLR')
  end subroutine mbdvdw_TLR

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This method solves the self consistent screen equation to solve for \bar{A}
  ! It does by solving the equation
  ! \bar{A}(i\omega) = \left(  A^{-1} + TSR \right)^{-1}
  ! Because A is diagonal, the inverse is trivial, and will be dealt with in the inner loops.
  ! TSR will be generated for each pair of atoms
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mbdvdw_SCS()
  implicit none

  ! Local variables
  real(dp), dimension(3, 3) :: TSR
  real(dp), dimension(3*nat, 3*nat) :: temp, temp1
  real(dp), dimension(:,:,:,:), allocatable:: dTSRdR
  real(dp), dimension(:,:,:), allocatable :: dTSRdv
  real(dp), dimension(3,3,3,3) :: dTSRdh
  ! type(dtdr_lut), dimension(:), allocatable :: my_dtsrdr, collected_dtsrdr, temp_dtsrdr
  ! type(dtdh_lut), dimension(:), allocatable :: my_dtsrdh, collected_dtsrdh, temp_dtsrdh
  real(dp), dimension(:,:,:), allocatable :: my_tsr, collected_tsr, temp_tsr
  real(dp), dimension(:,:,:,:,:), allocatable :: my_dtsrdr, collected_dtsrdr, temp_dtsrdr
  real(dp), dimension(:,:,:,:,:), allocatable :: my_dtsrdh, collected_dtsrdh, temp_dtsrdh
  real(dp), dimension(:,:,:,:),   allocatable :: my_dtsrdV, collected_dtsrdV, temp_dtsrdV

  ! Looping varlables
  integer :: p, q, i_idx, j_idx, i_index, j_index, i, s

  ! lapack variables
  integer :: errorflag
  integer, dimension(3*nat) :: IPIV
  real(dp), dimension(3*nat) :: WORK

  !MPI Vars
  integer :: ierr
  integer, dimension(MPI_STATUS_SIZE) :: rstatus
  integer, dimension(:), allocatable :: offsets, force_offsets

  ! Looping varlables
  integer :: counter, counter_amat, counter_f, counter_p
  integer :: cpuid, j, my_num_pairs
  real(dp) :: divisor

  if(.not.allocated(dTSRdr)) allocate(dTSRdr(3, 3, nat, 3))
  if(.not.allocated(dTSRdV)) allocate(dTSRdV(3, 3, nat))
  call start_clock('mbd_scs')

  call start_clock('mbd_pinit')
  call mbdvdw_para_init()
  call stop_clock('mbd_pinit')

  my_num_pairs = n_pairs(me_image+1)
  if(.not.allocated(offsets)) allocate(offsets(nproc_image))
  if(.not.allocated(force_offsets)) allocate(force_offsets(nproc_image))
  offsets(1) = 0
  force_offsets(1) = 0
  do counter = 2, nproc_image, 1
    offsets(counter) = offsets(counter-1) + n_pairs(counter-1)
    force_offsets(counter) = force_offsets(counter - 1) + n_comps(counter - 1)
    end do

  if(.not.allocated(my_tsr))           allocate(my_tsr(my_num_pairs, 3, 3)); my_tsr = 0.0_DP
  if(.not.allocated(collected_tsr))    allocate(collected_tsr(num_pairs, 3, 3)); collected_tsr = 0.0_DP

  if(vdw_self_consistent) then
    if(.not.allocated(my_dtsrdV))           allocate(my_dtsrdV(my_num_pairs, 3, 3, nat)); my_dtsrdV = 0.0_DP
    if(.not.allocated(collected_dtsrdV))    allocate(collected_dtsrdV(num_pairs, 3, 3, n_comps(me_image+1)));
    collected_dtsrdV = 0.0_DP
    end if

  if(do_forces) then
    if(.not.allocated(my_dtsrdh))           allocate(my_dtsrdh(my_num_pairs, 3, 3, 3, 3)); my_dtsrdh = 0.0_DP
    if(.not.allocated(collected_dtsrdh))    allocate(collected_dtsrdh(num_pairs, 3, 3, 3, 3)); collected_dtsrdh = 0.0_DP

    if(.not.allocated(my_dtsrdr))           allocate(my_dtsrdr(my_num_pairs, 3, 3, nat, 3)); my_dtsrdr = 0.0_DP
    if(.not.allocated(collected_dtsrdr))    allocate(collected_dtsrdr(num_pairs, 3, 3, n_comps(me_image+1), 3));
    collected_dtsrdr = 0.0_DP
    end if

  call start_clock('mbd_lut_tsr')
  counter = 1
  do i = 1, num_pairs, 1
    p = pairs_scs(i)%p
    q = pairs_scs(i)%q
    cpuid = pairs_scs(i)%cpu

    if(cpuid.eq.me_image) then
      tsr = 0.0_DP
      dtsrdr = 0.0_DP
      dtsrdh = 0.0_DP
      dtsrdv = 0.0_DP
      call mbdvdw_tsr(p, q, tsr, dtsrdr, dtsrdh, dtsrdv)
      my_tsr(counter,:,:) = tsr

      if(vdw_self_consistent) my_dtsrdv(counter,:,:,:) = dtsrdv

      if(do_forces) then
        my_dtsrdh(counter,:,:,:,:) = dtsrdh
        my_dtsrdr(counter,:,:,:,:) = dtsrdr
        end if
      counter = counter + 1
      end if
    end do
  call stop_clock('mbd_lut_tsr')

  !!!!!!!!!!!!!!!!!!
  ! Send and Receives for tsr
  !!!!!!!!!!!!!!!!!!
  do counter = 1, n_pairs(me_image+1), 1
    collected_tsr(counter + offsets(me_image+1), :, :) = my_tsr(counter, :, :)
    end do

  do i = 0, max_proc_forces-1, 1                      ! Outer loop is the loop for recvs
    do j = 0, max_proc_pairs-1, 1                   ! Loop for recvs
      if((j.eq.me_image).and.(i.ne.me_image)) then
        call mpi_send(my_tsr,&                  ! Send buffer
                my_num_pairs*3*3,&        ! Count of communicated data
                mpi_double_precision,&    ! Data type
                i,&                       ! Destination
                j,&                       ! Tag (who it was sent from)
                intra_image_comm,&        ! Communicator to send with
                ierr)                     ! Error flag
      else
        if((i.eq.me_image).and.(j.ne.me_image)) then
          if(.not.allocated(temp_tsr))    allocate(temp_tsr(n_pairs(j+1), 3, 3)); temp_tsr = 0.0_DP
          call mpi_recv(temp_tsr,&                    ! Recv buffer
                  n_pairs(j+1)*3*3,&            ! Count of communicated data
                  mpi_double_precision,&        ! Data type
                  j,&                           ! Source (who it was sent from)
                  j,&                           ! Tag (also who it was sent from)
                  intra_image_comm,&            ! Communicator to send with
                  rstatus,&                     ! Status flag
                  ierr)                         ! Error flag
          do counter = 1, n_pairs(j+1), 1
            collected_tsr(counter + offsets(j+1), :, :) = temp_tsr(counter, :, :)
            end do
          if(allocated(temp_tsr)) deallocate(temp_tsr)
          end if
        end if
      end do
    end do

  if(allocated(my_tsr)) deallocate(my_tsr)

  if(vdw_self_consistent) then
    !!!!!!!!!!!!!!!!!!
    ! Send and Receives for dtsrdV
    !!!!!!!!!!!!!!!!!!
    do counter_p = 1, n_pairs(me_image+1), 1
      do counter_f = 1, n_comps(me_image+1), 1
        collected_dtsrdV(offsets(me_image+1)+counter_p,:,:,counter_f) = &
                      my_dtsrdV(counter_p,:,:,force_offsets(me_image+1)+counter_f)
        end do
      end do

    do i = 0, max_proc_forces-1, 1                      ! Outer loop is the loop for recvs
      do j = 0, max_proc_pairs-1, 1                   ! Loop for recvs
        if((j.eq.me_image).and.(i.ne.me_image)) then
          if(.not.allocated(temp_dtsrdV))    allocate(temp_dtsrdV(n_pairs(j+1), 3, 3, n_comps(i+1)));
          temp_dtsrdV = 0.0_DP
          do counter_f = 1, n_comps(i+1), 1
            temp_dtsrdV(:, :, :, counter_f) = my_dtsrdV(:, :, :, force_offsets(i+1) + counter_f)
            end do
          call mpi_send(temp_dtsrdV,&                               ! Send buffer
                  n_pairs(j+1)*3*3*n_comps(i+1),& ! Count of communicated data
                  mpi_double_precision,&                   ! Data type
                  i,&                                      ! Destination
                  j,&                                      ! Tag (who it was sent from)
                  intra_image_comm,&                       ! Communicator to send with
                  ierr)                                    ! Error flag
          if(allocated(temp_dtsrdV))  deallocate(temp_dtsrdV)
        else
          if((i.eq.me_image).and.(j.ne.me_image)) then
            if(.not.allocated(temp_dtsrdV))    allocate(temp_dtsrdV(n_pairs(j+1), 3, 3, n_comps(i+1)));
            call mpi_recv(temp_dtsrdV,&                              ! Recv buffer
                    n_pairs(j+1)*3*3*n_comps(i+1),& ! Count of communicated data
                    mpi_double_precision,&                  ! Data type
                    j,&                                     ! Source (who it was sent from)
                    j,&                                     ! Tag (also who it was sent from)
                    intra_image_comm,&                      ! Communicator to send with
                    rstatus,&                               ! Status flag
                    ierr)                                   ! Error flag
            do counter = 1, n_pairs(j+1), 1
              collected_dtsrdV(counter+offsets(j+1),:,:,:) = temp_dtsrdV(counter,:,:,:)
              end do
            if(allocated(temp_dtsrdV))  deallocate(temp_dtsrdV)
            end if
          end if
        end do
      end do
    if(allocated(temp_dtsrdv)) deallocate(temp_dtsrdV)
    if(allocated(my_dtsrdv)) deallocate(my_dtsrdv)
    end if


  if(do_forces) then
    !!!!!!!!!!!!!!!!!!
    ! Send and Receives for dtsrdr
    !!!!!!!!!!!!!!!!!!
    do counter_p = 1, n_pairs(me_image+1), 1
      do counter_f = 1, n_comps(me_image+1), 1
        collected_dtsrdr(offsets(me_image+1)+counter_p,:,:,counter_f,:) = &
                      my_dtsrdr(counter_p,:,:,force_offsets(me_image+1)+counter_f,:)
        end do
      end do

    do i = 0, max_proc_forces-1, 1                      ! Outer loop is the loop for recvs
      do j = 0, max_proc_pairs-1, 1                   ! Loop for recvs
        if((j.eq.me_image).and.(i.ne.me_image)) then
          if(.not.allocated(temp_dtsrdr))    allocate(temp_dtsrdr(n_pairs(j+1), 3, 3, n_comps(i+1), 3));
          temp_dtsrdr = 0.0_DP
          do counter_f = 1, n_comps(i+1), 1
            temp_dtsrdr(:, :, :, counter_f, :) = my_dtsrdr(:, :, :, force_offsets(i+1) + counter_f, :)
            end do
          call mpi_send(temp_dtsrdr,&                               ! Send buffer
                  n_pairs(j+1)*3*3*n_comps(i+1)*3,& ! Count of communicated data
                  mpi_double_precision,&                   ! Data type
                  i,&                                      ! Destination
                  j,&                                      ! Tag (who it was sent from)
                  intra_image_comm,&                       ! Communicator to send with
                  ierr)                                    ! Error flag
          if(allocated(temp_dtsrdr))  deallocate(temp_dtsrdr)
        else
          if((i.eq.me_image).and.(j.ne.me_image)) then
            if(.not.allocated(temp_dtsrdr))    allocate(temp_dtsrdr(n_pairs(j+1), 3, 3, n_comps(i+1), 3));
            call mpi_recv(temp_dtsrdr,&                              ! Recv buffer
                    n_pairs(j+1)*3*3*n_comps(i+1)*3,& ! Count of communicated data
                    mpi_double_precision,&                  ! Data type
                    j,&                                     ! Source (who it was sent from)
                    j,&                                     ! Tag (also who it was sent from)
                    intra_image_comm,&                      ! Communicator to send with
                    rstatus,&                               ! Status flag
                    ierr)                                   ! Error flag
            do counter = 1, n_pairs(j+1), 1
              collected_dtsrdr(counter+offsets(j+1),:,:,:,:) = temp_dtsrdr(counter,:,:,:,:)
              end do
            if(allocated(temp_dtsrdr))  deallocate(temp_dtsrdr)
            end if
          end if
        end do
      end do

    !!!!!!!!!!!!!!!!!!
    ! Send and Receives for dtsrdh
    !!!!!!!!!!!!!!!!!!
    do counter = 1, n_pairs(me_image+1), 1
      collected_dtsrdh(counter + offsets(me_image+1), :, :, :, :) = my_dtsrdh(counter, :, :, :, :)
      end do

    do i = 0, max_proc_forces-1, 1                      ! Outer loop is the loop for recvs
      do j = 0, max_proc_pairs-1, 1                   ! Loop for recvs
        if((j.eq.me_image).and.(i.ne.me_image)) then
          call mpi_send(my_dtsrdh,&                  ! Send buffer
                  my_num_pairs*3*3*3*3,&        ! Count of communicated data
                  mpi_double_precision,&    ! Data type
                  i,&                       ! Destination
                  j,&                       ! Tag (who it was sent from)
                  intra_image_comm,&        ! Communicator to send with
                  ierr)                     ! Error flag
        else
          if((i.eq.me_image).and.(j.ne.me_image)) then
            if(.not.allocated(temp_dtsrdh))    allocate(temp_dtsrdh(n_pairs(j+1), 3, 3, 3, 3)); temp_dtsrdh=0.0_DP
            call mpi_recv(temp_dtsrdh,&                    ! Recv buffer
                    n_pairs(j+1)*3*3*3*3,&            ! Count of communicated data
                    mpi_double_precision,&        ! Data type
                    j,&                           ! Source (who it was sent from)
                    j,&                           ! Tag (also who it was sent from)
                    intra_image_comm,&            ! Communicator to send with
                    rstatus,&                     ! Status flag
                    ierr)                         ! Error flag
            do counter = 1, n_pairs(j+1), 1
              collected_dtsrdh(counter + offsets(j+1), :, :, :, :) = temp_dtsrdh(counter, :, :, :, :)
              end do
            if(allocated(temp_dtsrdh)) deallocate(temp_dtsrdh)
            end if
          end if
        end do
      end do
    if(allocated(my_dtsrdh)) deallocate(my_dtsrdh)
    if(allocated(dTSRdr)) deallocate(dTSRdr)
    if(allocated(temp_dtsrdr)) deallocate(temp_dtsrdr)
    if(allocated(my_dtsrdr)) deallocate(my_dtsrdr)
    end if

  call start_clock('mbd_a_mat')
  counter_amat = 1
  if(.not.allocated(A_matrix))        allocate(A_matrix(3*nat, 3*nat));                                   A_matrix = 0.0_DP
  if(vdw_self_consistent) then
    if(allocated(dTSRdV)) deallocate(dTSRdV)
    if(.not.allocated(dA_matrixdV))     allocate(dA_matrixdV(3*nat, 3*nat, n_comps(me_image+1)));       dA_matrixdV = 0.0_DP
    if(.not.allocated(dTSRdV))          allocate(dTSRdV(3, 3, n_comps(me_image+1)))
    end if

  if(do_forces) then
    if(.not.allocated(dA_matrixdR))     allocate(dA_matrixdR(3*nat, 3*nat, n_comps(me_image+1), 3));    dA_matrixdR = 0.0_DP
    if(.not.allocated(dA_matrixdh))     allocate(dA_matrixdh(3*nat, 3*nat, 3, 3));                      dA_matrixdh = 0.0_DP
    if(.not.allocated(dTSRdr))          allocate(dTSRdr(3, 3, n_comps(me_image+1), 3))
  end if

  counter = 1
  ! only loops over the upper triangle
  do p=1, nat, 1
    do q=p, nat, 1
      if(p.ne.q) then
        TSR = 0.0_DP
        dTSRdr = 0.0_DP
        dTSRdh = 0.0_DP
        dTSRdV = 0.0_DP
        !call mbdvdw_TSR(p, q, TSR, dTSRdR, dTSRdh)
        TSR = collected_tsr(counter, :, :)
        if(vdw_self_consistent) dTSRdv(:,:,:) = collected_dtsrdv(counter,:,:,:)
        if(do_forces) then
          dTSRdr = collected_dtsrdr(counter, :, :, :, :)
          dTSRdh = collected_dtsrdh(counter, :, :, :, :)
        end if
        counter = counter + 1
      end if
      ! Loops over the cartesian coordinates, fills the A matrix on the diagonals
      ! with the polarizability. And the short-range damping densely.
      do i_idx=1, 3, 1
        do j_idx=1, 3, 1
          ! Helper variables defined to help make the code more readable
          i_index = (3*p - 3 + i_idx)
          j_index = (3*q - 3 + j_idx)
          if(p.eq.q) then ! On the diagonal blocks
            if(i_idx.eq.j_idx) then
              ! Start by storing A^-1 in A
              A_matrix(i_index, j_index) = 1.0_DP/alpha_ts(p)
              ! then compute the derivative of A^-1
              if(do_forces.or.vdw_self_consistent) then
                counter_amat = 1
                do s=1, nat, 1
                  if(f_cpu_id(s).eq.me_image) then
                    divisor = 1.0_DP/alpha_ts(p)**2.0_DP
                    if(vdw_self_consistent) then
                      dA_matrixdV(i_index, j_index, counter_amat) = -dalpha_tsdV(p, s)*divisor
                      end if
                    if(do_forces) then
                      do i=1,3,1
                        dA_matrixdR(i_index, j_index, counter_amat, i) = -dalpha_tsdR(p, s, i)*divisor
                        if(s.le.3) dA_matrixdh(i_index, j_index, s, i) = -dalpha_tsdh(p, s, i)*divisor
                        end do
                      end if
                    counter_amat = counter_amat + 1
                    end if
                  end do
                end if
              end if
            end if
          if(p.ne.q) then
            ! for the off-diagonal blocks A=0 and we add TSR
            A_matrix(i_index, j_index) = A_matrix(i_index, j_index) + TSR(i_idx, j_idx)
            if(do_forces.or.vdw_self_consistent) then
              counter_amat = 1
              do s=1, nat, 1
                if(f_cpu_id(s).eq.me_image) then
                    if(vdw_self_consistent) then
                      dA_matrixdV(i_index,j_index,counter_amat)=&
                        dA_matrixdV(i_index,j_index,counter_amat)+dTSRdV(i_idx,j_idx,counter_amat)
                      dA_matrixdV(j_index,i_index,counter_amat)=dA_matrixdV(i_index,j_index,counter_amat)
                      end if
                    if(do_forces) then
                      do i=1,3,1
                        dA_matrixdR(i_index,j_index,counter_amat,i)=&
                        dA_matrixdR(i_index,j_index,counter_amat,i)+dTSRdR(i_idx,j_idx,counter_amat,i)
                        dA_matrixdR(j_index,i_index,counter_amat,i)=&
                              dA_matrixdR(i_index, j_index,counter_amat,i) !fill in symmetric elements
                        if(s.le.3) then
                          dA_matrixdh(i_index,j_index,s,i)=dA_matrixdh(i_index,j_index,s,i)&
                                              +dTSRdh(i_idx,j_idx,s,i)
                          dA_matrixdh(j_index,i_index,s,i)=dA_matrixdh(i_index, j_index, s, i)
                          end if
                        end do
                      end if
                  counter_amat = counter_amat + 1
                  end if
                end do
              end if
            end if
          ! Fills in the lower triangle of the matrix
          A_matrix(j_index, i_index) = A_matrix(i_index, j_index)
        end do
      end do
    end do
  end do
  call stop_clock('mbd_a_mat')
  ! Performs the LU Decomposition of the A_matrix in preparation for inverting it
  call DGETRF(3*nat, 3*nat, A_matrix, 3*nat, IPIV, errorflag)
    if(errorflag.ne.0) then
      write(stdout, '(3X, "ERROR IN LU DECOMPOSITION")')
      end if
   call DGETRI(3*nat, A_matrix, 3*nat, IPIV, WORK,3*nat,errorflag )
    if(errorflag.ne.0) then
      write(stdout, '(3X, "ERROR IN INVERSION")')
      end if
  ! At this point the A_matrix is now screened.  We complete the derivative of A by again
  ! applying the formula for the derivative of an inverse: \partial\bar{A} = -\bar{A} [\partial inv(\bar{A}) ] \bar{A}
  call start_clock('mbd_a_force')
  if(do_forces.or.vdw_self_consistent) then
    counter = 1
    do s=1, nat, 1
      if(f_cpu_id(s).eq.me_image) then
        if(vdw_self_consistent) then
          call dgemm('N', 'N', 3*nat, 3*nat, 3*nat, 1.0_DP, &
                dA_matrixdV(:, :, counter), 3*nat,&
                A_matrix, 3*nat, 0.0_DP, temp1, 3*nat)

          call dgemm('N', 'N', 3*nat, 3*nat, 3*nat, 1.0_DP, &
                A_matrix, 3*nat,&
                temp1, 3*nat, 0.0_DP, temp, 3*nat)
          dA_matrixdV(:, :, counter) = -temp
          end if

        if(do_forces) then
          do i=1, 3, 1
            call dgemm('N', 'N', 3*nat, 3*nat, 3*nat, 1.0_DP, &
                  dA_matrixdR(:, :, counter, i), 3*nat,&
                  A_matrix, 3*nat, 0.0_DP, temp1, 3*nat)

            call dgemm('N', 'N', 3*nat, 3*nat, 3*nat, 1.0_DP, &
                  A_matrix, 3*nat,&
                  temp1, 3*nat, 0.0_DP, temp, 3*nat)
            dA_matrixdR(:, :, counter, i) = -temp
            if(s.le.3) then
              call dgemm('N', 'N', 3*nat, 3*nat, 3*nat, 1.0_DP, &
                    dA_matrixdh(:, :, s, i), 3*nat,&
                    A_matrix, 3*nat, 0.0_DP, temp1, 3*nat)

              call dgemm('N', 'N', 3*nat, 3*nat, 3*nat, 1.0_DP, &
                    A_matrix, 3*nat,&
                    temp1, 3*nat, 0.0_DP, temp, 3*nat)
              dA_matrixdh(:, :, s, i) = -temp
              end if
            end do
          end if
        counter = counter + 1
        end if
      end do
    end if
  call stop_clock('mbd_a_force')

  call stop_clock('mbd_scs')
  return
  end subroutine mbdvdw_SCS

  subroutine mbdvdw_construct_sl()
    implicit none
    real(dp) :: alpha, R, phi

    ! Looping varlables
    integer :: p, q, i_idx, j_idx, k_idx, i_atom

    !write(stdout, '(3x, "About to construct the SL")')

    if(.not.mbdvdw_vacuum(1)) then
        sl_i = ceiling(supercell_cutoff / dsqrt(h_(1,1)**2.0_DP + h_(2,1)**2.0_DP + h_(3,1)**2.0_DP ))
      else
        sl_i = 1
      end if

    if(.not.mbdvdw_vacuum(2)) then
        sl_j= ceiling(supercell_cutoff / dsqrt(h_(1,2)**2.0_DP + h_(2,2)**2.0_DP + h_(3,2)**2.0_DP ))
      else
        sl_j = 1
      end if

    if(.not.mbdvdw_vacuum(3)) then
        sl_k = ceiling(supercell_cutoff / dsqrt(h_(1,3)**2.0_DP + h_(2,3)**2.0_DP + h_(3,3)**2.0_DP ))
      else
        sl_k = 1
      end if

    nat_sl = nat *sl_i*sl_j*sl_k

    if(.not.allocated(alpha_0_sl))      allocate(alpha_0_sl(nat_sl));                   alpha_0_sl = 0.0_DP
    if(.not.allocated(omega_scs_sl))    allocate(omega_scs_sl(nat_sl));                 omega_scs_sl = 0.0_DP
    if(.not.allocated(R_MBD_VdW_sl))    allocate(R_MBD_VdW_sl(nat_sl));                 R_MBD_VdW = 0.0_DP

    if(vdw_self_consistent) then
      if(.not.allocated(dalpha_0dV_sl))   allocate(dalpha_0dV_sl(nat_sl, nat));        dalpha_0dV_sl = 0.0_DP
      if(.not.allocated(domegadV_sl))     allocate(domegadV_sl(nat_sl, nat));          domegadV_sl = 0.0_DP
      if(.not.allocated(dR_MBD_VdWdV_sl)) allocate(dR_MBD_VdWdV_sl(nat_sl, nat));      dR_MBD_VdWdV_sl = 0.0_DP
      end if

    if(do_forces) then
      if(.not.allocated(dalpha_0dR_sl))   allocate(dalpha_0dR_sl(nat_sl, nat, 3));        dalpha_0dR_sl = 0.0_DP
      if(.not.allocated(dalpha_0dh_sl))   allocate(dalpha_0dh_sl(nat_sl, 3, 3));          dalpha_0dh_sl = 0.0_DP
      if(.not.allocated(domegadR_sl))     allocate(domegadR_sl(nat_sl, nat, 3));          domegadR_sl = 0.0_DP
      if(.not.allocated(domegadh_sl))     allocate(domegadh_sl(nat_sl, 3, 3));            domegadh_sl = 0.0_DP
      if(.not.allocated(dR_MBD_VdWdR_sl)) allocate(dR_MBD_VdWdR_sl(nat_sl, nat, 3));      dR_MBD_VdWdR_sl = 0.0_DP
      if(.not.allocated(dR_MBD_VdWdh_sl)) allocate(dR_MBD_VdWdh_sl(nat_sl, 3, 3));        dR_MBD_VdWdh_sl = 0.0_DP
    end if
    if(.not.allocated(tau_sl))          allocate(tau_sl(3, nat_sl)); tau_sl = 0.0_DP
    if(.not.allocated(tau_sl_s))          allocate(tau_sl_s(3, nat_sl)); tau_sl_s = 0.0_DP
    h_sl = 0.0_DP
    ainv_sl = 0.0_DP

    h_sl(:, 1) = h_(:, 1)*sl_i
    h_sl(:, 2) = h_(:, 2)*sl_j
    h_sl(:, 3) = h_(:, 3)*sl_k

    ainv_sl(:, 1) = ainv_(:, 1)/sl_i
    ainv_sl(:, 2) = ainv_(:, 2)/sl_j
    ainv_sl(:, 3) = ainv_(:, 3)/sl_k

    sl_mult(1) = sl_i
    sl_mult(2) = sl_j
    sl_mult(3) = sl_k

    q = 1
    do i_idx = 0, sl_i-1, 1
      do j_idx = 0, sl_j-1, 1
        do k_idx = 0, sl_k-1, 1
          do p=1, nat, 1
            phi = R_vdw_free(p)/(alpha_free(p))**(1.0_DP/3.0_DP)

            R_MBD_VdW_sl(q)=phi*alpha_0(p)**(1.0_DP/3.0_DP)
            alpha_0_sl(q) = alpha_0(p)
            omega_scs_sl(q) = omega_scs(p)

            if(vdw_self_consistent) then
              dalpha_0dV_sl(q, :) = dalpha_0dV(p, :)
              domegadV_sl(q, :) = domegadV(p, :)
              dR_MBD_VdWdV_sl(q,:) = phi/(3.0_DP*alpha_0(p)**(2.0_DP/3.0_DP))*dalpha_0dV(p, :)
              end if

            if(do_forces) then
              dalpha_0dR_sl(q, :, :) = dalpha_0dR(p, :, :)
              domegadR_sl(q, :, :) = domegadR(p, :, :)
              dR_MBD_VdWdR_sl(q,:,:) = phi/(3.0_DP*alpha_0(p)**(2.0_DP/3.0_DP))*dalpha_0dR(p, :, :)

              dalpha_0dh_sl(q, :, :) = dalpha_0dh(p, :, :)
              domegadh_sl(q, :, :) = domegadh(p, :, :)
              dR_MBD_VdWdh_sl(q,:,:) = phi/(3.0_DP*alpha_0(p)**(2.0_DP/3.0_DP))*dalpha_0dh(p, :, :)
              end if

            tau_sl(:, q) = tau(:, p) + i_idx*h_(:, 1) + j_idx*h_(:, 2) + k_idx*h_(:, 3)
            q = q + 1
            end do
          end do
        end do
      end do

    ! write(stdout, *), sl_i, sl_j, sl_k, nat_sl

    call mbdvdw_pbc(tau_sl, h_sl, ainv_sl, tau_sl_s)

    !write(stdout, '(3x, "Constructed the SL")')
  end subroutine mbdvdw_construct_sl

  subroutine mbdvdw_construct_hamiltonian()
    implicit none
    real(dp), dimension(3, 3) :: TLR
    real(dp), dimension(:, :, :, :), allocatable :: dTLRdr
    real(dp), dimension(3, 3, 3  , 3) :: dTLRdh
    real(dp), dimension(:, :, :), allocatable    :: dTLRdV
    real(dp), dimension(:,:,:), allocatable :: my_tlr, collected_tlr, temp_tlr
    real(dp), dimension(:,:,:,:,:), allocatable :: my_dtlrdr, collected_dtlrdr, temp_dtlrdr
    real(dp), dimension(:,:,:,:,:), allocatable :: my_dtlrdh, collected_dtlrdh, temp_dtlrdh
    real(dp), dimension(:,:,:,:),   allocatable :: my_dtlrdV, collected_dtlrdV, temp_dtlrdV
    real(dp) :: prefactor, pre1, pre2, der, aSqrt, omMult
    integer, parameter :: check = 60

    !MPI Vars
    integer :: ierr
    integer, dimension(MPI_STATUS_SIZE) :: rstatus
    integer, dimension(:), allocatable :: offsets, force_offsets

    ! Looping varlables
    integer :: p, q, i_idx, j_idx, i_index, j_index, i, s, counter
    integer :: cpuid, j, my_num_pairs, counter_c, counter_f, counter_p

    ! Supercell variables
    !write(stdout, '(3x, "building hamiltonian")')
    call start_clock('mbd_construct_hamiltonian')
    call mbdvdw_construct_sl()

    call start_clock('mbd_pinit_sl')
    call mbdvdw_para_init_sl()
    call stop_clock('mbd_pinit_sl')

    call start_clock('mbd_lut_tlr')
    my_num_pairs = n_pairs(me_image+1)
    if(.not.allocated(dTLRdr)) allocate(dTLRdr(3, 3, nat, 3))
    if(.not.allocated(dTLRdV)) allocate(dTLRdV(3, 3, nat))
    if(.not.allocated(offsets)) allocate(offsets(nproc_image))
    if(.not.allocated(force_offsets)) allocate(force_offsets(nproc_image))
    offsets(1) = 0
    force_offsets(1) = 0
    do counter = 2, nproc_image, 1
      offsets(counter) = offsets(counter-1) + n_pairs(counter-1)
      force_offsets(counter) = force_offsets(counter - 1) + n_comps(counter - 1)
      end do

    if(.not.allocated(my_tlr))           allocate(my_tlr(n_pairs(me_image+1), 3, 3)); my_tlr = 0.0_DP
    if(.not.allocated(collected_tlr))    allocate(collected_tlr(num_pairs, 3, 3)); collected_tlr = 0.0_DP

    if(do_forces) then
      if(.not.allocated(my_dtlrdh))           allocate(my_dtlrdh(n_pairs(me_image+1), 3, 3, 3, 3)); my_dtlrdh = 0.0_DP
      if(.not.allocated(collected_dtlrdh))    allocate(collected_dtlrdh(num_pairs, 3, 3, 3, 3)); collected_dtlrdh = 0.0_DP

      if(.not.allocated(my_dtlrdr))           allocate(my_dtlrdr(n_pairs(me_image+1), 3, 3, nat, 3)); my_dtlrdr = 0.0_DP
      if(.not.allocated(collected_dtlrdr))    allocate(collected_dtlrdr(num_pairs, 3, 3, n_comps(me_image+1), 3));
      collected_dtlrdr = 0.0_DP
      end if

    if(vdw_self_consistent) then
      if(.not.allocated(my_dtlrdv))           allocate(my_dtlrdv(n_pairs(me_image+1), 3, 3, nat)); my_dtlrdv = 0.0_DP
      if(.not.allocated(collected_dtlrdv))    allocate(collected_dtlrdv(num_pairs, 3, 3, n_comps(me_image+1)));
      collected_dtlrdv = 0.0_DP
      end if

    counter = 1
    do i = 1, num_pairs, 1
      p = unique_pairs(i)%p
      q = unique_pairs(i)%q
      cpuid = unique_pairs(i)%cpu

      if(cpuid.eq.me_image) then
        TLR = 0.0_DP
        if(do_forces) then
          dTLRdr = 0.0_DP
          dTLRdh = 0.0_DP
          end if
        call mbdvdw_TLR(p, q, TLR, dTLRdr, dTLRdh, dTLRdV)

        my_tlr(counter, :, :) = tlr
        if(vdw_self_consistent) my_dtlrdV(counter, :, :, :) = dTlrdV(:,:,:)
        if(do_forces) my_dtlrdh(counter, :, :, :, :) = dtlrdh(:,:,:,:)
        if(do_forces) my_dtlrdr(counter, :, :, :, :) = dTlrdr(:,:,:,:)
        counter = counter + 1
        end if
      end do
    call stop_clock('mbd_lut_tlr')

    call start_clock('mbd_cpq_comm')

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Send and Receives for tlr
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do counter = 1, n_pairs(me_image+1), 1
      collected_tlr(counter + offsets(me_image+1), :, :) = my_tlr(counter, :, :)
      end do

    do i = 0, max_proc_forces-1, 1                      ! Outer loop is the loop for recvs
      do j = 0, max_proc_pairs-1, 1                   ! Loop for recvs
        if((j.eq.me_image).and.(i.ne.me_image)) then
          call mpi_send(my_tlr,&                  ! Send buffer
                  my_num_pairs*3*3,&        ! Count of communicated data
                  mpi_double_precision,&    ! Data type
                  i,&                       ! Destination
                  j,&                       ! Tag (who it was sent from)
                  intra_image_comm,&        ! Communicator to send with
                  ierr)                     ! Error flag
        else
          if((i.eq.me_image).and.(j.ne.me_image)) then
            if(.not.allocated(temp_tlr))    allocate(temp_tlr(n_pairs(j+1), 3, 3)); temp_tlr = 0.0_DP
            call mpi_recv(temp_tlr,&                    ! Recv buffer
                    n_pairs(j+1)*3*3,&            ! Count of communicated data
                    mpi_double_precision,&        ! Data type
                    j,&                           ! Source (who it was sent from)
                    j,&                           ! Tag (also who it was sent from)
                    intra_image_comm,&            ! Communicator to send with
                    rstatus,&                     ! Status flag
                    ierr)                         ! Error flag
            do counter = 1, n_pairs(j+1), 1
              collected_tlr(counter + offsets(j+1), :, :) = temp_tlr(counter, :, :)
              end do
            if(allocated(temp_tlr)) deallocate(temp_tlr)
            end if
          end if
        end do
      end do
    if(allocated(my_tlr)) deallocate(my_tlr)

  !!!!!!!!!!!!!!!!!!
  ! Send and Receives for dtlrdV
  !!!!!!!!!!!!!!!!!!
  if(vdw_self_consistent) then
    do counter_p = 1, n_pairs(me_image+1), 1
      do counter_f = 1, n_comps(me_image+1), 1
        collected_dtlrdV(offsets(me_image+1)+counter_p,:,:,counter_f) = &
                      my_dtlrdV(counter_p,:,:,force_offsets(me_image+1)+counter_f)
        end do
      end do

    do i = 0, max_proc_forces-1, 1                      ! Outer loop is the loop for recvs
      do j = 0, max_proc_pairs-1, 1                   ! Loop for recvs
        if((j.eq.me_image).and.(i.ne.me_image)) then
          if(.not.allocated(temp_dtlrdV))    allocate(temp_dtlrdV(n_pairs(j+1), 3, 3, n_comps(i+1)));
          temp_dtlrdV = 0.0_DP
          do counter_f = 1, n_comps(i+1), 1
            temp_dtlrdV(:, :, :, counter_f) = my_dtlrdV(:, :, :, force_offsets(i+1) + counter_f)
            end do
          call mpi_send(temp_dtlrdV,&                               ! Send buffer
                  n_pairs(j+1)*3*3*n_comps(i+1),& ! Count of communicated data
                  mpi_double_precision,&                   ! Data type
                  i,&                                      ! Destination
                  j,&                                      ! Tag (who it was sent from)
                  intra_image_comm,&                       ! Communicator to send with
                  ierr)                                    ! Error flag
          if(allocated(temp_dtlrdV))  deallocate(temp_dtlrdV)
        else
          if((i.eq.me_image).and.(j.ne.me_image)) then
            if(.not.allocated(temp_dtlrdV))    allocate(temp_dtlrdV(n_pairs(j+1), 3, 3, n_comps(i+1)));
            call mpi_recv(temp_dtlrdV,&                              ! Recv buffer
                    n_pairs(j+1)*3*3*n_comps(i+1),& ! Count of communicated data
                    mpi_double_precision,&                  ! Data type
                    j,&                                     ! Source (who it was sent from)
                    j,&                                     ! Tag (also who it was sent from)
                    intra_image_comm,&                      ! Communicator to send with
                    rstatus,&                               ! Status flag
                    ierr)                                   ! Error flag
            do counter = 1, n_pairs(j+1), 1
              collected_dtlrdV(counter+offsets(j+1),:,:,:) = temp_dtlrdV(counter,:,:,:)
              end do
            if(allocated(temp_dtlrdV))  deallocate(temp_dtlrdV)
            end if
          end if
        end do
      end do
    if(allocated(temp_dtlrdv)) deallocate(temp_dtlrdV)
    if(allocated(my_dtlrdv)) deallocate(my_dtlrdv)
    if(allocated(dTLRdV))    deallocate(dTLRdV)
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Sends for dTLRdh
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(do_forces) then
      do counter = 1, n_pairs(me_image+1), 1
        collected_dtlrdh(counter + offsets(me_image+1), :, :, :, :) = my_dtlrdh(counter, :, :, :, :)
        end do

      do i = 0, max_proc_forces-1, 1                      ! Outer loop is the loop for recvs
        do j = 0, max_proc_pairs-1, 1                   ! Loop for recvs
          if((j.eq.me_image).and.(i.ne.me_image)) then
            call mpi_send(my_dtlrdh,&                  ! Send buffer
                    my_num_pairs*3*3*3*3,&        ! Count of communicated data
                    mpi_double_precision,&    ! Data type
                    i,&                       ! Destination
                    j,&                       ! Tag (who it was sent from)
                    intra_image_comm,&        ! Communicator to send with
                    ierr)                     ! Error flag
          else
            if((i.eq.me_image).and.(j.ne.me_image)) then
              if(.not.allocated(temp_dtlrdh))    allocate(temp_dtlrdh(n_pairs(j+1), 3, 3, 3, 3)); temp_dtlrdh=0.0_DP
              call mpi_recv(temp_dtlrdh,&                    ! Recv buffer
                      n_pairs(j+1)*3*3*3*3,&            ! Count of communicated data
                      mpi_double_precision,&        ! Data type
                      j,&                           ! Source (who it was sent from)
                      j,&                           ! Tag (also who it was sent from)
                      intra_image_comm,&            ! Communicator to send with
                      rstatus,&                     ! Status flag
                      ierr)                         ! Error flag
              do counter = 1, n_pairs(j+1), 1
                collected_dtlrdh(counter + offsets(j+1), :, :, :, :) = temp_dtlrdh(counter, :, :, :, :)
                end do
              if(allocated(temp_dtlrdh)) deallocate(temp_dtlrdh)
              end if
            end if
          end do
        end do
      if(allocated(my_dtlrdh)) deallocate(my_dtlrdh)

      !!!!!!!!!!!!!!!!!!
      !Send and Receives for dtlrdr
      !!!!!!!!!!!!!!!!!!

      do counter_p = 1, n_pairs(me_image+1), 1
        do counter_f = 1, n_comps(me_image+1), 1
          collected_dtlrdr(offsets(me_image+1)+counter_p, :, :, counter_f, :) = &
                        my_dtlrdr(counter_p, :, :, force_offsets(me_image+1) + counter_f, :)
          end do
        end do

      do i = 0, max_proc_forces-1, 1                      ! Outer loop is the loop for recvs
        do j = 0, max_proc_pairs-1, 1                   ! Loop for recvs
          if((j.eq.me_image).and.(i.ne.me_image)) then
            if(.not.allocated(temp_dtlrdr))    allocate(temp_dtlrdr(n_pairs(j+1), 3, 3, n_comps(i+1), 3));
            temp_dtlrdr = 0.0_DP
            do counter_f = 1, n_comps(i+1), 1
              temp_dtlrdr(:,:,:,counter_f,:) = my_dtlrdr(:,:,:,force_offsets(i+1)+counter_f,:)
              end do
            call mpi_send(temp_dtlrdr,&                               ! Send buffer
                    n_pairs(j+1)*3*3*n_comps(i+1)*3,& ! Count of communicated data
                    mpi_double_precision,&                   ! Data type
                    i,&                                      ! Destination
                    j,&                                      ! Tag (who it was sent from)
                    intra_image_comm,&                       ! Communicator to send with
                    ierr)                                    ! Error flag
            if(allocated(temp_dtlrdr))  deallocate(temp_dtlrdr)
          else
            if((i.eq.me_image).and.(j.ne.me_image)) then
              if(.not.allocated(temp_dtlrdr))    allocate(temp_dtlrdr(n_pairs(j+1), 3, 3, n_comps(i+1), 3));
              call mpi_recv(temp_dtlrdr,&                              ! Recv buffer
                      n_pairs(j+1)*3*3*n_comps(i+1)*3,& ! Count of communicated data
                      mpi_double_precision,&                  ! Data type
                      j,&                                     ! Source (who it was sent from)
                      j,&                                     ! Tag (also who it was sent from)
                      intra_image_comm,&                      ! Communicator to send with
                      rstatus,&                               ! Status flag
                      ierr)                                   ! Error flag
              do counter = 1, n_pairs(j+1), 1
                collected_dtlrdr(counter + offsets(j+1), :, :, :, :) = temp_dtlrdr(counter, :, :, :, :)
                end do
              if(allocated(temp_dtlrdr))  deallocate(temp_dtlrdr)
              end if
            end if
          end do
        end do
        if(allocated(my_dtlrdr)) deallocate(my_dtlrdr)
        if(allocated(dTLRdr)) deallocate(dtlrdr)
        end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Putting it all together!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call stop_clock('mbd_cpq_comm')
    ! only loops over the upper triangle
    if(.not.allocated(Cpq))             allocate(Cpq(3*nat_sl, 3*nat_sl));              Cpq = 0.0_DP
    call start_clock('mbd_build_c')
    if(vdw_self_consistent) then
      if(.not.allocated(dTLRdV))          allocate(dTLRdV(3, 3, n_comps(me_image+1)))
      if(.not.allocated(dCpqdV))          allocate(dCpqdV(3*nat_sl, 3*nat_sl, n_comps(me_image+1)));   dCpqdV = 0.0_DP
      end if
    if(do_forces) then
      if(.not.allocated(dTLRdr))          allocate(dTLRdr(3, 3, n_comps(me_image+1), 3))
      if(.not.allocated(dCpqdR))          allocate(dCpqdR(3*nat_sl, 3*nat_sl, n_comps(me_image+1), 3));   dCpqdR = 0.0_DP
      if(.not.allocated(dCpqdh))          allocate(dCpqdh(3*nat_sl, 3*nat_sl, 3  ,                 3));   dCpqdh = 0.0_DP
      end if
  counter = 1
  counter_c = 1

  do p=1, nat_sl, 1
    do q=p, nat_sl, 1
      if(p.ne.q) then
          TLR = 0.0_DP
          dTLRdr = 0.0_DP
          dTLRdh = 0.0_DP
          aSqrt = dsqrt(alpha_0_sl(p)*alpha_0_sl(q))
          omMult = omega_scs_sl(p)*omega_scs_sl(q)
          prefactor = omMult*aSqrt
          ! call mbdvdw_TLR(p, q, TLR, dTLRdr, dTLRdh, dTLRdV)
          TLR = collected_tlr(counter, :, :)
          if(do_forces) dTLRdr = collected_dtlrdr(counter, :, :, :, :)
          if(do_forces) dTLRdh = collected_dtlrdh(counter, :, :, :, :)
          if(vdw_self_consistent) dTLRdV = collected_dtlrdV(counter, :, :, :)
          counter = counter + 1
        end if
      do i_idx=1, 3, 1
        do j_idx=1, 3, 1
            i_index = (3*p - 3 + i_idx)
            j_index = (3*q - 3 + j_idx)
            if(p.eq.q) then
              ! Energy bit
              if(i_idx.eq.j_idx) then
                Cpq(i_index, j_index) = omega_scs_sl(p)**2.0_DP
                if(do_forces.or.vdw_self_consistent) then
                  counter_c = 1
                  do s=1, nat, 1
                    ! Self consistent derivative of diagonal
                    if((f_cpu_id(s).eq.me_image).and.(vdw_self_consistent)) then
                      dCpqdV(i_index, j_index, counter_c) = 2.0_DP*omega_scs_sl(p)*domegadV_sl(p, s)
                      end if
                    if(do_forces) then
                      do i=1, 3, 1
                        ! dR Derivative
                        if(f_cpu_id(s).eq.me_image) dCpqdR(i_index,j_index,counter_c,i)=2.0_DP*omega_scs_sl(p)*domegadR_sl(p,s,i)
                        ! dH Derivative
                        if(s.le.3) dCpqdh(i_index,j_index,s,i)=2.0_DP*omega_scs_sl(p)*domegadh_sl(p,s,i)
                        end do
                      end if
                      if(f_cpu_id(s).eq.me_image) counter_c = counter_c + 1
                    end do
                  end if
                end if
              else
              ! Energy bit
              Cpq(i_index, j_index) = Cpq(i_index, j_index) +  prefactor*TLR(i_idx, j_idx)
              Cpq(j_index, i_index) = Cpq(i_index, j_index)
              if(do_forces.or.vdw_self_consistent) then
                counter_c = 1
                do s=1, nat, 1
                  if((f_cpu_id(s).eq.me_image).and.(vdw_self_consistent)) then
                    ! Self consistent derivative of diagonal
                    pre1=(omega_scs_sl(p)*domegadV_sl(q,s)+omega_scs_sl(q)*domegadV_sl(p,s))*asqrt
                    pre2=omMult*(alpha_0_sl(q)*dalpha_0dV_sl(p,s)+alpha_0_sl(p)*dalpha_0dV_sl(q,s))
                    pre2=pre2/(2.0_DP*asqrt)
                    der = (pre1+pre2)*TLR(i_idx,j_idx)+prefactor*dTLRdV(i_idx,j_idx,i)
                    dCpqdV(i_index,j_index,counter_c)=dCpqdV(i_index,j_index,counter_c)+der
                    dCpqdV(j_index,i_index,counter_c)=dCpqdV(i_index,j_index,counter_c)
                    end if
                  if(do_forces) then
                    do i=1, 3, 1
                      ! dR Derivative
                      if((f_cpu_id(s).eq.me_image)) then
                        pre1=omega_scs_sl(p)*domegadR_sl(q,s,i)+omega_scs_sl(q)*domegadR_sl(p,s,i)
                        pre1=pre1*asqrt
                        pre2=omMult
                        pre2=pre2*(alpha_0_sl(q)*dalpha_0dR_sl(p,s,i)+alpha_0_sl(p)*dalpha_0dR_sl(q,s,i))
                        pre2=pre2/(2.0_DP*asqrt)
                        der = (pre1+pre2)*TLR(i_idx,j_idx)+prefactor*dTLRdR(i_idx,j_idx,counter_c,i)
                        dCpqdR(i_index,j_index,counter_c,i)=dCpqdR(j_index,i_index,counter_c,i)+der
                        dCpqdR(j_index,i_index,counter_c,i)=dCpqdR(i_index,j_index,counter_c,i)
                        end if
                      ! dH Derivative
                      if(s.le.3) then
                        pre1=omega_scs_sl(p)*domegadh_sl(q,s,i)+omega_scs_sl(q)*domegadh_sl(p,s,i)
                        pre1=pre1*asqrt
                        pre2=omMult
                        pre2=pre2*(alpha_0_sl(q)*dalpha_0dh_sl(p,s,i)+alpha_0_sl(p)*dalpha_0dh_sl(q,s,i))
                        pre2=pre2/(2.0_DP*asqrt)
                        der = (pre1+pre2)*TLR(i_idx, j_idx)+prefactor*dTLRdh(i_idx,j_idx,s,i)
                        dCpqdh(i_index,j_index,s,i)=dCpqdh(i_index,j_index,s,i) + der
                        dCpqdh(j_index, i_index, s, i) = dCpqdh(i_index, j_index, s, i)
                        end if
                      end do
                    end if
                  if(f_cpu_id(s).eq.me_image) counter_c = counter_c + 1
                  end do
                end if
              end if
          end do
        end do
      end do
    end do

  call stop_clock('mbd_build_c')
  if(allocated(collected_tlr)) deallocate(collected_tlr)
  if(allocated(collected_dtlrdr)) deallocate(collected_dtlrdr)
  if(allocated(collected_dtlrdh)) deallocate(collected_dtlrdh)
  if(allocated(collected_dtlrdV)) deallocate(collected_dtlrdv)
  call stop_clock('mbd_construct_hamiltonian')
  return
  end subroutine mbdvdw_construct_hamiltonian

  subroutine mbdvdw_write_check1()
  implicit none
  integer :: p, s, i
  integer, parameter :: e_check = 50
  integer, parameter :: de_dr_check = 51
  integer, parameter :: de_dh_check = 52
  integer, parameter :: e_nonint_check = 53
  integer, parameter :: de_nonint_dr_check = 54
  integer, parameter :: de_nonint_dh_check = 55
  integer, parameter :: e_int_check = 56
  integer, parameter :: de_int_dr_check = 57
  integer, parameter :: de_int_dh_check = 58
  integer, parameter :: alpha_0_check = 59
  integer, parameter :: dalpha_0_dr_check = 60
  integer, parameter :: dalpha_0_dh_check = 61
  integer, parameter :: omega_scs_check = 62
  integer, parameter :: domega_scs_dr_check = 63
  integer, parameter :: domega_scs_dh_check = 64

  if(mbd_debug_dh.or.mbd_debug_dh) then
    open (unit=omega_scs_check, file=trim(dir)//"omega_scs.tmp",action=act,status=stat,position=pos)
    open (unit=alpha_0_check,   file=trim(dir)//"alpha_0.tmp",action=act,status=stat,position=pos)
    do i=1,nat,1
      write(omega_scs_check, *) i, omega_scs(i)
      write(alpha_0_check,   *) i, alpha_0(i)
      end do
    close(omega_scs_check)
    close(alpha_0_check)
    end if

  if(mbd_debug_dr) then
    open (unit=dalpha_0_dr_check,   file=trim(dir)//"dalpha_0_dr.tmp",action=act,status=stat,position=pos)
    open (unit=domega_scs_dr_check, file=trim(dir)//"domega_scs_dr.tmp",action=act,status=stat,position=pos)
      do p=1,nat,1
        do s=1,nat,1
          do i=1,3,1
            write(dalpha_0_dr_check,   *) p, s, i, dalpha_0dr(p, s, i)
            write(domega_scs_dr_check, *) p, s, i, domegadr(p, s, i)
          end do
        end do
      end do
    close(dalpha_0_dr_check)
    close(domega_scs_dr_check)
    end if

  if(mbd_debug_dh) then
    open (unit=dalpha_0_dh_check,   file=trim(dir)//"dalpha_0_dh.tmp",action=act,status=stat,position=pos)
    open (unit=domega_scs_dh_check, file=trim(dir)//"domega_scs_dh.tmp",action=act,status=stat,position=pos)
      do p=1,nat,1
        do s=1,3,1
          do i=1,3,1
            write(dalpha_0_dh_check, *) p, s, i, dalpha_0dh(p, s, i)
            write(domega_scs_dh_check, *) p, s, i, domegadh(p, s, i)
          end do
        end do
      end do
    close(dalpha_0_dh_check)
    close(domega_scs_dh_check)
    end if

  end subroutine mbdvdw_write_check1

  subroutine mbdvdw_write_check2()
  implicit none
  integer :: s, i
  integer, parameter :: e_check = 50
  integer, parameter :: de_dr_check = 51
  integer, parameter :: de_dh_check = 52
  integer, parameter :: e_nonint_check = 53
  integer, parameter :: de_nonint_dr_check = 54
  integer, parameter :: de_nonint_dh_check = 55
  integer, parameter :: e_int_check = 56
  integer, parameter :: de_int_dr_check = 57
  integer, parameter :: de_int_dh_check = 58
  integer, parameter :: alpha_0_check = 59
  integer, parameter :: dalpha_0_dr_check = 60
  integer, parameter :: dalpha_0_dh_check = 61
  integer, parameter :: omega_scs_check = 62
  integer, parameter :: domega_scs_dr_check = 63
  integer, parameter :: domega_scs_dh_check = 64

  if(mbd_debug_dh.or.mbd_debug_dr) then
    open (unit=e_int_check, file=trim(dir)//"e_int.tmp",action=act,status=stat,position=pos)
    open (unit=e_nonint_check, file=trim(dir)//"e_nonint.tmp",action=act,status=stat,position=pos)

    write(e_int_check, *) interacting_energy
    write(e_nonint_check, *) non_interacting_energy
    close(e_nonint_check)
    close(e_int_check)
    end if
  if(mbd_debug_dh) then
    open (unit=de_int_dh_check, file=trim(dir)//"de_int_dh.tmp",action=act,status=stat,position=pos)
    open (unit=de_nonint_dh_check, file=trim(dir)//"de_nonint_dh.tmp",action=act,status=stat,position=pos)
    do i=1,3,1
      do s=1,3,1
        write(de_int_dh_check, *) i, s, int_HmbdVdW(i,s)
        write(de_nonint_dh_check, *) i, s, nonint_HmbdVdW(i,s)
        end do
      end do
    close(de_int_dh_check)
    close(de_nonint_dh_check)
    end if
  if(mbd_debug_dr) then
    open (unit=de_int_dr_check, file=trim(dir)//"de_int_dr.tmp",action=act,status=stat,position=pos)
    open (unit=de_nonint_dr_check, file=trim(dir)//"de_nonint_dr.tmp",action=act,status=stat,position=pos)
    do i=1,3,1
      do s=1,3,1
        write(de_int_dr_check, *) i, s, int_HmbdVdW(i,s)
        write(de_nonint_dr_check, *) i, s, nonint_HmbdVdW(i,s)
        end do
      end do
    close(de_int_dr_check)
    close(de_nonint_dr_check)
    end if

  end subroutine mbdvdw_write_check2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE mbdvdw_wfforce()
  IMPLICIT NONE
  INTEGER :: ia,ip,iq,off1, Ndim, nr1, nr2, nr3, iproc
  REAL(DP), DIMENSION(:), ALLOCATABLE :: UmbdvdWA

  nr1=dfftp%nr1; nr2=dfftp%nr2; nr3=dfftp%nr3
  nr1r=nr1/2; nr2r=nr2/2; nr3r=nr3/2
  IF(MOD(nr1,2).EQ.1) nr1r=(nr1+1)/2
  IF(MOD(nr2,2).EQ.1) nr2r=(nr2+1)/2
  IF(MOD(nr3,2).EQ.1) nr3r=(nr3+1)/2
  Ndim=MAX(nr1*nr2,dffts%npp(me_bgrp+1)*nr1*nr2)

  CALL start_clock('mbdvdw_wfforce')
  ALLOCATE(UmbdvdWA(nr1*nr2*nr3)); UmbdvdWA=0.0_DP
  DO iproc=1,nstates(me)
  ia=me+nproc_image*(iproc-1)
  DO iq=1,NsomegaA(ia)
    off1=somegaA(iq,1,iproc)+(somegaA(iq,2,iproc)-1)*nr1+(somegaA(iq,3,iproc)-1)*nr1*nr2    !global offset [nr1,nr2,nr3]
    UmbdvdWA(off1)=UmbdvdWA(off1)+Uprefactor(ia)*dveffAdn(iq,iproc)
  END DO
  END DO
  CALL mp_sum(UmbdvdWA,intra_image_comm)

  IF (dffts%npp(me_bgrp+1).NE.0) THEN
  DO ip=1,dffts%npp(me_bgrp+1)*nr1*nr2
    UmbdvdW(ip)=UmbdvdWA(ip+rdispls(me_bgrp+1))
  END DO
  END IF
  IF (ALLOCATED(UmbdvdWA))      DEALLOCATE(UmbdvdWA)
  CALL stop_clock('mbdvdw_wfforce')
  RETURN
  END SUBROUTINE mbdvdw_wfforce

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mbdvdw_calculate(tauin, rhor, at_x_alat, bg_trans_by_alat)
     implicit none
    ! IO Variables
    ! Atomic density
    real(dp), intent(in) :: rhor(:,:)
    ! Coordinates
    real(dp) :: tauin(3, nat)
    REAL(DP), INTENT(IN), OPTIONAL :: at_x_alat(3,3)
    REAL(DP), INTENT(IN), OPTIONAL :: bg_trans_by_alat(3,3)

    ! Local variables
    real(dp) :: toAng
    integer :: i_freq, p, s, i, i_idx, i_atom
    real(dp), dimension(nat) :: alpha_temp
    real(dp), dimension(nat, nat, 3) :: dalpha_tempdR
    real(dp), dimension(nat, 3  , 3) :: dalpha_tempdh
    real(dp), dimension(nat, nat   ) :: dalpha_tempdV
    real(dp), dimension(3,3) :: Htmp
    integer, parameter :: mbd_check=21
    integer, parameter :: fmbd_outfile=22
    integer :: ierr
    integer, dimension(MPI_STATUS_SIZE) :: rstatus
    logical :: alpha_zero_failed

    ! DEBUG
    integer, parameter :: e_check = 50
    integer, parameter :: de_dr_check = 51
    integer, parameter :: de_dh_check = 52
    integer, parameter :: e_nonint_check = 53
    integer, parameter :: de_nonint_dr_check = 54
    integer, parameter :: de_nonint_dh_check = 55
    integer, parameter :: e_int_check = 56
    integer, parameter :: de_int_dr_check = 57
    integer, parameter :: de_int_dh_check = 58
    integer, parameter :: alpha_0_check = 59
    integer, parameter :: dalpha_0_dr_check = 60
    integer, parameter :: dalpha_0_dh_check = 61
    integer, parameter :: omega_scs_check = 62
    integer, parameter :: domega_scs_dr_check = 63
    integer, parameter :: domega_scs_dh_check = 64

    call start_clock('mbd')

    if (present(at_x_alat)) then
      ! PW case
      h_(:,:) = at_x_alat(:,:)
      ainv_(:,:) = bg_trans_by_alat(:,:)
      else
      ! CP case
      h_(:,:) = h_(:,:)
      ainv_(:,:) = ainv_(:,:)
      end if

    EmbdvdW = 0.0_DP
    FmbdVdW = 0.0_DP
    HmbdVdW = 0.0_DP
    Umbdvdw = 0.0_DP
    alpha_zero_failed = .false.

    tau = tauin

    toAng = 0.52917721092_DP

    do_forces = mbd_vdw_forces
    ! write(stdout, *) do_forces, mbd_vdw_forces, vdw_self_consistent, mbd_conv_elec, mbd_first_step
    if(mbd_first_step) goto 10
    if((.not.vdw_self_consistent).and.(.not.mbd_conv_elec)) goto 10
    if(vdw_self_consistent.and.mbd_vdw_forces) then
      if(mbd_conv_elec) then
        do_forces = .true.
        else
        do_forces = .false.
        end if
      end if

    call start_clock('mbd_firstpinit')
    call mbdvdw_para_init()
    call stop_clock('mbd_firstpinit')

    ! Sets up the grid and computes the derivatives with respect to volume
    call start_clock('mbd_veff')
    call tsvdw_calculate(tau, rhor, at_x_alat, bg_trans_by_alat)
    if(mbd_debug_dh) call set_dveffdh()
    !  call set_dveffdr()
    call stop_clock('mbd_veff')

    call mbdvdw_zeroout()
    call start_clock('mbd_fparainit')
    call mbdvdw_para_init_forces()
    call stop_clock('mbd_fparainit')
    call start_clock('mbd_pbc')
    call mbdvdw_pbc(tau, h_, ainv_, tau_s)
    call stop_clock('mbd_pbc')

    ! print out FHI-aims test file in Ang
    open (unit=mbd_check, file=trim(dir)//"FHI-aims.xyz",action=act,status="replace")
    write(mbd_check, *) nat
    write(mbd_check, *)
    do p=1,nat,1
        write(mbd_check, '(A, F9.6, F9.6, F9.6, F9.6)') atm(ityp(p)), tau(1, p)*toAng, &
                                        tau(2, p)*toAng, tau(3, p)*toAng, VefftsvdW(p)/vfree(ityp(p))
    end do
    ! write(mbd_check, '("lattice_vector", F9.6,"  ", F9.6,"  ", F9.6)'), h_(1,1)*toAng, h_(1,2)*toAng, h_(1,3)*toAng
    ! write(mbd_check, '("lattice_vector", F9.6,"  ", F9.6,"  ", F9.6)'), h_(2,1)*toAng, h_(2,2)*toAng, h_(2,3)*toAng
    ! write(mbd_check, '("lattice_vector", F9.6,"  ", F9.6,"  ", F9.6)'), h_(3,1)*toAng, h_(3,2)*toAng, h_(3,3)*toAng
    close(mbd_check)

    ! Frequency dependent loop to perform the casimir polder integral to find the effective
    ! C6 coefficients. With these in hand, we'll toss this over to the energy expression

    call start_clock('mbd_loop')
    do i_freq = 1,npts,1
      ! write(stdout, '(3X, "omega =  ", F14.6)'), casimir_omega(i_freq)
      omega_to_print = casimir_omega(i_freq)
      alpha_temp = 0.0_DP
      dalpha_tempdR = 0.0_DP
      dalpha_tempdh = 0.0_DP
      ! Compute $R_{VdW}$
      call mbdvdw_effqts(casimir_omega(i_freq))
      ! Performs the self consistent screening
      call mbdvdw_SCS()
        ! Computes the contracted alphas
        if(i_freq.eq.1) then
          call mbdvdw_calculate_screened_pol(alpha_0, dalpha_0dR, dalpha_0dh, dalpha_0dV)
          if(me_image.eq.root_image) then
            do i_idx=1,nat,1
              if(alpha_0(i_idx).le.0.0) then
                write(stdout, '(3X, "Negative alpha! Aborting this iteration")')
                EmbdvdW = 0.0_DP
                Fmbdvdw = 0.0_DP
                Hmbdvdw = 0.0_DP
                alpha_zero_failed = .true.
                GOTO 20
                end if
              end do
              20 CONTINUE
              do i=1, nproc_image-1,1
                call mpi_send(alpha_zero_failed, 1, mpi_logical,&
                    i, 0, intra_image_comm,ierr)
                end do
              if(alpha_zero_failed) GOTO 10
            end if
          else
          alpha_temp = 0.0_DP
          call mbdvdw_calculate_screened_pol(alpha_temp, dalpha_tempdR, dalpha_tempdh, dalpha_tempdV)
          ! Performs the casimir polder integral
          do p=1, nat, 1
            omega_scs(p)=omega_scs(p)+(4.0_DP/pi)*casimir_omega_weight(i_freq)*&
                                    alpha_temp(p)**2.0_DP/(alpha_0(p)**2.0_DP)
            if(do_forces.or.vdw_self_consistent) then
              do s=1,nat,1
                if(f_cpu_id(s).eq.me_image) then
                  if(vdw_self_consistent) then
                    domegadV(p, s) =  domegadV(p, s) + casimir_omega_weight(i_freq)*(8.0_DP/pi)*&
                              (alpha_temp(p)*dalpha_tempdV(p, s) - 1.0_DP/alpha_0(p)*&
                                dalpha_0dV(p, s)*alpha_temp(p)**2.0_DP)/(alpha_0(p)**2.0)
                    end if
                  if(do_forces) then
                    do i=1,3,1
                      domegadR(p, s, i) =  domegadR(p, s, i) + casimir_omega_weight(i_freq)*(8.0_DP/pi)*&
                              (alpha_temp(p)*dalpha_tempdR(p, s, i) - 1.0_DP/alpha_0(p)*&
                                dalpha_0dR(p, s, i)*alpha_temp(p)**2.0_DP)/(alpha_0(p)**2.0)
                      if(s.le.3) then
                        domegadh(p, s, i) =  domegadh(p, s, i) + casimir_omega_weight(i_freq)*(8.0_DP/pi)*&
                                (alpha_temp(p)*dalpha_tempdh(p, s, i) - 1.0_DP/alpha_0(p)*&
                                  dalpha_0dh(p, s, i)*alpha_temp(p)**2.0_DP)/(alpha_0(p)**2.0)
                        end if
                      end do
                    end if
                  end if
                end do
              end if
            end do
          end if
      if((me_image.ne.root_image).and.(i_freq.eq.1)) then
        call mpi_recv(alpha_zero_failed,& ! Buffer
                1,& ! Count
                mpi_logical,& ! Type
                0,& ! Source
                0,& ! Tag
                intra_image_comm,& ! Communicator
                rstatus,& ! Status var
                ierr) ! Error flag
        if(alpha_zero_failed) then
          EmbdvdW = 0.0_DP
          Fmbdvdw = 0.0_DP
          Hmbdvdw = 0.0_DP
          GOTO 10
          end if
        end if
      end do
      call stop_clock('mbd_loop')
      ! Compute the noninteracting energy.
      ! This is done by simply summing over the effective frequencies.
      ! \sum_p \bar{\omega}_p
      ! Where \bar{\omega}_p is given by
      ! \bar{\omega}_p = \frac{4}{\pi} \frac{ \int_0^{\infty} \left[ \alpha_p(i \omega) \right]^2 d\omega  }{ \left[ \alpha_p^0 \right]^2 }.
    call mbdvdw_noninteracting_energy(non_interacting_energy, nonint_FmbdVdW, nonint_HmbdVdW, nonint_Uprefactor)
    if(vdw_self_consistent) call mp_sum(nonint_Uprefactor, intra_image_comm)
    if(mbd_vdw_forces) call mp_sum(nonint_FmbdVdW, intra_image_comm)
    if(mbd_vdw_forces) call mp_sum(nonint_HmbdVdW, intra_image_comm)

    if(mbd_debug_dh.or.mbd_debug_dr) call mbdvdw_write_check1()
    call mbdvdw_cleanup_postscs()
    !call mp_sum(alpha_0, intra_image_comm)
    !call mp_sum(omega_scs, intra_image_comm)

    if(vdw_self_consistent) then
      call mp_sum(dalpha_0dV, intra_image_comm)
      call mp_sum(domegadV, intra_image_comm)
      end if

    if(do_forces) then
      call mp_sum(dalpha_0dR, intra_image_comm)
      call mp_sum(domegadR, intra_image_comm)
      end if

    ! This method constructs C so that we can find its eigenvalues and sum them to find the interacting
    ! parts of the energy. This is done by constructing
    ! \mathcal{C}^{MBD}_{pq} = \delta_{pq} \bar{\omega}^2_p + \left(1 - \delta_pq \right) \bar{\omega}_p \bar{\omega}_q \sqrt{\bar{\alpha}_p^0 \bar{\alpha}_q^0} T_{pq}
    call mbdvdw_construct_hamiltonian()
    ! Computes the interacting energy by diagonalizing the MBD hamiltonian
    call mbdvdw_interacting_energy(interacting_energy, int_FmbdVdW, int_HmbdVdW, int_Uprefactor)

    if(vdw_self_consistent) then
      Uprefactor = 0.5_DP*int_Uprefactor - 1.5_DP*nonint_Uprefactor
      call mbdvdw_wfforce()
      end if

    if(me_image.eq.root_image) then
      EmbdvdW = 0.5_DP*interacting_energy - 1.5_DP*non_interacting_energy
      FmbdVdW = 0.5_DP*int_FmbdVdW - 1.5_DP*nonint_FmbdVdW
      HmbdVdW = 0.5_DP*int_HmbdVdW - 1.5_DP*nonint_HmbdVdW
      HmbdVdW = -1.0_DP*Hmbdvdw

      ! write(stdout, *), "MBD Energy is: ", EmbdvdW

      if(do_forces) then
        write(stdout, '(3X, "FmbdVdW")')
        write(stdout, '(3X, "DX       DY       DZ")')
        do s = 1, nat, 1
          write(stdout, '(3X, F14.10, F14.10, F14.10)') FmbdVdW(s, 1), FmbdVdW(s, 2), FmbdVdW(s, 3)
        end do

        Htmp = Hmbdvdw

        write(stdout, '(3X, "HmbdVdW")')
        write(stdout, '(3X, "DX       DY       DZ")')
        do s = 1, 3, 1
          write(stdout, '(3X, F14.10, F14.10, F14.10)') Htmp(s, 1), Htmp(s, 2), Htmp(s, 3)
        end do
        Htmp = 1.0_DP*MATMUL( HmbdvdW, transpose(h_) )/omega

      write(stdout, '(3X, "HtsVdW contracted with fractional coordinates")')
      write(stdout, '(3X, "DX       DY       DZ")')
      do s = 1, 3, 1
        write(stdout, '(3X, F14.10, F14.10, F14.10)') Htmp(s, 1), Htmp(s, 2), Htmp(s, 3)
      end do

        Htmp = 1.0_DP*MATMUL( HmbdvdW, h_ )/omega

        write(stdout, '(3X, "sigma_mbdVdW")')
        write(stdout, '(3X, "DX       DY       DZ")')
        do s = 1, 3, 1
          write(stdout, '(3X, F14.10, F14.10, F14.10)') Htmp(s, 1), Htmp(s, 2), Htmp(s, 3)
        end do

        end if

      if(mbd_debug_dh.or.mbd_debug_dr) call mbdvdw_write_check2()
      end if

    10 CONTINUE

    if(alpha_zero_failed) then
      Embdvdw = 0.0_DP
      Fmbdvdw = 0.0_DP
      Hmbdvdw = 0.0_DP
      Umbdvdw = 0.0_DP
      end if

    call MPI_Bcast(EmbdvdW, & ! buffer
             1,       & ! count
             mpi_double_precision, & ! MPI_Datatype datatype
             root_image, & !root
             intra_image_comm, & ! Communicator
             ierr)

    call MPI_Bcast(FmbdvdW, & ! buffer
             3*nat,   & ! count
             mpi_double_precision, & ! MPI_Datatype datatype
             root_image, & !root
             intra_image_comm, & ! Communicator
             ierr)

    call MPI_Bcast(HmbdvdW, & ! buffer
             3*3,     & ! count
             mpi_double_precision, & ! MPI_Datatype datatype
             root_image, & !root
             intra_image_comm, & ! Communicator
             ierr)


    call mbdvdw_cleanup()
    call tsvdw_cleanup()
    call tsvdw_cleanup_post_mbd()

    if(iverbosity.ge.1.and.mbd_conv_elec) then
    write(stdout, '(3X, "Done with MBD@rsSCS")')
    write(stdout, '(3X, "---------------------------")')
    end if

    call stop_clock('mbd')

  end subroutine mbdvdw_calculate

subroutine set_dveffdr()
implicit none

dveffdr(1, 1, 1) = -0.000227418632288_DP
dveffdr(1, 1, 2) = -0.000357282236816_DP
dveffdr(1, 1, 3) = 0.0016473119677_DP
dveffdr(1, 2, 1) = -3.25486414089_DP
dveffdr(1, 2, 2) = 4.76870534128e-06_DP
dveffdr(1, 2, 3) = 1.00075084622e-06_DP
dveffdr(1, 3, 1) = 3.2549824439_DP
dveffdr(1, 3, 2) = 0.000147752408957_DP
dveffdr(1, 3, 3) = 0.00017393131864_DP
dveffdr(2, 1, 1) = -0.574911184382_DP
dveffdr(2, 1, 2) = -0.000332414623576_DP
dveffdr(2, 1, 3) = 0.00109288209618_DP
dveffdr(2, 2, 1) = -0.63819758671_DP
dveffdr(2, 2, 2) = 6.40659119519e-05_DP
dveffdr(2, 2, 3) = -5.53277403454e-05_DP
dveffdr(2, 3, 1) = 1.20926191825_DP
dveffdr(2, 3, 2) = 0.000102710365261_DP
dveffdr(2, 3, 3) = 0.000144707089972_DP
dveffdr(3, 1, 1) = 0.574651499469_DP
dveffdr(3, 1, 2) = -0.000260664507041_DP
dveffdr(3, 1, 3) = 0.00101143311561_DP
dveffdr(3, 2, 1) = -1.20935563453_DP
dveffdr(3, 2, 2) = -3.92560081257e-05_DP
dveffdr(3, 2, 3) = 5.14700362787e-05_DP
dveffdr(3, 3, 1) = 0.638165561485_DP
dveffdr(3, 3, 2) = 0.000102865875411_DP
dveffdr(3, 3, 3) = 8.00748540168e-05_DP

end subroutine set_dveffdr

subroutine set_dveffdh()
implicit none

dveffdh(1, 1, 1) = -0.656265647065_DP
dveffdh(1, 1, 2) = 0.162840716052_DP
dveffdh(1, 1, 3) = 0.102944859057_DP
dveffdh(1, 2, 1) = 0.162687500554_DP
dveffdh(1, 2, 2) = -0.524452439217_DP
dveffdh(1, 2, 3) = -0.166438619598_DP
dveffdh(1, 3, 1) = 0.103020726315_DP
dveffdh(1, 3, 2) = -0.166367620479_DP
dveffdh(1, 3, 3) = -0.438694944172_DP
dveffdh(2, 1, 1) = -0.10994354915_DP
dveffdh(2, 1, 2) = 0.0669329833973_DP
dveffdh(2, 1, 3) = 0.0426923520115_DP
dveffdh(2, 2, 1) = 0.0667903860133_DP
dveffdh(2, 2, 2) = -0.260241509335_DP
dveffdh(2, 2, 3) = -0.0667801013765_DP
dveffdh(2, 3, 1) = 0.0427335394484_DP
dveffdh(2, 3, 2) = -0.0667409363987_DP
dveffdh(2, 3, 3) = -0.267242458017_DP
dveffdh(3, 1, 1) = -0.255929516787_DP
dveffdh(3, 1, 2) = 0.0650670407752_DP
dveffdh(3, 1, 3) = 0.038170853035_DP
dveffdh(3, 2, 1) = 0.0650823549224_DP
dveffdh(3, 2, 2) = -0.159106952966_DP
dveffdh(3, 2, 3) = -0.0697664638229_DP
dveffdh(3, 3, 1) = 0.0382107539494_DP
dveffdh(3, 3, 2) = -0.0697278678082_DP
dveffdh(3, 3, 3) = -0.228597121014_DP

end subroutine set_dveffdh

end module mbdvdw_module
