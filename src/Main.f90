!========================================================================================
!            Multiple Component Grand Canonical Monte Carlo Nucleation Code
!========================================================================================
!
! Authors:
!   - Original Developer: Troy Loeffler
!     Bin Chen Research Group, Louisiana State University (LSU)
!   - Extended by: Aliasghar Sepehri
!     For simulation of biomolecular condensate nucleation
!
!----------------------------------------------------------------------------------------
! PURPOSE:
!   This FORTRAN code performs Grand Canonical Monte Carlo (GCMC) simulations for 
!   multicomponent chemical systems to study nucleation phenomena. It is designed to
!   support both small molecule nucleation and biomolecular condensate formation.
!
!   Advanced Monte Carlo techniques implemented include:
!     - Aggregation-Volume-Bias Monte Carlo (AVBMC)
!     - Configurational-Bias Monte Carlo (CBMC)
!     - Free-Endpoints CBMC (FECBMC)

!----------------------------------------------------------------------------------------
! EXTENSIONS FOR BIOMOLECULAR CONDENSATES:
!   The code has been adapted to include coarse-grained one-bead-per-residue models
!   for proteins and RNA using various force fields.
!
!   Supported force fields (with references):
!
!   - Mpipi for proteins and RNA:
!     Joseph, J.A. et al., Nat Comput Sci 1, 732–743 (2021)
!     https://doi.org/10.1038/s43588-021-00155-3
!
!   - HPS-KR, HPS-KH:
!     Dignon, G.L. et al., PLoS Comput. Biol. 14, e1005941 (2018)
!
!   - HPS-cation-pi-i and HPS-cation-pi-ii:
!     Das, S. et al., Proc. Natl Acad. Sci. USA 117, 28795–28805 (2020)
!
!   - HPS-FB:
!     Dannenhoffer-Lafage, T. & Best, R.B., J. Phys. Chem. B 125, 4046–4056 (2021)
!
!   - HPS-TSCL-M2:
!     Tesei, G. et al., Proc. Natl Acad. Sci. USA 118, e2111696118 (2021)
!
!   - HPS-Urry:
!     Regy, R.M. et al., Protein Sci. 30, 1371–1379 (2021)
!
!   - HPS-Nucl:
!     Li, L. et al., J. Phys. Chem. Lett. 14, 1748–1755 (2023)
!
!----------------------------------------------------------------------------------------
! CLUSTER FORMATION CRITERIA:
!   - Original criteria: energy-based and distance-based (between first atoms)
!   - New: minimum inter-residue distance (see:
!     Zerze, G.H., J. Chem. Theory Comput. 2024, 20 (4), 1646–1655)
!
!----------------------------------------------------------------------------------------
! PHYSICAL ANALOGY:
!   - The simulation uses an implicit solvent model
!   - The system mimics vapor-liquid nucleation
!
!========================================================================================
!==============================================================
! Program: Nucleation
! Purpose: Main driver for Monte Carlo simulations of nucleation
!==============================================================
program Nucleation
  ! Step 0: Declare modules and dependencies for the serial Monte Carlo nucleation simulation
  ! - Includes precision, constants, simulation parameters, input parsing, coordinates, energy,
  !   clustering criteria, move types, histograms, WHAM, umbrella sampling, and analysis
  use VarPrecision, only: dp
  use Constants, only: pi
  use ParallelVar, only: myid, p_size, ierror, provided, nout
  use SimParameters, only: distCriteria, minDistCriteria, calcPressure, ncycle, ncycle2, nMolTypes, &
                           NMin, NMax, NPart, isActive, NTotal, maxMol, vmdAtoms, avbmc_vol, &
                           NHist, P_Avg, gas_dens, max_dist, max_dist_single, max_rot, temperature, &
                           Dist_Critr, Dist_Critr_sq, Eng_Critr, outputEngUnits, outputEConv, &
                           outFreq_Screen, outFreq_Traj, MolName, pressure, prevMoveAccepted
  use ScriptInput, only: Script_ReadParameters
  use Coords, only: MolArray, NeighborList, gasConfig
  use CoordinateFunctions
  use Forcefield, only: nAtoms, atomData, atomArray, nBonds, bondArray, nAngles, &
                        bendArray, nTorsional, torsArray, r_min
  use EnergyTables, only: ETable, NeiETable, E_Inter_T, E_NBond_T, E_Stretch_T, E_Bend_T, E_Torsion_T, &
                          AlphaTargetIns, AlphaTargetDel
  use EnergyPressurePointers, only: Detailed_ECalc, Detailed_PCalc
  use MoveTypeModule, only: nMoveTypes, mcMoveArray, moveProbability, movesAccepted, movesAttempt, moveName, &
                            avbmcUsed, cbmcUsed
!  use CBMC_Initialize, only: CBMC_CreateTopology
  use AcceptRates, only: acptTrans, acptRot, atmpTrans, atmpRot, acptSwapIn, atmpSwapIn, &
                         acptSwapOut, atmpSwapOut, angGen_accpt, angGen_atmp, dihedGen_accpt, &
                         dihedGen_atmp, distGen_accpt, distGen_atmp, clusterCritRej, nCBMCmax, &
                         nFECBMCmax, acptCBMC, atmpCBMC, acptFECBMC, atmpFECBMC
  use AVBMC_RejectionVar, only: totalRej, ovrlapRej, dbalRej, critriaRej, boundaryRej, &
                                totalRej_out, dbalRej_out, critriaRej_out, boundaryRej_out
  use CBMC_Variables, only: regrowType
  use MiscelaniousVars, only: CollectHistograms
  use WHAM_Functions, only: WHAM_AdjustHist, WHAM_Finalize
  use WHAM_Module
  use UmbrellaSamplingNew, only: useUmbrella, curUIndx, UmbrellaHistAdd, ScreenOutputUmbrella, &
                                 OutputUmbrellaHist, CheckInitialValues, energyAnalytics, &
                                 OutputUmbrellaAnalytics
  use AnalysisMain, only: useAnalysis, PostMoveAnalysis, OutputAnalysis
  use MPI
  implicit none

      
  ! Step 1: Declare variables for simulation control, energy, pressure, histograms, and MPI
  ! - Includes flags, counters, energy/pressure variables, random numbers, file formats, and arrays
  logical :: errRtn
  logical :: screenEcho
  integer(kind=8) :: iCycle, iMove
  integer :: i, j, seed, ios
  integer :: AllocateStatus
  integer :: nSel
  integer :: getBiasIndex
  real(dp) :: E_T, E_Final, E_Debug
  real(dp) :: E_Inter_Final
  real(dp) :: E_Bend_Final
  real(dp) :: E_Torsion_Final
  real(dp) :: E_Stretch_Final
  real(dp) :: E_NBond_Final
  real(dp) :: P_Final
  integer :: pressLimit
  integer :: curIndx, maxIndx
  real(dp) :: TimeStart, TimeFinish
  real(dp) :: grnd, ran_num
  real(dp) :: dist_limit, rot_limit
  real(dp) :: check, norm
  real(dp), allocatable :: P_avg_Sum(:)
  real(dp), allocatable :: NHist_Sum(:)
  character(len=100) :: format_string, fl_name, out1
  character(len=1500) :: outFormat1, outFormat2
  logical, allocatable, dimension(:,:) :: FinalNeighborList
  real(dp), allocatable, dimension(:) :: FinalETable
  real(dp), allocatable, dimension(:) :: FinalNeiETable
  integer :: status(MPI_STATUS_SIZE)
      
  ! Step 2: Initialize MPI with threading support
!      MPI CONTROL STATEMENTS.  This block initializes the MPI threads and assigns each process
!      a thread ID (myid) and collects the number of total processes (p_size)
  call MPI_INIT(ierror)
  ! - Sets up MPI environment with MPI_THREAD_FUNNELED, assigns process IDs (myid), and checks threading level
!  call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, provided, ierror)
!  if (provided < MPI_THREAD_FUNNELED) then
!    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)
!    if (myid == 0) then
!      print *, 'Error: MPI_THREAD_FUNNELED not supported, provided level: ', provided
!    end if
!    call MPI_FINALIZE(ierror)
!    stop
!  end if
  call MPI_COMM_SIZE(MPI_COMM_WORLD, p_size, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)  

  ! Step 3: Set up unique output file names for each MPI process
  ! - Generates file names with process ID padding (e.g., Debug_<myid>.txt) and opens debug/report files
  if (myid < 10) then
    format_string = "(A,I1,A)"
  elseif (myid < 100) then
    format_string = "(A,I2,A)"
  elseif (myid < 1000) then
    format_string = "(A,I3,A)"
  elseif (myid < 10000) then
    format_string = "(A,I4,A)"
  else
    format_string = "(A,I5,A)"
  endif
  nout = 100 + myid
  if (p_size == 1) then
    nout = 6
  endif
  write(fl_name, format_string) "Debug_", myid, ".txt"
  open(unit=34, file=trim(adjustl(fl_name)))
  write(fl_name, format_string) "Final_Report_", myid, ".txt"
  open(unit=35, file=trim(adjustl(fl_name)))
      
  ! Step 4: Read input parameters from script
  ! - Parses user-provided script to set simulation parameters and random seed
  call Script_ReadParameters(seed, screenEcho)

  ! Step 5: Initialize simulation parameters and counters
  ! - Sets up output units, move counters, energy, AVBMC geometry, and output formats
  if(p_size .eq. 1) then
    write(*,*) "Input Script Complete!"
  else
    write(nout,*) "Input Script Complete!"    
  endif
  write(34,*) ">> Finished Script_ReadParameters"
  if(screenEcho) then      
    if(myid .eq. 0) then
      nout = 6        
    endif
  endif         

  movesAccepted     = 0E0_dp
  movesAttempt      = 1E-30_dp
  distGen_accpt     = 0E0_dp
  angGen_accpt      = 0E0_dp
  dihedGen_accpt    = 0E0_dp
  acptRot           = 0E0_dp
  distGen_atmp      = 0E0_dp
  angGen_atmp       = 0E0_dp
  dihedGen_atmp     = 0E0_dp
  acptTrans         = 0E0_dp
  atmpTrans         = 1E-40_dp
  atmpRot           = 1E-40_dp
  acptCBMC          = 0E0_dp
  acptFECBMC        = 0E0_dp
  atmpCBMC          = 1E-30_dp
  atmpFECBMC        = 1E-30_dp
  acptSwapIn        = 0E0_dp
  acptSwapOut       = 0E0_dp
  atmpSwapIn        = 1E-30_dp
  atmpSwapOut       = 1E-30_dp
  clusterCritRej    = 0E0_dp
  totalRej          = 0E0_dp
  ovrlapRej         = 0E0_dp
  dbalRej           = 0E0_dp
  critriaRej        = 0E0_dp
  boundaryRej       = 0E0_dp
  totalRej_out      = 0E0_dp
  dbalRej_out       = 0E0_dp
  critriaRej_out    = 0E0_dp
  boundaryRej_out   = 0E0_dp
  max_dist          = 0.05E0_dp
  max_rot           = 0.05E0_dp * pi
  max_dist_single   = 0.01E0_dp
  dist_limit        = 2E0_dp
  rot_limit         = pi
  prevMoveAccepted  = .false.
  Dist_Critr_sq     = Dist_Critr * Dist_Critr
  avbmc_vol         = (4E0_dp / 3E0_dp) * pi * Dist_Critr**3

  write(34,*) ">> Initialized MC move counters and AVBMC geometry"
  write(out1,"(A,I2,A)") "(",(4 + nMoveTypes + nMolTypes),"(A))"
  write(outFormat1, out1) "(", "2x,I9", (",2x,I5", i=1,nMolTypes), ",E17.6,2x", (",2x,F6.2", i=1,nMoveTypes), ")"
  write(outFormat2, out1) "(", "2x,I9", (",2x,I5", i=1,nMolTypes), ",F17.4,2x", (",2x,F6.2", i=1,nMoveTypes), ")"
  write(34,*) ">> Initialized output format strings"
  E_T = 0E0_dp
  write(34,*) ">> Total energy E_T initialized"

  ! Step 6: Initialize random number generator
  ! - Seeds RNG using the base seed (single process)
  seed = p_size * seed + myid
  call sgrnd(seed)
  write(34,*) ">> RNG initialized with seed =", seed
  
  ! Step 7: Set up molecular topology and initial configuration
  ! - Creates CBMC topology, writes PSF file, and computes initial energy
  call CBMC_CreateTopology
  flush(34)
  close(34)
  if (myid .eq. 0) then
    call WritePSF
  endif
  call Detailed_ECalc(E_T, errRtn)
  if(calcPressure) then
    call Detailed_PCalc(pressure)
    pressLimit = 1
    do i = 1, nMolTypes
      pressLimit = pressLimit * (NMAX(i) + 1)
    enddo
    allocate( NHist(1:pressLimit) )
    allocate( NHist_Sum(1:pressLimit) )
    allocate( P_Avg(1:pressLimit) )
    allocate( P_Avg_Sum(1:pressLimit) )
    NHist = 0E0_dp
    NHist_Sum = 0E0_dp
    P_Avg = 0E0_dp
    P_Avg_Sum = 0E0_dp
  endif
  if (errRtn) then
    write(nout, '(A)') 'ERROR: Initialization failed.'
    stop
  end if

  ! Step 8: Output initial energy table and open trajectory file
  ! - Writes initial energy table for active molecules and opens trajectory file
  write(35, '(A)') 'Initial Energy Table:'
  do i = 1, maxMol
    if (isActive(i) .and. ETable(i) /= 0.0_dp) then
      write(35, '(I6, ES15.6)') i, ETable(i)
    end if
  end do
  write(35, *)
  write(fl_name, format_string) 'Traj', myid, '.xyz'
  open(unit=30, file=trim(fl_name), form='FORMATTED', status='REPLACE', iostat=ios)
  if (ios /= 0) then
    write(nout, '(A, I0)') 'ERROR: Failed to open trajectory file ', ios
    stop
  end if     

  ! Step 9: Write initial trajectory and initialize energy accumulators
  ! - Print Dummy frame to VMD so VMD will correctly display varying cluster sizes. 
  call InitialTrajOutput
  ! - Outputs initial trajectory frame and resets energy components
  call TrajOutput(iCycle, E_T) 

  ! Step 10: Output simulation parameters
  ! - Prints simulation, WHAM, clustering, energy, AVBMC, and thermodynamic parameters
  write(nout, '(A)') 'Simulation Parameters:'
  write(nout, '(A, I0)') 'MPI threads: ', p_size
  write(nout, '(A, I0)') 'Thread ID: ', myid
  write(nout, '(A, I0)') 'Random Seed: ', seed
  write(nout, '(A, I0)') 'Number of Molecule Types: ', nMolTypes
  write(nout, '(A, *(I0, 1X))') 'Number of Initial Particles: ', NPART
  write(nout, '(A, I0)') 'Minimum Number of Particles: ', NMIN
  write(nout, '(A, I0)') 'Maximum Number of Particles: ', NMAX
  write(nout, '(A, *(I0, 1X))') 'Number of Atoms per Molecule: ', nAtoms
  write(nout, '(A, I0)') 'Number of Cycles: ', ncycle
  write(nout, '(A, I0)') 'Number of Moves per Cycle: ', ncycle2
  write(nout, '(A, I0)') 'CBMC Regrow Type: ', regrowType
  if (useWham) then
    write(nout, '(A)') '----- WHAM Parameters -----'
    write(nout, '(A, I0)') 'Expected number of WHAM iterations: ', nWhamItter
    write(nout, *) 'Reference Bin: ', refSizeNumbers
    write(nout, '(A, I0)') 'Reference Bin Index: ', refBin
    write(nout, '(A, I0)') 'Number of Cycles before WHAM adjustment: ', intervalWHAM
    write(nout, '(A, ES10.3)') 'WHAM Convergence Tolerance: ', tolLimit
    write(nout, '(A, I0)') 'Maximum WHAM Iterations: ', maxSelfConsist
    write(nout, '(A, I0)') 'WHAM Number of Equilibration MC Cycles: ', equilInterval
    write(nout, '(A)') '---------------------------'
  end if
  if (distCriteria) then
    write(nout, '(A)') 'Criteria Type: Distance'
  else if (minDistCriteria) then
    write(nout, '(A)') 'Criteria Type: Minimum Distance'
  else
    write(nout, '(A)') 'Criteria Type: Energy'
  end if
  write(nout, '(A, ES10.3)') 'Distance Criteria: ', Dist_Critr
  write(nout, '(A, ES10.3)') 'Distance Criteria Squared: ', Dist_Critr_sq
  write(nout, '(A, ES10.3)') 'AVBMC Volume: ', avbmc_vol
  write(nout, '(A)') 'Energy Criteria:'
  do i = 1, nMolTypes
    write(nout, '(*(ES10.3, 1X))') (Eng_Critr(i, j), j=1, nMolTypes)
  end do
  write(nout, '(A)') 'AVBMC E_Bias Alpha:'
  write(nout, '(A)') 'AVBMC E_Bias Alpha Insertion:'
  do i = 1, nMolTypes
    write(nout, '(*(ES10.3, 1X))') (AlphaTargetIns(i, j), j=1, nMolTypes)
  end do
  write(nout, '(A)') 'AVBMC E_Bias Alpha Deletion:'
  do i = 1, nMolTypes
    write(nout, '(*(ES10.3, 1X))') (AlphaTargetDel(i, j), j=1, nMolTypes)
  end do
  write(nout, '(A, ES10.3)') 'Temperature: ', temperature
  write(nout, '(A, ES10.3)') 'Gas Phase Density: ', gas_dens
  write(nout, '(A, ES10.3, 1X, A)') 'Initial Energy: ', E_T/outputEConv, trim(outputEngUnits)
  write(nout, '(A, ES10.3, 1X, A)') 'Initial Energy (Per Molecule): ', &
    E_T/outputEConv/real(NTotal, dp), trim(outputEngUnits)

  ! Step 11: Start simulation and output initial statistics
  ! - Prints simulation start message, initial analysis, and umbrella sampling data
  write(nout, '(A)') '------------------------------------------------'
  write(nout, '(A)') '         Simulation Start!'
  write(nout, '(A)') '------------------------------------------------'
  write(nout, '(A)') 'Cycle #  Particles  Energy  Acceptance Rates'
  if(useAnalysis) then
    call PostMoveAnalysis
  endif
  if (abs(E_T) < 1.0E6_dp) then
    write(nout, outFormat2) 0, NPART, E_T, (100.0_dp * movesAccepted(j) / movesAttempt(j), j=1, nMoveTypes)
  else
    write(nout, outFormat1) 0, NPART, E_T, (100.0_dp * movesAccepted(j) / movesAttempt(j), j=1, nMoveTypes)
  end if
  if (useUmbrella) then
    call UmbrellaHistAdd(E_T)
    call CheckInitialValues
  end if
  flush(nout)
  flush(35)
  call CPU_TIME(TimeStart) 
!--------------------------------------------------------------------------------------------------      
  ! Step 12: Execute main Monte Carlo loop
  ! - Performs nCycle iterations with nCycle2 moves, updating energy, histograms, and parameters
  do iCycle = 1, nCycle
    ! Step 12.1: Perform moves within cycle
    do iMove = 1, nCycle2
      ! Select move type based on probability
      ran_num = grnd()
      nSel = 1
      do while (moveProbability(nSel) < ran_num)
        nSel = nSel + 1
      enddo
      call mcMoveArray(nSel)%moveFunction(E_T, movesAccepted(nSel), movesAttempt(nSel))


      ! Step 12.2: Update pressure histogram if enabled
      if (calcPressure) then
        curIndx = NPart(nMolTypes)
        maxIndx = NMAX(nMolTypes) + 1
        do i = 1, nMolTypes - 1
          curIndx = curIndx + maxIndx * NPart(nMolTypes - i)
          maxIndx = maxIndx * (NMAX(nMolTypes - i) + 1)
        enddo
        NHist(curIndx) = NHist(curIndx) + 1.0_dp
        P_Avg(curIndx) = P_Avg(curIndx) + pressure
      endif

      ! Step 12.3: Perform post-move analysis if enabled
      if (useAnalysis) then
        call PostMoveAnalysis
      endif

      ! Step 12.4: Update WHAM or umbrella histogram
      if (useWham) then
        if (mod(iCycle, intervalWham) > equilInterval .and. useUmbrella) then
          call UmbrellaHistAdd(E_T)
        endif
      elseif (useUmbrella) then
        call UmbrellaHistAdd(E_T)
      endif
    enddo ! End move loop

    ! Step 12.5: Adjust move parameters every 100 cycles
    if (mod(iCycle, 100) == 0) then
      do i = 1, nMolTypes
        if (NMAX(i) <= 0) cycle
        call AdjustMax(acptTrans(i), atmpTrans(i), max_dist(i), dist_limit)
        call AdjustMax(acptRot(i), atmpRot(i), max_rot(i), rot_limit)
      enddo
    endif    

    ! Step 12.6: Output results at specified intervals
    if (mod(iCycle, outFreq_Screen) == 0) then
      if (abs(E_T) < 1.0e6_dp) then
        write(nout, outFormat2) iCycle, NPart, E_T, &
                                (100.0_dp * movesAccepted(j) / movesAttempt(j), j=1, nMoveTypes)
      else
        write(nout, outFormat1) iCycle, NPart, E_T, &
                                (100.0_dp * movesAccepted(j) / movesAttempt(j), j=1, nMoveTypes)
      endif
      if (useUmbrella) call ScreenOutputUmbrella
      flush(nout)
    endif
    if (mod(iCycle, outFreq_Traj) == 0) then
      call TrajOutput(iCycle, E_T)
    endif

    ! Step 12.7: Adjust WHAM bias if enabled
    if (useWham .and. mod(iCycle, intervalWham) == 0) then
      call WHAM_AdjustHist
    endif
  enddo ! End cycle loop  
!--------------------------------------------------------------------------------------------------

  ! Step 13: Record simulation end time
  ! - Captures final CPU time for performance reporting
  call CPU_TIME(TimeFinish)
  write(nout, *) "--------------------------------------------"

  ! Step 14: Allocate and store final data for validation
  ! - Allocates arrays for final neighbor list and energy tables, copies current state
  allocate(FinalNeighborList(maxMol, maxMol), stat=AllocateStatus)
  allocate(FinalETable(maxMol), stat=AllocateStatus)
  allocate(FinalNeiETable(maxMol), stat=AllocateStatus)
  if (AllocateStatus /= 0) stop "Allocation failed for final validation arrays"
  FinalNeighborList = NeighborList
  FinalETable = ETable
  FinalNeiETable = NeiETable
  E_Inter_Final = E_Inter_T
  E_NBond_Final = E_NBond_T
  E_Stretch_Final = E_Stretch_T
  E_Bend_Final = E_Bend_T
  E_Torsion_Final = E_Torsion_T

  ! Step 15: Compute final energy and pressure
  ! - Recomputes energy and pressure using detailed functions for validation
  call Detailed_ECalc(E_Final, errRtn)
  if (errRtn) write(nout, *) "Warning: Detailed energy calculation error"
  if (calcPressure) then
    call Detailed_PCalc(P_Final)
  endif

  ! Step 16: Output final trajectory
  ! - Writes final configuration to trajectory file and closes it
  call TrajOutput(iCycle, E_T)
  close(30)
  
  ! Step 17: Write energy comparison to debug file
  ! - Compares cumulative and detailed energy components for validation
  write(35, *) "---------------------------------"
  write(35, *) "Cumulative Energy vs. Detailed Final Energy"
  write(35, *) "Inter:", E_Inter_Final, E_Inter_T
  write(35, *) "Intra Nonbonded:", E_NBond_Final, E_NBond_T
  write(35, *) "Stretch:", E_Stretch_Final, E_Stretch_T
  write(35, *) "Bend:", E_Bend_Final, E_Bend_T
  write(35, *) "Torsional:", E_Torsion_Final, E_Torsion_T
  
  ! Step 18: Output final simulation report
  ! - Prints simulation time, final energy, cluster size, and energy/pressure discrepancies
  write(nout, *) "Simulation Time:", TimeFinish - TimeStart
  write(nout, *) "Final Energy:", E_Final / outputEConv, outputEngUnits
  write(nout, *) "Final Energy (Per Molecule):", E_Final / outputEConv / real(NTotal, dp), outputEngUnits
  write(nout, *) "Final Cluster Size:", NPart
  if (abs(E_Final) > 0.0_dp) then
    check = abs((E_Final - E_T) / E_Final)
  else
    check = abs(E_Final - E_T)
  endif
  if (check > 1.0e-6_dp) then
    write(nout, *) "=========================================="
    write(nout, *) "Energy Disagreement"
    write(nout, *) "Cumulative Energy:", E_T / outputEConv, outputEngUnits
    write(nout, *) "Cumulative Energy (Per Molecule):", E_T / outputEConv / real(NTotal, dp), outputEngUnits
    write(nout, *) "=========================================="
    write(35, *) "=========================================="
    write(35, *) "Energy Disagreement Error:"
    write(35, *) "Cumulative Energy:", E_T / outputEConv, outputEngUnits
    write(35, *) "Final Energy:", E_Final / outputEConv, outputEngUnits
    write(35, *) "=========================================="
  endif
  if (calcPressure) then
    write(nout, *) "Final Pressure:", P_Final
    write(nout, *) "Final Pressure (Per Molecule):", P_Final / real(NTotal, dp)
    if (abs(P_Final) > 0.0_dp) then
      check = abs((P_Final - pressure) / P_Final)
    else
      check = abs(P_Final - pressure)
    endif
    if (check > 1.0e-6_dp) then
      write(nout, *) "=========================================="
      write(nout, *) "Pressure Disagreement"
      write(nout, *) "Cumulative Pressure:", pressure
      write(nout, *) "Cumulative Pressure (Per Molecule):", pressure / real(NTotal, dp)
      write(nout, *) "=========================================="
      write(35, *) "=========================================="
      write(35, *) "Pressure Disagreement Error:"
      write(35, *) "Cumulative Pressure:", pressure / outputEConv, outputEngUnits
      write(35, *) "Final Pressure:", P_Final / outputEConv, outputEngUnits
      write(35, *) "=========================================="
    endif
  endif

  ! Step 19: Output move parameters and acceptance rates
  ! - Prints displacement, rotation, and acceptance statistics for all moves, including AVBMC/CBMC
  write(nout, *) "Final Max Displacement:", (max_dist(j), j=1, nMolTypes)
  write(nout, *) "Final Max Rotation:", (max_rot(j), j=1, nMolTypes)
  do i = 1, nMoveTypes
    if (movesAttempt(i) > 0.0_dp) then
      write(nout, "(1x,A,1x,A,A,F8.2)") "Acceptance Rate", trim(adjustl(moveName(i))), ": ", &
                                        100.0_dp * movesAccepted(i) / movesAttempt(i)
    endif
  enddo
  if (any(atmpTrans > 0.0_dp)) then
    write(nout, *) "Acceptance Translate (Mol Type):", (100.0_dp * acptTrans(j) / atmpTrans(j), j=1, nMolTypes)
  endif
  if (any(atmpRot > 0.0_dp)) then
    write(nout, *) "Acceptance Rotate (Mol Type):", (100.0_dp * acptRot(j) / atmpRot(j), j=1, nMolTypes)
  endif
  if (distGen_atmp > 0.0_dp) then
    write(nout, *) "Distance Generation Success Rate:", 100.0_dp * distGen_accpt / distGen_atmp
  endif
  if (angGen_atmp > 0.0_dp) then
    write(nout, *) "Angle Generation Success Rate:", 100.0_dp * angGen_accpt / angGen_atmp
  endif
  if (dihedGen_atmp > 0.0_dp) then
    write(nout, *) "Dihedral Angle Generation Success Rate:", 100.0_dp * dihedGen_accpt / dihedGen_atmp
  endif
  if (avbmcUsed) then
    write(nout, *) "********** AVBMC Insertion Rejection Breakdown ********"
    write(nout, *) "Percent Rejected due to Overlap:", 100.0_dp * ovrlapRej / totalRej
    write(nout, *) "Percent Rejected due to Detailed Balance:", 100.0_dp * dbalRej / totalRej
    write(nout, *) "Percent Rejected due to Cluster Criteria:", 100.0_dp * critriaRej / totalRej
    write(nout, *) "Percent Rejected due to Boundary Condition:", 100.0_dp * boundaryRej / totalRej
    write(nout, *) "********** AVBMC Deletion Rejection Breakdown ********"
    write(nout, *) "Percent Rejected due to Detailed Balance:", 100.0_dp * dbalRej_out / totalRej_out
    write(nout, *) "Percent Rejected due to Cluster Criteria:", 100.0_dp * critriaRej_out / totalRej_out
    write(nout, *) "Percent Rejected due to Boundary Condition:", 100.0_dp * boundaryRej_out / totalRej_out
    write(nout, *) "**********"
    write(nout, *) "Acceptance Swap In (Mol Type):", (100.0_dp * acptSwapIn(j) / atmpSwapIn(j), j=1, nMolTypes)
    write(nout, *) "Acceptance Swap Out (Mol Type):", (100.0_dp * acptSwapOut(j) / atmpSwapOut(j), j=1, nMolTypes)
  endif
  if (cbmcUsed .and. minDistCriteria) then
    do i = 1, nCBMCmax
      write(nout, *) "Acceptance CBMC (Mol Type) of ", i, " atoms:", &
                     (100.0_dp * acptCBMC(j, i) / atmpCBMC(j, i), j=1, nMolTypes)
    enddo
    do i = 1, nFECBMCmax
      write(nout, *) "Acceptance FECBMC (Mol Type) of ", i, " atoms:", &
                     (100.0_dp * acptFECBMC(j, i) / atmpFECBMC(j, i), j=1, nMolTypes)
    enddo
  endif

  ! Step 20: Collect and output histograms
  ! - Synchronizes processes, reduces pressure histograms, and collects other histograms
  write(nout, *) "Waiting for all processes to finish..."
  call MPI_BARRIER(MPI_COMM_WORLD, ierror)
  write(nout, *) "Writing output..."
  if (calcPressure) then
    allocate(NHist_Sum(pressLimit), P_Avg_Sum(pressLimit))
    call MPI_REDUCE(NHist, NHist_Sum, pressLimit, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    call MPI_REDUCE(P_Avg, P_Avg_Sum, pressLimit, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    if (myid == 0) then
      open(unit=50, file="ClusterPressure.dat")
      do i = 1, pressLimit
        if (NHist_Sum(i) > 0.0_dp) then
          write(50, *) i, P_Avg_Sum(i) / (3.0_dp * NHist_Sum(i))
        endif
      enddo
      close(50)
    endif
    deallocate(NHist_Sum, P_Avg_Sum)
  endif
  call CollectHistograms
  write(nout, *) "Histograms Condensed..."

  ! Step 21: Output final configuration and analysis
  ! - Writes final VMD-compatible configuration and analysis results (root process)
  if (myid == 0) then
    call Output_VMD_Final
    if (useAnalysis) call OutputAnalysis
  endif
  if (useUmbrella) call OutputUmbrellaHist
  write(nout, *) "Histograms Outputted..."

  ! Step 22: Validate energy and neighbor lists
  ! - Compares final energy table and neighbor list with cumulative values for errors
  write(35, *) "Energy Table:"
  do i = 1, maxMol
    if (.not. isActive(i)) cycle
    if (FinalETable(i) /= 0.0_dp) then
      if (abs((FinalETable(i) - ETable(i)) / FinalETable(i)) > 1.0e-6_dp) then
        write(35, *) i, ETable(i) / outputEConv, FinalETable(i) / outputEConv, "<---ETable Error"
      else
        write(35, *) i, ETable(i) / outputEConv
      endif
    else
      if (ETable(i) /= 0.0_dp) then
        write(35, *) i, ETable(i) / outputEConv, FinalETable(i) / outputEConv, "<---ETable Error"
      else
        write(35, *) i, ETable(i) / outputEConv
      endif
    endif
  enddo
  write(35, *) "NeighborList:"
  do i = 1, maxMol
    if (.not. isActive(i)) cycle
    do j = 1, maxMol
      if (.not. isActive(j) .or. i == j) cycle
      if (FinalNeighborList(i, j) .neqv. NeighborList(i, j)) then
        write(35, *) "Error", i, j, NeighborList(i, j), FinalNeighborList(i, j)
      else
        write(35, *) i, j, NeighborList(i, j)
      endif
    enddo
  enddo

  ! Step 23: Finalize WHAM and analytics
  ! - Outputs WHAM results and umbrella analytics if enabled
  if (useWham) call WHAM_Finalize
  if (energyAnalytics) call OutputUmbrellaAnalytics

  ! Step 24: Cleanup and finalize MPI
  ! - Closes output files, synchronizes processes, and terminates MPI environment
  call MPI_BARRIER(MPI_COMM_WORLD, ierror)
  write(nout, *) "Finished!"
  close(nout)
  close(35)
  call MPI_BARRIER(MPI_COMM_WORLD, ierror)
  call MPI_FINALIZE(ierror)
end program Nucleation
!========================================================    
!     This function adjusts the maximum displacement for both the rotation and translational moves such that
!     a 50% acceptance rate is maitained through out the simulation.        
      subroutine AdjustMax(acc_x, atmp_x, max_x, limit)
      use VarPrecision
      implicit none
      real(dp), intent(in) :: acc_x,atmp_x,limit
      real(dp), intent(inout):: max_x
      
      if(atmp_x .lt. 0.5E0_dp) then
        return
      endif

      if(acc_x/atmp_x .gt. 0.5E0_dp) then
        if(max_x*1.01E0_dp .lt. limit) then
          max_x = max_x * 1.01E0_dp
        else 
          max_x = limit       
        endif
      else
        max_x = max_x * 0.99E0_dp
      endif

 
      end subroutine
!===========================================================
! Writes a PSF file describing the topology of the maximum system (all NMAX molecules).
! Includes atom, bond, angle, and torsion information for VMD visualization.
subroutine WritePSF
  use SimParameters, only: nMolTypes, NMAX, MolName
  use Forcefield, only: nAtoms, atomData, atomArray, nBonds, bondArray, nAngles, bendArray, nTorsional, torsArray
  implicit none
  integer :: iType, iMol, iAtom        ! Loop indices: molecule type, molecule, atom
  integer :: totalAtoms, totalBonds     ! Total counts for atoms and bonds
  integer :: totalAngles, totalTorsions ! Total counts for angles and torsions
  integer :: atmIdx, bondIdx            ! Counters for atoms and bonds
  integer :: angleIdx, torsIdx          ! Counters for angles and torsions
  integer :: baseIdx, iBond, iAngle, iTors ! Base index and loop indices for bonds, angles, torsions
  character(len=8) :: segName, resName, atmName, atmSymb ! Fields for PSF (8 chars for XPLOR)
  integer :: ios                       ! I/O status for error checking

  ! Calculate total counts for atoms, bonds, angles, and torsions
  totalAtoms = sum(NMAX(1:nMolTypes) * nAtoms(1:nMolTypes))
  totalBonds = sum(NMAX(1:nMolTypes) * nBonds(1:nMolTypes))
  totalAngles = sum(NMAX(1:nMolTypes) * nAngles(1:nMolTypes))
  totalTorsions = sum(NMAX(1:nMolTypes) * nTorsional(1:nMolTypes))

  ! Open PSF file for writing with error checking
  open(unit=301, file='system.psf', status='unknown', form='formatted', iostat=ios)
  if (ios /= 0) then
    write(*, *) 'Error opening system.psf, iostat = ', ios
    stop
  endif

  ! Write PSF header (XPLOR format)
  write(301, '(A)') 'PSF EXT'
  write(301, '(A)') ''
  write(301, '(I8,1X,A)') 1, '!NTITLE'
  write(301, '(A)') 'REMARKS Grand Canonical Nucleation Simulation - Coarse-Grained Biomolecules'
  write(301, '(A)') ''

  ! Write atom section
  write(301, '(I8,1X,A)') totalAtoms, '!NATOM'
  atmIdx = 0
  do iType = 1, nMolTypes
    ! Pad or truncate segment name to 8 characters
    segName = adjustl(MolName(iType))
    if (len_trim(segName) > 8) segName = segName(1:8)
    if (len_trim(segName) < 8) segName = segName(1:len_trim(segName)) // repeat(' ', 8 - len_trim(segName))
    ! Use MolName as residue name (modify if a specific resName field exists)
    resName = segName ! Replace with atomData%resName if available
    do iMol = 1, NMAX(iType)
      do iAtom = 1, nAtoms(iType)
        atmIdx = atmIdx + 1
        ! Pad or truncate atom name and symbol
        atmName = adjustl(atomData(atomArray(iType,iAtom))%atmName)
        if (len_trim(atmName) > 8) atmName = atmName(1:8)
        if (len_trim(atmName) < 8) atmName = atmName(1:len_trim(atmName)) // repeat(' ', 8 - len_trim(atmName))
        atmSymb = adjustl(atomData(atomArray(iType,iAtom))%Symb)
        if (len_trim(atmSymb) > 8) atmSymb = atmSymb(1:8)
        if (len_trim(atmSymb) < 8) atmSymb = atmSymb(1:len_trim(atmSymb)) // repeat(' ', 8 - len_trim(atmSymb))
        ! Write atom line
        write(301, '(I10,1X,A8,1X,I8,1X,A8,1X,A8,1X,A8,1X,I4,1X,F10.6,1X,F8.4,1X,I8)') &
          atmIdx, segName, iMol, resName, atmName, atmSymb, &
          atomArray(iType,iAtom), &
          atomData(atomArray(iType,iAtom))%q, &
          atomData(atomArray(iType,iAtom))%mass, 0
      enddo
    enddo
  enddo

  ! Write bond section
  write(301, '(A)') ''
  write(301, '(I8,1X,A)') totalBonds, '!NBOND'
  bondIdx = 0
  do iType = 1, nMolTypes
    do iMol = 1, NMAX(iType)
      baseIdx = sum(NMAX(1:iType-1) * nAtoms(1:iType-1)) + (iMol - 1) * nAtoms(iType)
      do iBond = 1, nBonds(iType)
        bondIdx = bondIdx + 1
        write(301, '(2I10)', advance='no') &
          baseIdx + bondArray(iType,iBond)%bondMembr(1), &
          baseIdx + bondArray(iType,iBond)%bondMembr(2)
        if (mod(bondIdx, 4) == 0 .or. bondIdx == totalBonds) write(301, *)
      enddo
    enddo
  enddo

  ! Write angle section
  write(301, '(A)') ''
  write(301, '(I8,1X,A)') totalAngles, '!NTHETA'
  angleIdx = 0
  do iType = 1, nMolTypes
    do iMol = 1, NMAX(iType)
      baseIdx = sum(NMAX(1:iType-1) * nAtoms(1:iType-1)) + (iMol - 1) * nAtoms(iType)
      do iAngle = 1, nAngles(iType)
        angleIdx = angleIdx + 1
        write(301, '(3I10)', advance='no') &
          baseIdx + bendArray(iType,iAngle)%bendMembr(1), &
          baseIdx + bendArray(iType,iAngle)%bendMembr(2), &
          baseIdx + bendArray(iType,iAngle)%bendMembr(3)
        if (mod(angleIdx, 3) == 0 .or. angleIdx == totalAngles) write(301, *)
      enddo
    enddo
  enddo

  ! Write torsion section
  write(301, '(A)') ''
  write(301, '(I8,1X,A)') totalTorsions, '!NPHI'
  torsIdx = 0
  do iType = 1, nMolTypes
    do iMol = 1, NMAX(iType)
      baseIdx = sum(NMAX(1:iType-1) * nAtoms(1:iType-1)) + (iMol - 1) * nAtoms(iType)
      do iTors = 1, nTorsional(iType)
        torsIdx = torsIdx + 1
        write(301, '(4I10)', advance='no') &
          baseIdx + torsArray(iType,iTors)%torsMembr(1), &
          baseIdx + torsArray(iType,iTors)%torsMembr(2), &
          baseIdx + torsArray(iType,iTors)%torsMembr(3), &
          baseIdx + torsArray(iType,iTors)%torsMembr(4)
        if (mod(torsIdx, 2) == 0 .or. torsIdx == totalTorsions) write(301, *)
      enddo
    enddo
  enddo

  ! Close PSF file
  close(301)
end subroutine WritePSF
!===========================================================
      subroutine TrajOutput(iCycle, E_T)
      use VarPrecision
      use SimParameters
      use Coords
      use Forcefield
      implicit none
      integer(kind=8), intent(in) :: iCycle
      real(dp), intent(in) :: E_T
      integer :: iType, iMol, iAtom
      integer ::  atmType
      integer :: cnt      
      real(dp) :: xcm,ycm,zcm


      xcm = 0E0_dp
      ycm = 0E0_dp
      zcm = 0E0_dp
      cnt = 0
      do iType = 1,nMolTypes
        do iMol = 1,NPART(iType)
          xcm = xcm + MolArray(iType)%mol(iMol)%x(1)
          ycm = ycm + MolArray(iType)%mol(iMol)%y(1)
          zcm = zcm + MolArray(iType)%mol(iMol)%z(1)
          cnt = cnt + 1
        enddo
      enddo
      
      xcm = xcm/real(cnt, dp)
      ycm = ycm/real(cnt, dp)
      zcm = zcm/real(cnt, dp)

      write(30,*) vmdAtoms
      write(30,*) NPART, E_T
      do iType = 1,nMolTypes
        do iMol = 1, NPART(iType)
          do iAtom = 1, nAtoms(iType)
            atmType = atomArray(iType,iAtom)
            write(30,*) atomData(atmType)%Symb, &
                        MolArray(iType)%mol(iMol)%x(iAtom)-xcm, &
                        MolArray(iType)%mol(iMol)%y(iAtom)-ycm, &
                        MolArray(iType)%mol(iMol)%z(iAtom)-zcm
          enddo
        enddo
        do iMol = NPART(iType)+1, NMAX(iType)
          do iAtom = 1, nAtoms(iType)
            atmType = atomArray(iType,iAtom)
            write(30,*) atomData(atmType)%Symb, 1E7_dp, 1E7_dp, 1E7_dp
          enddo
        enddo
      enddo


      end subroutine
!=============================================================================================
!     This subroutine prints a dummy frame to the Trajectory file.  This
!     is done so that when the trajectory is loaded into VMD it will properly
!     display the cluster configuration.
      subroutine InitialTrajOutput
      use VarPrecision
      use SimParameters      
      use Coords
      use Forcefield
      implicit none
      real(dp), parameter :: xOffset = 50E0_dp
      real(dp), parameter :: yOffset = 20E0_dp
      integer :: iType, iMol, iAtom
      integer :: atmType
      integer :: cnt      
      real(dp) :: xcm,ycm,zcm


      xcm = 0E0_dp
      ycm = 0E0_dp
      zcm = 0E0_dp
      cnt = 0
      do iType = 1,nMolTypes
        do iMol = 1,NMAX(iType)
          do iAtom = 1, nAtoms(iType)
            xcm = xcm + gasConfig(iType)%x(iAtom) + xOffset*( real(NMAX(iType),dp)/2E0_dp - real(iMol,dp) )
            ycm = ycm + gasConfig(iType)%y(iAtom) + yOffset*real(iType, dp)
            zcm = zcm + gasConfig(iType)%z(iAtom)
            cnt = cnt + 1
          enddo          
        enddo
      enddo

      xcm = xcm/real(cnt, dp)
      ycm = ycm/real(cnt, dp)
      zcm = zcm/real(cnt, dp)

      write(30,*) vmdAtoms
      write(30,*) NMAX
      do iType = 1,nMolTypes
        do iMol = 1, NMAX(iType)
          do iAtom = 1, nAtoms(iType)
            atmType = atomArray(iType,iAtom)
            write(30,*) atomData(atmType)%Symb, &
                      gasConfig(iType)%x(iAtom) + xOffset*( real(NMAX(iType),dp)/2E0_dp - real(iMol, dp) )-xcm, &
                      gasConfig(iType)%y(iAtom) + yOffset*real(iType, dp) - ycm, &
                      gasConfig(iType)%z(iAtom) - zcm
          enddo
        enddo
      enddo


      end subroutine
!========================================================            
! Subroutine: LowerCaseLine
! Purpose   : Converts all uppercase characters in a string to lowercase.
! Notes     : Uses ASCII values to check and convert characters.
!==================================================================
subroutine LowerCaseLine(line)
  implicit none
  character(len=*), intent(inout) :: line
  integer, parameter :: offset = ichar("a") - ichar("A")
  integer :: i, sizeLine
  integer :: curVal, newVal
  sizeLine = len(line)

  ! Loop through each character and convert uppercase to lowercase
  do i = 1, sizeLine
    curVal = ichar(line(i:i))
    if (curVal .le. ichar("Z")) then
      if (curVal .ge. ichar("A")) then
        newVal = curVal + offset
        line(i:i) = char(newVal)
      end if
    end if
  end do

end subroutine LowerCaseLine
!=============================================================================================      

