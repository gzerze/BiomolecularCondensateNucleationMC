!========================================================            
module ScriptInput
  integer, parameter :: maxLineLen = 500   
  contains
!========================================================            
!==============================================================
! Subroutine: Script_ReadParameters
! Purpose   : Read and parse user-provided simulation parameters
!             from a script/input file and initialize relevant
!             variables across modules.
! Author:
!   Aliasghar Sepehri
!==============================================================        
subroutine Script_ReadParameters(seed, screenEcho)
  use AnalysisMain, only: ScriptAnalysisInput
  use CoordinateFunctions, only: ReadInitialConfiguration, RecenterCoordinates, &
                                 ReadInitialGasPhase
  use ForceFieldInput, only: SetForcefieldType, ScriptForcefield, fieldTypeSet
  use MoveTypeModule, only: ScriptInput_MCMove
  use SimParameters, only: nMolTypes, Eng_Critr, temperature, ncycle
  use UmbrellaSamplingNew, only: useUmbrella, ScriptInput_Umbrella
  use WHAM_Functions, only: WHAM_Initialize, intervalWHAM, nWhamItter, useWHAM
  use EnergyTables, only: AlphaTargetIns, AlphaTargetDel
  use ParallelVar, only: nout
  implicit none
  !--------------------------------------------------------------
  ! Output arguments
  !--------------------------------------------------------------
  logical, intent(OUT)  :: screenEcho   ! Whether to print parsed values to screen
  integer, intent(OUT)  :: seed         ! Random number seed

  !--------------------------------------------------------------
  ! Internal variables
  !--------------------------------------------------------------
  integer :: i, ii, j                   ! Loop counters
  integer :: nArgs                      ! Number of command-line arguments
  integer :: iLine, lineStat            ! Line parsing and status flags
  integer :: nLines, nForceLines        ! Line counts for general and force field input
  integer :: lineBuffer, AllocateStat   ! Buffer size and allocation status

  !--------------------------------------------------------------
  ! Line storage arrays
  !--------------------------------------------------------------
  integer, allocatable :: lineNumber(:)       ! Line numbers for main input
  integer, allocatable :: ffLineNumber(:)     ! Line numbers for force field input
  character(len=maxLineLen), allocatable :: lineStore(:)        ! Input script lines
  character(len=maxLineLen), allocatable :: forcefieldStore(:)  ! Force field script lines

  !--------------------------------------------------------------
  ! Command parsing and file naming
  !--------------------------------------------------------------
  character(len=25) :: command, command2, dummy   ! Parsed command tokens
  character(len=50) :: fileName                   ! General input file name
  character(len=50) :: forcefieldFile             ! Force field input file name
  !=============================================================
  ! Parse command-line arguments to determine input file name
  !=============================================================
  nArgs = command_argument_count()

  if (nArgs > 1) then
    ! Too many arguments: terminate with an error message
    stop "This program only takes one argument"

  elseif (nArgs == 1) then
    ! One argument: treat it as the input script file name
    call get_command_argument(1, fileName)
    call LoadFile(lineStore, nLines, lineNumber, fileName)

  elseif (nArgs == 0) then
    ! No input file provided: terminate with an error message
    write(nout,*) "ERROR! No Input File has been given!"
    stop
  endif
!===============================================================
! Set default values for seed and screenEcho before parsing input
!===============================================================
call setDefaults(seed, screenEcho)

!===============================================================
! Begin processing each line from the input script
! - Skip lines buffered by multi-line commands (e.g., move definitions)
! - Skip blank or commented lines
! - Convert command to lowercase for case-insensitive matching
!===============================================================
  lineBuffer = 0
  do iLine = 1, nLines

    ! Skip lines if they are buffered (used by multi-line commands)
    if (lineBuffer > 0) then
      lineBuffer = lineBuffer - 1
      cycle
    end if

    ! Initialize line status flag
    lineStat = 0

    ! Extract first word/command from line and convert to lowercase
    call getCommand(lineStore(iLine), command, lineStat)
    call LowerCaseLine(command)

    ! Skip line if it's a comment or empty
    if (lineStat == 1) then
      cycle
    end if


    select case(trim(adjustl(command)))

!==============================
! Case: analysis block
!==============================
      case("analysis")
        call FindCommandBlock(iLine, lineStore, "end_analysis", lineBuffer)
        call ScriptAnalysisInput(lineStore(iLine:iLine+lineBuffer))

!==============================
! Case: set individual variable
!==============================
      case("set")
        call setVariable(lineStore(iLine), seed, screenEcho, lineStat)
        if (lineStat .eq. -1) then
            write(nout,"(A,2x,I10)") "ERROR! Unknown Variable Name on Line", lineNumber(iLine)
            write(nout,*) lineStore(iLine)
            stop
        endif

!==============================
! Case: movetypes block
!==============================
      case("movetypes")
        call FindCommandBlock(iLine, lineStore, "end_movetypes", lineBuffer)
        call ScriptInput_MCMove(lineStore(iLine:iLine+lineBuffer))

!==============================
! Case: umbrella sampling block
!==============================
      case("umbrella")
        call FindCommandBlock(iLine, lineStore, "end_umbrella", lineBuffer)
        call ScriptInput_Umbrella(lineStore(iLine:iLine+lineBuffer))

!==============================
! Case: specify force field type
!==============================
      case("forcefieldtype")
        read(lineStore(iLine),*) dummy, command2
        call LowerCaseLine(command2)
        call SetForcefieldType(command2)
        fieldTypeSet = .true.

!==============================
! Case: external force field file
!==============================
      case("forcefieldfile")
        read(lineStore(iLine),*) dummy, command2
        forcefieldFile = trim(adjustl(command2))
        call LoadFile(forcefieldStore, nForceLines, ffLineNumber, forcefieldFile)
        call ScriptForcefield(forcefieldStore)
        if (allocated(forcefieldStore)) then
            deallocate(forcefieldStore)
        endif

!==============================
! Case: boundary energy criteria matrix
!==============================
      case("boundary")
        if (nMolTypes .eq. 0) then
            write(nout,*) "INPUT ERROR! Boundary is called before the number of molecule types has been assigned"
            stop
        endif
        call FindCommandBlock(iLine, lineStore, "end_boundary", lineBuffer)
        allocate(Eng_Critr(1:nMolTypes,1:nMolTypes), stat = AllocateStat)
        ii = 0
        do i = iLine+1, iLine+lineBuffer-1
            ii = ii + 1
            read(lineStore(i),*) (Eng_Critr(ii,j), j = 1, nMolTypes)
        enddo

!==============================
! Case: AlphaTarget matrices
!==============================
      case("alphatargetins")
        if (nMolTypes .eq. 0) then
            write(nout,*) "INPUT ERROR! AlphaTargetIns is called before the number of molecule types has been assigned"
            stop
        endif
        call FindCommandBlock(iLine, lineStore, "end_alphatargetins", lineBuffer)
        if (lineBuffer - 1 .gt. nMolTypes) then
            write(nout,*) "ERROR! User has given the AlphaTargetIns command far too many lines."
            write(nout,*) "AlphaTargetIns takes a number of lines equal to the number of molecule types."
            stop
        endif
        ii = 0
        do i = iLine+1, iLine+lineBuffer-1
            ii = ii + 1
            read(lineStore(i),*) (AlphaTargetIns(ii,j), j = 1, nMolTypes)
        enddo

      case("alphatargetdel")
        if (nMolTypes .eq. 0) then
            write(nout,*) "INPUT ERROR! AlphaTargetDel is called before the number of molecule types has been assigned"
            stop
        endif
        call FindCommandBlock(iLine, lineStore, "end_alphatargetdel", lineBuffer)
        if (lineBuffer - 1 .gt. nMolTypes) then
            write(nout,*) "ERROR! User has given the AlphaTargetDel command far too many lines."
            write(nout,*) "AlphaTargetDel takes a number of lines equal to the number of molecule types."
            stop
        endif
        ii = 0
        do i = iLine+1, iLine+lineBuffer-1
            ii = ii + 1
            read(lineStore(i),*) (AlphaTargetDel(ii,j), j = 1, nMolTypes)
        enddo

!==============================
! Case: unknown/unrecognized command
!==============================
      case default
        write(nout,"(A,2x,I10)") "ERROR! Unknown Command on Line", lineNumber(iLine)
        write(nout,*) lineStore(iLine)
        stop

    end select
  enddo

  call ReadInitialConfiguration

  call RecenterCoordinates

  call ReadInitialGasPhase     

!--------------------------------------------------------------
! Initialize WHAM (Weighted Histogram Analysis Method) setup
!--------------------------------------------------------------
  if(useWHAM) then

  ! Ensure umbrella sampling is enabled when using WHAM
    if(.not. useUmbrella) then
      write(nout,*) "ERROR! The WHAM method cannot be used if no Umbrella sampling variables are given!"
      stop
    endif

  ! Calculate the total number of WHAM iterations to perform
  ! Based on total number of MC cycles and the interval at which WHAM is updated
    nWhamItter = ceiling(dble(ncycle) / dble(intervalWHAM))

  ! Initialize arrays and parameters for WHAM calculations
    call WHAM_Initialize

  endif

end subroutine Script_ReadParameters
!========================================================     
! Subroutine: setVariable
! Purpose  : Parses and processes a single "set" line from the input file.
!            Assigns simulation parameters based on recognized keywords.
!            Outputs all recognized settings to debug unit 34.
! Author:
!   Aliasghar Sepehri
!==========================================================================      
subroutine setVariable(line, seed, screenEcho, lineStat)
  use VarPrecision, only: dp
  use SimParameters, only: calcPressure, ncycle, ncycle2, nMolTypes, NPart, NPart_New, NMin, NMax, maxMol, &
       gas_dens, MolName, distCriteria, minDistCriteria, Dist_Critr, Dist_Critr_sq, &
       temperature, beta, softCutoff, outFreq_Screen, outFreq_Traj, multipleInput, outputEngUnits, &
       outputEConv, outputLenUnits, outputLenConv, engDefault, convEng, distDefault, convDist, &
       angDefaults, convAng, DielecConst, DebyeLength, kapa
  use CBMC_Variables, only: nRosenTrials, maxRosenTrial, nSwapInTrials, maxSwapInTrial, nSwapOutTrials, &
       maxSwapOutTrial
  use EnergyTables, only: E_gasIntra, AlphaTargetIns, AlphaTargetDel
  use Coords, only: rosenTrial
  use WHAM_Module, only: useWHAM, intervalWHAM, maxSelfConsist, whamEstInterval, equilInterval, tolLimit
  use AcceptRates, only: nCBMCmax, nFECBMCmax
  use Units, only: FindEngUnit, FindLengthUnit, FindAngularUnit
  use PairStorage, only: useDistStore
  implicit none

  character(len=maxLineLen), intent(in) :: line      ! input line from the script
  logical, intent(out) :: screenEcho                ! toggles screen output
  integer, intent(out) :: seed, lineStat            ! random seed, status flag

  character(len=30) :: dummy, command, format_string   ! placeholders for parsed tokens
  logical :: logicValue                             ! logical value parsed
  integer :: j                                      ! loop index
  integer :: intValue, AllocateStat                 ! parsed int and alloc status
  real(dp) :: realValue                             ! parsed real(dp) value
      
  lineStat  = 0

  ! Extract command name and lowercase it
  read(line,*) dummy, command
  call LowerCaseLine(command)
  write(34,*) 'Processing set command:', trim(command)

  !=======================================================================
  ! Process recognized commands
  !=======================================================================
  select case(trim(adjustl(command)))

    ! Set AVBMC radius distance and its square
    case("avbmc_distance")
      read(line,*) dummy, command, realValue
      Dist_Critr = realValue
      Dist_Critr_sq = realValue**2
      write(34,*) '  avbmc_distance =', realValue

    ! Toggle pressure calculation
    case("calcpressure")
      read(line,*) dummy, command, logicValue
      calcPressure = logicValue
      write(34,*) '  calcPressure =', logicValue

    ! Number of MC cycles
    case("cycles")
      read(line,*) dummy, command, realValue
      nCycle = nint(realValue)
      write(34,*) '  nCycle =', nCycle

    ! Use distance criteria for clustering
    case("distancecriteria")
      read(line,*) dummy, command, logicValue
      distCriteria = logicValue
      write(34,*) '  distCriteria =', logicValue

    ! Use minimum distance criteria for clustering
    case("minimumdistancecriteria")
      read(line,*) dummy, command, logicValue
      minDistCriteria = logicValue
      write(34,*) '  minDistCriteria =', logicValue

    ! Toggle distance storage
    case("distancestorage")
      read(line,*) dummy, command, logicValue
      useDistStore = logicValue
      write(34,*) '  useDistStore =', logicValue

    ! Input dilute densities and convert units
    case("gasdensity")
      if(.not. allocated(gas_dens)) then
        write(*,*) "INPUT ERROR! GasDensity is called before molecule types assigned"
        stop
      endif
      read(line,*) dummy, command, (gas_dens(j), j=1, nMolTypes)
      do j = 1, nMolTypes
        gas_dens(j) = gas_dens(j) / (convDist**3)
      enddo
      write(34,*) '  gas_dens assigned for all molecule types'

    ! Number of MC moves in one cycle
    case("moves")
      read(line,*) dummy, command, realValue
      ncycle2 = nint(realValue)
      write(34,*) '  ncycle2 =', ncycle2

    ! Number of molecule types and allocate arrays
    case("moleculetypes")
      read(line,*) dummy, command, realValue
      nMolTypes = nint(realValue)
      write(34,*) '  nMolTypes =', nMolTypes
      allocate(NPART(1:nMolTypes), STAT=AllocateStat)
      allocate(NPART_new(1:nMolTypes), STAT=AllocateStat)
      allocate(NMIN(1:nMolTypes), STAT=AllocateStat)
      allocate(NMAX(1:nMolTypes), STAT=AllocateStat)
      allocate(gas_dens(1:nMolTypes), STAT=AllocateStat)
      allocate(E_gasIntra(1:nMolTypes), STAT=AllocateStat)
      allocate(nRosenTrials(1:nMolTypes), STAT=AllocateStat)
      allocate(nSwapInTrials(1:nMolTypes), STAT=AllocateStat)
      allocate(nSwapOutTrials(1:nMolTypes), STAT=AllocateStat)
      allocate(MolName(1:nMolTypes), STAT=AllocateStat)
      NMIN = 0
      NMAX = 0
      gas_dens = 0E0_dp
      E_gasIntra = 0E0_dp
      nRosenTrials = 1
      nSwapInTrials = 1
      nSwapOutTrials = 1
      allocate(AlphaTargetIns(1:nMolTypes,1:nMolTypes), STAT=AllocateStat)
      AlphaTargetIns = 0E0_dp
      allocate(AlphaTargetDel(1:nMolTypes,1:nMolTypes), STAT=AllocateStat)
      AlphaTargetDel = 0E0_dp
      format_string = "(A,I1)"
      do j = 1, nMolTypes
        write(MolName(j), format_string) "mol", j
        MolName(j) = trim(adjustl(MolName(j)))
      enddo

    ! Molecule Name per type
    case("moleculename")
      if(.not. allocated(MolName)) then
        write(*,*) "INPUT ERROR! MolName called before molecule types assigned"
        stop
      endif
      read(line,*) dummy, command, (MolName(j), j=1, nMolTypes)
      do j = 1, nMolTypes
        MolName(j) = trim(adjustl(MolName(j)))
      enddo
      write(34,*) '  MolName =', (MolName(j), j=1, nMolTypes)

    ! Minimum molecules per type
    case("molmin")
      if(.not. allocated(NMIN)) then
        write(*,*) "INPUT ERROR! molmin called before molecule types assigned"
        stop
      endif
      read(line,*) dummy, command, (NMIN(j), j=1, nMolTypes)
      write(34,*) '  NMIN =', (NMIN(j), j=1, nMolTypes)

    ! Maximum molecules per type
    case("molmax")
      if(.not. allocated(NMAX)) then
        write(*,*) "INPUT ERROR! molmax called before molecule types assigned"
        stop
      endif
      read(line,*) dummy, command, (NMAX(j), j=1, nMolTypes)
      maxMol = sum(NMAX)
      write(34,*) '  NMAX =', (NMAX(j), j=1, nMolTypes)

    ! Number of trials for CBMC and FECBMC moves
    case("rosentrials")        
      if(.not. allocated(nRosenTrials)) then
        write(*,*) "INPUT ERROR! molmax is called before the number of molecule types has been assigned"
        stop
      endif
      read(line,*) dummy, command, (nRosenTrials(j), j=1, nMolTypes)
      maxRosenTrial = maxval(nRosenTrials)
      allocate(rosenTrial(1:maxRosenTrial))

    ! Number of trials for Insertion move
    case("swapintrials")        
      if(.not. allocated(nSwapInTrials)) then
        write(*,*) "INPUT ERROR! molmax is called before the number of molecule types has been assigned"
        stop
      endif
      read(line,*) dummy, command, (nSwapInTrials(j), j=1, nMolTypes)
      maxSwapInTrial = maxval(nSwapInTrials)
      write(34,*) '  TrilasIn =', (nSwapInTrials(j), j=1, nMolTypes)

    ! Number of trials for Deletion move
    case("swapouttrials")        
      if(.not. allocated(nSwapOutTrials)) then
        write(*,*) "INPUT ERROR! molmax is called before the number of molecule types has been assigned"
        stop
      endif
      read(line,*) dummy, command, (nSwapOutTrials(j), j=1, nMolTypes)
      maxSwapOutTrial = maxval(nSwapOutTrials)
      write(34,*) '  TrilasOut =', (nSwapInTrials(j), j=1, nMolTypes)

    ! Set temperature and compute beta
    case("temperature")
      read(line,*) dummy, command, realValue
      temperature = realValue
      beta = 1E0_dp / temperature
      write(34,*) '      temperature =', temperature, ', beta =', beta

    ! Toggle screen echo
    case("screenecho")
      read(line,*) dummy, command, logicValue
      screenEcho = logicValue
      write(34,*) '      screenEcho =', screenEcho

    ! Set random number generator seed
    case("rng_seed")
      read(line,*) dummy, command, intValue
      seed = intValue
      if (seed .lt. 0) then
        call system_clock(intValue)
        seed = mod(intValue, 10000)
      endif
      write(34,*) '      rng_seed =', seed

    ! Soft cutoff energy to reject a move immediately 
    case("softcutoff")
      read(line,*) dummy, command, realValue
      softCutoff = realValue
      write(34,*) '      softCutoff =', softCutoff

    ! Output frequency for screen
    case("screen_outfreq")
      read(line,*) dummy, command, intValue
      outFreq_Screen = intValue
      write(34,*) '      outFreq_Screen =', outFreq_Screen

    ! Output frequency for trajectory
    case("trajectory_outfreq")
      read(line,*) dummy, command, intValue
      outFreq_Traj = intValue
      write(34,*) '      outFreq_Traj =', outFreq_Traj

    ! Toggle multiple input config mode
    case("multipleinputconfig")
      read(line,*) dummy, command, logicValue
      multipleInput = logicValue
      write(34,*) '      multipleInput =', multipleInput

    ! Set energy output units
    case("out_energyunits")
      read(line,*) dummy, command, outputEngUnits
      outputEConv = FindEngUnit(outputEngUnits)
      write(34,*) '      outputEngUnits =', outputEngUnits

    ! Set distance output units
    case("out_distunits")
      read(line,*) dummy, command, outputLenUnits
      outputLenConv = FindLengthUnit(outputLenUnits)
      write(34,*) '      outputLenUnits =', outputLenUnits

    ! Set default energy unit
    case("defulat_energyunits")
      read(line,*) dummy, command, engDefault
      convEng = FindEngUnit(engDefault)
      write(34,*) '      engDefault =', engDefault

    ! Set default distance unit
    case("defulat_distunits")
      read(line,*) dummy, command, distDefault
      convDist = FindLengthUnit(distDefault)
      write(34,*) '      distDefault =', distDefault

    ! Set default angle unit
    case("defulat_angleunits")
      read(line,*) dummy, command, angDefaults
      convAng = FindAngularUnit(angDefaults)
      write(34,*) '      angDefaults =', angDefaults

    ! Toggle WHAM usage
    case("usewham")
      read(line,*) dummy, command, logicValue
      useWHAM = logicValue
      write(34,*) '      useWHAM =', useWHAM

    ! WHAM segment length
    case("whamseglength")
      read(line,*) dummy, command, intValue
      intervalWHAM = intValue
      write(34,*) '      intervalWHAM =', intervalWHAM

    ! WHAM maximum iterations
    case("whammaxiteration")
      read(line,*) dummy, command, intValue
      maxSelfConsist = intValue
      write(34,*) '      maxSelfConsist =', maxSelfConsist

    ! WHAM delta G replace interval
    case("whamdgreplace")
      read(line,*) dummy, command, intValue
      whamEstInterval = intValue
      write(34,*) '      whamEstInterval =', whamEstInterval

    ! WHAM equilibration cycles
    case("whamequilcycle")
      read(line,*) dummy, command, intValue
      equilInterval = intValue
      write(34,*) '      equilInterval =', equilInterval

    ! WHAM tolerance limit
    case("whamtol")
      read(line,*) dummy, command, realValue
      tolLimit = realValue
      write(34,*) '      tolLimit =', tolLimit

    ! CBMC maximum growing segments
    case("maxcbmcunits")
      read(line,*) dummy, command, intValue
      nCBMCmax = intValue
      write(34,*) '      nCBMCmax =', nCBMCmax

    ! FECBMC maximum growing segments
    case("maxfecbmcunits")
      read(line,*) dummy, command, intValue
      nFECBMCmax = intValue
      write(34,*) '      nFECBMCmax =', nFECBMCmax

    ! Debye length and inverse
    case("debyelength")
      read(line,*) dummy, command, realValue
      DebyeLength = realValue * convDist
      kapa = 1E0_dp / DebyeLength
      write(34,*) '      DebyeLength =', DebyeLength, ', kapa =', kapa

    ! Dielectric constant
    case("dielectricconstant")
      read(line,*) dummy, command, realValue
      DielecConst = realValue
      write(34,*) '      DielecConst =', DielecConst

    ! Catch unknown commands
    case default
      lineStat = -1
  end select

end subroutine setVariable
!========================================================     
! Subroutine: setDefaults
! Purpose:
!   Initializes all major simulation parameters to their default
!   values before reading user-defined settings from input script.
! Author:
!   Aliasghar Sepehri
!===============================================================
subroutine setDefaults(seed, screenEcho)
  use VarPrecision, only: dp
  use SimParameters, only: ncycle, ncycle2, maxMol, Dist_Critr, Dist_Critr_sq, temperature, beta, &
       DielecConst, softCutoff, outFreq_Screen, outFreq_Traj, multipleInput, outputEConv, &
       outputLenConv, convEng, convDist, convAng
  use CBMC_Variables, only: maxRosenTrial
  use WHAM_Module, only: useWHAM, intervalWHAM, maxSelfConsist, whamEstInterval, equilInterval, tolLimit
  use Units, only: FindEngUnit, FindLengthUnit, FindAngularUnit
  implicit none

  !==========================
  ! Output arguments
  !==========================
  logical, intent(out) :: screenEcho
  integer, intent(out) :: seed

  !==========================
  ! Random seed
  !==========================
  seed = -1

  !==========================
  ! General simulation parameters
  !==========================
  Dist_Critr = 0E0_dp         ! Cluster distance criterion
  Dist_Critr_sq = 0E0_dp
  nCycle = 0                  ! Number of MC cycles
  ncycle2 = 0
  maxMol = 0                  ! Max number of molecules

  maxRosenTrial = 1           ! CBMC trial count

  temperature = 0E0_dp
  beta = 0E0_dp               ! 1/(kB*T)
  DielecConst = 1E0_dp        ! Dielectric constant

  screenEcho = .true.         ! Echo inputs to screen
  softCutoff = 100E0_dp       ! Soft interaction cutoff

  outFreq_Screen = 1000       ! Screen output frequency
  outFreq_Traj   = 1000       ! Trajectory output frequency
  multipleInput  = .false.    ! Multiple input configuration definitions

  !==========================
  ! Unit conversion factors
  !==========================
  outputEConv   = FindEngUnit("kb")
  outputLenConv = FindLengthUnit("ang")
  convEng       = FindEngUnit("kb")
  convDist      = FindLengthUnit("ang")
  convAng       = FindAngularUnit("deg")

  !==========================
  ! WHAM (Weighted Histogram Analysis Method) options
  !==========================
  useWHAM           = .false.
  intervalWHAM      = 10000
  maxSelfConsist    = 100
  whamEstInterval   = -1
  equilInterval     = 10
  tolLimit          = 1E-5_dp    ! WHAM convergence tolerance

end subroutine setDefaults
!==============================================================
! Subroutine: LoadFile
! Purpose   : Reads a text file line by line into an array while
!             filtering out comments and blank lines. Also stores
!             corresponding line numbers from the original file.
! Author:
!   Aliasghar Sepehri
!========================================================            
subroutine LoadFile(lineArray, nLines, lineNumber, fileName)
  use ParallelVar, only: myid                       ! MPI rank
  use SimParameters, only: echoInput                ! Verbosity toggle
  implicit none

  !===========================
  ! Input / Output Variables
  !===========================
  character(len=maxLineLen), allocatable, intent(inout) :: lineArray(:)
  integer, allocatable, intent(inout) :: lineNumber(:)
  character(len=50), intent(in) :: fileName

  !===========================
  ! Local Variables
  !===========================
  character(len=maxLineLen), allocatable :: rawLines(:)   ! Temporary raw lines
  character(len=25) :: command                            ! Buffer for first token
  integer :: nLines                ! Number of parsed (valid) lines
  integer :: nRawLines             ! Total lines read from file
  integer :: i, iLine              ! Loop indices
  integer :: lineStat              ! Status for reading lines
  integer :: AllocateStat          ! Status of allocation
  integer :: InOutStat             ! I/O status flag
  !=============================================================
  ! Attempt to open the input file containing script commands
  !=============================================================
  open(unit=54, file=trim(adjustl(fileName)), status='OLD', iostat=InOutStat)

  ! If file open fails, print an error message and abort the run
  if (InOutStat .gt. 0) then
    if (myid .eq. 0) then
      write(*,*) "ERROR: The file specified could not be opened."
      write(*,*) "Please check if the file exists and the name is correct."
      write(*,*) "Attempted to open file: ", trim(adjustl(fileName))
    end if
    stop
  end if

  !===============================================================
  ! Count the number of lines in the input file to allocate storage
  !===============================================================
  nRawLines = 0
  do iLine = 1, nint(1d7)
    read(54, *, iostat = lineStat)
    if (lineStat < 0) exit
    nRawLines = nRawLines + 1
  end do
  rewind(54)

  !===============================================================
  ! Read in the input file line by line and optionally echo content
  !===============================================================
  allocate(rawLines(1:nRawLines), stat = AllocateStat)
  do iLine = 1, nRawLines
    read(54, "(A)") rawLines(iLine)
    if (echoInput) then
      write(35, *) rawLines(iLine)
    end if
  end do
  close(54)

  !===============================================================
  ! Identify valid (non-commented) lines for further processing
  !===============================================================
  nLines = 0
  do iLine = 1, nRawLines
    lineStat = 0
    call getCommand(rawLines(iLine), command, lineStat)
    if (lineStat == 0) then
      nLines = nLines + 1
    end if
  end do

  !===============================================================
  ! Store the valid lines and their original line numbers
  !===============================================================
  allocate(lineArray(1:nLines))
  allocate(lineNumber(1:nLines))
  i = 0
  do iLine = 1, nRawLines
    lineStat = 0
    call getCommand(rawLines(iLine), command, lineStat)
    if (lineStat == 0) then
      i = i + 1
      lineArray(i) = rawLines(iLine)
      lineNumber(i) = iLine
    end if
  end do

end subroutine LoadFile
!========================================================            
! Subroutine: GetCommand
! Purpose   : Extracts the first keyword/command from a line of input
!             and identifies whether the line is a comment or empty
! Arguments :
!   - line      : Full line of input text
!   - command   : Extracted command keyword (first word)
!   - lineStat  : Status flag (0 = valid line, 1 = comment or empty)
! Author:
!   Aliasghar Sepehri
!===============================================================
subroutine GetCommand(line, command, lineStat)
  use VarPrecision
  implicit none

  character(len=*), intent(in)  :: line
  character(len=25), intent(out) :: command
  integer, intent(out) :: lineStat

  integer :: i, sizeLine, lowerLim, upperLim

  sizeLine = len(line)
  lineStat = 0
  i = 1

  !===============================================================
  ! Skip leading spaces and detect if line is a comment
  !===============================================================
  do while (i <= sizeLine)
    if (ichar(line(i:i)) /= ichar(' ')) then
      if (ichar(line(i:i)) == ichar('#')) then
        lineStat = 1    ! Comment line
        return
      else
        exit            ! Found start of command
      end if
    end if
    i = i + 1
  end do

  !===============================================================
  ! If we reached the end, the line is empty
  !===============================================================
  if (i >= sizeLine) then
    lineStat = 1        ! Empty line
    return
  end if
  lowerLim = i

  !===============================================================
  ! Find the end of the first word (command)
  !===============================================================
  do while (i <= sizeLine)
    if (line(i:i) == ' ') exit
    i = i + 1
  end do
  upperLim = i - 1

  !===============================================================
  ! Extract the command from the line
  !===============================================================
  command = ''
  command = line(lowerLim:upperLim)

end subroutine GetCommand
!========================================================
!======================================================================
! Subroutine to find the extent of a block of lines between a command
! and its corresponding END statement in the input script.
! - iLine: Index of the opening command line
! - lineStore: Array of input script lines
! - endCommand: Keyword marking the end of the block (e.g., end_analysis)
! - lineBuffer: Number of lines in the block (including the END command)
! Author:
!   Aliasghar Sepehri
!======================================================================
subroutine FindCommandBlock(iLine, lineStore, endCommand, lineBuffer)
  implicit none
  integer, intent(in) :: iLine
  character(len=maxLineLen), intent(in) :: lineStore(:)   
  character(len=*), intent(in) :: endCommand   
  integer, intent(out) :: lineBuffer

  logical :: found
  integer :: i, lineStat, nLines
  character(len=35) :: command 

  ! Initialize search parameters
  command = " "
  nLines = size(lineStore)
  found = .false.

  !---------------------------
  ! Search for the endCommand
  !---------------------------
  do i = iLine + 1, nLines
    call GetCommand(lineStore(i), command, lineStat)
    call LowerCaseLine(command)

    ! Check if this line matches the expected END command
    if (trim(adjustl(command)) .eq. trim(adjustl(endCommand))) then
      lineBuffer = i - iLine
      found = .true.
      exit
    endif
  enddo

  !-----------------------------------
  ! Error handling if END not found
  !-----------------------------------
  if (.not. found) then
    write(*,*) "ERROR! A command block was opened in the input script, but no closing END statement found!"
    write(*,*) lineStore(iLine)
    stop
  endif

end subroutine FindCommandBlock
!========================================================            
      end module
