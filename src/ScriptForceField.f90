 !========================================================  
module ForceFieldInput
  use VarPrecision        ! Defines dp (double precision kind) and numerical constants
  use Units               ! Handles unit conversions for energy, distance, angle, etc.

  !======================
  ! Constants
  !======================
  integer, parameter :: maxLineLen = 500    ! Maximum character length per line read from file
  real(dp), parameter :: coulombConst = 1.671009770E5_dp  ! Coulomb constant in simulation units

  !======================
  ! Abstract Interfaces
  !======================
  ! These define expected subroutine signatures to be assigned dynamically 
  ! based on the selected force field type

  ! Interface for subroutines with no arguments (used for initialization tasks)
  abstract interface 
    subroutine CommonSub
      use VarPrecision
      implicit none
    end subroutine
  end interface 

  ! Interface for subroutines that take an array of lines from the force field file
  abstract interface 
    subroutine interSub(lineStore)
      implicit none
      character(len=*), intent(in) :: lineStore(:)
    end subroutine
  end interface 

  !======================
  ! Private Members
  !======================
  private                        ! Everything private by default unless made public explicitly
  logical :: fieldTypeSet = .false.  ! Tracks whether a force field type has been initialized

  !======================
  ! Function Pointers
  !======================
  ! These are set based on the chosen force field type. Each force field must provide
  ! its own implementations of these procedures, which are then assigned to the pointers below.

  procedure(CommonSub), pointer :: commonFunction => NULL()
  ! -> Allocates or initializes force-field-specific variables, e.g., potential tables

  procedure(CommonSub), pointer :: FFSpecificFlags => NULL()
  ! -> Optionally sets any force-field-specific toggles or switches used elsewhere

  procedure(interSub), pointer :: interFunction => NULL()
  ! -> Parses and assigns pairwise interaction parameters from the force field input block

  !======================
  ! Public Interfaces
  !======================
  public :: SetForcefieldType      ! Initializes the above pointers based on the type (e.g., HPS, Mpipi)
  public :: ScriptForcefield       ! Invokes interFunction to read parameters from the force field file
  public :: fieldTypeSet           ! Allows external access to check if the force field type has been set

!========================================================  
      contains
!========================================================  
subroutine SetForcefieldType(potenType)
  !=========================================================================
  ! Subroutine : SetForcefieldType
  ! Purpose    : Configure function pointers and data related to the selected
  !              force field type. This sets the core behavior of energy,
  !              pressure, boundary, and Rosenbluth calculations.
  ! Arguments  :
  !   potenType : character string specifying the selected force field
  !=========================================================================
  use SwapBoundary
  use ForceField,     only: ForceFieldName        ! Global name for selected force field
  use EnergyPressurePointers                       ! Interfaces for E/P calculation
  use ParallelVar,   only: nout                    ! Parallel output stream
  implicit none

  character(len=20) :: potenType              ! Input: name of selected force field


  select case(trim(adjustl(potenType)))
  case("lj_q")
    write(nout,*) "Forcefield Type: Standard Lennard-Jones w/ Electrostatic"
    ForceFieldName = "LJ_Q"

    ! Energy calculation function pointers
    Detailed_ECalc         => Detailed_ECalc_LJ_Q
    Shift_ECalc            => Shift_ECalc_LJ_Q
    SwapIn_ECalc           => SwapIn_ECalc_LJ_Q
    SwapOut_ECalc          => SwapOut_ECalc_LJ_Q

    ! Pressure calculation function pointers
    Detailed_PCalc         => Detailed_PCalc_LJ_Q
    Shift_PCalc            => Shift_PCalc_LJ_Q
    SwapIn_PCalc           => SwapIn_PCalc_LJ_Q
    SwapOut_PCalc          => SwapOut_PCalc_LJ_Q

    ! Rosenbluth scheme: molecule-level
    Rosen_Mol_New          => Rosen_Mol_New_LJ_Q
    Rosen_Mol_Old          => Rosen_Mol_Old_LJ_Q

    ! Rosenbluth scheme: atom-level (inter/intra/gas)
    Rosen_Atom_New             => Rosen_Atom_New_LJ_Q
    Rosen_Atom_Old             => Rosen_Atom_Old_LJ_Q
    Rosen_Atom_Intra_New       => Rosen_Atom_Intra_New_LJ_Q
    Rosen_Atom_Intra_Old       => Rosen_Atom_Intra_Old_LJ_Q
    Rosen_Atom_Intra_Gas_Old   => Rosen_Atom_Intra_Gas_Old_LJ_Q

    ! Rosenbluth scheme: long linear chains
    RosenLL_Atom_New           => RosenLL_Atom_New_LJ_Q
    RosenLL_Atom_Old           => RosenLL_Atom_Old_LJ_Q
    RosenLL_Atom_Intra_New     => RosenLL_Atom_Intra_New_LJ_Q
    RosenLL_Atom_Intra_Old     => RosenLL_Atom_Intra_Old_LJ_Q
    RosenLL_Atom_Intra_Gas_Old => RosenLL_Atom_Intra_Gas_Old_LJ_Q

    ! Interaction and boundary methods
    Quick_Nei_ECalc       => QuickNei_ECalc_Inter_LJ_Q
    Find_InterAtms        => Find_InterAtms_LJ_Q
    boundaryFunction      => Bound_MaxMin

    ! Force field-specific allocation and parameter reading
    commonFunction        => Allocate_LJ_Q
    FFSpecificFlags       => LJ_SetFlags
    interFunction         => Read_LJ_Q

    fieldTypeSet = .true.

  case("mpipi")
    write(nout,*) "Forcefield Type: Mpipi for Biomolecules (Proteins and RNA)"
    ForceFieldName = "Mpipi"

    ! Energy calculation function pointers
    Detailed_ECalc         => Detailed_ECalc_Mpipi
    Shift_ECalc            => Shift_ECalc_Mpipi
    SwapIn_ECalc           => SwapIn_ECalc_Mpipi
    SwapOut_ECalc          => SwapOut_ECalc_Mpipi

    ! Pressure calculation function pointers
    Detailed_PCalc         => Detailed_PCalc_Mpipi
    Shift_PCalc            => Shift_PCalc_Mpipi
    SwapIn_PCalc           => SwapIn_PCalc_Mpipi
    SwapOut_PCalc          => SwapOut_PCalc_Mpipi

    ! Rosenbluth scheme: molecule-level
    Rosen_Mol_New          => Rosen_Mol_New_Mpipi
    Rosen_Mol_Old          => Rosen_Mol_Old_Mpipi

    ! Rosenbluth scheme: atom-level (inter/intra/gas)
    Rosen_Atom_New             => Rosen_Atom_New_Mpipi
    Rosen_Atom_Old             => Rosen_Atom_Old_Mpipi
    Rosen_Atom_Intra_New       => Rosen_Atom_Intra_New_Mpipi
    Rosen_Atom_Intra_Old       => Rosen_Atom_Intra_Old_Mpipi
    Rosen_Atom_Intra_Gas_Old   => Rosen_Atom_Intra_Gas_Old_Mpipi

    ! Rosenbluth scheme: long linear chains
    RosenLL_Atom_New           => RosenLL_Atom_New_Mpipi
    RosenLL_Atom_Old           => RosenLL_Atom_Old_Mpipi
    RosenLL_Atom_Intra_New     => RosenLL_Atom_Intra_New_Mpipi
    RosenLL_Atom_Intra_Old     => RosenLL_Atom_Intra_Old_Mpipi
    RosenLL_Atom_Intra_Gas_Old => RosenLL_Atom_Intra_Gas_Old_Mpipi

    ! Interaction and neighbor list functions
    Quick_Nei_ECalc       => QuickNei_ECalc_Inter_Mpipi
    Find_InterAtms        => Find_InterAtms_Mpipi
    boundaryFunction      => Bound_MaxMin

    ! Force field-specific setup routines
    commonFunction        => Allocate_Mpipi
    FFSpecificFlags       => Mpipi_SetFlags
    interFunction         => Read_Mpipi

    fieldTypeSet = .true.


  case("hps_single")
    write(nout,*) "Forcefield Type: HPS_single for Biomolecules (Proteins)"
    ! HPS-Nucl: Single short-range interaction function
    ForceFieldName               = "HPS_single"

    ! Energy calculation function pointers
    Detailed_ECalc              => Detailed_ECalc_HPS_single
    Shift_ECalc                 => Shift_ECalc_HPS_single
    SwapIn_ECalc                => SwapIn_ECalc_HPS_single
    SwapOut_ECalc               => SwapOut_ECalc_HPS_single

    ! Pressure calculation function pointers
    Detailed_PCalc              => Detailed_PCalc_HPS_single
    Shift_PCalc                 => Shift_PCalc_HPS_single
    SwapIn_PCalc                => SwapIn_PCalc_HPS_single
    SwapOut_PCalc               => SwapOut_PCalc_HPS_single

    ! Rosenbluth scheme: molecule-level
    Rosen_Mol_New               => Rosen_Mol_New_HPS_single
    Rosen_Mol_Old               => Rosen_Mol_Old_HPS_single

    ! Rosenbluth scheme: atom-level (inter/intra/gas)
    Rosen_Atom_New              => Rosen_Atom_New_HPS_single
    Rosen_Atom_Old              => Rosen_Atom_Old_HPS_single
    Rosen_Atom_Intra_New        => Rosen_Atom_Intra_New_HPS_single
    Rosen_Atom_Intra_Old        => Rosen_Atom_Intra_Old_HPS_single
    Rosen_Atom_Intra_Gas_Old    => Rosen_Atom_Intra_Gas_Old_HPS_single

    ! Rosenbluth scheme: long linear chains
    RosenLL_Atom_New            => RosenLL_Atom_New_HPS_single
    RosenLL_Atom_Old            => RosenLL_Atom_Old_HPS_single
    RosenLL_Atom_Intra_New      => RosenLL_Atom_Intra_New_HPS_single
    RosenLL_Atom_Intra_Old      => RosenLL_Atom_Intra_Old_HPS_single
    RosenLL_Atom_Intra_Gas_Old  => RosenLL_Atom_Intra_Gas_Old_HPS_single

    ! Interaction and neighbor list functions
    Quick_Nei_ECalc             => QuickNei_ECalc_Inter_HPS_single
    Find_InterAtms              => Find_InterAtms_HPS_single
    boundaryFunction            => Bound_MaxMin

    ! Force field-specific setup routines
    commonFunction              => Allocate_HPS_single
    FFSpecificFlags             => HPS_single_SetFlags
    interFunction               => Read_HPS_single
    fieldTypeSet = .true.

  case("hps_piecewise")
    write(nout,*) "Forcefield Type: HPS_piecewise for Biomolecules (Proteins)"
    ! HPS-KR, KH, FB, TSCL-M2, Urry: Piecewise form, different parameter sets
    ForceFieldName               = "HPS_piecewise"

    ! Energy calculation function pointers
    Detailed_ECalc              => Detailed_ECalc_HPS_piecewise
    Shift_ECalc                 => Shift_ECalc_HPS_piecewise
    SwapIn_ECalc                => SwapIn_ECalc_HPS_piecewise
    SwapOut_ECalc               => SwapOut_ECalc_HPS_piecewise

    ! Pressure calculation function pointers
    Detailed_PCalc              => Detailed_PCalc_HPS_piecewise
    Shift_PCalc                 => Shift_PCalc_HPS_piecewise
    SwapIn_PCalc                => SwapIn_PCalc_HPS_piecewise
    SwapOut_PCalc               => SwapOut_PCalc_HPS_piecewise

    ! Rosenbluth scheme: molecule-level
    Rosen_Mol_New               => Rosen_Mol_New_HPS_piecewise
    Rosen_Mol_Old               => Rosen_Mol_Old_HPS_piecewise

    ! Rosenbluth scheme: atom-level (inter/intra/gas)
    Rosen_Atom_New              => Rosen_Atom_New_HPS_piecewise
    Rosen_Atom_Old              => Rosen_Atom_Old_HPS_piecewise
    Rosen_Atom_Intra_New        => Rosen_Atom_Intra_New_HPS_piecewise
    Rosen_Atom_Intra_Old        => Rosen_Atom_Intra_Old_HPS_piecewise
    Rosen_Atom_Intra_Gas_Old    => Rosen_Atom_Intra_Gas_Old_HPS_piecewise

    ! Rosenbluth scheme: long linear chains
    RosenLL_Atom_New            => RosenLL_Atom_New_HPS_piecewise
    RosenLL_Atom_Old            => RosenLL_Atom_Old_HPS_piecewise
    RosenLL_Atom_Intra_New      => RosenLL_Atom_Intra_New_HPS_piecewise
    RosenLL_Atom_Intra_Old      => RosenLL_Atom_Intra_Old_HPS_piecewise
    RosenLL_Atom_Intra_Gas_Old  => RosenLL_Atom_Intra_Gas_Old_HPS_piecewise

    ! Interaction and neighbor list functions
    Quick_Nei_ECalc             => QuickNei_ECalc_Inter_HPS_piecewise
    Find_InterAtms              => Find_InterAtms_HPS_piecewise
    boundaryFunction            => Bound_MaxMin

    ! Force field-specific setup routines
    commonFunction              => Allocate_HPS_piecewise
    FFSpecificFlags             => HPS_piecewise_SetFlags
    interFunction               => Read_HPS_piecewise
    fieldTypeSet = .true.

  case("hps_cation_pi")
    write(nout,*) "Forcefield Type: HPS_cation_pi for Biomolecules (Proteins)"
    ! HPS-cation-pi-i and -ii: Captures cation–π interactions explicitly
    ForceFieldName               = "HPS_cation_pi"

    ! Energy calculation function pointers
    Detailed_ECalc              => Detailed_ECalc_HPS_cation_pi
    Shift_ECalc                 => Shift_ECalc_HPS_cation_pi
    SwapIn_ECalc                => SwapIn_ECalc_HPS_cation_pi
    SwapOut_ECalc               => SwapOut_ECalc_HPS_cation_pi

    ! Pressure calculation function pointers
    Detailed_PCalc              => Detailed_PCalc_HPS_cation_pi
    Shift_PCalc                 => Shift_PCalc_HPS_cation_pi
    SwapIn_PCalc                => SwapIn_PCalc_HPS_cation_pi
    SwapOut_PCalc               => SwapOut_PCalc_HPS_cation_pi

    ! Rosenbluth scheme: molecule-level
    Rosen_Mol_New               => Rosen_Mol_New_HPS_cation_pi
    Rosen_Mol_Old               => Rosen_Mol_Old_HPS_cation_pi

    ! Rosenbluth scheme: atom-level (inter/intra/gas)
    Rosen_Atom_New              => Rosen_Atom_New_HPS_cation_pi
    Rosen_Atom_Old              => Rosen_Atom_Old_HPS_cation_pi
    Rosen_Atom_Intra_New        => Rosen_Atom_Intra_New_HPS_cation_pi
    Rosen_Atom_Intra_Old        => Rosen_Atom_Intra_Old_HPS_cation_pi
    Rosen_Atom_Intra_Gas_Old    => Rosen_Atom_Intra_Gas_Old_HPS_cation_pi

    ! Rosenbluth scheme: long linear chains
    RosenLL_Atom_New            => RosenLL_Atom_New_HPS_cation_pi
    RosenLL_Atom_Old            => RosenLL_Atom_Old_HPS_cation_pi
    RosenLL_Atom_Intra_New      => RosenLL_Atom_Intra_New_HPS_cation_pi
    RosenLL_Atom_Intra_Old      => RosenLL_Atom_Intra_Old_HPS_cation_pi
    RosenLL_Atom_Intra_Gas_Old  => RosenLL_Atom_Intra_Gas_Old_HPS_cation_pi

    ! Interaction and neighbor list functions
    Quick_Nei_ECalc             => QuickNei_ECalc_Inter_HPS_cation_pi
    Find_InterAtms              => Find_InterAtms_HPS_cation_pi
    boundaryFunction            => Bound_MaxMin

    ! Force field-specific setup routines
    commonFunction              => Allocate_HPS_cation_pi
    FFSpecificFlags             => HPS_cation_pi_SetFlags
    interFunction               => Read_HPS_cation_pi
    fieldTypeSet = .true.

  case("pedone")
!     ─────────────────────────────────────────────────────────────
!     Force field: Pedone
!     - Designed for ionic materials and oxides
!     - Includes detailed and shifted energy calculations
!     - No Rosenbluth scheme at atomic level
!     - Uses a specialized boundary condition to enforce charge neutrality
!     ─────────────────────────────────────────────────────────────
    write(nout,*) "Forcefield Type: Pedone"
    ForceFieldName = "Pedone"

!       Assign energy calculation subroutines
    Detailed_ECalc           => Detailed_ECalc_Pedone
    Shift_ECalc              => Shift_ECalc_Pedone
    SwapIn_ECalc             => SwapIn_ECalc_Pedone
    SwapOut_ECalc            => SwapOut_ECalc_Pedone

!       Assign Rosenbluth scheme at molecule level only
    Rosen_Mol_New            => Rosen_Mol_New_Pedone
    Rosen_Mol_Old            => Rosen_Mol_Old_Pedone

!       Assign fast interaction routine for neighbor list
    Quick_Nei_ECalc          => QuickNei_ECalc_Inter_Pedone

!       Use custom boundary condition for charge balancing
    boundaryFunction         => Bound_PedoneChargeBalance

!       Assign allocation, flags, and input reader for this force field
    commonFunction           => Allocate_Pedone
    FFSpecificFlags          => Pedone_SetFlags
    interFunction            => Read_Pedone

!       Confirm that a force field type has been set
    fieldTypeSet = .true.

  case("tersoff")
!     ─────────────────────────────────────────────────────────────
!     Force field: Tersoff
!     - Suitable for covalent materials (e.g., Si, C)
!     - Includes many-body interactions via bond-order potentials
!     - No atomic-level Rosenbluth scheme; only molecular-level regrowth
!     ─────────────────────────────────────────────────────────────
    write(nout,*) "Forcefield Type: Tersoff"
    ForceFieldName = "Tersoff"

!       Assign energy calculation routines
    Detailed_ECalc           => Detailed_ECalc_Tersoff
    Shift_ECalc              => Shift_ECalc_Tersoff
    SwapIn_ECalc             => SwapIn_ECalc_Tersoff
    SwapOut_ECalc            => SwapOut_ECalc_Tersoff

!       Only molecule-level Rosenbluth regrowth is supported
    Rosen_Mol_New            => Rosen_Mol_New_Tersoff
    Rosen_Mol_Old            => Rosen_Mol_Old_Tersoff

!       Assign fast interaction routine for neighbor list
    Quick_Nei_ECalc          => QuickNei_ECalc_Inter_Tersoff

!       Use general bounding condition
    boundaryFunction         => Bound_MaxMin

!       Assign allocation, flags, and parameter reader
    commonFunction           => Allocate_Tersoff
    FFSpecificFlags          => Tersoff_SetFlags
    interFunction            => Read_Tersoff

!       Confirm that a force field type has been set
    fieldTypeSet = .true.

  case("custompairwise")
!     ─────────────────────────────────────────────────────────────
!     Placeholder for user-defined pairwise force fields
!     - Function pointers must be manually set elsewhere
!     - Useful for testing or extending with new models
!     ─────────────────────────────────────────────────────────────
!     No assignment here; the user must set up the behavior manually.

  case default
!     ─────────────────────────────────────────────────────────────
!     Catch invalid or unknown force field types
!     - Prints error to output
!     - Stops execution to avoid undefined behavior
!     ─────────────────────────────────────────────────────────────
    write(nout,*) potenType, "Unknown potential type given in forcefield input"
    stop "Unknown potential type given in forcefield input"
  end select

end subroutine SetForcefieldType
!========================================================            
subroutine ScriptForcefield(lineStore)
  !=========================================================================
  ! Subroutine: ScriptForcefield
  ! Purpose   : Parses the force field file content stored in lineStore.
  !             The force field type and number of molecule types must
  !             already be set before calling this routine.
  !=========================================================================
  use CBMC_Variables, only: maxRosenTrial
  use Coords,          only: rosenTrial
  use ForceField,      only: atomData, atomArray, nonBondArray, bondArray, bendArray, &
                             torsArray, totalMass, nAtomTypes, nBondTypes, nAngleTypes, &
                             nTorsionalTypes, nAtoms
  use SimParameters,   only: nMolTypes, NMAX, vmdAtoms
  use VarPrecision,    only: dp
  use CoordinateFunctions, only: AllocateCoordinateArrays
  use ParallelVar,     only: nout
  implicit none

  character(len=maxLineLen), intent(in) :: lineStore(:)  ! Force field lines

  character(len=25) :: command                           ! Parsed command
  integer :: i, j, iLine, nLines                         ! Loop counters
  integer :: lineBuffer, lineStat                        ! Control flags
  integer :: AllocateStat                                ! Allocation status
  integer :: nAtomsMax, nBondsMax, nAnglesMax            ! Allocation sizes
  integer :: nTorsMax, nNonBondMax
      
  !=========================================================================
  ! Error checking before parsing:
  ! Ensure that a valid force field type has been set
  ! and that the number of molecule types is defined.
  !=========================================================================
  if(.not. fieldTypeSet) then
    write(nout,*) "ERROR! Forcefield input has been called before the forcefield type has been set!"
    stop
  endif

  if(nMolTypes .eq. 0) then
    write(nout,*) "ERROR! Forcefield input has been called before the number of molecule types have been defined"
    stop
  endif

  !=========================================================================
  ! Initialize counts of interaction types
  !=========================================================================
  nAtomTypes      = 0
  nBondTypes      = 0
  nAngleTypes     = 0
  nTorsionalTypes = 0

  !===============================================================
  ! Determine maximum required sizes for atom, bond, angle, etc.
  ! by scanning the forcefield input lines for each entity type.
  !===============================================================
  call FindMax(lineStore, "atoms",      nAtomsMax)
  write(34,*) "Maximum number of atoms per molecule type:", nAtomsMax

  call FindMax(lineStore, "nonbonded",  nNonBondMax)
  write(34,*) "Maximum number of nonbonded entries per molecule type:", nNonBondMax

  call FindMax(lineStore, "bonds",      nBondsMax)
  write(34,*) "Maximum number of bonds per molecule type:", nBondsMax

  call FindMax(lineStore, "angles",     nAnglesMax)
  write(34,*) "Maximum number of angles per molecule type:", nAnglesMax

  call FindMax(lineStore, "torsional",  nTorsMax)
  write(34,*) "Maximum number of torsions per molecule type:", nTorsMax


  !===============================================================
  ! Allocate arrays for storing forcefield parameters.
  ! Dimensions: [# molecule types, # max entries for each type]
  !===============================================================
  allocate(atomArray(1:nMolTypes,1:nAtomsMax),      STAT=AllocateStat)
  allocate(nonBondArray(1:nMolTypes,1:nNonBondMax), STAT=AllocateStat)
  allocate(bondArray(1:nMolTypes,1:nBondsMax),      STAT=AllocateStat)
  allocate(bendArray(1:nMolTypes,1:nAnglesMax),     STAT=AllocateStat)
  allocate(torsArray(1:nMolTypes,1:nTorsMax),       STAT=AllocateStat)

  !===============================================================
  ! Allocate Rosenbluth trial coordinate arrays (for CBMC/AVBMC)
  ! These arrays store trial x/y/z coordinates for each atom
  ! across a maximum number of trials.
  !===============================================================
  do i = 1, maxRosenTrial
    allocate(rosenTrial(i)%x(1:nAtomsMax))
    allocate(rosenTrial(i)%y(1:nAtomsMax))
    allocate(rosenTrial(i)%z(1:nAtomsMax))
  end do

  lineBuffer = 0
  nLines = size(lineStore)

  ! Loop through all lines in the force field input
  do iLine = 1, nLines

    ! If lineBuffer is positive, we are skipping previously processed lines
    if(lineBuffer .gt. 0) then
      lineBuffer = lineBuffer - 1
      cycle
    endif

    ! Extract and process the command on the current line
    call GetCommand(lineStore(iLine), command, lineStat)    
    if(lineStat .gt. 0) then
      cycle
    endif 

    call LowerCaseLine(command)

    ! Interpret command and dispatch appropriate subroutine
    select case( adjustl(trim(command)) )

      case("define")
        ! Handle 'define' block
        call FindCommandBlock(iLine, lineStore, "end_define", lineBuffer)
        call DefineForcefield(lineStore(iLine:iLine+lineBuffer) )

      case("create")
        ! Handle 'create' block
        call FindCommandBlock(iLine, lineStore, "end_create", lineBuffer)
        call CreateForcefield(lineStore(iLine:iLine+lineBuffer) )

      case("set")
        ! Handle single-line 'set' command
        call SetForcefieldParam(lineStore(iLine))

      case default
        ! Report unknown commands
        write(nout,*) "ERROR! Invalid command in Forcefield file on line:", iLine
        write(nout,*) lineStore(iLine)
        stop "INPUT ERROR!"

    end select

  enddo

  vmdAtoms = 0
  do i = 1,nMolTypes
    vmdAtoms = vmdAtoms + (NMAX(i)*nAtoms(i))
  enddo
      
  totalMass = 0d0
  do i = 1, nMolTypes
    do j = 1,nAtoms(i)
      totalMass(i) = totalMass(i) + atomData(atomArray(i,j))%mass
    enddo
  enddo
      
  call AllocateCoordinateArrays
  call FFSpecificFlags
   
end subroutine ScriptForcefield
!========================================================            
subroutine DefineForcefield(lineStore)
  !----------------------------------------------------------------------
  ! This subroutine parses a 'define' block from the forcefield script.
  ! It initializes parameters such as number of molecules, atoms, etc.,
  ! and stores topological and geometric definitions for each molecule type.
  !----------------------------------------------------------------------
  use ForceField,     only: atomData, bondData, bendData, torsData, nAtomTypes, &
                            nBondTypes, nAngleTypes, nTorsionalTypes
  use SimParameters,  only: beta, echoInput, convEng, convDist, convAng
  use VarPrecision,   only: dp
  use ParallelVar,    only: nout
  implicit none

  character(len=maxLineLen), intent(in) :: lineStore(:)
  character(len=30) :: dummy, defType
  integer :: i, j, iLine, nLines, nParam
  integer :: intValue, AllocateStat

  nLines = size(lineStore)

  write(34,*) "Entering subroutine DefineForcefield"
  write(34,*) "Total number of lines to parse in 'define' block:", nLines

  ! Read the definition type (e.g., atoms, bonds, etc.) from the first line
  read(lineStore(1),*) dummy, defType
  call LowerCaseLine(defType)

  write(34,*) "DefineForcefield: Parsing definition type =", trim(defType)
  write(34,*) "Number of lines in definition block =", nLines

  select case( trim(adjustl(defType)) )

    case("atomtypes")
    ! Check if atom types have already been defined
      if( allocated(atomData) ) then
        write(nout,*) "ERROR! AtomTypes already defined in the forcefield file"
        stop
      endif

    ! Read the number of atom types
      read(lineStore(1),*) dummy, defType, intValue
      nAtomTypes = intValue

      write(34,*) "DefineForcefield: Defining AtomTypes, number =", nAtomTypes

    ! Allocate memory for atom data
      ALLOCATE(atomData(1:nAtomTypes), STAT = AllocateStat)

    ! Call forcefield-specific allocation routine
      call commonFunction

    ! Call forcefield-specific interaction parameter reader
      call interFunction(lineStore)

    case("bondtypes")
    ! Read number of bond types
      read(lineStore(1),*) dummy, defType, intValue
      nBondTypes = intValue

      write(34,*) "DefineForcefield: Defining BondTypes, number =", nBondTypes

    ! Allocate array to hold bond parameters
      ALLOCATE(bondData(1:nBondTypes), STAT = AllocateStat)

      i = 0
      do iLine = 2, nLines-1
        i = i + 1

      ! Read name, force constant, and equilibrium distance for this bond type
        read(lineStore(iLine),*) bondData(i)%bondName, bondData(i)%k_eq, bondData(i)%r_eq

      ! Convert force constant and distance to internal units
        bondData(i)%k_eq = bondData(i)%k_eq * convEng
        bondData(i)%r_eq = bondData(i)%r_eq * convDist

      ! Compute Gaussian width of bond length distribution (r_sigma)
      ! and maximum allowed squared bond distance (rmax_sq) based on 5σ criterion
        if (bondData(i)%k_eq .gt. 0E0_dp) then
          bondData(i)%r_sigma = 1E0_dp / sqrt(beta * bondData(i)%k_eq)
          bondData(i)%rmax_sq = (bondData(i)%r_eq + 5E0_dp * bondData(i)%r_sigma)**2 
        else
          bondData(i)%r_sigma = 0E0_dp
          bondData(i)%rmax_sq = bondData(i)%r_eq
        endif

      ! Optionally echo input to file 35
        if(echoInput) then
          write(35,*) bondData(i)%bondName, bondData(i)%k_eq, bondData(i)%r_eq    
        endif        
      enddo     

    case("angletypes")
      ! Read number of angle types
      read(lineStore(1),*) dummy, defType, intValue
      nAngleTypes = intValue

      write(34,*) "DefineForcefield: Defining AngleTypes, number =", nAngleTypes

      ! Allocate array to hold angle parameters
      ALLOCATE(bendData(1:nAngleTypes), STAT = AllocateStat)

      i = 0
      do iLine = 2, nLines-1
        i = i + 1
        ! Read name, force constant, and equilibrium angle
        read(lineStore(iLine),*) bendData(i)%angleName, bendData(i)%k_eq, bendData(i)%ang_eq

        ! Optionally echo raw values from input
        if(echoInput) then
          write(35,*) bendData(i)%angleName, bendData(i)%k_eq, bendData(i)%ang_eq 
        endif          

        ! Convert to internal energy and angle units
        bendData(i)%k_eq = bendData(i)%k_eq * convEng
        bendData(i)%ang_eq = bendData(i)%ang_eq * convAng
      enddo
 
    case("torsiontypes")
      ! Read number of torsional types
      read(lineStore(1),*) dummy, defType, intValue
      nTorsionalTypes = intValue

      write(34,*) "DefineForcefield: Defining TorsionTypes, number =", nTorsionalTypes

      ! Allocate torsional parameter array
      ALLOCATE(torsData(1:nTorsionalTypes), STAT = AllocateStat)

      i = 0
      do iLine = 2, nLines - 1
        i = i + 1
        ! First parse the name and number of parameters
        read(lineStore(iLine),*) torsData(i)%torsName, nParam
        allocate(torsData(i)%a(1:nParam), STAT = AllocateStat)

        ! Read the torsion name, number, and the parameter list
        read(lineStore(iLine),*) torsData(i)%torsName, torsData(i)%nPara, (torsData(i)%a(j), j = 1, nParam)

        ! Optionally echo input
        if(echoInput) then
          write(35,*) torsData(i)%torsName, torsData(i)%nPara, (torsData(i)%a(j), j = 1, nParam)
        endif         

        ! Convert all parameters to internal energy units
        do j = 1, nParam
          torsData(i)%a(j) = torsData(i)%a(j) * convEng
        enddo
      enddo
            
    case("rmin")
      ! Handle RMin setting: ensures atom types have been defined
      if(nAtomTypes .eq. 0) then
        write(nout,*) "ERROR! RMin is called before the number of atom types"
        write(nout,*) "have been defined"
        stop
      endif
      write(34,*) "DefineForcefield: Calling SetRMin"
      call SetRMin(lineStore)

    case default 
      ! Catch unrecognized block types
      write(nout,*) "ERROR! Invalid command in Forcefield file on line:"
      write(nout,*) lineStore(1)
      stop "INPUT ERROR!"
    end select
    flush(34)
     
end subroutine DefineForcefield
!========================================================            
subroutine SetForcefieldParam(line)
  ! Sets unit conversion parameters for the force field using a keyword-based input line.
  use SimParameters,  only: convDist, convEng, convAng
  use Units,          only: FindEngUnit, FindLengthUnit, FindAngularUnit
  use VarPrecision,   only: dp
  implicit none
  character(len=maxLineLen), intent(in) :: line
  character(len=30) :: dummy, deftype, stringValue

  ! Parse the line to extract the keyword
  read(line,*) dummy, defType
  call LowerCaseLine(defType)

  ! Determine the type of unit to set and call corresponding conversion function
  select case(adjustl(trim(defType)))

    case("lenunits")
      read(line,*) dummy, defType, stringValue
      convDist = FindLengthUnit(stringValue)
      write(34,*) "Length unit set to:", trim(stringValue), "Conversion factor:", convDist

    case("engunits")
      read(line,*) dummy, defType, stringValue
      convEng = FindEngUnit(stringValue)
      write(34,*) "Energy unit set to:", trim(stringValue), "Conversion factor:", convEng

    case("angunits")
      read(line,*) dummy, defType, stringValue
      convAng = FindAngularUnit(stringValue)
      write(34,*) "Angular unit set to:", trim(stringValue), "Conversion factor:", convAng

  end select

end subroutine SetForcefieldParam
!========================================================            
!***************************************************************************
! Subroutine: CreateForcefield
! Purpose   : Defines each molecule type in the simulation.
!             For each molecule, it parses:
!             - number of atoms and atom types
!             - bond lists and bond members
!             - angle lists and angle members
!             - torsion lists and torsion members
!             - intra-nonbonded pairs
!             These definitions are typically read from the "create" block in the input file.
!***************************************************************************            
subroutine CreateForcefield(lineStore)
  use ForceField,     only: atomArray, bondArray, bendArray, torsArray, nonBondArray, &
                             nAtoms, nBonds, nAngles, nTorsional, nIntraNonBond
  use ParallelVar,    only: nout
  implicit none
  character(len=*), intent(in) :: lineStore(:)
  character(len=30) :: dummy      ! Temporary string for reading initial label
  character(len=25) :: command    ! To store the command read from each line
  integer :: j, iUnit, iLine, jLine, nLines, nMol
  integer :: intValue, lineStat, lineBuffer

  nLines = size(lineStore)        ! Total number of input lines
  read(lineStore(1),*) dummy, nMol  ! Read number of molecules from first line



  write(34,*) "Reading force field definitions for molecule type:", nMol

!     Loop through all remaining lines in the input block
  do iLine = 2, nLines - 1
!     Extract the command from the current line
    call GetCommand(lineStore(iLine), command, lineStat)
    if(lineStat .gt. 0) then
      cycle
    endif
    call LowerCaseLine(command)

!       Identify the section by command keyword
    select case(adjustl(trim(command)))
!       Handle atom definitions
      case("atoms")
!         Locate the block ending with 'end_atoms'
        call FindCommandBlock(iLine, lineStore, "end_atoms", lineBuffer)
        read(lineStore(iLine),*) dummy, intValue

!           Check consistency in atom count
        if(lineBuffer-1 .ne. intValue) then
          write(nout,*) "Error in Forcefield Def! The Number of Atoms specified is different than the number of lines given."
          write(nout,*) "Number of atoms specified:", intValue
          write(nout,*) "Number of lines in command block:", lineBuffer-1
          stop
        endif

!           Store number of atoms for this molecule type
        nAtoms(nMol) = intValue
        write(34,*) "Number of atoms in molecule", nMol, ":", intValue

!           Read each atom type into atomArray
        iUnit = 0
        do jLine = iLine+1, iLine+lineBuffer-1
          iUnit = iUnit + 1
          read(lineStore(jLine), *) atomArray(nMol, iUnit)
          write(34,*) "  Atom", iUnit, "Type:", atomArray(nMol, iUnit)
        enddo
      case("bonds")
!           Locate the block ending with 'end_bonds'
        call FindCommandBlock(iLine, lineStore, "end_bonds", lineBuffer)
        read(lineStore(iLine),*) dummy, intValue

!           Store number of bonds for this molecule type
        nBonds(nMol) = intValue
        write(34,*) "Number of bonds in molecule", nMol, ":", intValue

!           Read bond definitions into bondArray
        iUnit = 0
        do jLine = iLine+1, iLine+lineBuffer-1
          iUnit = iUnit + 1
          read(lineStore(jLine), *) bondArray(nMol, iUnit)%bondType, (bondArray(nMol, iUnit)%bondMembr(j), j=1,2)
          write(34,*) "  Bond", iUnit, "Type:", bondArray(nMol, iUnit)%bondType, &
                      "Members:", bondArray(nMol, iUnit)%bondMembr(1), bondArray(nMol, iUnit)%bondMembr(2)
        enddo

      case("angles")
!           Locate the block ending with 'end_angles'
        call FindCommandBlock(iLine, lineStore, "end_angles", lineBuffer)
        read(lineStore(iLine),*) dummy, intValue

!           Store number of angle interactions for this molecule type
        nAngles(nMol) = intValue
        write(34,*) "Number of angles in molecule", nMol, ":", intValue

!           Read angle definitions into bendArray
        iUnit = 0
        do jLine = iLine+1, iLine+lineBuffer-1
          iUnit = iUnit + 1
          read(lineStore(jLine), *) bendArray(nMol, iUnit)%bendType, (bendArray(nMol, iUnit)%bendMembr(j), j=1,3)
          write(34,*) "  Angle", iUnit, "Type:", bendArray(nMol, iUnit)%bendType, &
                      "Members:", bendArray(nMol, iUnit)%bendMembr(1), bendArray(nMol, iUnit)%bendMembr(2), bendArray(nMol, iUnit)%bendMembr(3)
        enddo

      case("torsional")
!           Locate the torsional block
        call FindCommandBlock(iLine, lineStore, "end_torsional", lineBuffer)
        read(lineStore(iLine),*) dummy, intValue

!           Store number of torsional interactions for this molecule type
        nTorsional(nMol) = intValue
        write(34,*) "Number of torsional interactions in molecule", nMol, ":", intValue

!           Read torsional interaction definitions
        iUnit = 0
        do jLine = iLine+1, iLine+lineBuffer-1
          iUnit = iUnit + 1
          read(lineStore(jLine), *) torsArray(nMol, iUnit)%torsType, (torsArray(nMol, iUnit)%torsMembr(j), j=1,4)
          write(34,*) "  Torsional", iUnit, "Type:", torsArray(nMol, iUnit)%torsType, &
                      "Members:", torsArray(nMol, iUnit)%torsMembr(1), torsArray(nMol, iUnit)%torsMembr(2), &
                         torsArray(nMol, iUnit)%torsMembr(3), torsArray(nMol, iUnit)%torsMembr(4)
        enddo

      case("nonbonded")
!           Locate the nonbonded intramolecular block
        call FindCommandBlock(iLine, lineStore, "end_nonbonded", lineBuffer)
        read(lineStore(iLine),*) dummy, intValue

!           Store number of intramolecular nonbonded pairs for this molecule type
        nIntraNonBond(nMol) = intValue
        write(34,*) "Number of intra-nonbonded pairs in molecule", nMol, ":", intValue

!           Read nonbonded pair definitions
        iUnit = 0
        do jLine = iLine+1, iLine+lineBuffer-1
          iUnit = iUnit + 1
          read(lineStore(jLine), *) (nonBondArray(nMol, iUnit)%nonMembr(j), j=1,2)
          write(34,*) "  Nonbonded pair", iUnit, ":", nonBondArray(nMol, iUnit)%nonMembr(1), &
                         nonBondArray(nMol, iUnit)%nonMembr(2)
        enddo

    end select
  enddo

end subroutine CreateForcefield
!========================================================
!  This subroutine searches forward in lineStore starting at iLine + 1
!  to find the position of the corresponding endCommand (e.g., "end_define").
!  The number of lines in the command block (including the endCommand line)
!  is returned in lineBuffer.
!
!  Arguments:
!    iLine       : index of the current line where the block begins
!    lineStore   : array of input lines to search through
!    endCommand  : expected closing command (e.g., "end_define")
!    lineBuffer  : output, number of lines until endCommand is found
!=======================================================================
subroutine FindCommandBlock(iLine, lineStore, endCommand, lineBuffer)
  implicit none

  integer, intent(in)               :: iLine
  character(len=*), intent(in)     :: lineStore(:)
  character(len=*), intent(in)     :: endCommand
  integer, intent(out)             :: lineBuffer

  logical                          :: found
  integer                          :: i, lineStat, nLines
  character(len=35)               :: dummy

  dummy   = " "
  found   = .false.
  nLines  = size(lineStore)

  ! Loop from the next line to the end to find the endCommand
  do i = iLine + 1, nLines
    call GetCommand(lineStore(i), dummy, lineStat)
    if (trim(adjustl(dummy)) .eq. trim(adjustl(endCommand))) then
      lineBuffer = i - iLine
      found = .true.
      exit
    endif
  enddo

  ! If endCommand not found, stop with an error
  if (.not. found) then
    write(*,*) "ERROR! A command block was opened in the input script, but no closing END statement found!"
    write(*,*) lineStore(iLine)
    stop
  endif

end subroutine FindCommandBlock
!========================================================            
!     This subroutine sets the r_min mixing rules and values from the forcefield input
subroutine SetRMin(lineStore)
  use VarPrecision,         only: dp
  use ForceField,           only: nAtomTypes, r_min, r_min_tab
  use ForceFieldFunctions,  only: MixRule, GeoMean_MixingFunc, Mean_MixingFunc, &
                                   Min_MixingFunc, Max_MixingFunc
  use SimParameters,        only: echoInput
  implicit none
  character(len=maxLineLen), intent(in) :: lineStore(:)  

  logical :: custom
  integer :: i, j, iLine, nLines
  integer :: indx1, indx2
  character(len=30) :: dummy, dummy2, mixingRule
  procedure (MixRule), pointer :: rmin_func => null()   
  real(dp) :: curVal

  nLines = size(lineStore)
  read(lineStore(1),*) dummy, dummy2, mixingRule

  write(34,*) "SetRMin: Mixing rule or mode selected =>", trim(mixingRule)

  custom = .false.
  select case(trim(adjustl(mixingRule)))
    case("geometric")
      rmin_func => GeoMean_MixingFunc
      write(34,*) "SetRMin: Using geometric mixing rule"
    case("average")
      rmin_func => Mean_MixingFunc
      write(34,*) "SetRMin: Using average mixing rule"
    case("min")
      rmin_func => Min_MixingFunc
      write(34,*) "SetRMin: Using min mixing rule"
    case("max")
      rmin_func => Max_MixingFunc
      write(34,*) "SetRMin: Using max mixing rule"
    case("custom")
      custom = .true. 
      write(34,*) "SetRMin: Using custom r_min table"
    case default
      write(34,*) "SetRMin: ERROR - invalid mixing rule type: ", trim(mixingRule)
      stop "Error! Invalid R_Min Mixing Rule Type"
  end select

  r_min = 0E0_dp
  r_min_tab = 0E0_dp
  if(custom) then           
    do iLine = 2, nLines-1
      read(lineStore(iLine), *) indx1, indx2, curVal 
      if(indx1 .gt. nAtomTypes) then
        write(34,*) "SetRMin: Skipping invalid indx1 =", indx1
        cycle
      endif
      if(indx2 .gt. nAtomTypes) then
        write(34,*) "SetRMin: Skipping invalid indx2 =", indx2
        cycle
      endif
      if(echoInput) then
        write(35,*) indx1, indx2, curVal
        flush(35)
      endif
      r_min_tab(indx1, indx2) = curVal
      r_min_tab(indx2, indx1) = curVal
      write(34,*) "SetRMin: Custom r_min(", indx1, ",", indx2, ") =", curVal
    enddo
    do i = 1,nAtomTypes
      do j = 1,nAtomTypes
        r_min_tab(i,j) = r_min_tab(i,j)**2
      enddo
    enddo
  else
    do iLine = 2, nLines-1
      read(lineStore(iLine), *) indx1, curVal 
      r_min(indx1) = curVal
      write(34,*) "SetRMin: r_min(", indx1, ") =", curVal
    enddo
    do i = 1,nAtomTypes
      do j = i,nAtomTypes
        r_min_tab(i,j) = rmin_func(r_min(i), r_min(j))**2
        r_min_tab(j,i) = r_min_tab(i,j)
        write(34,*) "SetRMin: Mixed r_min(", i, ",", j, ") = ", sqrt(r_min_tab(i,j))
      enddo
    enddo
  endif

  write(35,*) "Rmin Table:"
  do i = 1, nAtomTypes
    write(35,*) (sqrt(r_min_tab(i,j)), j= 1, nAtomTypes)
  enddo
  write(34,*) "SetRMin: Completed setting r_min table"

end subroutine SetRMin
!========================================================            
!  This subroutine extracts the first command from a line.
!  It skips leading spaces, ignores comment lines (starting with '#'),
!  and returns the first word (token) as 'command'.
!
!  Arguments:
!    line      : input string (a line from the input file)
!    command   : output command string (first word on the line)
!    lineStat  : status flag
!                = 1 → line is empty or comment
!                = 0 → valid command found
!===========================================================
subroutine GetCommand(line, command, lineStat)
  use VarPrecision
  implicit none

  character(len=*), intent(in)     :: line
  character(len=25), intent(out)   :: command
  integer, intent(out)             :: lineStat

  integer :: i, sizeLine, lowerLim, upperLim

  sizeLine = len(line)
  lineStat = 0
  command  = " "
  i = 1

  ! Find the first non-blank character in the line
  do while(i <= sizeLine)
    if (ichar(line(i:i)) /= ichar(' ')) then
      ! If it's a comment line, mark it and exit
      if (ichar(line(i:i)) == ichar('#')) then
        lineStat = 1
        return
      else
        exit
      endif
    endif
    i = i + 1
  enddo

  ! If only spaces or nothing is found, mark as empty
  if (i >= sizeLine) then
    lineStat = 1
    return
  endif

  lowerLim = i

  ! Continue scanning until a space or end of string
  do while(i <= sizeLine)
    if (line(i:i) == " ") then
      exit
    endif
    i = i + 1
  enddo

  upperLim = i

  ! Extract the command
  command = line(lowerLim:upperLim)

end subroutine GetCommand
!================================================================ 
subroutine Allocate_LJ_Q
  use SimParameters,        only: nMolTypes, max_dist, max_dist_single, max_rot
  use ForceField,           only: r_min_sq, totalMass, nAtoms, nIntraNonBond, nBonds, nAngles, &
                                  nTorsional, r_min, r_min_tab, nAtomTypes, nImproper
  use ForceFieldPara_LJ_Q,  only: ep_tab, sig_tab, q_tab
  use AcceptRates,          only: acptTrans, atmpTrans, acptRot, atmpRot, &
                                  acptSwapIn, atmpSwapIn, acptSwapOut, atmpSwapOut
  use VarPrecision,         only: dp
  implicit none
  integer :: AllocateStatus

  !---------------------------------------------------------------
  ! Allocate and initialize LJ_Q-specific arrays and parameters
  ! This subroutine is called by commonFunction pointer after
  ! reading 'atomtypes' block for LJ_Q force field.
  !---------------------------------------------------------------

  ! Allocate van der Waals and electrostatics parameter arrays
  ALLOCATE (r_min(1:nAtomTypes),             STAT = AllocateStatus)
  ALLOCATE (r_min_sq(1:nAtomTypes),          STAT = AllocateStatus)
  ALLOCATE (r_min_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (ep_tab(1:nAtomTypes,1:nAtomTypes),     STAT = AllocateStatus)
  ALLOCATE (sig_tab(1:nAtomTypes,1:nAtomTypes),    STAT = AllocateStatus)
  ALLOCATE (q_tab(1:nAtomTypes,1:nAtomTypes),      STAT = AllocateStatus)

  ! Initialize interaction tables to zero
  ep_tab      = 0E0_dp
  sig_tab     = 0E0_dp
  q_tab       = 0E0_dp
  r_min       = 0E0_dp
  r_min_sq    = 0E0_dp
  r_min_tab   = 0E0_dp

  ! Allocate molecule-type-level parameters
  ALLOCATE (totalMass(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nAtoms(1:nMolTypes),           STAT = AllocateStatus)
  ALLOCATE (nIntraNonBond(1:nMolTypes),    STAT = AllocateStatus)
  ALLOCATE (nBonds(1:nMolTypes),           STAT = AllocateStatus)
  ALLOCATE (nAngles(1:nMolTypes),          STAT = AllocateStatus)
  ALLOCATE (nTorsional(1:nMolTypes),       STAT = AllocateStatus)
  ALLOCATE (nImproper(1:nMolTypes),        STAT = AllocateStatus)

  ! Initialize molecule counts to 0
  nAtoms        = 0
  nIntraNonBond = 0
  nBonds        = 0
  nAngles       = 0
  nTorsional    = 0
  nImproper     = 0

  ! Allocate acceptance statistics arrays for each molecule type
  ALLOCATE (acptTrans(1:nMolTypes),     STAT = AllocateStatus)
  ALLOCATE (atmpTrans(1:nMolTypes),     STAT = AllocateStatus)
  ALLOCATE (acptRot(1:nMolTypes),       STAT = AllocateStatus)
  ALLOCATE (atmpRot(1:nMolTypes),       STAT = AllocateStatus)
  ALLOCATE (acptSwapIn(1:nMolTypes),    STAT = AllocateStatus)
  ALLOCATE (atmpSwapIn(1:nMolTypes),    STAT = AllocateStatus)
  ALLOCATE (acptSwapOut(1:nMolTypes),   STAT = AllocateStatus)
  ALLOCATE (atmpSwapOut(1:nMolTypes),   STAT = AllocateStatus)

  ! Allocate maximum displacement/rotation trackers
  ALLOCATE (max_dist(1:nMolTypes),         STAT = AllocateStatus)
  ALLOCATE (max_dist_single(1:nMolTypes),  STAT = AllocateStatus)
  ALLOCATE (max_rot(1:nMolTypes),          STAT = AllocateStatus)

  ! Terminate if memory allocation fails
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

end subroutine Allocate_LJ_Q
!================================================================ 
subroutine Allocate_Mpipi
  use SimParameters,         only: nMolTypes, max_dist, max_dist_single, max_rot
  use AcceptRates,           only: acptTrans, atmpTrans, acptRot, atmpRot, acptSwapIn, &
                                   atmpSwapIn, acptSwapOut, atmpSwapOut, acptCBMC, atmpCBMC, &
                                   acptFECBMC, atmpFECBMC, nCBMCmax, nFECBMCmax
  use ForceField,            only: totalMass, nAtoms, nIntraNonBond, nBonds, nAngles, &
                                   r_min_sq, nTorsional, r_min, r_min_tab, nAtomTypes, nImproper
  use ForceFieldPara_Mpipi,  only: sig_tab, q_tab, Mu_tab, sigsq_tab, &
                                   epAlpha_tab, q_Nonzero
  use VarPrecision,          only: dp
  implicit none
  integer :: AllocateStatus

  ALLOCATE (r_min(1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (r_min_sq(1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (r_min_tab(1:nAtomTypes, 1:nAtomTypes), STAT = AllocateStatus)

  ALLOCATE (sig_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (q_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (Mu_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (sigsq_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (epAlpha_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (q_Nonzero(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
     
  sig_tab      = 0E0_dp
  q_tab        = 0E0_dp
  Mu_tab       = 0
  sigsq_tab    = 0E0_dp
  epAlpha_tab  = 0E0_dp
  r_min        = 0E0_dp
  r_min_sq     = 0E0_dp
  r_min_tab    = 0E0_dp

  ALLOCATE (totalMass(1:nMolTypes), STAT = AllocateStatus)

  ALLOCATE (nAtoms(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nIntraNonBond(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nBonds(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nAngles(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nTorsional(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nImproper(1:nMolTypes), STAT = AllocateStatus)

  nAtoms        = 0
  nIntraNonBond = 0
  nBonds        = 0
  nAngles       = 0
  nTorsional    = 0
  nImproper     = 0

  ALLOCATE (acptTrans(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (atmpTrans(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (acptRot(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (atmpRot(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (acptSwapIn(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (atmpSwapIn(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (acptSwapOut(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (atmpSwapOut(1:nMolTypes), STAT = AllocateStatus)

  ALLOCATE (acptCBMC(1:nMolTypes,1:nCBMCmax), STAT = AllocateStatus)
  ALLOCATE (atmpCBMC(1:nMolTypes,1:nCBMCmax), STAT = AllocateStatus)
  ALLOCATE (acptFECBMC(1:nMolTypes,1:nFECBMCmax), STAT = AllocateStatus)
  ALLOCATE (atmpFECBMC(1:nMolTypes,1:nFECBMCmax), STAT = AllocateStatus)

  ALLOCATE (max_dist(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (max_dist_single(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (max_rot(1:nMolTypes), STAT = AllocateStatus)

  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
end subroutine Allocate_Mpipi
!================================================================ 
subroutine Allocate_HPS_single
  use SimParameters,              only: nMolTypes, max_dist, max_dist_single, max_rot
  use AcceptRates,               only: acptTrans, atmpTrans, acptRot, atmpRot, acptSwapIn, &
                            atmpSwapIn, acptSwapOut, atmpSwapOut, acptCBMC, atmpCBMC, &
                                       acptFECBMC, atmpFECBMC, nCBMCmax, nFECBMCmax
  use ForceField,                only: totalMass, nAtoms, nIntraNonBond, nBonds, nAngles, &
                          r_min_sq, nTorsional, r_min, r_min_tab, nAtomTypes, nImproper
  use ForceFieldPara_HPS_single, only: sig_tab, sigsq_tab, eps_tab, lambda_tab, &
                                       q_tab, cutoffNB_tab, cutoffNBsq_tab, q_Nonzero
  use VarPrecision,              only: dp
  implicit none
  integer :: AllocateStatus

  ALLOCATE (r_min(1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (r_min_sq(1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (r_min_tab(1:nAtomTypes, 1:nAtomTypes), STAT = AllocateStatus) 

  ALLOCATE (sig_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (sigsq_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (eps_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (lambda_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (q_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (cutoffNB_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (cutoffNBsq_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (q_Nonzero(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
     
  sig_tab        = 0.0_dp  
  sigsq_tab      = 0.0_dp 
  eps_tab        = 0.0_dp
  lambda_tab     = 0.0_dp
  q_tab          = 0.0_dp
  cutoffNB_tab   = 0.0_dp
  cutoffNBsq_tab = 0.0_dp

  r_min      = 0.0_dp
  r_min_sq   = 0.0_dp
  r_min_tab  = 0.0_dp

  ALLOCATE (totalMass(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nAtoms(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nIntraNonBond(1:nMolTypes), STAT = AllocateStatus)      
  ALLOCATE (nBonds(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nAngles(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nTorsional(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nImproper(1:nMolTypes), STAT = AllocateStatus)

  nAtoms        = 0
  nIntraNonBond = 0
  nBonds        = 0
  nAngles       = 0
  nTorsional    = 0
  nImproper     = 0

  ALLOCATE (acptTrans(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (atmpTrans(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (acptRot(1:nMolTypes),   STAT = AllocateStatus)
  ALLOCATE (atmpRot(1:nMolTypes),   STAT = AllocateStatus)
  ALLOCATE (acptSwapIn(1:nMolTypes),  STAT = AllocateStatus)
  ALLOCATE (atmpSwapIn(1:nMolTypes),  STAT = AllocateStatus)
  ALLOCATE (acptSwapOut(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (atmpSwapOut(1:nMolTypes), STAT = AllocateStatus)

  ALLOCATE (acptCBMC(1:nMolTypes,1:nCBMCmax), STAT = AllocateStatus)
  ALLOCATE (atmpCBMC(1:nMolTypes,1:nCBMCmax), STAT = AllocateStatus)
  ALLOCATE (acptFECBMC(1:nMolTypes,1:nFECBMCmax), STAT = AllocateStatus)
  ALLOCATE (atmpFECBMC(1:nMolTypes,1:nFECBMCmax), STAT = AllocateStatus)

  ALLOCATE (max_dist(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (max_dist_single(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (max_rot(1:nMolTypes), STAT = AllocateStatus) 

  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

end subroutine Allocate_HPS_single
!================================================================ 
subroutine Allocate_HPS_piecewise
  use SimParameters,           only: nMolTypes, max_dist, max_dist_single, max_rot
  use AcceptRates,            only: acptTrans, atmpTrans, acptRot, atmpRot, acptSwapIn, &
                              atmpSwapIn, acptSwapOut, atmpSwapOut, acptCBMC, atmpCBMC, &
                                acptFECBMC, atmpFECBMC, nCBMCmax, nFECBMCmax
  use ForceField,             only: totalMass, nAtoms, nIntraNonBond, nBonds, nAngles, &
                              r_min_sq, nTorsional, r_min, r_min_tab, nAtomTypes, nImproper
  use ForceFieldPara_HPS_piecewise, only: sig_tab, sigsq_tab, eps_tab, lambda_tab, &
                                          q_tab, cutoffNB_tab, cutoffNBsq_tab, q_Nonzero
  use VarPrecision,                only: dp
  implicit none
  integer :: AllocateStatus

  ALLOCATE (r_min(1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (r_min_sq(1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (r_min_tab(1:nAtomTypes, 1:nAtomTypes), STAT = AllocateStatus)

  ALLOCATE (sig_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (sigsq_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (eps_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (lambda_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (q_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (cutoffNB_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (cutoffNBsq_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (q_Nonzero(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
     
  sig_tab        = 0.0_dp
  sigsq_tab      = 0.0_dp
  eps_tab        = 0.0_dp
  lambda_tab     = 0.0_dp
  q_tab          = 0.0_dp
  cutoffNB_tab   = 0.0_dp
  cutoffNBsq_tab = 0.0_dp
  r_min          = 0.0_dp
  r_min_sq       = 0.0_dp
  r_min_tab      = 0.0_dp

  ALLOCATE (totalMass(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nAtoms(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nIntraNonBond(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nBonds(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nAngles(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nTorsional(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nImproper(1:nMolTypes), STAT = AllocateStatus)

  nAtoms        = 0
  nIntraNonBond = 0
  nBonds        = 0
  nAngles       = 0
  nTorsional    = 0
  nImproper     = 0

  ALLOCATE (acptTrans(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (atmpTrans(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (acptRot(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (atmpRot(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (acptSwapIn(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (atmpSwapIn(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (acptSwapOut(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (atmpSwapOut(1:nMolTypes), STAT = AllocateStatus)

  ALLOCATE (acptCBMC(1:nMolTypes,1:nCBMCmax), STAT = AllocateStatus)
  ALLOCATE (atmpCBMC(1:nMolTypes,1:nCBMCmax), STAT = AllocateStatus)
  ALLOCATE (acptFECBMC(1:nMolTypes,1:nFECBMCmax), STAT = AllocateStatus)
  ALLOCATE (atmpFECBMC(1:nMolTypes,1:nFECBMCmax), STAT = AllocateStatus)

  ALLOCATE (max_dist(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (max_dist_single(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (max_rot(1:nMolTypes), STAT = AllocateStatus)

  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

end subroutine Allocate_HPS_piecewise
!================================================================ 
subroutine Allocate_HPS_cation_pi
  use SimParameters,    only: nMolTypes, max_dist, max_dist_single, max_rot
  use AcceptRates,      only: acptTrans, atmpTrans, acptRot, atmpRot, acptSwapIn, &
                        atmpSwapIn, acptSwapOut, atmpSwapOut, acptCBMC, atmpCBMC, &
                        acptFECBMC, atmpFECBMC, nCBMCmax, nFECBMCmax
  use ForceField,       only: totalMass, nAtoms, nIntraNonBond, nBonds, nAngles, &
                        r_min_sq, nTorsional, r_min, r_min_tab, nAtomTypes, nImproper
  use ForceFieldPara_HPS_cation_pi, only: sig_tab, sigsq_tab, eps_tab, lambda_tab, &
                        q_tab, cutoffNB_tab, cutoffNBsq_tab, epsLJ_tab, q_Nonzero
  use VarPrecision,     only: dp
  implicit none
  integer :: AllocateStatus

  ALLOCATE (r_min(1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (r_min_sq(1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (r_min_tab(1:nAtomTypes, 1:nAtomTypes), STAT = AllocateStatus)

  ALLOCATE (sig_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (sigsq_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (eps_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (lambda_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (q_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (cutoffNB_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (cutoffNBsq_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (epsLJ_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
  ALLOCATE (q_Nonzero(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
     
  sig_tab = 0.0_dp
  sigsq_tab = 0.0_dp
  eps_tab = 0.0_dp
  lambda_tab = 0.0_dp
  q_tab = 0.0_dp
  cutoffNB_tab = 0.0_dp
  cutoffNBsq_tab = 0.0_dp
  epsLJ_tab = 0.0_dp

  r_min = 0.0_dp
  r_min_sq = 0.0_dp
  r_min_tab = 0.0_dp

  ALLOCATE (totalMass(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nAtoms(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nIntraNonBond(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nBonds(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nAngles(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nTorsional(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (nImproper(1:nMolTypes), STAT = AllocateStatus)

  nAtoms = 0
  nIntraNonBond = 0
  nBonds = 0
  nAngles = 0
  nTorsional = 0
  nImproper = 0

  ALLOCATE (acptTrans(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (atmpTrans(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (acptRot(1:nMolTypes),   STAT = AllocateStatus)
  ALLOCATE (atmpRot(1:nMolTypes),  STAT = AllocateStatus)
  ALLOCATE (acptSwapIn(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (atmpSwapIn(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (acptSwapOut(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (atmpSwapOut(1:nMolTypes), STAT = AllocateStatus)

  ALLOCATE (acptCBMC(1:nMolTypes,1:nCBMCmax), STAT = AllocateStatus)
  ALLOCATE (atmpCBMC(1:nMolTypes,1:nCBMCmax), STAT = AllocateStatus)
  ALLOCATE (acptFECBMC(1:nMolTypes,1:nFECBMCmax), STAT = AllocateStatus)
  ALLOCATE (atmpFECBMC(1:nMolTypes,1:nFECBMCmax), STAT = AllocateStatus)

  ALLOCATE (max_dist(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (max_dist_single(1:nMolTypes), STAT = AllocateStatus)
  ALLOCATE (max_rot(1:nMolTypes), STAT = AllocateStatus)

  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

end subroutine Allocate_HPS_cation_pi
!================================================================ 
      subroutine Allocate_Pedone
      use SimParameters
      use ForceField
      use ForceFieldPara_Pedone
      use AcceptRates
      implicit none
      integer :: AllocateStatus
      
      ALLOCATE (pedoneData(1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (atomData(1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (nAtoms(1:nMolTypes), STAT = AllocateStatus)

      nAtoms = 1

      ALLOCATE (r_min(1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (r_min_sq(1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (r_min_tab(1:nAtomTypes, 1:nAtomTypes), STAT = AllocateStatus) 
      ALLOCATE (bornRad(1:nAtomTypes), STAT = AllocateStatus)


      ALLOCATE (alpha_Tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (D_Tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (repul_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (rEq_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (q_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)

      ALLOCATE (totalMass(1:nAtomTypes), STAT = AllocateStatus)
     
      repul_tab = 0d0
      D_Tab = 0d0       
      alpha_Tab = 0d0   
      rEq_tab = 0d0
      q_tab = 0d0
      r_min_tab = 0d0

      ALLOCATE (acptTrans(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (atmpTrans(1:nMolTypes), STAT = AllocateStatus)
      acptTrans = 0E0_dp
      atmpTrans = 1E-40_dp


      ALLOCATE (acptRot(1:nMolTypes),   STAT = AllocateStatus)
      ALLOCATE (atmpRot(1:nMolTypes),  STAT = AllocateStatus)
      acptRot = 0E0_dp
      atmpRot = 1E-40_dp

      ALLOCATE (acptSwapIn(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (atmpSwapIn(1:nMolTypes), STAT = AllocateStatus)
      acptSwapIn = 0E0_dp
      atmpSwapIn = 1E-40_dp

      ALLOCATE (acptSwapOut(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (atmpSwapOut(1:nMolTypes), STAT = AllocateStatus)
      acptSwapOut = 0E0_dp
      atmpSwapOut = 1E-40_dp


      ALLOCATE (acptInSize(1:maxMol), STAT = AllocateStatus)
      ALLOCATE (atmpInSize(1:maxMol), STAT = AllocateStatus)

      acptInSize = 0E0_dp
      atmpInSize = 1E-100_dp

      ALLOCATE (max_dist(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (max_dist_single(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (max_rot(1:nMolTypes), STAT = AllocateStatus) 

      
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  
      end subroutine

!================================================================ 
      subroutine Allocate_Tersoff
      use SimParameters
      use ForceField
      use ForceFieldPara_Tersoff
      use AcceptRates
      implicit none
      integer :: AllocateStatus
      
      ALLOCATE (tersoffData(1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (atomData(1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (nAtoms(1:nMolTypes), STAT = AllocateStatus)

      nAtoms = 1

      ALLOCATE (r_min(1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (r_min_sq(1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (r_min_tab(1:nAtomTypes, 1:nAtomTypes), STAT = AllocateStatus) 
      ALLOCATE (totalMass(1:nAtomTypes), STAT = AllocateStatus)

      ALLOCATE (acptTrans(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (atmpTrans(1:nMolTypes), STAT = AllocateStatus)
      acptTrans = 0E0_dp
      atmpTrans = 1E-40_dp


      ALLOCATE (acptRot(1:nMolTypes),   STAT = AllocateStatus)
      ALLOCATE (atmpRot(1:nMolTypes),  STAT = AllocateStatus)
      acptRot = 0E0_dp
      atmpRot = 1E-40_dp

      ALLOCATE (acptSwapIn(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (atmpSwapIn(1:nMolTypes), STAT = AllocateStatus)
      acptSwapIn = 0E0_dp
      atmpSwapIn = 1E-40_dp

      ALLOCATE (acptSwapOut(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (atmpSwapOut(1:nMolTypes), STAT = AllocateStatus)
      acptSwapOut = 0E0_dp
      atmpSwapOut = 1E-40_dp


      ALLOCATE (acptInSize(1:maxMol), STAT = AllocateStatus)
      ALLOCATE (atmpInSize(1:maxMol), STAT = AllocateStatus)

      acptInSize = 0E0_dp
      atmpInSize = 1E-100_dp

      ALLOCATE (max_dist(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (max_dist_single(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (max_rot(1:nMolTypes), STAT = AllocateStatus) 

      
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  
      end subroutine
!===================================================================================
subroutine Read_LJ_Q(lineStore)
  use SimParameters,             only: convEng, convDist, echoInput, DielecConst
  use ForceField,               only: atomData, nAtomTypes
  use ForceFieldPara_LJ_Q,      only: ep_tab, sig_tab, q_tab, ep_func, sig_func
  use ForceFieldFunctions,      only: MixRule, GeoMean_MixingFunc, Mean_MixingFunc
  use VarPrecision,             only: dp
  implicit none
  character(len=maxLineLen), intent(in) :: lineStore(:)
  integer :: i, j, iLine, nLines

  nLines = size(lineStore)
  i = 0

  !---------------------------------------------------------------
  ! Read atom types and parameters: name, symbol, epsilon, sigma,
  ! charge, and mass. Convert to internal units.
  !---------------------------------------------------------------
  do iLine = 2, nLines - 1
    i = i + 1
    read(lineStore(iLine),*) atomData(i)%atmName, atomData(i)%Symb, atomData(i)%ep, &
                             atomData(i)%sig, atomData(i)%q, atomData(i)%mass
    atomData(i)%ep  = atomData(i)%ep * convEng
    atomData(i)%sig = atomData(i)%sig * convDist
  enddo

  ! Set mixing rules for LJ and charge parameters
  ep_func  => GeoMean_MixingFunc
  sig_func => Mean_MixingFunc

  !---------------------------------------------------------------
  ! Construct epsilon, sigma, and charge lookup tables between
  ! all atom type pairs using defined mixing rules.
  !---------------------------------------------------------------
  if(echoInput) then
    write(35,*) "---------------------------------------------"
    write(35,*) "Interaction Table"
    write(35,*) " i "," j ", " eps ", " sig ", " q "
  endif     

  do i = 1, nAtomTypes
    do j = i, nAtomTypes
      ep_tab(i,j)  = 4d0 * ep_func(atomData(i)%ep, atomData(j)%ep)
      sig_tab(i,j) = sig_func(atomData(i)%sig, atomData(j)%sig)**2
      q_tab(i,j)   = atomData(i)%q * atomData(j)%q * coulombConst / DielecConst

      ! Fill symmetric entries
      ep_tab(j,i)  = ep_tab(i,j)
      sig_tab(j,i) = sig_tab(i,j)
      q_tab(j,i)   = q_tab(i,j)

      ! Optional input echo
      if(echoInput) then
        write(35,*) i, j, ep_tab(i,j)/4d0, sqrt(sig_tab(i,j)), q_tab(i,j)
      endif
    enddo
  enddo

  if(echoInput) then
    write(35,*) "---------------------------------------------"
    flush(35)
  endif

end subroutine Read_LJ_Q
!===================================================================================
subroutine Read_Mpipi(lineStore)
  use SimParameters,        only: convEng, convDist, echoInput, outputEConv, &
                                  outputLenConv, DielecConst
  use ForceField,           only: atomData, nAtomTypes
  use ForceFieldPara_Mpipi, only: sig_tab, q_tab, Mu_tab, sigsq_tab, &
                                  epAlpha_tab, q_Nonzero, rcutElec, rcutElecsq
  use VarPrecision,         only: dp
  implicit none
  character(len=maxLineLen), intent(in) :: lineStore(:)
  integer :: i, j, k, nParameterSets, iLine, nLines, Nu, Mu
  real(dp) :: eps, sig, x, R, Alpha, cutElec

  nLines = size(lineStore)

  i = 0
  iLine = 1
  do k = 1, nAtomTypes
    i = i + 1
    iLine = iLine + 1
    read(lineStore(iLine),*) atomData(i)%atmName, atomData(i)%Symb, atomData(i)%ep, &
                             atomData(i)%sig, atomData(i)%q, atomData(i)%mass
    atomData(i)%ep = atomData(i)%ep * convEng
    atomData(i)%sig = atomData(i)%sig * convDist
  enddo

  if(echoInput) then
    write(35,*) "---------------------------------------------"
    write(35,*) "Interaction Table"
    write(35,*) " i "," j ", " eps ", " sig ", " q ", " R ", " Mu ", " Nu ", " Alpha ", " q_Nonzero "
  endif   

  nParameterSets = nAtomTypes * (nAtomTypes + 1) / 2
  q_Nonzero = .false.

  do k = 1, nParameterSets
    iLine = iLine + 1
    read(lineStore(iLine),*) i, j, eps, sig, Nu, Mu, cutElec
    eps = eps * convEng
    sig = sig * convDist
    cutElec = cutElec * convDist
    R = 3.0_dp * sig
    sigsq_tab(i,j) = sig * sig
    sig_tab(i,j) = sig
    Mu_tab(i,j) = Mu
    q_tab(i,j) = atomData(i)%q * atomData(j)%q * coulombConst / DielecConst
    if (abs(q_tab(i,j)) > 1.0E-2_dp) q_Nonzero(i,j) = .true.
    x = (R / sig)**(2 * Mu_tab(i,j))
    Alpha = 2.0_dp * Nu * x * (((2.0_dp * Nu + 1.0_dp) / (2.0_dp * Nu * (x - 1.0_dp)))**(2 * Nu + 1))
    epAlpha_tab(i,j) = eps * Alpha
    rcutElec = cutElec

    sigsq_tab(j,i) = sigsq_tab(i,j)
    sig_tab(j,i) = sig_tab(i,j)
    Mu_tab(j,i) = Mu_tab(i,j)
    epAlpha_tab(j,i) = epAlpha_tab(i,j)
    q_tab(j,i) = q_tab(i,j)
    q_Nonzero(j,i) = q_Nonzero(i,j)

    if(echoInput) then
      write(35,*) i, j, eps * outputEConv, sig * outputLenConv, q_tab(i,j), R * outputLenConv, &
                  Mu_tab(i,j), Nu, Alpha, q_Nonzero(i,j)
    endif
  enddo

  if(echoInput) then
    write(35,*) "---------------------------------------------"
    flush(35)
  endif

  rcutElecsq = rcutElec * rcutElec
end subroutine Read_Mpipi
!===================================================================================
subroutine Read_HPS_single(lineStore)
  use SimParameters,             only: convEng, convDist, echoInput, outputEConv, &
                                       outputLenConv, DielecConst
  use ForceField,                only: atomData, nAtomTypes
  use ForceFieldPara_HPS_single, only: sig_tab, sigsq_tab, eps_tab, lambda_tab, q_tab, &
                                       cutoffNB_tab, cutoffNBsq_tab, q_Nonzero, &
                                       rcutElec, rcutElecsq
  use VarPrecision,              only: dp
  implicit none
  character(len=maxLineLen), intent(in) :: lineStore(:)
  integer :: i, j, k, nParameterSets, iLine, nLines
  real(dp) :: eps, sig, lambda, cutNB, cutElec

  nLines = size(lineStore)

  i = 0
  iLine = 1
  do k = 1, nAtomTypes
    i = i + 1
    iLine = iLine + 1
    read(lineStore(iLine),*) atomData(i)%atmName, atomData(i)%Symb, atomData(i)%ep, &
                             atomData(i)%sig, atomData(i)%q, atomData(i)%mass
    atomData(i)%ep  = atomData(i)%ep * convEng
    atomData(i)%sig = atomData(i)%sig * convDist
  enddo

! Generate the look up tables for the inter molecular interactions
  if(echoInput) then
    write(35,*) "---------------------------------------------"
    write(35,*) "Interaction Table"
    write(35,*) " i "," j ", " eps ", " sig ", " q ", " lambda ", " rcutNB ", " rcutElec ", " q_Nonzero "
  endif 

  nParameterSets = nAtomTypes * (nAtomTypes+1) / 2
  q_Nonzero = .false.

  do k = 1, nParameterSets
    iLine = iLine + 1
    read(lineStore(iLine),*) i, j, eps, sig, lambda, cutNB, cutElec
    eps     = eps     * convEng
    sig     = sig     * convDist
    cutNB   = cutNB   * convDist
    cutElec = cutElec * convDist

    eps_tab(i,j)        = eps
    sigsq_tab(i,j)      = sig * sig
    sig_tab(i,j)        = sig
    lambda_tab(i,j)     = lambda
    cutoffNB_tab(i,j)   = cutNB
    cutoffNBsq_tab(i,j) = cutNB * cutNB
    rcutElec            = cutElec
    rcutElecsq          = rcutElec * rcutElec
    q_tab(i,j)          = atomData(i)%q * atomData(j)%q * coulombConst / DielecConst
    if (abs(q_tab(i,j)) > 1.0e-2_dp) q_Nonzero(i,j) = .true.

    sigsq_tab(j,i)      = sigsq_tab(i,j)
    sig_tab(j,i)        = sig_tab(i,j)
    eps_tab(j,i)        = eps_tab(i,j)
    lambda_tab(j,i)     = lambda_tab(i,j)
    cutoffNB_tab(j,i)   = cutoffNB_tab(i,j)
    cutoffNBsq_tab(j,i) = cutoffNBsq_tab(i,j)
    q_tab(j,i)          = q_tab(i,j)
    q_Nonzero(j,i)      = q_Nonzero(i,j)

    if(echoInput) then
      write(35,*) i, j, eps * outputEConv, sig * outputLenConv, q_tab(i,j), lambda, &
                  cutNB * outputLenConv, cutElec * outputLenConv, q_Nonzero(i,j)
    endif
  enddo

  if(echoInput) then
    write(35,*) "---------------------------------------------"
    flush(35)
  endif

end subroutine Read_HPS_single
!===================================================================================
subroutine Read_HPS_piecewise(lineStore)
  use SimParameters,                only: convEng, convDist, echoInput, outputEConv, &
                                          outputLenConv, DielecConst
  use ForceField,                  only: atomData, nAtomTypes
  use ForceFieldPara_HPS_piecewise, only: sig_tab, sigsq_tab, eps_tab, lambda_tab, q_tab, &
                                          cutoffNB_tab, cutoffNBsq_tab, q_Nonzero, &
                                          rcutElec, rcutElecsq
  use VarPrecision,                only: dp
  implicit none
  character(len=maxLineLen), intent(in) :: lineStore(:)
  integer :: i, j, k, nParameterSets, iLine, nLines
  real(dp) :: eps, sig, lambda, cutNB, cutElec

  nLines = size(lineStore)

  i = 0
  iLine = 1
  do k = 1, nAtomTypes
    i = i + 1
    iLine = iLine + 1
    read(lineStore(iLine),*) atomData(i)%atmName, atomData(i)%Symb, atomData(i)%ep, &
                             atomData(i)%sig, atomData(i)%q, atomData(i)%mass
    atomData(i)%ep  = atomData(i)%ep * convEng
    atomData(i)%sig = atomData(i)%sig * convDist
  enddo

  if(echoInput) then
    write(35,*) "---------------------------------------------"
    write(35,*) "Interaction Table"
    write(35,*) " i "," j ", " eps ", " sig ", " q ", " lambda ", " rcutNB ", " rcutElec ", " q_Nonzero "
  endif

  nParameterSets = nAtomTypes * (nAtomTypes + 1) / 2
  q_Nonzero = .false.
  do k = 1, nParameterSets
    iLine = iLine + 1
    read(lineStore(iLine),*) i, j, eps, sig, lambda, cutNB, cutElec
    eps     = eps * convEng
    sig     = sig * convDist
    cutNB   = cutNB * convDist
    cutElec = cutElec * convDist

    eps_tab(i,j)        = eps
    sigsq_tab(i,j)      = sig * sig
    sig_tab(i,j)        = sig
    lambda_tab(i,j)     = lambda
    cutoffNB_tab(i,j)   = cutNB
    cutoffNBsq_tab(i,j) = cutNB * cutNB
    rcutElec            = cutElec
    rcutElecsq          = rcutElec * rcutElec
    q_tab(i,j)          = atomData(i)%q * atomData(j)%q * coulombConst / DielecConst
    if (abs(q_tab(i,j)) .gt. 1.0e-2_dp) q_Nonzero(i,j) = .true.

    sigsq_tab(j,i)      = sigsq_tab(i,j)
    sig_tab(j,i)        = sig_tab(i,j)
    eps_tab(j,i)        = eps_tab(i,j)
    lambda_tab(j,i)     = lambda_tab(i,j)
    cutoffNB_tab(j,i)   = cutoffNB_tab(i,j)
    cutoffNBsq_tab(j,i) = cutoffNBsq_tab(i,j)
    q_tab(j,i)          = q_tab(i,j)
    q_Nonzero(j,i)      = q_Nonzero(i,j)

    if(echoInput) then
      write(35,*) i, j, eps * outputEConv, sig * outputLenConv, q_tab(i,j), lambda, &
                  cutNB * outputLenConv, cutElec * outputLenConv, q_Nonzero(i,j)
    endif
  enddo

  if(echoInput) then
    write(35,*) "---------------------------------------------"
    flush(35)
  endif

end subroutine Read_HPS_piecewise
!===================================================================================
subroutine Read_HPS_cation_pi(lineStore)
  use SimParameters,                  only: convEng, convDist, echoInput, outputEConv, &
                                            outputLenConv, DielecConst
  use ForceField,                    only: atomData, nAtomTypes
  use ForceFieldPara_HPS_cation_pi, only: sig_tab, sigsq_tab, eps_tab, lambda_tab, q_tab, &
                                            cutoffNB_tab, cutoffNBsq_tab, epsLJ_tab, q_Nonzero, &
                                            rcutElec, rcutElecsq
  use VarPrecision,                  only: dp

  implicit none
  character(len=maxLineLen), intent(in) :: lineStore(:)
  integer :: i, j, k, nParameterSets, iLine, nLines
  real(dp) :: eps, sig, lambda, cutNB, cutElec

  nLines = size(lineStore)
  i = 0
  iLine = 1

  do k = 1, nAtomTypes
    i = i + 1
    iLine = iLine + 1
    read(lineStore(iLine),*) atomData(i)%atmName, atomData(i)%Symb, atomData(i)%ep, &
                             atomData(i)%sig, atomData(i)%q, atomData(i)%mass
    atomData(i)%ep = atomData(i)%ep * convEng
    atomData(i)%sig = atomData(i)%sig * convDist
  enddo

  if(echoInput) then
    write(35,*) "---------------------------------------------"
    write(35,*) "Interaction Table"
    write(35,*) " i "," j ", " eps ", " sig ", " q ", " lambda ", " rcutNB ", " rcutElec ", " q_Nonzero "
  endif
     
  nParameterSets = nAtomTypes * (nAtomTypes+1) / 2
  q_Nonzero = .false.

  do k = 1, nParameterSets
    iLine = iLine + 1
    read(lineStore(iLine),*) i, j, eps, sig, lambda, cutNB, cutElec
    eps = eps * convEng
    sig = sig * convDist
    cutNB = cutNB * convDist
    cutElec = cutElec * convDist
    eps_tab(i,j) = eps
    sigsq_tab(i,j) = sig * sig
    sig_tab(i,j) = sig
    lambda_tab(i,j) = lambda
    cutoffNB_tab(i,j) = cutNB
    cutoffNBsq_tab(i,j) = cutNB * cutNB
    rcutElec = cutElec
    rcutElecsq = rcutElec * rcutElec
    q_tab(i,j) = atomData(i)%q * atomData(j)%q * coulombConst / DielecConst
    if (abs(q_tab(i,j)) .gt. 1E-2_dp) q_Nonzero(i,j) = .true.

    sigsq_tab(j,i) = sigsq_tab(i,j)
    sig_tab(j,i) = sig_tab(i,j)
    eps_tab(j,i) = eps_tab(i,j)
    lambda_tab(j,i) = lambda_tab(i,j)
    cutoffNB_tab(j,i) = cutoffNB_tab(i,j)
    cutoffNBsq_tab(j,i) = cutoffNBsq_tab(i,j)
    q_tab(j,i) = q_tab(i,j)
    q_Nonzero(j,i) = q_Nonzero(i,j)

    if(echoInput) then
      write(35,*) i,j, eps * outputEConv, sig * outputLenConv, q_tab(i,j), lambda, &
                  cutNB * outputLenConv, cutElec * outputLenConv, q_Nonzero(i,j)
    endif
  enddo

  if(echoInput) then
    write(35,*) "---------------------------------------------"
    flush(35)
  endif

  do k = 1, 24
    iLine = iLine + 1
    read(lineStore(iLine),*) i, j, eps, sig
    eps = eps * convEng
    sig = sig * convDist
    epsLJ_tab(i,j) = eps
    epsLJ_tab(j,i) = epsLJ_tab(i,j)
  enddo

end subroutine Read_HPS_cation_pi
!===================================================================================
      subroutine Read_Pedone(lineStore)
      use SimParameters
      use ForceField
      use ForceFieldPara_Pedone
      implicit none
      character(len=maxLineLen), intent(in) :: lineStore(:)
      integer :: i, j, iLine, nLines

      nLines = size(lineStore)

      bornRad = 0E0_dp

      i = 0
      do iLine = 2, nLines-1
        i = i + 1
!        write(*,*)  lineStore(iLine)
        read(lineStore(iLine),*) pedoneData(i)%atmName, pedoneData(i)%Symb, pedoneData(i)%repul, pedoneData(i)%rEq, &
                                 pedoneData(i)%alpha, pedoneData(i)%delta, pedoneData(i)%q, pedoneData(i)%mass, bornRad(i)
      enddo
      
      do i = 1, nAtomTypes
        if(echoInput) then
          write(35,*) pedoneData(i)%atmName, pedoneData(i)%Symb, pedoneData(i)%repul, pedoneData(i)%rEq, &
                    pedoneData(i)%alpha, pedoneData(i)%delta, pedoneData(i)%q, pedoneData(i)%mass, bornRad(i)
        endif
        pedoneData(i)%repul = pedoneData(i)%repul * convEng
        pedoneData(i)%delta = pedoneData(i)%delta * convEng
        pedoneData(i)%rEq = pedoneData(i)%rEq * convDist   
        bornRad(i) = bornRad(i) * convDist   
        atomData(i)%atmName = pedoneData(i)%atmName
        atomData(i)%Symb = pedoneData(i)%Symb
      enddo

      q_tab = 0d0
      do i = 1,nAtomTypes
        do j = i,nAtomTypes
          q_tab(i,j) = pedoneData(i)%q * pedoneData(j)%q * coulombConst
          q_tab(j,i) = q_tab(i,j)
          if(echoInput) then
             write(35,*) i, j, q_tab(i,j) 
           endif          
        enddo
      enddo
 
      repul_tab = 0d0
      alpha_Tab = 0d0
      rEq_tab = 0d0
      D_Tab = 0d0
      do i = 1,nAtomTypes
        repul_tab(i,1) = pedoneData(i)%repul
        alpha_Tab(i,1) = pedoneData(i)%alpha
        rEq_tab(i,1) = pedoneData(i)%rEq
        D_Tab(i,1) = pedoneData(i)%delta

        repul_tab(1,i) = repul_tab(i,1)
        alpha_Tab(1,i) = alpha_Tab(i,1)
        rEq_tab(1,i) = rEq_tab(i,1)
        D_Tab(1,i) = D_Tab(i,1)
         
      enddo

   
      if(echoInput) then
        do i = 1, nAtomTypes 
          do j = 1, nAtomTypes 
            write(35,*) i, j, repul_tab(i,j), alpha_Tab(i,j), rEq_tab(i,j), D_Tab(i,j), q_tab(i,j)
          enddo
        enddo
      endif 


      if(any(bornRad .ne. 0E0_dp)) then
        implcSolvent = .true.
      else
        implcSolvent = .false.
      endif

      end subroutine

!===================================================================================
      subroutine Read_Tersoff(lineStore)
      use SimParameters
      use ForceField
      use ForceFieldPara_Tersoff
      implicit none
      character(len=maxLineLen), intent(in) :: lineStore(:)
      integer :: i, iLine, nLines

      nLines = size(lineStore)

      i = 0
      do iLine = 2, nLines-1
        i = i + 1
!        write(*,*)  lineStore(iLine)
        read(lineStore(iLine),*) tersoffData(i)%atmName, tersoffData(i)%Symb, tersoffData(i)%A, & 
                                 tersoffData(i)%B, tersoffData(i)%c, tersoffData(i)%d, tersoffData(i)%n, &
                                 tersoffData(i)%lam1, tersoffData(i)%lam2, tersoffData(i)%h, tersoffData(i)%R, &
                                 tersoffData(i)%D2, tersoffData(i)%beta, tersoffData(i)%mass
      enddo


      do i = 1, nAtomTypes
        if(echoInput) then
          write(35,*) tersoffData(i)%atmName, tersoffData(i)%Symb, tersoffData(i)%A, & 
                                 tersoffData(i)%B, tersoffData(i)%c, tersoffData(i)%d, tersoffData(i)%n, &
                                 tersoffData(i)%lam1, tersoffData(i)%lam2, tersoffData(i)%h, tersoffData(i)%R, &
                                 tersoffData(i)%D2, tersoffData(i)%beta, tersoffData(i)%mass
        endif
        tersoffData(i)%A = tersoffData(i)%A * convEng
        tersoffData(i)%B = tersoffData(i)%B * convEng

        tersoffData(i)%lam1 = tersoffData(i)%lam1 / convDist
        tersoffData(i)%lam2 = tersoffData(i)%lam2 / convDist

        tersoffData(i)%R = tersoffData(i)%R * convDist   
        tersoffData(i)%D2 = tersoffData(i)%D2 * convDist

        atomData(i)%atmName = tersoffData(i)%atmName
        atomData(i)%Symb = tersoffData(i)%Symb
      enddo


      end subroutine
!===================================================================================
      subroutine LJ_SetFlags
      use SimParameters
      use Coords
      use ForceField
      use ForceFieldPara_LJ_Q
      use PairStorage, only: SetStorageFlags, rPair, useDistStore
      implicit none
      integer :: iType, iMol, iAtom
      integer :: jType, jMol, jAtom
      integer :: atmType1, atmType2, globIndx1, globIndx2
      real(dp) :: q, ep
  
      call IntegrateBendAngleProb

      if(useDistStore) then
       call SetStorageFlags(q_tab) 
         do iType = 1, nMolTypes
           do iMol = 1, NMAX(iType)
             do iAtom = 1, nAtoms(iType)
              atmType1 = atomArray(iType,iAtom)
              globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)            
              do jType = 1, nMolTypes
                do jMol = 1, NMAX(jType)
                  do jAtom = 1, nAtoms(jType)
                    atmType2 = atomArray(jType, jAtom)
                    globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
                    ep = ep_tab(atmType1,atmType2) 
                    q = q_tab(atmType1,atmType2)
                    rPair(globIndx1, globIndx2)%p%usePair = .true.
                    if(abs(q) .gt. 1E-16_dp) then
                      rPair(globIndx1, globIndx2)%p%storeRValue = .true.
                    endif
                    if(ep .eq. 0E0_dp) then
                      if(q .eq. 0E0_dp) then
                        rPair(globIndx1, globIndx2)%p%usePair = .false.

                      endif
                    endif 
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      endif

      end subroutine
!===================================================================================
      subroutine Mpipi_SetFlags
      use SimParameters
      use Coords
      use ForceField
      use ForceFieldPara_Mpipi
      use PairStorage, only: SetStorageFlags, rPair, useDistStore
      implicit none
      integer :: iType, iMol, iAtom
      integer :: jType, jMol, jAtom
      integer :: atmType1, atmType2, globIndx1, globIndx2
      real(dp) :: q, ep
  

      if(useDistStore) then
       call SetStorageFlags(q_tab) 
         do iType = 1, nMolTypes
           do iMol = 1, NMAX(iType)
             do iAtom = 1, nAtoms(iType)
              atmType1 = atomArray(iType,iAtom)
              globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)            
              do jType = 1, nMolTypes
                do jMol = 1, NMAX(jType)
                  do jAtom = 1, nAtoms(jType)
                    atmType2 = atomArray(jType, jAtom)
                    globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
                    ep = epAlpha_tab(atmType1,atmType2) 
                    q = q_tab(atmType1,atmType2)
                    rPair(globIndx1, globIndx2)%p%usePair = .true.
                    if(abs(q) .gt. 1E-16_dp) then
                      rPair(globIndx1, globIndx2)%p%storeRValue = .true.
                    endif
                    if(ep .eq. 0E0_dp) then
                      if(q .eq. 0E0_dp) then
                        rPair(globIndx1, globIndx2)%p%usePair = .false.

                      endif
                    endif 
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      endif

      end subroutine
!===================================================================================
      subroutine HPS_single_SetFlags
      use SimParameters
      use Coords
      use ForceField
      use ForceFieldPara_HPS_single
      use PairStorage, only: SetStorageFlags, rPair, useDistStore
      implicit none
      integer :: iType, iMol, iAtom
      integer :: jType, jMol, jAtom
      integer :: atmType1, atmType2, globIndx1, globIndx2
      real(dp) :: q, ep
  

      if(useDistStore) then
       call SetStorageFlags(q_tab) 
         do iType = 1, nMolTypes
           do iMol = 1, NMAX(iType)
             do iAtom = 1, nAtoms(iType)
              atmType1 = atomArray(iType,iAtom)
              globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)            
              do jType = 1, nMolTypes
                do jMol = 1, NMAX(jType)
                  do jAtom = 1, nAtoms(jType)
                    atmType2 = atomArray(jType, jAtom)
                    globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
                    ep = eps_tab(atmType1,atmType2) 
                    q = q_tab(atmType1,atmType2)
                    rPair(globIndx1, globIndx2)%p%usePair = .true.
                    if(abs(q) .gt. 1E-16_dp) then
                      rPair(globIndx1, globIndx2)%p%storeRValue = .true.
                    endif
                    if(ep .eq. 0E0_dp) then
                      if(q .eq. 0E0_dp) then
                        rPair(globIndx1, globIndx2)%p%usePair = .false.

                      endif
                    endif 
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      endif

      end subroutine
!===================================================================================
      subroutine HPS_piecewise_SetFlags
      use SimParameters
      use Coords
      use ForceField
      use ForceFieldPara_HPS_piecewise
      use PairStorage, only: SetStorageFlags, rPair, useDistStore
      implicit none
      integer :: iType, iMol, iAtom
      integer :: jType, jMol, jAtom
      integer :: atmType1, atmType2, globIndx1, globIndx2
      real(dp) :: q, ep
  

      if(useDistStore) then
       call SetStorageFlags(q_tab) 
         do iType = 1, nMolTypes
           do iMol = 1, NMAX(iType)
             do iAtom = 1, nAtoms(iType)
              atmType1 = atomArray(iType,iAtom)
              globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)            
              do jType = 1, nMolTypes
                do jMol = 1, NMAX(jType)
                  do jAtom = 1, nAtoms(jType)
                    atmType2 = atomArray(jType, jAtom)
                    globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
                    ep = eps_tab(atmType1,atmType2) 
                    q = q_tab(atmType1,atmType2)
                    rPair(globIndx1, globIndx2)%p%usePair = .true.
                    if(abs(q) .gt. 1E-16_dp) then
                      rPair(globIndx1, globIndx2)%p%storeRValue = .true.
                    endif
                    if(ep .eq. 0E0_dp) then
                      if(q .eq. 0E0_dp) then
                        rPair(globIndx1, globIndx2)%p%usePair = .false.

                      endif
                    endif 
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      endif

      end subroutine
!===================================================================================
      subroutine HPS_cation_pi_SetFlags
      use SimParameters
      use Coords
      use ForceField
      use ForceFieldPara_HPS_cation_pi
      use PairStorage, only: SetStorageFlags, rPair, useDistStore
      implicit none
      integer :: iType, iMol, iAtom
      integer :: jType, jMol, jAtom
      integer :: atmType1, atmType2, globIndx1, globIndx2
      real(dp) :: q, ep
  

      if(useDistStore) then
       call SetStorageFlags(q_tab) 
         do iType = 1, nMolTypes
           do iMol = 1, NMAX(iType)
             do iAtom = 1, nAtoms(iType)
              atmType1 = atomArray(iType,iAtom)
              globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)            
              do jType = 1, nMolTypes
                do jMol = 1, NMAX(jType)
                  do jAtom = 1, nAtoms(jType)
                    atmType2 = atomArray(jType, jAtom)
                    globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
                    ep = eps_tab(atmType1,atmType2) 
                    q = q_tab(atmType1,atmType2)
                    rPair(globIndx1, globIndx2)%p%usePair = .true.
                    if(abs(q) .gt. 1E-16_dp) then
                      rPair(globIndx1, globIndx2)%p%storeRValue = .true.
                    endif
                    if(ep .eq. 0E0_dp) then
                      if(q .eq. 0E0_dp) then
                        rPair(globIndx1, globIndx2)%p%usePair = .false.

                      endif
                    endif 
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      endif

      end subroutine
!===================================================================================
      subroutine Pedone_SetFlags
      use SimParameters
      use ForceField
      use ForceFieldPara_Pedone
      use PairStorage, only: SetStorageFlags, useDistStore
      implicit none

      if(useDistStore) then
        call SetStorageFlags(q_tab) 
      endif

      end subroutine

!===================================================================================
      subroutine Tersoff_SetFlags
      use SimParameters
      use ForceField
      use ForceFieldPara_Pedone
      use PairStorage, only: TurnOnAllStorageFlags, useDistStore
      implicit none

      if(useDistStore) then
        call TurnOnAllStorageFlags
      endif

      end subroutine
!===================================================================================
      subroutine FindMax(lineStore, targetCommand, commandMax)
      use SimParameters
      use ForceField
      use ForceFieldPara_LJ_Q
      implicit none
      character(len=maxLineLen), intent(in) :: lineStore(:)
      character(len=*), intent(in) :: targetCommand
      integer, intent(out) :: commandMax

      integer :: intValue
      character(len=25) :: curCommand, dummy

      integer :: iLine, nLines, lineStat

      nLines = size(lineStore)
      commandMax = 0
      do iLine = 1, nLines
        call GetCommand(lineStore(iLine), curCommand, lineStat)
        call LowerCaseLine(curCommand)
        if(trim(adjustl(curCommand)) .eq. trim(adjustl(targetCommand))) then
          read(lineStore(iLine), *) dummy, intValue
          commandMax = max(commandMax, intValue)
        endif
      enddo

      


      end subroutine

!================================================================================
      end module
