!======================================================
! Module: CoordinateTypes
! Purpose:
!   Defines data structures for representing atomic and molecular
!   coordinates, displacements, neighbor lists, and global indexing.
!
! Dependencies:
!   - Uses `VarPrecision` for consistent precision control.
!
! Author:
!   Aliasghar Sepehri
!======================================================
module CoordinateTypes
  use VarPrecision
  implicit none

  !------------------------------------------------------
  ! Pointer to a single double-precision real value
  type :: FloatPointer
    real(dp), pointer :: pnt
  end type FloatPointer

  !------------------------------------------------------
  ! Coordinates of atoms in a molecule (separate arrays for each axis)
  type :: SimpleMolCoords
    real(dp), allocatable :: x(:), y(:), z(:)
  end type SimpleMolCoords

  !------------------------------------------------------
  ! Coordinates of a single atom (structure of arrays alternative)
  type :: SimpleAtomCoords
    real(dp) :: x, y, z
  end type SimpleAtomCoords

  !------------------------------------------------------
  ! Global indexing structure to identify an atom across the system
  type :: GlobalAtomIndex
    integer(atomIntType) :: nType     ! Molecule type
    integer(atomIntType) :: nMol      ! Molecule index
    integer(atomIntType) :: nAtom     ! Atom index within molecule
    integer(atomIntType) :: atmType   ! Atom type (e.g., for force field lookup)
  end type GlobalAtomIndex

  !------------------------------------------------------
  ! A single molecule with coordinates and global atom indices
  type :: Molecule
    integer(atomIntType) :: indx                  ! Molecule index
    integer, allocatable :: globalIndx(:)         ! Global atom indices
    real(dp), allocatable :: x(:), y(:), z(:)     ! Atom positions
  end type Molecule

  !------------------------------------------------------
  ! Array wrapper for multiple Molecule objects
  type :: MolArrayType
    type(Molecule), allocatable :: mol(:)
  end type MolArrayType

  !------------------------------------------------------
  ! Pointers to molecule coordinates, used for referencing other copies
  type :: MolPointers
    integer(atomIntType) :: molType               ! Molecule type
    integer(atomIntType) :: molIndx               ! Molecule index
    type(FloatPointer), allocatable :: x(:), y(:), z(:) ! Pointers to positions
  end type MolPointers

  !------------------------------------------------------
  ! Information about a displacement of a single atom
  type :: Displacement
    integer(atomIntType) :: molType               ! Molecule type
    integer(atomIntType) :: atmIndx               ! Atom index within molecule
    integer(atomIntType) :: molIndx               ! Molecule index

    logical :: Displaced                          ! Whether atom was displaced

    real(dp) :: x_new, y_new, z_new               ! New trial coordinates
    real(dp), pointer :: x_old, y_old, z_old      ! Pointers to original coordinates
  end type Displacement

  !------------------------------------------------------
  ! Trial coordinates for all atoms in a molecule
  type :: TrialCoordinates
    integer(atomIntType) :: molType, molIndx
    real(dp), allocatable :: x(:), y(:), z(:)
  end type TrialCoordinates

! Types and variables for sparse neighbor list storage in a grand canonical ensemble
! nucleation simulation. Tracks atom pair indices for molecule pairs that are neighbors
! (NeighborList(i,j) = .true.) based on minDistCriteria for long biomolecules. Used for
! cluster detection and swap moves, optimized for memory efficiency.

! Neighbor details for a molecule pair, storing atom pair indices
  type :: NeighborDetails
    integer :: nPairs                    ! Number of atom pairs within minDistCriteria
    integer :: pairIndices(1:2000)  ! Atom indices [atm1_1, atm2_1, atm1_2, atm2_2, ...]
  end type NeighborDetails

! Molecule pair with neighbor details
  type :: MolPair
    type(NeighborDetails) :: details   ! Atom pair indices
  end type MolPair

end module CoordinateTypes
!==============================================================
! Module: AcceptRates
! Purpose:
!   Stores acceptance and attempt statistics for Monte Carlo moves,
!   including translation, rotation, insertion/removal, and CBMC-related moves.
!
!   These statistics are typically used for:
!     - Monitoring simulation performance
!     - Dynamically tuning move probabilities
!     - Reporting final acceptance ratios
!
! Dependencies:
!   - Uses `VarPrecision` for precision control
!
! Author:
!   Aliasghar Sepehri
!==============================================================
module AcceptRates
  use VarPrecision
  implicit none

  !----------------------------------------------------------
  ! Translation and rotation acceptance tracking
  real(dp), allocatable :: acptTrans(:)    ! Accepted translation moves
  real(dp), allocatable :: acptRot(:)      ! Accepted rotation moves
  real(dp), allocatable :: atmpTrans(:)    ! Attempted translation moves
  real(dp), allocatable :: atmpRot(:)      ! Attempted rotation moves

  !----------------------------------------------------------
  ! Grand Canonical Swap-In and Swap-Out (Insertion/Deletion)
  real(dp), allocatable :: acptSwapIn(:)   ! Accepted insertions
  real(dp), allocatable :: atmpSwapIn(:)   ! Attempted insertions
  real(dp), allocatable :: acptSwapOut(:)  ! Accepted deletions
  real(dp), allocatable :: atmpSwapOut(:)  ! Attempted deletions

  !----------------------------------------------------------
  ! Acceptance by inserted cluster size (for histogram-based biasing)
  real(dp), allocatable :: acptInSize(:)   ! Accepted insertions by cluster size
  real(dp), allocatable :: atmpInSize(:)   ! Attempted insertions by cluster size

  !----------------------------------------------------------
  ! Acceptance rates for angle/dihedral/distance generation (e.g., CBMC)
  real(dp) :: angGen_accpt     ! Accepted bond angle generation
  real(dp) :: angGen_atmp      ! Attempted bond angle generation
  real(dp) :: dihedGen_accpt   ! Accepted dihedral angle generation
  real(dp) :: dihedGen_atmp    ! Attempted dihedral angle generation
  real(dp) :: distGen_accpt    ! Accepted bond distance generation
  real(dp) :: distGen_atmp     ! Attempted bond distance generation

  !----------------------------------------------------------
  ! Rejection counter for cluster criterion violations
  real(dp) :: clusterCritRej   ! Rejections due to not satisfying cluster criterion

  !----------------------------------------------------------
  ! Maximum number of monomers allowed in CBMC and FECBMC chain growth
  integer :: nCBMCmax = 20     ! Max number of segments for CBMC
  integer :: nFECBMCmax = 10   ! Max number of segments for FECBMC

  !----------------------------------------------------------
  ! CBMC and FECBMC acceptance/attempt statistics (by chain length or molecule type)
  real(dp), allocatable :: acptCBMC(:,:)   ! Accepted CBMC moves
  real(dp), allocatable :: atmpCBMC(:,:)   ! Attempted CBMC moves

  real(dp), allocatable :: acptFECBMC(:,:) ! Accepted FECBMC moves
  real(dp), allocatable :: atmpFECBMC(:,:) ! Attempted FECBMC moves

end module AcceptRates
!===================================================================
! Module: AVBMC_RejectionVar
! Purpose:
!   Tracks rejection statistics for Aggregation-Volume-Bias Monte Carlo (AVBMC)
!   moves in nucleation simulations, especially for insertions and deletions
!   under species count constraints.
!
! Rejection Categories:
!   - ovrlapRej     : Overlap between particles (excluded volume violation)
!   - dbalRej       : Detailed balance violation
!   - critriaRej    : Cluster criterion not satisfied
!
!   Special handling for grand canonical insertion/deletion:
!     - boundaryRej     : Insertion rejected — exceeds max count for particle type
!     - boundaryRej_out : Deletion rejected — falls below min count for particle type
!
! Dependencies:
!   - Uses CoordinateTypes
!
! Author:
!   Aliasghar Sepehri
!===================================================================
module AVBMC_RejectionVar
  use CoordinateTypes
  implicit none

  !------------------------------------------------------
  ! Total rejections for AVBMC moves
  real(dp) :: totalRej             ! Total rejections (all causes)

  ! Specific rejection reasons for insertion moves
  real(dp) :: ovrlapRej            ! Rejected due to overlap
  real(dp) :: dbalRej              ! Rejected due to detailed balance
  real(dp) :: critriaRej           ! Rejected due to cluster criterion
  real(dp) :: boundaryRej          ! Rejected due to exceeding max particle count (insertion)

  ! Specific rejection reasons for deletion (outward) moves
  real(dp) :: totalRej_out         ! Total outward move rejections
  real(dp) :: dbalRej_out          ! Detailed balance violation (deletion)
  real(dp) :: critriaRej_out       ! Cluster criterion not satisfied (deletion)
  real(dp) :: boundaryRej_out      ! Rejected due to falling below min particle count (deletion)

end module AVBMC_RejectionVar
!======================================================
! Module: ParallelVar
! Purpose:
!   Stores global MPI variables used across all parallel routines.
!   These variables are essential for handling processor ranks,
!   communicator size, error flags, and MPI message tagging.
!
! Usage:
!   - Use this module with `use ParallelVar` in any MPI-parallel routine.
!   - Variables are set during MPI initialization and remain in scope.
!
! Author:
!   Aliasghar Sepehri
!======================================================
module ParallelVar
  implicit none

  !------------------------------------------------------
  ! MPI rank of the current process (0 to p_size - 1)
  integer :: myid

  ! Total number of processes in MPI_COMM_WORLD
  integer :: p_size

  ! Return flag from MPI calls
  integer :: ierror

  ! Indicates whether MPI was initialized with required OpenMP thread level
  integer :: provided

  ! Message tag used in MPI communication (can be varied for context)
  integer :: tag

  ! Output unit identifier for rank-specific files (e.g., nout = 100 + myid)
  integer :: nout

end module ParallelVar
!==============================================================
module Coords
  use CoordinateTypes
  implicit none

  !------------------------------------------------------
  ! Array of all molecules in the simulation
  ! Each entry contains the coordinates and global indices of atoms in one molecule
  type(MolArrayType), allocatable, target :: MolArray(:)

  !------------------------------------------------------
  ! List of global atom indices for the entire system
  type(GlobalAtomIndex), allocatable :: atomIndicies(:)

  !------------------------------------------------------
  ! Trial configurations for two new molecules (used in CBMC, AVBMC, etc.)
  type(TrialCoordinates) :: newMol      ! Trial molecule #1
  type(TrialCoordinates) :: newMol2     ! Trial molecule #2 
  !------------------------------------------------------
  ! Temporary coordinates for trial chain generations (Rosenbluth sampling)
  type(SimpleMolCoords), allocatable :: rosenTrial(:)

  !------------------------------------------------------
! Neighbor flags (existing, NeighborList(i,j) = .true. if molecules i,j are neighbors)
  logical, allocatable :: NeighborList(:,:)

  !------------------------------------------------------
! Sparse list of neighbor molecule pairs
  type(MolPair), allocatable :: NeighborPairs(:,:)
!      type(NeighVars2), allocatable :: NeighborDetails(:,:)
!      integer, allocatable :: NumNei(:)   
  !------------------------------------------------------
  ! Indexing and mapping utilities
  integer, allocatable :: typeList(:)      ! List of molecule types (used in trial generation)
  integer, allocatable :: subIndxList(:)   ! Sub-index within each molecule type

  !------------------------------------------------------
  ! Gas-phase configurations (e.g., for insertion move candidates or test volumes)
  type(SimpleMolCoords), allocatable :: gasConfig(:)

  !------------------------------------------------------
end module Coords
    
!==============================================================
! Module: CBMC_Variables
! Purpose:
!   Stores all topology-related and sampling variables required for
!   Configurational Bias Monte Carlo (CBMC) and Fixed-Endpoints CBMC (FECBMC)
!   simulations of flexible molecules (e.g., polymers, proteins).
!
!   Includes:
!     - Topological data (e.g., paths, bonding info)
!     - Sampling schedules and growth orders
!     - Dihedral angle biasing histograms
!     - CBMC swap-in/swap-out trial settings
!
! Dependencies:
!   - Constants (for 2π and other constants)
!   - VarPrecision (for dp and atomIntType)
!
! Author:
!   Aliasghar Sepehri
!===================================================================
module CBMC_Variables
  use Constants
  use VarPrecision
  implicit none

  !------------------------------------------------------
  ! Bonds associated with atoms in a molecule
  type :: BondNumber
    integer(atomIntType), allocatable :: atom(:)
  end type BondNumber

  !------------------------------------------------------
  ! Predefined paths and molecular graph info for CBMC growth
  type :: Pathing
    integer :: nTerminal         ! Number of terminal atoms
    integer :: nLinker           ! Number of linker atoms (internal chain atoms)
    integer :: nHub              ! Number of hub (branching) atoms
    integer(atomIntType), allocatable :: termAtoms(:)   ! List of terminal atoms
    integer(atomIntType), allocatable :: hubAtoms(:)    ! List of hub atoms
    integer :: nPaths                                    ! Number of available growth paths
    integer(atomIntType), allocatable :: path(:,:)       ! Atom indices along each path
    integer(atomIntType), allocatable :: pathMax(:)      ! Length of each path
  end type Pathing

  !------------------------------------------------------
  ! Dihedral angle biasing setup
  integer, parameter :: nDihBins = 2000                   ! Number of bins for dihedral histogram
  real(dp), parameter :: diBinSize = two_pi / real(nDihBins, dp)  ! Bin size in radians
  real(dp) :: startProb = 0.05_dp                         ! Default startup probability for biasing

  type :: DihedralAngle
    integer :: molType          ! Molecule type
    integer :: hubIndx          ! Hub atom index (branch point)
    integer :: dihedIndx        ! Dihedral angle index
    integer :: startBin         ! Initial bin index
    real(dp) :: accConst        ! Acceptance constant (normalization)
    real(dp) :: Hist(0:nDihBins)      ! Histogram of accepted angles
    real(dp) :: Prob(0:nDihBins)      ! Normalized biasing probabilities
    real(dp) :: Integral(0:nDihBins)  ! Cumulative integral for sampling
  end type DihedralAngle

  !------------------------------------------------------
  ! Trial counts for different Monte Carlo operations
  integer :: maxRosenTrial     ! Max Rosenbluth trials for standard CBMC
  integer :: maxSwapInTrial    ! Max CBMC trials for swap-in (insertion)
  integer :: maxSwapOutTrial   ! Max CBMC trials for swap-out (deletion)

  ! Number of trials per molecule type
  integer, allocatable :: nRosenTrials(:)
  integer, allocatable :: nSwapInTrials(:)
  integer, allocatable :: nSwapOutTrials(:)

  !------------------------------------------------------
  ! Regrowth strategy
  integer, allocatable :: regrowType(:)       ! Type of regrowth per molecule
  integer, allocatable :: regrowOrder(:,:)    ! Regrowth atom ordering
  real(dp), allocatable :: probTypeCBMC(:)    ! Regrowth move probabilities

  ! Topological data arrays
  type(BondNumber), allocatable :: topolArray(:)
  type(Pathing),    allocatable :: pathArray(:)

  ! Mapping atoms to paths
  integer, allocatable :: usedByPath(:,:)     ! (atom, pathID)
  integer, allocatable :: atomPathIndex(:,:)  ! Index of atom in path

  !------------------------------------------------------
  ! Dihedral angle control and storage
  integer :: totalDihed                        ! Total number of dihedrals in system
  type(DihedralAngle), allocatable :: dihedData(:)

  !------------------------------------------------------
  ! Predefined limits for growth branches and schedules
  integer, parameter :: maxBranches         = 6
  integer, parameter :: nRosenTwoBranch     = 5
  integer, parameter :: nRosenThreeBranch   = 1
  integer, parameter :: nRosenTorsion       = 100

  !------------------------------------------------------
  ! Growth schedule for molecule construction (CBMC regrowth logic)
  type :: GrowthSchedule
    integer :: GrowthSteps                     ! Total number of growth steps
    integer, allocatable :: GrowFrom(:)        ! Index of parent atom
    integer, allocatable :: GrowPrev(:)        ! Index of atom added before current
    integer, allocatable :: GrowNum(:)         ! Atom to grow at each step
    integer, allocatable :: TorNum(:)          ! Dihedral to use at each step
    integer, allocatable :: GrowList(:,:)      ! Candidate atoms to grow
    integer, allocatable :: TorList(:,:)       ! Candidate torsions
  end type GrowthSchedule

  ! Growth schedules for swap moves (1 per molecule type)
  type(GrowthSchedule), allocatable :: SwapGrowOrder(:)

end module CBMC_Variables
!==============================================================
module EnergyTables
  use VarPrecision
  implicit none

  !------------------------------------------------------
  ! Neighbor count of each molecule in the current cluster
  integer, allocatable :: neiCount(:)

  ! Total interaction energy of each molecule with the rest of the cluster
  ! Used to efficiently recompute energy changes during insert/delete moves
  real(dp), allocatable :: ETable(:)

  ! Total energy of each molecule type in the gas phase
  real(dp), allocatable :: E_gasIntra(:)

  ! Maximum interaction energy that each molecule has with its neighbors
  ! Used for AVBMC weighting and identifying strongly interacting partners
  real(dp), allocatable :: NeiETable(:)

  ! AVBMC-specific biasing weights, e.g., based on neighbor identity/distance
  ! Dimensions: (molecule type, cluster index) or AVBMC-specific interpretation
  real(dp), allocatable :: AlphaTargetIns(:,:)
  real(dp), allocatable :: AlphaTargetDel(:,:)

  !------------------------------------------------------
  ! Intermolecular energy components (entire cluster)

  ! Total interaction energy of the current configuration
  real(dp) :: E_Inter_T

  ! Change in total interaction energy due to a proposed MC move
  real(dp) :: E_Inter_Diff

  !------------------------------------------------------
  ! Additional energy components (only used if bonded terms are active)
  real(dp) :: E_NBond_T,   E_NBond_Diff     ! Non-bonded (e.g., LJ)
  real(dp) :: E_Stretch_T, E_Strch_Diff     ! Bond stretch
  real(dp) :: E_Bend_T,    E_Bend_Diff      ! Angle bending
  real(dp) :: E_Torsion_T, E_Tors_Diff      ! Dihedral torsion

end module EnergyTables
 !================================================================ 
! Module containing indexing functions for a Grand Canonical Monte Carlo nucleation
! simulation. Maps global molecule indices to molecule types and indices within their
! type, supporting systems with multiple molecule types (e.g., for long linear biomolecules
! with minDistCriteria). Used to select molecules and determine their type and index for
! Monte Carlo moves (e.g., SingleAtom_Translation, TemperatureMove).
module IndexingFunctions
  implicit none
  contains

  ! Maps a global molecule index (nMove, 1 to NTotal) to its molecule type (molType) and
  ! index within that type (molIndx). Uses NPart array (number of molecules per type).
  ! Returns molType and molIndx; stops with error if nMove is invalid.
  subroutine Get_MolIndex(nMove, NPart, molType, molIndx)
    implicit none
    integer, intent(in) :: nMove                ! Global molecule index (1 to NTotal)
    integer, intent(in) :: NPart(:)            ! Number of molecules per type (module-defined)
    integer, intent(out) :: molType            ! Molecule type index
    integer, intent(out) :: molIndx            ! Molecule index within type
    integer :: iType                           ! Loop index over molecule types
    integer :: nTypes                          ! Number of molecule types
    integer :: curLimit                        ! Cumulative molecule count

    ! Initialize variables
    nTypes = size(NPart)
    curLimit = 0
    molType = 1
    molIndx = 1

    ! Loop over molecule types to find the type and index
    do iType = 1, nTypes
      curLimit = curLimit + NPart(iType)
      if (nMove <= curLimit) then
        molIndx = nMove - (curLimit - NPart(iType))  ! Index within type
        molType = iType
        return
      endif
    enddo

    ! Error handling: nMove exceeds total molecules (NTotal)
    write(35, *) "Error in Get_MolIndex: Invalid nMove=", nMove, &
                 " molType=", molType, " molIndx=", molIndx
    stop "Get_MolIndex: nMove out of bounds"
  end subroutine Get_MolIndex

  ! Maps a global molecule index (nIndx) to its molecule type (nType) using NMAX array
  ! (maximum molecules per type). Returns nType; assumes nIndx is valid (1 to sum(NMAX)).
  ! If nIndx is invalid, returns nType=1 to avoid undefined behavior in pure context.
  pure subroutine Get_SubIndex(nIndx, nType, NMAX)
    implicit none
    integer, intent(in) :: nIndx               ! Global molecule index
    integer, intent(in) :: NMAX(:)            ! Maximum molecules per type (module-defined)
    integer, intent(out) :: nType              ! Molecule type index
    integer :: iType                          ! Loop index over molecule types
    integer :: nTypes                         ! Number of molecule types
    integer :: curLimit                       ! Cumulative molecule count

    ! Initialize variables
    nTypes = size(NMAX)
    curLimit = 0
    nType = 1  ! Default if nIndx is invalid

    ! Loop over molecule types to find the type
    do iType = 1, nTypes
      curLimit = curLimit + NMAX(iType)
      if (nIndx <= curLimit) then
        nType = iType
        return
      endif
    enddo
  end subroutine Get_SubIndex

  ! Returns the molecule type for a global molecule index (Indx) using NMAX array.
  ! Returns type index (1 to size(NMAX)) or 0 if Indx is invalid.
  pure integer function Get_MolType(Indx, NMAX) result(molType)
    implicit none
    integer, intent(in) :: Indx               ! Global molecule index
    integer, intent(in) :: NMAX(:)           ! Maximum molecules per type (module-defined)
    integer :: iType                         ! Loop index over molecule types
    integer :: nTypes                        ! Number of molecule types
    integer :: curLimit                      ! Cumulative molecule count

    ! Initialize variables
    nTypes = size(NMAX)
    curLimit = 0
    molType = 0  ! Default if Indx is invalid

    ! Loop over molecule types to find the type
    do iType = 1, nTypes
      curLimit = curLimit + NMAX(iType)
      if (Indx <= curLimit) then
        molType = iType
        return
      endif
    enddo
  end function Get_MolType

end module IndexingFunctions
!==============================================================
! Module: SimParameters
! Purpose:
!   Stores global parameters and flags that control the Monte Carlo
!   simulation of a single nucleating cluster in the grand canonical
!   ensemble, including:
!     - Cluster criteria and ensemble settings
!     - Particle number and type management
!     - Pressure, temperature, electrostatics
!     - WHAM/umbrella sampling
!     - Output formatting and units
!
! Dependencies:
!   - VarPrecision
!
! Author:
!   Aliasghar Sepehri
!===================================================================
module SimParameters
  use VarPrecision
  implicit none

  !------------------------------------------------------
  ! Input and simulation control flags
  logical, parameter :: echoInput      = .false.   ! Echo input file to output
  logical, parameter :: useScriptInput = .true.    ! Use script file for input

  logical :: distCriteria     = .false.            ! Use distance-based cluster criterion
  logical :: minDistCriteria  = .false.            ! Use minimum interatomic distance (Zerze 2024)

  ! For MPI: input contains multiple independent configurations
  logical :: multipleInput    = .false.

  ! Whether to calculate pressure using the Virial equation
  logical :: calcPressure     = .false.

  ! Whether to pair atom distances, distace squares, energy, etc.
!  logical :: useDistStore = .false.

  !------------------------------------------------------
  ! Monte Carlo cycle counters
  integer(kind=8) :: ncycle, ncycle2               ! MC cycle trackers

  !------------------------------------------------------
  ! Molecule type and count management
  integer(atomIntType) :: nMolTypes = 1            ! Number of molecule types

  integer, allocatable :: NMin(:), NMax(:)         ! Min/max per species (for GC constraints)
  integer, allocatable, target :: NPart(:)         ! Current number of particles per type
  integer, allocatable, target :: NPart_New(:)     ! Trial number after MC move

  logical, allocatable :: isActive(:)              ! Whether a molecule type is present in the cluster

  integer, target :: NTotal, NTotal_New            ! Total number of molecules (before/after trial)

  integer :: maxMol                                 ! Max number of molecules allowed in memory
  integer :: maxAtoms, vmdAtoms                     ! Atom bookkeeping for largest mol, VMD output

  !------------------------------------------------------
  ! AVBMC volume
  real(dp) :: avbmc_vol                             ! Sampling volume for AVBMC moves

  !------------------------------------------------------
  ! Umbrella sampling & WHAM data
  integer :: umbrellaLimit
  real(dp), allocatable :: NHist(:)                 ! Cluster size histogram
  real(dp), allocatable :: NBias(:)                 ! Biasing potential by cluster size
  real(dp), allocatable :: E_Avg(:)                 ! Energy average by cluster size

  !------------------------------------------------------
  ! Thermodynamic and electrostatic parameters
  real(dp), allocatable :: gas_dens(:)              ! Dilute-phase density for each molecule type

  real(dp), allocatable :: max_dist(:)              ! Max translation trial magnitude per type
  real(dp), allocatable :: max_dist_single(:)       ! Max per-atom displacement
  real(dp), allocatable :: max_rot(:)               ! Max rotation angle per molecule

  real(dp) :: beta                                   ! 1 / (k_B * T)
  real(dp), target :: temperature, tempNew           ! Current and trial temperatures

  real(dp) :: DielecConst                            ! Dielectric constant
  real(dp) :: DebyeLength                            ! Debye screening length
  real(dp) :: kapa                                   ! Inverse Debye length

  real(dp) :: global_r_min                           ! Shortest allowable atom-atom distance
  real(dp) :: softCutoff                             ! Soft interaction truncation 

  real(dp) :: Pressure, P_Diff                       ! Pressure (Virial-based) and change in pressure
  real(dp), allocatable :: P_Avg(:)                  ! Pressure as a function of cluster size

  !------------------------------------------------------
  ! Cluster formation criteria
  real(dp) :: Dist_Critr                             ! Cutoff distance for bonding
  real(dp) :: Dist_Critr_sq                          ! Precomputed squared version
  real(dp), allocatable :: Eng_Critr(:,:)            ! Energy-based clustering criteria
  real(dp) :: ECritMax                               ! Max energy threshold for cluster inclusion

  !------------------------------------------------------
  ! Output units and conversion factors
  character(len=10) :: outputEngUnits                ! Energy unit (e.g., kBT)
  character(len=10) :: outputLenUnits                ! Length unit (e.g., Å)
  character(len=10) :: engDefault, distDefault, angDefaults

  real(dp) :: outputEConv, outputLenConv             ! Conversion factors for formatted output
  real(dp) :: convEng, convDist, convAng             ! Internal conversion constants

  !------------------------------------------------------
  ! Output control: frequency of writing data
  integer :: outFreq_GCD    = 1000                   ! Energy output
  integer :: outFreq_Screen = 1000                   ! Screen print interval
  integer :: outFreq_Traj   = 1000                   ! Trajectory snapshot interval

  !------------------------------------------------------
  ! Molecule Names per type
  character(len=10), allocatable :: MolName(:)

  !------------------------------------------------------
  ! Previous Move is accpted
  logical :: prevMoveAccepted

end module SimParameters
!================================================================ 
      module UmbrellaFunctions

      contains

      pure integer function getBiasIndex(NPart,NMAX)
       implicit none
       integer,intent(in) :: NPart(:), NMAX(:)
       integer :: i, sizeN
       integer :: curIndx,maxIndx
        
       sizeN = size(NPART)
!       curIndx = 0
       curIndx = NPart(sizeN)
       maxIndx = NMAX(sizeN) + 1 
       do i=1,sizeN-1
         curIndx = curIndx + maxIndx*NPart(sizeN-i)
         maxIndx = maxIndx * (NMAX(sizeN-i) + 1)
       enddo
       getBiasIndex = curIndx + 1
       
      end function
!     ---------------------------------------------------------------
      pure integer function getNewBiasIndex(NPart,NMAX, increment)
       implicit none
       integer,intent(in) :: NPart(:), NMAX(:), increment(:)
       integer :: i, sizeN
       integer :: curIndx,maxIndx
        
       sizeN = size(NPART)
!       curIndx = 0
       curIndx = NPart(sizeN) + increment(sizeN)
       maxIndx = NMAX(sizeN) + 1 
       do i=1,sizeN-1
         curIndx = curIndx + maxIndx*(NPart(sizeN-i) + increment(sizeN-i))
         maxIndx = maxIndx * (NMAX(sizeN-i) + 1)
       enddo
       getNewBiasIndex = curIndx + 1
       
      end function      
    
!     ---------------------------------------------------------------      
      end module   
!======================================================
! Module: WHAM_Module
! Purpose:
!   Stores all WHAM-related data structures and parameters for computing
!   free energy profiles as a function of cluster size based on accumulated
!   biased/unbiased histograms from Monte Carlo sampling.
!
!   Includes:
!     - Self-consistent WHAM iteration control
!     - Histogram and bias storage
!     - Output toggles and convergence criteria
!
! Dependencies:
!   - VarPrecision (for real kind and integer precision)
!
! Author:
!   Aliasghar Sepehri
!===================================================================
module WHAM_Module
  use VarPrecision
  implicit none

  !------------------------------------------------------
  ! Output control
  logical, parameter :: WHAM_ExtensiveOutput = .false.  
  ! If true, output detailed WHAM data each iteration

  !------------------------------------------------------
  ! WHAM status and global settings
  logical :: useWHAM         = .false.  ! Enable/disable WHAM calculation
  integer :: intervalWHAM               ! MC cycles between WHAM updates
  integer :: maxSelfConsist             ! Max WHAM iterations per update
  real(dp) :: tolLimit                  ! Convergence tolerance (e.g., 1e-5)
  integer :: equilInterval              ! MC cycles before WHAM accumulation starts
  integer :: whamEstInterval            ! Interval between updates of bias estimates

  ! Reference bin for zero free energy 
  ! FreeEnergyEst(refBin) and NewBias(refBin) are set to zero; other values are relative to this
  integer :: refBin                     ! Index of reference cluster size
  real(dp), allocatable :: refSizeNumbers(:)  
  ! Cluster sizes used as reference for normalization (optional, continuous)

  !------------------------------------------------------
  ! WHAM iteration counters
  integer :: nWhamItter      ! Total WHAM iterations performed
  integer :: nCurWhamItter   ! WHAM iterations in current cycle

  !------------------------------------------------------
  ! Histogram and biasing arrays
  real(dp), allocatable :: WHAM_Numerator(:)        ! WHAM numerator per bin (accumulated)
  real(dp), allocatable :: WHAM_Denominator(:,:)    ! Denominator (per window, per bin)
  real(dp), allocatable :: HistStorage(:)           ! Accumulated total histogram

  real(dp), allocatable :: FreeEnergyEst(:)         ! F(x): WHAM-estimated free energy per bin
  real(dp), allocatable :: BiasStorage(:,:)         ! Bias applied during simulation (per cycle/bin)
  real(dp), allocatable :: NewBias(:)               ! Updated bias to be applied in next MC run

  real(dp), allocatable :: ProbArray(:)             ! Reconstructed probability P(x) from WHAM
  real(dp), allocatable :: LogProbArray(:)          ! Reconstructed logarithm of probability P(x)
  real(dp), allocatable :: TempHist(:)              ! Temporary histogram for current MC window

end module WHAM_Module
!================================================================   
