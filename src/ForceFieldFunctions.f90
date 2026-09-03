!==============================================================
!  This file contains a set of small functions that are either used in creating forcefield
!  tables such as mixing rules or contain functional forms that are used in the energy calculations
!  themselves. 
!==============================================================
!==============================================================
! Module: ForceFieldVariableType
! Purpose: Defines custom data types used to store atomistic and interaction
!          parameters for bonded and nonbonded potentials in force field models
!==============================================================
module ForceFieldVariableType
  use VarPrecision

  ! AtomDef stores atomic information: symbol, name, Lennard-Jones parameters, charge, and mass
  type AtomDef
    character(len=5) :: Symb           ! Chemical symbol (e.g., "C", "Na")
    character(len=20) :: atmName       ! Descriptive name
    real(dp) :: sig, ep                ! Lennard-Jones sigma and epsilon
    real(dp) :: q                      ! Atomic charge
    real(dp) :: mass                   ! Atomic mass (amu)
  end type
        
  ! BondDef stores harmonic bond parameters
  ! r_sigma = 1 / sqrt(beta * k_eq)
  ! rmax_sq = (r_eq + 5 * r_sigma)^2
  type BondDef
    character(len=20) :: bondName      ! Descriptive name of bond
    real(dp) :: k_eq                   ! Force constant
    real(dp) :: r_eq                   ! Equilibrium bond length
    real(dp) :: r_sigma                ! Bond width parameter for tabulation
    real(dp) :: rmax_sq                ! Squared maximum bond distance (for validation)
  end type

  ! Parameters for histogram-based bend angle sampling in CBMC
  integer, parameter :: nBendHistBins = 200
  real(dp), parameter :: bendBinWidth = 4e0 * atan(1e0) / real(nBendHistBins, dp)
  real(dp), parameter :: startProb = 0.05e0

  ! BendAngleDef stores harmonic bend angle parameters and CBMC sampling histogram
  type BendAngleDef
    character(len=20) :: angleName     ! Descriptive name of angle
    real(dp) :: k_eq                   ! Force constant
    real(dp) :: ang_eq                 ! Equilibrium angle (radians)
    integer :: startBin                ! Starting bin for histogram sampling
    real(dp) :: accptConstant          ! Precomputed acceptance weight
    real(dp) :: Prob(1:nBendHistBins)  ! Probability histogram
  end type

  ! Torsion angle definition (e.g., OPLS style with nPara parameters)
  type TorsionAngleDef
    character(len=20) :: torsName      ! Torsion name
    integer(kind=atomIntType) :: nPara ! Number of torsion parameters
    real(dp), allocatable :: a(:)      ! Array of parameters (e.g., V1, V2, ...)
  end type

  ! Defines indices for nonbonded intramolecular atom pairs
  type NonBondedIndex
    integer(kind=atomIntType) :: nonMembr(1:2)  ! Atom indices of nonbonded pair
  end type

  ! Defines indices and type for bonded atom pairs
  type BondIndex
    integer(kind=atomIntType) :: bondType       ! Bond type ID
    integer(kind=atomIntType) :: bondMembr(1:2)  ! Indices of atoms in the bond
  end type

  ! Bend angle definition based on atom triplet
  type BendAngleIndex
    integer(kind=atomIntType) :: bendType        ! Bend angle type ID
    integer(kind=atomIntType) :: bendMembr(1:3)   ! Indices of atoms forming angle
  end type

  ! Torsion angle defined by 4 atoms
  type TorsionAngleIndex
    integer(kind=atomIntType) :: TorsType        ! Torsion type ID
    integer(kind=atomIntType) :: torsMembr(1:4)   ! Atom indices for torsion
  end type

  ! Improper torsion angle, typically for planarity constraints
  type ImproperAngleIndex
    integer(kind=atomIntType) :: ImpropType       ! Improper type ID
    integer(kind=atomIntType) :: impropMembr(1:4)  ! Atom indices involved
  end type

end module ForceFieldVariableType
!==============================================================
! Module: ForceField
! Purpose: Stores force field definitions and atomistic topology for each molecule type
!==============================================================
module ForceField
  use ForceFieldVariableType
  use VarPrecision

  character(len=20) :: ForceFieldName             ! Name of the current force field

  ! Number of interaction types
  integer :: nAtomTypes = 0
  integer :: nBondTypes = 0
  integer :: nAngleTypes = 0
  integer :: nTorsionalTypes = 0
  integer :: nImproperTypes = 0

  ! Force field parameters
  type(AtomDef), allocatable :: atomData(:)       ! Atom parameters
  type(BondDef), allocatable :: bondData(:)       ! Bond parameters
  type(BendAngleDef), allocatable :: bendData(:)  ! Angle parameters
  type(TorsionAngleDef), allocatable :: torsData(:)    ! Torsion parameters
  type(TorsionAngleDef), allocatable :: impropData(:)  ! Improper torsion parameters

  ! Tabulated minimum interaction distances
  real(dp), allocatable :: r_min(:), r_min_sq(:)           ! Minimum distances and squared values
  real(dp), allocatable :: r_min_tab(:,:)                  ! Pairwise min distances by type

  ! Topology data per molecule type
  integer, allocatable :: nAtoms(:)                        ! Number of atoms per molecule
  integer, allocatable :: nIntraNonBond(:)                 ! Number of nonbonded pairs
  integer, allocatable :: nBonds(:)                        ! Number of bonds
  integer, allocatable :: nAngles(:)                       ! Number of bend angles
  integer, allocatable :: nTorsional(:)                    ! Number of torsions
  integer, allocatable :: nImproper(:)                     ! Number of impropers

  real(dp), allocatable :: totalMass(:)                    ! Total mass of molecule type

  ! Connectivity arrays
  integer(kind=atomIntType), allocatable :: atomArray(:,:)         ! Atom types per molecule
  type(NonBondedIndex), allocatable :: nonBondArray(:,:)          ! Intra nonbonded pairs
  type(BondIndex), allocatable :: bondArray(:,:)                  ! Bond indices and types
  type(BendAngleIndex), allocatable :: bendArray(:,:)            ! Angle indices and types
  type(TorsionAngleIndex), allocatable :: torsArray(:,:)         ! Torsion indices and types
  type(ImproperAngleIndex), allocatable :: impropArray(:,:)      ! Improper indices and types

end module ForceField
!==============================================================
! Module: ForceFieldFunctions
! Purpose: Defines reusable potential energy functions and mixing rules
!==============================================================
module ForceFieldFunctions

  !=========================== INTERFACES ============================

  interface
    ! Geometric/mean mixing rule interface
    real(dp) function MixRule(par1, par2)
      use VarPrecision
      real(dp), intent(in) :: par1, par2
    end function
  end interface

  interface
    ! Torsional energy function interface
    real(dp) function TorsEng(angle, coeffArray)
      use VarPrecision
      real(dp), intent(in) :: angle
      real(dp), intent(in) :: coeffArray(:)
    end function
  end interface

  !=========================== IMPLEMENTATIONS ============================

contains

  !------------------------------------------------------------------
  ! Geometric Mean Mixing Rule
  ! Returns sqrt(par1 * par2) unless one of the inputs is zero
  real(dp) function GeoMean_MixingFunc(par1, par2)
    use VarPrecision
    real(dp), intent(in) :: par1, par2

    if ((par1 .eq. 0) .or. (par2 .eq. 0)) then
      GeoMean_MixingFunc = 0e0
    else
      GeoMean_MixingFunc = sqrt(par1 * par2)
    endif
  end function

  !------------------------------------------------------------------
  ! Arithmetic Mean Mixing Rule
  ! Returns (par1 + par2)/2 unless one of the inputs is zero
  real(dp) function Mean_MixingFunc(par1, par2)
    use VarPrecision
    real(dp), intent(in) :: par1, par2

    if ((par1 .eq. 0e0) .or. (par2 .eq. 0e0)) then
      Mean_MixingFunc = 0e0
    else
      Mean_MixingFunc = 0.5e0 * (par1 + par2)
    endif
  end function

  !------------------------------------------------------------------
  ! Mean-Max Mixing Rule (custom scheme)
  ! Returns mean if both inputs > 0, else max of abs values
  real(dp) function MeanMax_MixingFunc(par1, par2)
    use VarPrecision
    real(dp), intent(in) :: par1, par2

    if ((par1 .gt. 0e0) .and. (par2 .gt. 0e0)) then
      MeanMax_MixingFunc = abs(0.5e0 * (par1 + par2))
    else
      MeanMax_MixingFunc = max(abs(par1), abs(par2))
    endif
  end function

  !------------------------------------------------------------------
  ! Minimum Mixing Rule
  real(dp) function Min_MixingFunc(par1, par2)
    use VarPrecision
    real(dp), intent(in) :: par1, par2

    Min_MixingFunc = min(par1, par2)
  end function

  !------------------------------------------------------------------
  ! Maximum Mixing Rule
  real(dp) function Max_MixingFunc(par1, par2)
    use VarPrecision
    real(dp), intent(in) :: par1, par2

    Max_MixingFunc = max(par1, par2)
  end function

  !------------------------------------------------------------------
  ! Harmonic Bond Potential
  ! Returns 0.5 * k_eq * (angle - ang_eq)^2
  pure real(dp) function Harmonic(angle, k_eq, ang_eq)
    use VarPrecision
    real(dp), intent(in) :: angle, k_eq, ang_eq

    Harmonic = 0.5e0 * k_eq * (angle - ang_eq)**2
  end function

  !------------------------------------------------------------------
  ! Harmonic Torsional Potential
  ! CoeffArray(1) = k_eq, CoeffArray(2) = angle_eq
  pure real(dp) function TorsHarmonic(angle, coeffArray)
    use VarPrecision
    real(dp), intent(in) :: angle
    real(dp), intent(in) :: coeffArray(:)

    TorsHarmonic = 0.5e0 * coeffArray(1) * (angle - coeffArray(2))**2
  end function

  !------------------------------------------------------------------
  ! Torsional Potential: Cosine raised to powers
  ! Returns energy using Σ coeff_i * cos(angle)^(i-1)
  pure real(dp) function Torsion_Cos_N(angle, coeffArray)
    use VarPrecision
    integer :: nCoeff, i
    real(dp), intent(in) :: angle
    real(dp), intent(in) :: coeffArray(:)
    real(dp) :: eng

    nCoeff = size(coeffArray)
    eng = coeffArray(1)
    do i = 2, nCoeff
      if (coeffArray(i) .ne. 0e0) then
        eng = eng + coeffArray(i) * cos(angle)**(i - 1)
      endif
    enddo
    Torsion_Cos_N = eng
  end function

  !------------------------------------------------------------------
  ! Torsional Potential: Cosine of integer multiples of angle
  ! Returns energy using Σ coeff_i * cos((i-1) * angle)
  pure real(dp) function Torsion_CosNx(angle, coeffArray)
    use VarPrecision
    integer :: nCoeff, i
    real(dp), intent(in) :: angle
    real(dp), intent(in) :: coeffArray(:)
    real(dp) :: eng

    nCoeff = size(coeffArray)
    eng = coeffArray(1)
    do i = 2, nCoeff
      if (coeffArray(i) .ne. 0e0) then
        eng = eng + coeffArray(i) * cos(dble(i - 1) * angle)
      endif
    enddo
    Torsion_CosNx = eng
  end function

  !------------------------------------------------------------------
  ! Trappe Torsional Potential
  ! Uses specific form with cos(x), cos(2x), cos(3x)
  pure real(dp) function Trappe_CosNx(angle, coeffArray)
    use VarPrecision
    integer :: nCoeff
    real(dp), intent(in) :: angle
    real(dp), intent(in) :: coeffArray(:)
    real(dp) :: eng

    nCoeff = size(coeffArray)
    eng = coeffArray(1)
    eng = eng + coeffArray(2) * (1e0 + cos(angle))
    eng = eng + coeffArray(3) * (1e0 - cos(2e0 * angle))
    eng = eng + coeffArray(4) * (1e0 + cos(3e0 * angle))
    Trappe_CosNx = eng
  end function

end module ForceFieldFunctions
!==============================================================
! Module: ForceFieldPara_LJ_Q
! Purpose: Stores mixing rule functions and parameter tables for Lennard-Jones + Electrostatic force field
!==============================================================
module ForceFieldPara_LJ_Q

  use ForceFieldFunctions          ! Provides mixing rule interfaces
  use ForceFieldVariableType       ! Type definitions for force field data
  use VarPrecision                 ! Precision definition (e.g., dp = double precision)

  !====================== PROCEDURE POINTERS =======================
  procedure(MixRule), pointer :: ep_func  => null()  ! Mixing rule for epsilon
  procedure(MixRule), pointer :: sig_func => null()  ! Mixing rule for sigma

  !====================== FORCE FIELD TABLES =======================
  real(dp), allocatable :: ep_tab(:,:)     ! Epsilon parameters between atom types
  real(dp), allocatable :: sig_tab(:,:)    ! Sigma parameters between atom types
  real(dp), allocatable :: q_tab(:,:)      ! Charge table between atom types

end module ForceFieldPara_LJ_Q
!==============================================================
! Module: ForceFieldPara_Mpipi
! Purpose: Stores interaction parameter tables for Mpipi coarse-grained force field
!==============================================================
module ForceFieldPara_Mpipi

  use ForceFieldFunctions         ! Provides mixing rules and energy functions
  use ForceFieldVariableType      ! Custom types used for force field data
  use VarPrecision                ! Precision specification (e.g., dp = double precision)

  !======================= GLOBAL CONSTANTS ========================
  real(dp) :: rcutElec            ! Real-space electrostatic cutoff
  real(dp) :: rcutElecsq          ! Square of electrostatic cutoff (for efficiency)

  !======================= FORCE FIELD TABLES ======================
  real(dp), allocatable :: sig_tab(:,:)       ! Sigma values between atom types
  real(dp), allocatable :: sigsq_tab(:,:)     ! Square of sigma
  real(dp), allocatable :: epAlpha_tab(:,:)   ! Scaled epsilon parameters
  real(dp), allocatable :: q_tab(:,:)         ! Partial charge matrix

  integer, allocatable :: Mu_tab(:,:)         ! Charge scaling matrix (multipole order, etc.)
  logical, allocatable :: q_Nonzero(:,:)      ! Flags if charges are nonzero

end module ForceFieldPara_Mpipi
!==============================================================
! Module: ForceFieldPara_HPS_single
! Purpose: Stores interaction parameters for the single-component HPS force field model
!==============================================================
module ForceFieldPara_HPS_single

  use ForceFieldFunctions         ! Provides generic mixing rule functions
  use ForceFieldVariableType      ! Custom types for force field variables
  use VarPrecision                ! Precision definition (e.g., dp for double precision)

  !======================= GLOBAL CONSTANTS ========================
  real(dp) :: rcutElec            ! Electrostatic cutoff distance
  real(dp) :: rcutElecsq          ! Square of electrostatic cutoff (for efficiency)

  !======================= FORCE FIELD PARAMETERS ==================
  real(dp), allocatable :: sig_tab(:,:)        ! Sigma values between atom types
  real(dp), allocatable :: sigsq_tab(:,:)      ! Square of sigma
  real(dp), allocatable :: eps_tab(:,:)        ! Epsilon values (interaction strength)
  real(dp), allocatable :: lambda_tab(:,:)     ! Lambda scaling factor
  real(dp), allocatable :: q_tab(:,:)          ! Partial charges

  real(dp), allocatable :: cutoffNB_tab(:,:)   ! Non-bonded cutoff distances
  real(dp), allocatable :: cutoffNBsq_tab(:,:) ! Squares of non-bonded cutoffs

  logical, allocatable :: q_Nonzero(:,:)       ! Flags for nonzero charges

end module ForceFieldPara_HPS_single
!==============================================================
module ForceFieldPara_HPS_piecewise

  use ForceFieldFunctions         ! Provides mixing rule and energy functions
  use ForceFieldVariableType      ! Type definitions for force field parameters
  use VarPrecision                ! Precision definition (e.g., dp = double precision)

  !======================= GLOBAL CONSTANTS ========================
  real(dp) :: rcutElec            ! Electrostatic interaction cutoff
  real(dp) :: rcutElecsq          ! Square of electrostatic cutoff (for performance)

  !======================= FORCE FIELD PARAMETERS ==================
  real(dp), allocatable :: sig_tab(:,:)         ! Sigma values between atom types
  real(dp), allocatable :: sigsq_tab(:,:)       ! Square of sigma
  real(dp), allocatable :: eps_tab(:,:)         ! Epsilon interaction strengths
  real(dp), allocatable :: lambda_tab(:,:)      ! Lambda interaction scaling factors
  real(dp), allocatable :: q_tab(:,:)           ! Partial charge matrix

  real(dp), allocatable :: cutoffNB_tab(:,:)    ! Non-bonded interaction cutoffs
  real(dp), allocatable :: cutoffNBsq_tab(:,:)  ! Square of non-bonded cutoffs

  logical, allocatable :: q_Nonzero(:,:)        ! Flags for nonzero charges

end module ForceFieldPara_HPS_piecewise
!==============================================================
module ForceFieldPara_HPS_cation_pi

  use ForceFieldFunctions         ! Provides mixing rules and energy functions
  use ForceFieldVariableType      ! Custom types for force field definitions
  use VarPrecision                ! Precision control (e.g., dp = double precision)

  !======================= GLOBAL CONSTANTS ========================
  real(dp) :: rcutElec            ! Electrostatic cutoff distance
  real(dp) :: rcutElecsq          ! Square of electrostatic cutoff

  !======================= FORCE FIELD PARAMETERS ==================
  real(dp), allocatable :: sig_tab(:,:)          ! Sigma values between atom types
  real(dp), allocatable :: sigsq_tab(:,:)        ! Square of sigma
  real(dp), allocatable :: eps_tab(:,:)          ! Epsilon interaction parameters
  real(dp), allocatable :: lambda_tab(:,:)       ! Lambda scaling factors
  real(dp), allocatable :: q_tab(:,:)            ! Charge matrix

  real(dp), allocatable :: cutoffNB_tab(:,:)     ! Non-bonded cutoff distances
  real(dp), allocatable :: cutoffNBsq_tab(:,:)   ! Squares of non-bonded cutoffs
  real(dp), allocatable :: epsLJ_tab(:,:)        ! Additional LJ epsilon table for cation–π interactions

  logical, allocatable :: q_Nonzero(:,:)         ! Flags indicating nonzero charge entries

end module ForceFieldPara_HPS_cation_pi
!==============================================================
module ForceFieldPara_Pedone

  use VarPrecision                ! Precision definition (e.g., dp = double precision)

  !======================= ATOM TYPE DEFINITION =====================
  type AtomDefPedone
    character(len=20) :: atmName   ! Atom name or label
    character(len=5)  :: Symb      ! Chemical symbol
    real(dp) :: repul              ! Repulsion parameter
    real(dp) :: rEq                ! Equilibrium bond length or vdW radius
    real(dp) :: q                  ! Partial charge
    real(dp) :: alpha              ! Polarizability
    real(dp) :: delta              ! Dispersion correction or parameter
    real(dp) :: mass               ! Atomic mass
  end type

  !======================= GLOBAL FLAGS =============================
  logical :: implcSolvent = .false.   ! Flag for using implicit solvent model

  !======================= FORCE FIELD TABLES =======================
  type(AtomDefPedone), allocatable :: pedoneData(:)     ! Atom-wise parameter table

  real(dp), allocatable :: repul_tab(:,:)   ! Repulsion parameters between atom types
  real(dp), allocatable :: rEq_tab(:,:)     ! Equilibrium distances between atom types
  real(dp), allocatable :: q_tab(:,:)       ! Partial charges
  real(dp), allocatable :: alpha_Tab(:,:)   ! Polarizabilities
  real(dp), allocatable :: D_Tab(:,:)       ! Dispersion/delta table
  real(dp), allocatable :: bornRad(:)       ! Born radii for implicit solvent calculations

end module ForceFieldPara_Pedone
!==============================================================
module ForceFieldPara_Tersoff

  use VarPrecision                ! Precision definition (e.g., dp = double precision)

  !======================= ATOM TYPE DEFINITION =====================
  type AtomDefTersoff
    character(len=20) :: atmName   ! Atom name or label
    character(len=5)  :: Symb      ! Chemical symbol
    real(dp) :: A                  ! Tersoff parameter A (repulsive term)
    real(dp) :: B                  ! Tersoff parameter B (attractive term)
    real(dp) :: c, d               ! Cutoff function parameters
    real(dp) :: n                  ! Exponential parameter for bond order
    real(dp) :: lam1, lam2         ! Lambda parameters for repulsion and attraction
    real(dp) :: h                  ! Bond angle term coefficient
    real(dp) :: R                  ! Cutoff distance
    real(dp) :: D2                 ! Width of cutoff region
    real(dp) :: beta               ! Bond order parameter
    real(dp) :: mass               ! Atomic mass
  end type

  !======================= FORCE FIELD DATA =========================
  type(AtomDefTersoff), allocatable :: tersoffData(:)   ! Atom-wise Tersoff parameter table

end module ForceFieldPara_Tersoff
!==============================================================
