!****************************************************************************************
!   HPS_piecewise (HPS-KR, KH, FB, TSCL-M2, Urry: Piecewise HPS LJ + Debye-Hückel
!      - u_LJ(r)           = 4 * ε * [ (σ / r)**12 − (σ / r)**6 ]
!      - U_HPS(r) =
!            u_LJ(r) + ε * (1 − λ)       , if r <= 2**(1/6) * σ
!            λ * u_LJ(r)                 , if r >  2**(1/6) * σ
!      - U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!
!   HPS_cation_pi (HPS-cation-pi-i and -ii): HPS_piecewise + LJ for specific π interactions
!      - U_total(r)        = U_HPS_piecewise(r) + U_extra_LJ(r) for cation-π pairs
!      - U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!*********************************************************************************************************************
      module E_Interface_HPS_cation_pi_Diststore
      use CoordinateTypes
!=============================================================================
      contains
!=============================================================================      
! Computes total system energy for the HPS_cation_pi force field in a grand canonical
! ensemble nucleation simulation with useDistStore=.true. Calls CalcAllDistPairs to store
! distances, followed by intermolecular, intramolecular non-bonded, bond stretching, and
! gas-phase intramolecular energy calculations. Applies distance or energy criteria to
! determine move rejection. Sets zero bending and torsional energies. Updates total energy
! (E_T) and rejection flag (rejMove). Compatible with single-cluster simulation and MPI.
subroutine Detailed_ECalc_HPS_cation_pi_DStore(E_T, rejMove)
  use InterEnergy_HPS_cation_pi_DistStore, only: Detailed_ECalc_Inter
  use IntraEnergy_HPS_cation_pi, only: Detailed_ECalc_IntraNonBonded, Detailed_ECalc_GasIntra
  use BondStretchFunctions, only: Detailed_ECalc_BondStretch
  use EnergyTables, only: E_Bend_T, E_Torsion_T
  use EnergyCriteria, only: Detailed_EnergyCriteria
  use DistanceCriteria_PairStore, only: Detailed_DistanceCriteria
  use SimParameters, only: maxMol, distCriteria, minDistCriteria
  use PairStorage, only: CalcAllDistPairs
  implicit none
  real(dp), intent(inout) :: E_T                   ! Total system energy
  logical, intent(inout) :: rejMove                ! Move rejection flag
  real(dp) :: PairList(1:maxMol, 1:maxMol)         ! Intermolecular pair data

  ! Initialize total energy
  E_T = 0.0_dp

  ! Compute and store all atom pair distances
  call CalcAllDistPairs

  ! Calculate intermolecular energies
  call Detailed_ECalc_Inter(E_T, PairList)

  ! Apply rejection criteria
  if (distCriteria .or. minDistCriteria) then
    call Detailed_DistanceCriteria(rejMove)
  else
    call Detailed_EnergyCriteria(PairList, rejMove)
  endif

  ! Calculate intramolecular and gas-phase energies
  call Detailed_ECalc_IntraNonBonded(E_T)
  call Detailed_ECalc_BondStretch(E_T)
  call Detailed_ECalc_GasIntra

  ! Set zero bending and torsional energies
  E_Bend_T = 0.0_dp
  E_T = E_T + E_Bend_T
  E_Torsion_T = 0.0_dp
  E_T = E_T + E_Torsion_T
end subroutine Detailed_ECalc_HPS_cation_pi_DStore
!=============================================================================      
! Calculates the energy change for a single molecule displacement in the HPS_cation_pi force field when
! useDistStore = .true. in a Grand Canonical Monte Carlo nucleation simulation. Computes
! intermolecular (E_Inter) and intramolecular (E_Intra) energy differences, updates pair lists,
! and checks cluster criteria. Populates NeighborDetailsNew for minDistCriteria when provided.
! Rejects moves if atoms overlap, energy exceeds softCutOff, or cluster criteria fail. When
! distCriteria or minDistCriteria is .true., PairList is not used since all distances are stored.
subroutine Shift_ECalc_HPS_cation_pi_DStore(E_Inter, E_Intra, disp, PairList, dETable, useIntra, rejMove, &
                                            useInter, NeighborDetailsNew)
  use SimParameters, only: distCriteria, minDistCriteria, beta, softCutOff, NTotal, maxMol
  use BendingFunctions, only: Shift_ECalc_Bending
  use BondStretchFunctions, only: Shift_ECalc_BondStretch
  use CBMC_Variables, only: regrowType
  use Coords, only: Displacement, MolArray
  use DistanceCriteria_PairStore, only: Shift_DistanceCriteria
  use EnergyCriteria, only: Shift_EnergyCriteria
  use EnergyTables, only: E_NBond_Diff, E_Strch_Diff, E_Bend_Diff, E_Tors_Diff, E_Inter_Diff
  use InterEnergy_HPS_cation_pi_DistStore, only: Shift_ECalc_Inter
  use IntraEnergy_HPS_cation_pi, only: Shift_ECalc_IntraNonBonded
  use TorsionalFunctions, only: Shift_ECalc_Torsional
  use PairStorage, only: CalcNewDistPairs, newDist
  use VarPrecision
  implicit none

  ! Input/Output variables
  logical, intent(in) :: useIntra(1:4)
  logical, intent(in), optional :: useInter
  logical, intent(inout) :: rejMove
  type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(1:maxMol)
  type(Displacement), intent(in) :: disp(:)
  real(dp), intent(inout) :: PairList(:), dETable(:)
  real(dp), intent(out) :: E_Inter, E_Intra

  ! Local variables
  logical :: interSwitch
  integer :: nDisp, nIndx
  real(dp) :: E_NonBond, E_Stretch, E_Bend, E_Torsion, E_Improper

  ! Step 1: Initialize variables
  nDisp = size(disp)
  nIndx = MolArray(disp(1)%molType)%mol(disp(1)%molIndx)%indx
  interSwitch = .true.
  if (present(useInter)) interSwitch = useInter
  rejMove = .false.
  E_Inter = 0.0_dp
  E_Intra = 0.0_dp
  E_NonBond = 0.0_dp
  E_Stretch = 0.0_dp
  E_Bend = 0.0_dp
  E_Torsion = 0.0_dp
  E_Improper = 0.0_dp
  dETable = 0.0_dp
  E_Inter_Diff = 0.0_dp
  E_NBond_Diff = 0.0_dp
  E_Strch_Diff = 0.0_dp
  E_Bend_Diff = 0.0_dp
  E_Tors_Diff = 0.0_dp

  ! Step 2: Calculate intermolecular energy and check for overlaps
  if (interSwitch .and. NTotal > 1) then
    call CalcNewDistPairs(disp, rejMove)
    if (rejMove) return

    if (present(NeighborDetailsNew)) then
      call Shift_ECalc_Inter(E_Inter, disp, newDist, PairList, dETable, rejMove, NeighborDetailsNew)
    else
      call Shift_ECalc_Inter(E_Inter, disp, newDist, PairList, dETable, rejMove)
    endif
    if (rejMove) return
    E_Inter_Diff = E_Inter

    if (E_Inter * beta > softCutOff) then
      rejMove = .true.
      return
    endif

    ! Step 3: Check cluster criteria
    if (distCriteria .or. minDistCriteria) then
      if (present(NeighborDetailsNew)) then
        call Shift_DistanceCriteria(nIndx, rejMove, NeighborDetailsNew)
      else
        call Shift_DistanceCriteria(nIndx, rejMove)
      endif
    else
      call Shift_EnergyCriteria(PairList, nIndx, rejMove)
    endif
    if (rejMove) return
  endif

  ! Step 4: Calculate intramolecular energy if applicable
  if (regrowType(disp(1)%molType) /= 0) then
    if (useIntra(1)) then
      call Shift_ECalc_IntraNonBonded(E_NonBond, disp)
      E_NBond_Diff = E_NonBond
    endif
    if (useIntra(2)) then
      call Shift_ECalc_BondStretch(E_Stretch, disp)
      E_Strch_Diff = E_Stretch
    endif
    if (useIntra(3)) then
      call Shift_ECalc_Bending(E_Bend, disp)
      E_Bend_Diff = E_Bend
    endif
    if (useIntra(4)) then
      call Shift_ECalc_Torsional(E_Torsion, disp)
      E_Tors_Diff = E_Torsion
    endif
    ! Note: Improper angle calculation is excluded
    E_Intra = E_NonBond + E_Stretch + E_Bend + E_Torsion + E_Improper
  endif

  ! Step 5: Verify E_Intra for rigid moves
  if (.not. any(useIntra) .and. abs(E_Intra) > 1.0E-6_dp) then
    write(35, *) "Error: Non-zero E_Intra=", E_Intra, " in Shift_ECalc_HPS_cation_pi_DStore with useIntra=.false."
    stop "Shift_ECalc_HPS_cation_pi_DStore: Invalid intramolecular energy"
  endif
end subroutine Shift_ECalc_HPS_cation_pi_DStore
!=============================================================================      
! Calculates the intermolecular and intramolecular energies for a newly inserted molecule in a Grand
! Canonical Monte Carlo nucleation simulation using the HPS_cation_pi force field with useDistStore=.true.
! Computes E_Inter (via NewMol_ECalc_Inter) and E_Intra (via NewMol_ECalc_IntraNonBonded,
! NewMol_ECalc_BondStretch, and optionally bending/torsional terms). Updates PairList and dETable
! for cluster criteria and NeighborDetailsNew (if present) for neighbor tracking. Checks for overlaps
! using CalcSwapInDistPairs. Sets rejMove = .true. if overlaps or invalid configurations are detected.
! Supports optional useInter flag to toggle intermolecular calculations. Conditionally passes
! NeighborDetailsNew to NewMol_ECalc_Inter based on its presence. Used in swap-in moves (e.g., AVBMC).
subroutine SwapIn_ECalc_HPS_cation_pi_DStore(E_Inter, E_Intra, PairList, dETable, rejMove, &
                                             useInter, NeighborDetailsNew)
  use InterEnergy_HPS_cation_pi_DistStore, only: NewMol_ECalc_Inter
  use IntraEnergy_HPS_cation_pi, only: NewMol_ECalc_IntraNonBonded
  use BondStretchFunctions, only: NewMol_ECalc_BondStretch
  use BendingFunctions, only: NewMol_ECalc_Bending
  use TorsionalFunctions, only: NewMol_ECalc_Torsional
  use EnergyTables, only: E_NBond_Diff, E_Strch_Diff, E_Bend_Diff, E_Tors_Diff, E_Inter_Diff
  use Coords, only: newMol
  use CBMC_Variables, only: regrowType
  use PairStorage, only: CalcSwapInDistPairs
  use VarPrecision
  implicit none

  ! Input/Output variables
  real(dp), intent(out) :: E_Inter             ! Intermolecular energy
  real(dp), intent(out) :: E_Intra             ! Intramolecular energy
  real(dp), intent(inout) :: PairList(:)       ! Pairwise interaction list (distances or energies)
  real(dp), intent(inout) :: dETable(:)        ! Energy difference table for cluster criteria
  logical, intent(out) :: rejMove               ! Flag for invalid configurations (e.g., overlaps)
  logical, intent(in), optional :: useInter     ! Flag to toggle intermolecular calculations
  type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)  ! Neighbor relationships

  ! Local variables
  logical :: interSwitch                       ! Local flag for intermolecular calculations
  real(dp) :: E_NonBond, E_Stretch, E_Bend     ! Intramolecular energy components
  real(dp) :: E_Torsion, E_Improper            ! Additional intramolecular energy components

  ! Step 1: Initialize variables
  ! Overview: Sets up energy accumulators, difference tables, and determines intermolecular switch.
  rejMove = .false.                            ! Initialize rejection flag
  E_Inter = 0.0_dp                             ! Initialize intermolecular energy
  E_Intra = 0.0_dp                             ! Intramolecular energy
  E_NonBond = 0.0_dp                           ! Initialize non-bonded intramolecular energy
  E_Stretch = 0.0_dp                           ! Initialize bond stretching energy
  E_Bend = 0.0_dp                              ! Initialize bending energy
  E_Torsion = 0.0_dp                           ! Initialize torsional energy
  E_Improper = 0.0_dp                          ! Initialize improper torsion energy
  E_Inter_Diff = 0.0_dp                        ! Initialize intermolecular energy difference
  E_NBond_Diff = 0.0_dp                        ! Initialize non-bonded energy difference
  E_Strch_Diff = 0.0_dp                        ! Initialize bond stretching energy difference
  E_Bend_Diff = 0.0_dp                         ! Initialize bending energy difference
  E_Tors_Diff = 0.0_dp                         ! Initialize torsional energy difference
  interSwitch = .true.                         ! Default: calculate intermolecular interactions
  if (present(useInter)) interSwitch = useInter ! Override with optional flag

  ! Step 2: Calculate intermolecular interactions
  ! Overview: Computes distances and intermolecular energies if interSwitch is true.
  if (interSwitch) then
    call CalcSwapInDistPairs(rejMove)          ! Check overlaps and store distances
    if (rejMove) return                        ! Return if overlap detected

    ! Conditionally call NewMol_ECalc_Inter based on NeighborDetailsNew presence
    if (present(NeighborDetailsNew)) then
      call NewMol_ECalc_Inter(E_Inter, PairList, dETable, rejMove, NeighborDetailsNew)
    else
      call NewMol_ECalc_Inter(E_Inter, PairList, dETable, rejMove)
    end if
    if (rejMove) return                        ! Return if invalid configuration detected
    E_Inter_Diff = E_Inter                     ! Store intermolecular energy difference
  end if

  ! Step 3: Calculate intramolecular interactions
  ! Overview: Computes non-bonded, bond stretching, and optionally bending/torsional energies.
  if (regrowType(newMol%molType) /= 0) then     ! Check if molecule requires intramolecular calculations
    call NewMol_ECalc_IntraNonBonded(E_NonBond) ! Non-bonded intramolecular energy
    call NewMol_ECalc_BondStretch(E_Stretch)   ! Bond stretching energy
    ! Note: Bending and torsional calculations are commented out in the original code
    ! call NewMol_ECalc_Bending(E_Bend)        ! Bond angle bending energy (disabled)
    ! call NewMol_ECalc_Torsional(E_Torsion)     ! Torsional energy (disabled)
    ! call NewMol_ECalc_Improper(E_Improper)   ! Improper torsion energy (disabled)
    E_Intra = E_NonBond + E_Stretch + E_Bend + E_Torsion + E_Improper
    E_NBond_Diff = E_NonBond                   ! Store non-bonded energy difference
    E_Strch_Diff = E_Stretch                   ! Store bond stretching energy difference
    E_Bend_Diff = E_Bend                       ! Store bending energy difference
    E_Tors_Diff = E_Torsion                   ! Store torsional energy difference
  end if
end subroutine SwapIn_ECalc_HPS_cation_pi_DStore
!=============================================================================
! Calculates the intermolecular and intramolecular energy of a molecule being removed in a Grand Canonical
! Monte Carlo swap-out move using the HPS_cation_pi force field with useDistStore=.true. Computes E_Inter
! (piecewise Hydropathy Scale Lennard-Jones with switching at (sigma/r)^6 = 0.5, cation-pi Lennard-Jones, and
! Debye-Hückel electrostatics using stored distances) and E_Intra (non-bonded, bond stretch, bending, torsional,
! and improper terms) for the specified molecule, updating dETable for cluster criteria. Optionally skips
! intermolecular interactions if useInter is present and .false. Used in AVBMC or similar moves.
subroutine SwapOut_ECalc_HPS_cation_pi_DStore(E_Inter, E_Intra, nType, nMol, dETable, useInter)
  use VarPrecision
  use CBMC_Variables, only: regrowType
  use EnergyTables, only: E_Inter_Diff, E_NBond_Diff, E_Strch_Diff, E_Bend_Diff, E_Tors_Diff
  use InterEnergy_HPS_cation_pi_DistStore, only: Mol_ECalc_Inter
  use IntraEnergy_HPS_cation_pi, only: Mol_ECalc_IntraNonBonded
  use BondStretchFunctions, only: Mol_ECalc_BondStretch
  implicit none

  ! Interface variables
  real(dp), intent(out) :: E_Inter, E_Intra   ! Inter- and intramolecular energies
  integer, intent(in) :: nType, nMol          ! Molecule type and index
  real(dp), intent(inout) :: dETable(:)       ! Energy difference table
  logical, intent(in), optional :: useInter   ! Flag to include intermolecular interactions

  ! Local variables
  logical :: interSwitch                      ! Determines if intermolecular interactions are computed
  real(dp) :: E_NonBond, E_Stretch, E_Bend    ! Intramolecular energy components
  real(dp) :: E_Torsion, E_Improper           ! Intramolecular energy components

  ! Step 1: Initialize variables
  ! Overview: Sets up energy accumulators and difference tables.
  E_Inter = 0.0_dp
  E_Intra = 0.0_dp
  E_NonBond = 0.0_dp
  E_Stretch = 0.0_dp
  E_Bend = 0.0_dp
  E_Torsion = 0.0_dp
  E_Improper = 0.0_dp
  E_Inter_Diff = 0.0_dp
  E_NBond_Diff = 0.0_dp
  E_Strch_Diff = 0.0_dp
  E_Bend_Diff = 0.0_dp
  E_Tors_Diff = 0.0_dp
  dETable = 0.0_dp
  interSwitch = .true.
  if (present(useInter)) interSwitch = useInter

  ! Step 2: Calculate intermolecular energy
  ! Overview: Computes piecewise HPS Lennard-Jones, cation-pi Lennard-Jones, and Debye-Hückel electrostatic energies using stored distances, negated for removal.
  if (interSwitch) then
    call Mol_ECalc_Inter(nType, nMol, dETable, E_Inter)
    E_Inter = -E_Inter
    E_Inter_Diff = E_Inter
  end if

  ! Step 3: Calculate intramolecular energy
  ! Overview: Computes non-bonded, bond stretch, bending, torsional, and improper energies for flexible molecules.
  if (regrowType(nType) /= 0) then
    call Mol_ECalc_IntraNonBonded(nType, nMol, E_NonBond)
    call Mol_ECalc_BondStretch(nType, nMol, E_Stretch)
    ! Note: Bending, torsional, and improper calculations are currently disabled.
    E_Intra = E_NonBond + E_Stretch + E_Bend + E_Torsion + E_Improper
    E_Intra = -E_Intra
    E_NBond_Diff = -E_NonBond
    E_Strch_Diff = -E_Stretch
    E_Bend_Diff = -E_Bend
    E_Tors_Diff = -E_Torsion
  end if
end subroutine SwapOut_ECalc_HPS_cation_pi_DStore
!=============================================================================
 
      
      end module
      
      
