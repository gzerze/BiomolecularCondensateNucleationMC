!****************************************************************************************
!   HPS_piecewise (HPS-KR, KH, FB, TSCL-M2, Urry: Piecewise HPS LJ + Debye-Hückel
!      - u_LJ(r)           = 4 * ε * [ (σ / r)**12 − (σ / r)**6 ]
!      - U_HPS(r) =
!            u_LJ(r) + ε * (1 − λ)       , if r <= 2**(1/6) * σ
!            λ * u_LJ(r)                 , if r >  2**(1/6) * σ
!      - U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!****************************************************************************************
      module E_Interface_HPS_piecewise
      use CoordinateTypes
!=============================================================================
      contains
!=============================================================================      
! Calculates the total energy of the initial configuration for the HPS_piecewise force field
! (piecewise hydrophobic-polar model for biomolecules) in a grand canonical ensemble
! nucleation simulation. Checks if the configuration forms a single cluster using one
! of three criteria: energy-based (interaction energy below threshold), distance-based
! (first atom distance below threshold), or minimum distance-based (any atom pair
! distance below threshold). Sets rejectMove to .true. if the cluster criterion fails.
! Uses stored pair distances if useDistStore=.true. for faster calculations.
subroutine Detailed_ECalc_HPS_piecewise(E_T, rejectMove)
  use InterEnergy_HPS_piecewise, only: Detailed_ECalc_Inter
  use IntraEnergy_HPS_piecewise, only: Detailed_ECalc_IntraNonBonded, Detailed_ECalc_GasIntra
  use BondStretchFunctions, only: Detailed_ECalc_BondStretch
  use EnergyTables, only: E_Bend_T, E_Torsion_T
  use DistanceCriteria, only: Detailed_DistanceCriteria
  use EnergyCriteria, only: Detailed_EnergyCriteria
  use SimParameters, only: maxMol, distCriteria, minDistCriteria
  use E_Interface_HPS_piecewise_Diststore, only: Detailed_ECalc_HPS_piecewise_DStore
  use PairStorage, only: useDistStore
  implicit none
  real(dp), intent(inout) :: E_T              ! Total energy of the configuration
  logical, intent(inout) :: rejectMove        ! Flag to reject configuration if not a single cluster
  real(dp) :: pairList(1:maxMol, 1:maxMol)    ! Pairwise interaction data (e.g., energies or distances)

  ! Use stored pair distances if enabled
  if (useDistStore) then
    call Detailed_ECalc_HPS_piecewise_DStore(E_T, rejectMove)
    return
  endif

  ! Initialize total energy
  E_T = 0.0_dp

  ! Calculate intermolecular energies for HPS_piecewise force field
  call Detailed_ECalc_Inter(E_T, pairList)

  ! Check cluster criteria
  if (distCriteria .or. minDistCriteria) then
    ! Use distance-based or minimum distance-based criterion
    call Detailed_DistanceCriteria(pairList, rejectMove)
  else
    ! Use energy-based criterion
    call Detailed_EnergyCriteria(pairList, rejectMove)
  endif

  ! Calculate intramolecular non-bonded and bonded energies
  call Detailed_ECalc_IntraNonBonded(E_T)
  call Detailed_ECalc_BondStretch(E_T)
  call Detailed_ECalc_GasIntra

  ! Add bending and torsional energy contributions (currently zero for HPS_piecewise)
  E_Bend_T = 0.0_dp
  E_T = E_T + E_Bend_T
  E_Torsion_T = 0.0_dp
  E_T = E_T + E_Torsion_T
end subroutine Detailed_ECalc_HPS_piecewise
!=============================================================================      
! Calculates the energy change and checks cluster criteria for a molecular displacement
! (e.g., translation, rotation) in a Grand Canonical Monte Carlo nucleation simulation
! using the HPS_piecewise force field. Computes intermolecular (E_Inter) and intramolecular (E_Intra)
! energy differences for a single molecule’s displacement, ensuring the single-cluster
! requirement (distCriteria or minDistCriteria for long linear biomolecules). Supports
! useDistStore=.true. (delegates to Shift_ECalc_HPS_piecewise_DStore) and .false. (uses PairList).
! Populates NeighborDetailsNew for minDistCriteria when provided (e.g., LLTranslation).
subroutine Shift_ECalc_HPS_piecewise(E_Inter, E_Intra, disp, PairList, dETable, useIntra, rejMove, &
                                     useInter, NeighborDetailsNew)
  use VarPrecision, only: dp
  use SimParameters, only: beta, softCutOff, distCriteria, minDistCriteria, NTotal
  use Coords, only: MolArray
  use CBMC_Variables, only: regrowType
  use EnergyTables, only: E_Inter_Diff, E_NBond_Diff, E_Strch_Diff, E_Bend_Diff, E_Tors_Diff
  use IntraEnergy_HPS_piecewise, only: Shift_ECalc_IntraNonBonded
  use BondStretchFunctions, only: Shift_ECalc_BondStretch
  use BendingFunctions, only: Shift_ECalc_Bending
  use TorsionalFunctions, only: Shift_ECalc_Torsional
  use EnergyCriteria, only: Shift_EnergyCriteria
  use DistanceCriteria, only: Shift_DistanceCriteria
  use InterEnergy_HPS_piecewise, only: Shift_ECalc_Inter
  use E_Interface_HPS_piecewise_Diststore, only: Shift_ECalc_HPS_piecewise_DStore
  use PairStorage, only: useDistStore
  implicit none

  ! Input/Output variables (interface-defined)
  real(dp), intent(out) :: E_Inter, E_Intra
  real(dp), intent(inout) :: PairList(:), dETable(:)
  logical, intent(in) :: useIntra(1:4)
  logical, intent(inout) :: rejMove
  type(Displacement), intent(in) :: disp(:)
  logical, intent(in), optional :: useInter
  type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)

  ! Local variables
  logical :: interSwitch                     ! Flag to include intermolecular interactions
  integer :: nDisp                           ! Number of displaced atoms
  integer :: nIndx                           ! Global molecule index
  real(dp) :: E_NonBond, E_Stretch, E_Bend   ! Intramolecular energy components
  real(dp) :: E_Torsion, E_Improper          ! Intramolecular energy components

  ! Step 1: Delegate to distance storage version if useDistStore=.true.
  if (useDistStore) then
    call Shift_ECalc_HPS_piecewise_DStore(E_Inter, E_Intra, disp, PairList, dETable, useIntra, rejMove, &
                                          useInter, NeighborDetailsNew)
    return
  endif

  ! Step 2: Initialize variables
  nDisp = size(disp)
  rejMove = .false.
  E_Inter = 0.0_dp
  E_Intra = 0.0_dp
  E_NonBond = 0.0_dp
  E_Stretch = 0.0_dp
  E_Bend = 0.0_dp
  E_Torsion = 0.0_dp
  E_Improper = 0.0_dp
  dETable = 0.0_dp
  PairList = 0.0_dp
  E_Inter_Diff = 0.0_dp
  E_NBond_Diff = 0.0_dp
  E_Strch_Diff = 0.0_dp
  E_Bend_Diff = 0.0_dp
  E_Tors_Diff = 0.0_dp
  interSwitch = .true.
  if (present(useInter)) interSwitch = useInter

  ! Step 3: Calculate intermolecular energy and check for overlaps
  if (interSwitch .and. NTotal > 1) then
    if (present(NeighborDetailsNew)) then
      call Shift_ECalc_Inter(E_Inter, disp, PairList, dETable, rejMove, NeighborDetailsNew)
    else
      call Shift_ECalc_Inter(E_Inter, disp, PairList, dETable, rejMove)
    endif
    if (rejMove) return
    E_Inter_Diff = E_Inter

    ! Reject move if energy exceeds soft cutoff
    if (E_Inter * beta > softCutOff) then
      rejMove = .true.
      return
    endif

    ! Step 4: Check cluster criteria
    nIndx = MolArray(disp(1)%molType)%mol(disp(1)%molIndx)%indx
    if (distCriteria .or. minDistCriteria) then
      if (present(NeighborDetailsNew)) then
        call Shift_DistanceCriteria(PairList, nIndx, rejMove, NeighborDetailsNew)
      else
        call Shift_DistanceCriteria(PairList, nIndx, rejMove)
      endif
    else
      call Shift_EnergyCriteria(PairList, nIndx, rejMove)
    endif
    if (rejMove) return
  endif

  ! Step 5: Calculate intramolecular energy if applicable
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

  ! Step 6: Verify E_Intra for rigid moves
  if (.not. any(useIntra) .and. abs(E_Intra) > 1.0E-6_dp) then
    write(35, *) "Error: Non-zero E_Intra=", E_Intra, " in Shift_ECalc_HPS_piecewise with useIntra=.false."
    stop "Shift_ECalc_HPS_piecewise: Invalid intramolecular energy"
  endif
end subroutine Shift_ECalc_HPS_piecewise
!=============================================================================      
! Computes the intra- and intermolecular energy contributions for inserting a new molecule into the
! cluster using the HPS_piecewise force field in a Grand Canonical Monte Carlo nucleation simulation,
! as part of LLAVBMC_EBias_Rosen_In. Conforms to SwapInEnergyInterface, calculating E_Inter
! (intermolecular energy with the cluster) and E_Intra (intramolecular energy of the new molecule).
! Updates PairList and dETable for energy and cluster criteria, and NeighborDetailsNew for
! minDistCriteria when provided. Delegates to SwapIn_ECalc_HPS_piecewise_DStore if
! useDistStore=.true. Supports the minimum-distance cluster criterion (minDistCriteria) and is
! extensible to LJ_Q, Mpipi, HPS_single, and HPS_cation_pi. The HPS_piecewise potential is defined as:
! u_LJ(r) = 4 * eps * [ (sigma/r)**12 - (sigma/r)**6 ]
! U_HPS(r) = u_LJ(r) + eps * (1 - lambda) if r <= 2**(1/6) * sigma, else lambda * u_LJ(r)
! U_Electrostatic(r) = (q1 * q2 / r) * exp(-kappa * r)
subroutine SwapIn_ECalc_HPS_piecewise(E_Inter, E_Intra, PairList, dETable, rejMove, useInter, NeighborDetailsNew)
  use VarPrecision, only: dp
  use SimParameters, only: beta, softCutOff
  use Coords, only: newMol
  use CBMC_Variables, only: regrowType
  USE CoordinateTypes, ONLY: NeighborDetails
  use InterEnergy_HPS_piecewise, only: NewMol_ECalc_Inter
  use IntraEnergy_HPS_piecewise, only: NewMol_ECalc_IntraNonBonded
  use BondStretchFunctions, only: NewMol_ECalc_BondStretch
  use EnergyTables, only: E_Inter_Diff, E_NBond_Diff, E_Strch_Diff, E_Bend_Diff, E_Tors_Diff
  use E_Interface_HPS_piecewise_Diststore, only: SwapIn_ECalc_HPS_piecewise_DStore
  use PairStorage, only: useDistStore
  implicit none

  ! Interface-defined inputs/outputs
  real(dp), intent(out) :: E_Inter              ! Intermolecular energy with cluster
  real(dp), intent(out) :: E_Intra              ! Intramolecular energy of new molecule
  logical, intent(out) :: rejMove               ! Flag for invalid configurations (e.g., overlaps)
  logical, intent(in), optional :: useInter      ! If .false., skip intermolecular calculation
  type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)  ! Neighbor relationships
  real(dp), intent(inout) :: PairList(:)        ! Pairwise interaction list
  real(dp), intent(inout) :: dETable(:)         ! Energy difference table for cluster criteria

  ! Local variables
  logical :: interSwitch                        ! Flag to include intermolecular interactions
  real(dp) :: E_NonBond                        ! Non-bonded intramolecular energy
  real(dp) :: E_Stretch                        ! Bond stretching energy
  real(dp) :: E_Bend                           ! Bending energy (currently unused)
  real(dp) :: E_Torsion                        ! Torsional energy (currently unused)
  real(dp) :: E_Improper                       ! Improper angle energy (currently unused)

  ! Step 1: Delegate to distance storage version if useDistStore=.true.
  ! Overview: If distance storage is enabled, use the optimized version to compute energies.
  if (useDistStore) then
    call SwapIn_ECalc_HPS_piecewise_DStore(E_Inter, E_Intra, PairList, dETable, rejMove, &
                                           useInter, NeighborDetailsNew)
    return
  end if

  ! Step 2: Initialize variables
  ! Overview: Sets up energy accumulators, PairList, dETable, and interSwitch based on useInter.
  rejMove = .false.                             ! Initialize rejection flag
  interSwitch = .true.                          ! Default to computing intermolecular energy
  if (present(useInter)) interSwitch = useInter  ! Override if useInter provided
  E_Inter = 0.0_dp                              ! Initialize intermolecular energy
  E_Intra = 0.0_dp                              ! Initialize intramolecular energy
  E_NonBond = 0.0_dp                            ! Initialize non-bonded energy
  E_Stretch = 0.0_dp                            ! Initialize bond stretching energy
  E_Bend = 0.0_dp                               ! Initialize bending energy
  E_Torsion = 0.0_dp                            ! Initialize torsional energy
  E_Improper = 0.0_dp                           ! Initialize improper angle energy
  PairList = 0.0_dp                             ! Clear pairwise interaction list
  dETable = 0.0_dp                              ! Clear energy difference table
  E_Inter_Diff = 0.0_dp                         ! Clear global intermolecular energy difference
  E_NBond_Diff = 0.0_dp                         ! Clear global non-bonded energy difference
  E_Strch_Diff = 0.0_dp                         ! Clear global bond stretching energy difference
  E_Bend_Diff = 0.0_dp                          ! Clear global bending energy difference
  E_Tors_Diff = 0.0_dp                          ! Clear global torsional energy difference

  ! Step 3: Calculate intermolecular energy if applicable
  ! Overview: Computes intermolecular energy between the new molecule and the cluster using the
  ! HPS_piecewise force field (piecewise U_HPS with Debye-Hückel electrostatics). Updates PairList
  ! and dETable for energy and cluster criteria. If NeighborDetailsNew is present, updates neighbor
  ! relationships for minDistCriteria; otherwise, uses default neighbor handling. Rejects the move
  ! if overlaps occur or energy exceeds soft cutoff.
  if (interSwitch) then
    if (present(NeighborDetailsNew)) then
      call NewMol_ECalc_Inter(E_Inter, PairList, dETable, rejMove, NeighborDetailsNew)
    else
      call NewMol_ECalc_Inter(E_Inter, PairList, dETable, rejMove)
    end if
    if (rejMove) return
    E_Inter_Diff = E_Inter                      ! Update global intermolecular energy difference
    ! Reject move if energy exceeds soft cutoff
    if (E_Inter * beta > softCutOff) then
      rejMove = .true.
      return
    end if
  end if

  ! Step 4: Calculate intramolecular energy for non-rigid molecules
  ! Overview: Computes intramolecular energy components (non-bonded, bond stretching) for the
  ! new molecule if it is not rigid (regrowType /= 0). Uses HPS_piecewise non-bonded interactions.
  if (regrowType(newMol%molType) /= 0) then
    call NewMol_ECalc_IntraNonBonded(E_NonBond)  ! Non-bonded intramolecular energy
    call NewMol_ECalc_BondStretch(E_Stretch)     ! Bond stretching energy
    ! Note: Bending, torsional, and improper angle calculations are disabled
    E_Intra = E_NonBond + E_Stretch + E_Bend + E_Torsion + E_Improper  ! Sum intramolecular energies
    E_NBond_Diff = E_NonBond                    ! Update global non-bonded energy difference
    E_Strch_Diff = E_Stretch                    ! Update global bond stretching energy difference
    E_Bend_Diff = E_Bend                        ! Update global bending energy difference
    E_Tors_Diff = E_Torsion                     ! Update global torsional energy difference
  end if
end subroutine SwapIn_ECalc_HPS_piecewise
!=============================================================================
! Calculates the intermolecular and intramolecular energy of a molecule being removed in a Grand Canonical
! Monte Carlo swap-out move using the HPS_piecewise force field with useDistStore=.false. Computes E_Inter
! (piecewise Hydropathy Scale Lennard-Jones with switching at (sigma/r)^6 = 0.5 and Debye-Hückel electrostatics)
! and E_Intra (non-bonded, bond stretch, bending, torsional, and improper terms) for the specified molecule,
! updating dETable for cluster criteria. Optionally skips intermolecular interactions if useInter is present
! and .false. Used in AVBMC or similar moves.
subroutine SwapOut_ECalc_HPS_piecewise(E_Inter, E_Intra, nType, nMol, dETable, useInter)
  use VarPrecision
  use PairStorage, only: useDistStore
  use CBMC_Variables, only: regrowType
  use EnergyTables, only: E_Inter_Diff, E_NBond_Diff, E_Strch_Diff, E_Bend_Diff, E_Tors_Diff
  use InterEnergy_HPS_piecewise, only: Mol_ECalc_Inter
  use IntraEnergy_HPS_piecewise, only: Mol_ECalc_IntraNonBonded
  use BondStretchFunctions, only: Mol_ECalc_BondStretch
  use E_Interface_HPS_piecewise_Diststore, only: SwapOut_ECalc_HPS_piecewise_DStore
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

  ! Step 1: Check for distance storage
  ! Overview: Redirects to DStore version if useDistStore is true.
  if (useDistStore) then
    call SwapOut_ECalc_HPS_piecewise_DStore(E_Inter, E_Intra, nType, nMol, dETable, useInter)
    return
  end if

  ! Step 2: Initialize variables
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

  ! Step 3: Calculate intermolecular energy
  ! Overview: Computes piecewise HPS Lennard-Jones and Debye-Hückel electrostatic energies, negated for removal.
  if (interSwitch) then
    call Mol_ECalc_Inter(nType, nMol, dETable, E_Inter)
    E_Inter = -E_Inter
    E_Inter_Diff = E_Inter
  end if

  ! Step 4: Calculate intramolecular energy
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
end subroutine SwapOut_ECalc_HPS_piecewise
!=============================================================================
    
      
      end module
      
      
