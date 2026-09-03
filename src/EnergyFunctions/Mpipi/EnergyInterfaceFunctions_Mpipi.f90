!****************************************************************************************

!****************************************************************************************
      module E_Interface_Mpipi
      use CoordinateTypes
!=============================================================================
      contains
!============================================================================= 
! Calculates the total energy of the initial configuration for the Mpipi force field
! (designed for biomolecules like proteins and RNA) in a grand canonical ensemble
! nucleation simulation. Checks if the configuration forms a single cluster using one
! of three criteria: energy-based (interaction energy below threshold), distance-based
! (first atom distance below threshold), or minimum distance-based (any atom pair
! distance below threshold). Sets rejectMove to .true. if the cluster criterion fails.
! Uses stored pair distances if useDistStore=.true. for faster calculations.
!   Mpipi: Wang–Frankel + Debye-Hückel Electrostatics
!      - U_WF(r)           = ε * α * [ ( (σ / r)**(2μ) − 1 ) − ( (3σ / r)**(2μ) − 1 ) ]**(2ν)
!      - U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)     
subroutine Detailed_ECalc_Mpipi(E_T,rejMove)

  ! Calculates total energy using Mpipi force field for the initial configuration.
  ! Applies distance or energy-based cluster rejection criterion.
  ! Uses stored distances if available (useDistStore = .true.).

  use InterEnergy_Mpipi, only: Detailed_ECalc_Inter
  use IntraEnergy_Mpipi, only: Detailed_ECalc_IntraNonBonded, Detailed_ECalc_GasIntra
  use BondStretchFunctions, only: Detailed_ECalc_BondStretch
  use EnergyTables, only: E_Bend_T, E_Torsion_T
  use DistanceCriteria, only: Detailed_DistanceCriteria
  use EnergyCriteria, only: Detailed_EnergyCriteria
  use SimParameters, only: maxMol, distCriteria, minDistCriteria
  use E_Interface_Mpipi_Diststore, only: Detailed_ECalc_Mpipi_DStore
  use PairStorage, only: useDistStore

  implicit none

  logical, intent(inout) :: rejMove
  real(dp), intent(inout) :: E_T
  real(dp) :: PairList(1:maxMol,1:maxMol)

  ! Use precomputed distance storage if enabled
  if(useDistStore) then
    call Detailed_ECalc_Mpipi_DStore(E_T, rejMove)
    return
  endif

  ! Compute intermolecular energy
  E_T = 0.0_dp
  call Detailed_ECalc_Inter(E_T, PairList)

  ! Apply cluster rejection criteria
  if(distCriteria .or. minDistCriteria) then
    call Detailed_DistanceCriteria(PairList, rejMove)
  else
    call Detailed_EnergyCriteria(PairList, rejMove)
  endif

  ! Compute intramolecular contributions
  call Detailed_ECalc_IntraNonBonded(E_T)
  call Detailed_ECalc_BondStretch(E_T)
  call Detailed_ECalc_GasIntra

  ! Initialize and add bending energy
  E_Bend_T = 0.0_dp
  E_T = E_T + E_Bend_T

  ! Initialize and add torsional energy
  E_Torsion_T = 0.0_dp
  E_T = E_T + E_Torsion_T

  ! Optional: Uncomment if improper angle energy is implemented
  ! call Detailed_ECalc_Improper(E_T)

end subroutine Detailed_ECalc_Mpipi
!=============================================================================      
!     This function contains the energy and cluster criteria functions for any move
!     where a molecule is moved within a cluster.  This function takes a set of displacement
!     vectors (the "disp" variable) and returns the change in energy for both the Intra- and
!     Inter-molecular components.  This function can be used for for moves that any number of
!     atoms in a given molecule, but it can not be used if more than one molecule changes
!     in a given move.  
! Calculates the energy change and checks cluster criteria for a molecular displacement
! (e.g., translation, rotation) in a Grand Canonical Monte Carlo nucleation simulation
! using the Mpipi force field. Computes intermolecular (E_Inter) and intramolecular (E_Intra)
! energy differences for a single molecule’s displacement, ensuring the single-cluster
! requirement (distCriteria or minDistCriteria for long linear biomolecules). Supports
! useDistStore=.true. (delegates to Shift_ECalc_Mpipi_DStore) and .false. (uses PairList).
! Populates NeighborDetailsNew for minDistCriteria when provided (e.g., LLTranslation).
subroutine Shift_ECalc_Mpipi(E_Inter, E_Intra, disp, PairList, dETable, useIntra, rejMove, &
                             useInter, NeighborDetailsNew)

  ! Calculates the energy change and checks cluster criteria for a molecular displacement
  ! using the Mpipi force field. Handles both intra- and inter-molecular components.
  ! Uses distance- or energy-based rejection and supports distance caching.

  use VarPrecision, only: dp
  use SimParameters, only: distCriteria, minDistCriteria, beta, softCutOff, NTotal
  use BondStretchFunctions, only: Shift_ECalc_BondStretch
  use BendingFunctions, only: Shift_ECalc_Bending
  use TorsionalFunctions, only: Shift_ECalc_Torsional
  use CBMC_Variables, only: regrowType
  use Coords, only: MolArray
  use DistanceCriteria, only: Shift_DistanceCriteria
  use EnergyCriteria, only: Shift_EnergyCriteria
  use EnergyTables, only: E_Inter_Diff, E_NBond_Diff, E_Strch_Diff, E_Bend_Diff, E_Tors_Diff
  use InterEnergy_Mpipi, only: Shift_ECalc_Inter
  use IntraEnergy_Mpipi, only: Shift_ECalc_IntraNonBonded
  use E_Interface_Mpipi_Diststore, only: Shift_ECalc_Mpipi_DStore
  use PairStorage, only: useDistStore

  implicit none

  logical, intent(in), optional :: useInter
  type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)
  logical, intent(in) :: useIntra(1:4)
  type(Displacement), intent(in) :: disp(:)
  logical, intent(inout) :: rejMove
  real(dp), intent(inout) :: dETable(:)
  real(dp), intent(out) :: E_Intra, E_Inter
  real(dp), intent(inout) :: PairList(:)

  logical :: interSwitch
  integer :: nIndx, nDisp
  real(dp) :: E_NonBond, E_Stretch, E_Bend
  real(dp) :: E_Torsion, E_Improper

  ! If precomputed distances are used, call specialized version and return
  if(useDistStore) then
    call Shift_ECalc_Mpipi_DStore(E_Inter, E_Intra, disp, PairList, dETable, useIntra, rejMove, &
                                   useInter, NeighborDetailsNew)
    return
  endif

  nDisp        = size(disp)
  rejMove      = .false.
  E_Inter      = 0.0_dp
  E_Intra      = 0.0_dp
  E_NonBond    = 0.0_dp
  E_Stretch    = 0.0_dp
  E_Bend       = 0.0_dp
  E_Torsion    = 0.0_dp
  E_Improper   = 0.0_dp
  dETable      = 0.0_dp

  E_Inter_Diff = 0.0_dp
  E_NBond_Diff = 0.0_dp
  E_Strch_Diff = 0.0_dp
  E_Bend_Diff  = 0.0_dp
  E_Tors_Diff  = 0.0_dp

  ! Determine if intermolecular energy calculation is enabled
  if(present(useInter)) then
    interSwitch = useInter
  else
    interSwitch = .true.
  endif

  ! Intermolecular energy + cluster rejection criteria
  if(interSwitch) then
    if(NTotal > 1) then
      dETable = 0.0_dp
      PairList = 0.0_dp

      if (present(NeighborDetailsNew)) then
        call Shift_ECalc_Inter(E_Inter, disp, PairList, dETable, rejMove, NeighborDetailsNew)
      else
        call Shift_ECalc_Inter(E_Inter, disp, PairList, dETable, rejMove)
      endif

      if(rejMove) return
      E_Inter_Diff = E_Inter

      if(E_Inter * beta > softCutOff) then
        rejMove = .true.
        return
      endif
    endif

    ! Apply cluster criteria based on energy or distance
    nIndx = MolArray(disp(1)%molType)%mol(disp(1)%molIndx)%indx
    if(distCriteria .or. minDistCriteria) then
      call Shift_DistanceCriteria(PairList, nIndx, rejMove)
    else
      call Shift_EnergyCriteria(PairList, nIndx, rejMove)
    endif

    if(rejMove) return
  endif

  ! Intra-molecular energy terms (conditional on regrowType and useIntra flags)
  if(regrowType(disp(1)%molType) /= 0) then
    if(useIntra(1)) then
      call Shift_ECalc_IntraNonBonded(E_NonBond, disp)
      E_NBond_Diff = E_NonBond
    endif
    if(useIntra(2)) then
      call Shift_ECalc_BondStretch(E_Stretch, disp)
      E_Strch_Diff = E_Stretch
    endif
    if(useIntra(3)) then
      call Shift_ECalc_Bending(E_Bend, disp)
      E_Bend_Diff = E_Bend
    endif
    if(useIntra(4)) then
      call Shift_ECalc_Torsional(E_Torsion, disp)
      E_Tors_Diff = E_Torsion
    endif
    ! call Shift_ECalc_Improper(E_Improper, disp)  ! Optional
    E_Intra = E_NonBond + E_Stretch + E_Bend + E_Torsion + E_Improper
  endif

end subroutine Shift_ECalc_Mpipi
!=============================================================================      
! Computes the intra- and intermolecular energy contributions for inserting a new molecule into the 
! cluster using the Mpipi force field in a Grand Canonical Monte Carlo nucleation simulation, as 
! part of LLAVBMC_EBias_Rosen_In. Conforms to SwapInEnergyInterface, calculating E_Inter (intermolecular 
! energy with the cluster) and E_Intra (intramolecular energy of the new molecule). Updates PairList 
! and dETable for energy and cluster criteria, and NeighborDetailsNew for minDistCriteria when provided.
! Delegates to SwapIn_ECalc_Mpipi_DStore if useDistStore=.true. Supports the minimum-distance cluster 
! criterion (minDistCriteria) and is extensible to LJ_Q, HPS_single, HPS_piecewise, and HPS_cation_pi.
subroutine SwapIn_ECalc_Mpipi(E_Inter, E_Intra, PairList, dETable, rejMove, useInter, NeighborDetailsNew)

  ! Computes intra- and intermolecular energy contributions for inserting a new molecule
  ! using the Mpipi force field. Supports useDistStore acceleration and energy/distance
  ! cluster criteria for Grand Canonical Monte Carlo nucleation simulations.

  use VarPrecision, only: dp
  use Coords, only: newMol
  use CBMC_Variables, only: regrowType
  use InterEnergy_Mpipi, only: NewMol_ECalc_Inter
  use IntraEnergy_Mpipi, only: NewMol_ECalc_IntraNonBonded
  use BondStretchFunctions, only: NewMol_ECalc_BondStretch
  use EnergyTables, only: E_Inter_Diff, E_NBond_Diff, E_Strch_Diff, E_Bend_Diff, E_Tors_Diff
  use E_Interface_Mpipi_Diststore, only: SwapIn_ECalc_Mpipi_DStore
  use PairStorage, only: useDistStore

  implicit none

  logical, intent(out) :: rejMove
  logical, intent(in), optional :: useInter
  type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)
  real(dp), intent(out) :: E_Inter, E_Intra
  real(dp), intent(inout) :: PairList(:), dETable(:)

  logical :: interSwitch
  real(dp) :: E_NonBond, E_Stretch, E_Bend
  real(dp) :: E_Torsion, E_Improper

  ! Use precomputed distance version if enabled
  if(useDistStore) then
    call SwapIn_ECalc_Mpipi_DStore(E_Inter, E_Intra, PairList, dETable, rejMove, &
                                    useInter, NeighborDetailsNew)
    return
  endif

  rejMove     = .false.
  E_Inter     = 0.0_dp
  E_Intra     = 0.0_dp
  E_NonBond   = 0.0_dp
  E_Stretch   = 0.0_dp
  E_Bend      = 0.0_dp
  E_Torsion   = 0.0_dp
  E_Improper  = 0.0_dp
  PairList    = 0.0_dp
  dETable     = 0.0_dp

  E_Inter_Diff = 0.0_dp
  E_NBond_Diff = 0.0_dp
  E_Strch_Diff = 0.0_dp
  E_Bend_Diff  = 0.0_dp
  E_Tors_Diff  = 0.0_dp

  ! Determine if intermolecular interactions should be evaluated
  if(present(useInter)) then
    interSwitch = useInter
  else
    interSwitch = .true.
  endif

  ! Compute intermolecular energy and reject if overlapping
  if(interSwitch) then
    if (present(NeighborDetailsNew)) then
      call NewMol_ECalc_Inter(E_Inter, PairList, dETable, rejMove, NeighborDetailsNew)
    else
      call NewMol_ECalc_Inter(E_Inter, PairList, dETable, rejMove)
    endif
    if(rejMove) return
  endif

  ! Compute intramolecular energy terms
  if(regrowType(newMol%molType) /= 0) then
    call NewMol_ECalc_IntraNonBonded(E_NonBond)
    call NewMol_ECalc_BondStretch(E_Stretch)
    ! call NewMol_ECalc_Bending(E_Bend)
    ! call NewMol_ECalc_Torsional(E_Torsion)
    ! call NewMol_ECalc_Improper(E_Improper)

    E_Intra = E_NonBond + E_Stretch + E_Bend + E_Torsion + E_Improper
    E_NBond_Diff = E_NonBond
    E_Strch_Diff = E_Stretch
    E_Bend_Diff  = E_Bend
    E_Tors_Diff  = E_Torsion
  endif

  if(interSwitch) then
    E_Inter_Diff = E_Inter
  endif

end subroutine SwapIn_ECalc_Mpipi
!=============================================================================
! Calculates the intermolecular and intramolecular energy of a molecule being removed in a Grand Canonical
! Monte Carlo swap-out move using the Mpipi force field with useDistStore=.false. Computes E_Inter
! (Wang-Frenkel and Debye-Hückel electrostatics) and E_Intra (non-bonded, bond stretch, bending, torsional,
! and improper terms) for the specified molecule, updating dETable for cluster criteria. Optionally skips
! intermolecular interactions if useInter is present and .false. Used in AVBMC or similar moves.
subroutine SwapOut_ECalc_Mpipi(E_Inter, E_Intra, nType, nMol, dETable, useInter)

  ! Calculates inter- and intramolecular energy of a molecule being removed
  ! using the Mpipi force field in a GCMC swap-out move. Optionally disables
  ! intermolecular terms. Supports distance cache via useDistStore.

  use VarPrecision, only: dp
  use PairStorage, only: useDistStore
  use CBMC_Variables, only: regrowType
  use InterEnergy_Mpipi, only: Mol_ECalc_Inter
  use IntraEnergy_Mpipi, only: Mol_ECalc_IntraNonBonded
  use BondStretchFunctions, only: Mol_ECalc_BondStretch
  use EnergyTables, only: E_Inter_Diff, E_NBond_Diff, E_Strch_Diff, E_Bend_Diff, E_Tors_Diff
  use E_Interface_Mpipi_Diststore, only: SwapOut_ECalc_Mpipi_DStore

  implicit none

  logical, intent(in), optional :: useInter
  real(dp), intent(out) :: E_Inter, E_Intra
  integer, intent(in) :: nType, nMol
  real(dp), intent(inout) :: dETable(:)

  logical :: interSwitch
  real(dp) :: E_NonBond, E_Stretch, E_Bend
  real(dp) :: E_Torsion, E_Improper

  ! Use cached version if enabled
  if(useDistStore) then
    call SwapOut_ECalc_Mpipi_DStore(E_Inter, E_Intra, nType, nMol, dETable, useInter)
    return
  endif

  E_Inter     = 0.0_dp
  E_Intra     = 0.0_dp
  E_NonBond   = 0.0_dp
  E_Stretch   = 0.0_dp
  E_Bend      = 0.0_dp
  E_Torsion   = 0.0_dp
  E_Improper  = 0.0_dp

  E_Inter_Diff = 0.0_dp
  E_NBond_Diff = 0.0_dp
  E_Strch_Diff = 0.0_dp
  E_Bend_Diff  = 0.0_dp
  E_Tors_Diff  = 0.0_dp
  dETable      = 0.0_dp

  if(present(useInter)) then
    interSwitch = useInter
  else
    interSwitch = .true.
  endif

  ! Intermolecular energy calculation
  if(interSwitch) then
    call Mol_ECalc_Inter(nType, nMol, dETable, E_Inter)
    E_Inter = -E_Inter
  endif

  ! Intramolecular energy calculation
  if(regrowType(nType) /= 0) then
    call Mol_ECalc_IntraNonBonded(nType, nMol, E_NonBond)
    call Mol_ECalc_BondStretch(nType, nMol, E_Stretch)
    ! call Mol_ECalc_Bending(nType, nMol, E_Bend)
    ! call Mol_ECalc_Torsional(nType, nMol, E_Torsion)
    ! call Mol_ECalc_Improper(nType, nMol, E_Improper)

    E_Intra = E_NonBond + E_Stretch + E_Bend + E_Torsion + E_Improper
    E_Intra = -E_Intra
    E_NBond_Diff = -E_NonBond
    E_Strch_Diff = -E_Stretch
    E_Bend_Diff  = -E_Bend
    E_Tors_Diff  = -E_Torsion
  endif

  if(interSwitch) then
    E_Inter_Diff = E_Inter
  endif

end subroutine SwapOut_ECalc_Mpipi
!=============================================================================
    
      
      end module
      
      
