!===========================================================================================
      module SimpleMCMoves_Module
      contains
!===========================================================================================
! Performs a single atom translation move in a Grand Canonical Monte Carlo nucleation
! simulation. Randomly selects an atom from an active molecule in the cluster, applies a
! small translational displacement, and evaluates the move based on energy and cluster
! criteria (distance, minimum distance, or energy-based). Updates coordinates, energies,
! and neighbor lists if accepted. Supports distCriteria (distance-based) and minDistCriteria
! (minimum distance for long linear biomolecules). Handles useDistStore=.true. (all pair
! distances stored) and .false. (uses PairList for distances).
subroutine SingleAtom_Translation(E_T, acc_x, atmp_x)
  USE VarPrecision, ONLY: dp, atomIntType
  USE CoordinateTypes, ONLY: displacement
  USE SimParameters, ONLY: NTotal, NPart, beta, calcPressure, pressure, &
       distCriteria, maxMol, max_dist_single, prevMoveAccepted, pressure
  USE Coords, ONLY: MolArray
  USE ForceField, ONLY: nAtoms
  USE CBMC_Variables, ONLY: regrowType
  USE IndexingFunctions, ONLY: Get_MolIndex
  USE EnergyTables, ONLY: ETable
  USE EnergyPressurePointers, ONLY: Shift_ECalc, Shift_PCalc
  USE EnergyCriteria, ONLY: NeighborUpdate
  USE DistanceCriteria, ONLY: NeighborUpdate_Distance
  USE PairStorage, ONLY: UpdateDistArray, useDistStore
  USE E_Interface_LJ_Q, ONLY: Update_SubEnergies
  implicit none

  ! Input/Output variables (module-defined names preserved)
  real(dp), intent(inout) :: E_T          ! Total system energy
  real(dp), intent(inout) :: acc_x        ! Number of accepted moves
  real(dp), intent(inout) :: atmp_x       ! Number of attempted moves

  ! Local variables
  logical, parameter :: useIntra(1:4) = [.true., .true., .true., .true.]  ! Include intramolecular energies
  logical :: rejectMove                   ! Flag to reject move if cluster criteria fail
  integer :: nType, nMol, nIndx          ! Molecule type, instance, and global index
  integer :: nMove                       ! Index of selected molecule (1 to NTotal)
  integer :: selectedAtom                ! Selected atom index within molecule
  real(dp) :: dx, dy, dz                 ! Displacement components
  real(dp) :: E_Inter, E_Intra           ! Inter- and intramolecular energy changes
  real(dp) :: E_Diff                     ! Total energy change (inter + intra)
  real(dp) :: PairList(1:maxMol)         ! Pairwise distances or energies for neighbors
  real(dp) :: dETable(1:maxMol)          ! Energy table changes for neighbors
  real(dp) :: P_Diff                     ! Pressure change (if calcPressure=.true.)
  type(displacement) :: disp(1:1)         ! Displacement array for single atom
  real(dp), external :: grnd

  ! Initialize move
  atmp_x = atmp_x + 1.0_dp
  rejectMove = .false.
  prevMoveAccepted = .false.

  ! Step 1: Select a random molecule from the cluster
  nMove = floor(NTotal * grnd() + 1.0_dp)
  call Get_MolIndex(nMove, NPart, nType, nMol)
  nIndx = MolArray(nType)%mol(nMol)%indx

  ! Step 2: Check if molecule can be regrown (skip if regrowType=0)
  if (regrowType(nType) == 0) return

  ! Step 3: Select a random atom within the molecule
  selectedAtom = floor(nAtoms(nType) * grnd() + 1.0_dp)

  ! Step 4: Generate random translational displacement
  ! max_dist_single(nType): Maximum displacement for the atom (type-specific)
  dx = max_dist_single(nType) * (2.0_dp * grnd() - 1.0_dp)
  dy = max_dist_single(nType) * (2.0_dp * grnd() - 1.0_dp)
  dz = max_dist_single(nType) * (2.0_dp * grnd() - 1.0_dp)

  ! Step 5: Set up displacement structure
  disp(1)%molType = int(nType, atomIntType)
  disp(1)%molIndx = int(nMol, atomIntType)
  disp(1)%atmIndx = int(selectedAtom, atomIntType)
  disp(1)%x_old => MolArray(nType)%mol(nMol)%x(selectedAtom)
  disp(1)%y_old => MolArray(nType)%mol(nMol)%y(selectedAtom)
  disp(1)%z_old => MolArray(nType)%mol(nMol)%z(selectedAtom)
  disp(1)%x_new = disp(1)%x_old + dx
  disp(1)%y_new = disp(1)%y_old + dy
  disp(1)%z_new = disp(1)%z_old + dz

  ! Step 6: Calculate energy change and check cluster criteria
  E_Inter = 0.0_dp
  E_Intra = 0.0_dp
  call Shift_ECalc(E_Inter, E_Intra, disp(1:1), PairList, dETable, useIntra, rejectMove)
  if (rejectMove) return  ! Reject if cluster criteria (e.g., minDistCriteria) fail
  E_Diff = E_Inter + E_Intra

  ! Step 7: Apply Metropolis criterion for acceptance
  if (E_Diff <= 0.0_dp .or. -beta * E_Diff > log(grnd())) then
    ! Move accepted: Update coordinates, energies, and related quantities
    disp(1)%x_old = disp(1)%x_new
    disp(1)%y_old = disp(1)%y_new
    disp(1)%z_old = disp(1)%z_new
    E_T = E_T + E_Diff
    ETable = ETable + dETable
    acc_x = acc_x + 1.0_dp
    prevMoveAccepted = .true.

    ! Update pressure if required
    if (calcPressure) then
      call Shift_PCalc(P_Diff, disp(1:1))
      pressure = pressure + P_Diff
    endif

    ! Update distance storage if useDistStore=.true.
    if (useDistStore) then
      call UpdateDistArray
    endif

    ! Update neighbor lists based on cluster criterion
    if (distCriteria) then
      call NeighborUpdate_Distance(PairList, nIndx)
    else
      call NeighborUpdate(PairList, nIndx)
    endif

    ! Update sub-energies (e.g., LJ and Coulomb contributions)
    call Update_SubEnergies
  endif
end subroutine SingleAtom_Translation
!===========================================================================================
! Performs a translation move for a small molecule in a Grand Canonical Monte Carlo nucleation
! simulation. Randomly selects a molecule from the cluster, applies a uniform translational
! displacement to all its atoms, and evaluates the move based on intermolecular energy and
! cluster criteria (distance, minimum distance, or energy-based). Updates coordinates, energies,
! and neighbor lists if accepted. Designed for rigid small molecules (useIntra=.false.) but
! supports minDistCriteria for long linear biomolecules via Shift_ECalc and NeighborUpdate.
! Handles useDistStore=.true. (all pair distances stored) and .false. (uses PairList).
subroutine Translation(E_T, acc_x, atmp_x)
  USE VarPrecision, ONLY: dp, atomIntType
  USE CoordinateTypes, ONLY: displacement
  USE SimParameters, ONLY: NTotal, NPart, beta, calcPressure, pressure, &
       distCriteria, maxMol, maxAtoms, max_dist, P_Diff, prevMoveAccepted, pressure
  USE Coords, ONLY: MolArray
  USE ForceField, ONLY: nAtoms
  USE AcceptRates, ONLY: acptTrans, atmpTrans
  USE IndexingFunctions, ONLY: Get_MolIndex
  USE EnergyTables, ONLY: ETable
  USE EnergyPressurePointers, ONLY: Shift_ECalc, Shift_PCalc
  USE EnergyCriteria, ONLY: NeighborUpdate
  USE DistanceCriteria, ONLY: NeighborUpdate_Distance
  USE PairStorage, ONLY: UpdateDistArray, useDistStore
  USE UmbrellaSamplingNew, ONLY: useUmbrella, GetUmbrellaBias_Disp
  USE E_Interface_LJ_Q, ONLY: Update_SubEnergies
  implicit none

  ! Input/Output variables (module-defined)
  real(dp), intent(inout) :: E_T, acc_x, atmp_x

  ! Local variables
  logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]  ! Exclude intramolecular energies
  logical :: rejMove                                                  ! Flag to reject move if criteria fail
  integer :: iAtom, nType, nMol, nIndx, nMove                         ! Atom, type, molecule, and indices
  real(dp) :: biasDiff, biasEnergy                                    ! Umbrella bias and acceptance energy
  real(dp) :: dx, dy, dz                                              ! Displacement components
  real(dp) :: E_Inter, E_Intra                                        ! Inter- and intramolecular energies
  real(dp) :: PairList(1:maxMol)                                      ! Pairwise distances or energies
  real(dp) :: dETable(1:maxMol)                                       ! Energy table changes
  type(displacement) :: disp(1:maxAtoms)                              ! Displacement array for atoms
  real(dp), external :: grnd

  ! Initialize move
  prevMoveAccepted = .false.
  if (NTotal == 1) return  ! Skip if only one molecule (no relative movement)

  ! Step 1: Select a random molecule from the cluster
  nMove = floor(NTotal * grnd() + 1.0_dp)
  call Get_MolIndex(nMove, NPart, nType, nMol)
  nIndx = MolArray(nType)%mol(nMol)%indx

  ! Step 2: Increment attempt counters
  atmp_x = atmp_x + 1.0_dp
  atmpTrans(nType) = atmpTrans(nType) + 1.0_dp

  ! Step 3: Generate random translational displacement
  dx = max_dist(nType) * (2.0_dp * grnd() - 1.0_dp)
  dy = max_dist(nType) * (2.0_dp * grnd() - 1.0_dp)
  dz = max_dist(nType) * (2.0_dp * grnd() - 1.0_dp)

  ! Step 4: Set up displacement for all atoms in the molecule
  do iAtom = 1, nAtoms(nType)
    disp(iAtom)%molType = int(nType, atomIntType)
    disp(iAtom)%molIndx = int(nMol, atomIntType)
    disp(iAtom)%atmIndx = int(iAtom, atomIntType)
    disp(iAtom)%Displaced = .true.
    disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
    disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
    disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
    disp(iAtom)%x_new = disp(iAtom)%x_old + dx
    disp(iAtom)%y_new = disp(iAtom)%y_old + dy
    disp(iAtom)%z_new = disp(iAtom)%z_old + dz
  enddo

  ! Step 5: Calculate energy change and check cluster criteria
  E_Inter = 0.0_dp
  E_Intra = 0.0_dp
  call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
  if (rejMove) return  ! Reject if cluster criteria (e.g., minDistCriteria) fail
  if (abs(E_Intra) > 1.0E-6_dp) then
    write(35, *) "Warning: Non-zero E_Intra=", E_Intra, " in Translation with useIntra=.false."
  endif

  ! Step 6: Handle umbrella sampling bias
  biasDiff = 0.0_dp
  if (useUmbrella) then
    call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
    if (rejMove) return  ! Reject if bias calculation fails
  endif
  biasEnergy = beta * E_Inter - biasDiff  ! Exclude E_Intra (rigid molecule)

  ! Step 7: Apply Metropolis criterion for acceptance
  if (biasEnergy <= 0.0_dp .or. -biasEnergy > log(grnd())) then
    ! Move accepted: Update coordinates, energies, and related quantities
    do iAtom = 1, nAtoms(nType)
      disp(iAtom)%x_old = disp(iAtom)%x_new
      disp(iAtom)%y_old = disp(iAtom)%y_new
      disp(iAtom)%z_old = disp(iAtom)%z_new
    enddo
    E_T = E_T + E_Inter
    ETable = ETable + dETable
    acc_x = acc_x + 1.0_dp
    acptTrans(nType) = acptTrans(nType) + 1.0_dp
    prevMoveAccepted = .true.

    ! Update pressure if required
    if (calcPressure) then
      call Shift_PCalc(P_Diff, disp(1:nAtoms(nType)))
      pressure = pressure + P_Diff
    endif

    ! Update distance storage if useDistStore=.true.
    if (useDistStore) then
      call UpdateDistArray
    endif

    ! Update neighbor lists based on cluster criterion
    if (distCriteria) then
      call NeighborUpdate_Distance(PairList, nIndx)
    else
      call NeighborUpdate(PairList, nIndx)
    endif

    ! Update sub-energies (e.g., LJ and Coulomb contributions)
    call Update_SubEnergies
  endif
end subroutine Translation
!===========================================================================================
! Performs a translation move for a long linear biomolecule in a Grand Canonical Monte Carlo
! nucleation simulation. Randomly selects a molecule from the cluster, applies a uniform
! translational displacement to all its atoms, and evaluates the move based on intermolecular
! energy and minimum distance criterion (minDistCriteria). Updates coordinates, energies,
! and neighbor lists (NeighborList, NeighborPairs) if accepted. Supports useDistStore=.true.
! (all pair distances stored) and .false. (uses PairList). Excludes intramolecular energies
! (useIntra=.false.) since translation does not alter relative atom positions within the molecule.
subroutine LLTranslation(E_T, acc_x, atmp_x)
  USE VarPrecision, ONLY: dp, atomIntType
  USE CoordinateTypes, ONLY: displacement, NeighborDetails
  USE SimParameters, ONLY: NTotal, NPart, beta, calcPressure, pressure, &
       distCriteria, minDistCriteria, maxMol, maxAtoms, max_dist, P_Diff, prevMoveAccepted, pressure
  USE Coords, ONLY: MolArray
  USE ForceField, ONLY: nAtoms
  USE AcceptRates, ONLY: acptTrans, atmpTrans
  USE IndexingFunctions, ONLY: Get_MolIndex
  USE EnergyTables, ONLY: ETable
  USE EnergyPressurePointers, ONLY: Shift_ECalc, Shift_PCalc
  USE EnergyCriteria, ONLY: NeighborUpdate
  USE DistanceCriteria, ONLY: NeighborUpdate_Distance
  USE PairStorage, ONLY: UpdateDistArray, useDistStore
  USE UmbrellaSamplingNew, ONLY: useUmbrella, GetUmbrellaBias_Disp
  USE E_Interface_LJ_Q, ONLY: Update_SubEnergies
  implicit none

  ! Input/Output variables (module-defined)
  real(dp), intent(inout) :: E_T, acc_x, atmp_x

  ! Local variables
  logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]  ! No intramolecular energies
  logical :: rejMove                                                  ! Flag to reject move if criteria fail
  integer :: iAtom, nType, nMol, nIndx, nMove, i                     ! Atom, type, molecule, and indices
  real(dp) :: biasDiff, biasEnergy                                    ! Umbrella bias and acceptance energy
  real(dp) :: dx, dy, dz                                              ! Displacement components
  real(dp) :: E_Inter, E_Intra                                        ! Inter- and intramolecular energies
  real(dp) :: PairList(1:maxMol)                                      ! Pairwise minimum distances
  real(dp) :: dETable(1:maxMol)                                       ! Energy table changes
  type(displacement) :: disp(1:maxAtoms)                              ! Displacement array for atoms
  type(NeighborDetails) :: NeighborDetailsNew(1:maxMol)                ! New neighbor details for minDistCriteria
  real(dp), external :: grnd

  ! Initialize move
  prevMoveAccepted = .false.
  if (NTotal == 1) return  ! Skip if only one molecule (no relative movement)

  ! Step 1: Allocate NeighborDetailsNew pairIndices
  do i = 1, maxMol
    NeighborDetailsNew(i)%nPairs = 0
    NeighborDetailsNew(i)%pairIndices = 0
  enddo

  ! Step 2: Select a random molecule from the cluster
  nMove = floor(NTotal * grnd() + 1.0_dp)
  call Get_MolIndex(nMove, NPart, nType, nMol)
  nIndx = MolArray(nType)%mol(nMol)%indx

  ! Step 3: Increment attempt counters
  atmp_x = atmp_x + 1.0_dp
  atmpTrans(nType) = atmpTrans(nType) + 1.0_dp

  ! Step 4: Generate random translational displacement
  dx = max_dist(nType) * (2.0_dp * grnd() - 1.0_dp)
  dy = max_dist(nType) * (2.0_dp * grnd() - 1.0_dp)
  dz = max_dist(nType) * (2.0_dp * grnd() - 1.0_dp)

  ! Step 5: Set up displacement for all atoms in the molecule
  do iAtom = 1, nAtoms(nType)
    disp(iAtom)%molType = int(nType, atomIntType)
    disp(iAtom)%molIndx = int(nMol, atomIntType)
    disp(iAtom)%atmIndx = int(iAtom, atomIntType)
    disp(iAtom)%Displaced = .true.
    disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
    disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
    disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
    disp(iAtom)%x_new = disp(iAtom)%x_old + dx
    disp(iAtom)%y_new = disp(iAtom)%y_old + dy
    disp(iAtom)%z_new = disp(iAtom)%z_old + dz
  enddo

  ! Step 6: Calculate energy change and check minDistCriteria
  E_Inter = 0.0_dp
  E_Intra = 0.0_dp
  call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, &
                   rejMove, .true., NeighborDetailsNew)
  if (rejMove) return  ! Reject if minDistCriteria or other criteria fail
  if (abs(E_Intra) > 1.0E-6_dp) then
    write(35, *) "Error: Non-zero E_Intra=", E_Intra, " in LLTranslation with useIntra=.false."
    stop "LLTranslation: Invalid intramolecular energy"
  endif

  ! Step 7: Handle umbrella sampling bias
  biasDiff = 0.0_dp
  if (useUmbrella) then
    call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
    if (rejMove) return  ! Reject if bias calculation fails
  endif
  biasEnergy = beta * E_Inter - biasDiff  ! Use only E_Inter (rigid translation)

  ! Step 8: Apply Metropolis criterion for acceptance
  if (biasEnergy <= 0.0_dp .or. -biasEnergy > log(grnd())) then
    ! Move accepted: Update coordinates, energies, and related quantities
    do iAtom = 1, nAtoms(nType)
      disp(iAtom)%x_old = disp(iAtom)%x_new
      disp(iAtom)%y_old = disp(iAtom)%y_new
      disp(iAtom)%z_old = disp(iAtom)%z_new
    enddo
    E_T = E_T + E_Inter
    ETable = ETable + dETable
    acc_x = acc_x + 1.0_dp
    acptTrans(nType) = acptTrans(nType) + 1.0_dp
    prevMoveAccepted = .true.

    ! Update pressure if required
    if (calcPressure) then
      call Shift_PCalc(P_Diff, disp(1:nAtoms(nType)))
      pressure = pressure + P_Diff
    endif

    ! Update distance storage if useDistStore=.true.
    if (useDistStore) then
      call UpdateDistArray
    endif

    ! Update neighbor lists based on cluster criterion
    if (distCriteria .or. minDistCriteria) then
      call NeighborUpdate_Distance(PairList, nIndx, NeighborDetailsNew)
    else
      call NeighborUpdate(PairList, nIndx)
    endif

    ! Update sub-energies (e.g., LJ and Coulomb contributions)
    call Update_SubEnergies
  endif


end subroutine LLTranslation
!===========================================================================================
! Performs a rotational move for a small molecule in a Grand Canonical Monte Carlo nucleation
! simulation. Randomly selects a molecule from the cluster and applies a rotation around one of
! three planes (xy, xz, or yz) with equal probability. Evaluates the move based on intermolecular
! energy and cluster criteria (distance or energy-based). Updates coordinates, energies, and neighbor
! lists if accepted. Designed for rigid small molecules (LJ_Q, Mpipi) with useIntra=.false. Supports
! useDistStore=.true. and minDistCriteria via Shift_ECalc and NeighborUpdate. Called for rotational
! moves in the simulation, complementing Translation.
subroutine Rotation(E_T, acc_x, atmp_x)
  use SimParameters, only: NTotal, prevMoveAccepted
  use VarPrecision, only: dp
  implicit none

  ! Input/Output variables
  real(dp), intent(inout) :: E_T    ! Total energy
  real(dp), intent(inout) :: acc_x  ! Acceptance counter
  real(dp), intent(inout) :: atmp_x ! Attempt counter

  ! Local variables
  real(dp) :: ran_num  ! Random number for selecting rotation type
  real(dp), external :: grnd

  ! Initialize move
  prevMoveAccepted = .false.
  if (NTotal == 1) return  ! Skip if only one molecule (no relative rotation)

  ! Select rotation type with equal probability (1/3 each)
  ran_num = grnd()
  if (ran_num < 1.0_dp / 3.0_dp) then
    call Rot_xy(E_T, acc_x, atmp_x)
  elseif (ran_num < 2.0_dp / 3.0_dp) then
    call Rot_xz(E_T, acc_x, atmp_x)
  else
    call Rot_yz(E_T, acc_x, atmp_x)
  endif
end subroutine Rotation
!=======================================================      
! Performs a rotational move around the xy-plane for a small molecule in a Grand Canonical Monte Carlo
! nucleation simulation. Randomly selects a molecule from the cluster, rotates all its atoms around the
! center of mass, and evaluates the move based on intermolecular energy and cluster criteria (distance
! or energy-based). Updates coordinates, energies, and neighbor lists if accepted. Designed for rigid
! small molecules (LJ_Q, Mpipi) with useIntra=.false. Supports useDistStore=.true. and minDistCriteria
! via Shift_ECalc and NeighborUpdate. Called by Rotation for rotational moves.
subroutine Rot_xy(E_T, acc_x, atmp_x)
  USE VarPrecision, ONLY: dp, atomIntType
  USE CoordinateTypes, ONLY: displacement
  USE SimParameters, ONLY: NTotal, NPart, maxAtoms, maxMol, max_rot, beta, &
       calcPressure, distCriteria, isActive, P_Diff, prevMoveAccepted, pressure
  USE Coords, ONLY: MolArray
  USE Forcefield, ONLY: nAtoms, atomArray, atomData, totalMass
  USE AcceptRates, ONLY: atmpRot, acptRot
  USE IndexingFunctions, ONLY: Get_MolIndex
  USE EnergyTables, ONLY: ETable
  USE EnergyPressurePointers, ONLY: Shift_ECalc, Shift_PCalc
  USE EnergyCriteria, ONLY: NeighborUpdate
  USE DistanceCriteria, ONLY: NeighborUpdate_Distance
  USE PairStorage, ONLY: UpdateDistArray, useDistStore
  USE UmbrellaSamplingNew, ONLY: GetUmbrellaBias_Disp, useUmbrella
  USE E_Interface_LJ_Q, ONLY: Update_SubEnergies
  implicit none

  ! Input/Output variables
  real(dp), intent(inout) :: E_T    ! Total energy
  real(dp), intent(inout) :: acc_x  ! Acceptance counter
  real(dp), intent(inout) :: atmp_x ! Attempt counter

  ! Local variables
  logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]  ! No intramolecular energies
  logical :: rejMove                                                  ! Flag to reject move if criteria fail
  integer :: iAtom, nMove, nType, nMol, nIndx, atmType                ! Indices
  real(dp) :: angle, c_term, s_term                                   ! Rotation angle and trigonometric terms
  real(dp) :: xcm, ycm, x_scale, y_scale                              ! Center of mass and scaling factors
  real(dp) :: E_Inter, E_Intra                                        ! Inter- and intramolecular energies
  real(dp) :: biasDiff, biasEnergy                                    ! Umbrella bias and acceptance energy
  real(dp) :: PairList(maxMol)                                        ! Pairwise distances or energies
  real(dp) :: dETable(maxMol)                                         ! Energy table changes
  type(displacement) :: disp(maxAtoms)                                ! Displacement array for atoms
  real(dp), external :: grnd

  ! Initialize move
  prevMoveAccepted = .false.
  if (NTotal == 1) return  ! Skip if only one molecule (no relative rotation)

  ! Step 1: Select a random molecule from the cluster
  nMove = floor(NTotal * grnd() + 1.0_dp)
  call Get_MolIndex(nMove, NPart, nType, nMol)
  nIndx = MolArray(nType)%mol(nMol)%indx
  if (nAtoms(nType) == 1 .or. .not. isActive(nIndx)) return  ! Skip single-atom or inactive molecules

  ! Step 2: Increment attempt counters
  atmp_x = atmp_x + 1.0_dp
  atmpRot(nType) = atmpRot(nType) + 1

  ! Step 3: Set up displacement array
  do iAtom = 1, nAtoms(nType)
    disp(iAtom)%molType = int(nType, atomIntType)
    disp(iAtom)%molIndx = nMol
    disp(iAtom)%atmIndx = iAtom
    disp(iAtom)%Displaced = .true.
    disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
    disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
    disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
  enddo

  ! Step 4: Generate random rotational displacement
  angle = max_rot(nType) * (2.0_dp * grnd() - 1.0_dp)
  c_term = cos(angle)
  s_term = sin(angle)

  ! Step 5: Calculate center of mass as pivot point
  xcm = 0.0_dp
  ycm = 0.0_dp
  do iAtom = 1, nAtoms(nType)
    atmType = atomArray(nType, iAtom)
    xcm = xcm + atomData(atmType)%mass * disp(iAtom)%x_old
    ycm = ycm + atomData(atmType)%mass * disp(iAtom)%y_old
  enddo
  xcm = xcm / totalMass(nType)
  ycm = ycm / totalMass(nType)

  ! Step 6: Apply rotation to all atoms
  do iAtom = 1, nAtoms(nType)
    x_scale = disp(iAtom)%x_old - xcm
    y_scale = disp(iAtom)%y_old - ycm
    disp(iAtom)%x_new = c_term * x_scale - s_term * y_scale + xcm
    disp(iAtom)%y_new = s_term * x_scale + c_term * y_scale + ycm
    disp(iAtom)%z_new = disp(iAtom)%z_old
  enddo

  ! Step 7: Calculate energy change and check cluster criteria
  E_Inter = 0.0_dp
  E_Intra = 0.0_dp
  call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
  if (rejMove) return
  if (abs(E_Intra) > 1.0E-6_dp) then
    write(35, *) "Warning: Non-zero E_Intra=", E_Intra, " in Rot_xy with useIntra=.false."
  endif

  ! Step 8: Handle umbrella sampling bias
  biasDiff = 0.0_dp
  if (useUmbrella) then
    call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
    if (rejMove) return
  endif
  biasEnergy = beta * E_Inter - biasDiff

  ! Step 9: Apply Metropolis criterion for acceptance
  if (biasEnergy <= 0.0_dp .or. -biasEnergy > log(grnd())) then
    ! Move accepted: Update coordinates, energies, and related quantities
    do iAtom = 1, nAtoms(nType)
      disp(iAtom)%x_old = disp(iAtom)%x_new
      disp(iAtom)%y_old = disp(iAtom)%y_new
      disp(iAtom)%z_old = disp(iAtom)%z_new
    enddo
    E_T = E_T + E_Inter
    ETable = ETable + dETable
    acc_x = acc_x + 1.0_dp
    acptRot(nType) = acptRot(nType) + 1
    prevMoveAccepted = .true.

    ! Update pressure if required
    if (calcPressure) then
      call Shift_PCalc(P_Diff, disp(1:nAtoms(nType)))
      pressure = pressure + P_Diff
    endif

    ! Update distance storage if useDistStore=.true.
    if (useDistStore) then
      call UpdateDistArray
    endif

    ! Update neighbor lists based on cluster criterion
    if (distCriteria) then
      call NeighborUpdate_Distance(PairList, nIndx)
    else
      call NeighborUpdate(PairList, nIndx)
    endif

    ! Update sub-energies
    call Update_SubEnergies
  endif
end subroutine Rot_xy
!=======================================================      
! Performs a rotational move around the xz-plane for a small molecule in a Grand Canonical Monte Carlo
! nucleation simulation. Randomly selects a molecule from the cluster, rotates all its atoms around the
! center of mass, and evaluates the move based on intermolecular energy and cluster criteria (distance
! or energy-based). Updates coordinates, energies, and neighbor lists if accepted. Designed for rigid
! small molecules (LJ_Q, Mpipi) with useIntra=.false. Supports useDistStore=.true. and minDistCriteria
! via Shift_ECalc and NeighborUpdate. Called by Rotation for rotational moves.
subroutine Rot_xz(E_T, acc_x, atmp_x)
  USE VarPrecision, ONLY: dp, atomIntType
  USE CoordinateTypes, ONLY: displacement
  USE SimParameters, ONLY: NTotal, NPart, maxAtoms, maxMol, max_rot, beta, &
       calcPressure, distCriteria, isActive, P_Diff, prevMoveAccepted, pressure
  USE Coords, ONLY: MolArray
  USE Forcefield, ONLY: nAtoms, atomArray, atomData, totalMass
  USE AcceptRates, ONLY: atmpRot, acptRot
  USE IndexingFunctions, ONLY: Get_MolIndex
  USE EnergyTables, ONLY: ETable
  USE EnergyPressurePointers, ONLY: Shift_ECalc, Shift_PCalc
  USE EnergyCriteria, ONLY: NeighborUpdate
  USE DistanceCriteria, ONLY: NeighborUpdate_Distance
  USE PairStorage, ONLY: UpdateDistArray, useDistStore
  USE UmbrellaSamplingNew, ONLY: GetUmbrellaBias_Disp, useUmbrella
  USE E_Interface_LJ_Q, ONLY: Update_SubEnergies
  implicit none

  ! Input/Output variables
  real(dp), intent(inout) :: E_T    ! Total energy
  real(dp), intent(inout) :: acc_x  ! Acceptance counter
  real(dp), intent(inout) :: atmp_x ! Attempt counter

  ! Local variables
  logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]  ! No intramolecular energies
  logical :: rejMove                                                  ! Flag to reject move if criteria fail
  integer :: iAtom, nMove, nType, nMol, nIndx, atmType                ! Indices
  real(dp) :: angle, c_term, s_term                                   ! Rotation angle and trigonometric terms
  real(dp) :: xcm, zcm, x_scale, z_scale                              ! Center of mass and scaling factors
  real(dp) :: E_Inter, E_Intra                                        ! Inter- and intramolecular energies
  real(dp) :: biasDiff, biasEnergy                                    ! Umbrella bias and acceptance energy
  real(dp) :: PairList(maxMol)                                        ! Pairwise distances or energies
  real(dp) :: dETable(maxMol)                                         ! Energy table changes
  type(displacement) :: disp(maxAtoms)                                ! Displacement array for atoms
  real(dp), external :: grnd

  ! Initialize move
  prevMoveAccepted = .false.
  if (NTotal == 1) return  ! Skip if only one molecule (no relative rotation)

  ! Step 1: Select a random molecule from the cluster
  nMove = floor(NTotal * grnd() + 1.0_dp)
  call Get_MolIndex(nMove, NPart, nType, nMol)
  nIndx = MolArray(nType)%mol(nMol)%indx
  if (nAtoms(nType) == 1 .or. .not. isActive(nIndx)) return  ! Skip single-atom or inactive molecules

  ! Step 2: Increment attempt counters
  atmp_x = atmp_x + 1.0_dp
  atmpRot(nType) = atmpRot(nType) + 1

  ! Step 3: Set up displacement array
  do iAtom = 1, nAtoms(nType)
    disp(iAtom)%molType = int(nType, atomIntType)
    disp(iAtom)%molIndx = nMol
    disp(iAtom)%atmIndx = iAtom
    disp(iAtom)%Displaced = .true.
    disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
    disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
    disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
  enddo

  ! Step 4: Generate random rotational displacement
  angle = max_rot(nType) * (2.0_dp * grnd() - 1.0_dp)
  c_term = cos(angle)
  s_term = sin(angle)

  ! Step 5: Calculate center of mass as pivot point
  xcm = 0.0_dp
  zcm = 0.0_dp
  do iAtom = 1, nAtoms(nType)
    atmType = atomArray(nType, iAtom)
    xcm = xcm + atomData(atmType)%mass * disp(iAtom)%x_old
    zcm = zcm + atomData(atmType)%mass * disp(iAtom)%z_old
  enddo
  xcm = xcm / totalMass(nType)
  zcm = zcm / totalMass(nType)

  ! Step 6: Apply rotation to all atoms
  do iAtom = 1, nAtoms(nType)
    x_scale = disp(iAtom)%x_old - xcm
    z_scale = disp(iAtom)%z_old - zcm
    disp(iAtom)%x_new = c_term * x_scale - s_term * z_scale + xcm
    disp(iAtom)%z_new = s_term * x_scale + c_term * z_scale + zcm
    disp(iAtom)%y_new = disp(iAtom)%y_old
  enddo

  ! Step 7: Calculate energy change and check cluster criteria
  E_Inter = 0.0_dp
  E_Intra = 0.0_dp
  call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
  if (rejMove) return
  if (abs(E_Intra) > 1.0E-6_dp) then
    write(35, *) "Warning: Non-zero E_Intra=", E_Intra, " in Rot_xz with useIntra=.false."
  endif

  ! Step 8: Handle umbrella sampling bias
  biasDiff = 0.0_dp
  if (useUmbrella) then
    call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
    if (rejMove) return
  endif
  biasEnergy = beta * E_Inter - biasDiff

  ! Step 9: Apply Metropolis criterion for acceptance
  if (biasEnergy <= 0.0_dp .or. -biasEnergy > log(grnd())) then
    ! Move accepted: Update coordinates, energies, and related quantities
    do iAtom = 1, nAtoms(nType)
      disp(iAtom)%x_old = disp(iAtom)%x_new
      disp(iAtom)%y_old = disp(iAtom)%y_new
      disp(iAtom)%z_old = disp(iAtom)%z_new
    enddo
    E_T = E_T + E_Inter
    ETable = ETable + dETable
    acc_x = acc_x + 1.0_dp
    acptRot(nType) = acptRot(nType) + 1
    prevMoveAccepted = .true.

    ! Update pressure if required
    if (calcPressure) then
      call Shift_PCalc(P_Diff, disp(1:nAtoms(nType)))
      pressure = pressure + P_Diff
    endif

    ! Update distance storage if useDistStore=.true.
    if (useDistStore) then
      call UpdateDistArray
    endif

    ! Update neighbor lists based on cluster criterion
    if (distCriteria) then
      call NeighborUpdate_Distance(PairList, nIndx)
    else
      call NeighborUpdate(PairList, nIndx)
    endif

    ! Update sub-energies
    call Update_SubEnergies
  endif
end subroutine Rot_xz
!=======================================================      
! Performs a rotational move around the yz-plane for a small molecule in a Grand Canonical Monte Carlo
! nucleation simulation. Randomly selects a molecule from the cluster, rotates all its atoms around the
! center of mass, and evaluates the move based on intermolecular energy and cluster criteria (distance
! or energy-based). Updates coordinates, energies, and neighbor lists if accepted. Designed for rigid
! small molecules (LJ_Q, Mpipi) with useIntra=.false. Supports useDistStore=.true. and minDistCriteria
! via Shift_ECalc and NeighborUpdate. Called by Rotation for rotational moves.
subroutine Rot_yz(E_T, acc_x, atmp_x)
  USE VarPrecision, ONLY: dp, atomIntType
  USE CoordinateTypes, ONLY: displacement
  USE SimParameters, ONLY: NTotal, NPart, maxAtoms, maxMol, max_rot, beta, &
       calcPressure, distCriteria, isActive, P_Diff, prevMoveAccepted, pressure
  USE Coords, ONLY: MolArray
  USE Forcefield, ONLY: nAtoms, atomArray, atomData, totalMass
  USE AcceptRates, ONLY: atmpRot, acptRot
  USE IndexingFunctions, ONLY: Get_MolIndex
  USE EnergyTables, ONLY: ETable
  USE EnergyPressurePointers, ONLY: Shift_ECalc, Shift_PCalc
  USE EnergyCriteria, ONLY: NeighborUpdate
  USE DistanceCriteria, ONLY: NeighborUpdate_Distance
  USE PairStorage, ONLY: UpdateDistArray, useDistStore
  USE UmbrellaSamplingNew, ONLY: GetUmbrellaBias_Disp, useUmbrella
  USE E_Interface_LJ_Q, ONLY: Update_SubEnergies
  implicit none

  ! Input/Output variables
  real(dp), intent(inout) :: E_T    ! Total energy
  real(dp), intent(inout) :: acc_x  ! Acceptance counter
  real(dp), intent(inout) :: atmp_x ! Attempt counter

  ! Local variables
  logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]  ! No intramolecular energies
  logical :: rejMove                                                  ! Flag to reject move if criteria fail
  integer :: iAtom, nMove, nType, nMol, nIndx, atmType                ! Indices
  real(dp) :: angle, c_term, s_term                                   ! Rotation angle and trigonometric terms
  real(dp) :: ycm, zcm, y_scale, z_scale                              ! Center of mass and scaling factors
  real(dp) :: E_Inter, E_Intra                                        ! Inter- and intramolecular energies
  real(dp) :: biasDiff, biasEnergy                                    ! Umbrella bias and acceptance energy
  real(dp) :: PairList(maxMol)                                        ! Pairwise distances or energies
  real(dp) :: dETable(maxMol)                                         ! Energy table changes
  type(displacement) :: disp(maxAtoms)                                ! Displacement array for atoms
  real(dp), external :: grnd

  ! Initialize move
  prevMoveAccepted = .false.
  if (NTotal == 1) return  ! Skip if only one molecule (no relative rotation)

  ! Step 1: Select a random molecule from the cluster
  nMove = floor(NTotal * grnd() + 1.0_dp)
  call Get_MolIndex(nMove, NPart, nType, nMol)
  nIndx = MolArray(nType)%mol(nMol)%indx
  if (nAtoms(nType) == 1 .or. .not. isActive(nIndx)) return  ! Skip single-atom or inactive molecules

  ! Step 2: Increment attempt counters
  atmp_x = atmp_x + 1.0_dp
  atmpRot(nType) = atmpRot(nType) + 1

  ! Step 3: Set up displacement array
  do iAtom = 1, nAtoms(nType)
    disp(iAtom)%molType = int(nType, atomIntType)
    disp(iAtom)%molIndx = nMol
    disp(iAtom)%atmIndx = iAtom
    disp(iAtom)%Displaced = .true.
    disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
    disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
    disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
  enddo

  ! Step 4: Generate random rotational displacement
  angle = max_rot(nType) * (2.0_dp * grnd() - 1.0_dp)
  c_term = cos(angle)
  s_term = sin(angle)

  ! Step 5: Calculate center of mass as pivot point
  ycm = 0.0_dp
  zcm = 0.0_dp
  do iAtom = 1, nAtoms(nType)
    atmType = atomArray(nType, iAtom)
    ycm = ycm + atomData(atmType)%mass * disp(iAtom)%y_old
    zcm = zcm + atomData(atmType)%mass * disp(iAtom)%z_old
  enddo
  ycm = ycm / totalMass(nType)
  zcm = zcm / totalMass(nType)

  ! Step 6: Apply rotation to all atoms
  do iAtom = 1, nAtoms(nType)
    y_scale = disp(iAtom)%y_old - ycm
    z_scale = disp(iAtom)%z_old - zcm
    disp(iAtom)%y_new = c_term * y_scale - s_term * z_scale + ycm
    disp(iAtom)%z_new = s_term * y_scale + c_term * z_scale + zcm
    disp(iAtom)%x_new = disp(iAtom)%x_old
  enddo

  ! Step 7: Calculate energy change and check cluster criteria
  E_Inter = 0.0_dp
  E_Intra = 0.0_dp
  call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
  if (rejMove) return
  if (abs(E_Intra) > 1.0E-6_dp) then
    write(35, *) "Warning: Non-zero E_Intra=", E_Intra, " in Rot_yz with useIntra=.false."
  endif

  ! Step 8: Handle umbrella sampling bias
  biasDiff = 0.0_dp
  if (useUmbrella) then
    call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
    if (rejMove) return
  endif
  biasEnergy = beta * E_Inter - biasDiff

  ! Step 9: Apply Metropolis criterion for acceptance
  if (biasEnergy <= 0.0_dp .or. -biasEnergy > log(grnd())) then
    ! Move accepted: Update coordinates, energies, and related quantities
    do iAtom = 1, nAtoms(nType)
      disp(iAtom)%x_old = disp(iAtom)%x_new
      disp(iAtom)%y_old = disp(iAtom)%y_new
      disp(iAtom)%z_old = disp(iAtom)%z_new
    enddo
    E_T = E_T + E_Inter
    ETable = ETable + dETable
    acc_x = acc_x + 1.0_dp
    acptRot(nType) = acptRot(nType) + 1
    prevMoveAccepted = .true.

    ! Update pressure if required
    if (calcPressure) then
      call Shift_PCalc(P_Diff, disp(1:nAtoms(nType)))
      pressure = pressure + P_Diff
    endif

    ! Update distance storage if useDistStore=.true.
    if (useDistStore) then
      call UpdateDistArray
    endif

    ! Update neighbor lists based on cluster criterion
    if (distCriteria) then
      call NeighborUpdate_Distance(PairList, nIndx)
    else
      call NeighborUpdate(PairList, nIndx)
    endif

    ! Update sub-energies
    call Update_SubEnergies
  endif
end subroutine Rot_yz
!=======================================================      
! Performs a rotational move for a long linear biomolecule in a Grand Canonical Monte Carlo nucleation
! simulation. Randomly selects a molecule from the cluster and applies a rotation around one of three
! planes (xy, xz, or yz) with equal probability. Evaluates the move based on intermolecular energy and
! minimum distance criterion (minDistCriteria). Updates coordinates, energies, and neighbor lists
! (NeighborList, NeighborPairs) if accepted. Designed for HPS force fields (HPS_single, HPS_piecewise,
! HPS_cation_pi) with useIntra=.false. Supports useDistStore=.true. via Shift_ECalc and
! NeighborUpdate_Distance. Called for rotational moves, complementing LLTranslation.
subroutine LLRotation(E_T, acc_x, atmp_x)
  use SimParameters, only: NTotal, prevMoveAccepted
  use VarPrecision, only: dp
  implicit none

  ! Input/Output variables
  real(dp), intent(inout) :: E_T    ! Total energy
  real(dp), intent(inout) :: acc_x  ! Acceptance counter
  real(dp), intent(inout) :: atmp_x ! Attempt counter

  ! Local variables
  real(dp) :: grnd, ran_num  ! Random number for selecting rotation type

  ! Initialize move
  prevMoveAccepted = .false.
  if (NTotal == 1) return  ! Skip if only one molecule (no relative rotation)

  ! Select rotation type with equal probability (1/3 each)
  ran_num = grnd()
  if (ran_num < 1.0_dp / 3.0_dp) then
    call LLRot_xy(E_T, acc_x, atmp_x)
  elseif (ran_num < 2.0_dp / 3.0_dp) then
    call LLRot_xz(E_T, acc_x, atmp_x)
  else
    call LLRot_yz(E_T, acc_x, atmp_x)
  endif
end subroutine LLRotation

!=======================================================      
! Performs a rotational move around the xy-plane for a long linear biomolecule in a Grand Canonical
! Monte Carlo nucleation simulation. Randomly selects a molecule from the cluster, rotates all its
! atoms around the center of mass, and evaluates the move based on intermolecular energy and minimum
! distance criterion (minDistCriteria). Updates coordinates, energies, and neighbor lists (NeighborList,
! NeighborPairs) if accepted. Designed for HPS force fields (HPS_single, HPS_piecewise, HPS_cation_pi)
! with useIntra=.false. Supports useDistStore=.true. via Shift_ECalc and NeighborUpdate_Distance.
! Called by LLRotation for rotational moves.
subroutine LLRot_xy(E_T, acc_x, atmp_x)
  USE VarPrecision, ONLY: dp, atomIntType
  USE CoordinateTypes, ONLY: displacement, NeighborDetails
  USE SimParameters, ONLY: NTotal, NPart, maxAtoms, maxMol, max_rot, beta, &
       calcPressure, distCriteria, minDistCriteria, &
       isActive, P_Diff, prevMoveAccepted, pressure
  USE Coords, ONLY: MolArray
  USE Forcefield, ONLY: nAtoms, atomArray, atomData, totalMass
  USE AcceptRates, ONLY: atmpRot, acptRot
  USE IndexingFunctions, ONLY: Get_MolIndex
  USE EnergyTables, ONLY: ETable
  USE EnergyPressurePointers, ONLY: Shift_ECalc, Shift_PCalc
  USE EnergyCriteria, ONLY: NeighborUpdate
  USE DistanceCriteria, ONLY: NeighborUpdate_Distance
  USE PairStorage, ONLY: UpdateDistArray, useDistStore
  USE UmbrellaSamplingNew, ONLY: GetUmbrellaBias_Disp, useUmbrella
  USE E_Interface_LJ_Q, ONLY: Update_SubEnergies
  implicit none

  ! Input/Output variables
  real(dp), intent(inout) :: E_T    ! Total energy
  real(dp), intent(inout) :: acc_x  ! Acceptance counter
  real(dp), intent(inout) :: atmp_x ! Attempt counter

  ! Local variables
  logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]  ! No intramolecular energies
  logical :: rejMove                                                  ! Flag to reject move if criteria fail
  integer :: iAtom, nMove, nType, nMol, nIndx, atmType, i             ! Indices
  real(dp) :: angle, c_term, s_term                                   ! Rotation angle and trigonometric terms
  real(dp) :: xcm, ycm, x_scale, y_scale                              ! Center of mass and scaling factors
  real(dp) :: E_Inter, E_Intra                                        ! Inter- and intramolecular energies
  real(dp) :: biasDiff, biasEnergy                                    ! Umbrella bias and acceptance energy
  real(dp) :: PairList(maxMol)                                        ! Pairwise minimum distances
  real(dp) :: dETable(maxMol)                                         ! Energy table changes
  type(displacement) :: disp(maxAtoms)                                ! Displacement array for atoms
  type(NeighborDetails) :: NeighborDetailsNew(1:maxMol)                  ! New neighbor details for minDistCriteria
  real(dp), external :: grnd

  ! Initialize move
  prevMoveAccepted = .false.
  if (NTotal == 1) return  ! Skip if only one molecule (no relative rotation)

  ! Step 1: Allocate NeighborDetailsNew pairIndices
  do i = 1, maxMol
    NeighborDetailsNew(i)%nPairs = 0
    NeighborDetailsNew(i)%pairIndices = 0
  enddo

  ! Step 2: Select a random molecule from the cluster
  nMove = floor(NTotal * grnd() + 1.0_dp)
  call Get_MolIndex(nMove, NPart, nType, nMol)
  nIndx = MolArray(nType)%mol(nMol)%indx
  if (nAtoms(nType) == 1 .or. .not. isActive(nIndx)) then
    return  ! Skip single-atom or inactive molecules
  endif

  ! Step 3: Increment attempt counters
  atmp_x = atmp_x + 1.0_dp
  atmpRot(nType) = atmpRot(nType) + 1

  ! Step 4: Set up displacement array
  do iAtom = 1, nAtoms(nType)
    disp(iAtom)%molType = int(nType, atomIntType)
    disp(iAtom)%molIndx = nMol
    disp(iAtom)%atmIndx = iAtom
    disp(iAtom)%Displaced = .true.
    disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
    disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
    disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
  enddo

  ! Step 5: Generate random rotational displacement
  angle = max_rot(nType) * (2.0_dp * grnd() - 1.0_dp)
  c_term = cos(angle)
  s_term = sin(angle)

  ! Step 6: Calculate center of mass as pivot point
  xcm = 0.0_dp
  ycm = 0.0_dp
  do iAtom = 1, nAtoms(nType)
    atmType = atomArray(nType, iAtom)
    xcm = xcm + atomData(atmType)%mass * disp(iAtom)%x_old
    ycm = ycm + atomData(atmType)%mass * disp(iAtom)%y_old
  enddo
  xcm = xcm / totalMass(nType)
  ycm = ycm / totalMass(nType)

  ! Step 7: Apply rotation to all atoms
  do iAtom = 1, nAtoms(nType)
    x_scale = disp(iAtom)%x_old - xcm
    y_scale = disp(iAtom)%y_old - ycm
    disp(iAtom)%x_new = c_term * x_scale - s_term * y_scale + xcm
    disp(iAtom)%y_new = s_term * x_scale + c_term * y_scale + ycm
    disp(iAtom)%z_new = disp(iAtom)%z_old
  enddo

  ! Step 8: Calculate energy change and check minDistCriteria
  E_Inter = 0.0_dp
  E_Intra = 0.0_dp
  call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, &
                   rejMove, .true., NeighborDetailsNew)
  if (rejMove) return
  if (abs(E_Intra) > 1.0E-6_dp) then
    write(35, *) "Error: Non-zero E_Intra=", E_Intra, " in LLRot_xy with useIntra=.false."
    return
  endif

  ! Step 9: Handle umbrella sampling bias
  biasDiff = 0.0_dp
  if (useUmbrella) then
    call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
    if (rejMove) return
  endif
  biasEnergy = beta * E_Inter - biasDiff

  ! Step 10: Apply Metropolis criterion for acceptance
  if (biasEnergy <= 0.0_dp .or. -biasEnergy > log(grnd())) then
    ! Move accepted: Update coordinates, energies, and related quantities
    do iAtom = 1, nAtoms(nType)
      disp(iAtom)%x_old = disp(iAtom)%x_new
      disp(iAtom)%y_old = disp(iAtom)%y_new
      disp(iAtom)%z_old = disp(iAtom)%z_new
    enddo
    E_T = E_T + E_Inter
    ETable = ETable + dETable
    acc_x = acc_x + 1.0_dp
    acptRot(nType) = acptRot(nType) + 1
    prevMoveAccepted = .true.

    ! Update pressure if required
    if (calcPressure) then
      call Shift_PCalc(P_Diff, disp(1:nAtoms(nType)))
      pressure = pressure + P_Diff
    endif

    ! Update distance storage if useDistStore=.true.
    if (useDistStore) then
      call UpdateDistArray
    endif

    ! Update neighbor lists based on cluster criterion
    if (distCriteria .or. minDistCriteria) then
      call NeighborUpdate_Distance(PairList, nIndx, NeighborDetailsNew)
    else
      call NeighborUpdate(PairList, nIndx)
    endif

    ! Update sub-energies
    call Update_SubEnergies
  endif

end subroutine LLRot_xy
!=======================================================      
! Performs a rotational move around the xz-plane for a long linear biomolecule in a Grand Canonical
! Monte Carlo nucleation simulation. Randomly selects a molecule from the cluster, rotates all its
! atoms around the center of mass, and evaluates the move based on intermolecular energy and minimum
! distance criterion (minDistCriteria). Updates coordinates, energies, and neighbor lists (NeighborList,
! NeighborPairs) if accepted. Designed for HPS force fields (HPS_single, HPS_piecewise, HPS_cation_pi)
! with useIntra=.false. Supports useDistStore=.true. via Shift_ECalc and NeighborUpdate_Distance.
! Called by LLRotation for rotational moves.
subroutine LLRot_xz(E_T, acc_x, atmp_x)
  USE VarPrecision, ONLY: dp, atomIntType
  USE CoordinateTypes, ONLY: displacement, NeighborDetails
  USE SimParameters, ONLY: NTotal, NPart, maxAtoms, maxMol, max_rot, beta, &
       calcPressure, distCriteria, minDistCriteria, &
       isActive, P_Diff, prevMoveAccepted, pressure
  USE Coords, ONLY: MolArray
  USE Forcefield, ONLY: nAtoms, atomArray, atomData, totalMass
  USE AcceptRates, ONLY: atmpRot, acptRot
  USE IndexingFunctions, ONLY: Get_MolIndex
  USE EnergyTables, ONLY: ETable
  USE EnergyPressurePointers, ONLY: Shift_ECalc, Shift_PCalc
  USE EnergyCriteria, ONLY: NeighborUpdate
  USE DistanceCriteria, ONLY: NeighborUpdate_Distance
  USE PairStorage, ONLY: UpdateDistArray, useDistStore
  USE UmbrellaSamplingNew, ONLY: GetUmbrellaBias_Disp, useUmbrella
  USE E_Interface_LJ_Q, ONLY: Update_SubEnergies
  implicit none

  ! Input/Output variables
  real(dp), intent(inout) :: E_T    ! Total energy
  real(dp), intent(inout) :: acc_x  ! Acceptance counter
  real(dp), intent(inout) :: atmp_x ! Attempt counter

  ! Local variables
  logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]  ! No intramolecular energies
  logical :: rejMove                                                  ! Flag to reject move if criteria fail
  integer :: iAtom, nMove, nType, nMol, nIndx, atmType, i             ! Indices
  real(dp) :: angle, c_term, s_term                                   ! Rotation angle and trigonometric terms
  real(dp) :: xcm, zcm, x_scale, z_scale                              ! Center of mass and scaling factors
  real(dp) :: E_Inter, E_Intra                                        ! Inter- and intramolecular energies
  real(dp) :: biasDiff, biasEnergy                                    ! Umbrella bias and acceptance energy
  real(dp) :: PairList(maxMol)                                        ! Pairwise minimum distances
  real(dp) :: dETable(maxMol)                                         ! Energy table changes
  type(displacement) :: disp(maxAtoms)                                ! Displacement array for atoms
  type(NeighborDetails) :: NeighborDetailsNew(1:maxMol)                  ! New neighbor details for minDistCriteria
  real(dp), external :: grnd

  ! Initialize move
  prevMoveAccepted = .false.
  if (NTotal == 1) return  ! Skip if only one molecule (no relative rotation)

  ! Step 1: Allocate NeighborDetailsNew pairIndices
  do i = 1, maxMol
    NeighborDetailsNew(i)%nPairs = 0
    NeighborDetailsNew(i)%pairIndices = 0
  enddo

  ! Step 2: Select a random molecule from the cluster
  nMove = floor(NTotal * grnd() + 1.0_dp)
  call Get_MolIndex(nMove, NPart, nType, nMol)
  nIndx = MolArray(nType)%mol(nMol)%indx
  if (nAtoms(nType) == 1 .or. .not. isActive(nIndx)) then
    return  ! Skip single-atom or inactive molecules
  endif

  ! Step 3: Increment attempt counters
  atmp_x = atmp_x + 1.0_dp
  atmpRot(nType) = atmpRot(nType) + 1

  ! Step 4: Set up displacement array
  do iAtom = 1, nAtoms(nType)
    disp(iAtom)%molType = int(nType, atomIntType)
    disp(iAtom)%molIndx = nMol
    disp(iAtom)%atmIndx = iAtom
    disp(iAtom)%Displaced = .true.
    disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
    disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
    disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
  enddo

  ! Step 5: Generate random rotational displacement
  angle = max_rot(nType) * (2.0_dp * grnd() - 1.0_dp)
  c_term = cos(angle)
  s_term = sin(angle)

  ! Step 6: Calculate center of mass as pivot point
  xcm = 0.0_dp
  zcm = 0.0_dp
  do iAtom = 1, nAtoms(nType)
    atmType = atomArray(nType, iAtom)
    xcm = xcm + atomData(atmType)%mass * disp(iAtom)%x_old
    zcm = zcm + atomData(atmType)%mass * disp(iAtom)%z_old
  enddo
  xcm = xcm / totalMass(nType)
  zcm = zcm / totalMass(nType)

  ! Step 7: Apply rotation to all atoms
  do iAtom = 1, nAtoms(nType)
    x_scale = disp(iAtom)%x_old - xcm
    z_scale = disp(iAtom)%z_old - zcm
    disp(iAtom)%x_new = c_term * x_scale - s_term * z_scale + xcm
    disp(iAtom)%z_new = s_term * x_scale + c_term * z_scale + zcm
    disp(iAtom)%y_new = disp(iAtom)%y_old
  enddo

  ! Step 8: Calculate energy change and check minDistCriteria
  E_Inter = 0.0_dp
  E_Intra = 0.0_dp
  call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, &
                   rejMove, .true., NeighborDetailsNew)
  if (rejMove) return
  if (abs(E_Intra) > 1.0E-6_dp) then
    write(35, *) "Error: Non-zero E_Intra=", E_Intra, " in LLRot_xz with useIntra=.false."
    return
  endif

  ! Step 9: Handle umbrella sampling bias
  biasDiff = 0.0_dp
  if (useUmbrella) then
    call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
    if (rejMove) return
  endif
  biasEnergy = beta * E_Inter - biasDiff

  ! Step 10: Apply Metropolis criterion for acceptance
  if (biasEnergy <= 0.0_dp .or. -biasEnergy > log(grnd())) then
    ! Move accepted: Update coordinates, energies, and related quantities
    do iAtom = 1, nAtoms(nType)
      disp(iAtom)%x_old = disp(iAtom)%x_new
      disp(iAtom)%y_old = disp(iAtom)%y_new
      disp(iAtom)%z_old = disp(iAtom)%z_new
    enddo
    E_T = E_T + E_Inter
    ETable = ETable + dETable
    acc_x = acc_x + 1.0_dp
    acptRot(nType) = acptRot(nType) + 1
    prevMoveAccepted = .true.

    ! Update pressure if required
    if (calcPressure) then
      call Shift_PCalc(P_Diff, disp(1:nAtoms(nType)))
      pressure = pressure + P_Diff
    endif

    ! Update distance storage if useDistStore=.true.
    if (useDistStore) then
      call UpdateDistArray
    endif

    ! Update neighbor lists based on cluster criterion
    if (distCriteria .or. minDistCriteria) then
      call NeighborUpdate_Distance(PairList, nIndx, NeighborDetailsNew)
    else
      call NeighborUpdate(PairList, nIndx)
    endif

    ! Update sub-energies
    call Update_SubEnergies
  endif

end subroutine LLRot_xz
!=======================================================      
! Performs a rotational move around the yz-plane for a long linear biomolecule in a Grand Canonical
! Monte Carlo nucleation simulation. Randomly selects a molecule from the cluster, rotates all its
! atoms around the center of mass, and evaluates the move based on intermolecular energy and minimum
! distance criterion (minDistCriteria). Updates coordinates, energies, and neighbor lists (NeighborList,
! NeighborPairs) if accepted. Designed for HPS force fields (HPS_single, HPS_piecewise, HPS_cation_pi)
! with useIntra=.false. Supports useDistStore=.true. via Shift_ECalc and NeighborUpdate_Distance.
! Called by LLRotation for rotational moves.
subroutine LLRot_yz(E_T, acc_x, atmp_x)
  USE VarPrecision, ONLY: dp, atomIntType
  USE CoordinateTypes, ONLY: displacement, NeighborDetails
  USE SimParameters, ONLY: NTotal, NPart, maxAtoms, maxMol, max_rot, beta, &
       calcPressure, distCriteria, minDistCriteria, &
       isActive, P_Diff, pressure, prevMoveAccepted
  USE Coords, ONLY: MolArray
  USE Forcefield, ONLY: nAtoms, atomArray, atomData, totalMass
  USE AcceptRates, ONLY: atmpRot, acptRot
  USE IndexingFunctions, ONLY: Get_MolIndex
  USE EnergyTables, ONLY: ETable
  USE EnergyPressurePointers, ONLY: Shift_ECalc, Shift_PCalc
  USE EnergyCriteria, ONLY: NeighborUpdate
  USE DistanceCriteria, ONLY: NeighborUpdate_Distance
  USE PairStorage, ONLY: UpdateDistArray, useDistStore
  USE UmbrellaSamplingNew, ONLY: GetUmbrellaBias_Disp, useUmbrella
  USE E_Interface_LJ_Q, ONLY: Update_SubEnergies
  implicit none

  ! Input/Output variables
  real(dp), intent(inout) :: E_T    ! Total energy
  real(dp), intent(inout) :: acc_x  ! Acceptance counter
  real(dp), intent(inout) :: atmp_x ! Attempt counter

  ! Local variables
  logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]  ! No intramolecular energies
  logical :: rejMove                                                  ! Flag to reject move if criteria fail
  integer :: iAtom, nMove, nType, nMol, nIndx, atmType, i             ! Indices
  real(dp) :: angle, c_term, s_term                                   ! Rotation angle and trigonometric terms
  real(dp) :: ycm, zcm, y_scale, z_scale                              ! Center of mass and scaling factors
  real(dp) :: E_Inter, E_Intra                                        ! Inter- and intramolecular energies
  real(dp) :: biasDiff, biasEnergy                                    ! Umbrella bias and acceptance energy
  real(dp) :: PairList(maxMol)                                        ! Pairwise minimum distances
  real(dp) :: dETable(maxMol)                                         ! Energy table changes
  type(displacement) :: disp(maxAtoms)                                ! Displacement array for atoms
  type(NeighborDetails) :: NeighborDetailsNew(1:maxMol)                  ! New neighbor details for minDistCriteria
  real(dp), external :: grnd

  ! Initialize move
  prevMoveAccepted = .false.
  if (NTotal == 1) return  ! Skip if only one molecule (no relative rotation)

  ! Step 1: Allocate NeighborDetailsNew pairIndices
  do i = 1, maxMol
    NeighborDetailsNew(i)%nPairs = 0
    NeighborDetailsNew(i)%pairIndices = 0
  enddo

  ! Step 2: Select a random molecule from the cluster
  nMove = floor(NTotal * grnd() + 1.0_dp)
  call Get_MolIndex(nMove, NPart, nType, nMol)
  nIndx = MolArray(nType)%mol(nMol)%indx
  if (nAtoms(nType) == 1 .or. .not. isActive(nIndx)) then
    return  ! Skip single-atom or inactive molecules
  endif

  ! Step 3: Increment attempt counters
  atmp_x = atmp_x + 1.0_dp
  atmpRot(nType) = atmpRot(nType) + 1

  ! Step 4: Set up displacement array
  do iAtom = 1, nAtoms(nType)
    disp(iAtom)%molType = int(nType, atomIntType)
    disp(iAtom)%molIndx = nMol
    disp(iAtom)%atmIndx = iAtom
    disp(iAtom)%Displaced = .true.
    disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
    disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
    disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
  enddo

  ! Step 5: Generate random rotational displacement
  angle = max_rot(nType) * (2.0_dp * grnd() - 1.0_dp)
  c_term = cos(angle)
  s_term = sin(angle)

  ! Step 6: Calculate center of mass as pivot point
  ycm = 0.0_dp
  zcm = 0.0_dp
  do iAtom = 1, nAtoms(nType)
    atmType = atomArray(nType, iAtom)
    ycm = ycm + atomData(atmType)%mass * disp(iAtom)%y_old
    zcm = zcm + atomData(atmType)%mass * disp(iAtom)%z_old
  enddo
  ycm = ycm / totalMass(nType)
  zcm = zcm / totalMass(nType)

  ! Step 7: Apply rotation to all atoms
  do iAtom = 1, nAtoms(nType)
    y_scale = disp(iAtom)%y_old - ycm
    z_scale = disp(iAtom)%z_old - zcm
    disp(iAtom)%y_new = c_term * y_scale - s_term * z_scale + ycm
    disp(iAtom)%z_new = s_term * y_scale + c_term * z_scale + zcm
    disp(iAtom)%x_new = disp(iAtom)%x_old
  enddo

  ! Step 8: Calculate energy change and check minDistCriteria
  E_Inter = 0.0_dp
  E_Intra = 0.0_dp
  call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, &
                   rejMove, .true., NeighborDetailsNew)
  if (rejMove) return
  if (abs(E_Intra) > 1.0E-6_dp) then
    write(35, *) "Error: Non-zero E_Intra=", E_Intra, " in LLRot_yz with useIntra=.false."
    return
  endif

  ! Step 9: Handle umbrella sampling bias
  biasDiff = 0.0_dp
  if (useUmbrella) then
    call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
    if (rejMove) return
  endif
  biasEnergy = beta * E_Inter - biasDiff

  ! Step 10: Apply Metropolis criterion for acceptance
  if (biasEnergy <= 0.0_dp .or. -biasEnergy > log(grnd())) then
    ! Move accepted: Update coordinates, energies, and related quantities
    do iAtom = 1, nAtoms(nType)
      disp(iAtom)%x_old = disp(iAtom)%x_new
      disp(iAtom)%y_old = disp(iAtom)%y_new
      disp(iAtom)%z_old = disp(iAtom)%z_new
    enddo
    E_T = E_T + E_Inter
    ETable = ETable + dETable
    acc_x = acc_x + 1.0_dp
    acptRot(nType) = acptRot(nType) + 1
    prevMoveAccepted = .true.

    ! Update pressure if required
    if (calcPressure) then
      call Shift_PCalc(P_Diff, disp(1:nAtoms(nType)))
      pressure = pressure + P_Diff
    endif

    ! Update distance storage if useDistStore=.true.
    if (useDistStore) then
      call UpdateDistArray
    endif

    ! Update neighbor lists based on cluster criterion
    if (distCriteria .or. minDistCriteria) then
      call NeighborUpdate_Distance(PairList, nIndx, NeighborDetailsNew)
    else
      call NeighborUpdate(PairList, nIndx)
    endif

    ! Update sub-energies
    call Update_SubEnergies
  endif

end subroutine LLRot_yz
!=======================================================      
! Performs a temperature change move in a Grand Canonical Monte Carlo nucleation
! simulation. Proposes a new system temperature within a small range, evaluates the
! move using a Metropolis-like acceptance criterion, and accounts for umbrella sampling
! biases if enabled. Updates temperature and beta (1/temperature) if accepted. Intended
! for systems with a single molecule (NTotal=1), as it does not affect molecular
! coordinates or cluster criteria (e.g., minDistCriteria for long linear biomolecules).
! Note: This is an experimental move and may not guarantee accurate results.
subroutine TemperatureMove(E_T, acc_x, atmp_x)
  USE VarPrecision, ONLY: dp
  USE SimParameters, ONLY: NTotal, temperature, beta
  USE UmbrellaSamplingNew, ONLY: useUmbrella, GetUmbrellaBias_Temperature
  implicit none

  ! Input/Output variables (module-defined names preserved)
  real(dp), intent(inout) :: E_T      ! Total system energy
  real(dp), intent(inout) :: acc_x    ! Number of accepted moves
  real(dp), intent(inout) :: atmp_x   ! Number of attempted moves

  ! Local variables
  real(dp), parameter :: power = 3.0_dp  ! Exponent for acceptance ratio (6/2 = 3)
  logical :: rejectMove                  ! Flag to reject move due to umbrella bias
  real(dp) :: TempNew                    ! Proposed new temperature
  real(dp) :: betaNew                    ! Proposed new inverse temperature (1/TempNew)
  real(dp) :: biasOld, biasNew           ! Umbrella biases for old and new temperatures
  real(dp) :: biasDiff                   ! Difference in umbrella biases (new - old)
  real(dp) :: biasEnergy                 ! Energy term for acceptance criterion
  real(dp), external :: grnd
  
  ! Step 1: Check if system has exactly one molecule
  ! Temperature move is only valid for NTotal=1, as it does not alter coordinates
  if (NTotal /= 1) return

  ! Step 2: Increment attempt counter and propose new temperature
  atmp_x = atmp_x + 1.0_dp
  TempNew = temperature + 15.0_dp * (2.0_dp * grnd() - 1.0_dp)  ! Random change in [-15, +15]
  betaNew = 1.0_dp / TempNew

  ! Step 3: Handle umbrella sampling bias if enabled
  rejectMove = .false.
  biasOld = 0.0_dp
  biasNew = 0.0_dp
  if (useUmbrella) then
    call GetUmbrellaBias_Temperature(biasDiff, rejectMove)
    if (rejectMove) return  ! Reject if umbrella bias calculation fails
  else
    biasDiff = 0.0_dp
  endif

  ! Step 4: Calculate acceptance energy term
  ! biasEnergy = (betaNew - beta) * E_T - biasDiff, where biasDiff = biasNew - biasOld
  biasEnergy = (betaNew - beta) * E_T - biasDiff

  ! Step 5: Apply Metropolis-like acceptance criterion
  ! Acceptance ratio: (T_new/T_old)^power * exp(-biasEnergy), where power = 3
  if ((TempNew / temperature)**power * exp(-biasEnergy) > grnd()) then
    acc_x = acc_x + 1.0_dp
    temperature = TempNew
    beta = betaNew
  endif
end subroutine TemperatureMove
!===========================================================================================
      end module
