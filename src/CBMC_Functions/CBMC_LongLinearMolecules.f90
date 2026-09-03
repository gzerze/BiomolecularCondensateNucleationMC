      module LLCBMC_Module
      contains
!=============================================================
! Performs a Configurational Bias Monte Carlo (CBMC) or Fixed Endpoint CBMC (FECBMC) move for long
! linear molecules (e.g., biomolecules) in a grand canonical ensemble nucleation simulation.
! Randomly selects a molecule and atom, chooses a regrowth direction, and decides between CBMC (for
! chain ends) or FECBMC (for internal sections). Regrows new and old conformations using
! LongChain_EndSections_ConfigGen or LongChain_InternalSections_ConfigGen, computing Rosenbluth
! weights (rosenProb_New, rosenProb_Old) and acceptance probability based on energy differences,
! umbrella sampling bias, and Jacobian factors (for FECBMC). Updates energies, positions, and
! neighbor lists if accepted. Supports minimum distance criterion and distance storage (useDistStore).
subroutine LLCBMC(E_T, acc_x, atmp_x)
  USE VarPrecision, ONLY: dp, atomIntType
  USE SimParameters, ONLY: maxAtoms, maxMol, NPART, &
       beta, calcPressure, distCriteria, minDistCriteria, pressure, prevMoveAccepted
  USE AcceptRates, ONLY: atmpCBMC, acptCBMC, atmpFECBMC, acptFECBMC, nCBMCmax, nFECBMCmax
  USE Coords, ONLY: MolArray, newMol
  USE CoordinateTypes, ONLY: NeighborDetails, Displacement
  USE CBMC_Variables, ONLY: pathArray, probTypeCBMC
  USE ForceField, ONLY: nAtoms
  USE EnergyPressurePointers, ONLY: Shift_ECalc, Shift_PCalc
  USE E_Interface_LJ_Q, ONLY: Update_SubEnergies
  USE EnergyTables, ONLY: ETable, E_NBond_Diff
  USE DistanceCriteria, ONLY: NeighborUpdate_Distance
  USE EnergyCriteria, ONLY: NeighborUpdate
  USE PairStorage, ONLY: UpdateDistArray, useDistStore
  USE UmbrellaSamplingNew, ONLY: useUmbrella, GetUmbrellaBias_Disp
  USE ParallelVar, ONLY: nout
  implicit none

  ! Input/Output: Total energy, acceptance counters
  real(dp), intent(inout) :: E_T                    ! Total system energy
  real(dp), intent(inout) :: acc_x                  ! Accepted move counter
  real(dp), intent(inout) :: atmp_x                 ! Attempted move counter

  ! Local variables
  logical, parameter :: useIntra(4) = [.true., .true., .false., .false.]  ! Intra-energy flags
  logical :: rejMove                                ! Rejection flag
  logical :: regrown(maxAtoms)                      ! Regrown status (false: regrow)
  logical :: regrowDirection                        ! True: positive direction, False: negative
  logical :: lCBMC, lFECBMC                         ! CBMC or FECBMC move flags
  integer :: iAtom, nDisp, nDirect, i               ! Loop index, displacement count, direction
  integer :: nType, nMol, nIndx                     ! Molecule type, index, global index
  integer :: nAtom                                  ! Selected atom for regrowth
  integer :: nGrow, nGrowStart, nGrowEnd            ! Number of segments, start/end path indices
  integer :: iGrowStart, iGrowEnd                   ! Start/end atom indices
  integer :: iPath                                  ! Path loop index
  real(dp) :: ranNum, sumInt                        ! Random number, cumulative sum
  real(dp) :: rosenProb_New, rosenProb_Old          ! New/old Rosenbluth weights
  real(dp) :: rosenRatio                            ! Ratio of old/new Rosenbluth weights
  real(dp) :: biasDiff                              ! Umbrella sampling bias difference
  real(dp) :: E_Inter, E_Intra                      ! Inter/intra-molecular energy differences
  real(dp) :: P_Diff                                ! Pressure difference
  real(dp) :: PairList(maxMol)                      ! Neighbor pair list
  real(dp) :: dETable(maxMol)                       ! Energy table differences
  type(displacement) :: disp(maxAtoms)              ! Displacements for energy calculation
  type(NeighborDetails) :: NeighborDetailsNew(1:maxMol)     ! Neighbor details (unused but declared)
  ! External function: Generates a uniform random number in (0,1) for trial selection
  real(dp), external :: grnd

  ! Initialize
  prevMoveAccepted = .false.
  rejMove = .false.
  atmp_x = atmp_x + 1.0_dp
  do i = 1, maxMol
    NeighborDetailsNew(i)%nPairs = 0
    NeighborDetailsNew(i)%pairIndices = 0
  enddo

  ! Select molecule type based on CBMC probabilities
  ranNum = grnd()
  sumInt = probTypeCBMC(1)
  nType = 1
  do while (sumInt < ranNum .and. nType < size(probTypeCBMC))
    nType = nType + 1
    sumInt = sumInt + probTypeCBMC(nType)
  enddo
  if (NPART(nType) == 0) return

  ! Select random molecule of chosen type
  nMol = floor(NPART(nType) * grnd() + 1.0_dp)
  nIndx = molArray(nType)%mol(nMol)%indx
  regrown = .true.

  if (pathArray(nType)%pathMax(1) < 3) then
    write(nout, *) "Error: Molecule type ", nType, " is too short for LLCBMC move, pathMax=", &
                pathArray(nType)%pathMax(1), " nCBMCmax=", nCBMCmax
    stop
  elseif (nCBMCmax > pathArray(nType)%pathMax(1)) then
  ! Select random atom for regrowth (2 <= nAtom <= nAtoms(nType)-1)
    nAtom = floor(grnd() * (pathArray(nType)%pathMax(1) - 2.0_dp)) + 2
  ! Choose regrowth direction
    if (grnd() > 0.5_dp) then
      regrowDirection = .true.
      nDirect = 1
    else
      regrowDirection = .false.
      nDirect = -1
    endif
  ! Find path index for selected atom
    nGrowStart = 0
    do iPath = 1, pathArray(nType)%pathMax(1)
      if (pathArray(nType)%path(1, iPath) == nAtom) then
        nGrowStart = iPath
        exit
      endif
    enddo
    ! CBMC move at chain end
    lCBMC = .true.
    nGrow = 0
    if (regrowDirection) then
      do iPath = nGrowStart + 1, pathArray(nType)%pathMax(1)
        iAtom = pathArray(nType)%path(1, iPath)
        regrown(iAtom) = .false.
        nGrow = nGrow + 1
      enddo
    else
      do iPath = nGrowStart - 1, 1, -1
        iAtom = pathArray(nType)%path(1, iPath)
        regrown(iAtom) = .false.
        nGrow = nGrow + 1
      enddo
    endif
    atmpCBMC(nType, nGrow) = atmpCBMC(nType, nGrow) + 1.0_dp

    ! Regrow new conformation
    call LongChain_EndSections_ConfigGen(nType, nMol, regrown, regrowDirection, nAtom, &
                                         rosenProb_New, rejMove)
    if (rejMove) return

    ! Regrow old conformation
    call LongChain_EndSections_ConfigGen_Reverse(nType, nMol, regrown, regrowDirection, &
                                                 nAtom, rosenProb_Old)
  else
  ! Select random atom for regrowth (3 <= nAtom <= nAtoms(nType)-3)
    nAtom = floor(grnd() * (pathArray(nType)%pathMax(1) - 5.0_dp)) + 3

  ! Choose regrowth direction
    if (grnd() > 0.5_dp) then
      regrowDirection = .true.
      nDirect = 1
    else
      regrowDirection = .false.
      nDirect = -1
    endif

  ! Find path index for selected atom
    nGrowStart = 0
    do iPath = 1, pathArray(nType)%pathMax(1)
      if (pathArray(nType)%path(1, iPath) == nAtom) then
        nGrowStart = iPath
        exit
      endif
    enddo

  ! Determine number of segments to regrow (2 <= nGrow <= nCBMCmax)
    nGrow = floor(grnd() * (nCBMCmax - 1.0_dp)) + 2
    nGrowEnd = nGrowStart + (nGrow + 1) * nDirect
  
  ! Decide between CBMC (chain end) or FECBMC (internal section)
    if (nGrowEnd > pathArray(nType)%pathMax(1) .or. nGrowEnd < 1) then
    ! CBMC move at chain end
      lCBMC = .true.
      nGrow = 0
      if (regrowDirection) then
        do iPath = nGrowStart + 1, pathArray(nType)%pathMax(1)
          iAtom = pathArray(nType)%path(1, iPath)
          regrown(iAtom) = .false.
          nGrow = nGrow + 1
        enddo
      else
        do iPath = nGrowStart - 1, 1, -1
          iAtom = pathArray(nType)%path(1, iPath)
          regrown(iAtom) = .false.
          nGrow = nGrow + 1
        enddo
      endif
      atmpCBMC(nType, nGrow) = atmpCBMC(nType, nGrow) + 1.0_dp

    ! Regrow new conformation
      call LongChain_EndSections_ConfigGen(nType, nMol, regrown, regrowDirection, nAtom, &
                                         rosenProb_New, rejMove)
      if (rejMove) return

    ! Regrow old conformation
      call LongChain_EndSections_ConfigGen_Reverse(nType, nMol, regrown, regrowDirection, &
                                                 nAtom, rosenProb_Old)
    else
    ! FECBMC move for internal section
      lCBMC = .false.
      if (nGrow .gt. nFECBMCmax) then
        lFECBMC = .false.
        do while (.not. lFECBMC)
          nGrow = floor(grnd() * (nFECBMCmax - 1.0_dp)) + 2
          nGrowEnd = nGrowStart + (nGrow + 1) * nDirect
          if (nGrowEnd <= pathArray(nType)%pathMax(1) .and. nGrowEnd >= 1) lFECBMC = .true.
        enddo
      endif
      atmpFECBMC(nType, nGrow) = atmpFECBMC(nType, nGrow) + 1.0_dp

    ! Mark intermediate atoms for regrowth
      do iPath = nGrowStart + nDirect, nGrowEnd - nDirect, nDirect
        iAtom = pathArray(nType)%path(1, iPath)
        regrown(iAtom) = .false.
      enddo
      iGrowStart = pathArray(nType)%path(1, nGrowStart)
      iGrowEnd = pathArray(nType)%path(1, nGrowEnd)

    ! Regrow new conformation
      call LongChain_InternalSections_ConfigGen(nType, nMol, regrown, iGrowStart, iGrowEnd, &
                                              rosenProb_New, rejMove)
      if (rejMove) return

    ! Regrow old conformation
      call LongChain_InternalSections_ConfigGen_Reverse(nType, nMol, regrown, iGrowStart, &
                                                      iGrowEnd, rosenProb_Old)
    endif
  endif

  ! Prepare displacements for energy calculation
  nDisp = 0
  do iAtom = 1, nAtoms(nType)
    nDisp = nDisp + 1
    disp(nDisp)%molType = int(nType, atomIntType)
    disp(nDisp)%molIndx = nMol
    disp(nDisp)%atmIndx = iAtom
    disp(nDisp)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
    disp(nDisp)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
    disp(nDisp)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
    if (.not. regrown(iAtom)) then
      disp(nDisp)%Displaced = .true.
      disp(nDisp)%x_new = newMol%x(iAtom)
      disp(nDisp)%y_new = newMol%y(iAtom)
      disp(nDisp)%z_new = newMol%z(iAtom)
    else
      disp(nDisp)%Displaced = .false.
      disp(nDisp)%x_new = disp(nDisp)%x_old
      disp(nDisp)%y_new = disp(nDisp)%y_old
      disp(nDisp)%z_new = disp(nDisp)%z_old
    endif
  enddo

  ! Compute energy differences
  E_Inter = 0.0_dp
  E_Intra = 0.0_dp
  call Shift_ECalc(E_Inter, E_Intra, disp(1:nDisp), PairList, dETable, useIntra, &
                   rejMove, .true., NeighborDetailsNew)
  if (rejMove) return

  ! Compute acceptance probability
  rosenRatio = rosenProb_Old / rosenProb_New
  biasDiff = 0.0_dp
  if (useUmbrella) then
    call GetUmbrellaBias_Disp(disp(1:nDisp), biasDiff, rejMove)
    if (rejMove) return
  endif

  if (log(rosenRatio) - beta * (E_Inter + E_NBond_Diff) + biasDiff > log(grnd())) then
    ! Update system state
    if (calcPressure) then
      call Shift_PCalc(P_Diff, disp(1:nDisp))
      pressure = pressure + P_Diff
    endif
    do iAtom = 1, nDisp
      disp(iAtom)%x_old = disp(iAtom)%x_new
      disp(iAtom)%y_old = disp(iAtom)%y_new
      disp(iAtom)%z_old = disp(iAtom)%z_new
    enddo
    E_T = E_T + E_Inter + E_Intra
    ETable = ETable + dETable
    acc_x = acc_x + 1.0_dp
    if (lCBMC) then
      acptCBMC(nType, nGrow) = acptCBMC(nType, nGrow) + 1.0_dp
    else
      acptFECBMC(nType, nGrow) = acptFECBMC(nType, nGrow) + 1.0_dp
    endif
    prevMoveAccepted = .true.

    ! Update distance storage and neighbor lists
    if (useDistStore) call UpdateDistArray
    if (distCriteria .or. minDistCriteria) then
      call NeighborUpdate_Distance(PairList, nIndx, NeighborDetailsNew)
    else
      call NeighborUpdate(PairList, nIndx)
    endif    
    call Update_SubEnergies
  endif

end subroutine LLCBMC
!============================================================
! Performs Configurational Bias Monte Carlo (CBMC) regrowth of a chain end for long linear
! biomolecules in a grand canonical ensemble nucleation simulation. Starting from a specified
! atom, regrows atoms in the specified direction (positive or negative) using the Rosenbluth
! scheme. Generates trial positions with a single bond type, flexible bending (sine distribution),
! and uniform azimuthal angles. Computes intermolecular and intramolecular non-bonded energies
! for trial selection. Updates the new molecule configuration (`newMol`) and calculates the
! Rosenbluth weight (`rosenRatio`). Rejects moves if trial positions overlap or are invalid.
subroutine LongChain_EndSections_ConfigGen(nType, nMol, regrownIn, regrowDirection, startAtom, rosenRatio, rejMove)
  USE VarPrecision, ONLY: dp, atomIntType
  USE SimParameters, ONLY: maxAtoms, beta
  USE Coords, ONLY: MolArray, newMol
  USE CoordinateTypes, ONLY: SimpleAtomCoords
  USE CoordinateFunctions, ONLY: GenerateGaussianBondLength, Generate_UnitCone
  USE ForceField, ONLY: bondData, bondArray, atomArray, nonBondArray, nIntraNonBond, nAtoms
  USE EnergyPressurePointers, ONLY: Find_InterAtms, RosenLL_Atom_New, RosenLL_Atom_Intra_New
  USE CBMC_Variables, ONLY: pathArray, nRosenTrials, maxRosenTrial
  USE PairStorage, ONLY: nTotalAtoms
  USE ParallelVar, ONLY: nout
  implicit none

  ! Input/Output variables
  integer, intent(in) :: nType                      ! Molecule type index
  integer, intent(in) :: nMol                       ! Molecule instance index
  logical, intent(in) :: regrownIn(maxAtoms)        ! Input regrown status of atoms
  logical, intent(in) :: regrowDirection            ! True: positive direction, False: negative
  integer, intent(in) :: startAtom                  ! Starting atom index for regrowth
  real(dp), intent(out) :: rosenRatio               ! Rosenbluth weight for new conformation
  logical, intent(out) :: rejMove                   ! Flag to reject move if invalid

  ! Local variables
  integer :: i, iAtom, jAtom, nIndx, iRosen                ! Loop indices, molecule global index
  integer :: Atm2, Atm3, Atm4                      ! Atom indices for chain geometry
  integer :: bondType                              ! Bond type index (single for molecule type)
  integer :: curPos                                ! Current path position
  integer :: Incrmt                                ! Increment direction (+1 or -1)
  integer :: totalRegrown                          ! Count of regrown atoms
  integer :: nIntraUnits                           ! Number of intramolecular non-bonded units
  integer :: IntraUnits(maxAtoms)                  ! Indices of intramolecular non-bonded atoms
  integer :: nInterAtoms                           ! Number of intermolecular interacting atoms
  integer :: InterAtoms(nTotalAtoms, 3)            ! Interacting atom indices
  integer(kind=atomIntType) :: atomType            ! Atom type for current atom
  logical :: regrown(maxAtoms)                     ! Regrown status of atoms
  logical :: overlap(maxRosenTrial)                ! Overlap flags for trial positions
  real(dp) :: k_bond, r_eq, r_sigma, rmax_sq       ! Bond parameters
  real(dp) :: r, r0, Prob, ProbRosen(1:maxRosenTrial)  ! Bond length and probability
  real(dp) :: bend_angle                           ! Bending angle for trial position
  real(dp) :: E_Trial(maxRosenTrial)              ! Trial energies
  real(dp) :: E_Trial_Intra                  ! Intramolecular trial energy
  real(dp) :: E_Min, rosenNorm, ranNum, sumInt     ! Rosenbluth selection variables
  integer :: nSel                                  ! Selected trial index
  type(SimpleAtomCoords) :: trialPos(maxRosenTrial) ! Trial positions
  type(SimpleAtomCoords) :: RefAtom, v2, v3        ! Reference atom and vectors for geometry
  ! External function: Generates a uniform random number in (0,1) for trial selection
  real(dp), external :: grnd

  ! Initialize variables
  rejMove = .false.
  rosenRatio = 1.0_dp
  regrown = regrownIn
  totalRegrown = 0
  nIndx = MolArray(nType)%mol(nMol)%indx

  ! Initialize new molecule configuration
  newMol%molType = nType
  newMol%x = 0.0_dp
  newMol%y = 0.0_dp
  newMol%z = 0.0_dp

  ! Copy coordinates of non-regrown atoms
  do i = 1, nAtoms(nType)
    if (regrown(i)) then
      newMol%x(i) = MolArray(nType)%mol(nMol)%x(i)
      newMol%y(i) = MolArray(nType)%mol(nMol)%y(i)
      newMol%z(i) = MolArray(nType)%mol(nMol)%z(i)
      totalRegrown = totalRegrown + 1
    endif
  enddo

  ! Set bond parameters (single bond type for molecule type)
  bondType = bondArray(nType, 1)%bondType
  k_bond = bondData(bondType)%k_eq
  r_eq = bondData(bondType)%r_eq
  r_sigma = bondData(bondType)%r_sigma
  rmax_sq = bondData(bondType)%rmax_sq

  ! Set regrowth direction increment
  Incrmt = merge(1, -1, regrowDirection)

  ! Find starting path position
  curPos = 0
  do i = 1, pathArray(nType)%pathMax(1)
    if (pathArray(nType)%path(1, i) == startAtom) then
      curPos = i
      exit
    endif
  enddo
  if (curPos == 0) then
    write(nout, *) "Error: Start atom ", startAtom, " not found in path for molecule type ", nType
    stop
  endif
  curPos = curPos + Incrmt

  ! Regrow atoms until all are regrown
  do while (any(.not. regrown))
    ! Exit if current position is out of bounds
    if (curPos < 1 .or. curPos > pathArray(nType)%pathMax(1)) exit

    ! Get atom indices for chain geometry
    Atm4 = pathArray(nType)%path(1, curPos)            ! Current atom to regrow
    Atm3 = pathArray(nType)%path(1, curPos - Incrmt)   ! Reference atom
    Atm2 = pathArray(nType)%path(1, curPos - 2*Incrmt) ! Second reference atom for bending

    ! Find intramolecular non-bonded atoms already regrown
    nIntraUnits = 0
    do i = 1, nIntraNonBond(nType)
      iAtom = nonBondArray(nType, i)%nonMembr(1)
      jAtom = nonBondArray(nType, i)%nonMembr(2)
      ! Ensure Atm4 is one of the pair
      if (iAtom /= Atm4 .and. jAtom /= Atm4) cycle
      ! Swap indices if needed to make iAtom = Atm4
      if (iAtom /= Atm4) then
        jAtom = iAtom
        iAtom = Atm4
      endif
      ! Include only if jAtom is already regrown
      if (.not. regrown(jAtom)) cycle
      nIntraUnits = nIntraUnits + 1
      IntraUnits(nIntraUnits) = jAtom
    enddo

    ! Get atom type and reference atom coordinates
    atomType = atomArray(nType, Atm4)
    RefAtom%x = newMol%x(Atm3)
    RefAtom%y = newMol%y(Atm3)
    RefAtom%z = newMol%z(Atm3)

    ! Find intermolecular interacting atoms within cutoff + buffer
    r0 = r_eq
    call Find_InterAtms(RefAtom, r0, atomType, nIndx, nInterAtoms, InterAtoms)

    ! Generate trial positions
    v2%x = newMol%x(Atm2) - newMol%x(Atm3)
    v2%y = newMol%y(Atm2) - newMol%y(Atm3)
    v2%z = newMol%z(Atm2) - newMol%z(Atm3)
    overlap = .false.
    E_Trial = 0.0_dp

    do iRosen = 1, nRosenTrials(nType)
      ! Generate bond length with Gaussian distribution
      call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)

      ! Generate bending angle from sine distribution: θ = acos(1 - 2*rand)
      bend_angle = acos(1.0_dp - 2.0_dp * grnd())

      ! Generate trial position using unit cone method
      call Generate_UnitCone(v2, r, bend_angle, v3)

      ! Compute trial coordinates
      trialPos(iRosen)%x = v3%x + newMol%x(Atm3)
      trialPos(iRosen)%y = v3%y + newMol%y(Atm3)
      trialPos(iRosen)%z = v3%z + newMol%z(Atm3)

      ! Calculate intermolecular non-bonded energy
      call RosenLL_Atom_New(nType, Atm4, trialPos(iRosen), nInterAtoms, InterAtoms, &
                            E_Trial(iRosen), overlap(iRosen))

      ! Calculate intramolecular non-bonded energy if no overlap
      if (.not. overlap(iRosen)) then
        call RosenLL_Atom_Intra_New(nType, Atm4, trialPos(iRosen), nIntraUnits, IntraUnits, &
                                    E_Trial_Intra, overlap(iRosen))
        E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_Intra
      endif
    enddo

    ! Inline ChooseRosenTrial: Select trial position based on Rosenbluth weights
    E_Min = minval(E_Trial)
    ProbRosen = 0.0_dp
    do iRosen = 1, nRosenTrials(nType)
      ProbRosen(iRosen) = exp(-beta * (E_Trial(iRosen) - E_Min))
    enddo

    ! Check for valid trials
    if (all(ProbRosen <= 0.0_dp)) then
      rejMove = .true.
      return
    endif

    ! Normalize and select trial
    rosenNorm = sum(ProbRosen)
    ranNum = grnd() * rosenNorm
    sumInt = ProbRosen(1)
    nSel = 1
    do while (sumInt < ranNum .and. nSel < nRosenTrials(nType))
      nSel = nSel + 1
      sumInt = sumInt + ProbRosen(nSel)
    enddo

    ! Check for overlap in selected trial
    if (overlap(nSel)) then
      rejMove = .true.
!  write(*,*) "CBMC Rejected"
      return
    endif

    ! Update Rosenbluth weight and store selected position
    rosenRatio = rosenRatio * ProbRosen(nSel) / rosenNorm
    regrown(Atm4) = .true.
    newMol%x(Atm4) = trialPos(nSel)%x
    newMol%y(Atm4) = trialPos(nSel)%y
    newMol%z(Atm4) = trialPos(nSel)%z

    ! Move to next atom
    curPos = curPos + Incrmt
    totalRegrown = totalRegrown + 1
  enddo

end subroutine LongChain_EndSections_ConfigGen
!============================================================
! Performs reverse Configurational Bias Monte Carlo (CBMC) regrowth of a chain end for long linear
! biomolecules in a grand canonical ensemble nucleation simulation. Starting from a specified atom,
! regrows atoms in the specified direction (positive or negative) to compute the Rosenbluth weight
! of the existing (old) conformation for detailed balance. Uses the current atom positions as the
! first trial and generates additional trials using a single bond type, flexible bending (sine
! distribution), and uniform azimuthal angles. Computes intermolecular and intramolecular non-bonded
! energies for trial selection. Returns the Rosenbluth weight (`rosenRatio`) without updating
! coordinates. Rejects moves if trial generation fails.
subroutine LongChain_EndSections_ConfigGen_Reverse(nType, nMol, regrownIn, regrowDirection, startAtom, rosenRatio)
  USE VarPrecision, ONLY: dp, atomIntType
  USE SimParameters, ONLY: maxAtoms, beta
  USE Coords, ONLY: MolArray
  USE CoordinateTypes, ONLY: SimpleAtomCoords
  USE CoordinateFunctions, ONLY: GenerateGaussianBondLength, Generate_UnitCone
  USE ForceField, ONLY: bondData, bondArray, atomArray, nonBondArray, nIntraNonBond, nAtoms
  USE EnergyPressurePointers, ONLY: Find_InterAtms, RosenLL_Atom_New, RosenLL_Atom_Old, RosenLL_Atom_Intra_Old
  USE CBMC_Variables, ONLY: pathArray, nRosenTrials, maxRosenTrial
  USE PairStorage, ONLY: nTotalAtoms
  USE ParallelVar, ONLY: nout
  implicit none

  ! Input/Output variables
  integer, intent(in) :: nType                      ! Molecule type index
  integer, intent(in) :: nMol                       ! Molecule instance index
  logical, intent(in) :: regrownIn(maxAtoms)        ! Input regrown status of atoms
  logical, intent(in) :: regrowDirection            ! True: positive direction, False: negative
  integer, intent(in) :: startAtom                  ! Starting atom index for regrowth
  real(dp), intent(out) :: rosenRatio              ! Rosenbluth weight for old conformation

  ! Local variables
  integer :: i, iRosen, iAtom, jAtom               ! Loop indices
  integer :: Atm2, Atm3, Atm4                      ! Atom indices for chain geometry
  integer :: bondType                              ! Bond type index (single for molecule type)
  integer :: curPos                                ! Current path position
  integer :: Incrmt                                ! Increment direction (+1 or -1)
  integer :: totalRegrown                          ! Count of regrown atoms
  integer :: nIntraUnits                           ! Number of intramolecular non-bonded units
  integer :: IntraUnits(maxAtoms)                  ! Indices of intramolecular non-bonded atoms
  integer :: nInterAtoms                           ! Number of intermolecular interacting atoms
  integer :: InterAtoms(nTotalAtoms, 3)            ! Interacting atom indices
  integer(kind=atomIntType) :: atomType            ! Atom type for current atom
  logical :: regrown(maxAtoms)                     ! Regrown status of atoms
  logical :: overlap(maxRosenTrial)                ! Overlap flags for trial positions
  real(dp) :: k_bond, r_eq, r_sigma, rmax_sq       ! Bond parameters
  real(dp) :: r, Prob, r0                              ! Bond length and probability
  real(dp) :: bend_angle                           ! Bending angle for trial position
  real(dp) :: E_Trial(maxRosenTrial)              ! Trial energies
  real(dp) :: E_Trial_Intra                  ! Intramolecular trial energy
  real(dp) :: E_Min, rosenNorm, ProbRosen(maxRosenTrial) ! Rosenbluth selection variables
  type(SimpleAtomCoords) :: trialPos(maxRosenTrial) ! Trial positions
  type(SimpleAtomCoords) :: RefAtom, v2, v3        ! Reference atom and vectors for geometry
  integer :: nIndx                                 ! Molecule global index
  ! External function: Generates a uniform random number in (0,1) for trial selection
  real(dp), external :: grnd

  ! Initialize variables
  rosenRatio = 1.0_dp
  regrown = regrownIn
  totalRegrown = 0
  nIndx = MolArray(nType)%mol(nMol)%indx

  ! Count initially regrown atoms
  do i = 1, nAtoms(nType)
    if (regrown(i)) totalRegrown = totalRegrown + 1
  enddo

  ! Set bond parameters (single bond type for molecule type)
  bondType = bondArray(nType, 1)%bondType
  k_bond = bondData(bondType)%k_eq
  r_eq = bondData(bondType)%r_eq
  r_sigma = bondData(bondType)%r_sigma
  rmax_sq = bondData(bondType)%rmax_sq

  ! Set regrowth direction increment
  Incrmt = merge(1, -1, regrowDirection)

  ! Find starting path position
  curPos = 0
  do i = 1, pathArray(nType)%pathMax(1)
    if (pathArray(nType)%path(1, i) == startAtom) then
      curPos = i
      exit
    endif
  enddo
  if (curPos == 0) then
    write(nout, *) "Error: Start atom ", startAtom, " not found in path for molecule type ", nType
    stop
  endif
  curPos = curPos + Incrmt

  ! Regrow atoms to compute Rosenbluth weight
  do while (any(.not. regrown))
    ! Exit if current position is out of bounds
    if (curPos < 1 .or. curPos > pathArray(nType)%pathMax(1)) exit

    ! Get atom indices for chain geometry
    Atm4 = pathArray(nType)%path(1, curPos)            ! Current atom to regrow
    Atm3 = pathArray(nType)%path(1, curPos - Incrmt)   ! Reference atom
    Atm2 = pathArray(nType)%path(1, curPos - 2*Incrmt) ! Second reference atom for bending

    ! Find intramolecular non-bonded atoms already regrown
    nIntraUnits = 0
    do i = 1, nIntraNonBond(nType)
      iAtom = nonBondArray(nType, i)%nonMembr(1)
      jAtom = nonBondArray(nType, i)%nonMembr(2)
      ! Ensure Atm4 is one of the pair
      if (iAtom /= Atm4 .and. jAtom /= Atm4) cycle
      ! Swap indices if needed to make iAtom = Atm4
      if (iAtom /= Atm4) then
        jAtom = iAtom
        iAtom = Atm4
      endif
      ! Include only if jAtom is already regrown
      if (.not. regrown(jAtom)) cycle
      nIntraUnits = nIntraUnits + 1
      IntraUnits(nIntraUnits) = jAtom
    enddo

    ! Get atom type and reference atom coordinates
    atomType = atomArray(nType, Atm4)
    RefAtom%x = MolArray(nType)%mol(nMol)%x(Atm3)
    RefAtom%y = MolArray(nType)%mol(nMol)%y(Atm3)
    RefAtom%z = MolArray(nType)%mol(nMol)%z(Atm3)

    ! Find intermolecular interacting atoms within cutoff + buffer
    r0 = r_eq
    call Find_InterAtms(RefAtom, r0, atomType, nIndx, nInterAtoms, InterAtoms)

    ! Generate trial positions
    v2%x = MolArray(nType)%mol(nMol)%x(Atm2) - MolArray(nType)%mol(nMol)%x(Atm3)
    v2%y = MolArray(nType)%mol(nMol)%y(Atm2) - MolArray(nType)%mol(nMol)%y(Atm3)
    v2%z = MolArray(nType)%mol(nMol)%z(Atm2) - MolArray(nType)%mol(nMol)%z(Atm3)
    overlap = .false.
    E_Trial = 0.0_dp

    ! First trial: existing atom position
    iRosen = 1
    trialPos(iRosen)%x = MolArray(nType)%mol(nMol)%x(Atm4)
    trialPos(iRosen)%y = MolArray(nType)%mol(nMol)%y(Atm4)
    trialPos(iRosen)%z = MolArray(nType)%mol(nMol)%z(Atm4)
    call RosenLL_Atom_Old(nType, nMol, Atm4, nInterAtoms, InterAtoms, E_Trial(iRosen))
    call RosenLL_Atom_Intra_Old(nType, nMol, Atm4, trialPos(iRosen), nIntraUnits, IntraUnits, E_Trial_Intra)
    E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_Intra

    ! Generate additional trials
    do iRosen = 2, nRosenTrials(nType)
      ! Generate bond length with Gaussian distribution
      call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)

      ! Generate bending angle from sine distribution: θ = acos(1 - 2*rand)
      bend_angle = acos(1.0_dp - 2.0_dp * grnd())

      ! Generate trial position using unit cone method
      call Generate_UnitCone(v2, r, bend_angle, v3)

      ! Compute trial coordinates
      trialPos(iRosen)%x = v3%x + MolArray(nType)%mol(nMol)%x(Atm3)
      trialPos(iRosen)%y = v3%y + MolArray(nType)%mol(nMol)%y(Atm3)
      trialPos(iRosen)%z = v3%z + MolArray(nType)%mol(nMol)%z(Atm3)

      ! Calculate intermolecular non-bonded energy
      call RosenLL_Atom_New(nType, Atm4, trialPos(iRosen), nInterAtoms, InterAtoms, &
                            E_Trial(iRosen), overlap(iRosen))

      ! Calculate intramolecular non-bonded energy if no overlap
      if (.not. overlap(iRosen)) then
        call RosenLL_Atom_Intra_Old(nType, nMol, Atm4, trialPos(iRosen), nIntraUnits, IntraUnits, &
                                    E_Trial_Intra)
        E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_Intra
      endif
    enddo

    ! Inline ChooseRosenTrial: Compute Rosenbluth weight for existing position
    E_Min = minval(E_Trial)
    ProbRosen = 0.0_dp
    do iRosen = 1, nRosenTrials(nType)
      ProbRosen(iRosen) = exp(-beta * (E_Trial(iRosen) - E_Min))
    enddo

    ! Check for valid trials
    if (all(ProbRosen <= 0.0_dp)) then
      write(nout, *) "Error: No valid trials in reverse regrowth for molecule ", nIndx, " of type ", nType
      stop
    endif

    ! Compute Rosenbluth weight (select first trial, existing position)
    rosenNorm = sum(ProbRosen)
    rosenRatio = rosenRatio * ProbRosen(1) / rosenNorm

    ! Mark atom as regrown and move to next
    regrown(Atm4) = .true.
    curPos = curPos + Incrmt
    totalRegrown = totalRegrown + 1
  enddo

end subroutine LongChain_EndSections_ConfigGen_Reverse
!============================================================
! Performs Fixed Endpoint Configurational Bias Monte Carlo (FECBMC) regrowth of an internal chain
! section for long linear biomolecules in a grand canonical ensemble nucleation simulation.
! Regrows atoms between fixed start and end atoms, ensuring geometric connectivity to the endpoint.
! Uses a single bond type, flexible bending (sine distribution for intermediate atoms, fixed for
! the final atom), and uniform azimuthal angles. Computes intermolecular and intramolecular
! non-bonded energies for trial selection, applies a Jacobian factor for constrained growth, and
! returns the Rosenbluth weight (`rosenRatio`). Updates the new molecule configuration (`newMol`).
! Rejects moves if trial positions overlap or geometric constraints are violated.
! For details of the method look at
! Escobedo, F. A.; De Pablo, J. J. Extended Continuum Configurational Bias Monte Carlo Methods for Simulation of Flexible Molecules. 
! J Chem Phys 1995, 102 (6), 2636–2652. https://doi.org/10.1063/1.468695.
subroutine LongChain_InternalSections_ConfigGen(nType, nMol, regrownIn, startAtom, endAtom, rosenRatio, rejMove)
  USE VarPrecision, ONLY: dp, atomIntType
  USE SimParameters, ONLY: maxAtoms, beta
  USE Coords, ONLY: MolArray, newMol
  USE CoordinateTypes, ONLY: SimpleAtomCoords
  USE CoordinateFunctions, ONLY: GenerateGaussianBondLength, Generate_UnitCone
  USE ForceField, ONLY: bondData, bondArray, atomArray, nonBondArray, nIntraNonBond, nAtoms
  USE EnergyPressurePointers, ONLY: Find_InterAtms, RosenLL_Atom_New, RosenLL_Atom_Intra_New
  USE CBMC_Variables, ONLY: pathArray, nRosenTrials, maxRosenTrial
  USE PairStorage, ONLY: nTotalAtoms
  USE ParallelVar, ONLY: nout
  implicit none

  ! Input/Output variables
  integer, intent(in) :: nType                      ! Molecule type index
  integer, intent(in) :: nMol                       ! Molecule instance index
  logical, intent(in) :: regrownIn(maxAtoms)        ! Input regrown status of atoms
  integer, intent(in) :: startAtom, endAtom         ! Fixed start and end atom indices
  real(dp), intent(out) :: rosenRatio              ! Rosenbluth weight for new conformation
  logical, intent(out) :: rejMove                   ! Flag to reject move if invalid

  ! Local variables
  integer :: i, iRosen, iAtom, jAtom               ! Loop indices
  integer :: Atm2, Atm3, Atm4                      ! Atom indices for chain geometry
  integer :: bondType                              ! Bond type index (single for molecule type)
  integer :: curPos, lastPos                       ! Current and last path positions
  integer :: Incrmt                                ! Increment direction (+1 or -1)
  integer :: totalRegrown                          ! Count of regrown atoms
  integer :: nIntraUnits                           ! Number of intramolecular non-bonded units
  integer :: IntraUnits(maxAtoms)                  ! Indices of intramolecular non-bonded atoms
  integer :: nInterAtoms                           ! Number of intermolecular interacting atoms
  integer :: InterAtoms(nTotalAtoms, 3)            ! Interacting atom indices
  integer :: nGrowUnits, nGrowBonds                ! Number of atoms and bonds to regrow
  integer :: iGrowUnits                            ! Current regrowth unit index
  integer(kind=atomIntType) :: atomType            ! Atom type for current atom
  logical :: regrown(maxAtoms)                     ! Regrown status of atoms
  logical :: overlap(maxRosenTrial)                ! Overlap flags for trial positions
  real(dp) :: k_bond, r_eq, r_sigma, rmax_sq       ! Bond parameters
  real(dp) :: r, Prob                              ! Bond length and probability
  real(dp) :: bend_angle, costhmax                 ! Bending angle and max cosine for constraint
  real(dp) :: r0End, r1End                         ! Distances to endpoint
  real(dp) :: r0                                   ! Buffer distance for intermolecular interactions
  real(dp) :: E_Trial(maxRosenTrial)              ! Trial energies
  real(dp) :: E_Trial_Intra                        ! Intramolecular trial energy
  real(dp) :: ranNum, sumInt, jacobian       ! Jacobian factor for constrained growth
  real(dp) :: GrowBonds(maxAtoms)                  ! Pre-generated bond lengths
  real(dp) :: E_Min, rosenNorm, ProbRosen(maxRosenTrial) ! Rosenbluth selection variables
  integer :: nSel                                  ! Selected trial index
  type(SimpleAtomCoords) :: trialPos(maxRosenTrial) ! Trial positions
  type(SimpleAtomCoords) :: RefAtom, v2, v3        ! Reference atom and vectors for geometry
  integer :: nIndx                                 ! Molecule global index
  ! External function: Generates a uniform random number in (0,1) for trial selection
  real(dp), external :: grnd

  ! Initialize variables
  rejMove = .false.
  rosenRatio = 1.0_dp
  jacobian = 1.0_dp
  regrown = regrownIn
  totalRegrown = 0
  nIndx = MolArray(nType)%mol(nMol)%indx

  ! Initialize new molecule configuration
  newMol%molType = nType
  newMol%x = 0.0_dp
  newMol%y = 0.0_dp
  newMol%z = 0.0_dp

  ! Copy coordinates of non-regrown atoms
  do i = 1, nAtoms(nType)
    if (regrown(i)) then
      newMol%x(i) = MolArray(nType)%mol(nMol)%x(i)
      newMol%y(i) = MolArray(nType)%mol(nMol)%y(i)
      newMol%z(i) = MolArray(nType)%mol(nMol)%z(i)
      totalRegrown = totalRegrown + 1
    endif
  enddo

  ! Set bond parameters (single bond type for molecule type)
  bondType = bondArray(nType, 1)%bondType
  k_bond = bondData(bondType)%k_eq
  r_eq = bondData(bondType)%r_eq
  r_sigma = bondData(bondType)%r_sigma
  rmax_sq = bondData(bondType)%rmax_sq

  ! Find path positions for start and end atoms
  curPos = 0
  lastPos = 0
  do i = 1, pathArray(nType)%pathMax(1)
    if (pathArray(nType)%path(1, i) == startAtom) curPos = i
    if (pathArray(nType)%path(1, i) == endAtom) lastPos = i
  enddo
  if (curPos == 0 .or. lastPos == 0) then
    write(nout, *) "Error: Start atom ", startAtom, " or end atom ", endAtom, &
                " not found in path for molecule type ", nType
    stop
  endif

  ! Determine regrowth direction and number of bonds/units
  if (lastPos > curPos) then
    Incrmt = 1
    nGrowBonds = lastPos - curPos
  else
    Incrmt = -1
    nGrowBonds = curPos - lastPos
  endif
  nGrowUnits = nGrowBonds - 1

  ! Generate bond lengths for all regrowth bonds
  do i = 1, nGrowBonds
    call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
    GrowBonds(i) = r
  enddo

  ! Regrow intermediate atoms
  curPos = curPos + Incrmt
  do iGrowUnits = 1, nGrowUnits - 1
    ! Exit if current position is out of bounds
    if (curPos < 1 .or. curPos > pathArray(nType)%pathMax(1)) exit

    ! Get atom indices for chain geometry
    Atm4 = pathArray(nType)%path(1, curPos)            ! Current atom to regrow
    Atm3 = pathArray(nType)%path(1, curPos - Incrmt)   ! Reference atom
    Atm2 = pathArray(nType)%path(1, curPos - 2*Incrmt) ! Second reference atom for bending

    ! Find intramolecular non-bonded atoms already regrown
    nIntraUnits = 0
    do i = 1, nIntraNonBond(nType)
      iAtom = nonBondArray(nType, i)%nonMembr(1)
      jAtom = nonBondArray(nType, i)%nonMembr(2)
      if (iAtom /= Atm4 .and. jAtom /= Atm4) cycle
      if (iAtom /= Atm4) then
        jAtom = iAtom
        iAtom = Atm4
      endif
      if (.not. regrown(jAtom)) cycle
      nIntraUnits = nIntraUnits + 1
      IntraUnits(nIntraUnits) = jAtom
    enddo

    ! Get atom type and reference atom coordinates
    atomType = atomArray(nType, Atm4)
    RefAtom%x = newMol%x(Atm3)
    RefAtom%y = newMol%y(Atm3)
    RefAtom%z = newMol%z(Atm3)

    ! Find intermolecular interacting atoms within cutoff + buffer
    r0 = r_eq
    call Find_InterAtms(RefAtom, r0, atomType, nIndx, nInterAtoms, InterAtoms)

    ! Compute geometric constraints
    r = GrowBonds(iGrowUnits)
    r1End = sum(GrowBonds(iGrowUnits + 1:nGrowUnits + 1))
    v2%x = newMol%x(endAtom) - newMol%x(Atm3)
    v2%y = newMol%y(endAtom) - newMol%y(Atm3)
    v2%z = newMol%z(endAtom) - newMol%z(Atm3)
    r0End = sqrt(v2%x * v2%x + v2%y * v2%y + v2%z * v2%z)

    ! Check triangle inequality for geometric feasibility
    if (r + r1End < r0End) then
      rejMove = .true.
      return
    endif

    ! Compute maximum cosine for constrained angle
    if (r1End > r + r0End) then
      costhmax = -1.0_dp
    else
      costhmax = (r * r + r0End * r0End - r1End * r1End) / (2.0_dp * r * r0End)
      if (costhmax <= -1.0_dp) costhmax = -1.0_dp
    endif
    jacobian = jacobian * (1.0_dp - costhmax)

    ! Generate trial positions
    overlap = .false.
    E_Trial = 0.0_dp
    do iRosen = 1, nRosenTrials(nType)
      ! Generate bending angle within constrained range
      bend_angle = acos(1.0_dp - (1.0_dp - costhmax) * grnd())
      call Generate_UnitCone(v2, r, bend_angle, v3)
      trialPos(iRosen)%x = v3%x + newMol%x(Atm3)
      trialPos(iRosen)%y = v3%y + newMol%y(Atm3)
      trialPos(iRosen)%z = v3%z + newMol%z(Atm3)

      ! Calculate intermolecular non-bonded energy
      call RosenLL_Atom_New(nType, Atm4, trialPos(iRosen), nInterAtoms, InterAtoms, &
                            E_Trial(iRosen), overlap(iRosen))

      ! Calculate intramolecular non-bonded energy if no overlap
      if (.not. overlap(iRosen)) then
        call RosenLL_Atom_Intra_New(nType, Atm4, trialPos(iRosen), nIntraUnits, IntraUnits, &
                                    E_Trial_Intra, overlap(iRosen))
        E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_Intra
      endif
    enddo

    ! Inline ChooseRosenTrial: Select trial position based on Rosenbluth weights
    E_Min = minval(E_Trial)
    ProbRosen = 0.0_dp
    do iRosen = 1, nRosenTrials(nType)
      ProbRosen(iRosen) = exp(-beta * (E_Trial(iRosen) - E_Min))
    enddo

    ! Check for valid trials
    if (all(ProbRosen <= 0.0_dp)) then
      rejMove = .true.
      return
    endif

    ! Normalize and select trial
    rosenNorm = sum(ProbRosen)
    ranNum = grnd() * rosenNorm
    sumInt = ProbRosen(1)
    nSel = 1
    do while (sumInt < ranNum .and. nSel < nRosenTrials(nType))
      nSel = nSel + 1
      sumInt = sumInt + ProbRosen(nSel)
    enddo

    ! Check for overlap in selected trial
    if (overlap(nSel)) then
      rejMove = .true.
      return
    endif

    ! Update Rosenbluth weight and store selected position
    rosenRatio = rosenRatio * ProbRosen(nSel) / rosenNorm
    regrown(Atm4) = .true.
    newMol%x(Atm4) = trialPos(nSel)%x
    newMol%y(Atm4) = trialPos(nSel)%y
    newMol%z(Atm4) = trialPos(nSel)%z

    ! Move to next atom
    curPos = curPos + Incrmt
    totalRegrown = totalRegrown + 1
  enddo

  ! Regrow final atom with fixed angle to reach endpoint
  if (curPos < 1 .or. curPos > pathArray(nType)%pathMax(1)) then
    write(nout, *) "Error: Final atom position out of bounds for molecule ", nIndx, " of type ", nType
    stop
  endif

  Atm4 = pathArray(nType)%path(1, curPos)
  Atm3 = pathArray(nType)%path(1, curPos - Incrmt)
  Atm2 = pathArray(nType)%path(1, curPos - 2*Incrmt)

  ! Find intramolecular non-bonded atoms
  nIntraUnits = 0
  do i = 1, nIntraNonBond(nType)
    iAtom = nonBondArray(nType, i)%nonMembr(1)
    jAtom = nonBondArray(nType, i)%nonMembr(2)
    if (iAtom /= Atm4 .and. jAtom /= Atm4) cycle
    if (iAtom /= Atm4) then
      jAtom = iAtom
      iAtom = Atm4
    endif
    if (.not. regrown(jAtom)) cycle
    nIntraUnits = nIntraUnits + 1
    IntraUnits(nIntraUnits) = jAtom
  enddo

  ! Get atom type and reference coordinates
  atomType = atomArray(nType, Atm4)
  RefAtom%x = newMol%x(Atm3)
  RefAtom%y = newMol%y(Atm3)
  RefAtom%z = newMol%z(Atm3)
  r0 = r_eq
  call Find_InterAtms(RefAtom, r0, atomType, nIndx, nInterAtoms, InterAtoms)

  ! Compute geometric constraints for final atom
  r = GrowBonds(nGrowUnits)
  r1End = GrowBonds(nGrowBonds)
  v2%x = newMol%x(endAtom) - newMol%x(Atm3)
  v2%y = newMol%y(endAtom) - newMol%y(Atm3)
  v2%z = newMol%z(endAtom) - newMol%z(Atm3)
  r0End = sqrt(v2%x * v2%x + v2%y * v2%y + v2%z * v2%z)

  ! Check triangle inequality
  if (r + r1End < r0End) then
    rejMove = .true.
    return
  endif

  ! Compute Jacobian and fixed bending angle
  jacobian = jacobian * (1.0_dp / (r * r1End * r0End))
  costhmax = (r * r + r0End * r0End - r1End * r1End) / (2.0_dp * r * r0End)
  bend_angle = acos(costhmax)

  ! Generate trial positions with fixed angle
  overlap = .false.
  E_Trial = 0.0_dp
  do iRosen = 1, nRosenTrials(nType)
    call Generate_UnitCone(v2, r, bend_angle, v3)
    trialPos(iRosen)%x = v3%x + newMol%x(Atm3)
    trialPos(iRosen)%y = v3%y + newMol%y(Atm3)
    trialPos(iRosen)%z = v3%z + newMol%z(Atm3)

    ! Calculate intermolecular non-bonded energy
    call RosenLL_Atom_New(nType, Atm4, trialPos(iRosen), nInterAtoms, InterAtoms, &
                          E_Trial(iRosen), overlap(iRosen))

    ! Calculate intramolecular non-bonded energy if no overlap
    if (.not. overlap(iRosen)) then
      call RosenLL_Atom_Intra_New(nType, Atm4, trialPos(iRosen), nIntraUnits, IntraUnits, &
                                  E_Trial_Intra, overlap(iRosen))
      E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_Intra
    endif
  enddo

  ! Inline ChooseRosenTrial: Select trial position based on Rosenbluth weights
  E_Min = minval(E_Trial)
  ProbRosen = 0.0_dp
  do iRosen = 1, nRosenTrials(nType)
    ProbRosen(iRosen) = exp(-beta * (E_Trial(iRosen) - E_Min))
  enddo

  ! Check for valid trials
  if (all(ProbRosen <= 0.0_dp)) then
    rejMove = .true.
    return
  endif

  ! Normalize and select trial
  rosenNorm = sum(ProbRosen)
  ranNum = grnd() * rosenNorm
  sumInt = ProbRosen(1)
  nSel = 1
  do while (sumInt < ranNum .and. nSel < nRosenTrials(nType))
    nSel = nSel + 1
    sumInt = sumInt + ProbRosen(nSel)
  enddo

  ! Check for overlap in selected trial
  if (overlap(nSel)) then
    rejMove = .true.
    return
  endif

  ! Update Rosenbluth weight and store selected position
  rosenRatio = rosenRatio * ProbRosen(nSel) / rosenNorm
  regrown(Atm4) = .true.
  newMol%x(Atm4) = trialPos(nSel)%x
  newMol%y(Atm4) = trialPos(nSel)%y
  newMol%z(Atm4) = trialPos(nSel)%z
  totalRegrown = totalRegrown + 1

  ! Verify all atoms are regrown
  if (any(.not. regrown)) then
    write(nout, *) "Error: Non-regrown atoms remain for molecule ", nIndx, " of type ", nType
    stop
  endif

  ! Apply Jacobian correction to Rosenbluth weight
  rosenRatio = rosenRatio / jacobian

end subroutine LongChain_InternalSections_ConfigGen
!============================================================
! Performs reverse Fixed Endpoint Configurational Bias Monte Carlo (FECBMC) regrowth of an internal
! chain section for long linear biomolecules in a grand canonical ensemble nucleation simulation.
! Regrows atoms between fixed start and end atoms to compute the Rosenbluth weight of the existing
! (old) conformation for detailed balance. Uses the current atom positions as the first trial and
! generates additional trials with a single bond type, flexible bending (sine distribution for
! intermediate atoms, fixed for the final atom), and uniform azimuthal angles. Computes
! intermolecular and intramolecular non-bonded energies, applies a Jacobian factor for constrained
! growth, and returns the Rosenbluth weight (`rosenRatio`) without updating coordinates. Terminates
! if geometric constraints are violated or non-regrown atoms remain.
! For details of the method look at
! Escobedo, F. A.; De Pablo, J. J. Extended Continuum Configurational Bias Monte Carlo Methods for Simulation of Flexible Molecules. 
! J Chem Phys 1995, 102 (6), 2636–2652. https://doi.org/10.1063/1.468695.
subroutine LongChain_InternalSections_ConfigGen_Reverse(nType, nMol, regrownIn, startAtom, endAtom, rosenRatio)
  USE VarPrecision, ONLY: dp, atomIntType
  USE SimParameters, ONLY: maxAtoms, beta
  USE Coords, ONLY: MolArray
  USE CoordinateTypes, ONLY: SimpleAtomCoords
  USE CoordinateFunctions, ONLY: Generate_UnitCone
  USE ForceField, ONLY: bondData, bondArray, atomArray, nonBondArray, nIntraNonBond, nAtoms
  USE EnergyPressurePointers, ONLY: Find_InterAtms, RosenLL_Atom_New, RosenLL_Atom_Old, RosenLL_Atom_Intra_Old
  USE CBMC_Variables, ONLY: pathArray, nRosenTrials, maxRosenTrial
  USE PairStorage, ONLY: nTotalAtoms
  USE ParallelVar, ONLY: nout
  implicit none

  ! Input/Output variables
  integer, intent(in) :: nType                      ! Molecule type index
  integer, intent(in) :: nMol                       ! Molecule instance index
  logical, intent(in) :: regrownIn(maxAtoms)        ! Input regrown status of atoms
  integer, intent(in) :: startAtom, endAtom         ! Fixed start and end atom indices
  real(dp), intent(out) :: rosenRatio              ! Rosenbluth weight for old conformation

  ! Local variables
  integer :: i, j, iRosen, iAtom, jAtom               ! Loop indices
  integer :: Atm2, Atm3, Atm4                      ! Atom indices for chain geometry
  integer :: bondType                              ! Bond type index (single for molecule type)
  integer :: curPos, lastPos                       ! Current and last path positions
  integer :: Incrmt                                ! Increment direction (+1 or -1)
  integer :: totalRegrown                          ! Count of regrown atoms
  integer :: nIntraUnits                           ! Number of intramolecular non-bonded units
  integer :: IntraUnits(maxAtoms)                  ! Indices of intramolecular non-bonded atoms
  integer :: nInterAtoms                           ! Number of intermolecular interacting atoms
  integer :: InterAtoms(nTotalAtoms, 3)            ! Interacting atom indices
  integer :: nGrowUnits, nGrowBonds                ! Number of atoms and bonds to regrow
  integer :: iGrowUnits                            ! Current regrowth unit index
  integer(kind=atomIntType) :: atomType            ! Atom type for current atom
  logical :: regrown(maxAtoms)                     ! Regrown status of atoms
  logical :: overlap(maxRosenTrial)                ! Overlap flags for trial positions
  real(dp) :: r, r0End, r1End                      ! Bond length, distances to endpoint
  real(dp) :: r0                                   ! Buffer distance for intermolecular interactions
  real(dp) :: bend_angle, costhmax                 ! Bending angle and max cosine for constraint
  real(dp) :: jacobian                             ! Jacobian factor for constrained growth
  real(dp) :: GrowBonds(maxAtoms)                  ! Bond lengths from existing configuration
  real(dp) :: E_Trial(maxRosenTrial)              ! Trial energies
  real(dp) :: E_Trial_Intra                  ! Intramolecular trial energy
  real(dp) :: E_Min, rosenNorm, ProbRosen(maxRosenTrial) ! Rosenbluth selection variables
  type(SimpleAtomCoords) :: trialPos(maxRosenTrial) ! Trial positions
  type(SimpleAtomCoords) :: RefAtom, v2, v3        ! Reference atom and vectors for geometry
  integer :: nIndx                                 ! Molecule global index
  real(dp) :: dx, dy, dz                           ! Distance components for bond length
  ! External function: Generates a uniform random number in (0,1) for trial selection
  real(dp), external :: grnd
  
  ! Initialize variables
  rosenRatio = 1.0_dp
  jacobian = 1.0_dp
  regrown = regrownIn
  totalRegrown = 0
  nIndx = MolArray(nType)%mol(nMol)%indx

  ! Set bond type (single bond type for molecule type)
  bondType = bondArray(nType, 1)%bondType

  ! Count initially regrown atoms
  do i = 1, nAtoms(nType)
    if (regrown(i)) totalRegrown = totalRegrown + 1
  enddo

  ! Find path positions for start and end atoms
  curPos = 0
  lastPos = 0
  do i = 1, pathArray(nType)%pathMax(1)
    if (pathArray(nType)%path(1, i) == startAtom) curPos = i
    if (pathArray(nType)%path(1, i) == endAtom) lastPos = i
  enddo
  if (curPos == 0 .or. lastPos == 0) then
    write(nout, *) "Error: Start atom ", startAtom, " or end atom ", endAtom, &
                " not found in path for molecule type ", nType
    stop
  endif

  ! Determine regrowth direction and number of bonds/units
  if (lastPos > curPos) then
    Incrmt = 1
    nGrowBonds = lastPos - curPos
  else
    Incrmt = -1
    nGrowBonds = curPos - lastPos
  endif
  nGrowUnits = nGrowBonds - 1

  ! Compute bond lengths from existing configuration
  i = 0
  do j = curPos, lastPos - Incrmt, Incrmt
    i = i + 1
    Atm3 = pathArray(nType)%path(1, j)
    Atm4 = pathArray(nType)%path(1, j + Incrmt)
    dx = MolArray(nType)%mol(nMol)%x(Atm4) - MolArray(nType)%mol(nMol)%x(Atm3)
    dy = MolArray(nType)%mol(nMol)%y(Atm4) - MolArray(nType)%mol(nMol)%y(Atm3)
    dz = MolArray(nType)%mol(nMol)%z(Atm4) - MolArray(nType)%mol(nMol)%z(Atm3)
    GrowBonds(i) = sqrt(dx * dx + dy * dy + dz * dz)
  enddo

  ! Regrow intermediate atoms
  curPos = curPos + Incrmt
  do iGrowUnits = 1, nGrowUnits - 1
    ! Exit if current position is out of bounds
    if (curPos < 1 .or. curPos > pathArray(nType)%pathMax(1)) then
      write(nout, *) "Error: Current position out of bounds for molecule ", nIndx, " of type ", nType
      stop
    endif

    ! Get atom indices for chain geometry
    Atm4 = pathArray(nType)%path(1, curPos)
    Atm3 = pathArray(nType)%path(1, curPos - Incrmt)
    Atm2 = pathArray(nType)%path(1, curPos - 2*Incrmt)

    ! Find intramolecular non-bonded atoms already regrown
    nIntraUnits = 0
    do i = 1, nIntraNonBond(nType)
      iAtom = nonBondArray(nType, i)%nonMembr(1)
      jAtom = nonBondArray(nType, i)%nonMembr(2)
      if (iAtom /= Atm4 .and. jAtom /= Atm4) cycle
      if (iAtom /= Atm4) then
        jAtom = iAtom
        iAtom = Atm4
      endif
      if (.not. regrown(jAtom)) cycle
      nIntraUnits = nIntraUnits + 1
      IntraUnits(nIntraUnits) = jAtom
    enddo

    ! Get atom type and reference atom coordinates
    atomType = atomArray(nType, Atm4)
    RefAtom%x = MolArray(nType)%mol(nMol)%x(Atm3)
    RefAtom%y = MolArray(nType)%mol(nMol)%y(Atm3)
    RefAtom%z = MolArray(nType)%mol(nMol)%z(Atm3)
    r0 = bondData(bondType)%r_eq
    call Find_InterAtms(RefAtom, r0, atomType, nIndx, nInterAtoms, InterAtoms)

    ! Compute geometric constraints
    r = GrowBonds(iGrowUnits)
    r1End = sum(GrowBonds(iGrowUnits + 1:nGrowUnits + 1))
    v2%x = MolArray(nType)%mol(nMol)%x(endAtom) - MolArray(nType)%mol(nMol)%x(Atm3)
    v2%y = MolArray(nType)%mol(nMol)%y(endAtom) - MolArray(nType)%mol(nMol)%y(Atm3)
    v2%z = MolArray(nType)%mol(nMol)%z(endAtom) - MolArray(nType)%mol(nMol)%z(Atm3)
    r0End = sqrt(v2%x * v2%x + v2%y * v2%y + v2%z * v2%z)

    ! Check triangle inequality for geometric feasibility
    if (r + r1End < r0End) then
      write(nout, *) "Error: Triangle inequality violated for molecule ", nIndx, " of type ", nType
      stop
    endif

    ! Compute maximum cosine for constrained angle
    if (r1End > r + r0End) then
      costhmax = -1.0_dp
    else
      costhmax = (r * r + r0End * r0End - r1End * r1End) / (2.0_dp * r * r0End)
      if (costhmax <= -1.0_dp) costhmax = -1.0_dp
    endif
    jacobian = jacobian * (1.0_dp - costhmax)

    ! Generate trial positions
    overlap = .false.
    E_Trial = 0.0_dp
    iRosen = 1
    trialPos(iRosen)%x = MolArray(nType)%mol(nMol)%x(Atm4)
    trialPos(iRosen)%y = MolArray(nType)%mol(nMol)%y(Atm4)
    trialPos(iRosen)%z = MolArray(nType)%mol(nMol)%z(Atm4)
    call RosenLL_Atom_Old(nType, nMol, Atm4, nInterAtoms, InterAtoms, E_Trial(iRosen))
    call RosenLL_Atom_Intra_Old(nType, nMol, Atm4, trialPos(iRosen), nIntraUnits, IntraUnits, E_Trial_Intra)
    E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_Intra

    do iRosen = 2, nRosenTrials(nType)
      bend_angle = acos(1.0_dp - (1.0_dp - costhmax) * grnd())
      call Generate_UnitCone(v2, r, bend_angle, v3)
      trialPos(iRosen)%x = v3%x + MolArray(nType)%mol(nMol)%x(Atm3)
      trialPos(iRosen)%y = v3%y + MolArray(nType)%mol(nMol)%y(Atm3)
      trialPos(iRosen)%z = v3%z + MolArray(nType)%mol(nMol)%z(Atm3)
      call RosenLL_Atom_New(nType, Atm4, trialPos(iRosen), nInterAtoms, InterAtoms, E_Trial(iRosen), overlap(iRosen))
      if (.not. overlap(iRosen)) then
        call RosenLL_Atom_Intra_Old(nType, nMol, Atm4, trialPos(iRosen), nIntraUnits, IntraUnits, E_Trial_Intra)
        E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_Intra
      endif
    enddo

    ! Inline ChooseRosenTrial: Compute Rosenbluth weight for existing position
    E_Min = minval(E_Trial)
    ProbRosen = 0.0_dp
    do iRosen = 1, nRosenTrials(nType)
      ProbRosen(iRosen) = exp(-beta * (E_Trial(iRosen) - E_Min))
    enddo

    ! Check for valid trials
    if (all(ProbRosen <= 0.0_dp)) then
      write(nout, *) "Error: No valid trials in reverse regrowth for molecule ", nIndx, " of type ", nType
      stop
    endif

    ! Compute Rosenbluth weight (select first trial, existing position)
    rosenNorm = sum(ProbRosen)
    rosenRatio = rosenRatio * ProbRosen(1) / rosenNorm

    ! Mark atom as regrown and move to next
    regrown(Atm4) = .true.
    curPos = curPos + Incrmt
    totalRegrown = totalRegrown + 1
  enddo

  ! Regrow final atom with fixed angle
  if (curPos < 1 .or. curPos > pathArray(nType)%pathMax(1)) then
    write(nout, *) "Error: Final atom position out of bounds for molecule ", nIndx, " of type ", nType
    stop
  endif

  Atm4 = pathArray(nType)%path(1, curPos)
  Atm3 = pathArray(nType)%path(1, curPos - Incrmt)
  Atm2 = pathArray(nType)%path(1, curPos - 2*Incrmt)

  ! Find intramolecular non-bonded atoms
  nIntraUnits = 0
  do i = 1, nIntraNonBond(nType)
    iAtom = nonBondArray(nType, i)%nonMembr(1)
    jAtom = nonBondArray(nType, i)%nonMembr(2)
    if (iAtom /= Atm4 .and. jAtom /= Atm4) cycle
    if (iAtom /= Atm4) then
      jAtom = iAtom
      iAtom = Atm4
    endif
    if (.not. regrown(jAtom)) cycle
    nIntraUnits = nIntraUnits + 1
    IntraUnits(nIntraUnits) = jAtom
  enddo

  ! Get atom type and reference coordinates
  atomType = atomArray(nType, Atm4)
  RefAtom%x = MolArray(nType)%mol(nMol)%x(Atm3)
  RefAtom%y = MolArray(nType)%mol(nMol)%y(Atm3)
  RefAtom%z = MolArray(nType)%mol(nMol)%z(Atm3)
  r0 = bondData(bondType)%r_eq
  call Find_InterAtms(RefAtom, r0, atomType, nIndx, nInterAtoms, InterAtoms)

  ! Compute geometric constraints for final atom
  r = GrowBonds(nGrowUnits)
  r1End = GrowBonds(nGrowBonds)
  v2%x = MolArray(nType)%mol(nMol)%x(endAtom) - MolArray(nType)%mol(nMol)%x(Atm3)
  v2%y = MolArray(nType)%mol(nMol)%y(endAtom) - MolArray(nType)%mol(nMol)%y(Atm3)
  v2%z = MolArray(nType)%mol(nMol)%z(endAtom) - MolArray(nType)%mol(nMol)%z(Atm3)
  r0End = sqrt(v2%x * v2%x + v2%y * v2%y + v2%z * v2%z)

  ! Check triangle inequality
  if (r + r1End < r0End) then
    write(nout, *) "Error: Triangle inequality violated for final atom in molecule ", nIndx, " of type ", nType
    stop
  endif

  ! Compute Jacobian and fixed bending angle
  jacobian = jacobian * (1.0_dp / (r * r1End * r0End))
  costhmax = (r * r + r0End * r0End - r1End * r1End) / (2.0_dp * r * r0End)
  bend_angle = acos(costhmax)

  ! Generate trial positions
  overlap = .false.
  E_Trial = 0.0_dp
  iRosen = 1
  trialPos(iRosen)%x = MolArray(nType)%mol(nMol)%x(Atm4)
  trialPos(iRosen)%y = MolArray(nType)%mol(nMol)%y(Atm4)
  trialPos(iRosen)%z = MolArray(nType)%mol(nMol)%z(Atm4)
  call RosenLL_Atom_Old(nType, nMol, Atm4, nInterAtoms, InterAtoms, E_Trial(iRosen))
  call RosenLL_Atom_Intra_Old(nType, nMol, Atm4, trialPos(iRosen), nIntraUnits, IntraUnits, E_Trial_Intra)
  E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_Intra

  do iRosen = 2, nRosenTrials(nType)
    call Generate_UnitCone(v2, r, bend_angle, v3)
    trialPos(iRosen)%x = v3%x + MolArray(nType)%mol(nMol)%x(Atm3)
    trialPos(iRosen)%y = v3%y + MolArray(nType)%mol(nMol)%y(Atm3)
    trialPos(iRosen)%z = v3%z + MolArray(nType)%mol(nMol)%z(Atm3)
    call RosenLL_Atom_New(nType, Atm4, trialPos(iRosen), nInterAtoms, InterAtoms, E_Trial(iRosen), overlap(iRosen))
    if (.not. overlap(iRosen)) then
      call RosenLL_Atom_Intra_Old(nType, nMol, Atm4, trialPos(iRosen), nIntraUnits, IntraUnits, E_Trial_Intra)
      E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_Intra
    endif
  enddo

  ! Inline ChooseRosenTrial: Compute Rosenbluth weight for existing position
  E_Min = minval(E_Trial)
  ProbRosen = 0.0_dp
  do iRosen = 1, nRosenTrials(nType)
    ProbRosen(iRosen) = exp(-beta * (E_Trial(iRosen) - E_Min))
  enddo

  ! Check for valid trials
  if (all(ProbRosen <= 0.0_dp)) then
    write(nout, *) "Error: No valid trials for final atom in reverse regrowth for molecule ", nIndx, " of type ", nType
    stop
  endif

  ! Compute Rosenbluth weight (select first trial, existing position)
  rosenNorm = sum(ProbRosen)
  rosenRatio = rosenRatio * ProbRosen(1) / rosenNorm

  ! Mark final atom as regrown
  regrown(Atm4) = .true.
  totalRegrown = totalRegrown + 1

  ! Verify all atoms are regrown
  if (any(.not. regrown)) then
    write(nout, *) "Error: Non-regrown atoms remain for molecule ", nIndx, " of type ", nType
    stop
  endif

  ! Apply Jacobian correction to Rosenbluth weight
  rosenRatio = rosenRatio / jacobian

end subroutine LongChain_InternalSections_ConfigGen_Reverse
!============================================================
end module LLCBMC_Module