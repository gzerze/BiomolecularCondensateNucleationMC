module LLAVBMC_CBMC
  implicit none
  public
contains
!=======================================================================
! Generates trial configurations for inserting a long linear biomolecule (e.g., protein or polymer) of 
! type nType in a grand canonical Monte Carlo nucleation simulation. The insertion occurs near a target 
! molecule (index nTarget, type nTargType), with segment nAtom placed close to segment nTargAtom of the 
! target. This subroutine is part of the Long Linear Aggregation-Volume-Bias Monte Carlo (LLAVBMC) 
! insertion move, called by LLAVBMC_EBias_Rosen_In. It builds the new molecule segment by segment, 
! starting with nAtom, and grows in a random direction (toward the N-terminus or C-terminus of the 
! chain). For each segment, it generates multiple trial positions (nSwapInTrials(nType)), evaluates 
! their energies, and selects one using the Rosenbluth scheme (a weighted random selection based on 
! energy). The output is the Rosenbluth weight ratio:
!   rosenRatio = [exp(-β * U_new) / W_new] / [exp(-β * U_old) / W_old]
! where U_new is the energy of the inserted configuration, W_new is the Rosenbluth weight (sum of 
! trial probabilities), U_old is the gas-phase energy, and W_old is the gas-phase Rosenbluth weight. 
! This ratio is used in the acceptance probability for the insertion move, ensuring detailed balance 
! with the deletion move (LLAVBMC_EBias_Rosen_Out). To handle very small probabilities for long 
! chains, logarithmic probabilities (LogrosenRatio, LogwReverse) are summed across segments, and 
! rosenRatio is computed as exp(LogrosenRatio - LogwReverse) for numerical stability. The minimum-
! distance cluster criterion (minDistCriteria) is supported indirectly via NeighborPairs for neighbor 
! identification. Compatible with all force fields (LJ_Q, Mpipi, HPS_single, HPS_piecewise, 
! HPS_cation_pi) and molecule types (rigid, small, linear, branched).
subroutine LongChain_RosenConfigGen(nType, nAtom, nTarget, nTargType, nTargAtom, rosenRatio, rejMove)
  use SimParameters,         only: maxAtoms, beta, Dist_Critr
  use CBMC_Variables,        only: nSwapInTrials, pathArray, maxSwapInTrial
  use Coords,                only: MolArray, newMol, subIndxList
  use ForceField,            only: atomArray, bondArray, bondData, nIntraNonBond, nonBondArray
  use VarPrecision,          only: dp
  USE CoordinateTypes, ONLY: SimpleAtomCoords
  use EnergyPressurePointers, only: Find_InterAtms, RosenLL_Atom_New, RosenLL_Atom_Intra_New
  use PairStorage,           only: nTotalAtoms
  use CoordinateFunctions,   only: Generate_UnitSphere, GenerateGaussianBondLength, Generate_UnitCone
  implicit none

  ! Inputs: Define where and how the new molecule is inserted
  integer, intent(in) :: nType                  ! Type of molecule to insert (e.g., protein type)
  integer, intent(in) :: nAtom                  ! Segment (residue) index to place near nTargAtom
  integer, intent(in) :: nTarget                ! Global index of the target molecule
  integer, intent(in) :: nTargType              ! Type of the target molecule
  integer, intent(in) :: nTargAtom              ! Segment index of the target molecule
  ! Outputs: Results of the configuration generation
  real(dp), intent(out) :: rosenRatio           ! Rosenbluth weight ratio (W_new / W_old) for acceptance
  logical, intent(out) :: rejMove               ! True if move is rejected (e.g., due to overlaps)

  ! Local variables: Manage regrowth, energies, and probabilities
  integer :: nTargetMol                         ! Type-specific index of nTarget for coordinate access
  integer :: iGrow, iAtom, jAtom, nSel   ! Indices for growth steps in the regrowth loop
  integer :: growPos, curPos, prePos            ! Segment indices: growing, current, and previous
  integer :: iRosen, i, j                       ! Loop indices for trials and non-bonded pairs
  integer :: atomType, bondType                 ! Atom type of growing segment; bond type for growth
  integer :: nInterAtoms                        ! Number of intermolecular atoms within cutoff
  integer :: InterAtoms(1:nTotalAtoms, 1:3)     ! Indices of intermolecular interacting atoms
  integer :: nIntraUnits                        ! Number of non-bonded intramolecular segments
  integer :: IntraUnits(1:maxAtoms)             ! Indices of non-bonded intramolecular segments
  integer :: CurAtoms(1:maxAtoms)               ! Current segment indices in regrowth path
  integer :: GrowAtoms(1:maxAtoms)              ! Growing segment indices in regrowth path
  integer :: PrevAtoms(1:maxAtoms)              ! Previous segment indices in regrowth path
  integer :: totalRegrown                       ! Count of successfully grown segments
  logical :: regrown(1:maxAtoms)                ! Tracks which segments have been grown
  logical :: overlap(1:maxSwapInTrial)          ! True if trial position overlaps with other atoms
  real(dp) :: E_Trial(1:maxSwapInTrial)         ! Intermolecular energies for trial positions
  real(dp) :: E_Trial_Intra                     ! Non-bonded intramolecular energy for a trial
  real(dp) :: ProbRosen(1:maxSwapInTrial)       ! Rosenbluth probabilities for trial positions
  real(dp) :: LogrosenRatio, LogwReverse        ! Log of forward and reverse Rosenbluth weights
  real(dp) :: k_bond, r_eq, r_sigma, rmax_sq    ! Bond parameters: stiffness, equilibrium length, etc.
  real(dp) :: r, r0, dx, dy, dz, bend_angle     ! Trial distance, reference distance, and coordinates
  real(dp) :: E_Max, rosenNorm, ranNum, sumInt  ! Max energy, normalization, random number, cumulative prob
  real(dp) :: Prob                              ! Probability from Gaussian bond length generation
  type(SimpleAtomCoords) :: trialPos(1:maxSwapInTrial)  ! Trial coordinates for each position
  type(SimpleAtomCoords) :: RefAtom, vPrev, vGrow      ! Reference atom, previous, and growing vectors

  ! External function: Generates a uniform random number in (0,1) for trial selection
  real(dp), external :: grnd

  ! Step 1: Initialize the new molecule and schedule segment regrowth
  ! Overview: Sets up the new molecule’s properties, clears its coordinates, and determines the order 
  ! of segment growth (starting from nAtom, randomly toward N- or C-terminus) using a predefined path.
  nTargetMol = subIndxList(nTarget)             ! Convert nTarget to type-specific index for coords
  newMol%molType = nType                        ! Set the molecule type for the new insertion
  newMol%x = 0.0_dp                             ! Initialize x-coordinates to zero
  newMol%y = 0.0_dp                             ! Initialize y-coordinates to zero
  newMol%z = 0.0_dp                             ! Initialize z-coordinates to zero
  regrown = .false.                             ! Mark all segments as not yet grown
  rejMove = .false.                             ! Initialize rejection flag as false
  totalRegrown = 0                              ! Initialize count of grown segments
  LogrosenRatio = 0.0_dp                        ! Initialize log of forward Rosenbluth weight
  call Schedule_LongChain_Growth(nType, nAtom, PrevAtoms, CurAtoms, GrowAtoms)  ! Set growth order

  ! Step 2: Place the first segment near the target’s nTargAtom
  ! Overview: Positions the first segment (nAtom) within a sphere of radius Dist_Critr around 
  ! nTargAtom. Generates nSwapInTrials(nType) trial positions, evaluates their intermolecular 
  ! energies, and selects one using the Rosenbluth scheme (weighted by exp(-β * E)).

  growPos = patharray(nType)%path(1, GrowAtoms(1))  ! Index of first segment to grow
  E_Trial = 0.0_dp                              ! Initialize trial energies
  overlap = .false.                             ! Initialize overlap flags
  atomType = atomArray(nType, growPos)          ! Get atom type of growing segment
  RefAtom%x = molArray(nTargType)%mol(nTargetMol)%x(nTargAtom)  ! Target segment x-coordinate
  RefAtom%y = molArray(nTargType)%mol(nTargetMol)%y(nTargAtom)  ! Target segment y-coordinate
  RefAtom%z = molArray(nTargType)%mol(nTargetMol)%z(nTargAtom)  ! Target segment z-coordinate
  r0 = Dist_Critr                               ! Search radius for interacting atoms
  call Find_InterAtms(RefAtom, r0, atomType, -1, nInterAtoms, InterAtoms)  ! Find nearby atoms
  do iRosen = 1, nSwapInTrials(nType)           ! Loop over trial positions
    r = Dist_Critr * grnd()**(1.0_dp / 3.0_dp)  ! Random distance within sphere (uniform volume)
    call Generate_UnitSphere(dx, dy, dz)        ! Random direction on unit sphere
    trialPos(iRosen)%x = r * dx + molArray(nTargType)%mol(nTargetMol)%x(nTargAtom)  ! Trial x
    trialPos(iRosen)%y = r * dy + molArray(nTargType)%mol(nTargetMol)%y(nTargAtom)  ! Trial y
    trialPos(iRosen)%z = r * dz + molArray(nTargType)%mol(nTargetMol)%z(nTargAtom)  ! Trial z
    call RosenLL_Atom_New(nType, growPos, trialPos(iRosen), nInterAtoms, InterAtoms, &
                          E_Trial(iRosen), overlap(iRosen))  ! Compute intermolecular energy
  end do
  E_Max = minval(E_Trial)                       ! Find minimum energy for numerical stability
  ProbRosen = 0.0_dp                            ! Initialize Rosenbluth probabilities
  do iRosen = 1, nSwapInTrials(nType)           ! Compute probabilities for each trial
    ProbRosen(iRosen) = exp(-beta * (E_Trial(iRosen) - E_Max))  ! Boltzmann factor
  end do
  if (all(ProbRosen <= 0.0_dp)) then            ! Check if all trials have zero probability
    rejMove = .true.                            ! Reject move if no valid trials
    return
  end if
  rosenNorm = sum(ProbRosen)                    ! Sum probabilities for normalization
  ranNum = grnd() * rosenNorm                   ! Random number for weighted selection
  sumInt = ProbRosen(1)                         ! Cumulative probability
  nSel = 1                                      ! Selected trial index
  do while (sumInt < ranNum .and. nSel < nSwapInTrials(nType))  ! Select trial
    nSel = nSel + 1
    sumInt = sumInt + ProbRosen(nSel)
  end do
  if (overlap(nSel)) then                       ! Check if selected trial overlaps
    rejMove = .true.                            ! Reject move if overlap
    return
  end if
  LogrosenRatio = LogrosenRatio + log(ProbRosen(nSel) / rosenNorm)  ! Update log probability
  regrown(growPos) = .true.                     ! Mark segment as grown
  totalRegrown = totalRegrown + 1               ! Increment grown segment count
  newMol%x(growPos) = trialPos(nSel)%x          ! Store selected x-coordinate
  newMol%y(growPos) = trialPos(nSel)%y          ! Store selected y-coordinate
  newMol%z(growPos) = trialPos(nSel)%z          ! Store selected z-coordinate

  ! Step 3: Grow the second segment using a Gaussian bond length
  ! Overview: Positions the second segment relative to the first, using a bond length drawn from a 
  ! Gaussian distribution and a random orientation. This ensures proper chain connectivity while 
  ! exploring possible configurations. Trials are evaluated for intermolecular energies only.
  if (patharray(nType)%pathmax(1) > 1) then
    bondType = bondArray(nType, 1)%bondType       ! Get bond type (assumed same for all segments)
    k_bond = bondData(bondType)%k_eq              ! Bond stiffness for positioning
    r_eq = bondData(bondType)%r_eq                ! Equilibrium bond length
    r_sigma = bondData(bondType)%r_sigma          ! Bond length standard deviation
    rmax_sq = bondData(bondType)%rmax_sq          ! Maximum squared bond length
    growPos = patharray(nType)%path(1, GrowAtoms(2))  ! Index of second segment to grow
    curPos = patharray(nType)%path(1, CurAtoms(2))    ! Index of current (first) segment
    E_Trial = 0.0_dp                              ! Initialize trial energies
    overlap = .false.                             ! Initialize overlap flags
    atomType = atomArray(nType, growPos)          ! Get atom type of growing segment
    RefAtom%x = newMol%x(curPos)                  ! Current segment x-coordinate
    RefAtom%y = newMol%y(curPos)                  ! Current segment y-coordinate
    RefAtom%z = newMol%z(curPos)                  ! Current segment z-coordinate
    r0 = bondData(bondType)%r_eq                  ! Reference distance for interacting atoms
    call Find_InterAtms(RefAtom, r0, atomType, -1, nInterAtoms, InterAtoms)  ! Find nearby atoms
    do iRosen = 1, nSwapInTrials(nType)           ! Loop over trial positions
      call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)  ! Random bond length
      call Generate_UnitSphere(dx, dy, dz)        ! Random direction
      trialPos(iRosen)%x = r * dx + newMol%x(curPos)  ! Trial x-coordinate
      trialPos(iRosen)%y = r * dy + newMol%y(curPos)  ! Trial y-coordinate
      trialPos(iRosen)%z = r * dz + newMol%z(curPos)  ! Trial z-coordinate
      call RosenLL_Atom_New(nType, growPos, trialPos(iRosen), nInterAtoms, InterAtoms, &
                          E_Trial(iRosen), overlap(iRosen))  ! Compute intermolecular energy
    end do
    E_Max = minval(E_Trial)                       ! Find minimum energy
    ProbRosen = 0.0_dp                            ! Initialize probabilities
    do iRosen = 1, nSwapInTrials(nType)           ! Compute probabilities
      ProbRosen(iRosen) = exp(-beta * (E_Trial(iRosen) - E_Max))
    end do
    if (all(ProbRosen <= 0.0_dp)) then            ! Check for zero probabilities
      rejMove = .true.
      return
    end if
    rosenNorm = sum(ProbRosen)                    ! Normalize probabilities
    ranNum = grnd() * rosenNorm                   ! Random number for selection
    sumInt = ProbRosen(1)                         ! Cumulative probability
    nSel = 1                                      ! Selected trial index
    do while (sumInt < ranNum .and. nSel < nSwapInTrials(nType))  ! Select trial
      nSel = nSel + 1
      sumInt = sumInt + ProbRosen(nSel)
    end do
    if (overlap(nSel)) then                       ! Check for overlap
      rejMove = .true.
      return
    end if
    LogrosenRatio = LogrosenRatio + log(ProbRosen(nSel) / rosenNorm)  ! Update log probability
    regrown(growPos) = .true.                     ! Mark segment as grown
    totalRegrown = totalRegrown + 1               ! Increment grown count
    newMol%x(growPos) = trialPos(nSel)%x          ! Store x-coordinate
    newMol%y(growPos) = trialPos(nSel)%y          ! Store y-coordinate
    newMol%z(growPos) = trialPos(nSel)%z          ! Store z-coordinate
  endif

  ! Step 4: Grow remaining segments with bond lengths and bending angles
  ! Overview: Grows all subsequent segments (3rd to last) along the chain, using Gaussian bond 
  ! lengths and bending angles relative to the previous segment’s orientation. Each segment’s trials 
  ! include both intermolecular and non-bonded intramolecular energies to account for interactions 
  ! within the growing chain. The Rosenbluth scheme selects the best trial for each segment.
  do iGrow = 3, patharray(nType)%pathmax(1)     ! Loop over remaining segments
    growPos = patharray(nType)%path(1, GrowAtoms(iGrow))  ! Growing segment index
    curPos = patharray(nType)%path(1, CurAtoms(iGrow))    ! Current segment index
    prePos = patharray(nType)%path(1, PrevAtoms(iGrow))   ! Previous segment index
    atomType = atomArray(nType, growPos)        ! Atom type of growing segment
    RefAtom%x = newMol%x(curPos)                ! Current segment x-coordinate
    RefAtom%y = newMol%y(curPos)                ! Current segment y-coordinate
    RefAtom%z = newMol%z(curPos)                ! Current segment z-coordinate
    r0 = bondData(bondType)%r_eq                ! Reference distance for interactions
    call Find_InterAtms(RefAtom, r0, atomType, -1, nInterAtoms, InterAtoms)  ! Find nearby atoms
    vPrev%x = newMol%x(prePos) - newMol%x(curPos)  ! Vector from current to previous segment
    vPrev%y = newMol%y(prePos) - newMol%y(curPos)  ! y-component of vector
    vPrev%z = newMol%z(prePos) - newMol%z(curPos)  ! z-component of vector
    E_Trial = 0.0_dp                            ! Initialize trial energies
    overlap = .false.                           ! Initialize overlap flags
    nIntraUnits = 0                             ! Initialize count of intramolecular units
    do i = 1, nIntraNonBond(nType)              ! Identify non-bonded intramolecular pairs
      if (all(nonBondArray(nType, i)%nonMembr /= growPos)) cycle  ! Skip if growPos not involved
      iAtom = nonBondArray(nType, i)%nonMembr(1)  ! First atom in non-bonded pair
      jAtom = nonBondArray(nType, i)%nonMembr(2)  ! Second atom in non-bonded pair
      if (iAtom /= growPos) then                ! Ensure iAtom is growPos
        j = iAtom
        iAtom = jAtom
        jAtom = j
      end if
      if (.not. regrown(jAtom)) cycle           ! Skip if paired atom not grown
      nIntraUnits = nIntraUnits + 1             ! Increment intramolecular unit count
      IntraUnits(nIntraUnits) = jAtom           ! Store paired atom index
    end do
    do iRosen = 1, nSwapInTrials(nType)         ! Generate trials for segment
      call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)  ! Bond length
      bend_angle = acos(1.0_dp - 2.0_dp * grnd())  ! Random bending angle (uniform in cos)
      call Generate_UnitCone(vPrev, r, bend_angle, vGrow)  ! Position relative to previous
      trialPos(iRosen)%x = vGrow%x + newMol%x(curPos)  ! Trial x-coordinate
      trialPos(iRosen)%y = vGrow%y + newMol%y(curPos)  ! Trial y-coordinate
      trialPos(iRosen)%z = vGrow%z + newMol%z(curPos)  ! Trial z-coordinate
      call RosenLL_Atom_New(nType, growPos, trialPos(iRosen), nInterAtoms, InterAtoms, &
                            E_Trial(iRosen), overlap(iRosen))  ! Intermolecular energy
      if (.not. overlap(iRosen)) then           ! Compute intramolecular energy if no overlap
        call RosenLL_Atom_Intra_New(nType, growPos, trialPos(iRosen), nIntraUnits, IntraUnits, &
                                    E_Trial_Intra, overlap(iRosen))  ! Intramolecular energy
        E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_Intra  ! Total trial energy
      end if
    end do
    E_Max = minval(E_Trial)                     ! Find minimum energy
    ProbRosen = 0.0_dp                          ! Initialize probabilities
    do iRosen = 1, nSwapInTrials(nType)         ! Compute probabilities
      ProbRosen(iRosen) = exp(-beta * (E_Trial(iRosen) - E_Max))
    end do
    if (all(ProbRosen <= 0.0_dp)) then          ! Check for zero probabilities
      rejMove = .true.
      return
    end if
    rosenNorm = sum(ProbRosen)                  ! Normalize probabilities
    ranNum = grnd() * rosenNorm                 ! Random number for selection
    sumInt = ProbRosen(1)                       ! Cumulative probability
    nSel = 1                                    ! Selected trial index
    do while (sumInt < ranNum .and. nSel < nSwapInTrials(nType))  ! Select trial
      nSel = nSel + 1
      sumInt = sumInt + ProbRosen(nSel)
    end do
    if (overlap(nSel)) then                     ! Check for overlap
      rejMove = .true.
      return
    end if
    LogrosenRatio = LogrosenRatio + log(ProbRosen(nSel) / rosenNorm)  ! Update log probability
    regrown(growPos) = .true.                   ! Mark segment as grown
    totalRegrown = totalRegrown + 1             ! Increment grown count
    newMol%x(growPos) = trialPos(nSel)%x        ! Store x-coordinate
    newMol%y(growPos) = trialPos(nSel)%y        ! Store y-coordinate
    newMol%z(growPos) = trialPos(nSel)%z        ! Store z-coordinate
  end do

  ! Step 5: Verify all segments have been grown
  ! Overview: Ensures the entire chain has been successfully constructed by checking the number of 
  ! grown segments against the expected total. If incomplete, the move is rejected.
  if (totalRegrown /= patharray(nType)%pathmax(1)) then
    write(*, *) 'Error: Chain type ', nType, ' has ', totalRegrown, ' grown segments, expected ', &
                patharray(nType)%pathmax(1)  ! Report incomplete growth
    rejMove = .true.                          ! Reject move
    return
  end if

  ! Step 6: Compute the reverse (gas-phase) Rosenbluth weight
  ! Overview: Calculates the Rosenbluth weight for the reverse move (deleting the molecule to the 
  ! gas phase) using the same growth scheme in isolation. The logarithmic difference between forward 
  ! and reverse weights gives the final rosenRatio, ensuring numerical stability.
  call LongChain_RosenConfigGen_GasPhase_Reverse(nType, PrevAtoms, CurAtoms, GrowAtoms, LogwReverse)  ! Compute reverse log weight
  rosenRatio = exp(LogrosenRatio - LogwReverse)  ! Compute final Rosenbluth ratio

end subroutine LongChain_RosenConfigGen
!=======================================================================
! Computes the Rosenbluth weight ratio for the deletion move (old to new state) in a grand canonical 
! Monte Carlo nucleation simulation, as part of LLAVBMC_EBias_Rosen_Out (Step 9). It regrows a 
! deleted molecule of type nType and molecule number nMol in the cluster (old state), starting from 
! segment nAtom, which is within Dist_Critr of segment nTargAtom of target molecule nTarget (type 
! nTargType). The ratio is calculated as rosenRatio = [exp(-β * U_old) / W_old] / [exp(-β * U_new) / W_new].
! - Old state: Regrows the molecule in the cluster, using the original conformation as the first trial 
!   and generating nSwapOutTrials(nType)-1 additional trials, always selecting the first trial.
! - New state: Calls LongChain_RosenConfigGen_GasPhase to regrow the molecule in the gas phase, 
!   computing E_GasPhase and LogwForward.
! Uses logarithmic probabilities (LogrosenRatio, LogwForward) for numerical stability. Compatible 
! with all force fields (Mpipi, LJ_Q, HPS_single, HPS_piecewise, HPS_cation_pi) via 
! RosenLL_Atom_Old/New/Intra_Old. Applies the minimum-distance cluster criterion (minDistCriteria).
subroutine LongChain_RosenConfigGen_Reverse(nType, nAtom, nMol, nTarget, nTargType, nTargAtom, E_GasPhase, rosenRatio)
  use SimParameters,         only: maxAtoms, beta, Dist_Critr
  use CBMC_Variables,        only: nSwapOutTrials, pathArray, maxSwapOutTrial
  use Coords,                only: MolArray, subIndxList
  use ForceField,            only: nIntraNonBond, nonBondArray, bondArray, bondData, atomArray
  use VarPrecision,          only: dp, atomIntType
  USE CoordinateTypes, ONLY: SimpleAtomCoords
  use EnergyPressurePointers, only: Find_InterAtms, RosenLL_Atom_New, RosenLL_Atom_Old, RosenLL_Atom_Intra_Old
  use PairStorage,           only: nTotalAtoms
  use CoordinateFunctions,   only: GenerateGaussianBondLength, Generate_UnitSphere, Generate_UnitCone
  implicit none

  ! Inputs: Define the molecule, starting segment, and target
  integer, intent(in) :: nType                  ! Type of deleted molecule (e.g., protein type)
  integer, intent(in) :: nAtom                  ! Starting segment index for regrowth
  integer, intent(in) :: nMol                   ! Molecule number of deleted molecule
  integer, intent(in) :: nTarget                ! Index of target molecule
  integer, intent(in) :: nTargType              ! Type of target molecule
  integer, intent(in) :: nTargAtom              ! Target segment index
  ! Outputs: Gas-phase energy and Rosenbluth ratio
  real(dp), intent(out) :: E_GasPhase           ! Gas-phase energy from LongChain_RosenConfigGen_GasPhase
  real(dp), intent(out) :: rosenRatio           ! Rosenbluth weight ratio: (exp(-β * U_old) / W_old) / (exp(-β * U_new) / W_new)

  ! Local variables: Manage regrowth and probability calculations
  integer(kind=atomIntType) :: atomType         ! Atom type of growing segment
  integer :: CurAtoms(1:maxAtoms), GrowAtoms(1:maxAtoms), PrevAtoms(1:maxAtoms)
  integer :: iGrow                              ! Growth step index
  integer :: growPos, curPos, prePos            ! Growing, current, and previous segment indices
  integer :: bondType                           ! Bond type for growth parameters
  integer :: totalRegrown                       ! Count of grown segments
  integer :: nIntraUnits                        ! Number of non-bonded intramolecular segments
  integer :: IntraUnits(1:maxAtoms)             ! Indices of non-bonded intramolecular segments
  integer :: nInterAtoms                        ! Number of intermolecular interacting atoms
  integer :: InterAtoms(1:nTotalAtoms,1:3)      ! Intermolecular interacting atoms (mol, type, atom)
  integer :: i, iRosen, iAtom, jAtom, j         ! Loop indices for trials and non-bonded pairs
  integer :: nIndx                              ! Index of deleted molecule
  integer :: nTargetMol                         ! Molecule number of target molecule
  logical :: regrown(1:maxAtoms)                ! Tracks grown segments
  logical :: overlap(1:maxSwapOutTrial)         ! Overlap flags for trials
  real(dp) :: E_Trial(1:maxSwapOutTrial)        ! Trial energies (inter- + intramolecular)
  real(dp) :: E_Trial_Intra                     ! Intramolecular energy for a trial
  real(dp) :: LogrosenRatio                     ! Log of Rosenbluth weight for old state
  real(dp) :: LogwForward                       ! Log of Rosenbluth weight for gas phase
  real(dp) :: k_bond, r_eq, r_sigma, rmax_sq    ! Bond parameters: stiffness, equilibrium length, etc.
  real(dp) :: r, r0, dx, dy, dz                 ! Trial bond length and spherical coordinates
  real(dp) :: bend_angle                        ! Bending angle for growth
  real(dp) :: E_Min, rosenNorm                  ! Minimum trial energy and probability normalization
  real(dp) :: Prob, ProbRosen(1:maxSwapOutTrial)! Probability from bond length generation
  type(SimpleAtomCoords) :: trialPos(1:maxSwapOutTrial)  ! Trial positions for segments
  type(SimpleAtomCoords) :: vPrev, vGrow        ! Vectors for previous and growing segments
  type(SimpleAtomCoords) :: RefAtom             ! Reference atom for intermolecular interactions

  ! External function: Generates a uniform random number in (0,1)
  real(dp), external :: grnd

  ! Step 1: Initialize regrowth and schedule growth path
  ! Overview: Sets up tracking arrays, maps target molecule index to molecule number, and defines 
  ! the growth order starting from nAtom to ensure detailed balance with the insertion move.
  regrown = .false.                             ! Mark all segments as not grown
  totalRegrown = 0                              ! Initialize count of grown segments
  LogrosenRatio = 0.0_dp                        ! Initialize log of Rosenbluth weight (old state)
  nTargetMol = subIndxList(nTarget)             ! Get molecule number of target molecule
  nIndx = molArray(nType)%mol(nMol)%indx        ! Get index of deleted molecule
  call Schedule_LongChain_Growth(nType, nAtom, PrevAtoms, CurAtoms, GrowAtoms)  ! Set growth order

  ! Step 2: Set bond parameters
  ! Overview: Retrieves bond parameters for generating trial positions, assuming the same bond type 
  ! for all segments (may need verification).


  ! Step 3: Regrow first segment relative to target atom
  ! Overview: Positions the first segment (nAtom) within Dist_Critr of nTargAtom, using the 
  ! original conformation as the first trial and generating additional trials randomly within a 
  ! sphere of radius Dist_Critr. Computes intermolecular energies.
  iGrow = 1                                     ! First growth step
  growPos = patharray(nType)%path(1, GrowAtoms(iGrow))  ! Index of first segment to grow
  atomType = atomArray(nType, growPos)          ! Atom type of growing segment
  RefAtom%x = molArray(nTargType)%mol(nTargetMol)%x(nTargAtom)  ! Target atom position
  RefAtom%y = molArray(nTargType)%mol(nTargetMol)%y(nTargAtom)
  RefAtom%z = molArray(nTargType)%mol(nTargetMol)%z(nTargAtom)
  r0 = Dist_Critr                               ! Cluster criterion distance
  call Find_InterAtms(RefAtom, r0, atomType, nIndx, nInterAtoms, InterAtoms)  ! Find interacting atoms
  E_Trial = 0.0_dp                              ! Initialize trial energies
  overlap = .false.                             ! Initialize overlap flags
  iRosen = 1                                    ! First trial (original conformation)
  trialPos(iRosen)%x = molArray(nType)%mol(nMol)%x(nAtom)  ! Original x-coordinate
  trialPos(iRosen)%y = molArray(nType)%mol(nMol)%y(nAtom)  ! Original y-coordinate
  trialPos(iRosen)%z = molArray(nType)%mol(nMol)%z(nAtom)  ! Original z-coordinate
  call RosenLL_Atom_Old(nType, nMol, nAtom, nInterAtoms, InterAtoms, E_Trial(iRosen))  ! Intermolecular energy
  do iRosen = 2, nSwapOutTrials(nType)          ! Generate additional trials
    r = Dist_Critr * grnd()**(1.0_dp/3.0_dp)    ! Random radius within Dist_Critr
    call Generate_UnitSphere(dx, dy, dz)        ! Random direction on unit sphere
    trialPos(iRosen)%x = r * dx + molArray(nTargType)%mol(nTargetMol)%x(nTargAtom)  ! Trial x-coordinate
    trialPos(iRosen)%y = r * dy + molArray(nTargType)%mol(nTargetMol)%y(nTargAtom)  ! Trial y-coordinate
    trialPos(iRosen)%z = r * dz + molArray(nTargType)%mol(nTargetMol)%z(nTargAtom)  ! Trial z-coordinate
    call RosenLL_Atom_New(nType, nAtom, trialPos(iRosen), nInterAtoms, InterAtoms, &
                          E_Trial(iRosen), overlap(iRosen))  ! Intermolecular energy and overlap
  end do
  E_Min = minval(E_Trial)                       ! Find minimum energy for stability
  ProbRosen = 0.0_dp                            ! Initialize probabilities
  do iRosen = 1, nSwapOutTrials(nType)          ! Compute Rosenbluth probabilities
    ProbRosen(iRosen) = exp(-beta * (E_Trial(iRosen) - E_Min))
  end do
  if (all(ProbRosen <= 0.0_dp)) then            ! Check for zero probabilities
    write(*, *) 'Error: Zero probabilities for segment ', growPos, ' in cluster regrowth'
    return
  end if
  rosenNorm = sum(ProbRosen)                    ! Normalize probabilities
  LogrosenRatio = LogrosenRatio + log(ProbRosen(1) / rosenNorm)  ! Update log probability (first trial)
  regrown(growPos) = .true.                     ! Mark segment as grown
  totalRegrown = totalRegrown + 1               ! Increment grown count

  ! Step 4: Regrow second segment with bond constraints
  ! Overview: Regrows the second segment relative to the first, using bond parameters and the 
  ! original conformation as the first trial. Computes intermolecular energies.
  if (patharray(nType)%pathmax(1) > 1) then
    bondType = bondArray(nType, 1)%bondType       ! Get bond type (assumed uniform across segments)
    k_bond = bondData(bondType)%k_eq              ! Bond stiffness
    r_eq = bondData(bondType)%r_eq                ! Equilibrium bond length
    r_sigma = bondData(bondType)%r_sigma          ! Bond length standard deviation
    rmax_sq = bondData(bondType)%rmax_sq          ! Maximum squared bond length
    iGrow = 2                                     ! Second growth step
    growPos = patharray(nType)%path(1, GrowAtoms(iGrow))  ! Index of second segment
    curPos = patharray(nType)%path(1, CurAtoms(iGrow))    ! Current segment index
    atomType = atomArray(nType, growPos)          ! Atom type of growing segment
    RefAtom%x = molArray(nType)%mol(nMol)%x(curPos)  ! Reference position (first segment)
    RefAtom%y = molArray(nType)%mol(nMol)%y(curPos)
    RefAtom%z = molArray(nType)%mol(nMol)%z(curPos)
    r0 = r_eq                                     ! Use equilibrium bond length
    call Find_InterAtms(RefAtom, r0, atomType, nIndx, nInterAtoms, InterAtoms)  ! Find interacting atoms
    E_Trial = 0.0_dp                              ! Initialize trial energies
    overlap = .false.                             ! Initialize overlap flags
    iRosen = 1                                    ! First trial (original conformation)
    trialPos(iRosen)%x = molArray(nType)%mol(nMol)%x(growPos)  ! Original x-coordinate
    trialPos(iRosen)%y = molArray(nType)%mol(nMol)%y(growPos)  ! Original y-coordinate
    trialPos(iRosen)%z = molArray(nType)%mol(nMol)%z(growPos)  ! Original z-coordinate
    call RosenLL_Atom_Old(nType, nMol, growPos, nInterAtoms, InterAtoms, E_Trial(iRosen))  ! Intermolecular energy
    if (overlap(iRosen)) then                     ! Check for overlap in original conformation
      write(*, *) 'Error: Chain ', nIndx, ' of type ', nType, ' has overlap at segment ', growPos
      return
    end if
    do iRosen = 2, nSwapOutTrials(nType)          ! Generate additional trials
      call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)  ! Random bond length
      call Generate_UnitSphere(dx, dy, dz)        ! Random direction
      trialPos(iRosen)%x = r * dx + molArray(nType)%mol(nMol)%x(curPos)  ! Trial x-coordinate
      trialPos(iRosen)%y = r * dy + molArray(nType)%mol(nMol)%y(curPos)  ! Trial y-coordinate
      trialPos(iRosen)%z = r * dz + molArray(nType)%mol(nMol)%z(curPos)  ! Trial z-coordinate
      call RosenLL_Atom_New(nType, growPos, trialPos(iRosen), nInterAtoms, InterAtoms, &
                          E_Trial(iRosen), overlap(iRosen))  ! Intermolecular energy and overlap
    end do
    E_Min = minval(E_Trial)                       ! Find minimum energy
    ProbRosen = 0.0_dp                            ! Initialize probabilities
    do iRosen = 1, nSwapOutTrials(nType)          ! Compute Rosenbluth probabilities
      ProbRosen(iRosen) = exp(-beta * (E_Trial(iRosen) - E_Min))
    end do
    if (all(ProbRosen <= 0.0_dp)) then            ! Check for zero probabilities
      write(*, *) 'Error: Zero probabilities for segment ', growPos, ' in cluster regrowth'
      return
    end if
    rosenNorm = sum(ProbRosen)                    ! Normalize probabilities
    LogrosenRatio = LogrosenRatio + log(ProbRosen(1) / rosenNorm)  ! Update log probability
    regrown(growPos) = .true.                     ! Mark segment as grown
    totalRegrown = totalRegrown + 1               ! Increment grown count
  endif

  ! Step 5: Regrow remaining segments with bond and angle constraints
  ! Overview: Regrows subsequent segments using bond lengths, bending angles, and both inter- and 
  ! intramolecular energies. Identifies non-bonded intramolecular pairs and computes energies for 
  ! all trials, selecting the original conformation.
  do iGrow = 3, patharray(nType)%pathmax(1)     ! Loop over remaining segments
    growPos = patharray(nType)%path(1, GrowAtoms(iGrow))  ! Growing segment index
    curPos = patharray(nType)%path(1, CurAtoms(iGrow))    ! Current segment index
    prePos = patharray(nType)%path(1, PrevAtoms(iGrow))   ! Previous segment index
    atomType = atomArray(nType, growPos)        ! Atom type of growing segment
    RefAtom%x = molArray(nType)%mol(nMol)%x(curPos)  ! Reference position (current segment)
    RefAtom%y = molArray(nType)%mol(nMol)%y(curPos)
    RefAtom%z = molArray(nType)%mol(nMol)%z(curPos)
    r0 = r_eq                                   ! Use equilibrium bond length
    call Find_InterAtms(RefAtom, r0, atomType, nIndx, nInterAtoms, InterAtoms)  ! Find interacting atoms
    nIntraUnits = 0                             ! Initialize intramolecular unit count
    do i = 1, nIntraNonBond(nType)              ! Identify non-bonded intramolecular pairs
      if (all(nonBondArray(nType, i)%nonMembr /= growPos)) cycle  ! Skip if growPos not involved
      iAtom = nonBondArray(nType, i)%nonMembr(1)  ! First atom in pair
      jAtom = nonBondArray(nType, i)%nonMembr(2)  ! Second atom in pair
      if (iAtom /= growPos) then                ! Ensure iAtom is growPos
        j = iAtom
        iAtom = jAtom
        jAtom = j
      end if
      if (.not. regrown(jAtom)) cycle           ! Skip if paired atom not grown
      nIntraUnits = nIntraUnits + 1             ! Increment unit count
      IntraUnits(nIntraUnits) = jAtom           ! Store paired atom index
    end do
    vPrev%x = molArray(nType)%mol(nMol)%x(prePos) - molArray(nType)%mol(nMol)%x(curPos)  ! Vector to previous
    vPrev%y = molArray(nType)%mol(nMol)%y(prePos) - molArray(nType)%mol(nMol)%y(curPos)
    vPrev%z = molArray(nType)%mol(nMol)%z(prePos) - molArray(nType)%mol(nMol)%z(curPos)
    E_Trial = 0.0_dp                            ! Initialize trial energies
    overlap = .false.                           ! Initialize overlap flags
    iRosen = 1                                  ! First trial (original conformation)
    trialPos(iRosen)%x = molArray(nType)%mol(nMol)%x(growPos)  ! Original x-coordinate
    trialPos(iRosen)%y = molArray(nType)%mol(nMol)%y(growPos)  ! Original y-coordinate
    trialPos(iRosen)%z = molArray(nType)%mol(nMol)%z(growPos)  ! Original z-coordinate
    call RosenLL_Atom_Old(nType, nMol, growPos, nInterAtoms, InterAtoms, E_Trial(iRosen))  ! Intermolecular energy
    call RosenLL_Atom_Intra_Old(nType, nMol, growPos, trialPos(iRosen), nIntraUnits, IntraUnits, &
                                E_Trial_Intra)  ! Intramolecular energy
    E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_Intra
    do iRosen = 2, nSwapOutTrials(nType)        ! Generate additional trials
      call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)  ! Random bond length
      bend_angle = acos(1.0_dp - 2.0_dp * grnd())  ! Random bending angle
      call Generate_UnitCone(vPrev, r, bend_angle, vGrow)  ! Position relative to previous
      trialPos(iRosen)%x = vGrow%x + molArray(nType)%mol(nMol)%x(curPos)  ! Trial x-coordinate
      trialPos(iRosen)%y = vGrow%y + molArray(nType)%mol(nMol)%y(curPos)  ! Trial y-coordinate
      trialPos(iRosen)%z = vGrow%z + molArray(nType)%mol(nMol)%z(curPos)  ! Trial z-coordinate
      call RosenLL_Atom_New(nType, growPos, trialPos(iRosen), nInterAtoms, InterAtoms, &
                            E_Trial(iRosen), overlap(iRosen))  ! Intermolecular energy and overlap
      if (.not. overlap(iRosen)) then           ! Compute intramolecular energy if no overlap
        call RosenLL_Atom_Intra_Old(nType, nMol, growPos, trialPos(iRosen), nIntraUnits, IntraUnits, &
                                    E_Trial_Intra)
        E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_Intra
      end if
    end do
    E_Min = minval(E_Trial)                     ! Find minimum energy
    ProbRosen = 0.0_dp                          ! Initialize probabilities
    do iRosen = 1, nSwapOutTrials(nType)        ! Compute Rosenbluth probabilities
      ProbRosen(iRosen) = exp(-beta * (E_Trial(iRosen) - E_Min))
    end do
    if (all(ProbRosen <= 0.0_dp)) then          ! Check for zero probabilities
      write(*, *) 'Error: Zero probabilities for segment ', growPos, ' in cluster regrowth'
      return
    end if
    rosenNorm = sum(ProbRosen)                  ! Normalize probabilities
    LogrosenRatio = LogrosenRatio + log(ProbRosen(1) / rosenNorm)  ! Update log probability
    regrown(growPos) = .true.                   ! Mark segment as grown
    totalRegrown = totalRegrown + 1             ! Increment grown count
  end do

  ! Step 6: Verify all segments grown
  ! Overview: Ensures the entire chain was regrown by checking the number of grown segments 
  ! against the expected total.
  if (totalRegrown /= patharray(nType)%pathmax(1)) then
    write(*, *) 'Error: Chain type ', nType, ' has ', totalRegrown, ' grown segments, expected ', &
                patharray(nType)%pathmax(1)
    return                                    ! Return gracefully
  end if

  ! Step 7: Compute gas-phase configuration
  ! Overview: Calls LongChain_RosenConfigGen_GasPhase to regrow the molecule in the gas phase 
  ! (new state), obtaining E_GasPhase and LogwForward.
  call LongChain_RosenConfigGen_GasPhase(nType, PrevAtoms, CurAtoms, GrowAtoms, E_GasPhase, LogwForward)

  ! Step 8: Calculate Rosenbluth weight ratio
  ! Overview: Combines LogrosenRatio (W_old) and LogwForward (W_new) to compute the final 
  ! rosenRatio, accounting for energy differences implicitly via trial energies.
  rosenRatio = exp(LogrosenRatio - LogwForward)  ! rosenRatio = W_old / W_new

end subroutine LongChain_RosenConfigGen_Reverse
!=======================================================================      
! Regrows a molecule of type nType in the gas phase (new conformation for the deletion move) in a 
! grand canonical Monte Carlo nucleation simulation, called by LongChain_RosenConfigGen_Reverse 
! within LLAVBMC_EBias_Rosen_Out. It computes the intramolecular energy (E_GasPhase) and the log 
! of the Rosenbluth weight (LogwForward, ln W_new) for the gas-phase configuration, contributing 
! to the Rosenbluth weight ratio: rosenRatio = [exp(-β * U_old) / W_old] / [exp(-β * U_new) / W_new].
! - Selects a random starting segment and schedules growth using Schedule_LongChain_Growth.
! - First segment is placed at gasConfig; second uses bond constraints; remaining segments use 
!   bond lengths, bending angles, and intramolecular energies via RosenLL_Atom_Intra_New.
! - Selects trial conformations for segments beyond the second based on Rosenbluth weights.
! Compatible with all force fields (Mpipi, LJ_Q, HPS_single, HPS_piecewise, HPS_cation_pi) via 
! RosenLL_Atom_Intra_New. Applies the minimum-distance cluster criterion (minDistCriteria) 
! indirectly via parent subroutines.
subroutine LongChain_RosenConfigGen_GasPhase(nType, PrevAtoms, CurAtoms, GrowAtoms, E_GasPhase, LogwForward)
  use SimParameters,         only: maxAtoms, beta
  use CBMC_Variables,        only: nSwapOutTrials, pathArray, maxSwapOutTrial
  use Coords,                only: newMol, gasConfig
  use ForceField,            only: nIntraNonBond, nonBondArray, nAtoms, bondArray, bondData
  use VarPrecision,          only: dp
  USE CoordinateTypes, ONLY: SimpleAtomCoords
  use EnergyPressurePointers, only: RosenLL_Atom_Intra_New
!  use CBMC_Utility,          only: Schedule_LongChain_Growth
  use CoordinateFunctions,   only: GenerateGaussianBondLength, Generate_UnitSphere, Generate_UnitCone
  implicit none

  ! Input: Molecule type
  integer, intent(in) :: nType                  ! Molecule type (e.g., protein type)
  integer, intent(in) :: CurAtoms(1:maxAtoms)               ! Current segments in growth path
  integer, intent(in) :: GrowAtoms(1:maxAtoms)              ! Growing segments in growth path
  integer, intent(in) :: PrevAtoms(1:maxAtoms)              ! Previous segments in growth path
  ! Outputs: Gas-phase energy and log Rosenbluth weight
  real(dp), intent(out) :: E_GasPhase           ! Intramolecular energy of gas-phase molecule
  real(dp), intent(out) :: LogwForward          ! Log of Rosenbluth weight (ln W_new)

  ! Local variables: Manage regrowth and probability calculations
!  integer :: nAtom                              ! Randomly selected starting segment
  integer :: bondType                           ! Bond type for growth parameters
  integer :: totalRegrown                       ! Count of grown segments
  integer :: iGrow                              ! Growth step index
  integer :: growPos, curPos, prePos            ! Growing, current, and previous segment indices
  integer :: nIntraUnits                        ! Number of non-bonded intramolecular segments
  integer :: IntraUnits(1:maxAtoms)             ! Indices of non-bonded intramolecular segments
  integer :: i, iRosen, nSel                    ! Loop indices for non-bonded pairs and trial selection
  integer :: iAtom, jAtom, j                    ! Indices for non-bonded pair assignments
  logical :: regrown(1:maxAtoms)                ! Tracks grown segments
  logical :: overlap(1:maxSwapOutTrial)         ! Overlap flags for trials
  logical :: rejMove                            ! Flag for rejecting move due to overlaps
  real(dp) :: E_Trial(1:maxSwapOutTrial)        ! Trial intramolecular energies
  real(dp) :: ProbRosen(1:maxSwapOutTrial)      ! Rosenbluth probabilities for trials
  real(dp) :: k_bond, r_eq, r_sigma, rmax_sq    ! Bond parameters: stiffness, equilibrium length, etc.
  real(dp) :: r, dx, dy, dz                     ! Trial bond length and spherical coordinates
  real(dp) :: bend_angle                        ! Bending angle for growth
  real(dp) :: E_Max, rosenNorm                  ! Maximum trial energy and probability normalization
  real(dp) :: ranNum, sumInt                    ! Random number and cumulative probability for trial selection
  real(dp) :: Prob                              ! Probability from bond length generation
  type(SimpleAtomCoords) :: trialPos(1:maxSwapOutTrial)  ! Trial positions for segments
  type(SimpleAtomCoords) :: vPrev, vGrow        ! Vectors for previous and growing segments

  ! External function: Generates a uniform random number in (0,1)
  real(dp), external :: grnd

  ! Step 1: Initialize regrowth
  ! Overview: Sets up tracking arrays, initializes newMol coordinates, selects a random starting 
  ! segment, and schedules the growth path to ensure detailed balance with the deletion move.
  regrown = .false.                             ! Mark all segments as not grown
  rejMove = .false.                             ! Initialize move rejection flag
  totalRegrown = 0                              ! Initialize count of grown segments
  E_GasPhase = 0.0_dp                           ! Initialize gas-phase energy
  LogwForward = 0.0_dp                          ! Initialize log of Rosenbluth weight
  newMol%molType = nType                        ! Set molecule type for newMol
  newMol%x = 0.0_dp                             ! Initialize coordinates
  newMol%y = 0.0_dp
  newMol%z = 0.0_dp
!  nAtom = floor(nAtoms(nType) * grnd() + 1.0_dp)  ! Randomly select starting segment
!  call Schedule_LongChain_Growth(nType, nAtom, PrevAtoms, CurAtoms, GrowAtoms)  ! Set growth order

  ! Step 2: Set bond parameters
  ! Overview: Retrieves bond parameters for generating trial positions, assuming the same bond type 
  ! for all segments (may need verification).


  ! Step 3: Place first segment
  ! Overview: Places the first segment at the position stored in gasConfig, contributing no energy 
  ! in the gas phase, and adjusts LogwForward for the uniform selection probability.
  iGrow = 1                                     ! First growth step
  growPos = patharray(nType)%path(1, GrowAtoms(iGrow))  ! Index of first segment to grow
  newMol%x(growPos) = gasConfig(nType)%x(growPos)  ! Set x-coordinate from gasConfig
  newMol%y(growPos) = gasConfig(nType)%y(growPos)  ! Set y-coordinate from gasConfig
  newMol%z(growPos) = gasConfig(nType)%z(growPos)  ! Set z-coordinate from gasConfig
  regrown(growPos) = .true.                     ! Mark segment as grown
  totalRegrown = totalRegrown + 1               ! Increment grown count
  LogwForward = LogwForward + log(1.0_dp / real(nSwapOutTrials(nType), dp))  ! Uniform trial probability

  ! Step 4: Grow second segment with bond constraints
  ! Overview: Grows the second segment relative to the first using a Gaussian bond length and 
  ! random orientation, contributing no energy, and adjusts LogwForward for trial selection.
  if (patharray(nType)%pathmax(1) > 1) then
    bondType = bondArray(nType, 1)%bondType       ! Get bond type (assumed uniform across segments)
    k_bond = bondData(bondType)%k_eq              ! Bond stiffness
    r_eq = bondData(bondType)%r_eq                ! Equilibrium bond length
    r_sigma = bondData(bondType)%r_sigma          ! Bond length standard deviation
    rmax_sq = bondData(bondType)%rmax_sq          ! Maximum squared bond length
    iGrow = 2                                     ! Second growth step
    growPos = patharray(nType)%path(1, GrowAtoms(iGrow))  ! Index of second segment
    curPos = patharray(nType)%path(1, CurAtoms(iGrow))    ! Current segment index
    call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)  ! Random bond length
    call Generate_UnitSphere(dx, dy, dz)          ! Random direction on unit sphere
    newMol%x(growPos) = r * dx + newMol%x(curPos)  ! Set x-coordinate
    newMol%y(growPos) = r * dy + newMol%y(curPos)  ! Set y-coordinate
    newMol%z(growPos) = r * dz + newMol%z(curPos)  ! Set z-coordinate
    regrown(growPos) = .true.                     ! Mark segment as grown
    totalRegrown = totalRegrown + 1               ! Increment grown count
    LogwForward = LogwForward + log(1.0_dp / real(nSwapOutTrials(nType), dp))  ! Uniform trial probability
  endif

  ! Step 5: Grow remaining segments with bond, angle, and intramolecular energy constraints
  ! Overview: Grows subsequent segments using bond lengths, bending angles, and intramolecular 
  ! energies (via RosenLL_Atom_Intra_New). Generates nSwapOutTrials(nType) trials, selects one 
  ! based on Rosenbluth weights, and updates E_GasPhase and LogwForward.
  do iGrow = 3, patharray(nType)%pathmax(1)     ! Loop over remaining segments
    growPos = patharray(nType)%path(1, GrowAtoms(iGrow))  ! Growing segment index
    curPos = patharray(nType)%path(1, CurAtoms(iGrow))    ! Current segment index
    prePos = patharray(nType)%path(1, PrevAtoms(iGrow))   ! Previous segment index
    vPrev%x = newMol%x(prePos) - newMol%x(curPos)  ! Vector to previous segment
    vPrev%y = newMol%y(prePos) - newMol%y(curPos)
    vPrev%z = newMol%z(prePos) - newMol%z(curPos)
    E_Trial = 0.0_dp                            ! Initialize trial energies
    overlap = .false.                           ! Initialize overlap flags
    nIntraUnits = 0                             ! Initialize intramolecular unit count
    do i = 1, nIntraNonBond(nType)              ! Identify non-bonded intramolecular pairs
      if (all(nonBondArray(nType, i)%nonMembr /= growPos)) cycle  ! Skip if growPos not involved
      iAtom = nonBondArray(nType, i)%nonMembr(1)  ! First atom in pair
      jAtom = nonBondArray(nType, i)%nonMembr(2)  ! Second atom in pair
      if (iAtom /= growPos) then                ! Ensure iAtom is growPos
        j = iAtom
        iAtom = jAtom
        jAtom = j
      end if
      if (.not. regrown(jAtom)) cycle           ! Skip if paired atom not grown
      nIntraUnits = nIntraUnits + 1             ! Increment unit count
      IntraUnits(nIntraUnits) = jAtom           ! Store paired atom index
    end do
    do iRosen = 1, nSwapOutTrials(nType)        ! Generate trials
      call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)  ! Random bond length
      bend_angle = acos(1.0_dp - 2.0_dp * grnd())  ! Random bending angle
      call Generate_UnitCone(vPrev, r, bend_angle, vGrow)  ! Position relative to previous
      trialPos(iRosen)%x = vGrow%x + newMol%x(curPos)  ! Trial x-coordinate
      trialPos(iRosen)%y = vGrow%y + newMol%y(curPos)  ! Trial y-coordinate
      trialPos(iRosen)%z = vGrow%z + newMol%z(curPos)  ! Trial z-coordinate
      call RosenLL_Atom_Intra_New(nType, growPos, trialPos(iRosen), nIntraUnits, IntraUnits, &
                                  E_Trial(iRosen), overlap(iRosen))  ! Intramolecular energy and overlap
    end do
    E_Max = minval(E_Trial)                     ! Find minimum energy for stability
    ProbRosen = 0.0_dp                          ! Initialize probabilities
    do iRosen = 1, nSwapOutTrials(nType)        ! Compute Rosenbluth probabilities
      ProbRosen(iRosen) = exp(-beta * (E_Trial(iRosen) - E_Max))
    end do
    if (all(ProbRosen <= 0.0_dp)) then          ! Check for zero probabilities
      rejMove = .true.
      LogwForward = huge(dp)                    ! Reject move
      return
    end if
    rosenNorm = sum(ProbRosen)                  ! Normalize probabilities
    ranNum = grnd() * rosenNorm                 ! Random number for trial selection
    sumInt = ProbRosen(1)                       ! Cumulative probability
    nSel = 1                                    ! Selected trial index
    do while (sumInt < ranNum .and. nSel < nSwapOutTrials(nType))  ! Select trial
      nSel = nSel + 1
      sumInt = sumInt + ProbRosen(nSel)
    end do
    if (overlap(nSel)) then                     ! Check for overlap in selected trial
      rejMove = .true.
      LogwForward = huge(dp)                    ! Reject move
      return
    end if
    E_GasPhase = E_GasPhase + E_Trial(nSel)     ! Accumulate selected trial energy
    LogwForward = LogwForward + log(ProbRosen(nSel) / rosenNorm)  ! Update log weight
    newMol%x(growPos) = trialPos(nSel)%x        ! Store selected coordinates
    newMol%y(growPos) = trialPos(nSel)%y
    newMol%z(growPos) = trialPos(nSel)%z
    regrown(growPos) = .true.                   ! Mark segment as grown
    totalRegrown = totalRegrown + 1             ! Increment grown count
  end do

  ! Step 6: Verify all segments grown
  ! Overview: Ensures the entire chain was regrown by checking the number of grown segments 
  ! against the expected total.
  if (totalRegrown /= patharray(nType)%pathmax(1)) then
    write(*, *) 'Error: Chain type ', nType, ' has ', totalRegrown, ' grown segments, expected ', &
                patharray(nType)%pathmax(1)
    return                                    ! Return gracefully
  end if

  ! Step 7: Finalize gas-phase energy
  ! Overview: Adjusts E_GasPhase by negating it (convention may need verification) to match the 
  ! expected sign for the deletion move.
  E_GasPhase = -E_GasPhase                      ! Negate energy (convention-specific)
end subroutine LongChain_RosenConfigGen_GasPhase
!=======================================================================      
! Computes the logarithmic Rosenbluth weight for the reverse move (new to old state) in a grand 
! canonical Monte Carlo nucleation simulation, regrowing a long linear biomolecule of type nType in 
! the gas phase, where only non-bonded intramolecular interactions are considered. Called by 
! LongChain_RosenConfigGen within LLAVBMC_EBias_Rosen_In, this subroutine calculates LogwReverse, 
! the log of the gas-phase Rosenbluth weight (W_old), used in the Rosenbluth ratio:
!   rosenRatio = [exp(-β * U_new) / W_new] / [exp(-β * U_old) / W_old]
! It starts from a randomly chosen segment (nAtom) and follows the growth path defined by 
! Schedule_LongChain_Growth. For each segment (from the third onward), nSwapInTrials(nType) trials 
! are generated: the first uses the stored gas-phase configuration (gasConfig), and the rest use 
! Gaussian bond lengths and bending angles. The first trial is always selected, but its probability 
! is computed relative to all trials. Non-bonded intramolecular energies are calculated via 
! RosenLL_Atom_Intra_Gas_Old. Logarithmic probabilities ensure numerical stability for long chains. 
! Compatible with all force fields (LJ_Q, Mpipi, HPS_single, HPS_piecewise, HPS_cation_pi) and the 
! minimum-distance cluster criterion (minDistCriteria).
subroutine LongChain_RosenConfigGen_GasPhase_Reverse(nType, PrevAtoms, CurAtoms, GrowAtoms, LogwReverse)
  use SimParameters,          only: maxAtoms, beta
  use ForceField,             only: nAtoms, nIntraNonBond, nonBondArray, bondArray, bondData
  use VarPrecision,           only: dp
  use EnergyPressurePointers, only: RosenLL_Atom_Intra_Gas_Old
  use CBMC_Variables,         only: pathArray, nSwapInTrials, maxSwapInTrial
  use Coords,                 only: gasConfig
  USE CoordinateTypes, ONLY: SimpleAtomCoords
  use CoordinateFunctions,    only: GenerateGaussianBondLength, Generate_UnitCone
  implicit none

  ! Inputs/Outputs
  integer, intent(in) :: nType                  ! Type of molecule to regrow (e.g., protein type)
  integer, intent(in) :: CurAtoms(1:maxAtoms)               ! Current segment indices in growth path
  integer, intent(in) :: GrowAtoms(1:maxAtoms)              ! Growing segment indices in growth path
  integer, intent(in) :: PrevAtoms(1:maxAtoms)              ! Previous segment indices in growth path
  real(dp), intent(out) :: LogwReverse          ! Log of gas-phase Rosenbluth weight (W_old)

  ! Local variables: Manage regrowth and probability calculations
  integer :: iGrow                              ! Growth step index
  integer :: growPos, curPos, prePos            ! Growing, current, and previous segment indices
!  integer :: nAtom                              ! Randomly chosen starting segment
  integer :: bondType                           ! Bond type for growth parameters
  integer :: totalRegrown                       ! Count of grown segments
  integer :: nIntraUnits                        ! Number of non-bonded intramolecular segments
  integer :: IntraUnits(1:maxAtoms)             ! Indices of non-bonded intramolecular segments
  integer :: iRosen, i, j                       ! Loop indices for trials and non-bonded pairs
  integer :: iAtom, jAtom                       ! Atom indices in non-bonded pairs
  logical :: regrown(1:maxAtoms)                ! Tracks grown segments
  logical :: overlap(1:maxSwapInTrial)          ! Overlap flags for trials (used for first trial)
  real(dp) :: E_Trial(1:maxSwapInTrial)         ! Trial intramolecular energies
  real(dp) :: ProbRosen(1:maxSwapInTrial)       ! Rosenbluth probabilities for trials
  real(dp) :: k_bond, r_eq, r_sigma, rmax_sq    ! Bond parameters: stiffness, equilibrium length, etc.
  real(dp) :: r, bend_angle                     ! Trial bond length and bending angle
  real(dp) :: E_Max, rosenNorm                  ! Maximum energy and probability normalization
  real(dp) :: Prob                              ! Probability from bond length generation
  type(SimpleAtomCoords) :: trialPos(1:maxSwapInTrial)  ! Trial positions for segments
  type(SimpleAtomCoords) :: vPrev, vGrow        ! Vectors for previous and growing segments

  ! External function: Generates a uniform random number in (0,1)
  real(dp), external :: grnd

  ! Step 1: Initialize regrowth and schedule growth path
  ! Overview: Sets up tracking arrays, selects a random starting segment, and defines the growth 
  ! order using the same scheduling as the forward move to ensure detailed balance.
  regrown = .false.                             ! Mark all segments as not grown
  totalRegrown = 0                              ! Initialize count of grown segments
  LogwReverse = 0.0_dp                          ! Initialize log of Rosenbluth weight
!  nAtom = floor(nAtoms(nType) * grnd() + 1.0_dp)  ! Randomly choose starting segment (1 to nAtoms)
!  call Schedule_LongChain_Growth(nType, nAtom, PrevAtoms, CurAtoms, GrowAtoms)  ! Set growth order

  ! Step 2: Set bond parameters
  ! Overview: Retrieves bond parameters for generating trial positions, assuming the same bond type 
  ! for all segments (may need verification).


  ! Step 3: Grow first segment (no interactions)
  ! Overview: Marks the first segment as grown with no energy calculations, as the gas phase has no 
  ! intermolecular interactions and minimal intramolecular effects for the first segment. Adjusts 
  ! LogwReverse for the trial probability.
  growPos = patharray(nType)%path(1, GrowAtoms(1))  ! Index of first segment to grow
  regrown(growPos) = .true.                     ! Mark segment as grown
  totalRegrown = totalRegrown + 1               ! Increment grown segment count
  LogwReverse = LogwReverse + log(1.0_dp / real(nSwapInTrials(nType), dp))  ! Account for trial choice

  ! Step 4: Grow second segment (no interactions)
  ! Overview: Marks the second segment as grown, again with no energy calculations, and updates 
  ! LogwReverse for the trial probability.
  if (patharray(nType)%pathmax(1) > 1) then
    bondType = bondArray(nType, 1)%bondType       ! Get bond type (assumed uniform across segments)
    k_bond = bondData(bondType)%k_eq              ! Bond stiffness
    r_eq = bondData(bondType)%r_eq                ! Equilibrium bond length
    r_sigma = bondData(bondType)%r_sigma          ! Bond length standard deviation
    rmax_sq = bondData(bondType)%rmax_sq          ! Maximum squared bond length
    growPos = patharray(nType)%path(1, GrowAtoms(2))  ! Index of second segment
    regrown(growPos) = .true.                     ! Mark segment as grown
    totalRegrown = totalRegrown + 1               ! Increment grown count
    LogwReverse = LogwReverse + log(1.0_dp / real(nSwapInTrials(nType), dp))  ! Update log probability
  endif

  ! Step 5: Grow remaining segments with intramolecular energies
  ! Overview: For each subsequent segment, generates nSwapInTrials(nType) trials, with the first 
  ! trial using the stored gas-phase configuration (gasConfig). Computes non-bonded intramolecular 
  ! energies for all trials, selects the first trial, and updates LogwReverse based on its probability.
  do iGrow = 3, patharray(nType)%pathmax(1)     ! Loop over remaining segments
    growPos = patharray(nType)%path(1, GrowAtoms(iGrow))  ! Growing segment index
    curPos = patharray(nType)%path(1, CurAtoms(iGrow))    ! Current segment index
    prePos = patharray(nType)%path(1, PrevAtoms(iGrow))   ! Previous segment index
    vPrev%x = gasConfig(nType)%x(prePos) - gasConfig(nType)%x(curPos)  ! Vector to previous
    vPrev%y = gasConfig(nType)%y(prePos) - gasConfig(nType)%y(curPos)  ! y-component
    vPrev%z = gasConfig(nType)%z(prePos) - gasConfig(nType)%z(curPos)  ! z-component
    E_Trial = 0.0_dp                            ! Initialize trial energies
    overlap = .false.                           ! Initialize overlap flags
    nIntraUnits = 0                             ! Initialize intramolecular unit count
    do i = 1, nIntraNonBond(nType)              ! Identify non-bonded intramolecular pairs
      if (all(nonBondArray(nType, i)%nonMembr /= growPos)) cycle  ! Skip if growPos not involved
      iAtom = nonBondArray(nType, i)%nonMembr(1)  ! First atom in pair
      jAtom = nonBondArray(nType, i)%nonMembr(2)  ! Second atom in pair
      if (iAtom /= growPos) then                ! Ensure iAtom is growPos
        j = iAtom
        iAtom = jAtom
        jAtom = j
      end if
      if (.not. regrown(jAtom)) cycle           ! Skip if paired atom not grown
      nIntraUnits = nIntraUnits + 1             ! Increment unit count
      IntraUnits(nIntraUnits) = jAtom           ! Store paired atom index
    end do
    ! First trial: Use stored gas-phase configuration
    iRosen = 1
    trialPos(iRosen)%x = gasConfig(nType)%x(growPos)  ! Gas-phase x-coordinate
    trialPos(iRosen)%y = gasConfig(nType)%y(growPos)  ! Gas-phase y-coordinate
    trialPos(iRosen)%z = gasConfig(nType)%z(growPos)  ! Gas-phase z-coordinate
    call RosenLL_Atom_Intra_Gas_Old(nType, growPos, trialPos(iRosen), nIntraUnits, IntraUnits, &
                                    E_Trial(iRosen))  ! Compute intramolecular energy
    if (overlap(iRosen)) then                   ! Check for self-overlap (should be rare)
      write(*, *) 'Error: Chain type ', nType, ' has self-overlap in gas phase at segment ', growPos
      return                                    ! Return gracefully
    end if
    ! Remaining trials: Generate new positions
    do iRosen = 2, nSwapInTrials(nType)         ! Generate additional trials
      call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)  ! Random bond length
      bend_angle = acos(1.0_dp - 2.0_dp * grnd())  ! Random bending angle
      call Generate_UnitCone(vPrev, r, bend_angle, vGrow)  ! Position relative to previous
      trialPos(iRosen)%x = vGrow%x + gasConfig(nType)%x(curPos)  ! Trial x-coordinate
      trialPos(iRosen)%y = vGrow%y + gasConfig(nType)%y(curPos)  ! Trial y-coordinate
      trialPos(iRosen)%z = vGrow%z + gasConfig(nType)%z(curPos)  ! Trial z-coordinate
      call RosenLL_Atom_Intra_Gas_Old(nType, growPos, trialPos(iRosen), nIntraUnits, IntraUnits, &
                                      E_Trial(iRosen))  ! Compute intramolecular energy
    end do
    E_Max = minval(E_Trial)                     ! Find minimum energy for stability
    ProbRosen = 0.0_dp                          ! Initialize probabilities
    do iRosen = 1, nSwapInTrials(nType)         ! Compute Rosenbluth probabilities
      ProbRosen(iRosen) = exp(-beta * (E_Trial(iRosen) - E_Max))
    end do
    if (all(ProbRosen <= 0.0_dp)) then          ! Check for zero probabilities
      write(*, *) 'Error: Zero probabilities for segment ', growPos, ' in gas phase'
      return
    end if
    rosenNorm = sum(ProbRosen)                  ! Normalize probabilities
    regrown(growPos) = .true.                   ! Mark segment as grown
    totalRegrown = totalRegrown + 1             ! Increment grown count
    LogwReverse = LogwReverse + log(ProbRosen(1) / rosenNorm)  ! Update log probability for first trial
  end do

  ! Step 6: Verify all segments grown
  ! Overview: Ensures the entire chain was regrown by checking the number of grown segments 
  ! against the expected total.
  if (totalRegrown /= patharray(nType)%pathmax(1)) then
    write(*, *) 'Error: Chain type ', nType, ' has ', totalRegrown, ' grown segments, expected ', &
                patharray(nType)%pathmax(1)
    return                                    ! Return gracefully
  end if
end subroutine LongChain_RosenConfigGen_GasPhase_Reverse
!=======================================================================
! Schedules the regrowth order for inserting a long linear biomolecule of type nType in a grand 
! canonical Monte Carlo nucleation simulation, starting from segment nAtom (e.g., a residue in a 
! protein). Called by LongChain_RosenConfigGen within LLAVBMC_EBias_Rosen_In, this subroutine 
! determines the sequence in which segments are grown, randomly choosing to proceed toward the 
! N-terminus or C-terminus of the chain. It populates three arrays for use in configuration 
! generation:
! - GrowAtoms: Indices of segments to grow at each step.
! - CurAtoms: Indices of the current (reference) segment for each growth step.
! - PrevAtoms: Indices of the previous segment for orientation (used for bending angles).
! The growth path is defined by patharray(nType)%path, which maps segment indices, and the total 
! number of segments is patharray(nType)%pathmax(1). The subroutine ensures all segments are 
! scheduled by tracking the minimum (imin) and maximum (imax) indices reached. Random direction 
! choices (50% chance each way) ensure diverse configurations, supporting the Rosenbluth scheme in 
! LongChain_RosenConfigGen. Compatible with all force fields (LJ_Q, Mpipi, HPS_single, 
! HPS_piecewise, HPS_cation_pi) and the minimum-distance cluster criterion (minDistCriteria).
subroutine Schedule_LongChain_Growth(nType, nAtom, PrevAtoms, CurAtoms, GrowAtoms)
  use SimParameters,   only: maxAtoms
  use CBMC_Variables,  only: pathArray
  use VarPrecision,    only: dp
  implicit none

  integer, intent(in)  :: nType
  integer, intent(in)  :: nAtom
  integer, intent(out) :: PrevAtoms(1:maxAtoms)
  integer, intent(out) :: CurAtoms(1:maxAtoms)
  integer, intent(out) :: GrowAtoms(1:maxAtoms)

  integer :: iAtom
  integer :: curPos
  integer :: imin, imax
  integer :: L, Nfree
  logical :: left_avail, right_avail
  real(dp) :: ranNum
  integer :: i

  real(dp), external :: grnd

  PrevAtoms = 0
  CurAtoms  = 0
  GrowAtoms = 0

  L = patharray(nType)%pathmax(1)

  ! Find starting segment position in path array
  curPos = 0
  do i = 1, L
    if (patharray(nType)%path(1,i) == nAtom) then
      curPos = i
      exit
    end if
  end do
  if (curPos == 0) then
    write(*,*) 'Error: Segment ', nAtom, ' not found in patharray for type ', nType
    return
  end if

  ! Step 1: first scheduled segment
  iAtom = 1
  GrowAtoms(iAtom) = curPos
  CurAtoms(iAtom)  = 0
  PrevAtoms(iAtom) = 0
  imin = curPos
  imax = curPos

  Nfree = 0

  ! Steps 2..L: random walk expansion
  do while (iAtom < L)

    left_avail  = (imin > 1)
    right_avail = (imax < L)

    if (.not.left_avail .and. .not.right_avail) then
      write(*,*) 'Error: No available growth direction, type=',nType,' iAtom=',iAtom,' imin=',imin,' imax=',imax
      return
    end if

    ! Count "free" steps (both choices available) BEFORE making the random choice
    if (left_avail .and. right_avail) Nfree = Nfree + 1

    ! Choose direction
    if (left_avail .and. right_avail) then
      ranNum = grnd()
      if (ranNum > 0.5_dp) then
        ! choose right
        imax = imax + 1
        iAtom = iAtom + 1
        GrowAtoms(iAtom) = imax
        CurAtoms(iAtom)  = imax - 1
        if (iAtom == 2) then
          PrevAtoms(iAtom) = 0
        else
          PrevAtoms(iAtom) = imax - 2
        end if
      else
        ! choose left
        imin = imin - 1
        iAtom = iAtom + 1
        GrowAtoms(iAtom) = imin
        CurAtoms(iAtom)  = imin + 1
        if (iAtom == 2) then
          PrevAtoms(iAtom) = 0
        else
          PrevAtoms(iAtom) = imin + 2
        end if
      end if

    else if (right_avail) then
      ! forced right
      imax = imax + 1
      iAtom = iAtom + 1
      GrowAtoms(iAtom) = imax
      CurAtoms(iAtom)  = imax - 1
      if (iAtom == 2) then
        PrevAtoms(iAtom) = 0
      else
        PrevAtoms(iAtom) = imax - 2
      end if

    else
      ! forced left
      imin = imin - 1
      iAtom = iAtom + 1
      GrowAtoms(iAtom) = imin
      CurAtoms(iAtom)  = imin + 1
      if (iAtom == 2) then
        PrevAtoms(iAtom) = 0
      else
        PrevAtoms(iAtom) = imin + 2
      end if
    end if

  end do

  ! Final sanity check
  if (imin /= 1 .or. imax /= L) then
    write(*,*) 'Error: Growth scheduling failed. imin=',imin,' imax=',imax,' expected 1 and ',L
    return
  end if

end subroutine Schedule_LongChain_Growth
!=======================================================================
end module LLAVBMC_CBMC