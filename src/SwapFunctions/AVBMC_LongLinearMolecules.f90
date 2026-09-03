!===================================================================================            
      module LLAVBMC_Module
      use VarPrecision
      contains
!===================================================================================    
! INSERTION AND DELETION MOVES IN CLUSTER SIMULATIONS (WITH SEGMENT REGROWTH)
!
! In simulations using the minimum-distance cluster criterion:
!
! - INSERTION:
!     * A random segment from the inserting molecule is placed within the
!       bound region (typically a sphere) of a randomly selected segment
!       from the target molecule.
!     * The rest of the molecule is regrown segment by segment.
!     * Growth can proceed toward either the N-terminus or C-terminus,
!       direction chosen randomly.
!
! - DELETION:
!     * A neighboring molecule of the target is selected for deletion.
!     * Among all N_pair segment pairs (one segment from target, one from
!       the neighbor), one is randomly selected.
!     * The regrowth of the deleting molecule is initiated from that pair
!       to calculate the generation probability of the old state.
!
! - DILUTE-PHASE REGROWTH:
!     * Required in both insertion and deletion:
!         - For insertion: regrow the old (dilute-phase) conformation.
!         - For deletion: regrow the new (dilute-phase) conformation.
!     * Only intramolecular non-bonded interactions are considered.
!
! ============================================================================
! INSERTION ACCEPTANCE PROBABILITY
!
! P_Acc^Ins = min(1,
!     [ (exp(β * ε_i,new) / Σ_j^{n+1} exp(β * ε_j,new)) *
!       (exp(β * U_i',new) / Σ_{j'}^{N_i,nei + 1} exp(β * U_j',new)) *
!       (1 / N_seg^Ins) * ρ * W_new * exp(-β * [f_bias(n+1) - f_bias(n)]) ]
!     /
!     [ (exp(α * U_i,old) / Σ_j^n exp(α * U_j,old)) *
!       (1 / V_in) * (1 / N_seg^Ins) * (1 / N_seg^Target) * W_old ]
! )
!
! Where:
!   - ε_i,new: energy of the highest-energy neighbor of insertion target
!   - U_i',new: energy of the selected neighbor for insertion
!   - W: Rosenbluth weight
!   - ρ: dilute-phase density = exp(βμ) / Λ^3
!   - Λ: thermal de Broglie wavelength = sqrt(h^2 / (2π m k_B T))
!   - V_in: volume of the bound region
!   - N_seg^Target: number of segments in the target molecule
!   - N_seg^Ins: number of segments in the inserted molecule
!
! ============================================================================
! DELETION ACCEPTANCE PROBABILITY
!
! P_Acc^Del = min(1,
!     [ (exp(α * U_i,new) / Σ_j^{n-1} exp(α * U_j,new)) *
!       (1 / V_in) * (1 / N_pair) * W_new *
!       exp(-β * [f_bias(n-1) - f_bias(n)]) ]
!     /
!     [ (exp(β * ε_i,old) / Σ_j^n exp(β * ε_j,old)) *
!       (exp(β * U_i',old) / Σ_{j'}^{N_i,nei} exp(β * U_j',old)) *
!       (1 / N_seg^Del) * ρ * W_old ]
! )
!
! Where:
!   - N_pair: number of neighbor segment pairs between target and deletee
!   - N_seg^Del: number of segments in the deleted molecule
!
! ============================================================================          
! Performs a swap move (insertion or deletion) for long linear biomolecules in a cluster simulation
! using the Long Linear AVBMC (LLAVBMC) move with segment regrowth. Randomly selects insertion
! (LLAVBMC_EBias_Rosen_In) or deletion (LLAVBMC_EBias_Rosen_Out) with equal probability (0.5).
! Updates total energy (E_T), acceptance count (acc_x), and attempt count (atmp_x). Uses the
! minimum-distance cluster criterion (minDistCriteria) for neighbor identification, where
! molecules are neighbors if any segment pair is within Dist_Critr. Compatible with all force
! fields (LJ_Q, Mpipi, HPS_single, HPS_piecewise, HPS_cation_pi).
subroutine LLAVBMC(E_T, acc_x, atmp_x)
  use SimParameters, only: prevMoveAccepted, maxMol
  use VarPrecision, only: dp
  implicit none

  ! Inputs/Outputs: Total energy, acceptance count, attempt count
  real(dp), intent(inout) :: E_T, acc_x, atmp_x

  ! External function for uniform random number in (0,1)
  real(dp), external :: grnd

  ! Reset previous move acceptance flag
  prevMoveAccepted = .false.

  ! Randomly choose insertion or deletion with equal probability
  if (grnd() < 0.5_dp) then
    call LLAVBMC_EBias_Rosen_In(E_T, maxMol, acc_x, atmp_x)
  else
    call LLAVBMC_EBias_Rosen_Out(E_T, maxMol, acc_x, atmp_x)
  end if
end subroutine LLAVBMC
!===================================================================================
! Performs an insertion move for a long linear biomolecule using the LLAVBMC move with segment
! regrowth. Selects a molecule type probabilistically, checks system capacity and boundary
! conditions, chooses a target molecule and atoms (likely segments, e.g., residues), and
! computes the biased selection probability. Places the selected atom within Dist_Critr of the
! target atom, regrows the molecule, checks cluster criteria, computes umbrella bias, calculates
! energy changes, evaluates reverse move probabilities, and determines acceptance. Updates
! system state (coordinates, energies, neighbor lists) if accepted. Uses minimum-distance
! cluster criterion (minDistCriteria) for neighbor identification. Compatible with all force
! fields (LJ_Q, Mpipi, HPS_single, HPS_piecewise, HPS_cation_pi).
subroutine LLAVBMC_EBias_Rosen_In(E_T, arrayMax, acc_x, atmp_x)
  use AVBMC_RejectionVar, only: boundaryRej, totalRej, ovrlapRej, dbalRej, critriaRej
  use AcceptRates, only: atmpSwapIn, atmpInSize, acptSwapIn, acptInSize, clusterCritRej
  use AVBMC_Module, only: swapProb          
  use SwapBoundary, only: boundaryFunction   
  use VarPrecision, only: dp
  use Coords, only: MolArray, newMol
  use SimParameters, only: NTotal, maxMol, nMolTypes, NPART, NPART_New, NTotal_New, &
                           distCriteria, minDistCriteria, beta, avbmc_vol, gas_dens, &
                           calcPressure, isActive, prevMoveAccepted, maxAtoms, pressure
  use ForceField, only: nAtoms
  use EnergyPressurePointers, only: SwapIn_ECalc, Quick_Nei_ECalc, SwapIn_PCalc
  use E_Interface_LJ_Q, only: Update_SubEnergies
  use UmbrellaSamplingNew, only: GetUmbrellaBias_SwapIn
  use DistanceCriteria, only: NeighborUpdate_SwapIn_Distance 
  use EnergyTables, only: ETable, E_gasIntra, E_NBond_Diff
  use EnergyCriteria, only: NeighborUpdate
  use NeighborTable, only: Insert_NewNeiETable, Insert_NewNeiETable_Distance_V2
  use PairStorage, only: UpdateDistArray, useDistStore
  use CoordinateTypes, only: NeighborDetails
  use LLAVBMC_CBMC, only: LongChain_RosenConfigGen
  implicit none

  ! Inputs/Outputs
  integer, intent(in) :: arrayMax               ! Maximum number of molecules (likely maxMol)
  real(dp), intent(inout) :: E_T                ! Total system energy
  real(dp), intent(inout) :: acc_x, atmp_x      ! Acceptance and attempt counters

  ! Local variables
  logical :: rejMove                            ! Rejection flag for boundary, overlap, or criteria
  integer :: i, iAtom                           ! Loop index for coordinate updates
  integer :: nType                              ! Type of molecule to insert
  integer :: nAtom                              ! Atom index of inserting molecule (likely segment)
  integer :: nTarget, nTargType, nTargMol       ! Target molecule instance, type, and global index
  integer :: nTargAtom                          ! Atom index of target molecule (likely segment)
  integer :: nIndx                              ! Global index for new molecule
  integer :: NDiff(1:nMolTypes)                 ! Change in molecule count by type
  real(dp) :: ProbTarg_In                       ! Biased selection probability for target
  real(dp) :: rosenRatio                        ! Rosenbluth weight ratio
  real(dp) :: biasDiff                          ! Umbrella bias difference
  real(dp) :: E_Inter, E_Intra                  ! Inter- and intramolecular energies
  real(dp) :: ProbTarg_Out                      ! Target selection probability for reverse move
  real(dp) :: ProbSel_Out                       ! Neighbor selection probability for reverse move
  real(dp) :: ranNum, sumInt                    ! Random number and cumulative probability
  real(dp) :: genProbRatio                      ! Probability ratio for acceptance
  real(dp) :: Boltzterm                         ! Boltzmann term for energy and bias
  real(dp) :: P_Diff                            ! Pressure difference for update
  real(dp) :: PairList(1:arrayMax)              ! Array for pair energies
  real(dp) :: dETable(1:arrayMax)               ! Energy table for differences
  real(dp) :: newNeiETable(1:arrayMax)          ! Neighbor energy table for new state
  type(NeighborDetails) :: NeighborDetailsNew(1:maxMol)  ! New neighbor details for minDistCriteria
  ! External function for uniform random number in (0,1)
  real(dp), external :: grnd

  ! Step 1: Check system kapasity
  if (NTotal == maxMol) then
    boundaryRej = boundaryRej + 1.0_dp
    totalRej = totalRej + 1.0_dp
    return
  end if

  ! Step 2: Select molecule type for insertion
  if (nMolTypes == 1) then
    nType = 1
  else
    ranNum = grnd()
    sumInt = swapProb(1)
    nType = 1
    do while (sumInt < ranNum .and. nType < nMolTypes)
      nType = nType + 1
      sumInt = sumInt + swapProb(nType)
    end do
  end if

  ! Step 3: Check boundary conditions
  NDiff = 0
  NDiff(nType) = 1
  rejMove = boundaryFunction(NPART, NDiff)
  if (rejMove) then
    boundaryRej = boundaryRej + 1.0_dp
    totalRej = totalRej + 1.0_dp
    return
  end if

  ! Step 4: Update attempt counters
  atmp_x = atmp_x + 1.0_dp
  atmpSwapIn(nType) = atmpSwapIn(nType) + 1.0_dp
  atmpInSize(NTotal) = atmpInSize(NTotal) + 1.0_dp

  ! Step 5: Select target molecule and atoms (likely segments), compute biased selection probability
  call LongLinearEBias_Insert_ChooseTarget(nType, nAtom, nTarget, nTargType, nTargMol, nTargAtom, ProbTarg_In)

  ! Step 6: Assign global index for the new molecule
  nIndx = molArray(nType)%mol(NPART(nType) + 1)%indx

  ! Step 7: Generate new molecule configuration and compute Rosenbluth weight ratio
  ! rosenRatio = [exp(-β * U_new) / W_new] / [exp(-β * U_old) / W_old]
  call LongChain_RosenConfigGen(nType, nAtom, nTarget, nTargType, nTargAtom, rosenRatio, rejMove)
  if (rejMove) then
    totalRej = totalRej + 1.0_dp
    ovrlapRej = ovrlapRej + 1.0_dp
    return
  end if

  ! Step 8: Check energy-based cluster criterion (if not using distance criteria)
  if (.not. distCriteria .and. .not. minDistCriteria) then
    rejMove = .false.
    call Quick_Nei_ECalc(nTargType, nTargMol, rejMove)
    if (rejMove) then
      totalRej = totalRej + 1.0_dp
      critriaRej = critriaRej + 1.0_dp
      clusterCritRej = clusterCritRej + 1.0_dp
      return
    end if
  end if

  ! Step 9: Calculate umbrella sampling bias for proposed insertion
  NPART_new = NPART + NDiff
  NTotal_New = NTotal + 1
  call GetUmbrellaBias_SwapIn(biasDiff, rejMove)
  if (rejMove) then
    boundaryRej = boundaryRej + 1.0_dp
    totalRej = totalRej + 1.0_dp
    return
  end if

  ! Allocate NeighborDetailsNew pairIndices
  do iAtom = 1, maxMol
    NeighborDetailsNew(iAtom)%nPairs = 0
    NeighborDetailsNew(iAtom)%pairIndices = 0
  enddo

  ! Step 10: Calculate energy differences for the inserted molecule
  E_Inter = 0.0_dp
  E_Intra = 0.0_dp
  call SwapIn_ECalc(E_Inter, E_Intra, PairList, dETable, rejMove, .true., NeighborDetailsNew)
  if (rejMove) then
    totalRej = totalRej + 1.0_dp
    ovrlapRej = ovrlapRej + 1.0_dp
    return
  end if


  ! Step 11: Calculate reverse move probability for deletion from the new state
  if (distCriteria .or. minDistCriteria) then
    call Insert_NewNeiETable_Distance_V2(nType, PairList, dETable, newNeiETable)
  else
    call Insert_NewNeiETable(nType, PairList, dETable, newNeiETable)
  end if
  call LongLinearEBias_Insert_ReverseProbTarget(nTarget, nType, newNeiETable, ProbTarg_Out)
  call LongLinearEBias_Insert_ReverseProbSel(nTarget, nType, dETable, ProbSel_Out)
       
  ! Step 12: Calculate acceptance probability and determine if move is accepted
  ! P_Acc^Ins = min(1, [ProbTarg_Out * ProbSel_Out * ρ * W_new * exp(-β * [f_bias(n+1) - f_bias(n)])] /
  !                     [(1 / V_in) * (1 / N_seg^Ins) * (1 / N_seg^Target) * W_old])
  genProbRatio = (ProbTarg_Out * ProbSel_Out * avbmc_vol * gas_dens(nType)) / (ProbTarg_In * rosenRatio)
  Boltzterm = exp(-beta * (E_Inter + E_NBond_Diff - E_gasIntra(nType)) + biasDiff)

  if (genProbRatio * Boltzterm > grnd()) then
    ! Update system state for accepted move
    if (calcPressure) then
      call SwapIn_PCalc(P_Diff)
      pressure = pressure + P_Diff
    end if
    acptSwapIn(nType) = acptSwapIn(nType) + 1.0_dp
    acptInSize(NTotal) = acptInSize(NTotal) + 1.0_dp
    ! Copy new molecule coordinates
    do i = 1, nAtoms(nType)
      molArray(nType)%mol(NPART(nType) + 1)%x(i) = newMol%x(i)
      molArray(nType)%mol(NPART(nType) + 1)%y(i) = newMol%y(i)
      molArray(nType)%mol(NPART(nType) + 1)%z(i) = newMol%z(i)
    end do
    E_T = E_T + E_Inter + E_Intra
    acc_x = acc_x + 1.0_dp
    isActive(nIndx) = .true.
    if (useDistStore) then
      call UpdateDistArray
    end if
    ! Update neighbor lists based on cluster criterion
    if (distCriteria .or. minDistCriteria) then
      call NeighborUpdate_SwapIn_Distance(PairList, nType, NeighborDetailsNew)
    else
      call NeighborUpdate(PairList, nIndx)
    end if
    NTotal = NTotal + 1
    ETable = ETable + dETable
    NPART(nType) = NPART(nType) + 1
    call Update_SubEnergies
    prevMoveAccepted = .true.
  else
    totalRej = totalRej + 1.0_dp
    dbalRej = dbalRej + 1.0_dp
  end if

end subroutine LLAVBMC_EBias_Rosen_In
!===================================================================================
! Handles the deletion move for a long linear biomolecule in a grand canonical Monte Carlo
! nucleation simulation. Selects a molecule type, checks system capacity and boundary
! conditions, and prepares to select a target molecule and neighbor for deletion.
! Computes the acceptance probability:
! P_Acc^Del = min(1,
!     [ (exp(α * U_i,new) / Σ_j^{n-1} exp(α * U_j,new)) * (1 / V_in) * (1 / N_pair) * W_new *
!       exp(-β * [f_bias(n-1) - f_bias(n)]) ] /
!     [ (exp(β * ε_i,old) / Σ_j^n exp(β * ε_j,old)) * (exp(β * U_i',old) / Σ_j'^(N_i,nei) exp(β * U_j',old)) *
!       (1 / N_seg^Del) * ρ * W_old ]
! )
! where N_pair is the number of neighbor segment pairs, and N_seg^Del is the number of
! segments in the deleted molecule. Updates system state (coordinates, energies, neighbor
! lists) if accepted. Uses minimum-distance cluster criterion (minDistCriteria) for neighbor
! identification. Compatible with all force fields (LJ_Q, Mpipi, HPS_single, HPS_piecewise,
! HPS_cation_pi) and molecule types (rigid, small, linear, branched).
subroutine LLAVBMC_EBias_Rosen_Out(E_T, arrayMax, acc_x, atmp_x)
  use AcceptRates, only: atmpSwapOut, acptSwapOut
  use AVBMC_RejectionVar, only: critriaRej_out, boundaryRej_out, totalRej_out, dbalRej_out
  use AVBMC_Module, only: swapProb              
  use SwapBoundary, only: boundaryFunction     
  use VarPrecision, only: dp
  use Coords, only: MolArray, newMol, subIndxList, gasConfig
  use SimParameters, only: NTotal, nMolTypes, NPART, avbmc_vol, gas_dens, beta, isActive, &
                           NPART_New, NTotal_New, calcPressure, pressure, P_Diff, prevMoveAccepted  
  use EnergyTables, only: ETable, E_gasIntra, E_NBond_Diff
  use ForceField, only: nAtoms
  use EnergyPressurePointers, only: SwapOut_ECalc, SwapOut_PCalc
  use E_Interface_LJ_Q, only: Update_SubEnergies
  use UmbrellaSamplingNew, only: GetUmbrellaBias_SwapOut
  use EnergyCriteria, only: SwapOut_ClusterCriteria
  use EnergyCriteria, only: NeighborUpdate_Delete
  use PairStorage, only: UpdateDistArray_SwapOut, useDistStore
  use CoordinateTypes, only: NeighborDetails
  use NeighborTable, only: Create_NeiETable
  use LLAVBMC_CBMC, only: LongChain_RosenConfigGen_Reverse
  implicit none

  ! Inputs/Outputs
  integer, intent(in) :: arrayMax               ! Maximum number of molecules (likely maxMol)
  real(dp), intent(inout) :: E_T                ! Total system energy
  real(dp), intent(inout) :: acc_x, atmp_x      ! Acceptance and attempt counters

  ! Local variables
  logical :: rejMove                            ! Rejection flag for boundary or criteria
  integer :: nType                              ! Type of molecule to delete
  integer :: nAtom                              ! Atom index of deleted molecule (likely segment)
  integer :: nTarget, nTargType, nTargMol       ! Target molecule instance, type, and global index
  integer :: nTargAtom                          ! Atom index of target molecule (likely segment)
  integer :: nIndx                              ! Global index of deleted molecule
  integer :: nSel                               ! Selected molecule index (to be used)
  integer :: NDiff(1:nMolTypes)                 ! Change in molecule count by type
  real(dp) :: ProbTargOut                       ! Target selection probability (old state)
  real(dp) :: ProbSel                           ! Neighbor selection probability (old state)
  real(dp) :: ProbTargIn                        ! Target selection probability (new state)
  real(dp) :: genProbRatio                      ! Probability ratio for acceptance
  real(dp) :: Boltzterm                         ! Boltzmann term for energy and bias
  real(dp) :: biasDiff                          ! Umbrella bias difference
  real(dp) :: E_Inter, E_Intra, E_GasPhase      ! Inter-, intra-, and gas-phase energies
  real(dp) :: rosenRatio                        ! Rosenbluth weight ratio
  real(dp) :: ranNum, sumInt                    ! Random number and cumulative probability
  real(dp) :: dETable(1:arrayMax)               ! Energy table for differences
  integer :: i, nMol
  ! External function for uniform random number in (0,1)
  real(dp), external :: grnd

  ! Step 1: Check system capacity
  if (NTotal == 1) then
    boundaryRej_out = boundaryRej_out + 1.0_dp
    totalRej_out = totalRej_out + 1.0_dp
    return
  end if

  ! Step 2: Select molecule type for deletion
  if (nMolTypes == 1) then
    nType = 1
  else
    ranNum = grnd()
    sumInt = swapProb(1)
    nType = 1
    do while (sumInt < ranNum .and. nType < nMolTypes)
      nType = nType + 1
      sumInt = sumInt + swapProb(nType)
    end do
  end if

  ! Step 3: Check boundary conditions
  NDiff = 0
  NDiff(nType) = -1
  rejMove = boundaryFunction(NPART, NDiff)
  if (rejMove) then
    boundaryRej_out = boundaryRej_out + 1.0_dp
    totalRej_out = totalRej_out + 1.0_dp
    return
  end if

  ! Step 4: Select target molecule and neighbor for deletion
  ! Uses: ProbTargOut = exp(β * ε_i,old) / Σ_j exp(β * ε_j,old)
  ! ProbSel = (exp(β * U_i',old) / Σ_j'^N_i,nei exp(β * U_j',old)) * (1 / nAtoms(nType)) / (1 / N_pair)
  call Create_NeiETable(nType)
  call LongLinearEBias_Remove_ChooseTarget(nType, nTarget, nTargType, nTargMol, ProbTargOut)
  call LongLinearEBias_Remove_ChooseNeighbor(nTarget, nTargAtom, nType, nSel, nAtom, ProbSel)

  ! Step 5: Get type-specific index for nSel
  nMol = subIndxList(nSel)

  ! Step 6: Increment attempt counters
  atmp_x = atmp_x + 1.0_dp
  atmpSwapOut(nType) = atmpSwapOut(nType) + 1.0_dp

  ! Step 7: Calculate umbrella sampling bias
  ! biasDiff = -β * [f_bias(n-1) - f_bias(n)]
  NPART_new = NPART + NDiff
  NTotal_New = NTotal - 1
  call GetUmbrellaBias_SwapOut(nType, nMol, biasDiff, rejMove)
  if (rejMove) then
    boundaryRej_out = boundaryRej_out + 1.0_dp
    totalRej_out = totalRej_out + 1.0_dp
    return
  end if

  ! Step 8: Check if deletion breaks the cluster
  rejMove = .false.
  call SwapOut_ClusterCriteria(nSel, rejMove)
  if (rejMove) then
    critriaRej_out = critriaRej_out + 1.0_dp
    totalRej_out = totalRej_out + 1.0_dp
    return
  end if

  ! Step 9: Calculate Rosenbluth weight ratio
  ! rosenRatio = (exp(-β * U_old) / W_old) / (exp(-β * U_new) / W_new)
  call LongChain_RosenConfigGen_Reverse(nType, nAtom, nMol, nTarget, nTargType, nTargAtom, E_GasPhase, rosenRatio)

  ! Step 10: Calculate energy differences for deletion
  E_Inter = 0.0_dp
  E_Intra = 0.0_dp
  call SwapOut_ECalc(E_Inter, E_Intra, nType, nMol, dETable, .true.)

  ! Step 11: Calculate reverse move target selection probability
  ! ProbTargIn = exp(α * U_i,new) / Σ_j^(n-1) exp(α * U_j,new) for nTarget
  call LongLinearEBias_Remove_ReverseProbTarget(nTarget, nSel, nType, dETable, ProbTargIn)

  ! Step 12: Calculate acceptance probability and decide move
  ! genProbRatio = (ProbTargIn * rosenRatio) / (ProbTargOut * ProbSel * V_in * ρ)
  ! Boltzterm = exp(-β * (E_Inter + E_NBond_Diff - E_GasPhase) + biasDiff)
  genProbRatio = (ProbTargIn * rosenRatio) / (ProbTargOut * ProbSel * avbmc_vol * gas_dens(nType))
  Boltzterm = exp(-beta * (E_Inter + E_NBond_Diff - E_GasPhase) + biasDiff)
  ! Step 13: Accept or reject move
  if (genProbRatio * Boltzterm > grnd()) then
    ! Update pressure if required
    if (calcPressure) then
      call SwapOut_PCalc(nType, nMol, P_Diff)
      pressure = pressure - P_Diff
    end if
    ! Increment acceptance counters
    acptSwapOut(nType) = acptSwapOut(nType) + 1.0_dp
    acc_x = acc_x + 1.0_dp
    ! Shift coordinates to fill deleted molecule’s slot
    molArray(nType)%mol(nMol)%x(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%x(1:nAtoms(nType))
    molArray(nType)%mol(nMol)%y(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%y(1:nAtoms(nType))
    molArray(nType)%mol(nMol)%z(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%z(1:nAtoms(nType))
    ! Update total energy
    E_T = E_T + E_Inter + E_Intra
    ! Update neighbor lists and distance storage
    nIndx = molArray(nType)%mol(nMol)%indx
    if (useDistStore) call UpdateDistArray_SwapOut(nType, nMol)
    call NeighborUpdate_Delete(nIndx, molArray(nType)%mol(NPART(nType))%indx)
    ! Deactivate last molecule
    isActive(molArray(nType)%mol(NPART(nType))%indx) = .false.
    ! Update energy table
    ETable = ETable - dETable
    ETable(nIndx) = ETable(molArray(nType)%mol(NPART(nType))%indx)
    ETable(molArray(nType)%mol(NPART(nType))%indx) = 0.0_dp
    ! Update gas-phase configuration
    do i = 1, nAtoms(nType)
      gasConfig(nType)%x(i) = newMol%x(i)
      gasConfig(nType)%y(i) = newMol%y(i)
      gasConfig(nType)%z(i) = newMol%z(i)
    end do
    E_gasIntra(nType) = -E_GasPhase
    ! Update molecule counts
    NPART(nType) = NPART(nType) - 1
    NTotal = NTotal - 1
    ! Update sub-energies and acceptance flag
    call Update_SubEnergies
    prevMoveAccepted = .true.
  else
    ! Increment rejection counters
    dbalRej_out = dbalRej_out + 1.0_dp
    totalRej_out = totalRej_out + 1.0_dp
  end if
end subroutine LLAVBMC_EBias_Rosen_Out
!===================================================================================
! Selects a target molecule for AVBMC insertion with biased probability in a grand canonical
! ensemble nucleation simulation. Prefers molecules with higher interaction energies (e.g.,
! cluster surface molecules) to improve insertion efficiency. Uses:
!   P_Target_Ins(i) = (exp(AlphaTargetIns * beta * (U_i - U_avg(type_i))) / Σ_j exp(AlphaTargetIns * beta * (U_j - U_avg(type_j)))) * (1/nAtoms(nTargType)) * (1/nAtoms(nInsType))
! where U_i is molecule i’s interaction energy (ETable), U_avg(type_i) is the average energy
! for its type, and AlphaTargetIns(nInsType, nTargType) is a user-defined bias parameter (e.g.,
! 0.15 * β). Outputs the target’s global index (nTarget), type (nTargType), instance (nTargMol),
! atom indices (nInsAtom, nTargAtom, likely segments), and selection probability (ProbSel).
! Compatible with all force fields (LJ_Q, Mpipi, HPS_single, HPS_piecewise, HPS_cation_pi).
subroutine LongLinearEBias_Insert_ChooseTarget(nInsType, nInsAtom, nTarget, nTargType, nTargMol, nTargAtom, ProbSel)
  use SimParameters,   only: maxMol, nMolTypes, NMAX, isActive, beta
  use EnergyTables,    only: ETable, AlphaTargetIns
  use Coords,          only: typeList
  use ForceField,      only: nAtoms
  use VarPrecision,    only: dp
  implicit none

  ! Inputs/Outputs
  integer, intent(in) :: nInsType               ! Type of molecule to insert
  integer, intent(out) :: nInsAtom              ! Atom index of inserting molecule (likely segment)
  integer, intent(out) :: nTarget               ! Global index of target molecule
  integer, intent(out) :: nTargType             ! Type of target molecule
  integer, intent(out) :: nTargMol              ! Instance index of target molecule
  integer, intent(out) :: nTargAtom             ! Atom index of target molecule (likely segment)
  real(dp), intent(out) :: ProbSel              ! Total selection probability

  ! Local variables
  integer :: i, iType                           ! Loop indices for molecules and types
  integer :: cnt(1:nMolTypes)                   ! Count of active molecules per type
  integer :: offset                             ! Offset for instance index calculation
  real(dp) :: avgE(1:nMolTypes)                 ! Average interaction energy per type
  real(dp) :: ProbTable(1:maxMol)               ! Unnormalized selection probabilities
  real(dp) :: norm                              ! Normalization factor for probabilities
  real(dp) :: ranNum, sumInt                    ! Random number and cumulative probability

  ! External function for uniform random number in (0,1)
  real(dp), external :: grnd

  ! Step 1: Initialize arrays
  ProbTable = 0.0_dp
  avgE = 0.0_dp
  cnt = 0

  ! Step 2: Compute average interaction energy per molecule type
  do i = 1, maxMol
    if (isActive(i)) then
      iType = typeList(i)
      avgE(iType) = avgE(iType) + ETable(i)
      cnt(iType) = cnt(iType) + 1
    end if
  end do
  do iType = 1, nMolTypes
    if (cnt(iType) > 0) then
      avgE(iType) = avgE(iType) / real(cnt(iType), dp)
    end if
  end do

  ! Step 3: Build unnormalized probability table using biased energies
  do i = 1, maxMol
    if (isActive(i)) then
      iType = typeList(i)
      ProbTable(i) = exp(AlphaTargetIns(nInsType, iType) * beta * (ETable(i) - avgE(iType)))
    end if
  end do

  ! Step 4: Normalize probabilities and select target molecule
  norm = sum(ProbTable(1:maxMol))
  ranNum = norm * grnd()
  sumInt = 0.0_dp
  nTarget = 0
  do while(sumInt .lt. ranNum)
    nTarget = nTarget + 1
    sumInt = sumInt + ProbTable(nTarget)
  enddo

  ! Step 5: Set target properties
  nTargType = typeList(nTarget)
  ProbSel = ProbTable(nTarget) / norm

  ! Step 6: Compute target instance within its type
  offset = 0
  do i = 1, nTargType - 1
    offset = offset + NMAX(i)
  end do
  nTargMol = nTarget - offset

  ! Step 7: Select atoms (likely segments) uniformly and adjust selection probability
  nTargAtom = floor(real(nAtoms(nTargType), dp) * grnd() + 1.0_dp)
  ProbSel = ProbSel * (1.0_dp / real(nAtoms(nTargType), dp))
  nInsAtom = floor(real(nAtoms(nInsType), dp) * grnd() + 1.0_dp)
!  ProbSel = ProbSel * (1.0_dp / real(nAtoms(nInsType), dp))
end subroutine LongLinearEBias_Insert_ChooseTarget
!=================================================================================  
! Computes the reverse move probability of selecting the target molecule for deletion
! in the new state (post-insertion) of a grand canonical ensemble nucleation simulation.
! Uses: P_Target_Del = exp(β * ε_nTarget,new) / Σ_j^(n+1) exp(β * ε_j,new), where
! ε_i,new is the highest neighbor interaction energy for molecule i (from newNeiETable).
! Supports minDistCriteria for biomolecules and all force fields (LJ_Q, Mpipi, HPS_single,
! HPS_piecewise, HPS_cation_pi). Called by LLAVBMC_EBias_Rosen_In to compute ProbTarg_Out.
subroutine LongLinearEBias_Insert_ReverseProbTarget(nTarget, nType, newNeiETable, ProbRev)
  use SimParameters,   only: maxMol, NPART, beta
  use EnergyTables,    only: neiCount, AlphaTargetDel
  use Coords,          only: MolArray, typeList
  use VarPrecision,    only: dp
  implicit none

  ! Inputs/Outputs
  integer, intent(in) :: nTarget                ! Global index of target molecule
  integer, intent(in) :: nType                  ! Type of inserted molecule
  real(dp), intent(in) :: newNeiETable(1:maxMol) ! Highest neighbor energy for each molecule
  real(dp), intent(out) :: ProbRev              ! Target deletion probability (ProbTarg_Out)

  ! Local variables
  integer :: i, iType                           ! Loop index
  integer :: nIndx                              ! Global index of inserted molecule
  real(dp) :: ProbTable(1:maxMol)               ! Unnormalized probabilities
  real(dp) :: norm                              ! Normalization factor
  real(dp) :: EMax                              ! Maximum neighbor energy for scaling

  ! Step 1: Compute global index of inserted molecule
  nIndx = MolArray(nType)%mol(NPART(nType) + 1)%indx  ! Index for consistency, unused here

  ! Step 2: Initialize probability table and find maximum energy
  ProbTable = 0.0_dp
  EMax = -huge(dp)
  do i = 1, maxMol
    if (neiCount(i) > 0) then
      EMax = max(EMax, newNeiETable(i))
    end if
  end do

  ! Step 3: Compute unnormalized probabilities
  ! ProbTable(i) = exp(β * (ε_i,new - EMax)) for active molecules with neighbors
  do i = 1, maxMol
    if (neiCount(i) > 0) then
      iType = typeList(i)
      if (AlphaTargetDel(nType,iType) * beta * (newNeiETable(i) - EMax) > -30.0_dp) then
        ProbTable(i) = exp(AlphaTargetDel(nType,iType) * beta * (newNeiETable(i) - EMax))
      end if
    end if
  end do

  ! Step 4: Normalize to compute target deletion probability
  norm = sum(ProbTable(1:maxMol))
  ProbRev = ProbTable(nTarget) / norm
end subroutine LongLinearEBias_Insert_ReverseProbTarget
!=================================================================================
! Computes the neighbor selection probability for the reverse move (deletion) in the new state
! (post-insertion) of a grand canonical ensemble nucleation simulation. Uses:
! P_Sel_Neighbor = (exp(β * U_i',new) / Σ_j'^(N_i,nei + 1) exp(β * U_j',new)) * (1 / nAtoms(nType)),
! where U_i',new is the interaction energy of the inserted molecule (nIndx) with the target
! (nTarget), and nAtoms(nType) is the number of segments. Since ETable(nIndx) = 0, the inserted
! molecule’s contribution simplifies to 1. Supports minDistCriteria for biomolecules and all
! force fields (LJ_Q, Mpipi, HPS_single, HPS_piecewise, HPS_cation_pi). Called by
! LLAVBMC_EBias_Rosen_In to compute ProbSel_Out.
subroutine LongLinearEBias_Insert_ReverseProbSel(nTarget, nType, dE, ProbRev)
  use SimParameters, only: maxMol, NPART, beta
  use Coords,        only: MolArray, NeighborList
  use EnergyTables,  only: ETable
  use ForceField,    only: nAtoms
  use VarPrecision,  only: dp
  implicit none

  ! Inputs/Outputs
  integer, intent(in) :: nTarget                ! Global index of target molecule
  integer, intent(in) :: nType                  ! Type of inserted molecule
  real(dp), intent(in) :: dE(1:maxMol)          ! Energy change for molecule i due to insertion
  real(dp), intent(out) :: ProbRev              ! Neighbor selection probability (ProbSel_Out)

  ! Local variables
  integer :: iMol                               ! Loop index for molecules of type nType
  integer :: iIndx, nIndx                       ! Global indices for neighbor and inserted molecule
  real(dp) :: norm                              ! Normalization factor for Boltzmann factors

  ! Step 1: Compute global index of inserted molecule
  nIndx = molArray(nType)%mol(NPART(nType) + 1)%indx

  ! Step 2: Initialize normalization
  norm = 0.0_dp

  ! Step 3: Sum Boltzmann factors for nTarget’s neighbors of type nType
  ! norm += exp(β * (ETable(j) + dE(j) - dE(nIndx))) for neighbors
  do iMol = 1, NPART(nType)
    iIndx = molArray(nType)%mol(iMol)%indx
    if (NeighborList(iIndx, nTarget)) then
      norm = norm + exp(beta * (ETable(iIndx) + dE(iIndx) - dE(nIndx)))
    end if
  end do

  ! Step 4: Add inserted molecule’s contribution
  ! Since ETable(nIndx) = 0, exp(β * (0 + dE(nIndx) - dE(nIndx))) = 1
  norm = norm + 1.0_dp

  ! Step 5: Compute neighbor selection probability
  ProbRev = (1.0_dp / norm)
!  ProbRev = (1.0_dp / norm) * (1.0_dp / real(nAtoms(nType), dp))
end subroutine LongLinearEBias_Insert_ReverseProbSel
!=================================================================================   
! Selects a target molecule for deletion in a grand canonical Monte Carlo nucleation simulation
! with probability P_Target = exp(β * ε_i,old) / Σ_j^n exp(β * ε_j,old), where ε_i,old is the
! highest interaction energy among molecule i’s neighbors of type nType (from NeiETable).
! Outputs nTarget (global index), nTargType (molecule type), nTargMol (type-specific index),
! and ProbTarget (selection probability, ProbTargOut). Called by LLAVBMC_EBias_Rosen_Out.
! Supports minDistCriteria for biomolecules and all force fields (LJ_Q, Mpipi, HPS_single,
! HPS_piecewise, HPS_cation_pi).
subroutine LongLinearEBias_Remove_ChooseTarget(nType, nTarget, nTargType, nTargMol, ProbTarget)
  use SimParameters, only: maxMol, beta, NMAX
  use EnergyTables,  only: neiCount, NeiETable, AlphaTargetDel
  use Coords,        only: typeList
  use VarPrecision,  only: dp
  implicit none

  ! Outputs
  integer, intent(in) :: nType
  integer, intent(out) :: nTarget               ! Global index of target molecule
  integer, intent(out) :: nTargType             ! Type of target molecule
  integer, intent(out) :: nTargMol              ! Type-specific index of target molecule
  real(dp), intent(out) :: ProbTarget           ! Target selection probability (ProbTargOut)

  ! Local variables
  integer :: i, iType                           ! Loop index
  real(dp) :: ProbTable(1:maxMol)               ! Unnormalized probabilities
  real(dp) :: norm                              ! Normalization factor
  real(dp) :: ranNum, sumInt                    ! Random number and cumulative probability

  ! External function for uniform random number in (0,1)
  real(dp), external :: grnd

  ! Step 1: Initialize probability table
  ! ProbTable(i) = exp(β * ε_i,old) for active molecules with neighbors, else 0
  ProbTable = 0.0_dp
  do i = 1, maxMol
    if (neiCount(i) > 0) then
      iType = typeList(i)
      ProbTable(i) = exp(AlphaTargetDel(nType,iType) * beta * NeiETable(i))
    end if
  end do

  ! Step 2: Normalize probabilities and select target molecule
  norm = sum(ProbTable(1:maxMol))
  ranNum = norm * grnd()
  sumInt = ProbTable(1)
  nTarget = 1
  do while (sumInt < ranNum .and. nTarget < maxMol)
    nTarget = nTarget + 1
    sumInt = sumInt + ProbTable(nTarget)
  end do

  ! Step 3: Set target type and selection probability
  nTargType = typeList(nTarget)
  ProbTarget = ProbTable(nTarget) / norm

  ! Step 4: Compute type-specific molecule index
  ! nTargMol = nTarget - offset of previous types
  nTargMol = 0
  do i = 1, nTargType - 1
    nTargMol = nTargMol + NMAX(i)
  end do
  nTargMol = nTarget - nTargMol
end subroutine LongLinearEBias_Remove_ChooseTarget
!=================================================================================  
! Selects a neighbor molecule of type nType to delete from nTarget’s neighbors in a grand
! canonical Monte Carlo nucleation simulation, with probability:
! P_Sel = (exp(β * U_i',old) / Σ_j'^N_i,nei exp(β * U_j',old)) * (1 / nAtoms(nType)) / (1 / N_pair),
! where U_i',old is the total interaction energy of molecule i' with nTarget (from ETable),
! nAtoms(nType) is N_seg^Del, and N_pair is the number of segment pairs within minDistCriteria.
! Outputs nSel (global index of molecule to delete), nAtom (segment of deleted molecule),
! nTargAtom (segment of target), and ProbSel (selection probability, ProbSel in
! LLAVBMC_EBias_Rosen_Out). Uses NeighborPairs for neighbor segment pairs.
! Supports all force fields (LJ_Q, Mpipi, HPS_single, HPS_piecewise, HPS_cation_pi).
subroutine LongLinearEBias_Remove_ChooseNeighbor(nTarget, nTargAtom, nType, nSel, nAtom, ProbSel)
  use SimParameters, only: maxMol, NPART, beta
  use Coords, only: MolArray, typeList, NeighborList, NeighborPairs
  use EnergyTables, only: ETable
  use ForceField, only: nAtoms
  use VarPrecision, only: dp
  implicit none

  ! Inputs/Outputs
  integer, intent(in) :: nTarget                ! Global index of target molecule
  integer, intent(in) :: nType                  ! Type of molecule to delete
  integer, intent(out) :: nTargAtom             ! Segment index of target molecule
  integer, intent(out) :: nSel                  ! Global index of molecule to delete
  integer, intent(out) :: nAtom                 ! Segment index of deleted molecule
  real(dp), intent(out) :: ProbSel              ! Neighbor selection probability (ProbSel)

  ! Local variables
  integer :: iMol, iIndx, iPair                 ! Loop index and molecule index
  integer :: nPairAtoms                         ! Number of neighbor segment pairs (N_pair)
  integer :: nSelPair                           ! Selected pair index
  integer :: nTargType                          ! Type of target molecule
  real(dp) :: ProbTable(1:maxMol)               ! Unnormalized probabilities for molecules
  real(dp), allocatable :: ProbAtom(:)          ! Unnormalized probabilities for segment pairs
  real(dp) :: norm                              ! Normalization for molecule selection
  real(dp) :: normPairs                         ! Normalization for pair selection
  real(dp) :: ranNum, sumInt                    ! Random number and cumulative probability

  ! External function for uniform random number in (0,1)
  real(dp), external :: grnd

  ! Step 1: Initialize probability table for neighbor molecules
  ProbTable = 0.0_dp
  nTargType = typeList(nTarget)
  do iMol = 1, NPART(nType)
    iIndx = MolArray(nType)%mol(iMol)%indx
    if (NeighborList(iIndx, nTarget)) then
      ProbTable(iIndx) = exp(beta * ETable(iIndx))
    end if
  end do

  ! Step 2: Select neighbor molecule using cumulative probability
  norm = sum(ProbTable(1:maxMol))
  ranNum = norm * grnd()
  sumInt = ProbTable(1)
  nSel = 1
  do while (sumInt < ranNum .and. nSel < maxMol)
    nSel = nSel + 1
    sumInt = sumInt + ProbTable(nSel)
  end do
  ProbSel = ProbTable(nSel) / norm

  ! Step 3: Get number of neighbor segment pairs from NeighborPairs
  nPairAtoms = NeighborPairs(nSel, nTarget)%details%nPairs
  if (nPairAtoms == 0) then
    write(*, *) 'Error: nSel = ', nSel, ' and nTarget = ', nTarget, ' are not neighbors in NeighborPairs'
    ProbSel = 0.0_dp
    nSel = 0
    nAtom = 0
    nTargAtom = 0
    stop
  end if

  ! Step 4: Select segment pair
  allocate(ProbAtom(1:nPairAtoms))
  ProbAtom = 0.0_dp
  if (nPairAtoms == 1) then
    ! Single pair: select directly
    nAtom = NeighborPairs(nSel, nTarget)%details%pairIndices(1)
    nTargAtom = NeighborPairs(nSel, nTarget)%details%pairIndices(2)
!    ProbSel = ProbSel * (1.0_dp / real(nAtoms(nType), dp))
  else
    ! Multiple pairs: select uniformly
    do iPair = 1, nPairAtoms
      ProbAtom(iPair) = 1.0_dp
    end do
    normPairs = sum(ProbAtom(1:nPairAtoms))
    ranNum = normPairs * grnd()
    sumInt = ProbAtom(1)
    nSelPair = 1
    do while (sumInt < ranNum .and. nSelPair < nPairAtoms)
      nSelPair = nSelPair + 1
      sumInt = sumInt + ProbAtom(nSelPair)
    end do
    nAtom = NeighborPairs(nSel, nTarget)%details%pairIndices(2 * nSelPair - 1)
    nTargAtom = NeighborPairs(nSel, nTarget)%details%pairIndices(2 * nSelPair)
!    ProbSel = ProbSel * (1.0_dp / real(nAtoms(nType), dp)) / (1.0_dp / normPairs)
  end if
  deallocate(ProbAtom)
end subroutine LongLinearEBias_Remove_ChooseNeighbor
!================================================================================= 
! Calculates the probability of selecting nTarget in the reverse insertion move for the new
! state (post-deletion of nSel, type nType) in a grand canonical Monte Carlo nucleation
! simulation, with P_TargIn = exp(α * U_i,new) / Σ_j^(n-1) exp(α * U_j,new), where
! U_i,new = ETable(i) - dE(i) - avgE(iType) and α = AlphaTargetIns(nType, iType). Outputs
! ProbTargIn for LLAVBMC_EBias_Rosen_Out’s deletion acceptance probability. Supports
! minDistCriteria for biomolecules and all force fields (LJ_Q, Mpipi, HPS_single,
! HPS_piecewise, HPS_cation_pi).
subroutine LongLinearEBias_Remove_ReverseProbTarget(nTarget, nSel, nType, dE, ProbTargIn)
  use SimParameters, only: NTotal, maxMol, nMolTypes, isActive, beta
  use EnergyTables,  only: ETable, AlphaTargetIns
  use Coords,        only: typeList
  use VarPrecision,  only: dp
  implicit none

  ! Inputs/Outputs
  integer, intent(in) :: nTarget                ! Global index of target molecule
  integer, intent(in) :: nSel                   ! Global index of molecule to delete
  integer, intent(in) :: nType                  ! Type of molecule to delete
  real(dp), intent(in) :: dE(1:maxMol)          ! Energy differences due to deletion
  real(dp), intent(out) :: ProbTargIn           ! Target selection probability in new state

  ! Local variables
  integer :: i, iType                           ! Loop index and molecule type
  integer :: cnt(1:nMolTypes)                   ! Count of molecules per type
  real(dp) :: avgE(1:nMolTypes)                 ! Average energy per type in new state
  real(dp) :: ProbTable(1:maxMol)               ! Unnormalized probabilities
  real(dp) :: norm                              ! Normalization factor

  ! Step 1: Handle single remaining molecule
  ! If NTotal = 2, only one molecule remains after deletion (n-1 = 1)
  if (NTotal == 2) then
    ProbTargIn = 1.0_dp
    return
  end if

  ! Step 2: Initialize arrays
  ProbTable = 0.0_dp
  avgE = 0.0_dp
  cnt = 0

  ! Step 3: Compute type-specific average energies
  ! U_i,new = ETable(i) - dE(i), excluding nSel
  do i = 1, maxMol
    if (isActive(i) .and. i /= nSel) then
      iType = typeList(i)
      avgE(iType) = avgE(iType) + ETable(i) - dE(i)
      cnt(iType) = cnt(iType) + 1
    end if
  end do
  do iType = 1, nMolTypes
    if (cnt(iType) /= 0) then
      avgE(iType) = avgE(iType) / real(cnt(iType), dp)
    end if
  end do

  ! Step 4: Compute unnormalized probabilities
  ! ProbTable(i) = exp(AlphaTargetIns(nType, iType) * beta * (U_i,new - avgE(iType)))
  do i = 1, maxMol
    if (isActive(i) .and. i /= nSel) then
      iType = typeList(i)
      ProbTable(i) = exp(AlphaTargetIns(nType, iType) * beta * (ETable(i) - dE(i) - avgE(iType)))
    end if
  end do

  ! Step 5: Normalize and set target probability
  norm = sum(ProbTable(1:maxMol))
  ProbTargIn = ProbTable(nTarget) / norm
end subroutine LongLinearEBias_Remove_ReverseProbTarget
!================================================================================= 
      end module
