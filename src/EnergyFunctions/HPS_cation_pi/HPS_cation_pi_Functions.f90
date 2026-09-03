!*********************************************************************************************************************
!     This file contains the energy functions that work for Hydrophobicity scale w/ electrostatic screening style forcefields
!     these functions are enclosed inside of the module "InterMolecularEnergy" so that
!     the energy functions can be freely exchanged from the simulation.
!     The prefix naming scheme implies the following:
!           Detailed - Complete energy calculation inteded for use at the beginning and end
!                      of the simulation.  This function is not inteded for use mid-simulation.
!             Shift  - Calculates the energy difference for any move that does not result
!                      in molecules being added or removed from the cluster. This function
!                      receives any number of Displacement vectors from the parent function as input.
!              Mol   - Calculates the energy for a molecule already present in the system. For
!                      use in moves such as particle deletion moves. 
!             NewMol - Calculates the energy for a molecule that has been freshly inserted into the system.
!                      Intended for use in Rosenbluth Sampling, Swap In, etc.
!           Exchange - Combines the Mol and New Mol routines for moves that simultaniously add and remove a particle at the same time.
!*********************************************************************************************************************
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
      module InterEnergy_HPS_cation_pi
      use VarPrecision
      contains
!====================================================================================== 
! Calculates intermolecular energies for the HPS_cation_pi force field in a grand canonical
! ensemble nucleation simulation, optimized for long biomolecules. Computes piecewise
! Lennard-Jones (HPS) with cation-pi interactions and Debye-Hückel electrostatics.
! Updates pairList with distances (distCriteria: first atom pair; minDistCriteria: minimum
! atom pair) or interaction energies (energy criterion). Populates NeighborPairs for
! minDistCriteria, storing atom pair indices for neighbor molecule pairs within Dist_Critr_sq.
! Stops if overlapping atoms are detected (rare, <1% probability). Assumes NeighborPairs
! is allocated with at least 2 * maxAtoms entries per molecule pair.
subroutine Detailed_ECalc_Inter(E_T, pairList)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: atomArray, nAtoms, r_min_tab
  use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, &
                                        sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElecsq
  use Coords, only: MolArray, NeighborPairs
  use SimParameters, only: nMolTypes, NPART, distCriteria, minDistCriteria, &
                           Dist_Critr_sq, kapa
  use EnergyTables, only: ETable, E_Inter_T
  use ParallelVar, only: nout
  implicit none

  ! Input/Output variables
  real(dp), intent(inout) :: E_T              ! Total intermolecular energy
  real(dp), intent(inout) :: pairList(:, :)   ! Pairwise data: distances or energies

  ! Local variables
  integer :: iType, jType, iMol, jMol, iAtom, jAtom, iIndx, jIndx, iPair
  integer :: jMolMin
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_sq, r
  real(dp) :: eps, epsLJ, sig_sq, q, lambda, rcutNBsq, rmin_ij
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS, E_Ele
  logical :: q_Nzero

  ! Initialize energy terms and pair list
  E_HPS = 0.0_dp
  E_Ele = 0.0_dp
  E_Inter_T = 0.0_dp
  pairList = huge(dp)
  if (.not. distCriteria .and. .not. minDistCriteria) pairList = 0.0_dp
  ETable = 0.0_dp

  ! Loop over molecule types and pairs
  do iType = 1, nMolTypes
    do jType = iType, nMolTypes
      ! Set starting jMol to avoid double-counting
      jMolMin = 1
      if (iType == jType) jMolMin = 2

      do iMol = 1, NPART(iType)
        do jMol = jMolMin, NPART(jType)
          if (iType == jType .and. iMol >= jMol) cycle
          iIndx = MolArray(iType)%mol(iMol)%indx
          jIndx = MolArray(jType)%mol(jMol)%indx

          ! Loop over atom pairs
          do iAtom = 1, nAtoms(iType)
            atmType1 = atomArray(iType, iAtom)
            do jAtom = 1, nAtoms(jType)
              atmType2 = atomArray(jType, jAtom)

              ! Retrieve force field parameters
              eps = eps_tab(atmType1, atmType2)
              epsLJ = epsLJ_tab(atmType1, atmType2)
              q_Nzero = q_Nonzero(atmType1, atmType2)
              sig_sq = sigsq_tab(atmType1, atmType2)
              lambda = lambda_tab(atmType1, atmType2)
              rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
              rmin_ij = r_min_tab(atmType1, atmType2)

              ! Calculate squared distance
              rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r_sq = rx * rx + ry * ry + rz * rz

              ! Check for overlap
              if (r_sq < rmin_ij) then
                write(nout, '(A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, F12.6)') &
                  "ERROR! Overlapping atoms detected: iType=", iType, &
                  ", iMol=", iMol, ", iAtom=", iAtom, ", jType=", jType, &
                  ", jMol=", jMol, ", jAtom=", jAtom, ", r_sq=", r_sq
                stop "ERROR! Overlapping atoms found in the current configuration!"
              endif

              ! Update pair list and neighbor data
              if (minDistCriteria) then
                if (r_sq < pairList(iIndx, jIndx)) pairList(iIndx, jIndx) = r_sq
                if (r_sq < pairList(jIndx, iIndx)) pairList(jIndx, iIndx) = r_sq
                if (r_sq < Dist_Critr_sq .and. allocated(NeighborPairs)) then
                  iPair = NeighborPairs(iIndx, jIndx)%details%nPairs + 1
                  NeighborPairs(iIndx, jIndx)%details%nPairs = iPair
                  NeighborPairs(iIndx, jIndx)%details%pairIndices(2 * iPair - 1) = iAtom
                  NeighborPairs(iIndx, jIndx)%details%pairIndices(2 * iPair) = jAtom
                  NeighborPairs(jIndx, iIndx)%details%nPairs = iPair
                  NeighborPairs(jIndx, iIndx)%details%pairIndices(2 * iPair - 1) = jAtom
                  NeighborPairs(jIndx, iIndx)%details%pairIndices(2 * iPair) = iAtom
                endif
              elseif (distCriteria .and. iAtom == 1 .and. jAtom == 1) then
                pairList(iIndx, jIndx) = r_sq
                pairList(jIndx, iIndx) = r_sq
              endif

              ! Calculate HPS energy (piecewise LJ with cation-pi)
              HPS = 0.0_dp
              if (r_sq < rcutNBsq) then
                x = sig_sq / r_sq
                x = x * x * x
                LJ = 4.0_dp * eps * x * (x - 1.0_dp)
                if (x < 0.5_dp) then
                  HPS = lambda * LJ
                else
                  HPS = LJ + eps * (1.0_dp - lambda)
                endif
                if (epsLJ > 1.0E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
              endif
              E_HPS = E_HPS + HPS

              ! Calculate electrostatic energy
              Ele = 0.0_dp
              if (q_Nzero .and. r_sq < rcutElecsq) then
                q = q_tab(atmType1, atmType2)
                r = sqrt(r_sq)
                if (kapa * r < 50.0_dp) then
                  Ele = q * exp(-kapa * r) / r
                endif
              endif
              E_Ele = E_Ele + Ele

              ! Update pair list and energy table for energy criterion
              if (.not. distCriteria .and. .not. minDistCriteria) then
                pairList(iIndx, jIndx) = pairList(iIndx, jIndx) + Ele + HPS
                pairList(jIndx, iIndx) = pairList(iIndx, jIndx)
              endif
              ETable(iIndx) = ETable(iIndx) + Ele + HPS
              ETable(jIndx) = ETable(jIndx) + Ele + HPS
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  ! Output energy contributions
  write(nout, *) "Hydrophobicity scale Energy: ", E_HPS
  write(nout, *) "Electrostatic Energy: ", E_Ele

  ! Update total energies
  E_Inter_T = E_Ele + E_HPS
  E_T = E_T + E_Inter_T
end subroutine Detailed_ECalc_Inter
!======================================================================================      
! Calculates the intermolecular energy change for a molecular displacement in a Grand Canonical
! Monte Carlo nucleation simulation using the HPS_cation_pi force field. Computes the energy
! difference (E_Trial = E_HPS + E_Ele) for displaced atoms interacting with stationary atoms,
! optimized for long linear biomolecules (sizeDisp typically >= 8). Updates PairList with
! distances (distCriteria: first atom pair; minDistCriteria: minimum atom pair) or energies
! (energy criterion), and dETable with energy changes. Populates NeighborDetailsNew with atom
! pair indices for minDistCriteria when provided. Rejects moves immediately if overlaps are
! detected (rare, <1% probability). Assumes NeighborDetailsNew is allocated with at least
! 2 * maxAtoms entries per molecule.
subroutine Shift_ECalc_Inter(E_Trial, disp, PairList, dETable, rejMove, NeighborDetailsNew)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: atomArray, nAtoms, r_min_tab
  use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, &
                                        sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElecsq
  use Coords, only: MolArray
  use CoordinateTypes, only: NeighborDetails, Displacement
  use SimParameters, only: distCriteria, minDistCriteria, Dist_Critr_sq, kapa, &
                           NPART, nMolTypes, maxMol, maxAtoms
  implicit none

  ! Input/Output variables
  real(dp), intent(out) :: E_Trial                    ! Energy change for the trial move
  type(Displacement), intent(in) :: disp(:)           ! Displaced atom data
  real(dp), intent(inout) :: PairList(:)              ! Pairwise data: distances or energies
  real(dp), intent(inout) :: dETable(:)               ! Energy difference table
  logical, intent(out) :: rejMove                     ! Flag to reject move if atoms overlap
  type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)  ! Neighbor atom pairs

  ! Local variables
  integer :: iType, jType, iMol, jMol, iAtom, jAtom, iDisp, iIndx, jIndx, iPair
  integer :: sizeDisp
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_new_sq, r_old_sq, r_new, r_old
  real(dp) :: eps, epsLJ, sig_sq, q, lambda, rcutNBsq, rmin_ij
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS, E_Ele
  logical :: q_Nzero

  ! Initialize energy terms and arrays
  sizeDisp = size(disp)
  E_HPS = 0.0_dp
  E_Ele = 0.0_dp
  E_Trial = 0.0_dp
  PairList = huge(dp)
  if (.not. distCriteria .and. .not. minDistCriteria) PairList = 0.0_dp
  dETable = 0.0_dp
  rejMove = .false.

  ! Get molecule type and index for displaced fragment
  iType = disp(1)%molType
  iMol = disp(1)%molIndx
  iIndx = MolArray(iType)%mol(iMol)%indx

  ! Loop over displaced atoms
  do iDisp = 1, sizeDisp
    iAtom = disp(iDisp)%atmIndx
    atmType1 = atomArray(iType, iAtom)

    ! Loop over stationary molecules and atoms
    do jType = 1, nMolTypes
      do jMol = 1, NPART(jType)
        if (iType == jType .and. iMol == jMol) cycle
        jIndx = MolArray(jType)%mol(jMol)%indx

        do jAtom = 1, nAtoms(jType)
          atmType2 = atomArray(jType, jAtom)

          ! Retrieve force field parameters
          eps = eps_tab(atmType2, atmType1)
          epsLJ = epsLJ_tab(atmType2, atmType1)
          q_Nzero = q_Nonzero(atmType2, atmType1)
          sig_sq = sigsq_tab(atmType2, atmType1)
          lambda = lambda_tab(atmType2, atmType1)
          rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
          rmin_ij = r_min_tab(atmType2, atmType1)

          ! Calculate distance for new position
          rx = disp(iDisp)%x_new - MolArray(jType)%mol(jMol)%x(jAtom)
          ry = disp(iDisp)%y_new - MolArray(jType)%mol(jMol)%y(jAtom)
          rz = disp(iDisp)%z_new - MolArray(jType)%mol(jMol)%z(jAtom)
          r_new_sq = rx * rx + ry * ry + rz * rz

          ! Check for overlap
          if (r_new_sq < rmin_ij) then
            rejMove = .true.
            return
          endif

          ! Update pair list and neighbor data
          if (minDistCriteria) then
            if (r_new_sq < PairList(jIndx)) PairList(jIndx) = r_new_sq
            if (r_new_sq < Dist_Critr_sq .and. present(NeighborDetailsNew)) then
              iPair = NeighborDetailsNew(jIndx)%nPairs + 1
              NeighborDetailsNew(jIndx)%nPairs = iPair
              NeighborDetailsNew(jIndx)%pairIndices(2 * iPair - 1) = iAtom
              NeighborDetailsNew(jIndx)%pairIndices(2 * iPair) = jAtom
            endif
          elseif (distCriteria .and. iAtom == 1 .and. jAtom == 1) then
            PairList(jIndx) = r_new_sq
          endif

          ! Calculate energy changes for displaced atoms
          if (disp(iDisp)%Displaced) then
            ! Calculate distance for old position
            rx = disp(iDisp)%x_old - MolArray(jType)%mol(jMol)%x(jAtom)
            ry = disp(iDisp)%y_old - MolArray(jType)%mol(jMol)%y(jAtom)
            rz = disp(iDisp)%z_old - MolArray(jType)%mol(jMol)%z(jAtom)
            r_old_sq = rx * rx + ry * ry + rz * rz

            ! Calculate HPS energy for new position
            HPS = 0.0_dp
            if (r_new_sq < rcutNBsq) then
              x = sig_sq / r_new_sq
              x = x * x * x
              LJ = 4.0_dp * eps * x * (x - 1.0_dp)
              if (x < 0.5_dp) then
                HPS = lambda * LJ
              else
                HPS = LJ + eps * (1.0_dp - lambda)
              endif
              if (epsLJ > 1.0E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
            endif
            E_HPS = E_HPS + HPS
            if (.not. distCriteria .and. .not. minDistCriteria) then
              PairList(jIndx) = PairList(jIndx) + HPS
            endif
            dETable(iIndx) = dETable(iIndx) + HPS
            dETable(jIndx) = dETable(jIndx) + HPS

            ! Subtract HPS energy for old position
            HPS = 0.0_dp
            if (r_old_sq < rcutNBsq) then
              x = sig_sq / r_old_sq
              x = x * x * x
              LJ = 4.0_dp * eps * x * (x - 1.0_dp)
              if (x < 0.5_dp) then
                HPS = lambda * LJ
              else
                HPS = LJ + eps * (1.0_dp - lambda)
              endif
              if (epsLJ > 1.0E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
            endif
            E_HPS = E_HPS - HPS
            dETable(iIndx) = dETable(iIndx) - HPS
            dETable(jIndx) = dETable(jIndx) - HPS

            ! Calculate electrostatic energy
            if (q_Nzero) then
              q = q_tab(atmType2, atmType1)
              ! New position
              Ele = 0.0_dp
              if (r_new_sq < rcutElecsq) then
                r_new = sqrt(r_new_sq)
                if (kapa * r_new < 50.0_dp) then
                  Ele = q * exp(-kapa * r_new) / r_new
                endif
              endif
              E_Ele = E_Ele + Ele
              if (.not. distCriteria .and. .not. minDistCriteria) then
                PairList(jIndx) = PairList(jIndx) + Ele
              endif
              dETable(iIndx) = dETable(iIndx) + Ele
              dETable(jIndx) = dETable(jIndx) + Ele

              ! Subtract old position
              Ele = 0.0_dp
              if (r_old_sq < rcutElecsq) then
                r_old = sqrt(r_old_sq)
                if (kapa * r_old < 50.0_dp) then
                  Ele = q * exp(-kapa * r_old) / r_old
                endif
              endif
              E_Ele = E_Ele - Ele
              dETable(iIndx) = dETable(iIndx) - Ele
              dETable(jIndx) = dETable(jIndx) - Ele
            endif
          endif
        enddo
      enddo
    enddo
  enddo

  ! Correct pair list for partial displacements
  if (.not. distCriteria .and. .not. minDistCriteria) then
    if (sizeDisp < nAtoms(iType)) then
      call Shift_PairList_Correct(disp, PairList)
    endif
  endif

  ! Compute total energy change
  E_Trial = E_HPS + E_Ele
end subroutine Shift_ECalc_Inter
!======================================================================================
! Corrects the PairList array for partial molecule moves in a Grand Canonical Monte Carlo
! nucleation simulation using the HPS_cation_pi force field. Adds intermolecular energy
! contributions (piecewise Lennard-Jones with cation-pi term and Debye-Hückel electrostatics)
! from undisplaced atoms in the moved molecule to PairList, used for energy-based cluster
! criteria when sizeDisp < nAtoms(iType) and neither distCriteria nor minDistCriteria is active.
! Called by Shift_ECalc_Inter for moves like single-atom displacements in long biomolecules.
pure subroutine Shift_PairList_Correct(disp, PairList)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: atomArray, nAtoms
  use CoordinateTypes, only: Displacement
  use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, &
                                        sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElecsq
  use Coords, only: MolArray
  use SimParameters, only: NPART, nMolTypes, kapa
  implicit none

  ! Input/Output variables
  type(Displacement), intent(in) :: disp(:)     ! Displaced atom data
  real(dp), intent(inout) :: PairList(:)       ! Pairwise energies

  ! Local variables
  integer :: iType, jType, iMol, jMol, iAtom, jAtom, jIndx
  integer :: sizeDisp
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_sq, r
  real(dp) :: eps, epsLJ, sig_sq, q, lambda, rcutNBsq
  real(dp) :: HPS, LJ, Ele, x
  logical :: q_Nzero

  ! Initialize variables
  sizeDisp = size(disp)
  iType = disp(1)%molType
  iMol = disp(1)%molIndx

  ! Loop over undisplaced atoms in the moved molecule
  do iAtom = 1, nAtoms(iType)
    if (any(disp%atmIndx == iAtom)) cycle
    atmType1 = atomArray(iType, iAtom)

    ! Loop over stationary molecules and atoms
    do jType = 1, nMolTypes
      do jMol = 1, NPART(jType)
        if (iType == jType .and. iMol == jMol) cycle
        jIndx = MolArray(jType)%mol(jMol)%indx

        do jAtom = 1, nAtoms(jType)
          atmType2 = atomArray(jType, jAtom)

          ! Retrieve force field parameters
          eps = eps_tab(atmType2, atmType1)
          epsLJ = epsLJ_tab(atmType2, atmType1)
          q_Nzero = q_Nonzero(atmType2, atmType1)
          sig_sq = sigsq_tab(atmType2, atmType1)
          lambda = lambda_tab(atmType2, atmType1)
          rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)

          ! Calculate squared distance
          rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
          ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
          rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
          r_sq = rx * rx + ry * ry + rz * rz

          ! Calculate HPS energy (piecewise LJ with cation-pi)
          HPS = 0.0_dp
          if (r_sq < rcutNBsq) then
            x = sig_sq / r_sq
            x = x * x * x
            LJ = 4.0_dp * eps * x * (x - 1.0_dp)
            if (x < 0.5_dp) then
              HPS = lambda * LJ
            else
              HPS = LJ + eps * (1.0_dp - lambda)
            endif
            if (epsLJ > 1.0E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
          endif
          PairList(jIndx) = PairList(jIndx) + HPS

          ! Calculate electrostatic energy
          Ele = 0.0_dp
          if (q_Nzero .and. r_sq < rcutElecsq) then
            q = q_tab(atmType2, atmType1)
            r = sqrt(r_sq)
            if (kapa * r < 50.0_dp) then
              Ele = q * exp(-kapa * r) / r
            endif
          endif
          PairList(jIndx) = PairList(jIndx) + Ele
        enddo
      enddo
    enddo
  enddo
end subroutine Shift_PairList_Correct
!======================================================================================      
! Computes the intermolecular energy of a specified molecule with all other molecules in a Grand
! Canonical Monte Carlo simulation using the HPS_cation_pi force field (useDistStore=.false.).
! Calculates E_Trial (piecewise Lennard-Jones with cation-pi and Debye-Hückel electrostatics) for
! molecule (iType, iMol), updating dETable for cluster criteria. Distances are computed on-the-fly.
! Optimized for long biomolecules with rare overlaps (<1%). Used in swap-out moves (e.g.,
! SwapOut_ECalc_HPS_cation_pi) and other energy calculations.
pure subroutine Mol_ECalc_Inter(iType, iMol, dETable, E_Trial)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: atomArray, nAtoms
  use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, &
                                        sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElecsq
  use Coords, only: MolArray
  use SimParameters, only: nMolTypes, NPART, kapa
  implicit none

  ! Input/Output variables
  integer, intent(in) :: iType, iMol          ! Molecule type and index
  real(dp), intent(out) :: E_Trial           ! Total intermolecular energy
  real(dp), intent(inout) :: dETable(:)      ! Energy difference table

  ! Local variables
  integer :: iAtom, jType, jMol, jAtom, iIndx, jIndx
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_sq, r
  real(dp) :: eps, epsLJ, sig_sq, q, lambda, rcutNBsq
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS, E_Ele
  logical :: q_Nzero

  ! Initialize energy terms and table
  E_HPS = 0.0_dp
  E_Ele = 0.0_dp
  E_Trial = 0.0_dp
  dETable = 0.0_dp

  ! Get molecule index
  iIndx = MolArray(iType)%mol(iMol)%indx

  ! Loop over atoms in the specified molecule
  do iAtom = 1, nAtoms(iType)
    atmType1 = atomArray(iType, iAtom)

    ! Loop over stationary molecules and atoms
    do jType = 1, nMolTypes
      do jMol = 1, NPART(jType)
        if (iType == jType .and. iMol == jMol) cycle
        jIndx = MolArray(jType)%mol(jMol)%indx

        do jAtom = 1, nAtoms(jType)
          atmType2 = atomArray(jType, jAtom)

          ! Retrieve force field parameters
          eps = eps_tab(atmType2, atmType1)
          epsLJ = epsLJ_tab(atmType2, atmType1)
          q_Nzero = q_Nonzero(atmType2, atmType1)
          sig_sq = sigsq_tab(atmType2, atmType1)
          lambda = lambda_tab(atmType2, atmType1)
          rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)

          ! Calculate squared distance
          rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
          ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
          rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
          r_sq = rx * rx + ry * ry + rz * rz

          ! Skip if distance exceeds cutoffs
          if (r_sq > rcutElecsq) cycle

          ! Calculate HPS energy (piecewise LJ with cation-pi)
          HPS = 0.0_dp
          if (r_sq < rcutNBsq) then
            x = sig_sq / r_sq
            x = x * x * x
            LJ = 4.0_dp * eps * x * (x - 1.0_dp)
            if (x < 0.5_dp) then
              HPS = lambda * LJ
            else
              HPS = LJ + eps * (1.0_dp - lambda)
            endif
            if (epsLJ > 1.0E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
          endif
          E_HPS = E_HPS + HPS
          dETable(iIndx) = dETable(iIndx) + HPS
          dETable(jIndx) = dETable(jIndx) + HPS

          ! Calculate electrostatic energy
          Ele = 0.0_dp
          if (q_Nzero) then
            q = q_tab(atmType2, atmType1)
            r = sqrt(r_sq)
            if (kapa * r < 50.0_dp) then
              Ele = q * exp(-kapa * r) / r
            endif
          endif
          E_Ele = E_Ele + Ele
          dETable(iIndx) = dETable(iIndx) + Ele
          dETable(jIndx) = dETable(jIndx) + Ele
        enddo
      enddo
    enddo
  enddo

  ! Compute total energy
  E_Trial = E_HPS + E_Ele
end subroutine Mol_ECalc_Inter
!======================================================================================      
subroutine NewMol_ECalc_Inter(E_Trial, PairList, dETable, rejMove, NeighborDetailsNew)
  ! Calculates intermolecular energy for a new molecule in a Grand Canonical Monte Carlo
  ! nucleation simulation using the HPS_cation_pi force field. Computes E_Trial (E_HPS + E_Ele),
  ! where E_HPS includes HPS_piecewise and cation-pi Lennard-Jones terms. Updates PairList with
  ! distances (for distCriteria/minDistCriteria) or energies, and populates NeighborDetailsNew
  ! with atom pair indices for minDistCriteria. Rejects moves if atomic overlaps occur.
  ! Supports LLAVBMC and useDistStore=.false.
  ! HPS_cation_pi potential:
  !   u_LJ(r) = 4 * eps * [ (sigma/r)**12 - (sigma/r)**6 ]
  !   U_HPS(r) = u_LJ + eps * (1 - lambda) if (sigma/r)^6 >= 0.5, else lambda * u_LJ
  !   U_CationPi(r) = 4 * epsLJ * [ (sigma/r)**12 - (sigma/r)**6 ] for cation-pi pairs

  use ForceField, only: r_min_tab, nAtoms, atomArray
  use ForceFieldPara_HPS_cation_pi, only: rcutElec, rcutElecsq, eps_tab, lambda_tab, sigsq_tab, q_tab, q_Nonzero, &
                                          cutoffNBsq_tab, epsLJ_tab
  use Coords, only: newMol, MolArray
  use SimParameters, only: nMolTypes, NPART, distCriteria, minDistCriteria, Dist_Critr_sq, kapa
  use CoordinateTypes, only: NeighborDetails
  use VarPrecision
  implicit none

  ! Input/Output variables
  real(dp), intent(out) :: E_Trial              ! Total intermolecular energy (E_HPS + E_Ele)
  real(dp), intent(inout) :: PairList(:)        ! Pairwise interaction list (distances or energies)
  real(dp), intent(inout) :: dETable(:)         ! Energy difference table for cluster criteria
  logical, intent(out) :: rejMove               ! Flag for invalid configurations (e.g., overlaps)
  type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)  ! Neighbor relationships

  ! Local variables
  integer :: iAtom, iIndx, jType, jIndx, jMol, jAtom, iPair
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: eps, sig_sq, q, lambda, epsLJ
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS, rcutNBsq, rmax
  real(dp) :: rmin_ij

  ! Initialize variables
  E_HPS = 0.0_dp
  E_Ele = 0.0_dp      
  E_Trial = 0.0_dp
  dETable = 0.0_dp
  PairList = huge(dp)      
  if ((.not. distCriteria) .and. (.not. minDistCriteria)) PairList = 0.0_dp
  rejMove = .false.
  
  ! Get index for the new molecule
  iIndx = molArray(newMol%molType)%mol(NPART(newMol%molType)+1)%indx

  ! Loop over atoms of the new molecule
  do iAtom = 1, nAtoms(newMol%molType)
    atmType1 = atomArray(newMol%molType, iAtom)
    ! Loop over all molecule types and their atoms
    do jType = 1, nMolTypes
      do jAtom = 1, nAtoms(jType)        
        atmType2 = atomArray(jType, jAtom)
        ! Retrieve force field parameters
        eps = eps_tab(atmType2, atmType1)
        epsLJ = epsLJ_tab(atmType2, atmType1)
        q_Nzero = q_Nonzero(atmType2, atmType1)
        sig_sq = sigsq_tab(atmType2, atmType1)  
        rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
        lambda = lambda_tab(atmType2, atmType1)
        rmin_ij = r_min_tab(atmType2, atmType1)
        ! Loop over existing molecules of type jType
        do jMol = 1, NPART(jType)
          ! Calculate distance components
          rx = abs(newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom))
          ry = abs(newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom))
          rz = abs(newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom))
          rmax = max(rx, ry, rz) 
          ! Check non-bonded cutoff
          if (rmax > rcutElec) cycle
          if (.not. q_Nzero) then
            rmax = rmax * rmax
            if (rmax > rcutNBsq) cycle
          endif
          ! Compute squared distance
          r = rx*rx + ry*ry + rz*rz
          if (r > rcutElecsq) cycle 
          ! Check for atomic overlaps
          if (r < rmin_ij) then
            rejMove = .true.
            return
          endif
          jIndx = molArray(jType)%mol(jMol)%indx  
          ! Update PairList and NeighborDetailsNew for distance criteria
          if (minDistCriteria) then
            if (r < PairList(jIndx)) PairList(jIndx) = r
            if (r < Dist_Critr_sq .and. present(NeighborDetailsNew)) then
              iPair = NeighborDetailsNew(jIndx)%nPairs + 1
              NeighborDetailsNew(jIndx)%nPairs = iPair
              NeighborDetailsNew(jIndx)%pairIndices(2 * iPair - 1) = iAtom
              NeighborDetailsNew(jIndx)%pairIndices(2 * iPair) = jAtom
            endif
          elseif (distCriteria) then             
            if (iAtom == 1) then
              if (jAtom == 1) then
                PairList(jIndx) = r
              endif
            endif
          endif              
          ! Calculate HPS and LJ energies
          HPS = 0.0_dp
          if (r < rcutNBsq) then
            x = sig_sq/r
            x = x * x * x
            LJ = 4.0_dp * eps * x * (x - 1.0_dp)
            if (x < 0.5_dp) then
              HPS = lambda * LJ
            else
              HPS = LJ + eps * (1.0_dp - lambda)
            endif
            if (epsLJ > 1.0E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
          endif
          E_HPS = E_HPS + HPS
          ! Update PairList and dETable for energy-based criteria
          if ((.not. distCriteria) .and. (.not. minDistCriteria)) then                
            PairList(jIndx) = PairList(jIndx) + HPS
          endif
          dETable(jIndx) = dETable(jIndx) + HPS
          dETable(iIndx) = dETable(iIndx) + HPS
          ! Calculate electrostatic energy
          if (q_Nzero) then
            q = q_tab(atmType2, atmType1)
            r = sqrt(r)
            Ele = q * exp(-kapa * r) / r
            E_Ele = E_Ele + Ele
            if ((.not. distCriteria) .and. (.not. minDistCriteria)) then                
              PairList(jIndx) = PairList(jIndx) + Ele
            endif
            dETable(jIndx) = dETable(jIndx) + Ele
            dETable(iIndx) = dETable(iIndx) + Ele
          endif
        enddo
      enddo
    enddo
  enddo
     
  ! Compute total trial energy
  E_Trial = E_HPS + E_Ele
      
end subroutine NewMol_ECalc_Inter
!======================================================================================      
subroutine QuickNei_ECalc_Inter_HPS_cation_pi(jType, jMol, rejMove)
  ! Evaluates intermolecular interaction energy between a new molecule (newMol) and a target
  ! molecule (jType, jMol) for energy-based cluster criterion in a grand canonical ensemble
  ! nucleation simulation. Uses HPS_cation_pi force field:
  !   u_LJ(r) = 4 * eps * [ (sigma/r)^12 - (sigma/r)^6 ]
  !   U_HPS(r) = u_LJ(r) + eps * (1 - lambda), if r <= 2^(1/6) * sigma
  !              lambda * u_LJ(r), if r > 2^(1/6) * sigma
  !   U_cation_pi(r) = 4 * epsLJ * [ (sigma/r)^12 - (sigma/r)^6 ], for specific pairs
  !   U_Electrostatic(r) = (q1 * q2 / r) * exp(-kappa * r)
  ! Computes HPS, cation-π, and electrostatic energies for all atom pairs. Rejects move if:
  !   1. Any pair distance < r_min_tab
  !   2. Total energy (E_Trial) > Eng_Critr
  ! Used in Monte Carlo moves (e.g., LLAVBMC) when distCriteria=.false.

  use VarPrecision, only: dp
  use ForceField, only: nAtoms, atomArray, r_min_tab
  use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, &
                                          cutoffNBsq_tab, rcutElecsq
  use Coords, only: MolArray, newMol
  use SimParameters, only: kapa, Eng_Critr
  implicit none

  ! Input/Output variables
  integer, intent(in) :: jType        ! Target molecule type index
  integer, intent(in) :: jMol         ! Target molecule instance index
  logical, intent(out) :: rejMove     ! Flag to reject the move if criterion fails

  ! Local variables
  integer :: iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: eps, sig_sq, q, lambda, epsLJ
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Trial, E_Ele, E_HPS, rcutNBsq
  real(dp) :: rmin_ij

  ! Initialize variables
  E_HPS = 0.0_dp
  E_Ele = 0.0_dp      
  E_Trial = 0.0_dp
  rejMove = .false.
    
  ! Loop over atoms of the new molecule
  do iAtom = 1, nAtoms(newMol%molType)
    atmType1 = atomArray(newMol%molType, iAtom)
    ! Loop over atoms of the target molecule
    do jAtom = 1, nAtoms(jType)        
      atmType2 = atomArray(jType, jAtom)
      ! Retrieve force field parameters
      rmin_ij = r_min_tab(atmType2, atmType1)
      ! Calculate distance components
      rx = newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
      ry = newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
      rz = newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
      r = rx*rx + ry*ry + rz*rz
      ! Check for atomic overlaps
      if (r < rmin_ij) then
        rejMove = .true.
        return
      endif          
      ! Retrieve additional force field parameters
      eps = eps_tab(atmType2, atmType1)
      epsLJ = epsLJ_tab(atmType2, atmType1)
      q_Nzero = q_Nonzero(atmType2, atmType1)
      sig_sq = sigsq_tab(atmType2, atmType1)  
      rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
      lambda = lambda_tab(atmType2, atmType1) 
      ! Calculate HPS and LJ energies
      HPS = 0.0_dp
      if (r < rcutNBsq) then
        x = sig_sq/r
        x = x * x * x
        LJ = 4.0_dp * eps * x * (x - 1.0_dp)
        if (x < 0.5_dp) then
          HPS = lambda * LJ
        else
          HPS = LJ + eps * (1.0_dp - lambda)
        endif
        if (epsLJ > 1.0E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
      endif
      E_HPS = E_HPS + HPS
      ! Calculate electrostatic energy
      if (q_Nzero) then    
        if (r < rcutElecsq) then
          q = q_tab(atmType1, atmType2)       
          r = sqrt(r)
          Ele = q * exp(-kapa * r) / r
          E_Ele = E_Ele + Ele
        endif
      endif
    enddo
  enddo
     
  ! Compute total trial energy
  E_Trial = E_HPS + E_Ele
  ! Check energy-based cluster criterion
  if (E_Trial > Eng_Critr(newMol%molType, jType)) then
    rejMove = .true.
  endif
      
end subroutine QuickNei_ECalc_Inter_HPS_cation_pi
!======================================================================================
      end module
      
       
