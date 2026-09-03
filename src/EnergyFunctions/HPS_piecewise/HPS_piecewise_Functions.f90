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
!*********************************************************************************************************************
      module InterEnergy_HPS_piecewise
      use VarPrecision
      contains
!====================================================================================== 
! Calculates intermolecular energies (piecewise hydrophobic-polar scale and electrostatic)
! for the HPS_piecewise force field in a grand canonical ensemble nucleation simulation.
! Updates pairList with distances (distCriteria: first atom pair; minDistCriteria: minimum
! atom pair) or interaction energies (energy criterion). Populates NeighborPairs for
! minDistCriteria, storing atom pair indices for neighbor molecule pairs within Dist_Critr_sq.
subroutine Detailed_ECalc_Inter(E_T, pairList)
  use ParallelVar, only: nout
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: atomArray, nAtoms, r_min_tab
  use ForceFieldPara_HPS_piecewise, only: eps_tab, sigsq_tab, cutoffNBsq_tab, lambda_tab, &
                                        q_Nonzero, q_tab, rcutElecsq
  use Coords, only: MolArray, NeighborPairs
  use SimParameters, only: distCriteria, minDistCriteria, Dist_Critr_sq, kapa, NPART, nMolTypes
  use EnergyTables, only: ETable, E_Inter_T
  implicit none

  ! Input/Output variables
  real(dp), intent(inout) :: E_T              ! Total energy (updated with intermolecular terms)
  real(dp), intent(inout) :: pairList(:, :)   ! Pairwise data: distances or energies

  ! Local variables
  integer :: iType, jType, iMol, jMol, iAtom, jAtom, iIndx, jIndx, globIndx1, globIndx2
  integer :: jMolMin, iPair
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_sq, r
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq, rmin_ij
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS_total, E_Ele_total
  logical :: q_Nzero

  ! Initialize energy terms and pair list
  E_HPS_total = 0E0_dp
  E_Ele_total = 0E0_dp
  E_Inter_T = 0E0_dp
  pairList = huge(dp)
  if (.not. distCriteria .and. .not. minDistCriteria) pairList = 0E0_dp
  ETable = 0E0_dp

  ! Calculate intermolecular energies for all molecule pairs
  do iType = 1, nMolTypes
    do iMol = 1, NPART(iType)
      iIndx = MolArray(iType)%mol(iMol)%indx
      do jType = iType, nMolTypes
        jMolMin = merge(iMol + 1, 1, iType == jType)
        do jMol = jMolMin, NPART(jType)
          jIndx = MolArray(jType)%mol(jMol)%indx
          do iAtom = 1, nAtoms(iType)
            atmType1 = atomArray(iType, iAtom)
            globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
            do jAtom = 1, nAtoms(jType)
              atmType2 = atomArray(jType, jAtom)
              globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)

              ! Retrieve force field parameters
              eps = eps_tab(atmType1, atmType2)
              q_Nzero = q_Nonzero(atmType1, atmType2)
              sig_sq = sigsq_tab(atmType1, atmType2)
              lambda = lambda_tab(atmType1, atmType2)
              rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
              rmin_ij = r_min_tab(atmType1, atmType2)

              ! Calculate distance
              rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r_sq = rx * rx + ry * ry + rz * rz

              ! Check for atomic overlap
              if (r_sq < rmin_ij) then
                stop "ERROR: Overlapping atoms found in the current configuration!"
              endif

              ! Update pair list and neighbor data for distance criteria
              if (minDistCriteria) then
                if (r_sq < pairList(iIndx, jIndx)) pairList(iIndx, jIndx) = r_sq
                if (r_sq < pairList(jIndx, iIndx)) pairList(jIndx, iIndx) = r_sq
                if (r_sq < Dist_Critr_sq) then
                  iPair = NeighborPairs(iIndx, jIndx)%details%nPairs + 1
                  NeighborPairs(iIndx, jIndx)%details%nPairs = iPair
                  NeighborPairs(iIndx, jIndx)%details%pairIndices(2 * iPair - 1) = iAtom
                  NeighborPairs(iIndx, jIndx)%details%pairIndices(2 * iPair) = jAtom
                  iPair = NeighborPairs(jIndx, iIndx)%details%nPairs + 1
                  NeighborPairs(jIndx, iIndx)%details%nPairs = iPair
                  NeighborPairs(jIndx, iIndx)%details%pairIndices(2 * iPair - 1) = jAtom
                  NeighborPairs(jIndx, iIndx)%details%pairIndices(2 * iPair) = iAtom
                endif
              elseif (distCriteria) then
                if (iAtom == 1 .and. jAtom == 1) then
                  pairList(iIndx, jIndx) = r_sq
                  pairList(jIndx, iIndx) = r_sq
                endif
              endif

              ! Calculate HPS energy
              HPS = 0E0_dp
              if (r_sq < rcutNBsq) then
                x = sig_sq / r_sq
                x = x * x * x
                LJ = 4.0_dp * eps * x * (x - 1.0_dp)
                if (x < 0.5_dp) then
                  HPS = lambda * LJ
                else
                  HPS = LJ + eps * (1.0_dp - lambda)
                endif
              endif
              E_HPS_total = E_HPS_total + HPS

              ! Calculate electrostatic energy
              Ele = 0E0_dp
              if (q_Nzero .and. r_sq < rcutElecsq) then
                q = q_tab(atmType1, atmType2)
                r = sqrt(r_sq)
                if (kapa * r < 50.0_dp) then  ! Prevent underflow
                  Ele = q * exp(-kapa * r) / r
                endif
              endif
              E_Ele_total = E_Ele_total + Ele

              ! Update energy tables for energy criterion
              if (.not. distCriteria .and. .not. minDistCriteria) then
                pairList(iIndx, jIndx) = pairList(iIndx, jIndx) + HPS + Ele
                pairList(jIndx, iIndx) = pairList(iIndx, jIndx)
              endif
              ETable(iIndx) = ETable(iIndx) + HPS + Ele
              ETable(jIndx) = ETable(jIndx) + HPS + Ele
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  ! Write energy results
  write(nout, '(A,F12.6)') "Hydrophobicity scale Energy:", E_HPS_total
  write(nout, '(A,F12.6)') "Electrostatic Energy:", E_Ele_total

  ! Update total energy
  E_Inter_T = E_HPS_total + E_Ele_total
  E_T = E_T + E_Inter_T

end subroutine Detailed_ECalc_Inter
!======================================================================================      
! Calculates the intermolecular energy change for a molecular displacement in a Grand Canonical
! Monte Carlo nucleation simulation using the HPS_piecewise force field. Computes the energy difference
! (E_Trial = E_HPS + E_Ele) for displaced atoms interacting with stationary atoms, updates pairList
! (distances for distCriteria/minDistCriteria, energies otherwise), and populates NeighborDetailsNew
! with atom pair indices for minDistCriteria when provided. Supports LLTranslation for long linear
! biomolecules and Translation without NeighborDetailsNew.
subroutine Shift_ECalc_Inter(E_Trial, disp, pairList, dETable, rejMove, NeighborDetailsNew)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: atomArray, nAtoms, r_min_tab
  use ForceFieldPara_HPS_piecewise, only: eps_tab, sigsq_tab, cutoffNBsq_tab, lambda_tab, &
                                        q_Nonzero, q_tab, rcutElecsq
  use Coords, only: MolArray
  use CoordinateTypes, only: NeighborDetails, Displacement
  use SimParameters, only: distCriteria, minDistCriteria, Dist_Critr_sq, kapa, NPART, nMolTypes, maxMol
  implicit none

  ! Input/Output variables
  real(dp), intent(out) :: E_Trial                    ! Energy change for the trial move
  real(dp), intent(inout) :: pairList(:)              ! Pairwise data: distances or energies
  real(dp), intent(inout) :: dETable(:)               ! Energy difference table
  logical, intent(out) :: rejMove                     ! Flag to reject move if atoms overlap
  type(Displacement), intent(in) :: disp(:)           ! Displaced atom data
  type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)  ! Neighbor atom pairs for minDistCriteria

  ! Local variables
  integer :: iType, jType, iMol, jMol, iAtom, jAtom, iDisp, iIndx, jIndx, iPair
  integer :: sizeDisp
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_new_sq, r_old_sq, r_new, r_old
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq, rmin_ij
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS_total, E_Ele_total
  logical :: q_Nzero

  ! Initialize energy terms and pair list
  sizeDisp = size(disp)
  E_HPS_total = 0E0_dp
  E_Ele_total = 0E0_dp
  E_Trial = 0E0_dp
  pairList = huge(dp)
  if (.not. distCriteria .and. .not. minDistCriteria) pairList = 0E0_dp
  dETable = 0E0_dp
  rejMove = .false.

  ! Get molecule type and index for displaced fragment
  iType = disp(1)%molType
  iMol = disp(1)%molIndx
  iIndx = MolArray(iType)%mol(iMol)%indx

  ! Calculate energy changes for displaced atoms
  do iDisp = 1, sizeDisp
    iAtom = disp(iDisp)%atmIndx
    atmType1 = atomArray(iType, iAtom)
    do jType = 1, nMolTypes
      do jMol = 1, NPART(jType)
        if (iType == jType .and. iMol == jMol) cycle
        jIndx = MolArray(jType)%mol(jMol)%indx
        do jAtom = 1, nAtoms(jType)
          atmType2 = atomArray(jType, jAtom)

          ! Retrieve force field parameters
          eps = eps_tab(atmType2, atmType1)
          q_Nzero = q_Nonzero(atmType2, atmType1)
          sig_sq = sigsq_tab(atmType2, atmType1)
          rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
          lambda = lambda_tab(atmType2, atmType1)
          rmin_ij = r_min_tab(atmType2, atmType1)

          ! Calculate distance for new position
          rx = disp(iDisp)%x_new - MolArray(jType)%mol(jMol)%x(jAtom)
          ry = disp(iDisp)%y_new - MolArray(jType)%mol(jMol)%y(jAtom)
          rz = disp(iDisp)%z_new - MolArray(jType)%mol(jMol)%z(jAtom)
          r_new_sq = rx * rx + ry * ry + rz * rz

          ! Reject move if atoms overlap
          if (r_new_sq < rmin_ij) then
            rejMove = .true.
            return
          endif

          ! Update pair list and neighbor data for distance criteria
          if (minDistCriteria) then
            if (r_new_sq < pairList(jIndx)) pairList(jIndx) = r_new_sq
            if (r_new_sq < Dist_Critr_sq .and. present(NeighborDetailsNew)) then
              iPair = NeighborDetailsNew(jIndx)%nPairs + 1
              NeighborDetailsNew(jIndx)%nPairs = iPair
              NeighborDetailsNew(jIndx)%pairIndices(2 * iPair - 1) = iAtom
              NeighborDetailsNew(jIndx)%pairIndices(2 * iPair) = jAtom
            endif
          else if (distCriteria .and. iAtom == 1 .and. jAtom == 1) then
            pairList(jIndx) = r_new_sq
          endif

          ! Calculate energy changes for displaced atoms
          if (disp(iDisp)%Displaced) then
            ! Calculate distance for old position
            rx = disp(iDisp)%x_old - MolArray(jType)%mol(jMol)%x(jAtom)
            ry = disp(iDisp)%y_old - MolArray(jType)%mol(jMol)%y(jAtom)
            rz = disp(iDisp)%z_old - MolArray(jType)%mol(jMol)%z(jAtom)
            r_old_sq = rx * rx + ry * ry + rz * rz

            ! Calculate HPS energy for new position
            HPS = 0E0_dp
            if (r_new_sq < rcutNBsq) then
              x = sig_sq / r_new_sq
              x = x * x * x
              LJ = 4.0_dp * eps * x * (x - 1.0_dp)
              if (x < 0.5_dp) then
                HPS = lambda * LJ
              else
                HPS = LJ + eps * (1.0_dp - lambda)
              endif
            endif
            E_HPS_total = E_HPS_total + HPS
            if (.not. distCriteria .and. .not. minDistCriteria) then
              pairList(jIndx) = pairList(jIndx) + HPS
            endif
            dETable(iIndx) = dETable(iIndx) + HPS
            dETable(jIndx) = dETable(jIndx) + HPS

            ! Subtract HPS energy for old position
            HPS = 0E0_dp
            if (r_old_sq < rcutNBsq) then
              x = sig_sq / r_old_sq
              x = x * x * x
              LJ = 4.0_dp * eps * x * (x - 1.0_dp)
              if (x < 0.5_dp) then
                HPS = lambda * LJ
              else
                HPS = LJ + eps * (1.0_dp - lambda)
              endif
            endif
            E_HPS_total = E_HPS_total - HPS
            dETable(iIndx) = dETable(iIndx) - HPS
            dETable(jIndx) = dETable(jIndx) - HPS

            ! Calculate electrostatic energy
            if (q_Nzero) then
              q = q_tab(atmType2, atmType1)
              ! New position
              Ele = 0E0_dp
              if (r_new_sq < rcutElecsq) then
                r_new = sqrt(r_new_sq)
                if (kapa * r_new < 50.0_dp) then  ! Prevent underflow
                  Ele = q * exp(-kapa * r_new) / r_new
                endif
              endif
              E_Ele_total = E_Ele_total + Ele
              if (.not. distCriteria .and. .not. minDistCriteria) then
                pairList(jIndx) = pairList(jIndx) + Ele
              endif
              dETable(iIndx) = dETable(iIndx) + Ele
              dETable(jIndx) = dETable(jIndx) + Ele

              ! Subtract old position
              Ele = 0E0_dp
              if (r_old_sq < rcutElecsq) then
                r_old = sqrt(r_old_sq)
                if (kapa * r_old < 50.0_dp) then  ! Prevent underflow
                  Ele = q * exp(-kapa * r_old) / r_old
                endif
              endif
              E_Ele_total = E_Ele_total - Ele
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
      call Shift_PairList_Correct(disp, pairList)
    endif
  endif

  ! Compute total energy change
  E_Trial = E_HPS_total + E_Ele_total

end subroutine Shift_ECalc_Inter
!======================================================================================
! Corrects the pairList array for partial molecule moves in a Grand Canonical Monte Carlo
! nucleation simulation using the HPS_piecewise force field. Adds intermolecular energy contributions
! (piecewise Lennard-Jones and electrostatic) from undisplaced atoms in the moved molecule to pairList,
! used for energy-based cluster criteria when sizeDisp < nAtoms(iType) and neither distCriteria
! nor minDistCriteria is active. Called by Shift_ECalc_Inter for moves like single-atom displacements.
pure subroutine Shift_PairList_Correct(disp, pairList)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: atomArray, nAtoms
  use ForceFieldPara_HPS_piecewise, only: eps_tab, sigsq_tab, cutoffNBsq_tab, lambda_tab, &
                                        q_Nonzero, q_tab, rcutElecsq
  use Coords, only: MolArray
  use CoordinateTypes, only: Displacement
  use SimParameters, only: kapa, nMolTypes, NPART
  implicit none

  ! Input/Output variables
  type(Displacement), intent(in) :: disp(:)      ! Displaced atom data
  real(dp), intent(inout) :: pairList(:)        ! Pairwise data: energies

  ! Local variables
  integer :: iType, jType, iMol, jMol, iAtom, jAtom, jIndx
  integer :: sizeDisp
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_sq, r
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  real(dp) :: HPS, LJ, Ele, x
  logical :: q_Nzero

  ! Initialize variables
  sizeDisp = size(disp)
  iType = disp(1)%molType
  iMol = disp(1)%molIndx

  ! Calculate energy contributions from undisplaced atoms
  do iAtom = 1, nAtoms(iType)
    if (any(disp%atmIndx == iAtom)) cycle
    atmType1 = atomArray(iType, iAtom)
    do jType = 1, nMolTypes
      do jAtom = 1, nAtoms(jType)
        atmType2 = atomArray(jType, jAtom)
        ! Retrieve force field parameters
        eps = eps_tab(atmType2, atmType1)
        q_Nzero = q_Nonzero(atmType2, atmType1)
        sig_sq = sigsq_tab(atmType2, atmType1)
        rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
        lambda = lambda_tab(atmType2, atmType1)
        do jMol = 1, NPART(jType)
          if (iType == jType .and. iMol == jMol) cycle
          jIndx = MolArray(jType)%mol(jMol)%indx

          ! Calculate distance
          rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
          ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
          rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
          r_sq = rx * rx + ry * ry + rz * rz

          ! Calculate HPS energy
          HPS = 0E0_dp
          if (r_sq < rcutNBsq) then
            x = sig_sq / r_sq
            x = x * x * x
            LJ = 4.0_dp * eps * x * (x - 1.0_dp)
            if (x < 0.5_dp) then
              HPS = lambda * LJ
            else
              HPS = LJ + eps * (1.0_dp - lambda)
            endif
          endif
          pairList(jIndx) = pairList(jIndx) + HPS

          ! Calculate electrostatic energy
          Ele = 0E0_dp
          if (q_Nzero .and. r_sq < rcutElecsq) then
            q = q_tab(atmType2, atmType1)
            r = sqrt(r_sq)
            if (kapa * r < 50.0_dp) then  ! Prevent underflow
              Ele = q * exp(-kapa * r) / r
            endif
            pairList(jIndx) = pairList(jIndx) + Ele
          endif
        enddo
      enddo
    enddo
  enddo

end subroutine Shift_PairList_Correct
!======================================================================================      
! Computes the intermolecular energy of a specified molecule with all other molecules in a Grand Canonical
! Monte Carlo simulation using the HPS_piecewise force field with useDistStore=.false. Calculates E_Trial
! (piecewise Lennard-Jones with switching at (sigma/r)^6 = 0.5 and Debye-Hückel electrostatics) for
! molecule (iType, iMol), updating dETable for cluster criteria. Distances are computed on-the-fly. Used in
! swap-out moves (e.g., SwapOut_ECalc_HPS_piecewise) and other energy calculations.
pure subroutine Mol_ECalc_Inter(iType, iMol, dETable, E_Trial)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: nAtoms, atomArray
  use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, &
                                        cutoffNBsq_tab, lambda_tab, rcutElecsq
  use Coords, only: MolArray
  use SimParameters, only: nMolTypes, NPART, kapa
  implicit none

  ! Input/Output variables
  integer, intent(in) :: iType, iMol                  ! Molecule type and index
  real(dp), intent(out) :: E_Trial                   ! Total intermolecular energy
  real(dp), intent(inout) :: dETable(:)              ! Energy difference table

  ! Local variables
  integer :: iAtom, jType, jAtom, jMol, iIndx, jIndx
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_sq, r
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS_total, E_Ele_total
  logical :: q_Nzero

  ! Initialize energy terms and energy table
  E_HPS_total = 0E0_dp
  E_Ele_total = 0E0_dp
  E_Trial = 0E0_dp
  dETable = 0E0_dp

  ! Get molecule index
  iIndx = MolArray(iType)%mol(iMol)%indx

  ! Calculate intermolecular energies
  do iAtom = 1, nAtoms(iType)
    atmType1 = atomArray(iType, iAtom)
    do jType = 1, nMolTypes
      do jAtom = 1, nAtoms(jType)
        atmType2 = atomArray(jType, jAtom)
        ! Retrieve force field parameters
        eps = eps_tab(atmType2, atmType1)
        q_Nzero = q_Nonzero(atmType2, atmType1)
        sig_sq = sigsq_tab(atmType2, atmType1)
        rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
        lambda = lambda_tab(atmType2, atmType1)
        do jMol = 1, NPART(jType)
          if (iType == jType .and. iMol == jMol) cycle
          jIndx = MolArray(jType)%mol(jMol)%indx

          ! Calculate distance
          rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
          ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
          rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
          r_sq = rx * rx + ry * ry + rz * rz
          if (r_sq > rcutElecsq) cycle

          ! Calculate HPS energy
          HPS = 0E0_dp
          if (r_sq < rcutNBsq) then
            x = sig_sq / r_sq
            x = x * x * x
            LJ = 4.0_dp * eps * x * (x - 1.0_dp)
            if (x < 0.5_dp) then
              HPS = lambda * LJ
            else
              HPS = LJ + eps * (1.0_dp - lambda)
            endif
          endif
          E_HPS_total = E_HPS_total + HPS
          dETable(iIndx) = dETable(iIndx) + HPS
          dETable(jIndx) = dETable(jIndx) + HPS

          ! Calculate electrostatic energy
          Ele = 0E0_dp
          if (q_Nzero) then
            q = q_tab(atmType2, atmType1)
            r = sqrt(r_sq)
            if (kapa * r < 50.0_dp) then  ! Prevent underflow
              Ele = q * exp(-kapa * r) / r
            endif
            E_Ele_total = E_Ele_total + Ele
            dETable(iIndx) = dETable(iIndx) + Ele
            dETable(jIndx) = dETable(jIndx) + Ele
          endif
        enddo
      enddo
    enddo
  enddo

  ! Compute total energy
  E_Trial = E_HPS_total + E_Ele_total

end subroutine Mol_ECalc_Inter
!======================================================================================      
! Calculates the intermolecular energy for a newly inserted molecule (newMol) interacting
! with the existing cluster in a Grand Canonical Monte Carlo nucleation simulation using
! the HPS_piecewise force field. Computes E_Trial (piecewise Lennard-Jones with switching
! at (sigma/r)^6 = 0.5 and Debye-Hückel electrostatics), updates pairList (distances for
! distCriteria/minDistCriteria, energies otherwise), and populates NeighborDetailsNew with
! atom pair indices for minDistCriteria when provided. Rejects the move if atomic overlaps
! occur (r < r_min_tab). Compatible with useDistStore=.false and LLAVBMC for long linear
! biomolecules.
pure subroutine NewMol_ECalc_Inter(E_Trial, pairList, dETable, rejMove, NeighborDetailsNew)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: r_min_tab, nAtoms, atomArray
  use ForceFieldPara_HPS_piecewise, only: rcutElec, rcutElecsq, eps_tab, lambda_tab, &
                                        sigsq_tab, q_tab, q_Nonzero, cutoffNBsq_tab
  use Coords, only: newMol, MolArray
  use CoordinateTypes, only: NeighborDetails
  use SimParameters, only: nMolTypes, NPART, distCriteria, minDistCriteria, Dist_Critr_sq, kapa
  implicit none

  ! Input/Output variables
  real(dp), intent(out) :: E_Trial                    ! Total intermolecular energy (E_HPS + E_Ele)
  real(dp), intent(inout) :: pairList(:)              ! Pairwise interaction list (distances or energies)
  real(dp), intent(inout) :: dETable(:)               ! Energy difference table for cluster criteria
  logical, intent(out) :: rejMove                     ! Flag for invalid configurations (e.g., overlaps)
  type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)  ! Neighbor relationships

  ! Local variables
  integer :: iAtom, iIndx, jType, jIndx, jMol, jAtom, iPair
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r, rmax, rmin_ij
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS_total, E_Ele_total
  logical :: q_Nzero

  ! Initialize energy terms and arrays
  E_HPS_total = 0E0_dp
  E_Ele_total = 0E0_dp
  E_Trial = 0E0_dp
  pairList = 1E7_dp
  if (.not. distCriteria .and. .not. minDistCriteria) pairList = 0E0_dp
  dETable = 0E0_dp
  rejMove = .false.

  ! Get index for the new molecule
  iIndx = MolArray(newMol%molType)%mol(NPART(newMol%molType) + 1)%indx

  ! Calculate intermolecular energies for the new molecule
  do iAtom = 1, nAtoms(newMol%molType)
    atmType1 = atomArray(newMol%molType, iAtom)
    do jType = 1, nMolTypes
      do jAtom = 1, nAtoms(jType)
        atmType2 = atomArray(jType, jAtom)
        ! Retrieve force field parameters
        eps = eps_tab(atmType2, atmType1)
        q_Nzero = q_Nonzero(atmType2, atmType1)
        sig_sq = sigsq_tab(atmType2, atmType1)
        rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
        lambda = lambda_tab(atmType2, atmType1)
        rmin_ij = r_min_tab(atmType2, atmType1)
        do jMol = 1, NPART(jType)
          ! Calculate distance
          rx = abs(newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom))
          ry = abs(newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom))
          rz = abs(newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom))
          rmax = max(rx, ry, rz)
          if (rmax > rcutElec) cycle
          if (.not. q_Nzero) then
            rmax = rmax * rmax
            if (rmax > rcutNBsq) cycle
          endif
          r = rx * rx + ry * ry + rz * rz
          if (r > rcutElecsq) cycle

          ! Check for atomic overlap
          if (r < rmin_ij) then
            rejMove = .true.
            return
          endif

          jIndx = MolArray(jType)%mol(jMol)%indx

          ! Update pair list and neighbor data for distance criteria
          if (minDistCriteria) then
            if (r < pairList(jIndx)) pairList(jIndx) = r
            if (r < Dist_Critr_sq .and. present(NeighborDetailsNew)) then
              iPair = NeighborDetailsNew(jIndx)%nPairs + 1
              NeighborDetailsNew(jIndx)%nPairs = iPair
              NeighborDetailsNew(jIndx)%pairIndices(2 * iPair - 1) = iAtom
              NeighborDetailsNew(jIndx)%pairIndices(2 * iPair) = jAtom
            endif
          elseif (distCriteria) then
            if (iAtom .eq. 1) then
              if (jAtom .eq. 1) then
                pairList(jIndx) = r
              endif
            endif
          endif

          ! Calculate HPS energy
!         if (eps .ne. 0E0_dp) then
          HPS = 0E0_dp
          if (r < rcutNBsq) then
            x = sig_sq / r
            x = x * x * x
            LJ = 4.0_dp * eps * x * (x - 1.0_dp)
            if (x < 0.5_dp) then
              HPS = lambda * LJ
            else
              HPS = LJ + eps * (1.0_dp - lambda)
            endif
          endif
          E_HPS_total = E_HPS_total + HPS
          if (.not. distCriteria .and. .not. minDistCriteria) then
            pairList(jIndx) = pairList(jIndx) + HPS
          endif
          dETable(jIndx) = dETable(jIndx) + HPS
          dETable(iIndx) = dETable(iIndx) + HPS
!         endif

          ! Calculate electrostatic energy
          if (q_Nzero) then
!           if (r < rcutElecsq) then
            q = q_tab(atmType2, atmType1)
            r = sqrt(r)
            Ele = 0E0_dp
            if (kapa * r < 50.0_dp) then  ! Prevent underflow
              Ele = q * exp(-kapa * r) / r
            endif
            E_Ele_total = E_Ele_total + Ele
            if (.not. distCriteria .and. .not. minDistCriteria) then
              pairList(jIndx) = pairList(jIndx) + Ele
            endif
            dETable(jIndx) = dETable(jIndx) + Ele
            dETable(iIndx) = dETable(iIndx) + Ele
!           endif
          endif
        enddo
      enddo
    enddo
  enddo

  ! Compute total energy
  E_Trial = E_HPS_total + E_Ele_total

end subroutine NewMol_ECalc_Inter
!======================================================================================      
! Evaluates the intermolecular interaction energy between a newly inserted molecule (newMol)
! and a target molecule (jType, jMol) to check the energy-based cluster criterion in a Grand
! Canonical Monte Carlo nucleation simulation using the HPS_piecewise force field. Computes
! E_Trial (piecewise Lennard-Jones with switching at (sigma/r)^6 = 0.5 and Debye-Hückel
! electrostatics) for all atom pairs between the molecules. Sets rejMove=.true. if:
!   1. Any pair distance is below r_min_tab (atomic overlap).
!   2. E_Trial exceeds Eng_Critr(newMol%molType, jType).
! Used in Monte Carlo moves (e.g., LLAVBMC) when distCriteria=.false. Does not modify
! NeighborList or NeighborPairs.
pure subroutine QuickNei_ECalc_Inter_HPS_piecewise(jType, jMol, rejMove)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: nAtoms, atomArray, r_min_tab
  use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, &
                                        lambda_tab, cutoffNBsq_tab, rcutElecsq
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
  real(dp) :: rx, ry, rz, r_sq, r
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq, rmin_ij
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS_total, E_Ele_total, E_Trial
  logical :: q_Nzero

  ! Initialize energy terms and rejection flag
  E_HPS_total = 0E0_dp
  E_Ele_total = 0E0_dp
  E_Trial = 0E0_dp
  rejMove = .false.

  ! Calculate intermolecular energies
  do iAtom = 1, nAtoms(newMol%molType)
    atmType1 = atomArray(newMol%molType, iAtom)
    do jAtom = 1, nAtoms(jType)
      atmType2 = atomArray(jType, jAtom)
      ! Retrieve force field parameters
      eps = eps_tab(atmType2, atmType1)
      q_Nzero = q_Nonzero(atmType2, atmType1)
      sig_sq = sigsq_tab(atmType2, atmType1)
      rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
      lambda = lambda_tab(atmType2, atmType1)
      rmin_ij = r_min_tab(atmType2, atmType1)

      ! Calculate distance
      rx = newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
      ry = newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
      rz = newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
      r_sq = rx * rx + ry * ry + rz * rz
      if (r_sq > rcutElecsq) cycle

      ! Check for atomic overlap
      if (r_sq < rmin_ij) then
        rejMove = .true.
        return
      endif

      ! Calculate HPS energy
      HPS = 0E0_dp
      if (r_sq < rcutNBsq) then
        x = sig_sq / r_sq
        x = x * x * x
        LJ = 4.0_dp * eps * x * (x - 1.0_dp)
        if (x < 0.5_dp) then
          HPS = lambda * LJ
        else
          HPS = LJ + eps * (1.0_dp - lambda)
        endif
      endif
      E_HPS_total = E_HPS_total + HPS

      ! Calculate electrostatic energy
      Ele = 0E0_dp
      if (q_Nzero) then
        q = q_tab(atmType2, atmType1)
        r = sqrt(r_sq)
        if (kapa * r < 50.0_dp) then  ! Prevent underflow
          Ele = q * exp(-kapa * r) / r
        endif
        E_Ele_total = E_Ele_total + Ele
      endif
    enddo
  enddo

  ! Compute total energy and check energy criterion
  E_Trial = E_HPS_total + E_Ele_total
  if (E_Trial > Eng_Critr(newMol%molType, jType)) then
    rejMove = .true.
  endif

end subroutine QuickNei_ECalc_Inter_HPS_piecewise
!======================================================================================
      end module
      
       
