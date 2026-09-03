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
!   HPS_single (HPS-Nucl): Hydropathy Scale Lennard-Jones + Debye-Hückel
!      - U_HPS(r)          = 4 * ε * λ * [ (σ / r)**12 − 2 * (σ / r)**6 ]
!      - U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!****************************************************************************************
      module InterEnergy_HPS_single
      use VarPrecision
      contains
!====================================================================================== 
! Calculates intermolecular energies (hydrophobic-polar scale and electrostatic) for the
! HPS_single force field in a grand canonical ensemble nucleation simulation. Updates
! PairList with distances (distCriteria: first atom pair; minDistCriteria: minimum atom
! pair) or interaction energies (energy criterion). Populates NeighborList and NeighborPairs
! for minDistCriteria, storing atom pair indices for neighbor molecule pairs within
! Dist_Critr_sq.
subroutine Detailed_ECalc_Inter(E_T, pairList)
  use ParallelVar, only: nout
  use ForceField, only: nAtoms, r_min_tab, atomArray
  use ForceFieldPara_HPS_single, only: eps_tab, q_Nonzero, sigsq_tab, lambda_tab, cutoffNBsq_tab, q_tab, rcutElecsq
  use Coords, only: MolArray, NeighborPairs
  use SimParameters, only: nMolTypes, NPART, distCriteria, minDistCriteria, Dist_Critr_sq, kapa
  use EnergyTables, only: ETable, E_Inter_T
  implicit none

  real(dp), intent(inout) :: E_T              ! Total energy (updated with intermolecular terms)
  real(dp), intent(inout) :: pairList(:,:)    ! Pairwise data: distances or energies

  integer :: iType, jType, iMol, jMol, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  integer :: iIndx, jIndx, globIndx1, globIndx2, jMolMin, iPair
  real(dp) :: rx, ry, rz, r
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq, rmin_ij
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS

  ! Initialize energy terms and pair list
  E_HPS = 0E0_dp
  E_Ele = 0E0_dp
  E_Inter_T = 0E0_dp
  pairList = huge(dp)
  if ((.not. distCriteria) .and. (.not. minDistCriteria)) pairList = 0E0_dp
  ETable = 0E0_dp

  ! Loop over molecule types and pairs
  do iType = 1, nMolTypes
    do jType = iType, nMolTypes
      do iMol = 1, NPART(iType)
        if (iType == jType) then
          jMolMin = iMol + 1
        else
          jMolMin = 1
        endif
        do jMol = jMolMin, NPART(jType)
          iIndx = MolArray(iType)%mol(iMol)%indx
          jIndx = MolArray(jType)%mol(jMol)%indx
          ! Loop over atoms in each molecule pair
          do iAtom = 1, nAtoms(iType)
            atmType1 = atomArray(iType, iAtom)
            globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
            do jAtom = 1, nAtoms(jType)
              atmType2 = atomArray(jType, jAtom)
              ! Retrieve force field parameters
              eps = eps_tab(atmType1, atmType2)
              q_Nzero = q_Nonzero(atmType1, atmType2)
              sig_sq = sigsq_tab(atmType1, atmType2)
              lambda = lambda_tab(atmType1, atmType2)
              rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
              rmin_ij = r_min_tab(atmType1, atmType2)
              globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)

              ! Calculate squared distance between atoms
              rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r = rx**2 + ry**2 + rz**2

              ! Update pair list and neighbor data for distance criteria
              if (minDistCriteria) then
                if (r < pairList(iIndx, jIndx)) pairList(iIndx, jIndx) = r
                if (r < pairList(jIndx, iIndx)) pairList(jIndx, iIndx) = r
                if (r < Dist_Critr_sq) then
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
                  pairList(iIndx, jIndx) = r
                  pairList(jIndx, iIndx) = r
                endif
              endif

              ! Check for overlapping atoms
              if (r < rmin_ij) then
                stop "ERROR! Overlapping atoms found in the current configuration!"
              endif

              ! Calculate HPS (hydrophobic-polar scale) energy
              HPS = 0E0_dp
              if (r < rcutNBsq) then
                x = sig_sq / r
                x = x * x * x
                LJ = 4.0_dp * eps * x * (x - 1.0_dp)
                HPS = lambda * LJ
              endif
              E_HPS = E_HPS + HPS

              ! Calculate electrostatic energy
              Ele = 0E0_dp
              if (q_Nzero) then
                if (r < rcutElecsq) then
                  q = q_tab(atmType1, atmType2)
                  r = sqrt(r)
                  Ele = q * exp(-kapa * r) / r
                endif
              endif
              E_Ele = E_Ele + Ele

              ! Update pair list for energy criterion
              if ((.not. distCriteria) .and. (.not. minDistCriteria)) then
                pairList(iIndx, jIndx) = pairList(iIndx, jIndx) + Ele + HPS
                pairList(jIndx, iIndx) = pairList(iIndx, jIndx)
              endif

              ! Update energy table for each molecule
              ETable(iIndx) = ETable(iIndx) + Ele + HPS
              ETable(jIndx) = ETable(jIndx) + Ele + HPS
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  ! Output energy components
  write(nout, *) "Hydrophobicity scale Energy:", E_HPS
  write(nout, *) "Electrostatic Energy:", E_Ele

  ! Update total and intermolecular energies
  E_T = E_T + E_Ele + E_HPS
  E_Inter_T = E_Ele + E_HPS

end subroutine Detailed_ECalc_Inter
!======================================================================================      
! Calculates the intermolecular energy change for a molecular displacement in a Grand Canonical
! Monte Carlo nucleation simulation using the HPS_single force field. Computes the energy difference
! (E_Trial = E_HPS + E_Ele) for displaced atoms interacting with stationary atoms, updates PairList
! (distances for distCriteria/minDistCriteria, energies otherwise), and populates NeighborDetailsNew
! with atom pair indices for minDistCriteria when provided. Supports LLTranslation for long linear
! biomolecules and Translation without NeighborDetailsNew.
subroutine Shift_ECalc_Inter(E_Trial, disp, pairList, dETable, rejMove, NeighborDetailsNew)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: atomArray, nAtoms, r_min_tab
  use ForceFieldPara_HPS_single, only: eps_tab, sigsq_tab, cutoffNBsq_tab, lambda_tab, q_Nonzero, q_tab, rcutElecsq
  use Coords, only: MolArray
  use CoordinateTypes, only: NeighborDetails, Displacement
  use SimParameters, only: distCriteria, minDistCriteria, Dist_Critr_sq, kapa, NPART, nMolTypes, maxMol
  implicit none

  ! Input/Output variables
  real(dp), intent(out) :: E_Trial                      ! Energy change for the trial move
  real(dp), intent(inout) :: pairList(:)                ! Pairwise data: distances or energies
  real(dp), intent(inout) :: dETable(:)                 ! Energy difference table
  logical, intent(out) :: rejMove                       ! Flag to reject move if atoms overlap
  type(Displacement), intent(in) :: disp(:)         ! Displaced atom data
  type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)  ! Neighbor atom pairs for minDistCriteria

  integer :: iType, jType, iMol, jMol, iAtom, jAtom, iDisp
  integer(kind=atomIntType) :: atmType1, atmType2, iIndx, jIndx
  integer :: sizeDisp, iPair
  real(dp) :: rx, ry, rz, r_new, r_old
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq, rmin_ij
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS

  ! Initialize energy terms and pair list
  sizeDisp = size(disp)
  E_HPS = 0E0_dp
  E_Ele = 0E0_dp
  E_Trial = 0E0_dp
  pairList = huge(dp)
  if ((.not. distCriteria) .and. (.not. minDistCriteria)) pairList = 0E0_dp
  dETable = 0E0_dp

  ! Get molecule type and index for displaced fragment
  iType = disp(1)%molType
  iMol = disp(1)%molIndx
  iIndx = MolArray(iType)%mol(iMol)%indx

  ! Calculate intermolecular energy changes for displaced atoms
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
          r_new = rx * rx + ry * ry + rz * rz

          ! Reject move if atoms overlap
          if (r_new < rmin_ij) then
            rejMove = .true.
            return
          endif


          ! Update pair list and neighbor data for distance criteria
          if (minDistCriteria) then
            if (r_new < pairList(jIndx)) pairList(jIndx) = r_new
            if (r_new < Dist_Critr_sq .and. present(NeighborDetailsNew)) then
              iPair = NeighborDetailsNew(jIndx)%nPairs + 1
              NeighborDetailsNew(jIndx)%nPairs = iPair
              NeighborDetailsNew(jIndx)%pairIndices(2 * iPair - 1) = iAtom
              NeighborDetailsNew(jIndx)%pairIndices(2 * iPair) = jAtom
            endif
          elseif (distCriteria) then
            if (iAtom == 1 .and. jAtom == 1) pairList(jIndx) = r_new
          endif

          ! Calculate energy changes for displaced atoms
          if (disp(iDisp)%Displaced) then
            ! Distance for old position
            rx = disp(iDisp)%x_old - MolArray(jType)%mol(jMol)%x(jAtom)
            ry = disp(iDisp)%y_old - MolArray(jType)%mol(jMol)%y(jAtom)
            rz = disp(iDisp)%z_old - MolArray(jType)%mol(jMol)%z(jAtom)
            r_old = rx * rx + ry * ry + rz * rz

            ! Calculate HPS energy for new position
            HPS = 0E0_dp
            if (r_new < rcutNBsq) then
              x = sig_sq / r_new
              x = x * x * x
              LJ = 4.0_dp * eps * x * (x - 1.0_dp)
              HPS = lambda * LJ
            endif
            E_HPS = E_HPS + HPS
            if ((.not. distCriteria) .and. (.not. minDistCriteria)) pairList(jIndx) = pairList(jIndx) + HPS

            dETable(iIndx) = dETable(iIndx) + HPS
            dETable(jIndx) = dETable(jIndx) + HPS

            ! Subtract HPS energy for old position
            HPS = 0E0_dp
            if (r_old < rcutNBsq) then
              x = sig_sq / r_old
              x = x * x * x
              LJ = 4.0_dp * eps * x * (x - 1.0_dp)
              HPS = lambda * LJ
            endif
            E_HPS = E_HPS - HPS
            dETable(iIndx) = dETable(iIndx) - HPS
            dETable(jIndx) = dETable(jIndx) - HPS

            ! Calculate electrostatic energy
            if (q_Nzero) then
              q = q_tab(atmType2, atmType1)
              ! New position
              Ele = 0E0_dp
              if (r_new < rcutElecsq) then
                r_new = sqrt(r_new)
                Ele = q * exp(-kapa * r_new) / r_new
              endif
              E_Ele = E_Ele + Ele
              if ((.not. distCriteria) .and. (.not. minDistCriteria)) pairList(jIndx) = pairList(jIndx) + Ele

              dETable(iIndx) = dETable(iIndx) + Ele
              dETable(jIndx) = dETable(jIndx) + Ele

              ! Subtract old position
              Ele = 0E0_dp
              if (r_old < rcutElecsq) then
                r_old = sqrt(r_old)
                Ele = q * exp(-kapa * r_old) / r_old
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
      call Shift_PairList_Correct(disp, pairList)
    endif
  endif

  ! Compute total energy change
  E_Trial = E_HPS + E_Ele

end subroutine Shift_ECalc_Inter
!======================================================================================
! Corrects the PairList array for partial molecule moves in a Grand Canonical Monte Carlo
! nucleation simulation using the HPS_single force field. Adds intermolecular energy contributions
! (scaled Lennard-Jones and electrostatic) from undisplaced atoms in the moved molecule to PairList,
! used for energy-based cluster criteria when sizeDisp < nAtoms(iType) and neither distCriteria
! nor minDistCriteria is active. Called by Shift_ECalc_Inter for moves like single-atom displacements.
pure subroutine Shift_PairList_Correct(disp, PairList)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: atomArray, nAtoms
  use CoordinateTypes, only: Displacement
  use ForceFieldPara_HPS_single, only: eps_tab, sigsq_tab, cutoffNBsq_tab, lambda_tab, q_Nonzero, q_tab, rcutElecsq
  use Coords, only: MolArray
  use SimParameters, only: kapa, nMolTypes, NPART
  implicit none

  type(Displacement), intent(in) :: disp(:)      ! Displaced atom data
  real(dp), intent(inout) :: PairList(:)        ! Pairwise energy data

  integer :: iType, jType, iMol, jMol, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2, jIndx
  integer :: sizeDisp
  real(dp) :: rx, ry, rz, r
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  logical :: q_Nzero
  real(dp) :: HPS, LJ, x, Ele

  ! Initialize variables
  sizeDisp = size(disp)
  iType = disp(1)%molType
  iMol = disp(1)%molIndx

  ! Loop over undisplaced atoms in the moved molecule
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
          ! Calculate distance between atoms
          rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
          ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
          rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
          r = rx * rx + ry * ry + rz * rz

          ! Calculate HPS (hydrophobic-polar scale) energy
          HPS = 0E0_dp
          if (r < rcutNBsq) then
            x = sig_sq / r
            x = x * x * x
            LJ = 4.0_dp * eps * x * (x - 1.0_dp)
            HPS = lambda * LJ
          endif
          pairList(jIndx) = pairList(jIndx) + HPS

          ! Calculate electrostatic energy
          if (q_Nzero .and. r < rcutElecsq) then
            q = q_tab(atmType2, atmType1)
            r = sqrt(r)
            Ele = q * exp(-kapa * r) / r
            pairList(jIndx) = pairList(jIndx) + Ele
          endif
        enddo
      enddo
    enddo
  enddo

end subroutine Shift_PairList_Correct
!======================================================================================      
! Computes the intermolecular energy of a specified molecule with all other molecules in a Grand Canonical
! Monte Carlo simulation using the HPS_single force field with useDistStore=.false. Calculates E_Trial
! (Hydropathy Scale Lennard-Jones scaled by lambda and Debye-Hückel electrostatics) for molecule (iType, iMol),
! updating dETable for cluster criteria. Distances are computed on-the-fly. Used in swap-out moves
! (e.g., SwapOut_ECalc_HPS_single) and other energy calculations.
pure subroutine Mol_ECalc_Inter(iType, iMol, dETable, E_Trial)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: nAtoms, atomArray
  use ForceFieldPara_HPS_single, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, cutoffNBsq_tab, lambda_tab, rcutElec, rcutElecsq
  use Coords, only: MolArray
  use SimParameters, only: nMolTypes, NPART, kapa
  implicit none

  integer, intent(in) :: iType, iMol             ! Molecule type and index
  real(dp), intent(out) :: E_Trial              ! Total intermolecular energy
  real(dp), intent(inout) :: dETable(:)         ! Energy difference table

  integer :: iAtom, jType, jMol, jAtom, iIndx, jIndx
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r, rmax
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS

  ! Initialize energy terms and table
  E_HPS = 0E0_dp
  E_Ele = 0E0_dp
  E_Trial = 0E0_dp
  dETable = 0E0_dp

  iIndx = MolArray(iType)%mol(iMol)%indx

  ! Loop over atoms in the specified molecule
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
          ! Calculate distance between atoms
          rx = abs(MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom))
          ry = abs(MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom))
          rz = abs(MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom))
          rmax = max(rx, ry, rz)
          if (rmax > rcutElec) cycle
          if (.not. q_Nzero .and. rmax * rmax > rcutNBsq) cycle
          r = rx * rx + ry * ry + rz * rz
          if (r > rcutElecsq) cycle

          ! Calculate HPS (hydrophobic-polar scale) energy
          HPS = 0E0_dp
          if (r < rcutNBsq) then
            x = sig_sq / r
            x = x * x * x
            LJ = 4.0_dp * eps * x * (x - 1.0_dp)
            HPS = lambda * LJ
          endif
          E_HPS = E_HPS + HPS
          dETable(iIndx) = dETable(iIndx) + HPS
          dETable(jIndx) = dETable(jIndx) + HPS

          ! Calculate electrostatic energy
          if (q_Nzero) then
            q = q_tab(atmType2, atmType1)
            r = sqrt(r)
            Ele = q * exp(-kapa * r) / r
            E_Ele = E_Ele + Ele
            dETable(iIndx) = dETable(iIndx) + Ele
            dETable(jIndx) = dETable(jIndx) + Ele
          endif
        enddo
      enddo
    enddo
  enddo

  ! Compute total intermolecular energy
  E_Trial = E_HPS + E_Ele

end subroutine Mol_ECalc_Inter
!======================================================================================      
! Calculates the intermolecular energy for a newly inserted molecule interacting with the existing
! cluster in a Grand Canonical Monte Carlo nucleation simulation using the HPS_single force field.
! Computes E_Trial (E_HPS + E_Ele) for the new molecule, updates PairList (distances for
! distCriteria/minDistCriteria, energies otherwise), and populates NeighborDetailsNew with atom
! pair indices for minDistCriteria when provided. Supports LLAVBMC for long linear biomolecules.
! Rejects the move if atomic overlaps occur or other criteria fail. Compatible with useDistStore=.false.
pure subroutine NewMol_ECalc_Inter(E_Trial, PairList, dETable, rejMove, NeighborDetailsNew)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: r_min_tab, nAtoms, atomArray
  use ForceFieldPara_HPS_single, only: eps_tab, lambda_tab, sigsq_tab, q_tab, q_Nonzero, cutoffNBsq_tab, rcutElec, rcutElecsq
  use Coords, only: newMol, MolArray
  use CoordinateTypes, only: NeighborDetails
  use SimParameters, only: nMolTypes, NPART, distCriteria, minDistCriteria, Dist_Critr_sq, kapa
  implicit none

  logical, intent(out) :: rejMove                    ! Flag to reject move if atoms overlap
  real(dp), intent(out) :: E_Trial                  ! Total intermolecular energy
  real(dp), intent(inout) :: PairList(:)            ! Pairwise data: distances or energies
  real(dp), intent(inout) :: dETable(:)             ! Energy difference table
  type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)  ! Neighbor atom pairs for minDistCriteria

  integer :: iAtom, iIndx, jType, jIndx, jMol, jAtom, iPair
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r, rmax, rmin_ij
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS

  ! Initialize energy terms and arrays
  E_HPS = 0E0_dp
  E_Ele = 0E0_dp
  E_Trial = 0E0_dp
  dETable = 0E0_dp
  PairList = 1E7_dp
  if (.not. distCriteria .and. .not. minDistCriteria) PairList = 0E0_dp
  rejMove = .false.

  ! Get index for the new molecule
  iIndx = MolArray(newMol%molType)%mol(NPART(newMol%molType) + 1)%indx

  ! Loop over atoms in the new molecule
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
          ! Calculate distance between atoms
          rx = abs(newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom))
          ry = abs(newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom))
          rz = abs(newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom))
          rmax = max(rx, ry, rz)
          if (rmax > rcutElec) cycle
          if (.not. q_Nzero .and. rmax * rmax > rcutNBsq) cycle
          r = rx * rx + ry * ry + rz * rz
          if (r > rcutElecsq) cycle

          ! Reject move if atoms overlap
          if (r < rmin_ij) then
            rejMove = .true.
            return
          endif

          jIndx = MolArray(jType)%mol(jMol)%indx
          ! Update pair list and neighbor data for distance criteria
          if (minDistCriteria) then
            if (r < PairList(jIndx)) PairList(jIndx) = r
            if (r < Dist_Critr_sq .and. present(NeighborDetailsNew)) then
              iPair = NeighborDetailsNew(jIndx)%nPairs + 1
              NeighborDetailsNew(jIndx)%nPairs = iPair
              NeighborDetailsNew(jIndx)%pairIndices(2 * iPair - 1) = iAtom
              NeighborDetailsNew(jIndx)%pairIndices(2 * iPair) = jAtom
            endif
          elseif (distCriteria) then
            if (iAtom == 1 .and. jAtom == 1) then
              PairList(jIndx) = r
            endif
          endif

          ! Calculate HPS (hydrophobic-polar scale) energy
          HPS = 0E0_dp
          if (r < rcutNBsq) then
            x = sig_sq / r
            x = x * x * x
            LJ = 4.0_dp * eps * x * (x - 1.0_dp)
            HPS = lambda * LJ
          endif
          E_HPS = E_HPS + HPS
          if (.not. distCriteria .and. .not. minDistCriteria) then
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
            if (.not. distCriteria .and. .not. minDistCriteria) then
              PairList(jIndx) = PairList(jIndx) + Ele
            endif
            dETable(jIndx) = dETable(jIndx) + Ele
            dETable(iIndx) = dETable(iIndx) + Ele
          endif
        enddo
      enddo
    enddo
  enddo

  ! Compute total intermolecular energy
  E_Trial = E_HPS + E_Ele

end subroutine NewMol_ECalc_Inter
!======================================================================================      
! Evaluates the intermolecular interaction energy between a newly inserted molecule (newMol)
! and a target molecule (jType, jMol) to check the energy-based cluster criterion in a grand
! canonical ensemble nucleation simulation. Uses the HPS_single force field:
!   U_HPS(r) = 4 * eps * lambda * [ (sigma/r)^12 - 2 * (sigma/r)^6 ]
!   U_Electrostatic(r) = (q1 * q2 / r) * exp(-kappa * r)
! Computes hydropathy-scaled Lennard-Jones (HPS) and Debye-Hückel electrostatic (Ele)
! contributions for all atom pairs between the molecules. Rejects the move (sets
! rejMove=.true.) if:
!   1. Any pair distance is below the minimum distance (r_min_tab).
!   2. The total interaction energy (E_Trial) exceeds the energy criterion (Eng_Critr).
! Used in Monte Carlo moves (e.g., LLAVBMC) when distCriteria=.false. Does not modify
! NeighborList or NeighborPairs, as it only evaluates energy for cluster checks.
! Supports minimum distance criterion for long linear biomolecules via r_min_tab check.
pure subroutine QuickNei_ECalc_Inter_HPS_single(jType, jMol, rejMove)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: nAtoms, atomArray, r_min_tab
  use ForceFieldPara_HPS_single, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElecsq
  use Coords, only: MolArray, newMol
  use SimParameters, only: kapa, Eng_Critr
  implicit none

  integer, intent(in) :: jType        ! Target molecule type index
  integer, intent(in) :: jMol         ! Target molecule instance index
  logical, intent(out) :: rejMove     ! Flag to reject move if criterion fails

  integer :: iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r, rmin_ij
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Trial, E_Ele, E_HPS

  ! Initialize energy terms and rejection flag
  E_HPS = 0E0_dp
  E_Ele = 0E0_dp
  E_Trial = 0E0_dp
  rejMove = .false.

  ! Loop over atoms in the new molecule and target molecule
  do iAtom = 1, nAtoms(newMol%molType)
    atmType1 = atomArray(newMol%molType, iAtom)
    do jAtom = 1, nAtoms(jType)
      atmType2 = atomArray(jType, jAtom)
      rmin_ij = r_min_tab(atmType2, atmType1)
      ! Calculate distance between atoms
      rx = newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
      ry = newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
      rz = newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
      r = rx * rx + ry * ry + rz * rz

      ! Reject move if atoms overlap
      if (r < rmin_ij) then
        rejMove = .true.
        return
      endif

      ! Retrieve force field parameters
      eps = eps_tab(atmType2, atmType1)
      q_Nzero = q_Nonzero(atmType2, atmType1)
      sig_sq = sigsq_tab(atmType2, atmType1)
      rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
      lambda = lambda_tab(atmType2, atmType1)

      ! Calculate HPS (hydrophobic-polar scale) energy
      HPS = 0E0_dp
      if (r < rcutNBsq) then
        x = sig_sq / r
        x = x * x * x
        LJ = 4.0_dp * eps * x * (x - 1.0_dp)
        HPS = lambda * LJ
      endif
      E_HPS = E_HPS + HPS

      ! Calculate electrostatic energy
      if (q_Nzero .and. r < rcutElecsq) then
        q = q_tab(atmType1, atmType2)
        r = sqrt(r)
        Ele = q * exp(-kapa * r) / r
        E_Ele = E_Ele + Ele
      endif
    enddo
  enddo

  ! Compute total intermolecular energy
  E_Trial = E_HPS + E_Ele

  ! Reject move if energy exceeds criterion
  if (E_Trial > Eng_Critr(newMol%molType, jType)) then
    rejMove = .true.
  endif

end subroutine QuickNei_ECalc_Inter_HPS_single
!======================================================================================
      end module
      
       
