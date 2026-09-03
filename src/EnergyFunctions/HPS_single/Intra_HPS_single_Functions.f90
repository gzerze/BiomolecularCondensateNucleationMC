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
!*********************************************************************************************************************
!   HPS_single (HPS-Nucl): Hydropathy Scale Lennard-Jones + Debye-Hückel
!      - U_HPS(r)          = 4 * ε * λ * [ (σ / r)**12 − 2 * (σ / r)**6 ]
!      - U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!****************************************************************************************
      module IntraEnergy_HPS_single
      use VarPrecision
      contains
!======================================================================================      
! Calculates intramolecular non-bonded energies for the HPS_single force field in a
! grand canonical ensemble nucleation simulation. Iterates over specified non-bonded
! atom pairs within each molecule, computing hydrophobic-polar and screened electrostatic
! interactions with cutoffs. Checks for atomic overlaps. Updates total system energy
! (E_T) and non-bonded energy (E_NBond_T). Used in Detailed_ECalc to compute system
! energy. Does not use NeighborList or NeighborPairs, as it handles intramolecular
! interactions, not intermolecular clustering.
subroutine Detailed_ECalc_IntraNonBonded(E_T)
  use VarPrecision, only: dp, atomIntType
  use ParallelVar, only: nout
  use ForceField, only: atomArray, nonBondArray, r_min_tab, nIntraNonBond
  use ForceFieldPara_HPS_single, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
  use Coords, only: MolArray
  use SimParameters, only: nMolTypes, NPART, kapa
  use EnergyTables, only: E_NBond_T
  implicit none

  real(dp), intent(inout) :: E_T          ! Total system energy

  integer :: iType, iMol, iPair, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r, rmax, rmin_ij
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS

  ! Initialize energy terms
  E_HPS = 0E0_dp
  E_Ele = 0E0_dp

  ! Loop over molecule types, instances, and non-bonded atom pairs
  do iType = 1, nMolTypes
    do iMol = 1, NPART(iType)
      do iPair = 1, nIntraNonBond(iType)
        iAtom = nonBondArray(iType, iPair)%nonMembr(1)
        jAtom = nonBondArray(iType, iPair)%nonMembr(2)
        atmType1 = atomArray(iType, iAtom)
        atmType2 = atomArray(iType, jAtom)
        ! Retrieve force field parameters
        eps = eps_tab(atmType1, atmType2)
        sig_sq = sigsq_tab(atmType1, atmType2)
        rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
        q_Nzero = q_Nonzero(atmType1, atmType2)
        lambda = lambda_tab(atmType1, atmType2)
        rmin_ij = r_min_tab(atmType1, atmType2)
        ! Calculate distance between atoms
        rx = abs(MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(iType)%mol(iMol)%x(jAtom))
        ry = abs(MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(iType)%mol(iMol)%y(jAtom))
        rz = abs(MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(iType)%mol(iMol)%z(jAtom))
        rmax = max(rx, ry, rz)
        if (rmax > rcutElec) cycle
        if (.not. q_Nzero .and. rmax * rmax > rcutNBsq) cycle
        r = rx * rx + ry * ry + rz * rz
        if (r > rcutElecsq) cycle

        ! Check for atomic overlaps
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
          q = q_tab(atmType1, atmType2)
          r = sqrt(r)
          Ele = q * exp(-kapa * r) / r
        endif
        E_Ele = E_Ele + Ele
      enddo
    enddo
  enddo

  ! Output energy components
  write(nout, *) "Intra Hydrophobicity scale Energy:", E_HPS
  write(nout, *) "Intra Electrostatic Energy:", E_Ele

  ! Update total and non-bonded energies
  E_T = E_T + E_Ele + E_HPS
  E_NBond_T = E_Ele + E_HPS

end subroutine Detailed_ECalc_IntraNonBonded
!======================================================================================  
! Calculates intramolecular non-bonded energies for a single molecule of each type in
! the gas phase using the HPS_single force field in a grand canonical ensemble nucleation
! simulation. Iterates over specified non-bonded atom pairs, computing hydrophobic-polar
! and screened electrostatic interactions with cutoffs. Checks for atomic overlaps.
! Stores results in E_gasIntra(iType). Does not use NeighborList or NeighborPairs, as it
! handles intramolecular interactions for one molecule per type.
subroutine Detailed_ECalc_GasIntra
  use VarPrecision, only: dp, atomIntType
  use ParallelVar, only: nout
  use ForceField, only: atomArray, nonBondArray, r_min_tab, nIntraNonBond
  use ForceFieldPara_HPS_single, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
  use Coords, only: gasConfig
  use SimParameters, only: nMolTypes, kapa
  use EnergyTables, only: E_gasIntra
  implicit none

  integer :: iType, iPair, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r, rmax, rmin_ij
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS

  ! Loop over molecule types
  do iType = 1, nMolTypes
    ! Initialize energy terms
    E_HPS = 0E0_dp
    E_Ele = 0E0_dp

    ! Loop over non-bonded atom pairs
    do iPair = 1, nIntraNonBond(iType)
      iAtom = nonBondArray(iType, iPair)%nonMembr(1)
      jAtom = nonBondArray(iType, iPair)%nonMembr(2)
      atmType1 = atomArray(iType, iAtom)
      atmType2 = atomArray(iType, jAtom)
      ! Retrieve force field parameters
      eps = eps_tab(atmType1, atmType2)
      sig_sq = sigsq_tab(atmType1, atmType2)
      rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
      q_Nzero = q_Nonzero(atmType1, atmType2)
      lambda = lambda_tab(atmType1, atmType2)
      rmin_ij = r_min_tab(atmType1, atmType2)
      ! Calculate distance between atoms
      rx = abs(gasConfig(iType)%x(iAtom) - gasConfig(iType)%x(jAtom))
      ry = abs(gasConfig(iType)%y(iAtom) - gasConfig(iType)%y(jAtom))
      rz = abs(gasConfig(iType)%z(iAtom) - gasConfig(iType)%z(jAtom))
      rmax = max(rx, ry, rz)
      if (rmax > rcutElec) cycle
      if (.not. q_Nzero .and. rmax * rmax > rcutNBsq) cycle
      r = rx * rx + ry * ry + rz * rz
      if (r > rcutElecsq) cycle

      ! Check for atomic overlaps
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
        q = q_tab(atmType1, atmType2)
        r = sqrt(r)
        Ele = q * exp(-kapa * r) / r
      endif
      E_Ele = E_Ele + Ele
    enddo

    ! Store total intramolecular energy for the molecule type
    E_gasIntra(iType) = E_Ele + E_HPS
    write(nout, *) "Type: ", iType, " Gas Intra Hydrophobicity scale Energy:", E_HPS
    write(nout, *) "Type: ", iType, " Gas Intra Electrostatic Energy:", E_Ele
  enddo

end subroutine Detailed_ECalc_GasIntra
!======================================================================================   
! Calculates the intramolecular non-bonded energy change for a molecular displacement in a Grand Canonical
! Monte Carlo nucleation simulation using the HPS_single force field. Computes the energy difference
! (E_Trial = E_HPS + E_Ele) for HPS (lambda-scaled Lennard-Jones) and electrostatic interactions between
! atom pairs listed in nonBondArray when at least one atom is displaced.
pure subroutine Shift_ECalc_IntraNonBonded(E_Trial, disp)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: atomArray, nonBondArray, nIntraNonBond
  use CoordinateTypes, only: Displacement
  use ForceFieldPara_HPS_single, only: eps_tab, sigsq_tab, cutoffNBsq_tab, lambda_tab, q_Nonzero, q_tab, rcutElecsq
  use Coords, only: MolArray
  use SimParameters, only: maxAtoms, kapa
  implicit none

  real(dp), intent(inout) :: E_Trial              ! Energy change for the trial move
  type(Displacement), intent(in) :: disp(:)       ! Displaced atom data

  logical :: changed(1:maxAtoms)
  integer :: dispIndx(1:maxAtoms)
  integer :: i, nDisp, iType, iMol, iPair, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: x1_New, y1_New, z1_New
  real(dp) :: x1_Old, y1_Old, z1_Old
  real(dp) :: x2_New, y2_New, z2_New
  real(dp) :: x2_Old, y2_Old, z2_Old
  real(dp) :: E_Ele, E_HPS

  ! Initialize displacement tracking and energy terms
  nDisp = size(disp)
  iType = disp(1)%molType
  iMol = disp(1)%molIndx
  changed = .false.
  dispIndx = 0
  E_HPS = 0E0_dp
  E_Ele = 0E0_dp

  ! Mark displaced atoms
  do i = 1, nDisp
    changed(disp(i)%atmIndx) = disp(i)%Displaced
    dispIndx(disp(i)%atmIndx) = i
  enddo

  ! Loop over non-bonded atom pairs
  do iPair = 1, nIntraNonBond(iType)
    iAtom = nonBondArray(iType, iPair)%nonMembr(1)
    jAtom = nonBondArray(iType, iPair)%nonMembr(2)
    if (changed(iAtom) .or. changed(jAtom)) then
      ! Retrieve force field parameters
      atmType1 = atomArray(iType, iAtom)
      atmType2 = atomArray(iType, jAtom)
      eps = eps_tab(atmType1, atmType2)
      sig_sq = sigsq_tab(atmType1, atmType2)
      rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
      q_Nzero = q_Nonzero(atmType1, atmType2)
      lambda = lambda_tab(atmType1, atmType2)

      ! Get coordinates for new and old positions
      x1_Old = merge(disp(dispIndx(iAtom))%x_old, MolArray(iType)%mol(iMol)%x(iAtom), changed(iAtom))
      y1_Old = merge(disp(dispIndx(iAtom))%y_old, MolArray(iType)%mol(iMol)%y(iAtom), changed(iAtom))
      z1_Old = merge(disp(dispIndx(iAtom))%z_old, MolArray(iType)%mol(iMol)%z(iAtom), changed(iAtom))
      x1_New = merge(disp(dispIndx(iAtom))%x_new, MolArray(iType)%mol(iMol)%x(iAtom), changed(iAtom))
      y1_New = merge(disp(dispIndx(iAtom))%y_new, MolArray(iType)%mol(iMol)%y(iAtom), changed(iAtom))
      z1_New = merge(disp(dispIndx(iAtom))%z_new, MolArray(iType)%mol(iMol)%z(iAtom), changed(iAtom))

      x2_Old = merge(disp(dispIndx(jAtom))%x_old, MolArray(iType)%mol(iMol)%x(jAtom), changed(jAtom))
      y2_Old = merge(disp(dispIndx(jAtom))%y_old, MolArray(iType)%mol(iMol)%y(jAtom), changed(jAtom))
      z2_Old = merge(disp(dispIndx(jAtom))%z_old, MolArray(iType)%mol(iMol)%z(jAtom), changed(jAtom))
      x2_New = merge(disp(dispIndx(jAtom))%x_new, MolArray(iType)%mol(iMol)%x(jAtom), changed(jAtom))
      y2_New = merge(disp(dispIndx(jAtom))%y_new, MolArray(iType)%mol(iMol)%y(jAtom), changed(jAtom))
      z2_New = merge(disp(dispIndx(jAtom))%z_new, MolArray(iType)%mol(iMol)%z(jAtom), changed(jAtom))

      ! Calculate new position energy
      rx = x2_New - x1_New
      ry = y2_New - y1_New
      rz = z2_New - z1_New
      r = rx * rx + ry * ry + rz * rz
      HPS = 0E0_dp
      if (r < rcutNBsq) then
        x = sig_sq / r
        x = x * x * x
        LJ = 4.0_dp * eps * x * (x - 1.0_dp)
        LJ = 4.0_dp * eps * x * (x - 1.0_dp)
        HPS = lambda * LJ
      endif
      E_HPS = E_HPS + HPS
      Ele = 0E0_dp
      if (q_Nzero .and. r < rcutElecsq) then
        q = q_tab(atmType1, atmType2)
        r = sqrt(r)
        Ele = q * exp(-kapa * r) / r
        E_Ele = E_Ele + Ele
      endif

      ! Subtract old position energy
      rx = x2_Old - x1_Old
      ry = y2_Old - y1_Old
      rz = z2_Old - z1_Old
      r = rx * rx + ry * ry + rz * rz
      HPS = 0E0_dp
      if (r < rcutNBsq) then
        x = sig_sq / r
        x = x * x * x
        LJ = 4.0_dp * eps * x * (x - 1.0_dp)
        HPS = lambda * LJ
      endif
      E_HPS = E_HPS - HPS
      if (q_Nzero .and. r < rcutElecsq) then
        q = q_tab(atmType1, atmType2)
        r = sqrt(r)
        Ele = q * exp(-kapa * r) / r
        E_Ele = E_Ele - Ele
      endif
    endif
  enddo

  ! Compute total energy change
  E_Trial = E_HPS + E_Ele

end subroutine Shift_ECalc_IntraNonBonded
!======================================================================================
! Computes the intramolecular non-bonded energy of a specified molecule in a Grand Canonical Monte Carlo
! simulation using the HPS_single force field. Calculates E_Trial (Hydropathy Scale Lennard-Jones scaled by
! lambda and Debye-Hückel electrostatics) for non-bonded atom pairs within molecule (iType, iMol). Used in
! swap-out moves (e.g., SwapOut_ECalc_HPS_single) and other energy calculations for flexible molecules.
pure subroutine Mol_ECalc_IntraNonBonded(iType, iMol, E_Trial)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: atomArray, nonBondArray, nIntraNonBond
  use ForceFieldPara_HPS_single, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, cutoffNBsq_tab, lambda_tab, rcutElec, rcutElecsq
  use Coords, only: MolArray
  use SimParameters, only: kapa
  implicit none

  integer, intent(in) :: iType, iMol     ! Molecule type and index
  real(dp), intent(out) :: E_Trial       ! Total intramolecular non-bonded energy

  integer :: iPair, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r, rmax
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS

  ! Initialize energy terms
  E_HPS = 0E0_dp
  E_Ele = 0E0_dp
  E_Trial = 0E0_dp

  ! Loop over non-bonded atom pairs
  do iPair = 1, nIntraNonBond(iType)
    iAtom = nonBondArray(iType, iPair)%nonMembr(1)
    jAtom = nonBondArray(iType, iPair)%nonMembr(2)
    ! Calculate distance between atoms
    rx = abs(MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(iType)%mol(iMol)%x(jAtom))
    ry = abs(MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(iType)%mol(iMol)%y(jAtom))
    rz = abs(MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(iType)%mol(iMol)%z(jAtom))
    rmax = max(rx, ry, rz)
    if (rmax > rcutElec) cycle

    ! Retrieve force field parameters
    atmType1 = atomArray(iType, iAtom)
    atmType2 = atomArray(iType, jAtom)
    eps = eps_tab(atmType1, atmType2)
    sig_sq = sigsq_tab(atmType1, atmType2)
    rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
    q_Nzero = q_Nonzero(atmType1, atmType2)
    lambda = lambda_tab(atmType1, atmType2)
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

    ! Calculate electrostatic energy
    if (q_Nzero) then
      q = q_tab(atmType1, atmType2)
      r = sqrt(r)
      Ele = q * exp(-kapa * r) / r
      E_Ele = E_Ele + Ele
    endif
  enddo

  ! Compute total intramolecular energy
  E_Trial = E_HPS + E_Ele

end subroutine Mol_ECalc_IntraNonBonded
!======================================================================================      
! Calculates the intramolecular non-bonded energy for a newly inserted molecule in a Grand Canonical
! Monte Carlo nucleation simulation using the HPS_single force field. Computes E_Trial (E_HPS + E_ELE)
! for non-bonded atom pairs within the molecule, as specified by nonBondArray, using the HPS
! Lennard-Jones potential and Debye-Hückel electrostatics. Applies non-bonded cutoff (rcutNBsq) and
! cutoff for electrostatics (rcutElecsq). Does not perform overlap checks. Used in swap-in moves (e.g., AVBMC).
pure subroutine NewMol_ECalc_IntraNonBonded(E_Trial)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: atomArray, nonBondArray, nIntraNonBond
  use ForceFieldPara_HPS_single, only: eps_tab, lambda_tab, sigsq_tab, q_Nonzero, q_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
  use Coords, only: newMol
  use SimParameters, only: kapa
  implicit none

  real(dp), intent(out) :: E_Trial       ! Total intramolecular non-bonded energy

  integer :: iType, iPair, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r, rmax
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS

  ! Initialize energy terms
  E_HPS = 0E0_dp
  E_Ele = 0E0_dp
  E_Trial = 0E0_dp
  iType = newMol%molType

  ! Loop over non-bonded atom pairs
  do iPair = 1, nIntraNonBond(iType)
    iAtom = nonBondArray(iType, iPair)%nonMembr(1)
    jAtom = nonBondArray(iType, iPair)%nonMembr(2)
    ! Calculate distance between atoms
    rx = abs(newMol%x(iAtom) - newMol%x(jAtom))
    ry = abs(newMol%y(iAtom) - newMol%y(jAtom))
    rz = abs(newMol%z(iAtom) - newMol%z(jAtom))
    rmax = max(rx, ry, rz)
    if (rmax > rcutElec) cycle

    ! Retrieve force field parameters
    atmType1 = atomArray(iType, iAtom)
    atmType2 = atomArray(iType, jAtom)
    eps = eps_tab(atmType1, atmType2)
    sig_sq = sigsq_tab(atmType1, atmType2)
    rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
    q_Nzero = q_Nonzero(atmType1, atmType2)
    lambda = lambda_tab(atmType1, atmType2)
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

    ! Calculate electrostatic energy
    if (q_Nzero) then
      q = q_tab(atmType1, atmType2)
      r = sqrt(r)
      Ele = q * exp(-kapa * r) / r
      E_Ele = E_Ele + Ele
    endif
  enddo

  ! Compute total intramolecular energy
  E_Trial = E_HPS + E_Ele

end subroutine NewMol_ECalc_IntraNonBonded
!======================================================================================   
      end module
      
      
