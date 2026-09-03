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
!   HPS_piecewise (HPS-KR, KH, FB, TSCL-M2, Urry: Piecewise HPS LJ + Debye-Hückel
!      - u_LJ(r)           = 4 * ε * [ (σ / r)**12 − (σ / r)**6 ]
!      - U_HPS(r) =
!            u_LJ(r) + ε * (1 − λ)       , if r <= 2**(1/6) * σ
!            λ * u_LJ(r)                 , if r >  2**(1/6) * σ
!      - U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!*********************************************************************************************************************
      module IntraEnergy_HPS_piecewise
      use VarPrecision
      contains
!======================================================================================      
! Calculates intramolecular non-bonded energies for the HPS_piecewise force field in a
! Grand Canonical Monte Carlo nucleation simulation. Iterates over specified non-bonded
! atom pairs within each molecule, computing piecewise hydrophobic-polar (HPS) and screened
! electrostatic (Ele) interactions with cutoffs. Checks for atomic overlaps, stopping if
! found. Updates total system energy (E_T) and non-bonded energy (E_NBond_T). Used in
! Detailed_ECalc to compute system energy.
subroutine Detailed_ECalc_IntraNonBonded(E_T)
  use VarPrecision, only: dp, atomIntType
  use ParallelVar, only: nout
  use ForceField, only: atomArray, nonBondArray, r_min_tab, nIntraNonBond
  use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, &
                                        lambda_tab, cutoffNBsq_tab, rcutElecsq
  use Coords, only: MolArray
  use SimParameters, only: nMolTypes, NPART, kapa
  use EnergyTables, only: E_NBond_T
  implicit none

  ! Input/Output variables
  real(dp), intent(inout) :: E_T                   ! Total system energy

  ! Local variables
  integer :: iType, iMol, iPair, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_sq, r, rmin_ij
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS_total, E_Ele_total
  logical :: q_Nzero

  ! Initialize energy terms
  E_HPS_total = 0E0_dp
  E_Ele_total = 0E0_dp

  ! Calculate intramolecular non-bonded energies
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

        ! Calculate distance
        rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(iType)%mol(iMol)%x(jAtom)
        ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(iType)%mol(iMol)%y(jAtom)
        rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(iType)%mol(iMol)%z(jAtom)
        r_sq = rx * rx + ry * ry + rz * rz
        if (r_sq > rcutElecsq) cycle

        ! Check for atomic overlap
        if (r_sq < rmin_ij) then
          stop "ERROR: Overlapping atoms found in the current configuration!"
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
          q = q_tab(atmType1, atmType2)
          r = sqrt(r_sq)
          if (kapa * r < 50.0_dp) then  ! Prevent underflow
            Ele = q * exp(-kapa * r) / r
          endif
        endif
        E_Ele_total = E_Ele_total + Ele
      enddo
    enddo
  enddo

  ! Update energy totals and output
  E_NBond_T = E_HPS_total + E_Ele_total
  E_T = E_T + E_NBond_T
  write(nout, '(A,F12.6)') "Intra Hydrophobicity Scale Energy:", E_HPS_total
  write(nout, '(A,F12.6)') "Intra Electrostatic Energy:", E_Ele_total

end subroutine Detailed_ECalc_IntraNonBonded
!======================================================================================  
! Calculates intramolecular non-bonded energies for a single molecule of each type in
! the gas phase using the HPS_piecewise force field in a Grand Canonical Monte Carlo
! nucleation simulation. Iterates over specified non-bonded atom pairs, computing
! piecewise hydrophobic-polar (HPS) and screened electrostatic (Ele) interactions with
! cutoffs. Checks for atomic overlaps, stopping if found. Stores results in E_gasIntra(iType).
! Does not use NeighborList or NeighborPairs, as it handles intramolecular interactions.
subroutine Detailed_ECalc_GasIntra
  use VarPrecision, only: dp, atomIntType
  use ParallelVar, only: nout
  use ForceField, only: atomArray, nonBondArray, nIntraNonBond, r_min_tab
  use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, &
                                        lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
  use Coords, only: gasConfig
  use SimParameters, only: nMolTypes, kapa
  use EnergyTables, only: E_gasIntra
  implicit none

  ! Local variables
  integer :: iType, iPair, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_sq, r, rmin_ij, rmax
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS_total, E_Ele_total
  logical :: q_Nzero

  ! Calculate intramolecular non-bonded energies for each molecule type
  do iType = 1, nMolTypes
    E_HPS_total = 0E0_dp
    E_Ele_total = 0E0_dp
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

      ! Calculate distance
      rx = abs(gasConfig(iType)%x(iAtom) - gasConfig(iType)%x(jAtom))
      ry = abs(gasConfig(iType)%y(iAtom) - gasConfig(iType)%y(jAtom))
      rz = abs(gasConfig(iType)%z(iAtom) - gasConfig(iType)%z(jAtom))
      rmax = max(rx, ry, rz)
      if (rmax > rcutElec) cycle
      if (.not. q_Nzero) then
        rmax = rmax * rmax
        if (rmax > rcutNBsq) cycle
      endif

      r_sq = rx * rx + ry * ry + rz * rz
      if (r_sq > rcutElecsq) cycle

      ! Check for atomic overlap
      if (r_sq < rmin_ij) then
        stop "ERROR: Overlapping atoms found in the current configuration!"
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
        q = q_tab(atmType1, atmType2)
        r = sqrt(r_sq)
        if (kapa * r < 50.0_dp) then  ! Prevent underflow
          Ele = q * exp(-kapa * r) / r
        endif
      endif
      E_Ele_total = E_Ele_total + Ele
    enddo

    ! Store and output energy for the molecule type
    E_gasIntra(iType) = E_HPS_total + E_Ele_total
    write(nout, *) "Type: ", iType, " Gas Intra Hydrophobicity Scale Energy:", E_HPS_total
    write(nout, *) "Type: ", iType, " Gas Intra Electrostatic Energy:", E_Ele_total
  enddo

end subroutine Detailed_ECalc_GasIntra
!======================================================================================   
! Calculates the intramolecular non-bonded energy change for a molecular displacement in a
! Grand Canonical Monte Carlo nucleation simulation using the HPS_piecewise force field.
! Computes the energy difference (E_Trial = E_HPS + E_Ele) for piecewise Lennard-Jones (HPS)
! and screened electrostatic interactions for non-bonded atom pairs listed in nonBondArray
! where at least one atom is displaced. Used in Monte Carlo displacement moves.
pure subroutine Shift_ECalc_IntraNonBonded(E_Trial, disp)
  use VarPrecision, only: dp, atomIntType
  use SimParameters, only: kapa, maxAtoms
  use ForceField, only: atomArray, nonBondArray, nIntraNonBond
  use CoordinateTypes, only: Displacement
  use ForceFieldPara_HPS_piecewise, only: rcutElecsq, eps_tab, sigsq_tab, cutoffNBsq_tab, &
                                        lambda_tab, q_Nonzero, q_tab
  implicit none

  ! Input/Output variables
  real(dp), intent(inout) :: E_Trial                ! Energy difference (HPS + Ele)
  type(Displacement), intent(in) :: disp(:)         ! Displaced atom data

  ! Local variables
  integer :: i, nDisp, iType, iMol, iPair, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_new_sq, r_old_sq, r
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS_total, E_Ele_total
  real(dp) :: x1_New, y1_New, z1_New
  real(dp) :: x1_Old, y1_Old, z1_Old
  real(dp) :: x2_New, y2_New, z2_New
  real(dp) :: x2_Old, y2_Old, z2_Old
  logical :: q_Nzero
  logical :: changed(1:maxAtoms)
  integer :: dispIndx(1:maxAtoms)

  ! Initialize arrays and energy terms
  nDisp = size(disp)
  iType = disp(1)%molType
  iMol = disp(1)%molIndx
  changed = .false.
  dispIndx = 0
  E_HPS_total = 0E0_dp
  E_Ele_total = 0E0_dp

  ! Mark displaced atoms
  do i = 1, nDisp
    changed(disp(i)%atmIndx) = disp(i)%Displaced
    dispIndx(disp(i)%atmIndx) = disp(i)%atmIndx
  enddo

  ! Calculate energy changes for non-bonded pairs with displaced atoms
  do iPair = 1, nIntraNonBond(iType)
    iAtom = nonBondArray(iType, iPair)%nonMembr(1)
    jAtom = nonBondArray(iType, iPair)%nonMembr(2)
    if (changed(iAtom) .or. changed(jAtom)) then
      atmType1 = atomArray(iType, iAtom)
      atmType2 = atomArray(iType, jAtom)

      ! Retrieve force field parameters
      eps = eps_tab(atmType1, atmType2)
      sig_sq = sigsq_tab(atmType1, atmType2)
      rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
      q_Nzero = q_Nonzero(atmType1, atmType2)
      lambda = lambda_tab(atmType1, atmType2)
      q = q_tab(atmType1, atmType2)

      ! Get coordinates for new position
      x1_New = disp(dispIndx(iAtom))%x_new
      y1_New = disp(dispIndx(iAtom))%y_new
      z1_New = disp(dispIndx(iAtom))%z_new
      x2_New = disp(dispIndx(jAtom))%x_new
      y2_New = disp(dispIndx(jAtom))%y_new
      z2_New = disp(dispIndx(jAtom))%z_new

      ! Calculate distance for new position
      rx = x2_New - x1_New
      ry = y2_New - y1_New
      rz = z2_New - z1_New
      r_new_sq = rx * rx + ry * ry + rz * rz

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

      ! Calculate electrostatic energy for new position
      Ele = 0E0_dp
      if (q_Nzero .and. r_new_sq < rcutElecsq) then
        r = sqrt(r_new_sq)
        if (kapa * r < 50.0_dp) then  ! Prevent underflow
          Ele = q * exp(-kapa * r) / r
        endif
      endif
      E_Ele_total = E_Ele_total + Ele

      ! Get coordinates for old position
      x1_Old = disp(dispIndx(iAtom))%x_old
      y1_Old = disp(dispIndx(iAtom))%y_old
      z1_Old = disp(dispIndx(iAtom))%z_old
      x2_Old = disp(dispIndx(jAtom))%x_old
      y2_Old = disp(dispIndx(jAtom))%y_old
      z2_Old = disp(dispIndx(jAtom))%z_old

      ! Calculate distance for old position
      rx = x2_Old - x1_Old
      ry = y2_Old - y1_Old
      rz = z2_Old - z1_Old
      r_old_sq = rx * rx + ry * ry + rz * rz

      ! Calculate HPS energy for old position
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

      ! Calculate electrostatic energy for old position
      Ele = 0E0_dp
      if (q_Nzero .and. r_old_sq < rcutElecsq) then
        r = sqrt(r_old_sq)
        if (kapa * r < 50.0_dp) then  ! Prevent underflow
          Ele = q * exp(-kapa * r) / r
        endif
      endif
      E_Ele_total = E_Ele_total - Ele
    endif
  enddo

  ! Update total energy difference
  E_Trial = E_HPS_total + E_Ele_total

end subroutine Shift_ECalc_IntraNonBonded
!======================================================================================
! Computes the intramolecular non-bonded energy of a specified molecule in a Grand Canonical
! Monte Carlo simulation using the HPS_piecewise force field. Calculates E_Trial (piecewise
! Lennard-Jones with switching at (sigma/r)^6 = 0.5 and Debye-Hückel electrostatics) for
! non-bonded atom pairs within molecule (iType, iMol). Used in swap-out moves (e.g.,
! SwapOut_ECalc_HPS_piecewise) for flexible molecules.
pure subroutine Mol_ECalc_IntraNonBonded(iType, iMol, E_Trial)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: atomArray, nonBondArray, nIntraNonBond
  use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, &
                                        cutoffNBsq_tab, lambda_tab, rcutElec, rcutElecsq
  use Coords, only: MolArray
  use SimParameters, only: kapa
  implicit none

  ! Input/Output variables
  integer, intent(in) :: iType, iMol      ! Molecule type and index
  real(dp), intent(out) :: E_Trial        ! Intramolecular non-bonded energy

  ! Local variables
  integer :: iPair, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_sq, r, rmax
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS_total, E_Ele_total
  logical :: q_Nzero

  ! Initialize energy terms
  E_HPS_total = 0E0_dp
  E_Ele_total = 0E0_dp
  E_Trial = 0E0_dp

  ! Calculate intramolecular non-bonded energies
  do iPair = 1, nIntraNonBond(iType)
    iAtom = nonBondArray(iType, iPair)%nonMembr(1)
    jAtom = nonBondArray(iType, iPair)%nonMembr(2)
    ! Calculate distance
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
    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif

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

    ! Calculate electrostatic energy
    Ele = 0E0_dp
    if (q_Nzero) then
      q = q_tab(atmType1, atmType2)
      r = sqrt(r_sq)
      if (kapa * r < 50.0_dp) then  ! Prevent underflow
        Ele = q * exp(-kapa * r) / r
      endif
      E_Ele_total = E_Ele_total + Ele
    endif
  enddo

  ! Compute total energy
  E_Trial = E_HPS_total + E_Ele_total

end subroutine Mol_ECalc_IntraNonBonded
!======================================================================================      
! Calculates the intramolecular non-bonded energy for a newly inserted molecule in a
! Grand Canonical Monte Carlo nucleation simulation using the HPS_piecewise force field.
! Computes E_Trial (piecewise Lennard-Jones with switching at (sigma/r)^6 = 0.5 and
! Debye-Hückel electrostatics) for non-bonded atom pairs within newMol, as defined by
! nonBondArray. Applies non-bonded (rcutNBsq) and electrostatic (rcutElecsq) cutoffs.
! Used in swap-in moves (e.g., AVBMC). No overlap checks are performed.
pure subroutine NewMol_ECalc_IntraNonBonded(E_Trial)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: nIntraNonBond, nonBondArray, atomArray
  use ForceFieldPara_HPS_piecewise, only: rcutElec, rcutElecsq, eps_tab, lambda_tab, &
                                        sigsq_tab, q_Nonzero, q_tab, cutoffNBsq_tab
  use Coords, only: newMol
  use SimParameters, only: kapa
  implicit none

  ! Input/Output variables
  real(dp), intent(out) :: E_Trial        ! Intramolecular non-bonded energy

  ! Local variables
  integer :: iType, iPair, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_sq, r, rmax
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS_total, E_Ele_total
  logical :: q_Nzero

  ! Initialize energy terms
  E_HPS_total = 0E0_dp
  E_Ele_total = 0E0_dp
  E_Trial = 0E0_dp
  iType = newMol%molType

  ! Calculate intramolecular non-bonded energies
  do iPair = 1, nIntraNonBond(iType)
    iAtom = nonBondArray(iType, iPair)%nonMembr(1)
    jAtom = nonBondArray(iType, iPair)%nonMembr(2)
    ! Calculate distance
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
    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif

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

    ! Calculate electrostatic energy
    Ele = 0E0_dp
    if (q_Nzero) then
      q = q_tab(atmType1, atmType2)
      r = sqrt(r_sq)
      if (kapa * r < 50.0_dp) then  ! Prevent underflow
        Ele = q * exp(-kapa * r) / r
      endif
      E_Ele_total = E_Ele_total + Ele
    endif
  enddo

  ! Compute total energy
  E_Trial = E_HPS_total + E_Ele_total

end subroutine NewMol_ECalc_IntraNonBonded
!======================================================================================   
      end module
      
      
