!*********************************************************************************************************************
!     This file contains the energy functions that work for Wang–Frenkel w/ Columbic style forcefields
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
!   Mpipi: Wang–Frankel + Debye-Hückel Electrostatics
!      U_WF(r)           = ε * α * [ ( (σ / r)**(2μ) − 1 ) − ( (3σ / r)**(2μ) − 1 ) ]**(2ν)
!      U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!*********************************************************************************************************************
      module IntraEnergy_Mpipi
      use VarPrecision
      contains
!====================================================================================== 
pure function WF_Func(r_sq, epAlpha, sig_sq, Mu) result(WF)

  ! Computes the Wang–Frenkel potential function value given squared distance and parameters.
  ! Uses integer Mu for generalization; optimized for Mu = 2 or 3 with precomputed constants.

  implicit none
  real(dp), intent(in) :: r_sq, epAlpha, sig_sq
  integer, intent(in) :: Mu
  real(dp) :: WF, x, y, k

  x = sig_sq / r_sq
  if(Mu == 2) then
    x = x * x
    k = 81.0_dp  ! 9^2
  elseif(Mu == 3) then
    x = x * x * x
    k = 729.0_dp  ! 9^3
  else
    x = x ** Mu
    k = 9.0_dp ** Mu
  endif

  y = k * x - 1.0_dp
  WF = epAlpha * (x - 1.0_dp) * y * y

end function WF_Func
!======================================================================================      
! Calculates intramolecular non-bonded energies for the Mpipi force field in a grand
! canonical ensemble nucleation simulation. Iterates over specified non-bonded atom
! pairs within each molecule, computing Wang-Frenkel and screened electrostatic
! interactions with cutoffs. Checks for atomic overlaps. Updates total system energy
! (E_T) and non-bonded energy (E_NBond_T). Used in Detailed_ECalc to compute system
! energy. Does not use NeighborList or NeighborPairs, as it handles intramolecular
! interactions, not intermolecular clustering.
subroutine Detailed_ECalc_IntraNonBonded(E_T)

  ! Calculates intramolecular non-bonded energy using the Mpipi force field.
  ! Iterates over non-bonded atom pairs in each molecule and computes WF + electrostatics.
  ! Rejects overlaps and updates total and non-bonded energy accumulators.

  use ParallelVar, only: nout
  use ForceFieldPara_Mpipi, only: epAlpha_tab, q_Nonzero, q_tab, sigsq_tab, Mu_tab, rcutElec, rcutElecsq
  use Coords, only: MolArray
  use SimParameters, only: nMolTypes, NPART, kapa
  use ForceField, only: nIntraNonBond, r_min_tab, atomArray, nonBondArray
  use EnergyTables, only: E_NBond_T

  implicit none

  real(dp), intent(inout) :: E_T

  integer :: iType, iMol, iPair, iAtom, jAtom, Mu
  integer :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r, rmin_ij
  real(dp) :: epAlpha, sig_sq, q
  logical :: q_Nzero
  real(dp) :: WF, Ele, E_Ele, E_WF
  real(dp) :: rcutNBsq, rmax

  E_WF = 0.0_dp
  E_Ele = 0.0_dp

  do iType = 1, nMolTypes
    do iMol = 1, NPART(iType)
      do iPair = 1, nIntraNonBond(iType)

        iAtom = nonBondArray(iType, iPair)%nonMembr(1)
        jAtom = nonBondArray(iType, iPair)%nonMembr(2)
        atmType1 = atomArray(iType, iAtom)
        atmType2 = atomArray(iType, jAtom)

        rx = abs(MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(iType)%mol(iMol)%x(jAtom))
        ry = abs(MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(iType)%mol(iMol)%y(jAtom))
        rz = abs(MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(iType)%mol(iMol)%z(jAtom))
        rmax = max(rx, ry, rz)
        if (rmax > rcutElec) cycle

        epAlpha = epAlpha_tab(atmType1, atmType2)
        sig_sq  = sigsq_tab(atmType1, atmType2)
        rcutNBsq = 9.0_dp * sig_sq
        q_Nzero = q_Nonzero(atmType1, atmType2)
        Mu = Mu_tab(atmType1, atmType2)
        rmin_ij = r_min_tab(atmType1, atmType2)

        if (.not. q_Nzero) then
          rmax = rmax * rmax
          if (rmax > rcutNBsq) cycle
        endif

        r = rx**2 + ry**2 + rz**2
        if (r > rcutElecsq) cycle
        if (r < rmin_ij) then
          stop "ERROR! Overlapping atoms found in the current configuration!"
        endif

        WF = 0.0_dp
        if (r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
        E_WF = E_WF + WF

        Ele = 0.0_dp
        if (q_Nzero) then
          q = q_tab(atmType1, atmType2)
          r = sqrt(r)
          Ele = q * exp(-kapa * r) / r
        endif
        E_Ele = E_Ele + Ele

      enddo
    enddo
  enddo

  write(nout,*) "Intra Wang–Frenkel Energy:", E_WF
  write(nout,*) "Intra Electrostatic Energy:", E_Ele

  E_T = E_T + E_Ele + E_WF
  E_NBond_T = E_Ele + E_WF

end subroutine Detailed_ECalc_IntraNonBonded
!======================================================================================  
! Calculates intramolecular non-bonded energies for a single molecule of each type in
! the gas phase using the Mpipi force field in a grand canonical ensemble nucleation
! simulation. Iterates over specified non-bonded atom pairs, computing Wang-Frenkel and
! screened electrostatic interactions with cutoffs. Checks for atomic overlaps. Stores
! results in E_gasIntra(iType). Does not use NeighborList or NeighborPairs, as it handles
! intramolecular interactions for one molecule per type.
subroutine Detailed_ECalc_GasIntra

  ! Calculates intramolecular non-bonded energies for a single gas-phase molecule of each type.
  ! Uses the Mpipi force field and updates E_gasIntra(iType) for use in reference state corrections.

  use ParallelVar, only: nout
  use ForceFieldPara_Mpipi, only: epAlpha_tab, q_Nonzero, q_tab, sigsq_tab, Mu_tab, &
                                  rcutElec, rcutElecsq 
  use ForceField, only: r_min_tab, nIntraNonBond, atomArray, nonBondArray
  use Coords, only: gasConfig
  use SimParameters, only: nMolTypes, kapa
  use EnergyTables, only: E_gasIntra

  implicit none

  integer :: iType, iPair, iAtom, jAtom, Mu
  integer :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r, rmin_ij
  real(dp) :: epAlpha, sig_sq, q
  logical :: q_Nzero
  real(dp) :: WF, Ele
  real(dp) :: E_Ele, E_WF, rcutNBsq, rmax

  do iType = 1, nMolTypes
    E_WF = 0.0_dp
    E_Ele = 0.0_dp

    do iPair = 1, nIntraNonBond(iType)
      iAtom = nonBondArray(iType, iPair)%nonMembr(1)
      jAtom = nonBondArray(iType, iPair)%nonMembr(2)
      atmType1 = atomArray(iType, iAtom)
      atmType2 = atomArray(iType, jAtom)

      rx = abs(gasConfig(iType)%x(iAtom) - gasConfig(iType)%x(jAtom))
      ry = abs(gasConfig(iType)%y(iAtom) - gasConfig(iType)%y(jAtom))
      rz = abs(gasConfig(iType)%z(iAtom) - gasConfig(iType)%z(jAtom))
      rmax = max(rx, ry, rz)
      if (rmax > rcutElec) cycle

      epAlpha = epAlpha_tab(atmType1, atmType2)
      sig_sq  = sigsq_tab(atmType1, atmType2)
      rcutNBsq = 9.0_dp * sig_sq
      q_Nzero = q_Nonzero(atmType1, atmType2)
      Mu = Mu_tab(atmType1, atmType2)
      rmin_ij = r_min_tab(atmType1, atmType2)

      if (.not. q_Nzero) then
        rmax = rmax * rmax
        if (rmax > rcutNBsq) cycle
      endif

      r = rx**2 + ry**2 + rz**2
      if (r > rcutElecsq) cycle
      if (r < rmin_ij) then
        stop "ERROR! Overlapping atoms found in the current configuration!"
      endif

      WF = 0.0_dp
      if (r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
      E_WF = E_WF + WF

      Ele = 0.0_dp
      if (q_Nzero) then
        q = q_tab(atmType1, atmType2)
        r = sqrt(r)
        Ele = q * exp(-kapa * r) / r
      endif
      E_Ele = E_Ele + Ele

    enddo

    E_gasIntra(iType) = E_Ele + E_WF
    write(nout,*) "Type:", iType, "Gas Intra Wang–Frenkel Energy:", E_WF
    write(nout,*) "Type:", iType, "Gas Intra Electrostatic Energy:", E_Ele

  enddo

end subroutine Detailed_ECalc_GasIntra
!======================================================================================   
! Calculates the intramolecular non-bonded energy change for a molecular displacement in a Grand Canonical
! Monte Carlo nucleation simulation using the Mpipi force field. Computes the energy difference
! (E_Trial = E_WF + E_Ele) for Wang-Frenkel and electrostatic interactions between atom pairs listed in
! nonBondArray when at least one atom is displaced.
pure subroutine Shift_ECalc_IntraNonBonded(E_Trial, disp)

  ! Calculates the intramolecular non-bonded energy change for a molecular displacement
  ! in a Grand Canonical Monte Carlo simulation using the Mpipi force field.

  use VarPrecision, only: dp, atomIntType
  use SimParameters, only: maxAtoms, kapa
  use CoordinateTypes, only: Displacement
  use ForceField, only: atomArray, nonBondArray, nIntraNonBond
  use ForceFieldPara_Mpipi, only: epAlpha_tab, sigsq_tab, q_tab, q_Nonzero, Mu_tab, rcutElecsq

  implicit none

  real(dp), intent(inout) :: E_Trial
  type(Displacement), intent(in) :: disp(:)

  logical :: changed(1:maxAtoms)
  integer :: dispIndx(1:maxAtoms)
  integer :: i, nDisp
  integer :: iType, iMol, iPair, iAtom, jAtom, Mu
  integer :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: epAlpha, sig_sq, q
  logical :: q_Nzero
  real(dp) :: WF, Ele
  real(dp) :: E_Ele, E_WF, rcutNBsq
  real(dp) :: x1_New, y1_New, z1_New
  real(dp) :: x1_Old, y1_Old, z1_Old
  real(dp) :: x2_New, y2_New, z2_New
  real(dp) :: x2_Old, y2_Old, z2_Old

  nDisp = size(disp)
  iType = disp(1)%molType
  iMol  = disp(1)%molIndx
  changed = .false.
  E_WF = 0.0_dp
  E_Ele = 0.0_dp

  do i = 1, nDisp
    changed(disp(i)%atmIndx) = disp(i)%Displaced
    dispIndx(disp(i)%atmIndx) = disp(i)%atmIndx
  enddo

  do iPair = 1, nIntraNonBond(iType)
    iAtom = nonBondArray(iType, iPair)%nonMembr(1)
    jAtom = nonBondArray(iType, iPair)%nonMembr(2)

    if (changed(iAtom) .or. changed(jAtom)) then
      atmType1 = atomArray(iType, iAtom)
      atmType2 = atomArray(iType, jAtom)

      epAlpha = epAlpha_tab(atmType1, atmType2)
      sig_sq  = sigsq_tab(atmType1, atmType2)
      rcutNBsq = 9.0_dp * sig_sq
      q_Nzero = q_Nonzero(atmType1, atmType2)
      Mu = Mu_tab(atmType1, atmType2)

      x1_Old = disp(dispIndx(iAtom))%x_old
      y1_Old = disp(dispIndx(iAtom))%y_old
      z1_Old = disp(dispIndx(iAtom))%z_old
      x1_New = disp(dispIndx(iAtom))%x_new
      y1_New = disp(dispIndx(iAtom))%y_new
      z1_New = disp(dispIndx(iAtom))%z_new

      x2_Old = disp(dispIndx(jAtom))%x_old
      y2_Old = disp(dispIndx(jAtom))%y_old
      z2_Old = disp(dispIndx(jAtom))%z_old
      x2_New = disp(dispIndx(jAtom))%x_new
      y2_New = disp(dispIndx(jAtom))%y_new
      z2_New = disp(dispIndx(jAtom))%z_new

      ! New configuration
      rx = x2_New - x1_New
      ry = y2_New - y1_New
      rz = z2_New - z1_New
      r = rx*rx + ry*ry + rz*rz
      WF = 0.0_dp
      if (r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
      E_WF = E_WF + WF

      Ele = 0.0_dp
      if (q_Nzero) then
        if (r < rcutElecsq) then
          q = q_tab(atmType1, atmType2)
          r = sqrt(r)
          Ele = q * exp(-kapa * r) / r
          E_Ele = E_Ele + Ele
        endif
      endif

      ! Old configuration
      rx = x2_Old - x1_Old
      ry = y2_Old - y1_Old
      rz = z2_Old - z1_Old
      r = rx*rx + ry*ry + rz*rz
      WF = 0.0_dp
      if (r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
      E_WF = E_WF - WF

      if (q_Nzero) then
        if (r < rcutElecsq) then
          q = q_tab(atmType1, atmType2)
          r = sqrt(r)
          Ele = q * exp(-kapa * r) / r
          E_Ele = E_Ele - Ele
        endif
      endif
    endif
  enddo

  E_Trial = E_WF + E_Ele

end subroutine Shift_ECalc_IntraNonBonded
!======================================================================================
! Computes the intramolecular non-bonded energy of a specified molecule in a Grand Canonical Monte Carlo
! simulation using the Mpipi force field. Calculates E_Trial (Wang-Frenkel potential and Debye-Hückel
! electrostatics) for non-bonded atom pairs within molecule (iType, iMol). Used in swap-out moves
! (e.g., SwapOut_ECalc_Mpipi) and other energy calculations for flexible molecules.
pure subroutine Mol_ECalc_IntraNonBonded(iType, iMol, E_Trial)

  ! Computes the intramolecular non-bonded energy of a specified molecule in a Grand Canonical
  ! Monte Carlo simulation using the Mpipi force field. Calculates E_Trial (Wang-Frenkel potential
  ! and Debye-Hückel electrostatics) for non-bonded atom pairs within molecule (iType, iMol).
  ! Used in swap-out moves (e.g., SwapOut_ECalc_Mpipi) and other energy calculations for flexible molecules.

  use VarPrecision
  use ForceField, only: atomArray, nonBondArray, nIntraNonBond
  use ForceFieldPara_Mpipi, only: epAlpha_tab, q_Nonzero, q_tab, sigsq_tab, Mu_tab, rcutElec, rcutElecsq
  use Coords, only: MolArray
  use SimParameters, only: kapa

  implicit none

  integer, intent(in) :: iType, iMol
  real(dp), intent(out) :: E_Trial

  integer :: iPair, iAtom, jAtom, Mu
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: epAlpha, sig_sq, q
  logical :: q_Nzero
  real(dp) :: WF, Ele
  real(dp) :: E_Ele, E_WF, rcutNBsq, rmax

  E_WF = 0.0_dp
  E_Ele = 0.0_dp
  E_Trial = 0.0_dp

  do iPair = 1, nIntraNonBond(iType)
    iAtom = nonBondArray(iType, iPair)%nonMembr(1)
    jAtom = nonBondArray(iType, iPair)%nonMembr(2)

    rx = abs(MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(iType)%mol(iMol)%x(jAtom))
    ry = abs(MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(iType)%mol(iMol)%y(jAtom))
    rz = abs(MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(iType)%mol(iMol)%z(jAtom))
    rmax = max(rx, ry, rz)
    if (rmax > rcutElec) cycle

    atmType1 = atomArray(iType, iAtom)
    atmType2 = atomArray(iType, jAtom)
    epAlpha = epAlpha_tab(atmType1, atmType2)
    sig_sq = sigsq_tab(atmType1, atmType2)
    rcutNBsq = 9.0_dp * sig_sq
    q_Nzero = q_Nonzero(atmType1, atmType2)
    Mu = Mu_tab(atmType1, atmType2)

    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif

    r = rx*rx + ry*ry + rz*rz
    if (r > rcutElecsq) cycle

    WF = 0.0_dp
    if (r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
    E_WF = E_WF + WF

    if (q_Nzero) then
      q = q_tab(atmType1, atmType2)
      r = sqrt(r)
      Ele = q * exp(-kapa * r) / r
      E_Ele = E_Ele + Ele
    endif

  enddo

  E_Trial = E_WF + E_Ele

end subroutine Mol_ECalc_IntraNonBonded
!======================================================================================      
! Calculates the intramolecular non-bonded energy for a newly inserted molecule in a Grand Canonical
! Monte Carlo nucleation simulation using the Mpipi force field. Computes E_Trial (E_WF + E_Ele) for
! non-bonded atom pairs within the molecule, as defined by nonBondArray, using the Wang-Frenkel
! potential and Debye-Hückel electrostatics. Applies non-bonded cutoff (rcutNBsq = 9 * sigma^2) and
! electrostatic cutoff (rcutElecsq). No overlap checks are performed. Used in swap-in moves (e.g., AVBMC).
subroutine NewMol_ECalc_IntraNonBonded(E_Trial)

  ! Calculates the intramolecular non-bonded energy for a newly inserted molecule in a Grand Canonical
  ! Monte Carlo nucleation simulation using the Mpipi force field. Computes E_Trial (E_WF + E_Ele) for
  ! non-bonded atom pairs within the molecule, as defined by nonBondArray, using the Wang-Frenkel
  ! potential and Debye-Hückel electrostatics. Applies non-bonded cutoff (rcutNBsq = 9 * sigma^2) and
  ! electrostatic cutoff (rcutElecsq). No overlap checks are performed. Used in swap-in moves (e.g., AVBMC).

  use ForceFieldPara_Mpipi, only: epAlpha_tab, q_Nonzero, sigsq_tab, Mu_tab, q_tab, rcutElec, rcutElecsq
  use Coords, only: newMol
  use SimParameters, only: kapa
  use ForceField, only: nIntraNonBond, nonBondArray, atomArray
  use VarPrecision

  implicit none

  real(dp), intent(out) :: E_Trial

  integer :: iType, iPair, iAtom, jAtom, Mu
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: epAlpha, sig_sq, q
  logical :: q_Nzero
  real(dp) :: WF, Ele
  real(dp) :: E_Ele, E_WF, rcutNBsq, rmax

  E_WF = 0.0_dp
  E_Ele = 0.0_dp
  E_Trial = 0.0_dp

  iType = newMol%molType

  do iPair = 1, nIntraNonBond(iType)
    iAtom = nonBondArray(iType, iPair)%nonMembr(1)
    jAtom = nonBondArray(iType, iPair)%nonMembr(2)

    rx = abs(newMol%x(iAtom) - newMol%x(jAtom))
    ry = abs(newMol%y(iAtom) - newMol%y(jAtom))
    rz = abs(newMol%z(iAtom) - newMol%z(jAtom))
    rmax = max(rx, ry, rz)
    if (rmax > rcutElec) cycle

    atmType1 = atomArray(iType, iAtom)
    atmType2 = atomArray(iType, jAtom)

    epAlpha   = epAlpha_tab(atmType1, atmType2)
    sig_sq    = sigsq_tab(atmType1, atmType2)
    rcutNBsq  = 9.0_dp * sig_sq
    q_Nzero   = q_Nonzero(atmType1, atmType2)
    Mu        = Mu_tab(atmType1, atmType2)

    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif

    r = rx*rx + ry*ry + rz*rz
    if (r > rcutElecsq) cycle

    WF = 0.0_dp
    if (r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
    E_WF = E_WF + WF

    if (q_Nzero) then
      q = q_tab(atmType1, atmType2)
      r = sqrt(r)
      Ele = q * exp(-kapa * r) / r
      E_Ele = E_Ele + Ele
    endif

  enddo

  E_Trial = E_WF + E_Ele

end subroutine NewMol_ECalc_IntraNonBonded
!======================================================================================   
      end module
      
      
