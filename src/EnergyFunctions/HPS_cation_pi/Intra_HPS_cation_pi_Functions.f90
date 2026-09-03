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
!
!   HPS_cation_pi (HPS-cation-pi-i and -ii): HPS_piecewise + LJ for specific π interactions
!      - U_total(r)        = U_HPS_piecewise(r) + U_extra_LJ(r) for cation-π pairs
!      - U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!*********************************************************************************************************************
      module IntraEnergy_HPS_cation_pi
      use VarPrecision
      contains
!======================================================================================      
subroutine Detailed_ECalc_IntraNonBonded(E_T)
  ! Calculates intramolecular non-bonded energies for the HPS_cation_pi force field in a
  ! grand canonical ensemble nucleation simulation. Computes piecewise hydrophobic-polar
  ! (HPS), cation-pi, and screened electrostatic energies for specified non-bonded atom pairs
  ! within each molecule. Checks for atomic overlaps and stops if found. Updates total system
  ! energy (E_T) and non-bonded energy (E_NBond_T). Used in Detailed_ECalc for system energy
  ! computation. Does not use NeighborList or NeighborPairs, as it handles intramolecular
  ! interactions only.
  ! HPS_cation_pi potential:
  !   u_LJ(r) = 4 * eps * [ (sigma/r)^12 - (sigma/r)^6 ]
  !   U_HPS(r) = u_LJ(r) + eps * (1 - lambda), if r <= 2^(1/6) * sigma
  !              lambda * u_LJ(r), if r > 2^(1/6) * sigma
  !   U_cation_pi(r) = 4 * epsLJ * [ (sigma/r)^12 - (sigma/r)^6 ], for specific pairs
  !   U_Electrostatic(r) = (q1 * q2 / r) * exp(-kappa * r)

  use ParallelVar, only: nout
  use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, &
                                          lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
  use ForceField, only: r_min_tab, nIntraNonBond, atomArray, nonBondArray
  use Coords, only: MolArray
  use SimParameters, only: nMolTypes, NPART, kapa
  use EnergyTables, only: E_NBond_T
  implicit none

  ! Input/Output variables
  real(dp), intent(inout) :: E_T                   ! Total system energy

  ! Local variables
  integer :: iType, iMol, iPair, iAtom, jAtom
  integer :: atmType1, atmType2      
  real(dp) :: rx, ry, rz, r, rmin_ij
  real(dp) :: eps, sig_sq, q, lambda, epsLJ
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x, E_Ele, E_HPS
  real(dp) :: rcutNBsq, rmax

  ! Initialize energy variables
  E_HPS = 0.0_dp
  E_Ele = 0.0_dp      

  ! Loop over molecule types
  do iType = 1, nMolTypes
    ! Loop over molecules of current type
    do iMol = 1, NPART(iType)
      ! Loop over non-bonded atom pairs
      do iPair = 1, nIntraNonBond(iType)
        iAtom = nonBondArray(iType, iPair)%nonMembr(1)
        jAtom = nonBondArray(iType, iPair)%nonMembr(2)
        atmType1 = atomArray(iType, iAtom)
        atmType2 = atomArray(iType, jAtom)
        ! Calculate distance components
        rx = abs(MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(iType)%mol(iMol)%x(jAtom))
        ry = abs(MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(iType)%mol(iMol)%y(jAtom))
        rz = abs(MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(iType)%mol(iMol)%z(jAtom))
        rmax = max(rx, ry, rz)
        ! Check electrostatic cutoff
        if (rmax > rcutElec) cycle
        ! Retrieve force field parameters
        eps = eps_tab(atmType1, atmType2)
        epsLJ = epsLJ_tab(atmType1, atmType2)
        sig_sq = sigsq_tab(atmType1, atmType2)
        rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
        q_Nzero = q_Nonzero(atmType1, atmType2)
        lambda = lambda_tab(atmType1, atmType2)
        rmin_ij = r_min_tab(atmType1, atmType2) 
        ! Check non-bonded cutoff for non-charged pairs
        if (.not. q_Nzero) then
          rmax = rmax * rmax
          if (rmax > rcutNBsq) cycle
        endif
        ! Compute squared distance
        r = rx**2 + ry**2 + rz**2
        if (r > rcutElecsq) cycle
        ! Check for atomic overlaps
        if (r < rmin_ij) then
          stop "ERROR! Overlapping atoms found in the current configuration!"
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
        ! Calculate electrostatic energy
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
  ! Output energy contributions
  write(nout, *) "Intra Hydrophobicity scale Energy:", E_HPS
  write(nout, *) "Intra Electrostatic Energy:", E_Ele
  ! Update total system and non-bonded energies
  E_T = E_T + E_Ele + E_HPS    
  E_NBond_T = E_Ele + E_HPS  
      
end subroutine Detailed_ECalc_IntraNonBonded
!======================================================================================  
subroutine Detailed_ECalc_GasIntra
  ! Calculates intramolecular non-bonded energies for a single molecule of each type in
  ! the gas phase using the HPS_cation_pi force field in a grand canonical ensemble
  ! nucleation simulation. Computes piecewise hydrophobic-polar (HPS), cation-pi, and
  ! screened electrostatic energies for specified non-bonded atom pairs. Checks for
  ! atomic overlaps and stops if found. Stores results in E_gasIntra(iType). Does not
  ! use NeighborList or NeighborPairs, as it handles intramolecular interactions only.
  ! HPS_cation_pi potential:
  !   u_LJ(r) = 4 * eps * [ (sigma/r)^12 - (sigma/r)^6 ]
  !   U_HPS(r) = u_LJ(r) + eps * (1 - lambda), if r <= 2^(1/6) * sigma
  !              lambda * u_LJ(r), if r > 2^(1/6) * sigma
  !   U_cation_pi(r) = 4 * epsLJ * [ (sigma/r)^12 - (sigma/r)^6 ], for specific pairs
  !   U_Electrostatic(r) = (q1 * q2 / r) * exp(-kappa * r)

  use ParallelVar, only: nout
  use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, &
                                          lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
  use ForceField, only: r_min_tab, atomArray, nonBondArray, nIntraNonBond
  use Coords, only: gasConfig
  use SimParameters, only: nMolTypes, kapa
  use EnergyTables, only: E_gasIntra
  implicit none

  ! Local variables
  integer :: iType, iPair, iAtom, jAtom
  integer :: atmType1, atmType2      
  real(dp) :: rx, ry, rz, r, rmin_ij
  real(dp) :: eps, sig_sq, q, lambda, epsLJ
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS, rcutNBsq, rmax

  ! Loop over molecule types
  do iType = 1, nMolTypes
    ! Initialize energy variables
    E_HPS = 0.0_dp
    E_Ele = 0.0_dp  
    ! Loop over non-bonded atom pairs
    do iPair = 1, nIntraNonBond(iType)
      iAtom = nonBondArray(iType, iPair)%nonMembr(1)
      jAtom = nonBondArray(iType, iPair)%nonMembr(2)
      atmType1 = atomArray(iType, iAtom)
      atmType2 = atomArray(iType, jAtom)
      ! Calculate distance components
      rx = abs(gasConfig(iType)%x(iAtom) - gasConfig(iType)%x(jAtom))
      ry = abs(gasConfig(iType)%y(iAtom) - gasConfig(iType)%y(jAtom))
      rz = abs(gasConfig(iType)%z(iAtom) - gasConfig(iType)%z(jAtom))
      rmax = max(rx, ry, rz) 
      ! Check electrostatic cutoff
      if (rmax > rcutElec) cycle
      ! Retrieve force field parameters
      eps = eps_tab(atmType1, atmType2)
      epsLJ = epsLJ_tab(atmType1, atmType2)
      sig_sq = sigsq_tab(atmType1, atmType2)
      rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
      q_Nzero = q_Nonzero(atmType1, atmType2)
      lambda = lambda_tab(atmType1, atmType2)
      rmin_ij = r_min_tab(atmType1, atmType2)
      ! Check non-bonded cutoff for non-charged pairs
      if (.not. q_Nzero) then
        rmax = rmax * rmax
        if (rmax > rcutNBsq) cycle
      endif
      ! Compute squared distance
      r = rx**2 + ry**2 + rz**2
      if (r > rcutElecsq) cycle
      ! Check for atomic overlaps
      if (r < rmin_ij) then
        stop "ERROR! Overlapping atoms found in the current configuration!"
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
      ! Calculate electrostatic energy
      Ele = 0.0_dp
      if (q_Nzero) then
        q = q_tab(atmType1, atmType2)
        r = sqrt(r)
        Ele = q * exp(-kapa * r) / r
      endif
      E_Ele = E_Ele + Ele
    enddo
    ! Store total intramolecular energy for the molecule type
    E_gasIntra(iType) = E_Ele + E_HPS
    ! Output energy contributions
    write(nout, *) "Type: ", iType, " Gas Intra Hydrophobicity scale Energy:", E_HPS
    write(nout, *) "Type: ", iType, " Gas Intra Electrostatic Energy:", E_Ele
  enddo
      
end subroutine Detailed_ECalc_GasIntra
!======================================================================================   
pure subroutine Shift_ECalc_IntraNonBonded(E_Trial, disp)
  ! Calculates the intramolecular non-bonded energy change for a molecular displacement in a
  ! Grand Canonical Monte Carlo nucleation simulation using the HPS_cation_pi force field.
  ! Computes the energy difference (E_Trial = E_HPS + E_Ele) for HPS (piecewise Lennard-Jones
  ! with cation-pi term) and electrostatic interactions between atom pairs listed in
  ! nonBondArray when at least one atom is displaced.
  ! HPS_cation_pi potential:
  !   u_LJ(r) = 4 * eps * [ (sigma/r)^12 - (sigma/r)^6 ]
  !   U_HPS(r) = u_LJ(r) + eps * (1 - lambda), if r <= 2^(1/6) * sigma
  !              lambda * u_LJ(r), if r > 2^(1/6) * sigma
  !   U_cation_pi(r) = 4 * epsLJ * [ (sigma/r)^12 - (sigma/r)^6 ], for specific pairs
  !   U_Electrostatic(r) = (q1 * q2 / r) * exp(-kappa * r)

  use VarPrecision, only: dp, atomIntType
  use SimParameters, only: kapa, maxAtoms
  use CoordinateTypes, only: Displacement
  use ForceField, only: atomArray, nonBondArray, nIntraNonBond
  use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, sigsq_tab, cutoffNBsq_tab, &
                                          q_tab, lambda_tab, q_Nonzero, rcutElecsq
  implicit none

  ! Input/Output variables
  real(dp), intent(inout) :: E_Trial
  type(Displacement), intent(in) :: disp(:)

  ! Local variables
  logical :: changed(1:maxAtoms)
  integer :: dispIndx(1:maxAtoms)      
  integer :: i, nDisp
  integer :: iType, iMol, iPair, iAtom, jAtom
  integer :: atmType1, atmType2      
  real(dp) :: rx, ry, rz, r
  real(dp) :: eps, sig_sq, q, lambda, epsLJ
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS, rcutNBsq
  real(dp) :: x1_New, y1_New, z1_New
  real(dp) :: x1_Old, y1_Old, z1_Old
  real(dp) :: x2_New, y2_New, z2_New
  real(dp) :: x2_Old, y2_Old, z2_Old

  ! Initialize variables
  nDisp = size(disp)
  iType = disp(1)%molType
  iMol = disp(1)%molIndx
  changed = .false.
  E_HPS = 0.0_dp
  E_Ele = 0.0_dp  

  ! Mark displaced atoms
  do i = 1, nDisp
    changed(disp(i)%atmIndx) = disp(i)%Displaced
    dispIndx(disp(i)%atmIndx) = disp(i)%atmIndx
  enddo

  ! Loop over non-bonded atom pairs
  do iPair = 1, nIntraNonBond(iType)
    iAtom = nonBondArray(iType, iPair)%nonMembr(1)
    jAtom = nonBondArray(iType, iPair)%nonMembr(2)
    ! Process pairs where at least one atom is displaced
    if (changed(iAtom) .or. changed(jAtom)) then
      atmType1 = atomArray(iType, iAtom)
      atmType2 = atomArray(iType, jAtom)
      ! Retrieve force field parameters
      eps = eps_tab(atmType1, atmType2)
      epsLJ = epsLJ_tab(atmType1, atmType2)
      sig_sq = sigsq_tab(atmType1, atmType2)
      rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
      q_Nzero = q_Nonzero(atmType1, atmType2)
      lambda = lambda_tab(atmType1, atmType2)
      ! Get old and new coordinates
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
      ! Calculate new configuration energy
      rx = x2_New - x1_New
      ry = y2_New - y1_New
      rz = z2_New - z1_New
      r = rx*rx + ry*ry + rz*rz
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
      Ele = 0.0_dp
      if (q_Nzero) then
        if (r < rcutElecsq) then
          q = q_tab(atmType1, atmType2)
          r = sqrt(r)
          Ele = q * exp(-kapa * r) / r
          E_Ele = E_Ele + Ele
        endif
      endif
      ! Calculate old configuration energy and subtract
      rx = x2_Old - x1_Old
      ry = y2_Old - y1_Old
      rz = z2_Old - z1_Old
      r = rx*rx + ry*ry + rz*rz
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
      E_HPS = E_HPS - HPS
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

  ! Compute total energy difference
  E_Trial = E_HPS + E_Ele
      
end subroutine Shift_ECalc_IntraNonBonded
!======================================================================================
pure subroutine Mol_ECalc_IntraNonBonded(iType, iMol, E_Trial)
  ! Computes the intramolecular non-bonded energy of a specified molecule in a Grand Canonical
  ! Monte Carlo simulation using the HPS_cation_pi force field. Calculates E_Trial (piecewise
  ! Hydropathy Scale Lennard-Jones with switching at (sigma/r)^6 = 0.5, cation-pi Lennard-Jones,
  ! and Debye-Hückel electrostatics) for non-bonded atom pairs within molecule (iType, iMol).
  ! Used in swap-out moves (e.g., SwapOut_ECalc_HPS_cation_pi) and other energy calculations
  ! for flexible molecules.
  ! HPS_cation_pi potential:
  !   u_LJ(r) = 4 * eps * [ (sigma/r)^12 - (sigma/r)^6 ]
  !   U_HPS(r) = u_LJ(r) + eps * (1 - lambda), if r <= 2^(1/6) * sigma
  !              lambda * u_LJ(r), if r > 2^(1/6) * sigma
  !   U_cation_pi(r) = 4 * epsLJ * [ (sigma/r)^12 - (sigma/r)^6 ], for specific pairs
  !   U_Electrostatic(r) = (q1 * q2 / r) * exp(-kappa * r)

  use VarPrecision
  use ForceField, only: atomArray, nonBondArray, nIntraNonBond
  use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, &
                                          cutoffNBsq_tab, lambda_tab, rcutElec, rcutElecsq
  use Coords, only: MolArray
  use SimParameters, only: kapa
  implicit none

  ! Input/output variables
  integer, intent(in) :: iType, iMol          ! Molecule type and index
  real(dp), intent(out) :: E_Trial           ! Total intramolecular non-bonded energy

  ! Local variables
  integer :: iPair, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: eps, sig_sq, q, lambda, epsLJ
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS, rcutNBsq, rmax

  ! Initialize energy variables
  E_HPS = 0.0_dp
  E_Ele = 0.0_dp      
  E_Trial = 0.0_dp

  ! Loop over non-bonded atom pairs
  do iPair = 1, nIntraNonBond(iType)
    iAtom = nonBondArray(iType, iPair)%nonMembr(1)
    jAtom = nonBondArray(iType, iPair)%nonMembr(2)
    ! Calculate distance components
    rx = abs(MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(iType)%mol(iMol)%x(jAtom))
    ry = abs(MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(iType)%mol(iMol)%y(jAtom))
    rz = abs(MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(iType)%mol(iMol)%z(jAtom))
    rmax = max(rx, ry, rz) 
    ! Check electrostatic cutoff
    if (rmax > rcutElec) cycle
    ! Retrieve force field parameters
    atmType1 = atomArray(iType, iAtom)
    atmType2 = atomArray(iType, jAtom)
    eps = eps_tab(atmType1, atmType2)
    epsLJ = epsLJ_tab(atmType1, atmType2)
    sig_sq = sigsq_tab(atmType1, atmType2)
    rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
    q_Nzero = q_Nonzero(atmType1, atmType2)
    lambda = lambda_tab(atmType1, atmType2)
    ! Check non-bonded cutoff for non-charged pairs
    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif
    ! Compute squared distance
    r = rx*rx + ry*ry + rz*rz
    if (r > rcutElecsq) cycle
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
pure subroutine NewMol_ECalc_IntraNonBonded(E_Trial)
  ! Calculates the intramolecular non-bonded energy for a newly inserted molecule in a Grand
  ! Canonical Monte Carlo nucleation simulation using the HPS_cation_pi force field. Computes
  ! E_Trial (E_HPS + E_Ele) for non-bonded atom pairs within the molecule, as defined by
  ! nonBondArray, using piecewise HPS Lennard-Jones, cation-pi Lennard-Jones, and Debye-Hückel
  ! electrostatics. Applies non-bonded cutoff (rcutNBsq) and electrostatic cutoff (rcutElecsq).
  ! No overlap checks are performed. Used in swap-in moves (e.g., AVBMC).
  ! HPS_cation_pi potential:
  !   u_LJ(r) = 4 * eps * [ (sigma/r)^12 - (sigma/r)^6 ]
  !   U_HPS(r) = u_LJ(r) + eps * (1 - lambda), if r <= 2^(1/6) * sigma
  !              lambda * u_LJ(r), if r > 2^(1/6) * sigma
  !   U_cation_pi(r) = 4 * epsLJ * [ (sigma/r)^12 - (sigma/r)^6 ], for specific pairs
  !   U_Electrostatic(r) = (q1 * q2 / r) * exp(-kappa * r)

  use ForceFieldPara_HPS_cation_pi, only: rcutElec, rcutElecsq, eps_tab, lambda_tab, sigsq_tab, &
                                          q_Nonzero, q_tab, cutoffNBsq_tab, epsLJ_tab
  use Coords, only: newMol
  use ForceField, only: atomArray, nIntraNonBond, nonBondArray 
  use SimParameters, only: kapa
  use VarPrecision
  implicit none

  ! Input/Output variables
  real(dp), intent(out) :: E_Trial              ! Total intramolecular energy (E_HPS + E_Ele)

  ! Local variables
  integer :: iType, iPair, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: eps, sig_sq, q, lambda, epsLJ
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS, rcutNBsq, rmax

  ! Initialize energy variables
  E_HPS = 0.0_dp
  E_Ele = 0.0_dp      
  E_Trial = 0.0_dp
  iType = newMol%molType

  ! Loop over non-bonded atom pairs
  do iPair = 1, nIntraNonBond(iType)
    iAtom = nonBondArray(iType, iPair)%nonMembr(1)
    jAtom = nonBondArray(iType, iPair)%nonMembr(2)
    ! Calculate distance components
    rx = abs(newMol%x(iAtom) - newMol%x(jAtom))
    ry = abs(newMol%y(iAtom) - newMol%y(jAtom))
    rz = abs(newMol%z(iAtom) - newMol%z(jAtom))
    rmax = max(rx, ry, rz) 
    ! Check electrostatic cutoff
    if (rmax > rcutElec) cycle
    ! Retrieve force field parameters
    atmType1 = atomArray(iType, iAtom)
    atmType2 = atomArray(iType, jAtom)
    eps = eps_tab(atmType1, atmType2)
    epsLJ = epsLJ_tab(atmType1, atmType2)
    sig_sq = sigsq_tab(atmType1, atmType2)
    rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
    q_Nzero = q_Nonzero(atmType1, atmType2)
    lambda = lambda_tab(atmType1, atmType2)
    ! Check non-bonded cutoff for non-charged pairs
    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif
    ! Compute squared distance
    r = rx*rx + ry*ry + rz*rz
    if (r > rcutElecsq) cycle
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
      
      
