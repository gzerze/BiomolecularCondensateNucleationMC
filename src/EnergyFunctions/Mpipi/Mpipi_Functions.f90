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
!           Exchange - Combines the Mol and New Mol routines for moves that simultaniously add and remove a particle at the same time.
!   Mpipi: Wang–Frankel + Debye-Hückel Electrostatics
!      U_WF(r)           = ε * α * [ ( (σ / r)**(2μ) − 1 ) − ( (3σ / r)**(2μ) − 1 ) ]**(2ν)
!      U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!*********************************************************************************************************************
      module InterEnergy_Mpipi
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
! Calculates intermolecular energies (Wang-Frenkel and electrostatic) for the Mpipi
! force field in a grand canonical ensemble nucleation simulation. Updates PairList
! with distances (distCriteria: first atom pair; minDistCriteria: minimum atom pair)
! or interaction energies (energy criterion). Populates NeighborList and NeighborPairs
! for minDistCriteria, storing atom pair indices for neighbor molecule pairs.
subroutine Detailed_ECalc_Inter(E_T, pairList)

  ! Calculates intermolecular energy for the Mpipi force field.
  ! Populates PairList and NeighborPairs for cluster rejection criteria.

  use ParallelVar, only: nout
  use ForceField, only: r_min_tab, atomArray, nAtoms
  use ForceFieldPara_Mpipi, only: epAlpha_tab, q_Nonzero, sigsq_tab, Mu_tab, q_tab, rcutElecsq
  use Coords, only: MolArray, NeighborPairs
  use SimParameters, only: nMolTypes, NPART, distCriteria, minDistCriteria, Dist_Critr_sq, kapa
  use EnergyTables, only: ETable, E_Inter_T
  use VarPrecision, only: dp, atomIntType

  implicit none

  real(dp), intent(inout) :: E_T
  real(dp), intent(inout) :: PairList(:,:)

  integer :: iType, jType, iMol, jMol, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  integer :: Mu, iIndx, jIndx, globIndx1, globIndx2, jMolMin, iPair
  real(dp) :: rx, ry, rz, r
  real(dp) :: epAlpha, sig_sq, q
  logical :: q_Nzero
  real(dp) :: WF, Ele
  real(dp) :: E_Ele, E_WF, rcutNBsq
  real(dp) :: rmin_ij

  E_WF = 0.0_dp
  E_Ele = 0.0_dp
  E_Inter_T = 0.0_dp

  PairList = huge(dp)
  if ((.not. distCriteria) .and. (.not. minDistCriteria)) PairList = 0.0_dp

  ETable = 0.0_dp

  do iType = 1, nMolTypes
    do jType = iType, nMolTypes
      do iMol = 1, NPART(iType)
        if(iType == jType) then
          jMolMin = iMol + 1
        else
          jMolMin = 1
        endif
        do jMol = jMolMin, NPART(jType)
          iIndx = MolArray(iType)%mol(iMol)%indx
          jIndx = MolArray(jType)%mol(jMol)%indx
          do iAtom = 1, nAtoms(iType)
            atmType1 = atomArray(iType, iAtom)
            globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
            do jAtom = 1, nAtoms(jType)
              atmType2 = atomArray(jType, jAtom)
              epAlpha = epAlpha_tab(atmType1, atmType2)
              q_Nzero = q_Nonzero(atmType1, atmType2)
              sig_sq = sigsq_tab(atmType1, atmType2)
              rcutNBsq = 9.0_dp * sig_sq
              Mu = Mu_tab(atmType1, atmType2)
              rmin_ij = r_min_tab(atmType1, atmType2)
              globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)

              rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r = rx**2 + ry**2 + rz**2

              ! Cluster criteria: update PairList and NeighborPairs
              if(minDistCriteria) then
                if(r < PairList(iIndx, jIndx)) PairList(iIndx, jIndx) = r
                if(r < PairList(jIndx, iIndx)) PairList(jIndx, iIndx) = r
                if(r < Dist_Critr_sq) then
                  iPair = NeighborPairs(iIndx, jIndx)%details%nPairs + 1
                  NeighborPairs(iIndx, jIndx)%details%nPairs = iPair
                  NeighborPairs(iIndx, jIndx)%details%pairIndices(2 * iPair - 1) = iAtom
                  NeighborPairs(iIndx, jIndx)%details%pairIndices(2 * iPair) = jAtom
                  iPair = NeighborPairs(jIndx, iIndx)%details%nPairs + 1
                  NeighborPairs(jIndx, iIndx)%details%nPairs = iPair
                  NeighborPairs(jIndx, iIndx)%details%pairIndices(2 * iPair - 1) = jAtom
                  NeighborPairs(jIndx, iIndx)%details%pairIndices(2 * iPair) = iAtom
                endif
              elseif(distCriteria) then
                if(iAtom == 1 .and. jAtom == 1) then
                  PairList(iIndx, jIndx) = r
                  PairList(jIndx, iIndx) = r
                endif
              endif

              ! Overlap check
              if(r < rmin_ij) then
                stop "ERROR! Overlapping atoms found in the current configuration!"
              endif

              WF = 0.0_dp
              if(r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
              E_WF = E_WF + WF

              Ele = 0.0_dp
              if(q_Nzero) then
                if(r < rcutElecsq) then
                  q = q_tab(atmType1, atmType2)
                  r = sqrt(r)
                  Ele = q * exp(-kapa * r) / r
                endif
              endif
              E_Ele = E_Ele + Ele

              if((.not. distCriteria) .and. (.not. minDistCriteria)) then
                PairList(iIndx, jIndx) = PairList(iIndx, jIndx) + Ele + WF
                PairList(jIndx, iIndx) = PairList(iIndx, jIndx)
              endif

              ETable(iIndx) = ETable(iIndx) + Ele + WF
              ETable(jIndx) = ETable(jIndx) + Ele + WF

            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  write(nout,*) "Wang–Frenkel Energy:", E_WF
  write(nout,*) "Electrostatic Energy:", E_Ele

  E_T = E_T + E_Ele + E_WF
  E_Inter_T = E_Ele + E_WF

end subroutine Detailed_ECalc_Inter
!======================================================================================  
! Calculates the intermolecular energy change for a molecular displacement in a Grand Canonical
! Monte Carlo nucleation simulation using the Mpipi force field. Computes the energy difference
! (E_Trial = E_WF + E_Ele) for displaced atoms interacting with stationary atoms, updates PairList
! (distances for distCriteria/minDistCriteria, energies otherwise), and populates NeighborDetailsNew
! with atom pair indices for minDistCriteria when provided. Supports LLTranslation for long linear
! biomolecules and Translation without NeighborDetailsNew.    
subroutine Shift_ECalc_Inter(E_Trial, disp, PairList, dETable, rejMove, NeighborDetailsNew)
  use VarPrecision, only: dp, atomIntType
  use ForceField, only: atomArray, r_min_tab, nAtoms
  use ForceFieldPara_Mpipi, only: epAlpha_tab, sigsq_tab, Mu_tab, q_tab, q_Nonzero, rcutElecsq
  use Coords, only: MolArray, Displacement, NeighborDetails
  use SimParameters, only: distCriteria, minDistCriteria, maxMol, nMolTypes, NPART, Dist_Critr_sq, kapa

  implicit none
  type(Displacement), intent(in) :: disp(:)
  real(dp), intent(out) :: E_Trial
  real(dp), intent(inout) :: PairList(:), dETable(:)
  logical, intent(out) :: rejMove
  type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)

  integer :: iDisp, jType, jMol, jAtom
  integer :: iType, iMol, iAtom, iIndx, jIndx
  integer :: sizeDisp, Mu, iPair
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_new, r_old
  real(dp) :: rmin_ij, sig_sq, epAlpha, q, WF, Ele
  real(dp) :: E_WF, E_Ele, rcutNBsq
  logical :: q_Nzero

  E_WF    = 0.0_dp
  E_Ele   = 0.0_dp
  E_Trial = 0.0_dp
  rejMove = .false.
  sizeDisp = size(disp)

  PairList = huge(dp)
  if ((.not. distCriteria) .and. (.not. minDistCriteria)) PairList = 0.0_dp
  dETable  = 0.0_dp

  iType = disp(1)%molType
  iMol  = disp(1)%molIndx
  iIndx = MolArray(iType)%mol(iMol)%indx

  do iDisp = 1, sizeDisp
    iAtom    = disp(iDisp)%atmIndx
    atmType1 = atomArray(iType, iAtom)

    do jType = 1, nMolTypes
      do jMol = 1, NPART(jType)
        if (jType == iType .and. jMol == iMol) cycle
        jIndx = MolArray(jType)%mol(jMol)%indx

        do jAtom = 1, nAtoms(jType)
          atmType2 = atomArray(jType, jAtom)
          epAlpha  = epAlpha_tab(atmType2, atmType1)
          sig_sq   = sigsq_tab(atmType2, atmType1)
          Mu       = Mu_tab(atmType2, atmType1)
          rmin_ij  = r_min_tab(atmType2, atmType1)
          q_Nzero  = q_Nonzero(atmType2, atmType1)
          rcutNBsq = 9.0_dp * sig_sq

          ! New distance
          rx = disp(iDisp)%x_new - MolArray(jType)%mol(jMol)%x(jAtom)
          ry = disp(iDisp)%y_new - MolArray(jType)%mol(jMol)%y(jAtom)
          rz = disp(iDisp)%z_new - MolArray(jType)%mol(jMol)%z(jAtom)
          r_new = rx*rx + ry*ry + rz*rz

          if (r_new < rmin_ij) then
            rejMove = .true.
            return
          endif

          ! Neighbor tracking
          if (minDistCriteria) then
            if (r_new < PairList(jIndx)) PairList(jIndx) = r_new
            if (r_new < Dist_Critr_sq .and. present(NeighborDetailsNew)) then
              iPair = NeighborDetailsNew(jIndx)%nPairs + 1
              NeighborDetailsNew(jIndx)%nPairs = iPair
              NeighborDetailsNew(jIndx)%pairIndices(2*iPair - 1) = iAtom
              NeighborDetailsNew(jIndx)%pairIndices(2*iPair)     = jAtom
            endif
          elseif (distCriteria) then
            if (iAtom == 1 .and. jAtom == 1) PairList(jIndx) = r_new
          endif

          if (disp(iDisp)%Displaced) then
            ! Old distance
            rx = disp(iDisp)%x_old - MolArray(jType)%mol(jMol)%x(jAtom)
            ry = disp(iDisp)%y_old - MolArray(jType)%mol(jMol)%y(jAtom)
            rz = disp(iDisp)%z_old - MolArray(jType)%mol(jMol)%z(jAtom)
            r_old = rx*rx + ry*ry + rz*rz

            ! WF energy
            WF = 0.0_dp
            if (r_new < rcutNBsq) WF = WF_Func(r_new, epAlpha, sig_sq, Mu)
            E_WF = E_WF + WF
            dETable(iIndx) = dETable(iIndx) + WF
            dETable(jIndx) = dETable(jIndx) + WF
            if (.not. distCriteria .and. .not. minDistCriteria) PairList(jIndx) = PairList(jIndx) + WF

            WF = 0.0_dp
            if (r_old < rcutNBsq) WF = WF_Func(r_old, epAlpha, sig_sq, Mu)
            E_WF = E_WF - WF
            dETable(iIndx) = dETable(iIndx) - WF
            dETable(jIndx) = dETable(jIndx) - WF

            ! Electrostatics
            if (q_Nzero) then
              q = q_tab(atmType2, atmType1)

              Ele = 0.0_dp
              if (r_new < rcutElecsq) then
                r_new = sqrt(r_new)
                Ele = q * exp(-kapa * r_new) / r_new
              endif
              E_Ele = E_Ele + Ele
              dETable(iIndx) = dETable(iIndx) + Ele
              dETable(jIndx) = dETable(jIndx) + Ele
              if (.not. distCriteria .and. .not. minDistCriteria) PairList(jIndx) = PairList(jIndx) + Ele

              Ele = 0.0_dp
              if (r_old < rcutElecsq) then
                r_old = sqrt(r_old)
                Ele = q * exp(-kapa * r_old) / r_old
              endif
              E_Ele = E_Ele - Ele
              dETable(iIndx) = dETable(iIndx) - Ele
              dETable(jIndx) = dETable(jIndx) - Ele
            end if
          end if
        end do
      end do
    end do
  end do

  ! Optional correction if only a subset of atoms was displaced
  if (.not. distCriteria .and. .not. minDistCriteria) then
    if (sizeDisp < nAtoms(iType)) then
      call Shift_PairList_Correct(disp, PairList)
    end if
  end if

  E_Trial = E_WF + E_Ele
end subroutine Shift_ECalc_Inter
!======================================================================================
! Corrects the PairList array for partial molecule moves in a Grand Canonical Monte Carlo
! nucleation simulation using the Mpipi force field. Adds intermolecular energy contributions
! (Wang-Frenkel and electrostatic) from undisplaced atoms in the moved molecule to PairList,
! used for energy-based cluster criteria when sizeDisp < nAtoms(iType) and neither distCriteria
! nor minDistCriteria is active. Called by Shift_ECalc_Inter for moves like single-atom displacements.
pure subroutine Shift_PairList_Correct(disp, PairList)

  ! Corrects PairList for partial molecule moves by including interactions from
  ! non-displaced atoms with the rest of the system, using the Mpipi force field.

  use VarPrecision, only: dp, atomIntType
  use ForceField, only: atomArray, nAtoms
  use ForceFieldPara_Mpipi, only: epAlpha_tab, sigsq_tab, Mu_tab, q_tab, q_Nonzero, rcutElecsq
  use Coords, only: MolArray, Displacement
  use SimParameters, only: nMolTypes, NPART, kapa

  implicit none

  type(Displacement), intent(in) :: disp(:)
  real(dp), intent(inout) :: PairList(:)

  integer :: iType, jType, iMol, jMol, iAtom, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2, jIndx
  integer :: sizeDisp, Mu
  real(dp) :: rx, ry, rz, r
  real(dp) :: epAlpha, sig_sq, q
  logical :: q_Nzero
  real(dp) :: WF, Ele, rcutNBsq

  sizeDisp = size(disp)
  iType = disp(1)%molType
  iMol  = disp(1)%molIndx

  do iAtom = 1, nAtoms(iType)
    if(any(disp%atmIndx == iAtom)) cycle
    atmType1 = atomArray(iType, iAtom)

    do jType = 1, nMolTypes
      do jAtom = 1, nAtoms(jType)
        atmType2 = atomArray(jType, jAtom)
        epAlpha  = epAlpha_tab(atmType2, atmType1)
        q_Nzero  = q_Nonzero(atmType2, atmType1)
        sig_sq   = sigsq_tab(atmType2, atmType1)
        rcutNBsq = 9.0_dp * sig_sq
        Mu       = Mu_tab(atmType2, atmType1)

        do jMol = 1, NPART(jType)
          if(iType == jType .and. iMol == jMol) cycle
          jIndx = MolArray(jType)%mol(jMol)%indx

          ! Distance for the new position
          rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
          ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
          rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
          r  = rx*rx + ry*ry + rz*rz

          ! Wang–Frenkel interaction
          WF = 0.0_dp
          if(r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
          PairList(jIndx) = PairList(jIndx) + WF

          ! Electrostatics
          if(q_Nzero) then
            if(r < rcutElecsq) then
              q = q_tab(atmType2, atmType1)
              r = sqrt(r)
              Ele = q * exp(-kapa * r) / r
              PairList(jIndx) = PairList(jIndx) + Ele
            endif
          endif

        enddo
      enddo
    enddo
  enddo

end subroutine Shift_PairList_Correct
!======================================================================================      
! Computes the intermolecular energy of a specified molecule with all other molecules in a Grand Canonical
! Monte Carlo simulation using the Mpipi force field with useDistStore=.false. Calculates E_Trial (Wang-Frenkel
! potential and Debye-Hückel electrostatics) for molecule (iType, iMol), updating dETable for cluster criteria.
! Distances are computed on-the-fly. Used in swap-out moves (e.g., SwapOut_ECalc_Mpipi) and other energy calculations.
pure subroutine Mol_ECalc_Inter(iType, iMol, dETable, E_Trial)

  ! Computes the intermolecular energy (Wang-Frenkel + electrostatics)
  ! of molecule (iType, iMol) with all others using the Mpipi force field.
  ! Used in swap-out and trial energy calculations.

  use VarPrecision
  use ForceField, only: nAtoms, atomArray
  use ForceFieldPara_Mpipi, only: epAlpha_tab, q_Nonzero, q_tab, sigsq_tab, Mu_tab, rcutElec, rcutElecsq
  use Coords, only: MolArray
  use SimParameters, only: nMolTypes, NPART, kapa

  implicit none

  integer, intent(in) :: iType, iMol
  real(dp), intent(out) :: E_Trial
  real(dp), intent(inout) :: dETable(:)

  integer :: iAtom, iIndx, jType, jIndx, jMol, jAtom, Mu
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: epAlpha, sig_sq, q
  logical :: q_Nzero
  real(dp) :: WF, Ele
  real(dp) :: E_Ele, E_WF, rcutNBsq, rmax

  E_WF = 0.0_dp
  E_Ele = 0.0_dp
  E_Trial = 0.0_dp
  dETable = 0.0_dp

  iIndx = MolArray(iType)%mol(iMol)%indx

  do iAtom = 1, nAtoms(iType)
    atmType1 = atomArray(iType, iAtom)

    do jType = 1, nMolTypes
      do jAtom = 1, nAtoms(jType)
        atmType2 = atomArray(jType, jAtom)
        epAlpha  = epAlpha_tab(atmType2, atmType1)
        q_Nzero  = q_Nonzero(atmType2, atmType1)
        sig_sq   = sigsq_tab(atmType2, atmType1)
        rcutNBsq = 9.0_dp * sig_sq
        Mu       = Mu_tab(atmType2, atmType1)

        do jMol = 1, NPART(jType)
          if(iType == jType .and. iMol == jMol) cycle
          jIndx = MolArray(jType)%mol(jMol)%indx

          ! Intermolecular distance
          rx = abs(MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom))
          ry = abs(MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom))
          rz = abs(MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom))
          rmax = max(rx, ry, rz)
          if(rmax > rcutElec) cycle
          if(.not. q_Nzero) then
            rmax = rmax * rmax
            if(rmax > rcutNBsq) cycle
          endif

          r = rx*rx + ry*ry + rz*rz
          if(r > rcutElecsq) cycle

          ! Wang-Frenkel
          WF = 0.0_dp
          if(r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
          E_WF = E_WF + WF
          dETable(iIndx) = dETable(iIndx) + WF
          dETable(jIndx) = dETable(jIndx) + WF

          ! Electrostatics
          if(q_Nzero) then
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

  E_Trial = E_WF + E_Ele

end subroutine Mol_ECalc_Inter
!======================================================================================      
! Calculates the intermolecular energy for a newly inserted molecule interacting with the existing
! cluster in a Grand Canonical Monte Carlo nucleation simulation using the Mpipi force field.
! Computes E_Trial (E_WF + E_Ele) for the new molecule, updates PairList (distances for
! distCriteria/minDistCriteria, energies otherwise), and populates NeighborDetailsNew with atom
! pair indices for minDistCriteria when provided. Supports LLAVBMC for long linear biomolecules.
! Rejects the move if atomic overlaps occur or other criteria fail. Compatible with useDistStore=.false.
pure subroutine NewMol_ECalc_Inter(E_Trial, PairList, dETable, rejMove, NeighborDetailsNew)

  ! Computes intermolecular energy of the inserted molecule using the Mpipi force field.
  ! Updates PairList and dETable; supports dist/minDist criteria and NeighborDetailsNew.
  ! Rejects if overlaps occur. Used when useDistStore = .false.

  use ForceField, only: r_min_tab, nAtoms, atomArray
  use ForceFieldPara_Mpipi, only: rcutElec, rcutElecsq, epAlpha_tab, q_Nonzero, sigsq_tab, Mu_tab, q_tab
  use Coords, only: newMol, MolArray
  use SimParameters, only: nMolTypes, NPART, distCriteria, minDistCriteria, Dist_Critr_sq, &
                           kapa
  use CoordinateTypes, only: NeighborDetails
  use VarPrecision

  implicit none

  logical, intent(out) :: rejMove
  real(dp), intent(out) :: E_Trial
  real(dp), intent(inout) :: PairList(:), dETable(:)
  type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)

  integer :: iAtom, iIndx, jType, jIndx, jMol, jAtom, Mu, iPair
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: epAlpha, sig_sq, q
  logical :: q_Nzero
  real(dp) :: WF, Ele
  real(dp) :: E_Ele, E_WF, rcutNBsq, rmax, rmin_ij

  E_WF = 0.0_dp
  E_Ele = 0.0_dp
  E_Trial = 0.0_dp
  dETable = 0.0_dp
  PairList = huge(dp)
  if ((.not. distCriteria) .and. (.not. minDistCriteria)) PairList = 0.0_dp

  rejMove = .false.

  iIndx = MolArray(newMol%molType)%mol(NPART(newMol%molType)+1)%indx

  do iAtom = 1, nAtoms(newMol%molType)
    atmType1 = atomArray(newMol%molType, iAtom)

    do jType = 1, nMolTypes
      do jAtom = 1, nAtoms(jType)
        atmType2 = atomArray(jType, jAtom)
        epAlpha  = epAlpha_tab(atmType2, atmType1)
        q_Nzero  = q_Nonzero(atmType2, atmType1)
        sig_sq   = sigsq_tab(atmType2, atmType1)
        rcutNBsq = 9.0_dp * sig_sq
        Mu       = Mu_tab(atmType2, atmType1)
        rmin_ij  = r_min_tab(atmType2, atmType1)

        do jMol = 1, NPART(jType)
          rx = abs(newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom))
          ry = abs(newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom))
          rz = abs(newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom))
          rmax = max(rx, ry, rz)

          if (rmax > rcutElec) cycle
          if (.not. q_Nzero) then
            rmax = rmax * rmax
            if (rmax > rcutNBsq) cycle
          endif

          r = rx*rx + ry*ry + rz*rz
          if (r > rcutElecsq) cycle

          if (r < rmin_ij) then
            rejMove = .true.
            return
          endif

          jIndx = MolArray(jType)%mol(jMol)%indx

          ! Handle distance-based cluster criteria
          if(minDistCriteria) then
            if(r < PairList(jIndx)) PairList(jIndx) = r
            if(r < Dist_Critr_sq .and. present(NeighborDetailsNew)) then
              iPair = NeighborDetailsNew(jIndx)%nPairs + 1
              NeighborDetailsNew(jIndx)%nPairs = iPair
              NeighborDetailsNew(jIndx)%pairIndices(2 * iPair - 1) = iAtom
              NeighborDetailsNew(jIndx)%pairIndices(2 * iPair)     = jAtom
            endif
          elseif(distCriteria) then
            if(iAtom == 1 .and. jAtom == 1) then
              PairList(jIndx) = r
            endif
          endif

          ! Wang–Frenkel potential
          WF = 0.0_dp
          if (r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
          E_WF = E_WF + WF
          if ((.not. distCriteria) .and. (.not. minDistCriteria)) PairList(jIndx) = PairList(jIndx) + WF
          dETable(jIndx) = dETable(jIndx) + WF
          dETable(iIndx) = dETable(iIndx) + WF

          ! Electrostatics
          if(q_Nzero) then
            q = q_tab(atmType2, atmType1)
            r = sqrt(r)
            Ele = q * exp(-kapa * r) / r
            E_Ele = E_Ele + Ele
            if ((.not. distCriteria) .and. (.not. minDistCriteria)) PairList(jIndx) = PairList(jIndx) + Ele
            dETable(jIndx) = dETable(jIndx) + Ele
            dETable(iIndx) = dETable(iIndx) + Ele
          endif

        enddo
      enddo
    enddo
  enddo

  E_Trial = E_WF + E_Ele

end subroutine NewMol_ECalc_Inter
!======================================================================================      
! Evaluates the intermolecular interaction energy between a newly inserted molecule (newMol)
! and a target molecule (jType, jMol) to check the energy-based cluster criterion in a grand
! canonical ensemble nucleation simulation. Uses the Mpipi force field:
!   U_WF(r) = eps * alpha * [((sigma/r)^(2*mu) - 1) - ((3*sigma/r)^(2*mu) - 1)]^(2*nu)
!   U_Electrostatic(r) = (q1 * q2 / r) * exp(-kappa * r)
! Computes Wang-Frankel (WF) and Debye-Hückel electrostatic (Ele) contributions for all atom
! pairs between the molecules. Rejects the move (sets rejMove=.true.) if:
!   1. Any pair distance is below the minimum distance (r_min_tab).
!   2. The total interaction energy (E_Trial) exceeds the energy criterion (Eng_Critr).
! Used in Monte Carlo moves (e.g., LLAVBMC) when distCriteria=.false. Does not modify
! NeighborList or NeighborPairs, as it only evaluates energy for cluster checks.
! Supports minimum distance criterion for long linear biomolecules via r_min_tab check.
subroutine QuickNei_ECalc_Inter_Mpipi(jType, jMol, rejMove)

  ! Evaluates intermolecular energy between newMol and (jType, jMol) using Mpipi force field.
  ! Rejects move if atoms overlap or total interaction energy exceeds Eng_Critr.

  use VarPrecision, only: dp
  use ForceField, only: nAtoms, atomArray, r_min_tab
  use ForceFieldPara_Mpipi, only: epAlpha_tab, q_Nonzero, q_tab, sigsq_tab, Mu_tab, rcutElecsq
  use Coords, only: MolArray, newMol
  use SimParameters, only: kapa, Eng_Critr

  implicit none

  integer, intent(in) :: jType, jMol
  logical, intent(out) :: rejMove

  integer :: iAtom, jAtom, Mu
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: epAlpha, sig_sq, q
  logical :: q_Nzero
  real(dp) :: WF, Ele
  real(dp) :: E_Trial, E_Ele, E_WF, rcutNBsq
  real(dp) :: rmin_ij

  E_WF = 0.0_dp
  E_Ele = 0.0_dp
  E_Trial = 0.0_dp
  rejMove = .false.

  do iAtom = 1, nAtoms(newMol%molType)
    atmType1 = atomArray(newMol%molType, iAtom)
    do jAtom = 1, nAtoms(jType)
      atmType2 = atomArray(jType, jAtom)
      rmin_ij = r_min_tab(atmType2, atmType1)

      rx = newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
      ry = newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
      rz = newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
      r = rx*rx + ry*ry + rz*rz

      if(r < rmin_ij) then
        rejMove = .true.
        return
      endif

      epAlpha = epAlpha_tab(atmType2, atmType1)
      q_Nzero = q_Nonzero(atmType2, atmType1)
      sig_sq = sigsq_tab(atmType2, atmType1)
      rcutNBsq = 9.0_dp * sig_sq
      Mu = Mu_tab(atmType2, atmType1)

      WF = 0.0_dp
      if (r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
      E_WF = E_WF + WF

      if(q_Nzero) then
        if (r < rcutElecsq) then
          q = q_tab(atmType1, atmType2)
          r = sqrt(r)
          Ele = q * exp(-kapa * r) / r
          E_Ele = E_Ele + Ele
        endif
      endif

    enddo
  enddo

  E_Trial = E_WF + E_Ele

  if(E_Trial > Eng_Critr(newMol%molType, jType)) then
    rejMove = .true.
  endif

end subroutine QuickNei_ECalc_Inter_Mpipi
!======================================================================================
      end module
      
       
