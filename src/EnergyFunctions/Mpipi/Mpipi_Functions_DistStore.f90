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
!*********************************************************************************************************************
      module InterEnergy_Mpipi_DistStore
      use VarPrecision

      contains
!======================================================================================      
      pure function WF_Func(r_sq, epAlpha, sig_sq, Mu) result(WF)
      implicit none
      real(dp), intent(in) :: r_sq, epAlpha, sig_sq
      integer, intent(in) :: Mu
      real(dp) :: WF, x, y, k  
 
      x = sig_sq/r_sq
      if(Mu .eq. 2) then
        x = x*x
        k = 81.0_dp  ! 9^2
      elseif(Mu .eq. 3) then
        x = x*x*x
        k = 729.0_dp  ! 9^3
      else
        x = x**Mu
        k = 9.0_dp ** Mu
      endif
      y = k*x-1E0_dp
      WF = epAlpha * (x-1E0_dp) * y * y  

      end function
!====================================================================================== 
      subroutine Detailed_ECalc_Inter(E_T, PairList)
      use ParallelVar, only: nout
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_Mpipi
      use Coords, only: MolArray, NeighborPairs
      use SimParameters, only: nMolTypes, NPART, distCriteria,minDistCriteria, kapa, Dist_Critr_sq,maxMol
      use EnergyTables, only: ETable, E_Inter_T
      use PairStorage, only: rPair
      implicit none
      real(dp), intent(inOut) :: E_T
      real(dp), intent(inOut) :: PairList(:,:)
      integer :: iType,jType,iMol,jMol,iAtom,jAtom
      integer(kind=atomIntType) :: atmType1,atmType2      
      integer :: iIndx, jIndx, globIndx1, globIndx2, jMolMin, Mu, iPair
      real(dp) :: r
      real(dp) :: epAlpha,sig_sq,q
      logical :: q_Nzero
      real(dp) :: WF, Ele
      real(dp) :: E_Ele,E_WF, rcutNBsq    



      E_WF = 0E0_dp
      E_Ele = 0E0_dp
      E_Inter_T = 0E0_dp
      PairList = 1E7_dp      
      if ((.not. distCriteria) .and. (.not. minDistCriteria)) PairList = 0E0_dp
      ETable = 0E0_dp

      do iType = 1,nMolTypes
        do jType = iType, nMolTypes
          do iMol=1,NPART(iType)
           if(iType .eq. jType) then
             jMolMin = iMol+1
           else
             jMolMin = 1        
           endif
           do jMol = jMolMin,NPART(jType)
             iIndx = MolArray(iType)%mol(iMol)%indx
             jIndx = MolArray(jType)%mol(jMol)%indx  
             do iAtom = 1,nAtoms(iType)
               atmType1 = atomArray(iType,iAtom)
               globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
               do jAtom = 1,nAtoms(jType)        
                 atmType2 = atomArray(jType,jAtom)
                 epAlpha = epAlpha_tab(atmType1,atmType2)
                 q_Nzero = q_Nonzero(atmType1,atmType2)
                 sig_sq = sigsq_tab(atmType1,atmType2)  
                 rcutNBsq = 9E0_dp * sig_sq
                 Mu = Mu_tab(atmType1,atmType2)          
                 globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)

                 r = rPair(globIndx1, globIndx2)%p%r_sq

                 if(minDistCriteria) then
                   if(r .lt. PairList(iIndx, jIndx))  PairList(iIndx, jIndx) = r
                   if(r .lt. PairList(jIndx, iIndx))  PairList(jIndx, iIndx) = r
                   if(r .lt. Dist_Critr_sq) then
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
                   if(iAtom .eq. 1) then
                     if(jAtom .eq. 1) then
                       PairList(iIndx, jIndx) = r
                       PairList(jIndx, iIndx) = PairList(iIndx,jIndx)                    
                     endif
                   endif
                 endif
                 WF = 0E0_dp
                 if (r .lt. rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
                 E_WF = E_WF + WF
                 Ele = 0E0_dp
                 if(q_Nzero) then
                  if (r .lt. rcutElecsq) then
                   q = q_tab(atmType1,atmType2)
                   r = sqrt(r)
                   Ele = q * exp(-kapa * r) / r
                  endif
                 endif
                 E_Ele = E_Ele + Ele
                 rPair(globIndx1, globIndx2)%p%E_Pair = Ele + WF
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
      write(nout,*) "Eletrostatic Energy:", E_Ele

!      write(35,*) "Pair List:"
!      do iMol=1,maxMol
!        write(35,*) iMol, PairList(iMol)
!      enddo


!      do iAtom = 1, size(distStorage) - 1
!        write(35,*) distStorage(iAtom)%indx1, distStorage(iAtom)%indx2, distStorage(iAtom)%r_sq, distStorage(iAtom)%E_Pair
!      enddo
!      flush(35)
      
      E_T = E_T + E_Ele + E_WF    
      E_Inter_T = E_Ele + E_WF   
      
      end subroutine
!======================================================================================      
      pure subroutine Shift_ECalc_Inter(E_Trial,disp,newDist, PairList,dETable,rejMove, NeighborDetailsNew)

      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_Mpipi, only: epAlpha_tab, q_Nonzero, q_tab, sigsq_tab, Mu_tab, rcutElecsq
      use CoordinateTypes, only: NeighborDetails
      use Coords, only: Displacement,  atomIndicies, molArray
      use SimParameters, only: distCriteria, minDistCriteria, kapa, Dist_Critr_sq,maxMol
      use PairStorage, only: distStorage, rPair, DistArrayNew, nNewDist, oldIndxArray

      implicit none
      
      type(Displacement), intent(in) :: disp(:)  
      type(DistArrayNew), intent(inout) :: newDist(:)
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      logical, intent(out) :: rejMove
      type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)
      
      integer :: iType,jType,iMol,jMol,iAtom,jAtom, iPair, kPair
      integer(kind=atomIntType) :: atmType1,atmType2,iIndx,jIndx
      integer :: sizeDisp
!      integer, pointer :: oldIndx 
      integer :: gloIndx1, gloIndx2, Mu
      real(dp) :: r_new, r_new_sq   
      real(dp) :: epAlpha,sig_sq,q
      logical :: q_Nzero
      real(dp) :: WF, Ele, E_New,E_PairOld, E_Old
      real(dp) :: E_Ele,E_WF, rcutNBsq


      E_WF = 0E0_dp
      E_Ele = 0E0_dp
      E_Trial = 0E0_dp
      E_Old = 0E0_dp
      PairList = 1E7_dp      
      if ((.not. distCriteria) .and. (.not. minDistCriteria)) PairList = 0E0_dp
      rejMove = .false.

!      dETable = 0E0
!      if(NTotal .eq. 1) return
      iType = disp(1)%molType
      iMol = disp(1)%molIndx
      iIndx = MolArray(iType)%mol(iMol)%indx
      
!      !This section calculates the Intermolecular interaction between the atoms that
!      !have been modified in this trial move with the atoms that have remained stationary

      do iPair = 1, nNewDist
        gloIndx2 = newDist(iPair)%indx2
        gloIndx1 = newDist(iPair)%indx1
        if(.not. rPair(gloIndx1, gloIndx2)%p%usePair) then
          cycle
        endif
        jType = atomIndicies(gloIndx2)%nType
        jMol  = atomIndicies(gloIndx2)%nMol
        jIndx = MolArray(jType)%mol(jMol)%indx
!        if(jIndx .ne. iIndx) then
!          gloIndx1 = newDist(iPair)%indx1
          jAtom = atomIndicies(gloIndx2)%nAtom
          iAtom = atomIndicies(gloIndx1)%nAtom
       
          atmType1 = atomArray(iType,iAtom)
          atmType2 = atomArray(jType,jAtom)

          epAlpha = epAlpha_tab(atmType2, atmType1)
          q_Nzero = q_Nonzero(atmType2, atmType1)
          sig_sq = sigsq_tab(atmType2, atmType1)  
          rcutNBsq = 9E0_dp * sig_sq
          Mu = Mu_tab(atmType2, atmType1)  
         
          r_new_sq = newDist(iPair)%r_sq
          if(minDistCriteria) then
            if(r_new_sq .lt. PairList(jIndx))  PairList(jIndx) = r_new_sq
            if(r_new_sq .lt. Dist_Critr_sq .and. present(NeighborDetailsNew)) then
              kPair = NeighborDetailsNew(jIndx)%nPairs + 1
              NeighborDetailsNew(jIndx)%nPairs = kPair
              NeighborDetailsNew(jIndx)%pairIndices(2 * kPair - 1) = iAtom
              NeighborDetailsNew(jIndx)%pairIndices(2 * kPair) = jAtom
            endif
          elseif(distCriteria) then
            if(iAtom .eq. 1) then
              if(jAtom .eq. 1) then
                PairList(jIndx) = r_new_sq
              endif
            endif
          endif
          WF = 0E0_dp
          if (r_new_sq .lt. rcutNBsq) WF = WF_Func(r_new_sq, epAlpha, sig_sq, Mu)
          E_WF = E_WF + WF
          Ele = 0E0_dp
          if(q_Nzero) then
           if (r_new_sq .lt. rcutElecsq) then
            q = q_tab(atmType2, atmType1)
            r_new = newDist(iPair)%r
            Ele = q * exp(-kapa * r_new) / r_new
            E_Ele = E_Ele + Ele
           endif
          endif
          E_New = Ele + WF
          if((.not. distCriteria) .and. (.not. minDistCriteria)) then                
            PairList(jIndx) = PairList(jIndx) + E_New 
          endif
          newDist(iPair)%E_Pair = E_New 
          E_PairOld = distStorage(oldIndxArray(iPair))%E_Pair
!          E_PairOld = rPair(gloIndx1, gloIndx2)%p%E_Pair
          dETable(iIndx) = dETable(iIndx) + E_New  - E_PairOld
          dETable(jIndx) = dETable(jIndx) + E_New  - E_PairOld  
          E_Old = E_Old + E_PairOld
!        endif
      enddo




      sizeDisp = size(disp)
      if(.not. distCriteria) then   
       if(.not. minDistCriteria) then      
        if(sizeDisp .lt. nAtoms(iType)) then
          call Shift_PairList_Correct(disp, PairList)
        endif
       endif
      endif
     
      E_Trial = E_WF + E_Ele - E_Old
      
      
      end subroutine
!======================================================================================
      pure subroutine Shift_PairList_Correct(disp, PairList)
      use ForceField, only: nAtoms
      use Coords, only: Displacement, molArray
      use SimParameters, only: nMolTypes, NPART
      use PairStorage, only: rPair, DistArrayNew
      implicit none
      
      type(Displacement), intent(in) :: disp(:)      
      real(dp), intent(inout) :: PairList(:)
      
      integer :: iType,jType,iMol,jMol,iAtom,jAtom
      integer(kind=atomIntType) :: jIndx, iIndx
      integer :: sizeDisp 
      integer :: gloIndx1, gloIndx2
      real(dp) :: E_Pair

      sizeDisp = size(disp)
      iType = disp(1)%molType
      iMol = disp(1)%molIndx
      iIndx = molArray(iType)%mol(iMol)%indx

      do iAtom=1,nAtoms(iType)
        if(any(disp%atmIndx .eq. iAtom)) cycle
        gloIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            do jMol=1, NPART(jType)
              jIndx = MolArray(jType)%mol(jMol)%indx
              if(iIndx .ne. jIndx) then
                gloIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
                E_Pair = rPair(gloIndx1, gloIndx2) % p % E_Pair    
                PairList(jIndx) = PairList(jIndx) + E_Pair
              endif
            enddo
          enddo
        enddo
      enddo

      end subroutine      
!======================================================================================      
      pure subroutine Mol_ECalc_Inter(iType, iMol, dETable, E_Trial)
      use ForceField, only: nAtoms
      use Coords, only: Displacement, molArray
      use SimParameters, only: nMolTypes, NPART
      use PairStorage, only: rPair
      implicit none
      integer, intent(in) :: iType, iMol     
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: dETable(:)
      
      integer :: iAtom,iIndx,jType,jIndx,jMol,jAtom
      integer  :: gloIndx1, gloIndx2
      real(dp) :: E_Pair

      E_Trial = 0E0_dp
      dETable = 0E0_dp
      iIndx = MolArray(iType)%mol(iMol)%indx

      do iAtom = 1,nAtoms(iType)
        gloIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom) 
        do jType = 1, nMolTypes
          do jMol=1, NPART(jType)
            jIndx = MolArray(jType)%mol(jMol)%indx  
            if(iIndx .eq. jIndx) then
              cycle
            endif
            do jAtom = 1,nAtoms(jType)        
              gloIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
              if(.not. rPair(gloIndx1, gloIndx2)%p%usePair) then
                cycle
              endif 
              E_Pair = rPair(gloIndx1, gloIndx2)%p%E_Pair
              E_Trial = E_Trial + E_Pair
              dETable(iIndx) = dETable(iIndx) + E_Pair
              dETable(jIndx) = dETable(jIndx) + E_Pair
            enddo
          enddo
        enddo
      enddo

  
      
      end subroutine
!======================================================================================      
      subroutine NewMol_ECalc_Inter(E_Trial, PairList, dETable, rejMove, NeighborDetailsNew)
      use ForceField, only: atomArray
      use ForceFieldPara_Mpipi, only: epAlpha_tab, q_Nonzero, q_tab, sigsq_tab, Mu_tab, rcutElecsq
      use CoordinateTypes, only: NeighborDetails
      use Coords, only: Displacement,  atomIndicies, molArray, newmol
      use SimParameters, only: distCriteria, NPART, minDistCriteria, kapa, Dist_Critr_sq,maxMol
      use PairStorage, only: rPair, newDist, nNewDist
      implicit none
      logical, intent(out) :: rejMove
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)
      
      integer :: iPair
      integer :: iType, iMol, iAtom, iIndx, jType, jIndx, jMol, jAtom, Mu, kPair
      integer(kind=atomIntType) :: atmType1,atmType2
      integer :: gloIndx1, gloIndx2
      real(dp) :: r_sq, r
      real(dp) :: epAlpha,sig_sq,q
      logical :: q_Nzero
      real(dp) :: WF, Ele
      real(dp) :: E_Ele,E_WF, rcutNBsq


      E_WF = 0E0_dp
      E_Ele = 0E0_dp
      E_Trial = 0E0_dp
      dETable = 0E0_dp
      PairList = 1E7      
      if ((.not. distCriteria) .and. (.not. minDistCriteria)) PairList = 0E0

      rejMove = .false.
      
      iType = newMol%molType
      iMol = NPART(iType)+1
      iIndx = molArray(iType)%mol(iMol)%indx
      do iPair = 1, nNewDist
        gloIndx1 = newDist(iPair)%indx1
        gloIndx2 = newDist(iPair)%indx2
        if(.not. rPair(gloIndx1, gloIndx2)%p%usePair) then
          cycle
        endif
        jType = atomIndicies(gloIndx2)%nType
        jMol = atomIndicies(gloIndx2)%nMol
        jIndx = MolArray(jType)%mol(jMol)%indx
        if(iIndx .ne. jIndx) then
          iAtom = atomIndicies(gloIndx1)%nAtom
          jAtom = atomIndicies(gloIndx2)%nAtom

          atmType1 = atomArray(iType, iAtom)
          atmType2 = atomArray(jType, jAtom)

          epAlpha = epAlpha_tab(atmType2, atmType1)
          q_Nzero = q_Nonzero(atmType2, atmType1)
          sig_sq = sigsq_tab(atmType2, atmType1)  
          rcutNBsq = 9E0_dp * sig_sq
          Mu = Mu_tab(atmType2, atmType1)  

          WF = 0E0_dp
          r_sq = newDist(iPair)%r_sq
          if (r_sq .lt. rcutNBsq) WF = WF_Func(r_sq, epAlpha, sig_sq, Mu)
          E_WF = E_WF + WF
          if(minDistCriteria) then
            if(r_sq .lt. PairList(jIndx))  PairList(jIndx) = r_sq
            if(r_sq .lt. Dist_Critr_sq .and. present(NeighborDetailsNew)) then
              kPair = NeighborDetailsNew(jIndx)%nPairs + 1
              NeighborDetailsNew(jIndx)%nPairs = kPair
              NeighborDetailsNew(jIndx)%pairIndices(2 * kPair - 1) = iAtom
              NeighborDetailsNew(jIndx)%pairIndices(2 * kPair) = jAtom
            endif
          elseif(distCriteria) then
            if(iAtom .eq. 1) then
              if(jAtom .eq. 1) then
                PairList(jIndx) = r_sq
              endif
            endif
          else
            PairList(jIndx) = PairList(jIndx) + WF
          endif

          Ele = 0E0_dp
          if(q_Nzero) then
           if (r_sq .lt. rcutElecsq) then
            q = q_tab(atmType2, atmType1)
            r = newDist(iPair)%r
            Ele = q * exp(-kapa * r) / r
            E_Ele = E_Ele + Ele
            if((.not. distCriteria) .and. (.not. minDistCriteria)) then 
              PairList(jIndx) = PairList(jIndx) + Ele
            endif
           endif
          endif
          dETable(iIndx) = dETable(iIndx) + WF + Ele
          dETable(jIndx) = dETable(jIndx) + WF + Ele
          newDist(iPair)%E_Pair = WF + Ele
        endif 
      enddo
     
      E_Trial = E_WF + E_Ele
      
      
      end subroutine    
!======================================================================================      
      pure subroutine Exchange_ECalc_Inter(E_Trial, nType, nMol, PairList, dETable, rejMove)
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      implicit none
      logical, intent(out) :: rejMove
      integer, intent(in) :: nType, nMol
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      
      integer :: iAtom, newIndx, jType, jIndx, jMol, jAtom
      integer :: iIndx2
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij

      E_LJ = 0E0_dp
      E_Ele = 0E0_dp 
      E_Trial = 0E0_dp
      dETable = 0E0_dp
      PairList = 0E0_dp
      rejMove = .false.
      
      newIndx = molArray(newMol%molType)%mol(NPART(newMol%molType)+1)%indx
      iIndx2 = molArray(nType)%mol(nMol)%indx

       !Calculate the energy of the molecule that is entering the cluster

      do iAtom = 1,nAtoms(newMol%molType)
        atmType1 = atomArray(newMol%molType,iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType,jAtom)
            ep = ep_tab(atmType2,atmType1)
            q = q_tab(atmType2,atmType1)
            if(q .eq. 0.0E0_dp) then
              if(ep .eq. 0.0E0_dp) then
                cycle
              endif
            endif
            sig_sq = sig_tab(atmType2,atmType1)
            rmin_ij = r_min_tab(atmType2,atmType1)
            do jMol = 1,NPART(jType)
              if(jMol .eq. nMol) then
                if(nType .eq. jType) then
                  cycle
                endif
              endif
              jIndx = molArray(jType)%mol(jMol)%indx              
              
              rx = newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r = rx*rx + ry*ry + rz*rz
              if(r .lt. rmin_ij) then
                rejMove = .true.
                return
              endif
              if(distCriteria) then              
                if(iAtom .eq. 1) then
                  if(jAtom .eq. 1) then
                    PairList(jIndx) = r
                  endif
                endif
              endif              
              LJ = 0E0
              Ele = 0E0
              if(ep .ne. 0E0_dp) then
                LJ = (sig_sq/r)
                LJ = LJ * LJ * LJ              
                LJ = ep * LJ * (LJ-1E0_dp)                
                E_LJ = E_LJ + LJ
                if(.not. distCriteria) then                
                  PairList(jIndx) = PairList(jIndx) + LJ
                endif
                dETable(jIndx) = dETable(jIndx) + LJ
                dETable(newIndx) = dETable(newIndx) + LJ
              endif
              if(q .ne. 0E0_dp) then
                r = sqrt(r)
                Ele = q / r
                E_Ele = E_Ele + Ele
                if(.not. distCriteria) then                
                  PairList(jIndx) = PairList(jIndx) + Ele
                endif
                dETable(jIndx) = dETable(jIndx) + Ele
                dETable(newIndx) = dETable(newIndx) + Ele
              endif
            enddo
          enddo
        enddo
      enddo

       !Calculate the energy of the molecule that is exiting the cluster
   
      do iAtom = 1,nAtoms(nType)
        atmType1 = atomArray(nType, iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType,jAtom)
            ep = ep_tab(atmType2,atmType1)
            q = q_tab(atmType2,atmType1)
            if(q .eq. 0E0_dp) then
              if(ep .eq. 0E0_dp) then
                cycle
              endif
            endif
            sig_sq = sig_tab(atmType2,atmType1)
            do jMol=1,NPART(jType)
              if(nMol .eq. jMol) then
                if(nType .eq. jType) then
                  cycle
                endif
              endif
              jIndx = MolArray(jType)%mol(jMol)%indx               
              rx = MolArray(nType)%mol(nMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = MolArray(nType)%mol(nMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = MolArray(nType)%mol(nMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r = rx*rx + ry*ry + rz*rz
              if(ep .ne. 0E0_dp) then
                LJ = (sig_sq/r)
                LJ = LJ * LJ * LJ              
                LJ = ep * LJ * (LJ-1E0_dp)                
                E_LJ = E_LJ - LJ
                dETable(iIndx2) = dETable(iIndx2) - LJ
                dETable(jIndx) = dETable(jIndx) - LJ
              endif
              if(q .ne. 0E0_dp) then            
                r = sqrt(r)
                Ele = q / r
                E_Ele = E_Ele - Ele
                dETable(iIndx2) = dETable(iIndx2) - Ele
                dETable(jIndx) = dETable(jIndx) - Ele                
              endif
            enddo
          enddo
        enddo
      enddo
     

     
      E_Trial = E_LJ + E_Ele
      
      
      end subroutine    
!======================================================================================      
      subroutine QuickNei_ECalc_Inter_Mpipi(jType, jMol, rejMove)
      use ForceField, only: r_min_tab, atomArray, nAtoms
      use ForceFieldPara_Mpipi, only: epAlpha_tab, q_Nonzero, sigsq_tab, Mu_tab, q_tab, rcutElecsq
      use Coords, only: newMol, molArray
      use SimParameters, only: Eng_Critr, kapa
      implicit none
      integer, intent(in) :: jType, jMol     
      logical, intent(out) :: rejMove
      
      integer :: iAtom,jAtom,Mu
      integer(kind=atomIntType)  :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: epAlpha,sig_sq,q
      logical :: q_Nzero
      real(dp) :: WF, Ele
      real(dp) :: E_Trial,E_Ele,E_WF, rcutNBsq
      real(dp) :: rmin_ij

      E_WF = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      rejMove = .false.
    
      do iAtom = 1,nAtoms(newMol%molType)
        atmType1 = atomArray(newMol%molType, iAtom)
        do jAtom = 1,nAtoms(jType)        
          atmType2 = atomArray(jType, jAtom)
          rmin_ij = r_min_tab(atmType2, atmType1)
          rx = newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
          ry = newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
          rz = newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
          r = rx*rx + ry*ry + rz*rz

          if(r .lt. rmin_ij) then
            rejMove = .true.
            return
          endif          
          epAlpha = epAlpha_tab(atmType2, atmType1)
          q_Nzero = q_Nonzero(atmType2, atmType1)
          sig_sq = sigsq_tab(atmType2, atmType1)  
          rcutNBsq = 9E0 * sig_sq
          Mu = Mu_tab(atmType2, atmType1) 
          WF = 0E0
          if (r .lt. rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)  
          E_WF = E_WF + WF

          if(q_Nzero) then    
           if (r .lt. rcutElecsq) then
            q = q_tab(atmType1,atmType2)       
            r = sqrt(r)
            Ele = q * exp(-kapa * r) / r
            E_Ele = E_Ele + Ele
           endif
          endif
        enddo
      enddo
     
      E_Trial = E_WF + E_Ele

      if( E_Trial .gt. Eng_Critr(newMol%molType,jType) ) then
        rejMove = .true.
      endif
!      write(35,*) newMol%molType, jType, E_Trial, Eng_Critr(newMol%molType,jType), rejMove

      
      end subroutine
!======================================================================================
      end module
      
       
