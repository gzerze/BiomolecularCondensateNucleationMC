!*********************************************************************************************************************
!     This file contains the energy functions that work for Lennard-Jones w/ Columbic style forcefields
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
      module InterEnergy_LJ_Electro
      use VarPrecision

      real(dp), parameter :: lj_Cut = 12.0
      real(dp), parameter :: lj_Cut_sq = lj_Cut**2

      real(dp), parameter :: q_Cut = 24.0
      real(dp), parameter :: q_Cut_sq = q_Cut**2
      contains
!======================================================================================      
      pure function LJ_Func(r_sq, ep, sig) result(LJ)
      implicit none
      real(dp), intent(in) :: r_sq, ep, sig
      real(dp) :: LJ  
 
      LJ = (sig/r_sq)
      LJ = LJ * LJ * LJ
      LJ = ep * LJ * (LJ-1E0_dp)  

      end function
!======================================================================================      
      pure function Ele_Func(r, q) result(Ele)
      implicit none
      real(dp), intent(in) :: r, q
      real(dp) :: Ele
 
!      r = sqrt(r_sq)
      Ele = q/r

      end function
!======================================================================================      
      subroutine Detailed_ECalc_Inter(E_T, PairList)
!      use ParallelVar
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_LJ_Q, only: ep_tab, q_tab, sig_tab
      use Coords, only: MolArray
      use SimParameters, only: nMolTypes, NPART, distCriteria
      use EnergyTables, only: ETable, E_Inter_T
      use PairStorage, only: rPair
      implicit none
      real(dp), intent(inOut) :: E_T
      real(dp), intent(inOut) :: PairList(:,:)
      integer :: iType,jType,iMol,jMol,iAtom,jAtom
      integer(kind=atomIntType) :: atmType1,atmType2      
      integer :: iIndx, jIndx, globIndx1, globIndx2, jMolMin
      real(dp) :: r_sq, r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ    



      E_LJ = 0E0_dp
      E_Ele = 0E0_dp
      E_Inter_T = 0E0_dp
      PairList = 0E0_dp
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
                 ep = ep_tab(atmType1,atmType2)
                 q = q_tab(atmType1,atmType2)
                 sig_sq = sig_tab(atmType1,atmType2)          
                 globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)

                 r_sq = rPair(globIndx1, globIndx2)%p%r_sq
                 if(distCriteria) then
                   if(iAtom .eq. 1) then
                     if(jAtom .eq. 1) then
                       PairList(iIndx, jIndx) = r_sq
                       PairList(jIndx, iIndx) = PairList(iIndx,jIndx)                    
                     endif
                   endif
                 endif
                 if(r_sq .lt. lj_cut_sq) then
                   LJ = LJ_Func(r_sq, ep, sig_sq)             
                   E_LJ = E_LJ + LJ
                 endif                   
!                 r = rPair(globIndx1, globIndx2)%p%r
                 if(r_sq .lt. q_cut_sq) then
                   r = sqrt(r_sq)
                   Ele = q / r
!                   Ele = Ele_Func(r, q)
                   E_Ele = E_Ele + Ele
                 endif
                 rPair(globIndx1, globIndx2)%p%E_Pair = Ele + LJ
                 if(.not. distCriteria) then
                   PairList(iIndx, jIndx) = PairList(iIndx, jIndx) + Ele + LJ
                   PairList(jIndx, iIndx) = PairList(iIndx, jIndx)
                 endif
                 ETable(iIndx) = ETable(iIndx) + Ele + LJ
                 ETable(jIndx) = ETable(jIndx) + Ele + LJ              
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      
!      write(nout,*) "Lennard-Jones Energy:", E_LJ
!      write(nout,*) "Eletrostatic Energy:", E_Ele

!      write(35,*) "Pair List:"
!      do iMol=1,maxMol
!        write(35,*) iMol, PairList(iMol)
!      enddo


!      do iAtom = 1, size(distStorage) - 1
!        write(35,*) distStorage(iAtom)%indx1, distStorage(iAtom)%indx2, distStorage(iAtom)%r_sq, distStorage(iAtom)%E_Pair
!      enddo
!      flush(35)
      
      E_T = E_T + E_Ele + E_LJ    
      E_Inter_T = E_Ele + E_LJ   
      
      end subroutine
!======================================================================================      
      pure subroutine Shift_ECalc_Inter(E_Trial,disp,newDist, PairList,dETable,rejMove)

      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_LJ_Q, only: ep_tab, q_tab, sig_tab
      use Coords, only: Displacement,  atomIndicies, molArray
      use SimParameters, only: distCriteria
      use PairStorage, only: distStorage, rPair, DistArrayNew, nNewDist, oldIndxArray

      implicit none
      
      type(Displacement), intent(in) :: disp(:)  
      type(DistArrayNew), intent(inout) :: newDist(:)
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      logical, intent(out) :: rejMove
      
      integer :: iType,jType,iMol,jMol,iAtom,jAtom, iPair
      integer(kind=atomIntType) :: atmType1,atmType2,iIndx,jIndx
      integer :: sizeDisp
!      integer, pointer :: oldIndx 
      integer :: gloIndx1, gloIndx2
      real(dp) :: r_new, r_new_sq   
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele, E_New,E_PairOld, E_Old
      real(dp) :: E_Ele,E_LJ


      E_LJ = 0E0_dp
      E_Ele = 0E0_dp
      E_Trial = 0E0_dp
      E_Old = 0E0_dp
      PairList = 0E0_dp
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
        jType = atomIndicies(gloIndx2)%nType
        jMol  = atomIndicies(gloIndx2)%nMol
        jIndx = MolArray(jType)%mol(jMol)%indx
!        if(jIndx .ne. iIndx) then
          gloIndx1 = newDist(iPair)%indx1
          jAtom = atomIndicies(gloIndx2)%nAtom
          iAtom = atomIndicies(gloIndx1)%nAtom
       
          atmType1 = atomArray(iType,iAtom)
          atmType2 = atomArray(jType,jAtom)

          ep = ep_tab(atmType2, atmType1)
          q  = q_tab(atmType2, atmType1)
         
!          if(distCriteria) then
!            if(iAtom .eq. 1) then
!              if(jAtom .eq. 1) then
!                PairList(jIndx) = newDist(iPair)%r_sq
!              endif
!            endif
!          endif

!          LJ = 0E0
!          Ele = 0E0
          LJ = 0E0_dp
          if(ep .ne. 0E0_dp) then
            r_new_sq = newDist(iPair)%r_sq
            if(r_new_sq .lt. lj_cut_sq) then
              sig_sq = sig_tab(atmType2,atmType1)
              LJ = LJ_Func(r_new_sq, ep, sig_sq)             
              E_LJ = E_LJ + LJ
            endif
          endif
          if(q .ne. 0E0_dp) then
            r_new = newDist(iPair)%r
            Ele = q/r_new
!            Ele = Ele_Func(r_new, q)                
            E_Ele = E_Ele + Ele

!            dETable(iIndx) = dETable(iIndx) + Ele
!            dETable(jIndx) = dETable(jIndx) + Ele
!            newDist(iPair)%E_Pair = newDist(iPair)%E_Pair + Ele
          else
            Ele = 0E0_dp
          endif
          E_New = Ele + LJ
          if(.not. distCriteria) then                
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
        if(sizeDisp .lt. nAtoms(iType)) then
          call Shift_PairList_Correct(disp, PairList)
        endif
      endif
     
      E_Trial = E_LJ + E_Ele - E_Old
      
      
      end subroutine
!======================================================================================
      pure subroutine Shift_PairList_Correct(disp, PairList)
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_LJ_Q, only: ep_tab, q_tab, sig_tab
      use Coords, only: Displacement,  atomIndicies, molArray
      use SimParameters, only: distCriteria, nMolTypes, NPART
      use PairStorage, only: distStorage, rPair, DistArrayNew, nNewDist, oldIndxArray
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
      use ForceField, only: nAtoms, atomArray
      use Coords, only: Displacement,  atomIndicies, molArray
      use SimParameters, only: distCriteria, nMolTypes, NPART
      use PairStorage, only: distStorage, rPair
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
      subroutine NewMol_ECalc_Inter(E_Trial, PairList, dETable, rejMove)
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_LJ_Q, only: ep_tab, q_tab, sig_tab
      use Coords, only: Displacement,  atomIndicies, molArray, newmol
      use SimParameters, only: distCriteria, nMolTypes, NPART
      use PairStorage, only: distStorage, rPair, newDist, nNewDist
      implicit none
      logical, intent(out) :: rejMove
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      
      integer :: iPair
      integer :: iType, iMol, iAtom, iIndx, jType, jIndx, jMol, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      integer :: gloIndx1, gloIndx2
      real(dp) :: r_sq, r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele

      real(dp) :: E_Ele,E_LJ

      E_LJ = 0E0_dp
      E_Ele = 0E0_dp
      E_Trial = 0E0_dp
      dETable = 0E0_dp
      PairList = 0E0_dp
      rejMove = .false.
      
      iType = newMol%molType
      iMol = NPART(iType)+1
      iIndx = molArray(iType)%mol(iMol)%indx
      do iPair = 1, nNewDist
        gloIndx1 = newDist(iPair)%indx1
        gloIndx2 = newDist(iPair)%indx2

        jType = atomIndicies(gloIndx2)%nType
        jMol = atomIndicies(gloIndx2)%nMol
        jIndx = MolArray(jType)%mol(jMol)%indx
        if(iIndx .ne. jIndx) then
          iAtom = atomIndicies(gloIndx1)%nAtom
          jAtom = atomIndicies(gloIndx2)%nAtom

          atmType1 = atomArray(iType, iAtom)
          atmType2 = atomArray(jType, jAtom)

          ep = ep_tab(atmType2, atmType1)
          q = q_tab(atmType2, atmType1)

!          LJ = 0d0
!          Ele = 0d0
          LJ = 0E0_dp
          if(ep .ne. 0E0_dp) then
            r_sq = newDist(iPair)%r_sq
            if(r_sq .lt. lj_cut_sq) then
              sig_sq = sig_tab(atmType2,atmType1)
              LJ = LJ_Func(r_sq, ep, sig_sq)             
              E_LJ = E_LJ + LJ
              if(.not. distCriteria) then
                PairList(jIndx) = PairList(jIndx) + LJ
              endif
            endif
          endif

          if(q .ne. 0E0_dp) then
            r = newDist(iPair)%r
            Ele = q/r
!            Ele = Ele_Func(r_sq, q)                
            E_Ele = E_Ele + Ele
            if(.not. distCriteria) then                
              PairList(jIndx) = PairList(jIndx) + Ele
            endif
          else
            Ele = 0E0_dp
          endif
          dETable(iIndx) = dETable(iIndx) + LJ + Ele
          dETable(jIndx) = dETable(jIndx) + LJ + Ele
          newDist(iPair)%E_Pair = LJ + Ele
        endif 
      enddo
     
      E_Trial = E_LJ + E_Ele
      
      
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
      subroutine QuickNei_ECalc_Inter_LJ_Q(jType, jMol, rejMove)
      use ForceField, only: r_min_tab, atomArray, nAtoms
      use ForceFieldPara_LJ_Q, only: ep_tab, q_tab, sig_tab
      use Coords, only: newMol, molArray
      use SimParameters, only: Eng_Critr
      implicit none
      integer, intent(in) :: jType, jMol     
      logical, intent(out) :: rejMove
      
      integer :: iAtom,jAtom
      integer(kind=atomIntType)  :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Trial,E_Ele,E_LJ
      real(dp) :: rmin_ij

      E_LJ = 0E0
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
          ep = ep_tab(atmType2, atmType1)
          q = q_tab(atmType2, atmType1)
          if(ep .ne. 0E0) then
            sig_sq = sig_tab(atmType2, atmType1)
            LJ = LJ_Func(r, ep, sig_sq)
            E_LJ = E_LJ + LJ
          endif
          if(q .ne. 0E0) then   
            r = sqrt(r)         
            Ele = Ele_Func(r, q)
            E_Ele = E_Ele + Ele
          endif
        enddo
      enddo
     
      E_Trial = E_LJ + E_Ele

      if( E_Trial .gt. Eng_Critr(newMol%molType,jType) ) then
        rejMove = .true.
      endif
!      write(35,*) newMol%molType, jType, E_Trial, Eng_Critr(newMol%molType,jType), rejMove

      
      end subroutine
!======================================================================================
      end module
      
       
