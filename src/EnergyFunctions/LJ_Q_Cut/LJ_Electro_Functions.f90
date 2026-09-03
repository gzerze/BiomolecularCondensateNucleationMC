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
      contains
!======================================================================================      
      subroutine Detailed_ECalc_Inter(E_T, PairList)
      use ParallelVar
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      use EnergyTables
      use PairStorage, only: rPair, distStorage, nTotalAtoms
      implicit none
      real(dp), intent(inOut) :: E_T
      real(dp), intent(inOut) :: PairList(:,:)
      integer :: iType,jType,iMol,jMol,iAtom,jAtom
      integer(kind=atomIntType) :: atmType1,atmType2      
      integer :: iIndx, jIndx, globIndx1, globIndx2, jMolMin
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij      

      E_LJ = 0E0
      E_Ele = 0E0
      E_Inter_T = 0E0
      PairList = 0E0      
      ETable = 0E0
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
                 rmin_ij = r_min_tab(atmType1,atmType2)          
                 globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)

                 rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
                 ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
                 rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom) 
                 r = rx**2 + ry**2 + rz**2
                 if(distCriteria) then
                   if(iAtom .eq. 1) then
                     if(jAtom .eq. 1) then
                       PairList(iIndx, jIndx) = r
                       PairList(jIndx, iIndx) = PairList(iIndx,jIndx)                    
                     endif
                   endif
                 endif
                 if(r .lt. rmin_ij) then
                   stop "ERROR! Overlaping atoms found in the current configuration!"
                 endif 
                 rPair(globIndx1, globIndx2)%p%r_sq = r
                 LJ = (sig_sq/r)**3
                 LJ = ep * LJ * (LJ-1E0)              
                 E_LJ = E_LJ + LJ
              
                 r = sqrt(r)
                 Ele = q/r
                 E_Ele = E_Ele + Ele
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
      
      write(nout,*) "Lennard-Jones Energy:", E_LJ
      write(nout,*) "Eletrostatic Energy:", E_Ele

!      write(35,*) "Pair List:"
!      do iMol=1,maxMol
!        write(35,*) iMol, PairList(iMol)
!      enddo

      do iAtom = 1, size(distStorage) - 1
        write(35,*) distStorage(iAtom)%indx1, distStorage(iAtom)%indx2, distStorage(iAtom)%r_sq, distStorage(iAtom)%E_Pair
      enddo
      
      E_T = E_T + E_Ele + E_LJ    
      E_Inter_T = E_Ele + E_LJ   
      
      end subroutine
!======================================================================================      
      subroutine Shift_ECalc_Inter(E_Trial,disp, PairList,dETable,rejMove)
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      implicit none
      
      type(Displacement), intent(in) :: disp(:)      
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      logical, intent(out) :: rejMove

      
      integer :: iType,jType,iMol,jMol,iAtom,jAtom,iDisp
      integer(kind=atomIntType) :: atmType1,atmType2,iIndx,jIndx
      integer :: sizeDisp 
      integer :: gloIndx1, gloIndx2
      real(dp) :: rx,ry,rz
      real(dp) :: r_new, r_old
      real(dp) :: r_min1_sq      
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij    
      real(dp) :: time_r, time_LJ, time_Ele
      real(dp) :: cnt_r, cnt_LJ, cnt_Ele
      real(dp) :: time1, time2

      sizeDisp = size(disp)
      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      PairList = 0E0      

      dETable = 0E0
!      if(NTotal .eq. 1) return
      iType = disp(1)%molType
      iMol = disp(1)%molIndx
      iIndx = MolArray(iType)%mol(iMol)%indx
      
!      !This section calculates the Intermolecular interaction between the atoms that
!      !have been modified in this trial move with the atoms that have remained stationary

      do iDisp=1,sizeDisp
        iAtom = disp(iDisp)%atmIndx
        atmType1 = atomArray(iType,iAtom)
        gloIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType,jAtom)
            ep = ep_tab(atmType2, atmType1)
            q = q_tab(atmType2, atmType1)
            rmin_ij = r_min_tab(atmType2, atmType1)

!            if(q .eq. 0E0) then
!              if(ep .eq. 0E0) then
!                if(rmin_ij .eq. 0E0) then
!                  cycle
!                endif
!              endif
!            endif
            sig_sq = sig_tab(atmType2,atmType1)
            do jMol=1,NPART(jType)
              if(iType .eq. jType) then
                if(iMol .eq. jMol) then
                  cycle
                endif
              endif  
              gloIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
!               Distance for the New position
              rx = disp(iDisp)%x_new - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = disp(iDisp)%y_new - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = disp(iDisp)%z_new - MolArray(jType)%mol(jMol)%z(jAtom)
              r_new = rx*rx + ry*ry + rz*rz
!             If r_new is less than r_min reject the move.              
              if(r_new .lt. rmin_ij) then
                 rejMove = .true.
                 return
              endif    
              if(distCriteria) then
                if(iAtom .eq. 1) then
                  if(jAtom .eq. 1) then
                    PairList(jIndx) = r_new
                  endif
                endif
              endif

!             Distance for the Old position
              rx = disp(iDisp)%x_old - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = disp(iDisp)%y_old - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = disp(iDisp)%z_old - MolArray(jType)%mol(jMol)%z(jAtom)
              r_old = rx*rx + ry*ry + rz*rz              

              jIndx = MolArray(jType)%mol(jMol)%indx
!             Check to see if there is a non-zero Lennard-Jones parmaeter. If so calculate
!             the Lennard-Jones energy           
              if(ep .ne. 0E0) then
                LJ = (sig_sq/r_new)
                LJ = LJ * LJ * LJ              
                LJ = ep * LJ * (LJ-1E0)                
                E_LJ = E_LJ + LJ
                if(.not. distCriteria) then
                  PairList(jIndx) = PairList(jIndx) + LJ
                endif
                dETable(iIndx) = dETable(iIndx) + LJ
                dETable(jIndx) = dETable(jIndx) + LJ
                
                LJ = (sig_sq/r_old)
                LJ = LJ * LJ * LJ
                LJ = ep * LJ * (LJ-1E0)                
                E_LJ = E_LJ - LJ
                dETable(iIndx) = dETable(iIndx) - LJ
                dETable(jIndx) = dETable(jIndx) - LJ                                
              endif
!             Check to see if there is a non-zero Electrostatic parmaeter. If so calculate
!             the electrostatic energy              
              if(q .ne. 0E0) then
                r_new = sqrt(r_new)
                Ele = q / r_new
                E_Ele = E_Ele + Ele
                if(.not. distCriteria) then                
                  PairList(jIndx) = PairList(jIndx) + Ele
                endif
                dETable(iIndx) = dETable(iIndx) + Ele
                dETable(jIndx) = dETable(jIndx) + Ele
                
                r_old = sqrt(r_old)
                Ele = q / r_old
                E_Ele = E_Ele - Ele
                dETable(iIndx) = dETable(iIndx) - Ele
                dETable(jIndx) = dETable(jIndx) - Ele                                
              endif
            enddo
          enddo
        enddo
      enddo

     
      if(.not. distCriteria) then      
        if(sizeDisp .lt. nAtoms(iType)) then
          call Shift_PairList_Correct(disp, PairList)
        endif
      endif
     
      E_Trial = E_LJ + E_Ele
      
      
      end subroutine
!======================================================================================
      pure subroutine Shift_PairList_Correct(disp, PairList)
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      implicit none
      
      type(Displacement), intent(in) :: disp(:)      
      real(dp), intent(inout) :: PairList(:)
      
      integer :: iType,jType,iMol,jMol,iAtom,jAtom
      integer(kind=atomIntType) :: atmType1,atmType2, jIndx
      integer :: sizeDisp 
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele

      sizeDisp = size(disp)
      iType = disp(1)%molType
      iMol = disp(1)%molIndx

      do iAtom=1,nAtoms(iType)
        if(any(disp%atmIndx .eq. iAtom)) cycle
        atmType1 = atomArray(iType,iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType,jAtom)
            ep = ep_tab(atmType2,atmType1)
            q = q_tab(atmType2,atmType1)
            if(q .eq. 0E0) then
              if(ep .eq. 0E0) then
                cycle
              endif
            endif
            sig_sq = sig_tab(atmType2,atmType1)
            do jMol=1, NPART(jType)
              if(iType .eq. jType) then
                if(iMol .eq. jMol) then
                  cycle
                endif
              endif  
              jIndx = MolArray(jType)%mol(jMol)%indx              
!             Distance for the New position
              rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r = rx*rx + ry*ry + rz*rz

!             Check to see if there is a non-zero Lennard-Jones parmaeter. If so calculate
!             the Lennard-Jones energy           
              if(ep .ne. 0E0) then
                LJ = (sig_sq/r)
                LJ = LJ * LJ * LJ              
                LJ = ep * LJ * (LJ - 1E0)                
                PairList(jIndx) = PairList(jIndx) + LJ
              endif
!             Check to see if there is a non-zero Electrostatic parmaeter. If so calculate
!             the electrostatic energy              
              if(q .ne. 0E0) then
                r = sqrt(r)
                Ele = q / r
                PairList(jIndx) = PairList(jIndx) + Ele
              endif
            enddo
          enddo
        enddo
      enddo
      
      end subroutine      
!======================================================================================      
      pure subroutine Mol_ECalc_Inter(iType, iMol, dETable, E_Trial)
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      implicit none
      integer, intent(in) :: iType, iMol     
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: dETable(:)
      
      integer :: iAtom,iIndx,jType,jIndx,jMol,jAtom
      integer(kind=atomIntType)  :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ

      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      dETable = 0E0
      
      iIndx = MolArray(iType)%mol(iMol)%indx

   
      do iAtom = 1,nAtoms(iType)
        atmType1 = atomArray(iType,iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType,jAtom)
            ep = ep_tab(atmType2,atmType1)
            q = q_tab(atmType2,atmType1)
            if(q .eq. 0E0) then
              if(ep .eq. 0E0) then
                cycle
              endif
            endif
            sig_sq = sig_tab(atmType2,atmType1)

            do jMol=1,NPART(jType)
              if(iType .eq. jType) then
                if(iMol .eq. jMol) then
                  cycle
                endif
              endif
              jIndx = MolArray(jType)%mol(jMol)%indx               
!             New Energy Calculation
              rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r = rx*rx + ry*ry + rz*rz
              if(ep .ne. 0E0) then
                LJ = (sig_sq/r)
                LJ = LJ * LJ * LJ              
                LJ = ep * LJ * (LJ-1E0)                
                E_LJ = E_LJ + LJ
                dETable(iIndx) = dETable(iIndx) + LJ
                dETable(jIndx) = dETable(jIndx) + LJ
              endif
              if(q .ne. 0E0) then            
                r = sqrt(r)
                Ele = q / r
                E_Ele = E_Ele + Ele
                dETable(iIndx) = dETable(iIndx) + Ele
                dETable(jIndx) = dETable(jIndx) + Ele                
              endif
            enddo
          enddo
        enddo
      enddo
     
      E_Trial = E_LJ + E_Ele
      
      
      end subroutine
!======================================================================================      
      pure subroutine NewMol_ECalc_Inter(E_Trial, PairList, dETable, rejMove)
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      implicit none
      logical, intent(out) :: rejMove
      real(dp), intent(out) :: E_Trial

      real(dp), intent(inout) :: PairList(:), dETable(:)
      
      integer :: iAtom, iIndx, jType, jIndx, jMol, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele

      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij

      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      dETable = 0E0
      PairList = 0E0
      rejMove = .false.
      
      iIndx = molArray(newMol%molType)%mol(NPART(newMol%molType)+1)%indx
  
      do iAtom = 1,nAtoms(newMol%molType)
        atmType1 = atomArray(newMol%molType,iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType,jAtom)
            ep = ep_tab(atmType2,atmType1)
            q = q_tab(atmType2,atmType1)
            sig_sq = sig_tab(atmType2,atmType1)
            rmin_ij = r_min_tab(atmType2,atmType1)
            do jMol = 1,NPART(jType)
              rx = newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r = rx*rx + ry*ry + rz*rz
              if(r .lt. rmin_ij) then
                rejMove = .true.
                return
              endif
              jIndx = molArray(jType)%mol(jMol)%indx  
              if(distCriteria) then              
                if(iAtom .eq. 1) then
                  if(jAtom .eq. 1) then
                    PairList(jIndx) = r
                  endif
                endif
              endif              
              if(ep .ne. 0E0) then
                LJ = (sig_sq/r)
                LJ = LJ * LJ * LJ              
                LJ = ep * LJ * (LJ-1E0)                
                E_LJ = E_LJ + LJ
                if(.not. distCriteria) then                
                  PairList(jIndx) = PairList(jIndx) + LJ
                endif
                dETable(jIndx) = dETable(jIndx) + LJ
                dETable(iIndx) = dETable(iIndx) + LJ
              endif
              if(q .ne. 0E0) then
                r = sqrt(r)
                Ele = q / r
                E_Ele = E_Ele + Ele
                if(.not. distCriteria) then                
                  PairList(jIndx) = PairList(jIndx) + Ele
                endif
                dETable(jIndx) = dETable(jIndx) + Ele
                dETable(iIndx) = dETable(iIndx) + Ele
              endif
            enddo
          enddo
        enddo
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

      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      dETable = 0E0
      PairList = 0E0
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
            if(q .eq. 0.0E0) then
              if(ep .eq. 0.0E0) then
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
              if(ep .ne. 0E0) then
                LJ = (sig_sq/r)
                LJ = LJ * LJ * LJ              
                LJ = ep * LJ * (LJ-1E0)                
                E_LJ = E_LJ + LJ
                if(.not. distCriteria) then                
                  PairList(jIndx) = PairList(jIndx) + LJ
                endif
                dETable(jIndx) = dETable(jIndx) + LJ
                dETable(newIndx) = dETable(newIndx) + LJ
              endif
              if(q .ne. 0E0) then
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
            if(q .eq. 0E0) then
              if(ep .eq. 0E0) then
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
              if(ep .ne. 0E0) then
                LJ = (sig_sq/r)
                LJ = LJ * LJ * LJ              
                LJ = ep * LJ * (LJ-1E0)                
                E_LJ = E_LJ - LJ
                dETable(iIndx2) = dETable(iIndx2) - LJ
                dETable(jIndx) = dETable(jIndx) - LJ
              endif
              if(q .ne. 0E0) then            
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
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
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
          sig_sq = sig_tab(atmType2, atmType1)
          ep = ep_tab(atmType2, atmType1)
          q = q_tab(atmType2, atmType1)
          if(ep .ne. 0E0) then
            LJ = (sig_sq/r)
            LJ = LJ * LJ * LJ              
            LJ = ep * LJ * (LJ - 1E0)            
            E_LJ = E_LJ + LJ
          endif
          if(q .ne. 0E0) then            
            r = sqrt(r)
            Ele = q / r
            E_Ele = E_Ele + Ele
          endif
        enddo
      enddo
     
      E_Trial = E_LJ + E_Ele

      if( E_Trial .gt. Eng_Critr(newMol%molType,jType) ) then
        rejMove = .true.
      endif

      
      end subroutine
!======================================================================================
      end module
      
       
