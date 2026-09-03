      module Rosenbluth_Functions_LJ_Q
      use VarPrecision
      use InterEnergy_LJ_Electro, only: lj_cut_sq, q_cut_sq
      contains
!======================================================================================================
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which regrows an entire molecule for the given trial.
      pure subroutine Rosen_Mol_New_LJ_Q(nRosen, nType, included,  E_Trial, overlap)
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_LJ_Q, only: q_tab, ep_tab, sig_tab
      use Coords, only: rosenTrial, MolArray
      use SimParameters, only: nMolTypes, NPART
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nRosen      
      
      logical, intent(out) :: overlap
      real(dp), intent(out) :: E_Trial
      
      integer :: iAtom, jType, jIndx, jMol, jAtom
      integer(kind = atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij

      E_Trial = 0E0
      overlap = .false.
      
      E_LJ = 0E0
      E_Ele = 0E0      

      do iAtom = 1,nAtoms(nType)
        atmType1 = atomArray(nType, iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType, jAtom)
            ep = ep_tab(atmType2, atmType1)
            q = q_tab(atmType2, atmType1)
            sig_sq = sig_tab(atmType2, atmType1)
            rmin_ij = r_min_tab(atmType2, atmType1)
            do jMol = 1,NPART(jType)
              jIndx = molArray(jType)%mol(jMol)%indx              
              if(included(jIndx) .eqv. .false.) then
                cycle
              endif              

              rx = rosenTrial(nRosen)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = rosenTrial(nRosen)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = rosenTrial(nRosen)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r = rx*rx + ry*ry + rz*rz
              if(r .lt. rmin_ij) then
                E_Trial = huge(dp)
                overlap = .true.
                return
              endif

              if(ep .ne. 0E0) then
                if(r .lt. lj_cut_sq) then
                  LJ = (sig_sq / r)
                  LJ = LJ * LJ * LJ              
                  LJ = ep * LJ * (LJ - 1E0)                
                  E_LJ = E_LJ + LJ
                endif
              endif
              if(q .ne. 0E0) then
!                Ele = 0E0
                r = sqrt(r)
                Ele = q / r
                E_Ele = E_Ele + Ele
              endif
            enddo
          enddo
        enddo
      enddo
!      write(*,*) E_Trial
      E_Trial = E_LJ + E_Ele
      
      end subroutine 
!======================================================================================================
      pure subroutine Rosen_Mol_Old_LJ_Q(mol_x, mol_y, mol_z, nType, included,  E_Trial)
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_LJ_Q, only: ep_tab, sig_tab, q_tab
      use Coords, only: MolArray
      use SimParameters, only: NPART, nMolTypes
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType
      real(dp), intent(in) :: mol_x(:), mol_y(:), mol_z(:)
      real(dp), intent(out) :: E_Trial
      
      integer :: iAtom, jType, jIndx, jMol, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij

      E_Trial = 0E0
      E_LJ = 0E0
      E_Ele = 0E0      


      do iAtom = 1,nAtoms(nType)
        atmType1 = atomArray(nType, iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType, jAtom)
            ep = ep_tab(atmType2, atmType1)
            q = q_tab(atmType2, atmType1)
            rmin_ij = r_min_tab(atmType2, atmType1)
            sig_sq = sig_tab(atmType2, atmType1)

            do jMol = 1,NPART(jType)
              jIndx = molArray(jType)%mol(jMol)%indx              
              if(included(jIndx) .eqv. .false.) then
                cycle
              endif         
              rx = mol_x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = mol_y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = mol_z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r = rx*rx + ry*ry + rz*rz
              if(r .lt. rmin_ij) then
                E_Trial = huge(dp)
                return
              endif
              if(ep .ne. 0E0) then
                if(r .lt. lj_cut_sq) then
                  LJ = (sig_sq / r)
                  LJ = LJ * LJ * LJ              
                  LJ = ep * LJ * (LJ - 1E0)                
                  E_LJ = E_LJ + LJ
                endif
              endif
              if(q .ne. 0E0) then
!                Ele = 0E0
                r = sqrt(r)
                Ele = q / r
                E_Ele = E_Ele + Ele
              endif
            enddo
          enddo
        enddo
      enddo
     
      E_Trial = E_LJ + E_Ele
      
      end subroutine 
!======================================================================================================
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which each atom is regrown sequentially for the given trial.
      pure subroutine Rosen_Atom_New_LJ_Q(nType, nAtom, trialPos, included,  E_Trial, overlap)
      use Coords, only: MolArray, SimpleAtomCoords
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_LJ_Q, only: q_tab, sig_tab, ep_tab
      use SimParameters, only: nMolTypes, NPART
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      
      logical, intent(inout) :: overlap
      real(dp), intent(out) :: E_Trial
      
      integer :: jType, jIndx, jMol, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij

      
      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      overlap = .false.

      atmType1 = atomArray(nType, nAtom)
      do jType = 1, nMolTypes
        do jMol = 1,NPART(jType)
          jIndx = molArray(jType)%mol(jMol)%indx              
          if(included(jIndx) .eqv. .false.) then
            cycle
          endif      
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType, jAtom)
            rmin_ij = r_min_tab(atmType2, atmType1)
            rx = trialPos%x - MolArray(jType)%mol(jMol)%x(jAtom)
            ry = trialPos%y - MolArray(jType)%mol(jMol)%y(jAtom)
            rz = trialPos%z - MolArray(jType)%mol(jMol)%z(jAtom)
            r = rx*rx + ry*ry + rz*rz
            if(r .lt. rmin_ij) then
              E_Trial = huge(dp)
              overlap = .true.
              return
            endif
            ep = ep_tab(atmType2, atmType1)
            q = q_tab(atmType2, atmType1)
            sig_sq = sig_tab(atmType2, atmType1)   
            if(ep .ne. 0E0) then
              if(r .lt. lj_cut_sq) then
                LJ = (sig_sq / r)
                LJ = LJ * LJ * LJ              
                LJ = ep * LJ * (LJ - 1E0)                
                E_LJ = E_LJ + LJ
              endif
            endif
            if(q .ne. 0E0) then
!              Ele = 0E0
              r = sqrt(r)
              Ele = q / r
              E_Ele = E_Ele + Ele
            endif
          enddo
        enddo
      enddo
     
      E_Trial = E_LJ + E_Ele
      
      end subroutine 
!======================================================================================================
      pure subroutine Rosen_Atom_Old_LJ_Q(nType, nMol, nAtom, included, E_Trial)
      use Coords, only: MolArray
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_LJ_Q, only: q_tab, ep_tab, sig_tab
      use SimParameters, only: nMolTypes, NPART
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nAtom, nMol
      real(dp), intent(out) :: E_Trial
      
      integer :: jType, jIndx, jMol, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij

      
      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0

      atmType1 = atomArray(nType, nAtom)
      do jType = 1, nMolTypes
        do jMol = 1,NPART(jType)
          jIndx = molArray(jType)%mol(jMol)%indx              
          if(included(jIndx) .eqv. .false.) then
            cycle
          endif              
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType, jAtom)
            rmin_ij = r_min_tab(atmType2, atmType1)
            rx = MolArray(nType)%mol(nMol)%x(nAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
            ry = MolArray(nType)%mol(nMol)%y(nAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
            rz = MolArray(nType)%mol(nMol)%z(nAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
            r = rx*rx + ry*ry + rz*rz
            if(r .lt. rmin_ij) then
              E_Trial = huge(dp)
              return
            endif
            ep = ep_tab(atmType2, atmType1)
            q = q_tab(atmType2, atmType1)

            if(ep .ne. 0E0) then
              if(r .lt. lj_cut_sq) then
                sig_sq = sig_tab(atmType2, atmType1)
                LJ = (sig_sq / r)
                LJ = LJ * LJ * LJ              
                LJ = ep * LJ * (LJ - 1E0)                
                E_LJ = E_LJ + LJ
              endif
            endif
            if(q .ne. 0E0) then
!              Ele = 0E0
              r = sqrt(r)
              Ele = q / r
              E_Ele = E_Ele + Ele
            endif
          enddo
        enddo
      enddo
     
      E_Trial = E_LJ + E_Ele
      
      end subroutine 
!======================================================================================================
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which each atom is regrown sequentially for the given trial.
      pure subroutine Rosen_Atom_Intra_New_LJ_Q(nType, nAtom, trialPos, regrown,  E_Trial, overlap)
      use Coords, only: newMol, SimpleAtomCoords
      use ForceField, only: atomArray, nonBondArray, nIntraNonBond, r_min_tab
      use ForceFieldPara_LJ_Q, only: ep_tab, sig_tab, q_tab
      implicit none
      
      logical, intent(in) :: regrown(:)
      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      
      logical, intent(inout) :: overlap
      real(dp), intent(out) :: E_Trial
      
      integer :: temp, iPair, iAtom, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ,rmin_ij

      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      overlap = .false.

      do iPair = 1,nIntraNonBond(nType)
        if(all(nonBondArray(nType,iPair)%nonMembr .ne. nAtom)) then
          cycle
        endif
        iAtom = nonBondArray(nType,iPair)%nonMembr(1)
        jAtom = nonBondArray(nType,iPair)%nonMembr(2)
         
        if(iAtom .ne. nAtom) then
          temp = iAtom
          iAtom = jAtom
          jAtom = temp
        endif

        if(regrown(jAtom) .eqv. .false.) then
          cycle
        endif

        atmType1 = atomArray(nType,iAtom)
        atmType2 = atomArray(nType,jAtom)
        ep = ep_tab(atmType1,atmType2)
        sig_sq = sig_tab(atmType1,atmType2)
        q = q_tab(atmType1,atmType2)
        rmin_ij = r_min_tab(atmType1,atmType2)
        
        rx = trialPos%x - newMol%x(jAtom)
        ry = trialPos%y - newMol%y(jAtom)
        rz = trialPos%z - newMol%z(jAtom) 
        r = rx*rx + ry*ry + rz*rz
        if(r .lt. rmin_ij) then
          E_Trial = huge(dp)
          overlap = .true.
          return
        endif
        if(ep .ne. 0E0) then
          if(r .lt. lj_cut_sq) then
            LJ = (sig_sq/r)
            LJ = LJ * LJ * LJ
            LJ = ep * LJ * (LJ-1E0)
            E_LJ = E_LJ + LJ 
          endif
        endif
        if(q .ne. 0E0) then            
          r = sqrt(r)
          Ele = q/r
          E_Ele = E_Ele + Ele            
        endif
      enddo
     
      E_Trial = E_LJ + E_Ele
      
      
      end subroutine 
!======================================================================================================
      pure subroutine Rosen_Atom_Intra_Old_LJ_Q(nType, nMol, nAtom, trialPos, regrown, E_Trial)
      use Coords, only: MolArray, SimpleAtomCoords
      use ForceField, only: atomArray, nIntraNonBond, nonBondArray, r_min_tab
      use ForceFieldPara_LJ_Q, only: ep_tab, sig_tab, q_tab
!      use SimParameters
      implicit none
      
      logical, intent(in) :: regrown(:)
      integer, intent(in) :: nType, nAtom, nMol
      type(SimpleAtomCoords), intent(in) :: trialPos
      real(dp), intent(out) :: E_Trial
      
      integer :: temp, iPair, iAtom, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ,rmin_ij

      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      
      do iPair = 1,nIntraNonBond(nType)
        if(all(nonBondArray(nType,iPair)%nonMembr .ne. nAtom)) then
          cycle
        endif
        iAtom = nonBondArray(nType, iPair)%nonMembr(1)
        jAtom = nonBondArray(nType, iPair)%nonMembr(2)
         
        if(iAtom .ne. nAtom) then
          temp = iAtom
          iAtom = jAtom
          jAtom = temp
        endif

        if(regrown(jAtom) .eqv. .false.) then
          cycle
        endif

        atmType1 = atomArray(nType, iAtom)
        atmType2 = atomArray(nType, jAtom)
        ep = ep_tab(atmType1,atmType2)
        q = q_tab(atmType1,atmType2)
        rmin_ij = r_min_tab(atmType1,atmType2)
        
        rx = trialPos%x - MolArray(nType)%mol(nMol)%x(jAtom)
        ry = trialPos%y - MolArray(nType)%mol(nMol)%y(jAtom)
        rz = trialPos%z - MolArray(nType)%mol(nMol)%z(jAtom)
        r = rx*rx + ry*ry + rz*rz
        if(r .lt. rmin_ij) then
          E_Trial = huge(dp)
          return
        endif

        if(ep .ne. 0E0) then
          if(r .lt. lj_cut_sq) then
            sig_sq = sig_tab(atmType1,atmType2)
            LJ = (sig_sq/r)
            LJ = LJ * LJ * LJ
            LJ = ep * LJ * (LJ-1E0)
            E_LJ = E_LJ + LJ
          endif
        endif
        if(q .ne. 0E0) then            
          r = sqrt(r)
          Ele = q/r
          E_Ele = E_Ele + Ele            
        endif
      enddo
      E_Trial = E_LJ + E_Ele
  
      end subroutine 
!=====================================================================================  
      pure subroutine Rosen_Atom_Intra_Gas_Old_LJ_Q(nType, nAtom, trialPos, regrown, E_Trial)
      use Coords, only: gasConfig, SimpleAtomCoords
      use ForceField, only: atomArray, nIntraNonBond, nonBondArray, r_min_tab
      use ForceFieldPara_LJ_Q, only: ep_tab, sig_tab, q_tab
!      use SimParameters
      implicit none
      
      logical, intent(in) :: regrown(:)
      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      real(dp), intent(out) :: E_Trial
      
      integer :: temp, iPair, iAtom, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ,rmin_ij

      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      
      do iPair = 1,nIntraNonBond(nType)
        if(all(nonBondArray(nType,iPair)%nonMembr .ne. nAtom)) then
          cycle
        endif
        iAtom = nonBondArray(nType, iPair)%nonMembr(1)
        jAtom = nonBondArray(nType, iPair)%nonMembr(2)
         
        if(iAtom .ne. nAtom) then
          temp = iAtom
          iAtom = jAtom
          jAtom = temp
        endif

        if(regrown(jAtom) .eqv. .false.) then
          cycle
        endif

        atmType1 = atomArray(nType, iAtom)
        atmType2 = atomArray(nType, jAtom)
        ep = ep_tab(atmType1,atmType2)
        q = q_tab(atmType1,atmType2)
        rmin_ij = r_min_tab(atmType1,atmType2)
        
        rx = trialPos%x - gasConfig(nType)%x(jAtom)
        ry = trialPos%y - gasConfig(nType)%y(jAtom)
        rz = trialPos%z - gasConfig(nType)%z(jAtom)
        r = rx*rx + ry*ry + rz*rz
        if(r .lt. rmin_ij) then
          E_Trial = huge(dp)
          return
        endif

        if(ep .ne. 0E0) then
          if(r .lt. lj_cut_sq) then
            sig_sq = sig_tab(atmType1,atmType2)
            LJ = (sig_sq/r)
            LJ = LJ * LJ * LJ
            LJ = ep * LJ * (LJ-1E0)
            E_LJ = E_LJ + LJ
          endif
        endif
        if(q .ne. 0E0) then            
          r = sqrt(r)
          Ele = q/r
          E_Ele = E_Ele + Ele            
        endif
      enddo
      E_Trial = E_LJ + E_Ele
  
      end subroutine 
!=====================================================================================   
      pure subroutine RosenLL_Atom_New_LJ_Q(nType, nAtom, trialPos, nInterAtoms, InterAtoms, E_Trial, overlap)
      use Coords, only: MolArray, SimpleAtomCoords
      use ForceField, only: r_min_tab, atomArray
      use ForceFieldPara_LJ_Q, only: q_tab, sig_tab, ep_tab
      implicit none

      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      integer, intent(in) :: nInterAtoms
      integer, intent(in) :: InterAtoms(:,:)
      logical, intent(inout) :: overlap
      real(dp), intent(out) :: E_Trial

      integer :: iInterAtoms, jType, jMol, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij

      
      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      overlap = .false.

      atmType1 = atomArray(nType, nAtom)
      do iInterAtoms = 1, nInterAtoms
        jType = InterAtoms(iInterAtoms,1)
        jMol = InterAtoms(iInterAtoms,2) 
        jAtom = InterAtoms(iInterAtoms,3)
        atmType2 = atomArray(jType, jAtom)
        rmin_ij = r_min_tab(atmType2, atmType1)
        rx = trialPos%x - MolArray(jType)%mol(jMol)%x(jAtom)
        ry = trialPos%y - MolArray(jType)%mol(jMol)%y(jAtom)
        rz = trialPos%z - MolArray(jType)%mol(jMol)%z(jAtom)
        r = rx*rx + ry*ry + rz*rz
        if(r .lt. rmin_ij) then
          E_Trial = huge(dp)
          overlap = .true.
          return
        endif
        ep = ep_tab(atmType2, atmType1)
        q = q_tab(atmType2, atmType1)
        sig_sq = sig_tab(atmType2, atmType1)   
        if(ep .ne. 0E0) then
          if(r .lt. lj_cut_sq) then
            LJ = (sig_sq / r)
            LJ = LJ * LJ * LJ              
            LJ = ep * LJ * (LJ - 1E0)                
            E_LJ = E_LJ + LJ
          endif
        endif
        if(q .ne. 0E0) then
!              Ele = 0E0
          r = sqrt(r)
          Ele = q / r
          E_Ele = E_Ele + Ele
        endif
      enddo

      E_Trial = E_LJ + E_Ele

      end subroutine 
!=====================================================================================   
      pure subroutine RosenLL_Atom_Old_LJ_Q(nType, nMol, nAtom, nInterAtoms, InterAtoms, E_Trial)
      use Coords, only: MolArray
      use ForceField, only: r_min_tab, atomArray
      use ForceFieldPara_LJ_Q, only: q_tab, ep_tab, sig_tab
      implicit none

      integer, intent(in) :: nType, nMol, nAtom
      integer, intent(in) :: nInterAtoms
      integer, intent(in) :: InterAtoms(:,:)
      real(dp), intent(out) :: E_Trial

      integer :: iInterAtoms, jType, jMol, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij

      
      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0

      atmType1 = atomArray(nType, nAtom)
      do iInterAtoms = 1, nInterAtoms
        jType = InterAtoms(iInterAtoms,1)
        jMol = InterAtoms(iInterAtoms,2) 
        jAtom = InterAtoms(iInterAtoms,3)
        atmType2 = atomArray(jType, jAtom)
        rmin_ij = r_min_tab(atmType2, atmType1)
        rx = MolArray(nType)%mol(nMol)%x(nAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
        ry = MolArray(nType)%mol(nMol)%y(nAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
        rz = MolArray(nType)%mol(nMol)%z(nAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
        r = rx*rx + ry*ry + rz*rz
        if(r .lt. rmin_ij) then
          E_Trial = huge(dp)
          return
        endif
        ep = ep_tab(atmType2, atmType1)
        q = q_tab(atmType2, atmType1)
        if(ep .ne. 0E0) then
          if(r .lt. lj_cut_sq) then
            sig_sq = sig_tab(atmType2, atmType1)
            LJ = (sig_sq / r)
            LJ = LJ * LJ * LJ              
            LJ = ep * LJ * (LJ - 1E0)                
            E_LJ = E_LJ + LJ
          endif
        endif
        if(q .ne. 0E0) then
!              Ele = 0E0
          r = sqrt(r)
          Ele = q / r
          E_Ele = E_Ele + Ele
        endif
      enddo

      E_Trial = E_LJ + E_Ele

      end subroutine 
!=====================================================================================   
      pure subroutine RosenLL_Atom_Intra_New_LJ_Q(nType, nAtom, trialPos, nIntraUnits, IntraUnits, E_Trial, overlap)
      use Coords, only: newMol, SimpleAtomCoords
      use ForceField, only: atomArray, r_min_tab
      use ForceFieldPara_LJ_Q, only: ep_tab, sig_tab, q_tab
      implicit none

      integer, intent(in) :: IntraUnits(:)
      integer, intent(in) :: nType, nAtom, nIntraUnits
      type(SimpleAtomCoords), intent(in) :: trialPos
      logical, intent(inout) :: overlap
      real(dp), intent(out) :: E_Trial

      integer :: iPair, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ,rmin_ij

      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      overlap = .false.

      atmType1 = atomArray(nType,nAtom)
      do iPair = 1,nIntraUnits
        jAtom = IntraUnits(iPair)
        atmType2 = atomArray(nType,jAtom)
        ep = ep_tab(atmType1,atmType2)
        sig_sq = sig_tab(atmType1,atmType2)
        q = q_tab(atmType1,atmType2)
        rmin_ij = r_min_tab(atmType1,atmType2)
        
        rx = trialPos%x - newMol%x(jAtom)
        ry = trialPos%y - newMol%y(jAtom)
        rz = trialPos%z - newMol%z(jAtom) 
        r = rx*rx + ry*ry + rz*rz
        if(r .lt. rmin_ij) then
          E_Trial = huge(dp)
          overlap = .true.
          return
        endif
        if(ep .ne. 0E0) then
          if(r .lt. lj_cut_sq) then
            LJ = (sig_sq/r)
            LJ = LJ * LJ * LJ
            LJ = ep * LJ * (LJ-1E0)
            E_LJ = E_LJ + LJ 
          endif
        endif
        if(q .ne. 0E0) then            
          r = sqrt(r)
          Ele = q/r
          E_Ele = E_Ele + Ele            
        endif
      enddo
      
      E_Trial = E_LJ + E_Ele

      end subroutine 
!=====================================================================================  
      pure subroutine RosenLL_Atom_Intra_Old_LJ_Q(nType, nMol, nAtom, trialPos, nIntraUnits, IntraUnits, E_Trial)
      use Coords, only: MolArray, SimpleAtomCoords
      use ForceField, only: atomArray, r_min_tab
      use ForceFieldPara_LJ_Q, only: ep_tab, sig_tab, q_tab
!      use SimParameters
      implicit none

      integer, intent(in) :: IntraUnits(:)
      integer, intent(in) :: nType, nMol, nAtom, nIntraUnits
      type(SimpleAtomCoords), intent(in) :: trialPos      
      real(dp), intent(out) :: E_Trial

      integer :: iPair, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ,rmin_ij

      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      atmType1 = atomArray(nType,nAtom)
      do iPair = 1,nIntraUnits
        jAtom = IntraUnits(iPair)
        atmType2 = atomArray(nType, jAtom)
        ep = ep_tab(atmType1,atmType2)
        q = q_tab(atmType1,atmType2)
        rmin_ij = r_min_tab(atmType1,atmType2)
        
        rx = trialPos%x - MolArray(nType)%mol(nMol)%x(jAtom)
        ry = trialPos%y - MolArray(nType)%mol(nMol)%y(jAtom)
        rz = trialPos%z - MolArray(nType)%mol(nMol)%z(jAtom)
        r = rx*rx + ry*ry + rz*rz
        if(r .lt. rmin_ij) then
          E_Trial = huge(dp)
          return
        endif

        if(ep .ne. 0E0) then
          if(r .lt. lj_cut_sq) then
            sig_sq = sig_tab(atmType1,atmType2)
            LJ = (sig_sq/r)
            LJ = LJ * LJ * LJ
            LJ = ep * LJ * (LJ-1E0)
            E_LJ = E_LJ + LJ
          endif
        endif
        if(q .ne. 0E0) then            
          r = sqrt(r)
          Ele = q/r
          E_Ele = E_Ele + Ele            
        endif
      enddo

      E_Trial = E_LJ + E_Ele

      end subroutine 
!=====================================================================================   
      pure subroutine RosenLL_Atom_Intra_Gas_Old_LJ_Q(nType, nAtom, trialPos, nIntraUnits, IntraUnits,  E_Trial)
      use Coords, only: gasConfig, SimpleAtomCoords
      use ForceField, only: atomArray, nIntraNonBond, r_min_tab
      use ForceFieldPara_LJ_Q, only: ep_tab, sig_tab, q_tab
!      use SimParameters
      implicit none

      integer, intent(in) :: IntraUnits(:)
      integer, intent(in) :: nType, nAtom, nIntraUnits
      type(SimpleAtomCoords), intent(in) :: trialPos      
      real(dp), intent(out) :: E_Trial

      integer :: iPair, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ,rmin_ij

      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      atmType1 = atomArray(nType,nAtom)
      do iPair = 1,nIntraUnits
        jAtom = IntraUnits(iPair)
        atmType2 = atomArray(nType, jAtom)
        ep = ep_tab(atmType1,atmType2)
        q = q_tab(atmType1,atmType2)
        rmin_ij = r_min_tab(atmType1,atmType2)
        
        rx = trialPos%x - gasConfig(nType)%x(jAtom)
        ry = trialPos%y - gasConfig(nType)%y(jAtom)
        rz = trialPos%z - gasConfig(nType)%z(jAtom)
        r = rx*rx + ry*ry + rz*rz
        if(r .lt. rmin_ij) then
          E_Trial = huge(dp)
          return
        endif

        if(ep .ne. 0E0) then
          if(r .lt. lj_cut_sq) then
            sig_sq = sig_tab(atmType1,atmType2)
            LJ = (sig_sq/r)
            LJ = LJ * LJ * LJ
            LJ = ep * LJ * (LJ-1E0)
            E_LJ = E_LJ + LJ
          endif
        endif
        if(q .ne. 0E0) then            
          r = sqrt(r)
          Ele = q/r
          E_Ele = E_Ele + Ele            
        endif
      enddo


      E_Trial = E_LJ + E_Ele

      end subroutine 
!=====================================================================================   
!=====================================================================================   
!=====================================================================================   
      pure subroutine Find_InterAtms_LJ_Q(RefAtm, r0, atmType1, nIndx, nInterAtoms, InterAtoms)
      use Coords, only: MolArray, SimpleAtomCoords
      use ForceField, only: nAtoms
      use SimParameters, only: nMolTypes, NPART
      implicit none 

      type(SimpleAtomCoords), intent(in) :: RefAtm
      real(dp), intent(in) :: r0
      integer(kind=atomIntType), intent(in) :: atmType1
      integer, intent(in) :: nIndx
      integer, intent(inout) :: nInterAtoms
      integer, intent(out) :: InterAtoms(:,:)

      integer :: jType, jMol, jAtom, jIndx

      nInterAtoms = 0
      do jType = 1, nMolTypes
        do jMol = 1,NPART(jType) 
          jIndx = molArray(jType)%mol(jMol)%indx     
          if (jIndx .eq. nIndx) cycle
          do jAtom = 1,nAtoms(jType)    
            nInterAtoms = nInterAtoms + 1
            InterAtoms(nInterAtoms,1) = jType
            InterAtoms(nInterAtoms,2) = jMol
            InterAtoms(nInterAtoms,3) = jAtom
          enddo
        enddo
      enddo

      end subroutine 
!======================================================================================================
      end module 


