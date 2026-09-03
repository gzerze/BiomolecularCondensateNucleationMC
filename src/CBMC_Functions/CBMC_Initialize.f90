!=========================================================
!     The purpose of this function is to generate the data structures
!     that will be used in the Configuration Bias Monte Carlo regrowth
!     functions.
!     The terminology used in this segment is as followed
!           Terminal: An atom that has one bond associated with it. These atoms are 
!                     located at the end of a molecular chain.
!             Linker: An atom that has two bonds associated with it.  These atoms are
!                     segments found in the middle of a long chain
!                Hub: An atom that has three or more bonds assocaited with it. These atoms
!                     are where the molecule splits into multiple branches.
      subroutine CBMC_CreateTopology
      use Forcefield
      use SimParameters
      use Coords
      use CBMC_Variables
      use ParallelVar
      implicit none
      integer :: i, iBin, iType, iAtom, iBond
      integer :: AllocateStatus
      integer :: cnt,curAtom
      integer :: nTerminal, nLinker, nHub      
      

      allocate(regrowType(1:nMolTypes), STAT = AllocateStatus)
      select case(trim(adjustl(ForceFieldName)))
      case("Pedone") 
        regrowType = 0
        return
      case("Tersoff") 
        regrowType = 0
        return
      case default
        regrowType = -1
      end select


      do iType = 1, nMolTypes
        if(nAtoms(iType) .eq. 1) then
          regrowType(iType) = 0
        endif
      enddo
      if(all(regrowType .eq. 0) )then
  ! Still allocate minimal pathArray so any accidental CBMC call won't crash
        allocate(pathArray(1:nMolTypes), STAT=AllocateStatus)
        do iType = 1, nMolTypes
          pathArray(iType)%nPaths = 1
          allocate(pathArray(iType)%path(1:1,1:nAtoms(iType)), STAT=AllocateStatus)
          allocate(pathArray(iType)%pathMax(1:1), STAT=AllocateStatus)
          pathArray(iType)%path = 1
          pathArray(iType)%pathMax(1) = 1
        enddo
        return
      endif


!      Allocate the Topology arrays      
      allocate(topolArray(1:nMolTypes),STAT = AllocateStatus)
      do iType = 1,nMolTypes
        allocate(topolArray(iType)%atom(1:nAtoms(iType)),STAT = AllocateStatus)                  
      enddo
!      allocate(regrowType(1:nMolTypes), STAT = AllocateStatus)
      allocate(regrowOrder(1:nMolTypes, 1:maxAtoms), STAT = AllocateStatus)
      allocate(pathArray(1:nMolTypes), STAT = AllocateStatus)
      allocate(usedByPath(1:nMolTypes,1:maxAtoms), STAT = AllocateStatus)
      allocate(atomPathIndex(1:nMolTypes,1:maxAtoms), STAT = AllocateStatus)
      allocate(probTypeCBMC(1:nMolTypes), STAT = AllocateStatus)
      allocate(SwapGrowOrder(1:nMolTypes), STAT = AllocateStatus)
	  


	  
!     Allocate Arrays for Growing Branched Molecules
      do iType = 1, nMolTypes
!       GrowFrom(i) is the atom number where growing at step "i" occurs from.
        allocate(SwapGrowOrder(iType)%GrowFrom(1:maxAtoms), STAT = AllocateStatus)
!       GrowPrev(i) is the atom number that has already been grown and connected to GrowFrom(i) in step "i"
        allocate(SwapGrowOrder(iType)%GrowPrev(1:maxAtoms), STAT = AllocateStatus)
!       GrowNum(i) is the number of atoms that are supposed to be grown in step "i"
        allocate(SwapGrowOrder(iType)%GrowNum(1:maxAtoms), STAT = AllocateStatus)
!       GrowList(i,j) is the jth atom that is grown in ith step (1 <= j <= GrowNum(i))
        allocate(SwapGrowOrder(iType)%GrowList(1:maxAtoms,1:maxBranches), STAT = AllocateStatus)
!       TorNum(i) is the number of atoms that are connetcted to GrowPrev(i), except GrowFrom(i), that are used to 
!       calculate torsional energies for growing atoms in step "i"
        allocate(SwapGrowOrder(iType)%TorNum(1:maxAtoms), STAT = AllocateStatus)
!       TorList(i,j) is the jth atom that is connected to GrowPrev(i) in ith step (1 <= j <= TorNum(i))
        allocate(SwapGrowOrder(iType)%TorList(1:maxAtoms,1:maxBranches), STAT = AllocateStatus)
      enddo
      
!      Search through each molecule and determine the number of bonds each atom is
!      a part of. This is done in order to classify each atom as either a Terminal atom(numBonds = 1),
!      Linker Atom(numBonds = 2), or Hub Atom (numBonds >= 3).  This information is used to classify
!      which type of CBMC algorithm will be required to regrow each molecule.
      totalDihed = 0
      do iType = 1, nMolTypes
        do iAtom = 1, nAtoms(iType)
          cnt = 0
          topolArray(iType)%atom(iAtom) = 0
          do iBond = 1, nBonds(iType)
            if(bondArray(iType,iBond)%bondMembr(1) .eq. iAtom) then
              cnt = cnt + 1
              cycle
            endif
            if(bondArray(iType,iBond)%bondMembr(2) .eq. iAtom) then
              cnt = cnt + 1
              cycle
            endif            
          enddo 
          topolArray(iType)%atom(iAtom) = int(cnt, atomIntType)
        enddo
      enddo

!      Using the data that was previously calcualted classify each atom as terminal, linker, or hub and then tabulate
!      the data.
      do iType = 1, nMolTypes
        nTerminal = 0
        nHub = 0 
        nLinker = 0 
        do iAtom = 1, nAtoms(iType)         
           if(topolArray(iType)%atom(iAtom) .le. 1) then
             nTerminal = nTerminal + 1
           elseif(topolArray(iType)%atom(iAtom) .eq. 2) then
             nLinker = nLinker + 1
           elseif(topolArray(iType)%atom(iAtom) .ge. 3) then
             nHub = nHub + 1
!             Count the number of dihedral angles that are used.
             if(topolArray(iType)%atom(iAtom) .eq. 3) then
               totalDihed = totalDihed + 1
             elseif(topolArray(iType)%atom(iAtom) .eq. 4) then
               totalDihed = totalDihed + 2
             endif
           endif
        enddo
        pathArray(iType)%nTerminal = nTerminal        
        pathArray(iType)%nLinker = nLinker
        pathArray(iType)%nHub = nHub
        if(nTerminal .ne. 0) then
           allocate(pathArray(iType)%termAtoms(1:nTerminal),STAT = AllocateStatus)
        else
           write(nout,*) "Molecule is Cyclic, not yet supported"           
           stop "Molecule is Cyclic, not yet supported"
        endif
        if(nHub .ne. 0) then
           allocate(pathArray(iType)%hubAtoms(1:nHub),STAT = AllocateStatus)        
        endif        


!        This block stores the atom indicies for each Terminal and Hub atom for quick access
!        later in the program.
        nTerminal = 0
        nHub = 0
        nLinker = 0
        do iAtom = 1, nAtoms(iType)         
           if(topolArray(iType)%atom(iAtom) .le. 1) then
             nTerminal = nTerminal + 1
             pathArray(iType)%termAtoms(nTerminal) = iAtom
           elseif(topolArray(iType)%atom(iAtom) .ge. 3) then
             nHub = nHub + 1
             pathArray(iType)%hubAtoms(nHub) = iAtom
           endif
        enddo
      enddo


      allocate(dihedData(1:totalDihed), STAT = AllocateStatus)
      cnt = 0
      write(35,*) "--------------------------------------------"
      write(35,*) "Dihedral Angles"
      do iType = 1, nMolTypes
        do iAtom = 1, nAtoms(iType)
          if(topolArray(iType)%atom(iAtom) .eq. 3) then
            cnt = cnt + 1
            dihedData(cnt)%molType = iType
            dihedData(cnt)%hubIndx = iAtom
            dihedData(cnt)%dihedIndx = 1
            write(35,*) iType, iAtom, 1
          elseif(topolArray(iType)%atom(iAtom) .eq. 4) then
            cnt = cnt + 1
            dihedData(cnt)%molType = iType
            dihedData(cnt)%hubIndx = iAtom
            dihedData(cnt)%dihedIndx = 1
            write(35,*) iType, iAtom, 1
            cnt = cnt + 1
            dihedData(cnt)%molType = iType
            dihedData(cnt)%hubIndx = iAtom
            dihedData(cnt)%dihedIndx = 2
            write(35,*) iType, iAtom, 2
          endif
        enddo
      enddo
      
      do i = 1, totalDihed
        dihedData(i)%Prob = 1E0 / nDihBins
        dihedData(i)%Hist = 0E0
        dihedData(i)%Integral(0) = dihedData(i)%Prob(0)
        do iBin = 1, nDihBins
          dihedData(i)%Integral(iBin) = dihedData(i)%Integral(iBin-1) + dihedData(i)%Prob(iBin)
        enddo
        dihedData(i)%accConst = nDihBins*diBinSize
      enddo

      call CBMC_FindPathways

      usedByPath = 0
      do iType = 1, nMolTypes
        do i = 1, pathArray(iType)%nPaths
          do iAtom = 1, nAtoms(iType)
            curAtom = pathArray(iType)%path(i,iAtom)
            if(curAtom .eq. 0) then
              cycle
            endif
            if(topolArray(iType)%atom(curAtom) .le. 2) then
              usedByPath(iType,curAtom) = i
              atomPathIndex(iType,curAtom) = iAtom
            else
              usedByPath(iType,curAtom) = -1
              atomPathIndex(iType,curAtom) = -1
            endif
          enddo
        enddo
      enddo
      


      call DetermineRegrowthType
      call DetermineRegrowOrder

      
      
      end subroutine
!=========================================================
!      The purpose of this subroutine is to divide each molecule into chain segments.  A segment is the sequence 
!      of atoms that start with either a hub or a terminal atom and end with another terminal or hub atom.  
!      For instance in the molecule 2-methyl-butane [CH3-CH2-CH(CH3)-CH3] one segment would
!      be defined as CH3-CH2-CH since CH3 is a terminal atom (numBonds = 1) and CH is a Hub atom
!      (numBonds >= 3).  This information is used by the CBMC algorithm to map out which algorithms
!      need to be used in order to efficiently regrow each atom.  For instance a Hub atom is much more
!      difficult to regrow than a Linker (numBonds = 2) since there are multiple bending angles associated
!      with this regrowth.
      subroutine CBMC_FindPathways
      use Forcefield
      use SimParameters
      use Coords
      use CBMC_Variables
      implicit none
      integer :: cnt
      integer :: iType,iAtom, iPath
      integer :: nPaths,curBonds      
      integer, allocatable :: tempArray(:,:)
      integer :: atmPrev, atmCur, atmNext
      integer :: AllocateStatus      
	  
      integer :: nTerminal,nBranch,iBranch,nRemainPaths
      integer :: atmFirst,atmLast
      integer, allocatable :: OnePath(:),tempPathMax(:),BondedAtoms(:)


      allocate(tempArray(1:50,1:maxAtoms))
      allocate(OnePath(1:maxAtoms))
      allocate(tempPathMax(1:50))
      allocate(BondedAtoms(1:maxBranches))


      do iType = 1, nMolTypes
         tempArray = 0
!     For a monoatomic molecule, there is only one path with one member
         if (nAtoms(iType) .eq. 1) then
            nPaths = 1
            tempArray(1,1:1) = 1
            goto 56
         endif
!     Starting with the first Terminal
         nTerminal = 1
!     There is no previous atom
         atmPrev = 0
         atmCur = pathArray(iType)%termAtoms(nTerminal)
!     Find the atom connected to atmCur
         call getNextAtom(iType,atmPrev,atmCur,atmNext)
         nPaths = 1
         nRemainPaths = 1
         iPath = 1
         tempArray(nPaths,1) = atmCur
         tempArray(nPaths,2) = atmNext
         atmPrev = atmCur
         atmCur = atmNext
         cnt = 2
!     Start forming all paths
         do while (nRemainPaths .gt. 0)
            curBonds = topolArray(iType)%atom(atmCur)
!     Find all the linkers on the path until a terminal or hud is found (the end of the path)
            do while (curBonds .eq. 2)
               call getNextAtom(iType,atmPrev,atmCur,atmNext)
               atmPrev = atmCur
               atmCur = atmNext
               cnt = cnt + 1
               tempArray(iPath,cnt) = atmCur
               curBonds = topolArray(iType)%atom(atmCur)
            enddo
!     Number of atoms in iPath
            tempPathMax(iPath) = cnt
            nRemainPaths = nRemainPaths - 1 + curBonds - 1
            if (curBonds .gt. 2) then
               nBranch = curBonds - 1
!     Find all branches that are connected to "atmCur" except "atmPrev"
               call getBondedAtoms(iType,atmPrev,atmCur,nBranch,BondedAtoms)
               do iBranch = 1,nBranch
                  nPaths = nPaths + 1
!     Allocate the first two atms in the path
                  tempArray(nPaths,1) = atmCur
                  tempArray(nPaths,2) = BondedAtoms(iBranch)
               enddo
            endif
            iPath = iPath + 1
            atmPrev = tempArray(iPath,1)
            atmCur = tempArray(iPath,2)
            cnt = 2
         enddo
         iPath = iPath - 1
         if (iPath .ne. nPaths) then
            write(*,*) "Failed to trace paths", iPath, nPaths
            stop
         endif
56       pathArray(iType)%nPaths = nPaths
         allocate(pathArray(iType)%path(1:nPaths,1:nAtoms(iType)), STAT = AllocateStatus)
         allocate(pathArray(iType)%pathMax(1:nPaths), STAT = AllocateStatus)
         do iPath = 1,nPaths
            if (nAtoms(iType) .eq. 1) then
              pathArray(iType)%pathMax(iPath) = 1
            else
              pathArray(iType)%pathMax(iPath) = tempPathMax(iPath)
            endif
            if (nAtoms(iType) .eq. 1) then
              atmFirst = 1
              atmLast  = 1
            else
              atmFirst = tempArray(iPath,1)
              atmLast  = tempArray(iPath,tempPathMax(iPath))
            endif
!            atmFirst = tempArray(iPath,1)
!            atmLast = tempArray(iPath,tempPathMax(iPath))
!     For those paths that has one Terminal and one hub, the Terminal must be the first atom and the hub must be the last atom
            if (topolArray(iType)%atom(atmFirst) .gt. 2) then
               if (topolArray(iType)%atom(atmLast) .eq. 1) then
                  OnePath = 0
                  do iAtom = 1,tempPathMax(iPath)
                     OnePath(tempPathMax(iPath) - iAtom + 1) = tempArray(iPath,iAtom)
                  enddo
                  do iAtom = 1,tempPathMax(iPath)
                     tempArray(iPath,iAtom) = OnePath(iAtom)
                  enddo
               endif
            endif
         enddo
         do iPath = 1, nPaths
           pathArray(iType)%path(iPath,1:nAtoms(iType)) = tempArray(iPath,1:nAtoms(iType))
         enddo
      enddo

!      write(35,*) "------------------------------------------------------------------"
!      do iType = 1, nMolTypes
!        write(35,*) "iType:", iType
!        write(35,*) "Number of Paths:", pathArray(iType)%nPaths
!        do iPath = 1, pathArray(iType)%nPaths
!          write(35,*) pathArray(iType)%path(iPath,1:pathArray(iType)%pathMax(iPath)) 
!        enddo
!        write(35,*)
!      enddo
!      write(35,*) "------------------------------------------------------------------"

      end subroutine
!=========================================================
      subroutine getNextAtom(iType,atmPrev,atmCur,atmNext)
      use Forcefield
      use SimParameters
      use Coords
      implicit none
      
      integer, intent(in) :: iType,atmPrev,atmCur
      integer, intent(out) :: atmNext
      
      logical :: atmFound
      integer :: iBond
      
      

      do iBond = 1, nBonds(iType)
         atmFound = .false.      
         if(any(bondArray(iType,iBond)%bondMembr .eq. atmCur)) then
           atmFound = .true.     
         endif
         if(atmFound) then
           if(all(bondArray(iType,iBond)%bondMembr .ne. atmPrev) ) then
              if(bondArray(iType,iBond)%bondMembr(1) .eq. atmCur) then
                atmNext = bondArray(iType,iBond)%bondMembr(2)
!                bondUsed = iBond
                return
              else
                atmNext = bondArray(iType,iBond)%bondMembr(1)
!                bondUsed = iBond
                return              
              endif
           endif
         endif
      enddo
      
      write(*,*) "Error in getNextAtom subroutine!!" 
      write(*,*) "Dead End! Unable to find next atom in the chain"
      write(*,*) "Mol Type:", iType
      write(*,*) "Previous Atom:", atmPrev
      write(*,*) "Current Atom:", atmCur
      stop
      
      end subroutine
!===============================================================================
      subroutine getBondedAtoms(iType,atmPrev,atmCur,nBranch,BondedAtoms)
      use Forcefield
      use SimParameters
      use Coords
      use CBMC_Variables
      implicit none
	  
      integer, intent(in) :: iType,atmPrev,atmCur,nBranch
      
      logical :: atmFound
      integer :: iBond,BondedAtoms(1:maxBranches),iBranch
	  
	  
      iBranch = 0
      BondedAtoms = 0
	  
      do iBond = 1, nBonds(iType)
         atmFound = .false.      
         if(any(bondArray(iType,iBond)%bondMembr .eq. atmCur)) then
           atmFound = .true.     
         endif
         if(atmFound) then
           if(all(bondArray(iType,iBond)%bondMembr .ne. atmPrev) ) then
              iBranch = iBranch + 1
              if(bondArray(iType,iBond)%bondMembr(1) .eq. atmCur) then
                BondedAtoms(iBranch) = bondArray(iType,iBond)%bondMembr(2)
                if (iBranch .eq. nBranch) return
              else
                BondedAtoms(iBranch) = bondArray(iType,iBond)%bondMembr(1)
                if (iBranch .eq. nBranch) return              
              endif
           endif
         endif
      enddo
	  
	  
	  
      write(*,*) "Error in getBondedAtoms subroutine!!" 
      write(*,*) "Dead End! Unable to find all branches"
      write(*,*) "Mol Type:", iType
      write(*,*) "Previous Atom:", atmPrev
      write(*,*) "Current Atom:", atmCur
      stop
      
      end subroutine
!===============================================================================
!     This subroutine determines the type of regrowth function to be used in
!     a CBMC or AVBMC regrowth.  Each molecule is assigned an integer based on
!     the molecule's flexbility and size.  The integers correspond to these molecule types
!           0: Completely Ridgid Molecule
!           1: Simple Flexible Molecule
!           2: Straight Chain Molecule
!           3: Branched Chain Molecule
!
      subroutine DetermineRegrowthType
      use Forcefield
      use SimParameters
      use Coords
      use CBMC_Variables
      implicit none
      
      logical :: bondRidgid, bendRidgid
      integer :: iType, i
      integer :: xType
      real(dp) :: xForceConst, norm
      
      regrowType = -1
      
      
      do iType = 1, nMolTypes
         if(nAtoms(iType) .eq. 1) then
           regrowType(iType) = 0
           cycle
         endif
         bondRidgid = .true.
         do i = 1, nBonds(iType)
           xType = bondArray(iType,i)%bondType
           xForceConst = bondData(xType)%k_eq
           if(xForceConst .ne. 0d0) then
             bondRidgid = .false.
             exit
           endif
         enddo
         bendRidgid = .true.
         do i = 1, nAngles(iType)
           xType = bendArray(iType,i)%bendType
           xForceConst = bendData(xType)%k_eq
           if(xForceConst .ne. 0d0) then
             bendRidgid = .false.
             exit
           endif
         enddo
         if((bondRidgid .eqv. .true.) .and. (bendRidgid .eqv. .true.)) then
           if(nTorsional(iType) .eq. 0) then
             regrowType(iType) = 0
           endif
         else
           if(nTorsional(iType) .eq. 0) then
             regrowType(iType) = 1
           else
             if(pathArray(iType)%nHub .eq. 0) then
               regrowType(iType) = 2
             else
               regrowType(iType) = 3             
             endif
           endif
         endif
         write(35,*) "Regrow Type:",regrowType(iType)
      enddo

      probTypeCBMC = 0d0
      do iType = 1, nMolTypes
        if(regrowType(iType) .ge. 1) then
          probTypeCBMC(iType) = 1d0
        endif
      enddo
      norm = sum(probTypeCBMC)
      if( any(probTypeCBMC .ne. 0d0) ) then
        do iType = 1, nMolTypes
          probTypeCBMC(iType) = probTypeCBMC(iType)/norm
        enddo
      endif

      
      end subroutine  
!===============================================================================
!   
      subroutine DetermineRegrowOrder
      use Forcefield
      use SimParameters
      use Coords
      use CBMC_Variables
      implicit none
      
      integer :: iType, iAtom, Atom1Loc
      integer :: cnt
      integer :: iPath
      
      regrowOrder = 0
      
      do iType = 1, nMolTypes
        select case(regrowType(iType))
        case(0)
          continue

        case(1)
          continue

        case(2)
          regrowOrder(iType,1) = 1
          do iAtom = 1, pathArray(iType)%pathMax(1)
            if(pathArray(iType)%path(1,iAtom) .eq. 1)then
              Atom1Loc = iAtom
              exit
            endif
          enddo
          write(35,*) pathArray(iType)%path(1,:)

          cnt = 1
          if(dble(Atom1Loc) .ge. dble(pathArray(iType)%pathMax(1)+1)/2d0) then
            do iAtom = Atom1Loc+1,pathArray(iType)%pathMax(1)
              cnt = cnt + 1
              regrowOrder(iType,cnt) = pathArray(iType)%path(1,iAtom) 
            enddo
            do iAtom = 1, Atom1Loc-1
              cnt = cnt + 1
              regrowOrder(iType,cnt) = pathArray(iType)%path(1,iAtom) 
            enddo
          else
            do iAtom = 1, Atom1Loc-1
              cnt = cnt + 1
              regrowOrder(iType,cnt) = pathArray(iType)%path(1,iAtom) 
            enddo
            do iAtom = Atom1Loc+1,pathArray(iType)%pathMax(1)
              cnt = cnt + 1
              regrowOrder(iType,cnt) = pathArray(iType)%path(1,iAtom) 
            enddo 
          endif
        case (3)
          SwapGrowOrder(iType)%GrowFrom = 0
          SwapGrowOrder(iType)%GrowPrev = 0
          SwapGrowOrder(iType)%GrowNum = 0
          SwapGrowOrder(iType)%GrowList = 0
          SwapGrowOrder(iType)%TorNum = 0
          SwapGrowOrder(iType)%TorList = 0
          if (topolArray(iType)%atom(1) .eq. 1) then
             do iPath = 1,pathArray(iType)%nPaths
                if (pathArray(iType)%path(iPath,1) .eq. 1) then
                   if (pathArray(iType)%path(iPath,2) .ne. 2) then
                      write(*,*) "Error! The first and the second atoms in type", iType, "must be connected"
                      stop
                   endif
                   exit
                endif
             enddo
             SwapGrowOrder(iType)%GrowPrev(1) = 1
             SwapGrowOrder(iType)%GrowFrom(1) = 2
          elseif (topolArray(iType)%atom(2) .eq. 1) then
             do iPath = 1,pathArray(iType)%nPaths
                if (pathArray(iType)%path(iPath,1) .eq. 2) then
                   if (pathArray(iType)%path(iPath,2) .ne. 1) then
                      write(*,*) "Error! The first and the second atoms in type", iType, "must be connected"
                      stop
                   endif
                   exit
                endif
             enddo
             SwapGrowOrder(iType)%GrowPrev(1) = 2
             SwapGrowOrder(iType)%GrowFrom(1) = 1
          else
            write(*,*) iType, regrowType(iType)
            write(*,*) "Error in DetermineRegrowOrder function!"
            write(*,*) "Invalid Regrow Type!"
            write(*,*) "The first or the second atom must be a terminal!"
            stop 
          endif
		  
          call Schedule_BranchedMol_Swap(iType)
		  
        case default
          write(*,*) iType, regrowType(iType)
          write(*,*) "Error in DetermineRegrowOrder function!"
          write(*,*) "Invalid Regrow Type!"
          stop 


        end select

   
      enddo

      
      end subroutine  
!=========================================================      
      
      subroutine Schedule_BranchedMol_Swap(nType)
      use SimParameters
      use CBMC_Variables
      use CBMC_Utility, only: getRandomBondedAtoms
      implicit none
	  
      integer, intent(in) :: nType
      logical :: LGrowing
      integer :: atmCur, atmPrev, atmNext, curBonds, preBonds, TotalGrowing, OuterNum, OuterAtoms(1:maxAtoms)
      integer :: OuterPrev(1:maxAtoms), iBond, iu, iufrom, iCBunit, BondedAtoms(1:maxAtoms), OuterTry, nGrow, iGrow
      real(dp) :: grnd
	  
      nGrow = 1
      atmCur = SwapGrowOrder(nType)%GrowFrom(nGrow)
      atmPrev = SwapGrowOrder(nType)%GrowPrev(nGrow)
      SwapGrowOrder(nType)%TorNum(nGrow) = 0
      curBonds = topolArray(nType)%atom(atmCur)
	  
      select case(curBonds)
      case(1) ! A Terminal has been selected, an Error happened
         write(*,*) "Error! This molecule of type", nType, "is Diatomic NOT Branched"
         stop
      case(2) ! A Linker is Growing
         SwapGrowOrder(nType)%GrowNum(nGrow) = curBonds - 1
         iGrow = SwapGrowOrder(nType)%GrowNum(nGrow)
         call getNextAtom(nType,atmPrev,atmCur,atmNext)
         SwapGrowOrder(nType)%GrowList(nGrow,iGrow) = atmNext
      case default ! A Hub is Growing
         SwapGrowOrder(nType)%GrowNum(nGrow) = curBonds - 1
         iGrow = SwapGrowOrder(nType)%GrowNum(nGrow)
         LGrowing = .true.
         call getRandomBondedAtoms(nType,atmCur,atmPrev,iGrow,BondedAtoms,LGrowing)
         do iBond = 1,iGrow
            SwapGrowOrder(nType)%GrowList(nGrow,iBond) = BondedAtoms(iBond)
         enddo
      end select
      TotalGrowing = 0 
      OuterNum = 0
      iufrom = SwapGrowOrder(nType)%GrowFrom(nGrow)
      do iCBunit = 1,SwapGrowOrder(nType)%GrowNum(nGrow)
         TotalGrowing = TotalGrowing + 1
         iu = SwapGrowOrder(nType)%GrowList(nGrow,iCBunit)
         iBond = topolArray(nType)%atom(iu)
         if (iBond .gt. 1) then
            OuterNum = OuterNum + 1
            OuterAtoms(OuterNum) = iu
            OuterPrev(OuterNum) = iufrom
         endif
      enddo
      do while (OuterNum .gt. 0)
         OuterTry = floor(grnd()*OuterNum + 1d0)
         nGrow = nGrow + 1
         atmCur = OuterAtoms(OuterTry)
         atmPrev = OuterPrev(OuterTry)
         SwapGrowOrder(nType)%GrowFrom(nGrow) = atmCur
         SwapGrowOrder(nType)%GrowPrev(nGrow) = atmPrev
         SwapGrowOrder(nType)%GrowNum(nGrow) = topolArray(nType)%atom(atmCur) - 1
         iGrow = SwapGrowOrder(nType)%GrowNum(nGrow)
         LGrowing = .true.
         call getRandomBondedAtoms(nType,atmCur,atmPrev,iGrow,BondedAtoms,LGrowing)
         do iBond = 1,iGrow
            SwapGrowOrder(nType)%GrowList(nGrow,iBond) = BondedAtoms(iBond)
         enddo
         preBonds = topolArray(nType)%atom(atmPrev)
         SwapGrowOrder(nType)%TorNum(nGrow) = preBonds - 1
         iGrow = SwapGrowOrder(nType)%TorNum(nGrow)
         if (iGrow .gt. 0) then
            LGrowing = .false.
            call getRandomBondedAtoms(nType,atmPrev,atmCur,iGrow,BondedAtoms,LGrowing)
            do iBond = 1,iGrow
               SwapGrowOrder(nType)%TorList(nGrow,iBond) = BondedAtoms(iBond)
            enddo
         endif
! Update list of "Outer" beads: Remove bead that was just grown from Outer list
         OuterAtoms(OuterTry) = OuterAtoms(OuterNum)
         OuterPrev(OuterTry) = OuterPrev(OuterNum)
         OuterNum = OuterNum - 1
! Add the new beads if they have more to be grown
         iufrom = atmCur
         do iCBunit = 1,SwapGrowOrder(nType)%GrowNum(nGrow)
            TotalGrowing = TotalGrowing + 1
            iu = SwapGrowOrder(nType)%GrowList(nGrow,iCBunit)
            iBond = topolArray(nType)%atom(iu)
            if (iBond .gt. 1) then
               OuterNum = OuterNum + 1
               OuterAtoms(OuterNum) = iu
               OuterPrev(OuterNum) = iufrom
            endif
         enddo
      enddo
      SwapGrowOrder(nType)%GrowthSteps = nGrow
      
      end subroutine  
	  
!=========================================================      
      
