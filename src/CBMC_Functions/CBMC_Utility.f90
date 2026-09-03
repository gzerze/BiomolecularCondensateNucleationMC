      module CBMC_Utility
      use VarPrecision

      real(dp), parameter :: rosenCut_sq = 4d0**2
      contains
!=======================================================================
!     The purpose of this subroutine is to create a list of molecules
!     that will be included in the intermolecular component of the Rosenbluth
!     weight.  
      subroutine Rosen_CreateSubset(nTarget, included)
      use Coords
      use Constants
      use ForceField
      use SimParameters
      implicit none
      integer, intent(in) :: nTarget
      logical, intent(out) :: included(1:maxMol)
      
      logical:: inc_firstPass(1:maxMol)
      integer :: jIndx
      included = .false.
      inc_firstPass = .false.
      
!      Loop 1: Adds the neighbors of the target molecule to the list      
      do jIndx = 1, maxMol
         if(.not. isActive(jIndx)) then
           cycle
         endif
         if(NeighborList(nTarget, jIndx)) then
!            inc_firstPass(jIndx) = .true.
            included(jIndx) = .true.         
         endif
      enddo

!      nType = typeList(nTarget)
!      nMol = subIndxList(nTarget)
!      globIndx1 = molArray(nType)%mol(nMol)%globalIndx(1)
!      do jIndx = 1, maxMol
!         if(.not. isActive(jIndx)) then
!           cycle
!         endif
!         jType = typeList(jIndx)
!         jMol = subIndxList(jIndx)
!         globIndx2 = molArray(jType)%mol(jMol)%globalIndx(1)
!         if(rPair(globIndx1, globIndx2)%p%r_sq .lt. rosenCut_sq) then
!           included(jIndx) = .true.         
!         endif
!      enddo


!      Loop 2: Adds the neighbors of the molecules added in the previous step to the list      
!      do jIndx = 1, maxMol
!       if(.not. isActive(jIndx)) then
!         cycle
!       endif
!       if(inc_firstPass(jIndx) .eqv. .true.) then
!         do kIndx = 1, maxMol
!           if(.not. isActive(kIndx)) then
!             cycle
!           endif           
!           if(NeighborList(jIndx, kIndx)) then
!             included(kIndx) = .true.         
!           endif
!         enddo
!       endif
!      enddo

      included(nTarget) = .true.
      
      
      end subroutine 
!=======================================================================
!     This function is similar to the Rosen_CreateSubset function.  The
!     only difference being this is intended for calculating the reverse probability. Thus it
!     will excluse nIndx from the subset.
      subroutine Rosen_CreateSubset_Reverse(nTarget, nIndx ,included)
      use Coords
      use ForceField
      use Constants
      use SimParameters
      implicit none
      integer, intent(in) :: nTarget, nIndx
      logical, intent(out) :: included(1:maxMol)

      logical:: inc_firstPass(1:maxMol)
      integer :: jIndx
      
      included = .false.
      inc_firstPass = .false.
      
!      Loop 1: Adds the neighbors of the target molecule to the list      
      do jIndx = 1, maxMol
         if(NeighborList(nTarget, jIndx)) then
            if(.not. isActive(jIndx)) then
              cycle
            endif   
            if(jIndx .eq. nIndx) then
              cycle
            endif                   
            included(jIndx) = .true.         
         endif
      enddo


!      nType = typeList(nTarget)
!      nMol = subIndxList(nTarget)
!      globIndx1 = molArray(nType)%mol(nMol)%globalIndx(1)
!      do jIndx = 1, maxMol
!         if(.not. isActive(jIndx)) then
!           cycle
!         endif
!         if(jIndx .eq. nIndx) then
!           cycle
!         endif 
!         jType = typeList(jIndx)
!         jMol = subIndxList(jIndx)
!         globIndx2 = molArray(jType)%mol(jMol)%globalIndx(1)
!         if(rPair(globIndx1, globIndx2)%p%r_sq .lt. rosenCut_sq) then
!           included(jIndx) = .true.         
!         endif
!      enddo




!      Loop 2: Adds the neighbors of the molecules added in the previous step to the list      
!      do jIndx = 1, maxMol
!       if(inc_firstPass(jIndx) .eqv. .true.) then
!         if(.not. isActive(jIndx)) then
!           cycle
!         endif       
!         do kIndx = 1, maxMol
!           if(NeighborList(jIndx, kIndx)) then
!             if(.not. isActive(kIndx)) then
!               cycle
!             endif           
!             if(kIndx .eq. nIndx) then
!               cycle
!             endif              
!             included(kIndx) = .true.         
!           endif
!         enddo
!       endif
!      enddo
      
      included(nTarget) = .true.
      
      end subroutine 


!=======================================================================
      subroutine Schedule_BranchedMol_Growth(nType,nPath,nAtomLoc,nAtom,nGrow,regrown,GrowFrom,GrowPrev, &
                                               GrowNum,GrowList,TorNum,TorList)
      use SimParameters
      use CBMC_Variables
      implicit none
      
      integer, intent(in) :: nType, nPath, nAtomLoc, nAtom
      integer :: nGrow
      integer :: GrowFrom(1:maxAtoms),GrowPrev(1:maxAtoms),GrowNum(1:maxAtoms),TorNum(1:maxAtoms)
      integer :: TorList(1:maxAtoms,1:maxBranches),GrowList(1:maxAtoms,1:maxBranches)
      logical :: regrown(1:maxAtoms)
      logical :: LGrowing
      integer :: atmCur, atmPrev, atmNext, curBonds, preBonds, TotalGrowing, OuterNum, OuterAtoms(1:maxAtoms), OuterPrev(1:maxAtoms)
      integer :: Incrmt, iBond, iu, iufrom, iCBunit, BondedAtoms(1:maxBranches), OuterTry
      real(dp) :: grnd
      
      GrowFrom = 0
      GrowPrev = 0
      GrowNum = 0
      TorNum = 0
      TorList = 0
      GrowList = 0
      nGrow = 1
      atmCur = nAtom
      GrowFrom(nGrow) = atmCur
      curBonds = topolArray(nType)%atom(atmCur)
      select case(curBonds)
      case(1) ! A Terminal has been selected, the whole molecule is regrown
         GrowPrev(nGrow) = 0
         TorNum(nGrow) = 0
         GrowNum(nGrow) = curBonds
         atmNext = pathArray(nType)%path(nPath,nAtomLoc + 1)
         GrowList(nGrow,GrowNum(nGrow)) = atmNext
         regrown(atmNext) = .false.
      case(2) ! A Linker has been selected, regrow to left or right randomly
         if (grnd() .gt. 0.5d0) then
            Incrmt = 1
         else
            Incrmt = -1
         endif
         atmPrev = pathArray(nType)%path(nPath,nAtomLoc - Incrmt)
         GrowPrev(nGrow) = atmPrev
         GrowNum(nGrow) = curBonds - 1
         atmNext = pathArray(nType)%path(nPath,nAtomLoc + Incrmt)
         GrowList(nGrow,GrowNum(nGrow)) = atmNext
         regrown(atmNext) = .false.
         preBonds = topolArray(nType)%atom(atmPrev)
         TorNum(nGrow) = preBonds - 1
         if (TorNum(nGrow) .gt. 0) then
            LGrowing = .false.
            call getRandomBondedAtoms(nType,atmPrev,atmCur,TorNum(nGrow),BondedAtoms,LGrowing)
            do iBond = 1,TorNum(nGrow)
               TorList(nGrow,iBond) = BondedAtoms(iBond)
            enddo
         endif
      case default ! A Hub has been selected, the selected path is kept and regrow the other paths of the hub
         if (nAtomLoc .eq. 1) then
            Incrmt = 1
         else
            Incrmt = -1
         endif  
         atmPrev = pathArray(nType)%path(nPath,nAtomLoc + Incrmt)
         GrowPrev(nGrow) = atmPrev
         GrowNum(nGrow) = curBonds - 1
         LGrowing = .true.
         call getRandomBondedAtoms(nType,atmCur,atmPrev,GrowNum(nGrow),BondedAtoms,LGrowing)
         do iBond = 1,GrowNum(nGrow)
            GrowList(nGrow,iBond) = BondedAtoms(iBond)
            regrown(BondedAtoms(iBond)) = .false.
         enddo
         preBonds = topolArray(nType)%atom(atmPrev)
         TorNum(nGrow) = preBonds - 1
         if (TorNum(nGrow) .gt. 0) then
            LGrowing = .false.
            call getRandomBondedAtoms(nType,atmPrev,atmCur,TorNum(nGrow),BondedAtoms,LGrowing)
            do iBond = 1,TorNum(nGrow)
               TorList(nGrow,iBond) = BondedAtoms(iBond)
            enddo
         endif
      end select
      TotalGrowing = 0 
      OuterNum = 0
      iufrom = GrowFrom(nGrow)
      do iCBunit = 1,GrowNum(nGrow)
         TotalGrowing = TotalGrowing + 1
         iu = GrowList(nGrow,iCBunit)
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
         GrowFrom(nGrow) = atmCur
         GrowPrev(nGrow) = atmPrev
         GrowNum(nGrow) = topolArray(nType)%atom(atmCur) - 1
         LGrowing = .true.
         call getRandomBondedAtoms(nType,atmCur,atmPrev,GrowNum(nGrow),BondedAtoms,LGrowing)
         do iBond = 1,GrowNum(nGrow)
            GrowList(nGrow,iBond) = BondedAtoms(iBond)
            regrown(BondedAtoms(iBond)) = .false.
         enddo
         preBonds = topolArray(nType)%atom(atmPrev)
         TorNum(nGrow) = preBonds - 1
         if (TorNum(nGrow) .gt. 0) then
            LGrowing = .false.
            call getRandomBondedAtoms(nType,atmPrev,atmCur,TorNum(nGrow),BondedAtoms,LGrowing)
            do iBond = 1,TorNum(nGrow)
               TorList(nGrow,iBond) = BondedAtoms(iBond)
            enddo
         endif
! Update list of "Outer" beads: Remove bead that was just grown from Outer list
         OuterAtoms(OuterTry) = OuterAtoms(OuterNum)
         OuterPrev(OuterTry) = OuterPrev(OuterNum)
         OuterNum = OuterNum - 1
! Add the new beads if they have more to be grown
         iufrom = atmCur
         do iCBunit = 1,GrowNum(nGrow)
            TotalGrowing = TotalGrowing + 1
            iu = GrowList(nGrow,iCBunit)
            iBond = topolArray(nType)%atom(iu)
            if (iBond .gt. 1) then
               OuterNum = OuterNum + 1
               OuterAtoms(OuterNum) = iu
               OuterPrev(OuterNum) = iufrom
            endif
         enddo
      enddo

      
      end subroutine

!=======================================================================
      subroutine FindAtomsFromPath(nType, regrown, pathNum, Atm4, Atm1, Atm2, Atm3)
      use Forcefield
      use SimParameters
      use Coords
      use CBMC_Variables
      implicit none
      integer, intent(in) :: Atm4, pathNum, nType
      logical, intent(in) :: regrown(1:maxAtoms)
      integer, intent(out) ::  Atm1, Atm2, Atm3

      integer :: i, atm4Pos
      integer :: atm4plus1, atm4minus1

      do i = 1, pathArray(nType)%pathMax(pathNum)
        if(pathArray(nType)%path(pathNum,i) .eq. Atm4 ) then
          atm4Pos = i
          exit
        endif
      enddo
  
      atm4plus1 = atm4Pos + 1
      atm4minus1 = atm4Pos - 1

      write(35,*) atm4plus1, atm4minus1
      write(35,*) "regrown:", regrown(:)

      if(atm4minus1 .lt. 1) then
        Atm3 = pathArray(nType)%path(pathNum,atm4Pos+1)
        Atm2 = pathArray(nType)%path(pathNum,atm4Pos+2)
        Atm1 = pathArray(nType)%path(pathNum,atm4Pos+3)
        return
      endif

      if(atm4plus1 .gt. pathArray(nType)%pathMax(pathNum)) then
        Atm3 = pathArray(nType)%path(pathNum,atm4Pos-1)
        Atm2 = pathArray(nType)%path(pathNum,atm4Pos-2)
        Atm1 = pathArray(nType)%path(pathNum,atm4Pos-3)
        return
      endif

      if( regrown(pathArray(nType)%path(pathNum,atm4plus1)) ) then
        Atm3 = pathArray(nType)%path(pathNum,atm4Pos+1)
        Atm2 = pathArray(nType)%path(pathNum,atm4Pos+2)
        Atm1 = pathArray(nType)%path(pathNum,atm4Pos+3)
        return
      endif


      if( regrown(pathArray(nType)%path(pathNum,atm4minus1)) ) then
        Atm3 = pathArray(nType)%path(pathNum,atm4Pos-1)
        Atm2 = pathArray(nType)%path(pathNum,atm4Pos-2)
        Atm1 = pathArray(nType)%path(pathNum,atm4Pos-3)
        return
      endif

      write(*,*) "ERROR in FindAtomsFromPath subroutine! Could not find the torsional members"
      write(*,*) "from the stored pathway for atom:", Atm4
      stop
      
      end subroutine 
!=======================================================================
      subroutine FindNextBondFromPath(nType, regrown, pathNum, Atm2, Atm1)
      use Forcefield
      use SimParameters
      use Coords
      use CBMC_Variables
      implicit none
      integer, intent(in) :: Atm2, pathNum, nType
      logical, intent(in) :: regrown(1:maxAtoms)
      integer, intent(out) ::  Atm1

      integer :: i, atm2Pos
      integer :: atm2plus1, atm2minus1

      do i = 1, pathArray(nType)%pathMax(pathNum)
        if(pathArray(nType)%path(pathNum,i) .eq. Atm2 ) then
          atm2Pos = i
          exit
        endif
      enddo
  
      atm2plus1 = atm2Pos + 1
      atm2minus1 = atm2Pos - 1

      if(atm2minus1 .lt. 1) then
        Atm1 = pathArray(nType)%path(pathNum,atm2Pos+1)
        return
      endif

      if(atm2plus1 .gt. pathArray(nType)%pathMax(pathNum)) then
        Atm1 = pathArray(nType)%path(pathNum,atm2Pos-1)
        return
      endif

      if( regrown(pathArray(nType)%path(pathNum,atm2plus1)) ) then
        Atm1 = pathArray(nType)%path(pathNum,atm2Pos+1)
        return
      endif


      if(.not. regrown(pathArray(nType)%path(pathNum,atm2minus1))) then
        Atm1 = pathArray(nType)%path(pathNum,atm2Pos-1)
        return
      endif

      write(*,*) "ERROR in FindNextBondFromPath subroutine! Could not find the bond members"
      write(*,*) "from the stored pathway for atom:", Atm2
      stop
      
      end subroutine 

!=======================================================================
      subroutine getRandomBondedAtoms(iType,atm1,atm2,BondedNum,BondedAtoms,LGrowing)
      use SimParameters
      use CBMC_Variables
      implicit none
      integer, intent(in) :: iType, atm1, atm2, BondedNum
      logical, intent(in) :: LGrowing
      integer :: BondedAtoms(1:maxBranches), iBond, iLast, iPath, iAtom, atmConnect

      
      BondedAtoms = 0
      iBond = 0
      select case(BondedNum)
      case(1) ! atm1 is a Linker
         do iPath = 1,pathArray(iType)%nPaths
            iLast = pathArray(iType)%pathMax(iPath)
            do iAtom = 2,iLast - 1
               if (atm1 .eq. pathArray(iType)%path(iPath,iAtom)) then
                  if (atm2 .eq. pathArray(iType)%path(iPath,iAtom+1)) then
                     BondedAtoms(1) = pathArray(iType)%path(iPath,iAtom-1)
                  else
                     BondedAtoms(1) = pathArray(iType)%path(iPath,iAtom+1)
                  endif
                  return
               endif
            enddo
         enddo
      case default ! atm1 is a Hub
         do iPath = 1,pathArray(iType)%nPaths
            if (atm1 .eq. pathArray(iType)%path(iPath,1)) then
               atmConnect = pathArray(iType)%path(iPath,2)
               if (atm2 .ne. atmConnect) then
                  iBond = iBond + 1
                  BondedAtoms(iBond) = atmConnect
                  if (iBond .eq. BondedNum) goto 58
               endif
            endif
            iLast = pathArray(iType)%pathMax(iPath)
            if (atm1 .eq. pathArray(iType)%path(iPath,iLast)) then
               atmConnect = pathArray(iType)%path(iPath,iLast - 1)
               if (atm2 .ne. atmConnect) then
                  iBond = iBond + 1
                  BondedAtoms(iBond) = atmConnect
                  if (iBond .eq. BondedNum) goto 58
               endif
            endif
         enddo
      end select
      write(*,*) "Error! Bonded Atom Not Found"
      stop
58    if (LGrowing) call RandomSort(BondedAtoms,BondedNum)
      
      end subroutine
!=======================================================================
      subroutine RandomSort(BondedAtoms,nBranch)
      use SimParameters
      use CBMC_Variables
      implicit none
      integer :: BondedAtoms(1:maxBranches),nBranch
      integer :: temp_store(1:maxBranches),zz,iBranch,itry,iut
      real(dp) :: grnd
      temp_store = 0
      do iBranch = 1,nBranch
         temp_store(iBranch) = BondedAtoms(iBranch)
      enddo
      iBranch = 0
      do zz = nBranch,1,-1
         itry = floor(grnd()*zz + 1d0)
         iut = temp_store(itry)
         iBranch = iBranch + 1
         BondedAtoms(iBranch) = iut
         temp_store(itry) = temp_store(zz)
      enddo
      if (iBranch .ne. nBranch) then
         write(*,*) "Error! Bonded Atom Not Sorted"
         stop
      endif
      end subroutine
!=======================================================================
      subroutine FindNextAngleFromPath(nType, regrown, pathNum, Atm3, Atm1, Atm2)
      use Forcefield
      use SimParameters
      use Coords
      use CBMC_Variables
      implicit none
      integer, intent(in) :: Atm3, pathNum, nType
      logical, intent(in) :: regrown(1:maxAtoms)
      integer, intent(out) ::  Atm1, Atm2

      integer :: i, atm3Pos
      integer :: atm3plus1, atm3minus1

      do i = 1, pathArray(nType)%pathMax(pathNum)
        if(pathArray(nType)%path(pathNum,i) .eq. Atm3 ) then
          atm3Pos = i
          exit
        endif
      enddo
  
      atm3plus1 = atm3Pos + 1
      atm3minus1 = atm3Pos - 1

      if(atm3minus1 .lt. 1) then
        Atm2 = pathArray(nType)%path(pathNum,atm3Pos+1)
        Atm1 = pathArray(nType)%path(pathNum,atm3Pos+2)
        return
      endif

      if(atm3plus1 .gt. pathArray(nType)%pathMax(pathNum)) then
        Atm2 = pathArray(nType)%path(pathNum,atm3Pos-1)
        Atm1 = pathArray(nType)%path(pathNum,atm3Pos-2)
        return
      endif

      if( regrown(pathArray(nType)%path(pathNum,atm3plus1)) ) then
        Atm2 = pathArray(nType)%path(pathNum,atm3Pos+1)
        Atm1 = pathArray(nType)%path(pathNum,atm3Pos+2)
        return
      endif


      if(.not. regrown(pathArray(nType)%path(pathNum,atm3minus1))) then
        Atm2 = pathArray(nType)%path(pathNum,atm3Pos-1)
        Atm1 = pathArray(nType)%path(pathNum,atm3Pos-2)
        return
      endif

      write(*,*) "ERROR in FindNextAngleFromPath subroutine! Could not find the torsional members"
      write(*,*) "from the stored pathway for atom:", Atm3
      stop
      
      end subroutine 
!=======================================================================
      subroutine ChooseRosenTrial(atomNumber, nType, trialPos, regrown, E_Trial, overlap, rosenRatio, rejMove)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions
      use CBMC_Variables
      implicit none
      logical, intent(in) :: overlap(:)
      integer, intent(in) :: atomNumber, nType
      type(SimpleAtomCoords), intent(in) :: trialPos(:)
      logical, intent(inout) :: regrown(:)
      real(dp), intent(inout) :: E_Trial(:)
      real(dp), intent(inout) :: rosenRatio

 
      logical,intent(out) :: rejMove
      integer :: nSel, iRosen
      real(dp) :: E_Min, ProbRosen(1:maxRosenTrial), rosenNorm, grnd
      real(dp) :: ranNum, sumInt

      rejMove = .false.
      E_Min = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
         ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Min))         
      enddo
      if(all(ProbRosen .le. 0d0)) then
        rejMove = .true.
        return
      endif
      rosenNorm = sum(ProbRosen)
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbRosen(nSel)
      enddo
      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif
      rosenRatio = rosenRatio*ProbRosen(nSel)/rosenNorm
      regrown(atomNumber) = .true.
      newMol%x(atomNumber) = trialPos(nSel)%x 
      newMol%y(atomNumber) = trialPos(nSel)%y 
      newMol%z(atomNumber) = trialPos(nSel)%z

      end subroutine
!=======================================================================
      subroutine ChooseRosenTrial_Branched(nToGrow, atmGrow, nType, trialPos_Branched, regrown, E_Trial, wBending, &
                                           wTorsion, overlap, rosenRatio, rejMove)
      use SimParameters
      use Coords
      use ForceField
      use Constants
      use CBMC_Variables
      implicit none
      logical, intent(in) :: overlap(:)
      integer, intent(in) :: nToGrow, atmGrow(:), nType
      real(dp), intent(in) :: wBending(:), wTorsion(:)
      type(SimpleAtomCoords), intent(in) :: trialPos_Branched(:,:)
      logical, intent(inout) :: regrown(:)
      real(dp), intent(inout) :: E_Trial(:)
      real(dp), intent(inout) :: rosenRatio 
      logical,intent(out) :: rejMove
      
      integer :: nSel, iRosen, iAtom
      real(dp) :: E_Min, ProbRosen(1:maxRosenTrial), rosenNorm, grnd
      real(dp) :: ranNum, sumInt
      
      rejMove = .false.
      E_Min = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
         ProbRosen(iRosen) = wBending(iRosen) * wTorsion(iRosen) * exp(-beta*(E_Trial(iRosen)-E_Min))         
      enddo
      if(all(ProbRosen .le. 0d0)) then
        rejMove = .true.
        return
      endif
      rosenNorm = sum(ProbRosen)
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbRosen(nSel)
      enddo
      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif
      rosenRatio = rosenRatio*ProbRosen(nSel)/rosenNorm
      do iAtom = 1, nToGrow
        regrown(atmGrow(iAtom)) = .true.
        newMol%x(atmGrow(iAtom)) = trialPos_Branched(nSel,iAtom)%x 
        newMol%y(atmGrow(iAtom)) = trialPos_Branched(nSel,iAtom)%y 
        newMol%z(atmGrow(iAtom)) = trialPos_Branched(nSel,iAtom)%z
      enddo
      end subroutine
!=======================================================================
      end module
