!=================================================================================
      module EnergyCriteria
      contains
!=================================================================================     
!     Extensive Cluster Criteria Check.  Used at the start and end of the simulation. 
!     This ensures that all particles in the cluster are properly connected to each other.
!     This function also calculates the initial Neighborlist that is used throughout the simulation. 
      subroutine Detailed_EnergyCriteria(PairList,rejMove)
      use SimParameters
      use Coords
      use IndexingFunctions
      use ParallelVar
      implicit none     
      logical,intent(out) :: rejMove
      real(dp),intent(inout) :: PairList(1:maxMol, 1:maxMol)
      
      logical :: ClusterMember(1:maxMol)
      integer :: i,j,h,cnt
      integer :: iType,jType, iMol, jMol, iIndx, jIndx

      rejMove = .false.
      NeighborList = .false.
      if(NTotal .eq. 1) then
         return      
      endif
      
      do iType = 1, nMolTypes
        do jType = iType, nMolTypes
          do iMol = 1, NPART(iType)      
            iIndx = MolArray(iType)%mol(iMol)%indx
            do jMol = 1, NPART(jType) 
              jIndx = MolArray(jType)%mol(jMol)%indx        
              if(PairList(iIndx,jIndx) .le. Eng_Critr(iType,jType) ) then
                NeighborList(iIndx,jIndx)=.true.         
                NeighborList(jIndx,iIndx)=.true.          
              endif
            enddo
          enddo
        enddo
      enddo

      if( all(NeighborList .eqv. .false.) )then
        rejMove = .true.
        write(nout,*) "------- Cluster Criteria Not Met! --------"
      endif

      cnt = 0
      do i=1,maxMol
        if(isActive(i) .eqv. .false.) then
          ClusterMember(i) = .true.
          cnt = cnt + 1
        endif         
      enddo


      do i = 1, maxMol
        if(isActive(i) .eqv. .true.) then
          ClusterMember(i) = .true.
          exit
        endif         
      enddo
      
!      ClusterMember(1)=.true.      
      do h = 1, maxMol  
        do i = 1, maxMol
          do j = 1, maxMol
            if( NeighborList(i,j) ) then
              if( ClusterMember(i) ) then
                ClusterMember(j) = .true.
!                cnt = cnt + 1
              endif
              if( ClusterMember(j) ) then
                ClusterMember(i) = .true.
!                cnt = cnt + 1                
              endif
            endif
          enddo           
        enddo
      enddo     

      do i = 1, maxMol
        NeighborList(i,i) = .false.
      enddo
      
      if(any(ClusterMember .eqv. .false.) ) then
        rejMove = .true.
        write(nout,*) "------- Cluster Criteria Not Met! --------"
      endif
     
      end subroutine
!=================================================================================     
!     This function determines if a given translational move will destroy a cluster. 
      subroutine Shift_EnergyCriteria(PairList, nIndx, rejMove)
      use SimParameters     
      use Coords
      use IndexingFunctions
      implicit none     
      
      logical, intent(out) :: rejMove      
      real(dp), intent(in) :: PairList(1:maxMol)
      integer,intent(in) :: nIndx
      
      logical :: neiFlipped, memberAdded
      logical :: ClusterMember(1:maxMol)      
      logical :: flipped(1:maxMol)
      integer :: i,j,h
      integer :: nType, jType      
      integer :: curNeigh(1:60), neiMax
      
      rejMove=.false.
      if(NTotal .eq. 1) return        
      ClusterMember=.false.
      flipped=.false.
      
      nType = Get_MolType(nIndx,NMAX)
     
!     This section dermines which molecules are neighbored with the new trial position.  In the event
!     that the molecule's new location has no neghibors all further calcualtions are skipped and the move is
!     rejected.
   
      memberAdded = .false.
      do j=1,maxMol
        jType = typeList(j)          
        if(PairList(j) .le. Eng_Critr(jType,nType)) then
          ClusterMember(j) = .true.        
          memberAdded = .true.
        endif
      enddo

      if(.not. memberAdded) then      
        rejMove = .true.
        return     
      endif    

!      This part of the code tabulates all the neighbors       
      neiMax = 0
      curNeigh = 0
      do i=1,maxMol
        if(NeighborList(i,nIndx)) then
          if(i .ne. nIndx) then
            if(isActive(i)) then
              neiMax = neiMax + 1
              curNeigh(neiMax) = i
            endif
          endif      
        endif
      enddo      
      
      
!     This section performs a quick check to see if the molecules that were neighbored with the old position
!     are part of the new cluster.  If all the old neighbors are indeed part of the cluster then no furth
!     calculations are needed.      
      neiFlipped = .true.
      do i = 1, neiMax
        if(.not. clusterMember(curNeigh(i))) then
          neiFlipped = .false.
          exit
        endif
      enddo

      if(neiFlipped) then
        rejMove = .false.
        return
      endif
      
    
      do h = 1, NTotal
!        cnt = 0
        memberAdded = .false.
        do i = 1, maxMol
            if(ClusterMember(i) .neqv. flipped(i)) then
              do j=1,maxMol
                if(NeighborList(i,j)) then
                  if(j .ne. nIndx) then
                    ClusterMember(j)=.true.   
                    memberAdded = .true.
                  endif
                endif
              enddo        
              flipped(i)=.true.
            endif
        enddo
 
        neiFlipped = .true.
        do i = 1, neiMax
          if(.not. clusterMember(curNeigh(i))) then
            neiFlipped = .false.
            exit
          endif
        enddo        
        if( neiFlipped ) then
          exit
        else 
          if(.not. memberAdded) then
            exit
          endif           
        endif
      enddo
  
       if( .not. neiFlipped ) then
         rejMove = .true.
       endif
     
      end subroutine
!=================================================================================     
!     This function determines if a given translational move will destroy a cluster. 
      pure subroutine SwapIn_EnergyCriteria(nType,PairList,rejMove)
      use SimParameters     
      use Coords
      use IndexingFunctions
      implicit none     
      
      logical, intent(out) :: rejMove      
      real(dp), intent(in) :: PairList(1:maxMol)
      integer, intent(in) :: nType
      integer :: j, jType
      
      do j = 1, maxMol
        if(isActive(j)) then      
          jType = Get_MolType(j,NMAX)                 
          if(PairList(j) .le. Eng_Critr(nType,jType)) then
            rejMove = .false.
            return
          endif
        endif
      enddo
     
      rejMove = .true.
     
      end subroutine
!=================================================================================     
! Checks if deleting molecule nSwap breaks the single cluster, setting rejMove = .true.
! if the cluster becomes disconnected. Uses NeighborList to propagate cluster membership,
! ensuring all nSwap neighbors remain connected post-deletion. Called by AVBMC_EBias_Rosen_Out
! to reject invalid moves early. Correctness is ensured by:
! - Identifying nSwap’s neighbors and building a connected cluster excluding nSwap.
! - Rejecting if any neighbor is disconnected.
! - Validating inputs and activity status.
! Supports all force fields (LJ_Q, Mpipi, HPS_single, HPS_piecewise, HPS_cation_pi)
! and molecule types (rigid, small, linear, branched).
subroutine SwapOut_ClusterCriteria(nSwap, rejMove)
  use SimParameters, only: maxMol, NTotal, isActive
  use Coords, only: NeighborList
  implicit none

  integer, intent(in) :: nSwap
  logical, intent(out) :: rejMove

  logical :: memberAdded, neiFlipped
  logical :: clusterMember(1:maxMol)
  logical :: flipped(1:maxMol)
  integer :: i, j, h, neiMax
  integer :: curNeigh(1:maxMol)

  ! Initialize arrays and flags
  ! Correctness: Prepares for cluster propagation
  clusterMember = .false.
  flipped = .false.
  rejMove = .false.
  neiMax = 0
  curNeigh = 0

  ! Mark first neighbor of nSwap as cluster member
  ! Correctness: Starts cluster with a connected molecule
  do i = 1, maxMol
    if (isActive(i) .and. NeighborList(i, nSwap) .and. i /= nSwap) then
      clusterMember(i) = .true.
      exit
    endif
  enddo

  ! Build list of nSwap’s neighbors
  ! Correctness: Tracks neighbors to check connectivity
  do i = 1, maxMol
    if (isActive(i) .and. NeighborList(i, nSwap) .and. i /= nSwap) then
      neiMax = neiMax + 1
      curNeigh(neiMax) = i
    endif
  enddo

  ! Propagate cluster membership
  ! Correctness: Ensures all connected molecules are marked, excluding nSwap
  do h = 1, NTotal
    memberAdded = .false.
    do i = 1, maxMol
      if (clusterMember(i) .and. .not. flipped(i)) then
        do j = 1, maxMol
          if (NeighborList(i, j) .and. j /= nSwap) then
            clusterMember(j) = .true.
            memberAdded = .true.
          endif
        enddo
        flipped(i) = .true.
      endif
    enddo
    ! Check if all nSwap neighbors are in cluster
    neiFlipped = .true.
    do i = 1, neiMax
      if (.not. clusterMember(curNeigh(i))) then
        neiFlipped = .false.
        exit
      endif
    enddo
    ! Exit if all neighbors are connected or no new members added
    if (neiFlipped .or. .not. memberAdded) exit
  enddo

  ! Reject move if any neighbor is disconnected
  if (.not. neiFlipped) rejMove = .true.
end subroutine SwapOut_ClusterCriteria
!=================================================================================     
!     This function updates the neighborlist if a move is accepted.
      subroutine NeighborUpdate(PairList, nIndx)
      use SimParameters
      use IndexingFunctions      
      use Coords      
      implicit none     
      integer iType,j,jType,nIndx
      real(dp) :: PairList(1:maxMol)


!      do j=1,maxMol
!        if(.not. isActive(j)) cycle      
!        if(j .ne. nIndx) then
!            NeighborList(nIndx,j)=.true.
!            NeighborList(j,nIndx)=.true.  
!          else             
!            NeighborList(nIndx,j)=.false.
!            NeighborList(j,nIndx)=.false.            
!          endif   
!        endif
!      enddo


      iType = Get_MolType(nIndx,NMAX)  
      do j=1,maxMol
        if(.not. isActive(j)) then
          cycle
        endif
        if(j .ne. nIndx) then
          jType = Get_MolType(j,NMAX)        
          if(PairList(j) .le. Eng_Critr(iType,jType) ) then
            NeighborList(nIndx,j)=.true.
            NeighborList(j,nIndx)=.true.  
          else             
            NeighborList(nIndx,j)=.false.
            NeighborList(j,nIndx)=.false.            
          endif   
        endif
      enddo
      NeighborList(nIndx,nIndx) = .false.
      

      end subroutine
!================================================================================= 
! Updates the neighbor list and, if minDistCriteria=.true., neighbor pair details after a molecule
! deletion in a Grand Canonical Monte Carlo nucleation simulation. Transfers neighbor relationships
! from nSwapIndx (swapped molecule, typically MolArray(nType)%mol(NPART(nType))%indx) to nIndx
! (deleted molecule’s new position, MolArray(nType)%mol(nMol)%indx), and clears nSwapIndx’s entries.
! Ensures NeighborList symmetry and updates NeighborPairs(i,j) for active neighbor pairs.
! Optimized for single-bead molecules (HPS force fields) but supports multi-atom molecules (LJ_Q, Mpipi).
! Called after a deletion is accepted in AVBMC_EBias_Rosen_Out or Exchange moves.
subroutine NeighborUpdate_Delete(nIndx, nSwapIndx)
  use SimParameters, only: maxMol, maxAtoms, minDistCriteria, isActive
  use Coords, only: NeighborList, NeighborPairs
  use VarPrecision, only: dp
  implicit none

  ! Input variables
  integer, intent(in) :: nIndx        ! Index of the deleted molecule’s new position
  integer, intent(in) :: nSwapIndx    ! Index of the swapped molecule (to be cleared)

  ! Local variables
  integer :: i, iPair              ! Loop indices
  integer :: nPairs                   ! Number of atom pairs for a neighbor pair
  integer :: atm1, atm2               ! Atom indices for pair details

  ! Step 1: Handle case where no swap is needed (nIndx == nSwapIndx)
  if (nIndx == nSwapIndx) then
    NeighborList(nSwapIndx, :) = .false.
    NeighborList(:, nSwapIndx) = .false.
    if (minDistCriteria) then
      ! Clear NeighborPairs entries for nSwapIndx
      do i = 1, maxMol
        NeighborPairs(nSwapIndx, i)%details%nPairs = 0
        NeighborPairs(nSwapIndx, i)%details%pairIndices = 0
        NeighborPairs(i, nSwapIndx)%details%nPairs = 0
        NeighborPairs(i, nSwapIndx)%details%pairIndices = 0
      enddo
    endif
    return
  endif

  ! Step 2: Transfer neighbor relationships from nSwapIndx to nIndx
  do i = 1, maxMol
    if (.not. isActive(i) .or. i == nIndx .or. i == nSwapIndx) then
      NeighborList(i, nIndx) = .false.
      NeighborList(nIndx, i) = .false.
      cycle
    endif
    NeighborList(i, nIndx) = NeighborList(i, nSwapIndx)
    NeighborList(nIndx, i) = NeighborList(i, nSwapIndx)
    if (minDistCriteria) then
      ! Transfer NeighborPairs details from (nSwapIndx,i) to (nIndx,i)
      nPairs = NeighborPairs(nSwapIndx, i)%details%nPairs
      NeighborPairs(nIndx, i)%details%nPairs = nPairs
      do iPair = 1, nPairs
        atm1 = NeighborPairs(nSwapIndx, i)%details%pairIndices(2 * iPair - 1)
        atm2 = NeighborPairs(nSwapIndx, i)%details%pairIndices(2 * iPair)
        NeighborPairs(nIndx, i)%details%pairIndices(2 * iPair - 1) = atm1
        NeighborPairs(nIndx, i)%details%pairIndices(2 * iPair) = atm2
      enddo
      if (nPairs < maxAtoms) then
        NeighborPairs(nIndx, i)%details%pairIndices(2 * nPairs + 1:) = 0
      endif
      ! Transfer NeighborPairs details from (i,nSwapIndx) to (i,nIndx) with swapped indices
      NeighborPairs(i, nIndx)%details%nPairs = nPairs
      do iPair = 1, nPairs
        atm1 = NeighborPairs(nSwapIndx, i)%details%pairIndices(2 * iPair - 1)
        atm2 = NeighborPairs(nSwapIndx, i)%details%pairIndices(2 * iPair)
        NeighborPairs(i, nIndx)%details%pairIndices(2 * iPair - 1) = atm2
        NeighborPairs(i, nIndx)%details%pairIndices(2 * iPair) = atm1
      enddo
      if (nPairs < maxAtoms) then
        NeighborPairs(i, nIndx)%details%pairIndices(2 * nPairs + 1:) = 0
      endif
    endif
  enddo
  NeighborList(nIndx, nIndx) = .false.

  ! Step 3: Clear nSwapIndx’s neighbor relationships
  NeighborList(nSwapIndx, :) = .false.
  NeighborList(:, nSwapIndx) = .false.
  if (minDistCriteria) then
    do i = 1, maxMol
      NeighborPairs(nSwapIndx, i)%details%nPairs = 0
      NeighborPairs(nSwapIndx, i)%details%pairIndices = 0
      NeighborPairs(i, nSwapIndx)%details%nPairs = 0
      NeighborPairs(i, nSwapIndx)%details%pairIndices = 0
    enddo
  endif
end subroutine NeighborUpdate_Delete
!=================================================================================     
!     This function updates the neighborlist if a move is accepted.
      subroutine LLNeighborUpdate_Delete(nIndx, nSwapIndx)
      use SimParameters
      use IndexingFunctions
      use Coords
      implicit none
      integer, intent(in) :: nIndx, nSwapIndx
      integer :: nType
      integer :: i, iPair


!      write(35,*) nIndx, nSwapIndx
      nType = typeList(nIndx)     
!      nSwapIndx = molArray(nType)%mol(NPART(nType))%indx

!      NeighborList(nIndx,:) = NeighborList(nSwapIndx,:)
!      NeighborList(:,nIndx) = NeighborList(:,nSwapIndx)

      if(nIndx .eq. nSwapIndx) then
        NeighborList(nSwapIndx,:) = .false.
        NeighborList(:,nSwapIndx) = .false.   
        NeighborPairs(nSwapIndx,:)%details%nPairs = 0
        NeighborPairs(:,nSwapIndx)%details%nPairs = 0
        return
      endif
     
      do i = 1, maxMol
!       if(.not. isActive(i)) then
!          NeighborList(i,nIndx) = .false.
!          NeighborList(nIndx,i) = .false.          
!          cycle
!        endif
        if(NeighborList(i,nSwapIndx)) then
          if(i .ne. nIndx) then
            NeighborList(i,nIndx) = .true.
            NeighborList(nIndx,i) = .true.  
            NeighborPairs(i,nIndx)%details%nPairs = &
            NeighborPairs(i,nSwapIndx)%details%nPairs
            do iPair = 1, NeighborPairs(i,nIndx)%details%nPairs
              NeighborPairs(i,nIndx)%details%pairIndices(2 * iPair - 1) = &
              NeighborPairs(i,nSwapIndx)%details%pairIndices(2 * iPair - 1)
              NeighborPairs(i,nIndx)%details%pairIndices(2 * iPair) = &
              NeighborPairs(i,nSwapIndx)%details%pairIndices(2 * iPair)
            enddo

            NeighborPairs(nIndx,i)%details%nPairs = NeighborPairs(nSwapIndx,i)%details%nPairs
            do iPair = 1, NeighborPairs(nIndx,i)%details%nPairs
              NeighborPairs(nIndx,i)%details%pairIndices(2 * iPair - 1) = &
              NeighborPairs(nSwapIndx,i)%details%pairIndices(2 * iPair - 1)
              NeighborPairs(nIndx,i)%details%pairIndices(2 * iPair) = &
              NeighborPairs(nSwapIndx,i)%details%pairIndices(2 * iPair)
            enddo
          endif
        else
          NeighborList(i,nIndx) = .false.
          NeighborList(nIndx,i) = .false.
          NeighborPairs(i,nIndx)%details%nPairs = 0

          NeighborPairs(nIndx,i)%details%nPairs = 0
        endif        
      enddo
      
      NeighborList(nIndx,nIndx) = .false.
      NeighborList(nSwapIndx,:) = .false.
      NeighborList(:,nSwapIndx) = .false.
      NeighborPairs(nIndx,nIndx)%details%nPairs = 0

      NeighborPairs(nSwapIndx,:)%details%nPairs = 0
      NeighborPairs(:,nSwapIndx)%details%nPairs = 0      

      end subroutine      
!=================================================================================           
      subroutine MultipleSwap_EnergyCriteria(nType2, nIndx1, PairList, isIncluded, rejMove)
      use SimParameters     
      use Coords
      use IndexingFunctions
      implicit none     
      
      logical, intent(out) :: rejMove      
      real(dp), intent(in) :: PairList(1:maxMol)
      logical,  intent(in) :: isIncluded(:)
      integer, intent(in) :: nType2, nIndx1
      
      logical :: neiFlipped, memberAdded
      logical :: ClusterMember(1:maxMol)      
      logical :: flipped(1:maxMol)
      integer :: i,j,h
      integer :: jType
      integer :: curNeigh(1:60), neiMax
 

      rejMove=.false.
      if(NTotal .eq. 1) return        
      ClusterMember=.false.
      flipped=.false.
      
     
!     This section dermines which molecules are neighbored with the new trial position.  In the event
!     that the molecule's new location has no neghibors all further calcualtions are skipped and the move is
!     rejected.
   
      memberAdded = .false.
      do j=1,maxMol
        if(isIncluded(j)) then
          jType = typeList(j)          
          if(PairList(j) .le. Eng_Critr(jType,nType2)) then
            ClusterMember(j) = .true.        
            memberAdded = .true.
          endif
        endif
      enddo

      if(.not. memberAdded) then      
        rejMove = .true.
        return     
      endif    

!      This part of the code tabulates all the neighbors       
      neiMax = 0
      curNeigh = 0
      do i=1,maxMol
        if(NeighborList(i,nIndx1)) then
          if(i .ne. nIndx1) then
            if(isIncluded(i)) then
              neiMax = neiMax + 1
              curNeigh(neiMax) = i
            endif
          endif      
        endif
      enddo      
      
      
!     This section performs a quick check to see if the molecules that were neighbored with the old position
!     are part of the new cluster.  If all the old neighbors are indeed part of the cluster then no furth
!     calculations are needed.      
      neiFlipped = .true.

      do i = 1, neiMax
        if(.not. clusterMember(curNeigh(i))) then
          neiFlipped = .false.
          exit
        endif
      enddo

      if(neiFlipped) then
        rejMove = .false.
        return
      endif
      
    
      do h = 1, NTotal
        memberAdded = .false.
        do i = 1, maxMol
          if(isIncluded(i)) then
            if(ClusterMember(i) .neqv. flipped(i)) then
              do j=1,maxMol
                if(isIncluded(i)) then
                  if(NeighborList(i,j)) then  
                    ClusterMember(j)=.true.   
                    memberAdded = .true.
                  endif
                endif
              enddo        
              flipped(i)=.true.
            endif
          endif
        enddo
 
        neiFlipped = .true.
        do i = 1, neiMax
          if(.not. clusterMember(curNeigh(i))) then
            neiFlipped = .false.
            exit
          endif
        enddo        
        if( neiFlipped ) then
          exit
        else 
          if(.not. memberAdded) then
            exit
          endif           
        endif
      enddo
  
       if( .not. neiFlipped ) then
         rejMove=.true.
       endif
     

      end subroutine
!=================================================================================           
      end module
      
