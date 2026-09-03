!=================================================================================
      module DistanceCriteria

!      logical, allocatable :: ClusterMember(:), flipped(:) 
      contains
!=================================================================================     
!     Extensive Cluster Criteria Check.  Used at the start and end of the simulation. 
!     This ensures that all particles in the cluster are properly connected to each other.
!     This function also calculates the initial Neighborlist that is used throughout the simulation. 
      subroutine Detailed_DistanceCriteria(PairList, rejMove)
      use Coords
      use IndexingFunctions
      use ParallelVar
      use SimParameters
      implicit none     
      logical, intent(out) :: rejMove
      real(dp), intent(inout) :: PairList(:, :)
      
!      logical :: ClusterMember(1:maxMol)
      integer :: h,cnt
      integer :: iType,jType, iMol, jMol, iIndx, jIndx
      integer :: globIndx1, globIndx2
      logical :: ClusterMember(1:maxMol)      


      rejMove = .false.
      NeighborList = .false.
      if(NTotal .eq. 1) then
         return      
      endif


      do iType = 1, nMolTypes
        do iMol = 1 ,NPART(iType)      
          globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(1)
          iIndx = MolArray(iType)%mol(iMol)%indx
          do jType = 1, nMolTypes
            do jMol = 1 ,NPART(jType) 
              globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
              jIndx = MolArray(jType)%mol(jMol)%indx
              if(jIndx .eq. iIndx) then
                cycle
              endif        
              if(PairList(iIndx, jIndx) .lt. Dist_Critr_sq ) then
                NeighborList(iIndx, jIndx) = .true.         
                NeighborList(jIndx, iIndx) = .true.   
              else
                NeighborList(iIndx, jIndx) = .false.         
                NeighborList(jIndx, iIndx) = .false.   
              endif
            enddo
          enddo
        enddo
      enddo

      if( all(NeighborList .eqv. .false.) ) then
        write(nout,*) "------- Cluster Criteria Not Met!!!! --------"
        rejMove = .true.
        return
      endif



          
      do iType = 1, nMolTypes
        if(NPART(iType) .gt. 0) then
          iIndx = MolArray(iType)%mol(1)%indx
          ClusterMember(iIndx) = .true.
          exit
        endif
      enddo

      do h = 1, maxMol  
        do iIndx = 1, maxMol
          do jIndx = 1, maxMol
            if( NeighborList(iIndx, jIndx) ) then
              if( ClusterMember(iIndx) ) then
                ClusterMember(jIndx) = .true.
!                cnt = cnt + 1
              endif
              if( ClusterMember(jIndx) ) then
                ClusterMember(iIndx) = .true.
!                cnt = cnt + 1                
              endif
            endif
          enddo           
        enddo
      enddo     

      cnt = 0
      do iType = 1, nMolTypes
        if(NPART(iType) .lt. NMAX(iType)) then
          do iMol = NPART(iType)+1, NMAX(iType) 
            iIndx = MolArray(iType)%mol(iMol)%indx
            ClusterMember(iIndx) = .true.
            cnt = cnt + 1
          enddo
        endif   
      enddo
      
      if( any(ClusterMember .eqv. .false.) ) then
        rejMove = .true.
        write(nout,*) "------- Cluster Criteria Not Met!!!! --------"
      endif
     
      end subroutine
!=================================================================================     
! Determines if a translational move destroys a single cluster in a Grand Canonical Monte Carlo
! nucleation simulation. Compares neighbors of the molecule’s old and new positions using PairList
! (distances for new position) and NeighborList (old position neighbors). Rejects the move if the
! new position’s cluster doesn’t include all old neighbors. Supports distCriteria (atom 1 distances)
! and minDistCriteria (all atom distances). NeighborDetailsNew is optional for future minDistCriteria
! enhancements (e.g., using pairIndices for neighbor pairs).
subroutine Shift_DistanceCriteria(PairList, nIndx, rejMove, NeighborDetailsNew)
  USE VarPrecision, ONLY: dp
  USE SimParameters, ONLY: maxMol, NTotal, isActive, Dist_Critr_sq                           
  USE Coords, ONLY: molArray, NeighborList, typeList, subIndxList
  USE CoordinateTypes, ONLY: NeighborDetails
  implicit none

  ! Input/Output variables
  real(dp), intent(in) :: PairList(:)
  integer, intent(in) :: nIndx
  logical, intent(out) :: rejMove
  type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)

  ! Local variables
  logical :: neiFlipped, memberAdded
  logical :: ClusterMember(1:maxMol)      ! Tracks molecules in the new cluster
  logical :: flipped(1:maxMol)            ! Tracks processed molecules
  integer :: iIndx, jIndx, h, i
  integer :: nType, nMol, globIndx1
  integer :: neiMax                       ! Number of old position neighbors
  integer, allocatable :: curNeigh(:)     ! Indices of old position neighbors

  ! Step 1: Initialize variables
  rejMove = .false.
  if (NTotal == 1) return
  ClusterMember = .false.
  flipped = .false.
  neiMax = 0

  ! Allocate curNeigh dynamically based on maxMol
  allocate(curNeigh(1:maxMol))
  curNeigh = 0

  ! Get molecule info
  nType = typeList(nIndx)
  nMol = subIndxList(nIndx)
  globIndx1 = molArray(nType)%mol(nMol)%globalIndx(1)

  ! Step 2: Identify neighbors of the new trial position
  memberAdded = .false.
  do jIndx = 1, maxMol
    if (.not. isActive(jIndx)) cycle
    if (nIndx == jIndx) cycle
    if (PairList(jIndx) < Dist_Critr_sq) then
      ClusterMember(jIndx) = .true.
      memberAdded = .true.
    endif
  enddo

  ! Reject move if no neighbors in new position
  if (.not. memberAdded) then
    rejMove = .true.
    deallocate(curNeigh)
    return
  endif

  ! Step 3: Tabulate neighbors of the old position
  do iIndx = 1, maxMol
    if (NeighborList(iIndx, nIndx)) then
      if (iIndx /= nIndx .and. isActive(iIndx)) then
        neiMax = neiMax + 1
        curNeigh(neiMax) = iIndx
      endif
    endif
  enddo

  ! Step 4: Check if all old neighbors are in the new cluster
  neiFlipped = .true.
  do i = 1, neiMax
    if (.not. ClusterMember(curNeigh(i))) then
      neiFlipped = .false.
      exit
    endif
  enddo

  if (neiFlipped) then
    rejMove = .false.
    deallocate(curNeigh)
    return
  endif

  ! Step 5: Expand the cluster to include all connected molecules
  do h = 1, NTotal
    memberAdded = .false.
    do iIndx = 1, maxMol
      if (ClusterMember(iIndx) .neqv. flipped(iIndx)) then
        do jIndx = 1, maxMol
          if (NeighborList(iIndx, jIndx)) then
            if (jIndx /= nIndx) then
              ClusterMember(jIndx) = .true.
              memberAdded = .true.
            endif
          endif
        enddo
        flipped(iIndx) = .true.
      endif
    enddo

    neiFlipped = .true.
    do i = 1, neiMax
      if (.not. ClusterMember(curNeigh(i))) then
        neiFlipped = .false.
        exit
      endif
    enddo
    if (neiFlipped .or. .not. memberAdded) then
      exit
    endif
  enddo

  ! Step 6: Reject move if old neighbors are not in the final cluster
  if (.not. neiFlipped) rejMove = .true.

  deallocate(curNeigh)
end subroutine Shift_DistanceCriteria
!=================================================================================     
!     This function determines if removing a particle from the cluster will result in the destruction of the cluster criteria. 
      subroutine SwapOut_DistanceCriteria(nSwap, rejMove)
      use SimParameters     
      use Coords
      use IndexingFunctions
      implicit none     
      
      logical, intent(out) :: rejMove
      integer, intent(inout) :: nSwap
      integer :: iIndx, jIndx, h, cnt
      logical :: ClusterMember(1:maxMol)      
      logical :: flipped(1:maxMol)


      rejMove = .false.
      ClusterMember = .false.
      flipped = .false.

      do iIndx = 1, maxMol
        if( .not. isActive(iIndx) ) then      
          ClusterMember(iIndx) = .true.
          flipped(iIndx) = .true.
        endif
      enddo

!      In order to initialize the cluster criteria search, a single particle must be chosen as the starting point.
      do iIndx = 1, maxMol
        if( isActive(iIndx) ) then      
          if( nSwap .ne. iIndx ) then
            ClusterMember(iIndx) = .true.
            exit
          endif
        endif
      enddo

      cnt = 0
      do iIndx = 1, maxMol
        if( isActive(iIndx) .eqv. .false. ) then
          ClusterMember(iIndx) = .true.      
          flipped(iIndx) = .true.         
          cnt = cnt + 1
        endif
      enddo
      
      do h = 1, maxMol
        do iIndx = 1, maxMol
          if( iIndx .ne. nSwap ) then     
            if( ClusterMember(iIndx) .neqv. flipped(iIndx) ) then
              do jIndx = 1, maxMol
                if( NeighborList(iIndx, jIndx) ) then
                  ClusterMember(jIndx) = .true. 
                  cnt = cnt + 1
                endif
              enddo
              flipped(iIndx) = .true.
            endif
          endif
        enddo
        if(cnt .eq. maxMol-1) then
          exit
        endif         
      enddo
  
  
      ClusterMember(nSwap) = .true.
       
      if( any(ClusterMember .eqv. .false.) ) then
        rejMove = .true.
      endif      
     
      end subroutine
!=================================================================================     
! Updates the neighbor list and optionally neighbor pair details for a moved molecule based on the
! distance cluster criterion (distCriteria or minDistCriteria) after an accepted move in a Grand Canonical
! Monte Carlo nucleation simulation. Sets NeighborList(nIndx, j) and NeighborList(j, nIndx) to .true.
! if the minimum squared distance (PairList(j)) is less than Dist_Critr_sq, otherwise sets to .false.
! If NeighborDetailsNew is present, updates NeighborPairs(i,j) and NeighborPairs(j,i) with atom pair indices
! in pre-allocated pairIndices arrays, maintaining symmetry. Handles useDistStore=.true. by calling
! NeighborUpdate_Distance_PairStore. Optimized for single-bead molecules (HPS force fields) but supports
! multi-atom molecules (LJ_Q, Mpipi). Called after a move is accepted in Translation or LLTranslation moves.
subroutine NeighborUpdate_Distance(PairList, nIndx, NeighborDetailsNew)
  use SimParameters, only: maxMol, Dist_Critr_sq, isActive, maxAtoms
  use DistanceCriteria_PairStore, only: NeighborUpdate_Distance_PairStore
  use VarPrecision, only: dp
  use CoordinateTypes, only: NeighborDetails
  use Coords, only: NeighborPairs, NeighborList
  USE PairStorage, ONLY: useDistStore
  implicit none

  ! Input variables
  integer, intent(in) :: nIndx                                    ! Index of the moved molecule
  real(dp), intent(in) :: PairList(maxMol)                       ! Min squared distance between nIndx and other molecules
  type(NeighborDetails), intent(in), optional :: NeighborDetailsNew(maxMol)  ! New atom pair details (for minDistCriteria)

  ! Local variables
  integer :: j, iPair, nPairs                                    ! Loop indices and pair counts
  integer :: atm1, atm2                                          ! Atom indices for pair details
  logical :: isNeighbor                                          ! Neighbor status flag

  ! Step 1: Handle useDistStore case
  if (useDistStore) then
    if (present(NeighborDetailsNew)) then
      call NeighborUpdate_Distance_PairStore(nIndx, NeighborDetailsNew)
    else
      call NeighborUpdate_Distance_PairStore(nIndx)
    endif
    return
  endif

  ! Step 2: Update NeighborList for all molecules
  do j = 1, maxMol
    if (.not. isActive(j)) cycle
    if (j == nIndx) then
      NeighborList(nIndx, j) = .false.
      NeighborList(j, nIndx) = .false.
      if (present(NeighborDetailsNew)) then
        NeighborPairs(nIndx, j)%details%nPairs = 0
        NeighborPairs(nIndx, j)%details%pairIndices = 0
        NeighborPairs(j, nIndx)%details%nPairs = 0
        NeighborPairs(j, nIndx)%details%pairIndices = 0
      endif
      cycle
    endif

    isNeighbor = PairList(j) < Dist_Critr_sq
    NeighborList(nIndx, j) = isNeighbor
    NeighborList(j, nIndx) = isNeighbor

    ! Step 3: Update NeighborPairs if NeighborDetailsNew is present (minDistCriteria)
    if (present(NeighborDetailsNew) .and. isNeighbor) then
      nPairs = NeighborDetailsNew(j)%nPairs
      NeighborPairs(nIndx, j)%details%nPairs = nPairs
      NeighborPairs(j, nIndx)%details%nPairs = nPairs
      ! Update pairIndices for (nIndx,j)
      do iPair = 1, nPairs
        atm1 = NeighborDetailsNew(j)%pairIndices(2 * iPair - 1)
        atm2 = NeighborDetailsNew(j)%pairIndices(2 * iPair)
        NeighborPairs(nIndx, j)%details%pairIndices(2 * iPair - 1) = atm1
        NeighborPairs(nIndx, j)%details%pairIndices(2 * iPair) = atm2
        ! Mirror to (j,nIndx) with swapped atom indices
        NeighborPairs(j, nIndx)%details%pairIndices(2 * iPair - 1) = atm2
        NeighborPairs(j, nIndx)%details%pairIndices(2 * iPair) = atm1
      enddo
      ! Clear remaining indices to avoid stale data
      if (nPairs < maxAtoms) then
        NeighborPairs(nIndx, j)%details%pairIndices(2 * nPairs + 1:) = 0
        NeighborPairs(j, nIndx)%details%pairIndices(2 * nPairs + 1:) = 0
      endif
    elseif (present(NeighborDetailsNew)) then
      ! Clear NeighborPairs for non-neighbor pairs
      NeighborPairs(nIndx, j)%details%nPairs = 0
      NeighborPairs(nIndx, j)%details%pairIndices = 0
      NeighborPairs(j, nIndx)%details%nPairs = 0
      NeighborPairs(j, nIndx)%details%pairIndices = 0
    endif
  enddo
end subroutine NeighborUpdate_Distance
!=================================================================================     
!     This function updates the neighborlist if a move is accepted.
      subroutine LLNeighborUpdate_Distance(PairList, nIndx, NeighborDetailsNew)
      use Coords      
      use IndexingFunctions      
      use SimParameters
      implicit none     
      integer, intent(in) :: nIndx
      real(dp), intent(in) :: PairList(:)
      type(NeighborDetails), intent(in) :: NeighborDetailsNew(:)
      integer :: nType, nMol, jIndx, iPair

      nType = typeList(nIndx)
      nMol = subIndxList(nIndx)
      do jIndx=1,maxMol
        if(.not. isActive(jIndx)) then
          cycle
        endif
        if(jIndx .ne. nIndx) then  
          if( PairList(jIndx)  .lt.  Dist_Critr_sq ) then
            NeighborList(nIndx, jIndx) = .true.
            NeighborList(jIndx, nIndx) = .true.  
          else             
            NeighborList(nIndx, jIndx) = .false.
            NeighborList(jIndx, nIndx) = .false.            
          endif   
        endif
      enddo
!      NeighborList(nIndx,nIndx) = .false.
      do jIndx=1,maxMol
        if(.not. isActive(jIndx)) then
          cycle
        endif
        if(jIndx .ne. nIndx) then  
          NeighborPairs(nIndx, jIndx)%details%nPairs = NeighborDetailsNew(jIndx)%nPairs
          NeighborPairs(jIndx, nIndx)%details%nPairs = NeighborDetailsNew(jIndx)%nPairs
          do iPair = 1, NeighborDetailsNew(jIndx)%nPairs
            NeighborPairs(nIndx, jIndx)%details%pairIndices(2 * iPair - 1) = &
               NeighborDetailsNew(jIndx)%pairIndices(2 * iPair - 1)        
            NeighborPairs(nIndx, jIndx)%details%pairIndices(2 * iPair) = &
               NeighborDetailsNew(jIndx)%pairIndices(2 * iPair)
            NeighborPairs(jIndx, nIndx)%details%pairIndices(2 * iPair - 1) = &
               NeighborDetailsNew(jIndx)%pairIndices(2 * iPair)        
            NeighborPairs(jIndx, nIndx)%details%pairIndices(2 * iPair) = &
               NeighborDetailsNew(jIndx)%pairIndices(2 * iPair - 1)
          enddo   
        endif
      enddo



      end subroutine
!=================================================================================     
!     This function updates the neighborlist if a move is accepted.
      subroutine LLNeighborUpdate_SwapIn_Distance(PairList, nType, NeighborDetailsNew)
      use Coords      
      use IndexingFunctions      
      use SimParameters
      implicit none     
      integer, intent(in) :: nType
      real(dp), intent(in) :: PairList(:)
      type(NeighborDetails), intent(in) :: NeighborDetailsNew(:)
      integer ::  nIndx, nMol, jMol, jType, jIndx, iPair

      nMol = NPART(nType) + 1
      nIndx = MolArray(nType)%mol(nMol)%indx 
!      write(35,*) nType, nMol, nIndx , globIndx1
      do jType = 1, nMolTypes
        do jMol = 1, NPART(jType)
          jIndx = MolArray(jType)%mol(jMol)%indx
          if(jIndx .ne. nIndx) then  
            if( PairList(jIndx)  .lt.  Dist_Critr_sq ) then

              NeighborList(nIndx, jIndx) = .true.
              NeighborList(jIndx, nIndx) = .true.  
            else             
              NeighborList(nIndx, jIndx) = .false.
              NeighborList(jIndx, nIndx) = .false.            
            endif   
          endif
        enddo
      enddo

!      NeighborList(nIndx,nIndx) = .false.
      do jIndx=1,maxMol
        if(.not. isActive(jIndx)) then
          cycle
        endif
        if(jIndx .ne. nIndx) then  
          NeighborPairs(nIndx, jIndx)%details%nPairs = NeighborDetailsNew(jIndx)%nPairs
          NeighborPairs(jIndx, nIndx)%details%nPairs = NeighborDetailsNew(jIndx)%nPairs
          do iPair = 1, NeighborDetailsNew(jIndx)%nPairs
            NeighborPairs(nIndx, jIndx)%details%pairIndices(2 * iPair - 1) = &
               NeighborDetailsNew(jIndx)%pairIndices(2 * iPair - 1)        
            NeighborPairs(nIndx, jIndx)%details%pairIndices(2 * iPair) = &
               NeighborDetailsNew(jIndx)%pairIndices(2 * iPair)
            NeighborPairs(jIndx, nIndx)%details%pairIndices(2 * iPair - 1) = &
               NeighborDetailsNew(jIndx)%pairIndices(2 * iPair)        
            NeighborPairs(jIndx, nIndx)%details%pairIndices(2 * iPair) = &
               NeighborDetailsNew(jIndx)%pairIndices(2 * iPair - 1)
          enddo   
        endif
      enddo

      end subroutine
!=================================================================================  
! Updates the neighbor list and optionally neighbor pair details for a newly inserted molecule
! based on the distance cluster criterion (distCriteria or minDistCriteria) after an accepted
! insertion in a Grand Canonical Monte Carlo nucleation simulation. Sets NeighborList(nIndx, j)
! and NeighborList(j, nIndx) to .true. if the minimum squared distance (PairList(j)) is less
! than Dist_Critr_sq, otherwise sets to .false. If NeighborDetailsNew is present, updates
! NeighborPairs(i,j) and NeighborPairs(j,i) with atom pair indices in pre-allocated pairIndices
! arrays, maintaining symmetry. Handles useDistStore=.true. by calling
! NeighborUpdate_SwapIn_Distance_PairStore. Optimized for single-bead molecules (HPS force fields)
! but supports multi-atom molecules (LJ_Q, Mpipi). Called after a molecule insertion is accepted
! in Exchange or LLAVBMC moves.
subroutine NeighborUpdate_SwapIn_Distance(PairList, nType, NeighborDetailsNew)
  use SimParameters, only: maxMol, nMolTypes, NPART, Dist_Critr_sq, isActive, maxAtoms
  use DistanceCriteria_PairStore, only: NeighborUpdate_SwapIn_Distance_PairStore
  use VarPrecision, only: dp
  use CoordinateTypes, only: NeighborDetails
  use Coords, only: MolArray, NeighborPairs, NeighborList
  USE PairStorage, ONLY: useDistStore
  implicit none

  ! Input variables
  integer, intent(in) :: nType                                    ! Type index of the inserted molecule
  real(dp), intent(in) :: PairList(maxMol)                       ! Min squared distance between inserted molecule and others
  type(NeighborDetails), intent(in), optional :: NeighborDetailsNew(maxMol)  ! New atom pair details (for minDistCriteria)

  ! Local variables
  integer :: nMol, nIndx                                         ! Molecule instance and global index of inserted molecule
  integer :: jType, jMol, jIndx                                  ! Loop indices for molecules
  integer :: iPair, nPairs                                       ! Pair indices and counts for NeighborPairs
  integer :: atm1, atm2                                          ! Atom indices for pair details
  logical :: isNeighbor                                          ! Neighbor status flag

  ! Step 1: Handle useDistStore case
  if (useDistStore) then
    if (present(NeighborDetailsNew)) then
      call NeighborUpdate_SwapIn_Distance_PairStore(nType, NeighborDetailsNew)
    else
      call NeighborUpdate_SwapIn_Distance_PairStore(nType)
    endif
    return
  endif

  ! Step 2: Get index of the inserted molecule
  nMol = NPART(nType) + 1
  nIndx = MolArray(nType)%mol(nMol)%indx

  ! Step 3: Update NeighborList and NeighborPairs for all existing molecules
  do jType = 1, nMolTypes
    do jMol = 1, NPART(jType)
      jIndx = MolArray(jType)%mol(jMol)%indx
      if (.not. isActive(jIndx) .or. jIndx == nIndx) then
        NeighborList(nIndx, jIndx) = .false.
        NeighborList(jIndx, nIndx) = .false.
        if (present(NeighborDetailsNew)) then
          NeighborPairs(nIndx, jIndx)%details%nPairs = 0
          NeighborPairs(nIndx, jIndx)%details%pairIndices = 0
          NeighborPairs(jIndx, nIndx)%details%nPairs = 0
          NeighborPairs(jIndx, nIndx)%details%pairIndices = 0
        endif
        cycle
      endif

      isNeighbor = PairList(jIndx) < Dist_Critr_sq
      NeighborList(nIndx, jIndx) = isNeighbor
      NeighborList(jIndx, nIndx) = isNeighbor

      ! Update NeighborPairs if NeighborDetailsNew is present (minDistCriteria)
      if (present(NeighborDetailsNew)) then
        if (isNeighbor) then
          nPairs = NeighborDetailsNew(jIndx)%nPairs
          NeighborPairs(nIndx, jIndx)%details%nPairs = nPairs
          NeighborPairs(jIndx, nIndx)%details%nPairs = nPairs
          ! Update pairIndices for (nIndx,jIndx)
          do iPair = 1, nPairs
            atm1 = NeighborDetailsNew(jIndx)%pairIndices(2 * iPair - 1)
            atm2 = NeighborDetailsNew(jIndx)%pairIndices(2 * iPair)
            NeighborPairs(nIndx, jIndx)%details%pairIndices(2 * iPair - 1) = atm1
            NeighborPairs(nIndx, jIndx)%details%pairIndices(2 * iPair) = atm2
            ! Mirror to (jIndx,nIndx) with swapped atom indices
            NeighborPairs(jIndx, nIndx)%details%pairIndices(2 * iPair - 1) = atm2
            NeighborPairs(jIndx, nIndx)%details%pairIndices(2 * iPair) = atm1
          enddo
          ! Clear remaining indices to avoid stale data
          if (nPairs < maxAtoms) then
            NeighborPairs(nIndx, jIndx)%details%pairIndices(2 * nPairs + 1:) = 0
            NeighborPairs(jIndx, nIndx)%details%pairIndices(2 * nPairs + 1:) = 0
          endif
        else
          ! Clear NeighborPairs for non-neighbor pairs
          NeighborPairs(nIndx, jIndx)%details%nPairs = 0
          NeighborPairs(nIndx, jIndx)%details%pairIndices = 0
          NeighborPairs(jIndx, nIndx)%details%nPairs = 0
          NeighborPairs(jIndx, nIndx)%details%pairIndices = 0
        endif
      endif
    enddo
  enddo
  NeighborList(nIndx, nIndx) = .false.
end subroutine NeighborUpdate_SwapIn_Distance
!=================================================================================            
      end module
      
