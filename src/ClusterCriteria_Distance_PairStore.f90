!=================================================================================
      module DistanceCriteria_PairStore

!      logical, allocatable :: ClusterMember(:), flipped(:) 
      contains
!=================================================================================     
!     Extensive Cluster Criteria Check.  Used at the start and end of the simulation. 
!     This ensures that all particles in the cluster are properly connected to each other.
!     This function also calculates the initial Neighborlist that is used throughout the simulation. 
      subroutine Detailed_DistanceCriteria(rejMove)
      use Coords
      use IndexingFunctions
      use ParallelVar
      use PairStorage
      use SimParameters
      implicit none     
      logical, intent(out) :: rejMove
!      real(dp), intent(inout) :: PairList(:, :)
      
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

      
!      do iType=1,nMolTypes
!      do jType=iType,nMolTypes
!        do iMol = 1 ,NPART(iType)      
!        iIndx = MolArray(iType)%mol(iMol)%indx
!        do jMol = 1 ,NPART(jType) 
!          jIndx = MolArray(jType)%mol(jMol)%indx        
!          if(PairList(iIndx,jIndx) .le. Dist_Critr_sq ) then
!            NeighborList(iIndx,jIndx)=.true.         
!            NeighborList(jIndx,iIndx)=.true.          
!          endif
!        enddo
!        enddo
!      enddo
!      enddo

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
              if(rPair(globIndx1, globIndx2) % p % r_sq .lt. Dist_Critr_sq ) then
!                write(*,*) globIndx1, globIndx2, rPair(globIndx1, globIndx2) % p % r_sq, rPair(globIndx1, globIndx2) % p % r
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
! Determines if a translational move destroys a single cluster in a grand canonical ensemble
! nucleation simulation when useDistStore=.true. Compares neighbors of the molecule’s old position
! (via NeighborList) and new position (via rPairNew for distCriteria or NeighborDetailsNew for
! minDistCriteria). Rejects the move if the new position’s cluster does not include all old neighbors.
! Supports distCriteria (distance between first atoms) and minDistCriteria (minimum distance between
! any atom pair). NeighborDetailsNew is required for minDistCriteria to provide new neighbor pairs.
! Compatible with all force fields (Mpipi, LJ_Q, HPS_single, HPS_piecewise, HPS_cation_pi).
subroutine Shift_DistanceCriteria(nIndx, rejMove, NeighborDetailsNew)
  use Coords, only: MolArray, NeighborList, typeList, subIndxList
  use PairStorage, only: rPairNew
  use SimParameters, only: maxMol, NTotal, isActive, Dist_Critr_sq, distCriteria, minDistCriteria
  use VarPrecision, only: dp
  use CoordinateTypes, only: NeighborDetails
  implicit none

  ! Input/Output variables
  integer, intent(in) :: nIndx                      ! Global index of moved molecule
  logical, intent(out) :: rejMove                   ! Flag to reject move if cluster breaks
  type(NeighborDetails), intent(in), optional :: NeighborDetailsNew(:) ! New neighbor pair details (minDistCriteria)

  ! Local variables
  logical :: neiFlipped                             ! True if all old neighbors are in new cluster
  logical :: memberAdded                            ! True if new members added to cluster
  logical :: ClusterMember(1:maxMol)                ! Tracks molecules in new cluster
  logical :: flipped(1:maxMol)                      ! Tracks processed molecules in propagation
  integer :: iIndx, jIndx, h, i                     ! Loop indices
  integer :: nType, nMol, jType, jMol, iMol, iType  ! Molecule type and instance indices
  integer :: globIndx1, globIndx2                   ! Global atom indices
  integer :: iAtom, jAtom                           ! Atom indices for neighbor pairs
  integer :: neiMax                                 ! Number of old position neighbors
  integer, allocatable :: curNeigh(:)               ! Indices of old position neighbors
  integer :: pairIdx, allocStat                     ! Neighbor pair index and allocation status

  ! Step 1: Initialize variables
  rejMove = .false.
  if (NTotal == 1) return ! Single molecule case is always valid
  ClusterMember = .false.
  flipped = .false.
  neiMax = 0

  ! Allocate curNeigh dynamically
  allocate(curNeigh(1:maxMol), stat=allocStat)
  if (allocStat /= 0) stop 'Error allocating curNeigh in Shift_DistanceCriteria'
  curNeigh = 0

  ! Get moved molecule info
  nType = typeList(nIndx)
  nMol = subIndxList(nIndx)
  globIndx1 = MolArray(nType)%mol(nMol)%globalIndx(1)

  ! Step 2: Identify neighbors of the new trial position
  memberAdded = .false.
  if (distCriteria) then
    ! Use rPairNew for first atom distances
    do jIndx = 1, maxMol
      if (.not. isActive(jIndx)) cycle
      if (nIndx == jIndx) cycle
      jType = typeList(jIndx)
      jMol = subIndxList(jIndx)
      globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
      if (rPairNew(globIndx1, globIndx2)%p%r_sq < Dist_Critr_sq) then
        ClusterMember(jIndx) = .true.
        memberAdded = .true.
      endif
    enddo
  elseif (minDistCriteria) then
    ! Use NeighborDetailsNew for minimum distance criterion
    if (.not. present(NeighborDetailsNew)) stop 'Error: NeighborDetailsNew required for minDistCriteria'
    do pairIdx = 1, size(NeighborDetailsNew)
      if (NeighborDetailsNew(pairIdx)%nPairs > 0) then
        do i = 1, NeighborDetailsNew(pairIdx)%nPairs
          iAtom = NeighborDetailsNew(pairIdx)%pairIndices(2 * i - 1)  ! First atom index
          jAtom = NeighborDetailsNew(pairIdx)%pairIndices(2 * i)      ! Second atom index
          ! Map atom indices to molecule indices
          do iIndx = 1, maxMol
            if (.not. isActive(iIndx)) cycle
            iType = typeList(iIndx)
            iMol = subIndxList(iIndx)
            if (MolArray(iType)%mol(iMol)%globalIndx(1) == iAtom .or. &
                MolArray(iType)%mol(iMol)%globalIndx(1) == jAtom) then
              if (iIndx == nIndx) then
                ! Find the other molecule in the pair
                do jIndx = 1, maxMol
                  if (.not. isActive(jIndx) .or. jIndx == iIndx) cycle
                  jType = typeList(jIndx)
                  jMol = subIndxList(jIndx)
                  if (MolArray(jType)%mol(jMol)%globalIndx(1) == iAtom .or. &
                      MolArray(jType)%mol(jMol)%globalIndx(1) == jAtom) then
                    ClusterMember(jIndx) = .true.
                    memberAdded = .true.
                  endif
                enddo
              else
                ClusterMember(iIndx) = .true.
                memberAdded = .true.
              endif
            endif
          enddo
        enddo
      endif
    enddo
  else
    stop 'Error: Neither distCriteria nor minDistCriteria specified'
  endif

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

  ! Step 4: Quick check if all old neighbors are in the new cluster
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
    if (neiFlipped .or. .not. memberAdded) exit
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
! distance cluster criterion (distCriteria or minDistCriteria) after an accepted move when
! useDistStore=.true. in a Grand Canonical Monte Carlo nucleation simulation. Sets NeighborList(nIndx, j)
! and NeighborList(j, nIndx) to .true. if the minimum squared distance between any atom pair of molecules
! nIndx and j in newDist is less than Dist_Critr_sq, otherwise sets to .false. If NeighborDetailsNew is
! present, updates NeighborPairs(i,j) and NeighborPairs(j,i) with atom pair indices in pre-allocated
! pairIndices arrays, maintaining symmetry. Optimized for single-bead molecules (HPS force fields) but
! supports multi-atom molecules (LJ_Q, Mpipi). Called by NeighborUpdate_Distance after a move is accepted
! in Translation or LLTranslation.
subroutine NeighborUpdate_Distance_PairStore(nIndx, NeighborDetailsNew)
  use SimParameters, only: maxMol, maxAtoms, Dist_Critr_sq, isActive
  use Coords, only: atomIndicies, NeighborList, NeighborPairs
  use PairStorage, only: newDist, nNewDist
  use VarPrecision, only: dp
  use CoordinateTypes, only: NeighborDetails
  implicit none

  ! Input variables
  integer, intent(in) :: nIndx                                    ! Index of the moved molecule
  type(NeighborDetails), intent(in), optional :: NeighborDetailsNew(maxMol)  ! New atom pair details (for minDistCriteria)

  ! Local variables
  integer :: j, iPair, nPairs                                    ! Loop indices and pair counts
  integer :: jMol, atm1, atm2                                    ! Molecule and atom indices
  logical :: isNeighbor(maxMol)                                  ! Neighbor status flags
  real(dp) :: minDistSq(maxMol)                                  ! Minimum squared distances

  ! Step 1: Initialize neighbor status and minimum distances
  isNeighbor = .false.
  minDistSq = huge(1.0_dp)
  NeighborList(nIndx, :) = .false.
  NeighborList(:, nIndx) = .false.

  ! Step 2: Find neighbors using newDist
  do iPair = 1, nNewDist
    ! Get molecule index for the second atom
    jMol = atomIndicies(newDist(iPair)%indx2)%nMol
    if (jMol == nIndx .or. .not. isActive(jMol)) cycle

    ! Update minimum squared distance
    if (newDist(iPair)%r_sq < minDistSq(jMol)) then
      minDistSq(jMol) = newDist(iPair)%r_sq
      isNeighbor(jMol) = newDist(iPair)%r_sq < Dist_Critr_sq
    endif
  enddo

  ! Step 3: Update NeighborList and NeighborPairs
  do j = 1, maxMol
    if (.not. isActive(j) .or. j == nIndx) then
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

    NeighborList(nIndx, j) = isNeighbor(j)
    NeighborList(j, nIndx) = isNeighbor(j)

    ! Update NeighborPairs if NeighborDetailsNew is present (minDistCriteria)
    if (present(NeighborDetailsNew)) then
      if (isNeighbor(j)) then
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
      else
        ! Clear NeighborPairs for non-neighbor pairs
        NeighborPairs(nIndx, j)%details%nPairs = 0
        NeighborPairs(nIndx, j)%details%pairIndices = 0
        NeighborPairs(j, nIndx)%details%nPairs = 0
        NeighborPairs(j, nIndx)%details%pairIndices = 0
      endif
    endif
  enddo
end subroutine NeighborUpdate_Distance_PairStore
!=================================================================================     
! Updates the neighbor list and optionally neighbor pair details for a newly inserted molecule
! based on the distance cluster criterion (distCriteria or minDistCriteria) after an accepted
! insertion when useDistStore=.true. in a Grand Canonical Monte Carlo nucleation simulation.
! Sets NeighborList(nIndx, j) and NeighborList(j, nIndx) to .true. if the minimum squared
! distance between any atom pair of molecules nIndx and j in newDist is less than Dist_Critr_sq,
! otherwise sets to .false. If NeighborDetailsNew is present, updates NeighborPairs(i,j) and
! NeighborPairs(j,i) with atom pair indices in pre-allocated pairIndices arrays, maintaining
! symmetry. Optimized for single-bead molecules (HPS force fields) but supports multi-atom
! molecules (LJ_Q, Mpipi). Called by NeighborUpdate_SwapIn_Distance after a molecule insertion
! is accepted in Exchange or LLAVBMC moves.
subroutine NeighborUpdate_SwapIn_Distance_PairStore(nType, NeighborDetailsNew)
  use SimParameters, only: maxMol, maxAtoms, nMolTypes, NPART, Dist_Critr_sq, isActive
  use Coords, only: MolArray, atomIndicies, NeighborList, NeighborPairs
  use PairStorage, only: newDist, nNewDist
  use VarPrecision, only: dp
  use CoordinateTypes, only: NeighborDetails
  implicit none

  ! Input variables
  integer, intent(in) :: nType                                    ! Type index of the inserted molecule
  type(NeighborDetails), intent(in), optional :: NeighborDetailsNew(maxMol)  ! New atom pair details (for minDistCriteria)

  ! Local variables
  integer :: nMol, nIndx                                         ! Molecule instance and global index of inserted molecule
  integer :: jType, jMol, jIndx, iPair                        ! Loop indices for molecules and pairs
  integer :: atm1, atm2                                          ! Atom indices for pair details
  integer :: nPairs                                              ! Pair counts for NeighborPairs
  logical :: isNeighbor(maxMol)                                  ! Neighbor status flags
  real(dp) :: minDistSq(maxMol)                                  ! Minimum squared distances

  ! Step 1: Get index of the inserted molecule
  nMol = NPART(nType) + 1
  nIndx = MolArray(nType)%mol(nMol)%indx

  ! Step 2: Initialize neighbor status and minimum distances
  isNeighbor = .false.
  minDistSq = huge(1.0_dp)
  NeighborList(nIndx, :) = .false.
  NeighborList(:, nIndx) = .false.

  ! Step 3: Find neighbors using newDist
  do iPair = 1, nNewDist
    ! Get molecule index for the second atom
    jMol = atomIndicies(newDist(iPair)%indx2)%nMol
    if (jMol == nIndx .or. .not. isActive(jMol)) cycle

    ! Update minimum squared distance
    if (newDist(iPair)%r_sq < minDistSq(jMol)) then
      minDistSq(jMol) = newDist(iPair)%r_sq
      isNeighbor(jMol) = newDist(iPair)%r_sq < Dist_Critr_sq
    endif
  enddo

  ! Step 4: Update NeighborList and NeighborPairs
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

      NeighborList(nIndx, jIndx) = isNeighbor(jIndx)
      NeighborList(jIndx, nIndx) = isNeighbor(jIndx)

      ! Update NeighborPairs if NeighborDetailsNew is present (minDistCriteria)
      if (present(NeighborDetailsNew)) then
        if (isNeighbor(jIndx)) then
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
end subroutine NeighborUpdate_SwapIn_Distance_PairStore
!=================================================================================           
      end module
      
