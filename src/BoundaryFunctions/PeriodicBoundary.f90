!=================================================================================
      module PeriodicBoundary
      contains
!=================================================================================     
!     Extensive Cluster Criteria Check.  Used at the start and end of the simulation. 
!     This ensures that all particles in the cluster are properly connected to each other.
!     This function also calculates the initial Neighborlist that is used throughout the simulation. 
      subroutine Detailed_PeriodicBoundary(rejMove)
      use SimParameters
      use Coords
      use IndexingFunctions
      use ParallelVar
      implicit none     
      logical,intent(out) :: rejMove
      real(dp),intent(inout) :: PairList(1:maxMol,1:maxMol)
      
      logical :: ClusterMember(1:maxMol)
      integer :: i,j,h,cnt
      integer :: iType,jType, iMol, jMol, iIndx, jIndx


     
      end subroutine
!=================================================================================     
!     This function determines if a given translational move will destroy a cluster. 
      subroutine Shift_PeriodicBoundary(nIndx, rejMove)
      use SimParameters     
      use Coords
      use IndexingFunctions
      implicit none     
      
      logical, intent(out) :: rejMove      
      real(dp), intent(in) :: PairList(1:maxMol)
      integer,intent(in) :: nIndx
      
     
      end subroutine
!=================================================================================     
!     This function determines if a given translational move will destroy a cluster. 
      pure subroutine SwapIn_PeriodicBoundary(nType, rejMove)
      use SimParameters     
      use Coords
      use IndexingFunctions
      implicit none     
      
      logical, intent(out) :: rejMove      
      real(dp), intent(in) :: PairList(1:maxMol)
      integer, intent(in) :: nType



      end subroutine
!=================================================================================     
!     This function determines if a given translational move will destroy a cluster. 
      subroutine SwapOut_PeriodicBoundary(nSwap, rejMove)
      use SimParameters     
      use Coords
      use IndexingFunctions
      implicit none     
      
      logical, intent(out) :: rejMove
      integer, intent(inout) :: nSwap



     
      end subroutine
!=================================================================================     
!     This function updates the neighborlist if a move is accepted.
      subroutine NeighborUpdate(nIndx)
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
!     This function updates the neighborlist if a move is accepted.
      subroutine NeighborUpdate_Delete(nIndx)
      use SimParameters
      use IndexingFunctions
      use Coords
      implicit none
      integer, intent(in) :: nIndx
      integer :: nType, nSwapIndx
      integer :: i

      nType = typeList(nIndx)     
      nSwapIndx = molArray(nType)%mol(NPART(nType))%indx

!      NeighborList(nIndx,:) = NeighborList(nSwapIndx,:)
!      NeighborList(:,nIndx) = NeighborList(:,nSwapIndx)

      if(nIndx .eq. nSwapIndx) then
        NeighborList(nSwapIndx,:) = .false.
        NeighborList(:,nSwapIndx) = .false.       
        return
      endif
     
      do i = 1, maxMol
        if(NeighborList(i,nSwapIndx)) then
          if(i .ne. nIndx) then
            NeighborList(i,nIndx) = .true.
            NeighborList(nIndx,i) = .true.            
          endif
        else
          NeighborList(i,nIndx) = .false.
          NeighborList(nIndx,i) = .false.
        endif        
      enddo
      
      NeighborList(nIndx,nIndx) = .false.

      NeighborList(nSwapIndx,:) = .false.
      NeighborList(:,nSwapIndx) = .false.      

      end subroutine      
!=================================================================================           
      subroutine MultipleSwap_PeriodicBoundary(nType2, nIndx1, PairList, isIncluded, rejMove)
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
      integer :: jlowerIndx      
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
      
