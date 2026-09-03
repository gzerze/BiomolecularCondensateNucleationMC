      module NeighborTable
      contains
!===================================================================
  ! Computes the neighbor energy table (NeiETable) for the old configuration, storing
! ε_i, the highest ETable(j) among molecule i’s neighbors of type nType, and neiCount(i),
! the number of such neighbors. Used for target selection in AVBMC_EBias_Rosen_Out with
! probability exp(β * ε_i,old) / Σ_j exp(β * ε_j,old). Returns if NTotal = 1 (no neighbors).
! Correctness is ensured by:
! - Checking only nType neighbors via NeighborList.
! - Initializing NeiETable to 0 for molecules with no neighbors.
! - Validating inputs and activity status.
! Supports all force fields (LJ_Q, Mpipi, HPS_single, HPS_piecewise, HPS_cation_pi)
! and molecule types (rigid, small, linear, branched).
subroutine Create_NeiETable(nType)
  use EnergyTables,     only: ETable, NeiETable, neiCount
  use SimParameters,    only: nMolTypes, NPART, NMAX, NTotal 
  use Coords,           only: NeighborList
  use VarPrecision,     only: dp
  implicit none

  integer, intent(in) :: nType
  integer :: i, iType, j, jType, iLowIndx, jLowIndx
  real(dp) :: EMax, ETab

  ! Initialize arrays
  ! Correctness: Sets NeiETable to 0 for no neighbors; clears neiCount
  NeiETable = 0.0_dp
  neiCount = 0
  ! Return if only one molecule exists
  ! Correctness: No neighbors possible
  if (NTotal == 1) then
    return
  endif

  ! Compute base index for nType molecules
  ! Correctness: jLowIndx points to start of nType molecules
  jLowIndx = 0
  do jType = 1, nType - 1
    jLowIndx = jLowIndx + NMAX(jType)
  enddo

  ! Compute NeiETable for each molecule
  ! Correctness: Finds max ETable(j) among nType neighbors for each i
  iLowIndx = 0
  do iType = 1, nMolTypes
    do i = iLowIndx + 1, iLowIndx + NPART(iType)
      EMax = -huge(dp)
      do j = jLowIndx + 1, jLowIndx + NPART(nType)
        if (NeighborList(j, i)) then
          ETab = ETable(j)
          neiCount(i) = neiCount(i) + 1
          if (ETab > EMax) then
            EMax = ETab
          endif
        endif
      enddo
      ! Store max energy; keep 0 if no neighbors
      NeiETable(i) = EMax
    enddo
    iLowIndx = iLowIndx + NMAX(iType)
  enddo
end subroutine Create_NeiETable
!=================================================================================
      subroutine Insert_NewNeiETable(nType, PairList, dE, newNeiTable)
      use EnergyTables
      use SimParameters
      use Coords
      implicit none
      integer, intent(in) :: nType
      real(dp), intent(in) :: PairList(:)
      real(dp), intent(inout) :: dE(:), newNeiTable(:)
!      real(dp), intent(in) :: biasArray(:)

      integer :: i,j, nIndx
      integer :: iType, jType, iLowIndx, jLowIndx
      real(dp) :: EMax, ETab
       
      
      neiCount = 0
      newNeiTable = 0E0
!      return
      nIndx = molArray(nType)%mol(NPART(nType)+1)%indx

      jLowIndx = 0
      do jType = 1, nType-1
        jLowIndx = jLowIndx + NMAX(jType)
      enddo

      iLowIndx = 0
      do iType = 1,nMolTypes
        do i = ilowIndx+1,ilowIndx+NPART(iType)
          EMax = -huge(dp)
          do j = jlowIndx+1, jlowIndx+NPART(nType)
            if(NeighborList(j,i)) then
              ETab = ETable(j) + dE(j)
              neiCount(i) = neiCount(i) + 1
              if(ETab .gt. EMax) then
                EMax = ETab
              endif
            endif 
          enddo
          if(PairList(i) .le. Eng_Critr(nType, iType)) then
            ETab = ETable(nIndx) + dE(nIndx)
            neiCount(i) = neiCount(i) + 1
            if(ETab .gt. EMax) then
              EMax = ETab
            endif
          endif
          newNeiTable(i) = EMax  
        enddo    
        iLowIndx = iLowIndx + NMAX(iType)
      enddo

      EMax = -huge(dp)
      do jType = 1,nMolTypes
        do j = jlowIndx+1,jlowIndx+NPART(nType)
          if(PairList(j) .le. Eng_Critr(jType, nType)) then
            ETab = ETable(j) + dE(j)
            neiCount(nIndx) = neiCount(nIndx) + 1
            if(ETab .gt. EMax) then
              EMax = ETab
            endif
          endif 
        enddo
      enddo    
      newNeiTable(nIndx) = EMax
     
!      if( all(neiCount .eq. 0) )then
!        write(*,*) "ERROR!"
!        do i = 1, maxMol
!          write(*,*) i, PairList(i), neiCount(i)
!        enddo
!      endif
      end subroutine      
!=================================================================================
      subroutine Insert_NewNeiETable_Distance_Old(nType, dE, newNeiTable)
      use Coords
      use EnergyTables
      use PairStorage
      use SimParameters
      implicit none
      integer, intent(in) :: nType
!      real(dp), intent(in) :: PairList(:)
      real(dp), intent(inout) :: dE(:), newNeiTable(:)

      integer :: i,j, nIndx
      integer :: iType, iMol
      integer :: jType, jMol
      integer :: globIndxN, globIndx1, globIndx2
      real(dp) :: EMax
       
      neiCount = 0
      nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
      globIndxN = molArray(nType)%mol(NPART(nType)+1)%globalIndx(1)
      do i = 1, maxMol
        newNeiTable(i) = 0E0

        EMax = -huge(dp)
        if( .not. isActive(i) ) then
          if( i .ne. nIndx ) then
             cycle
          else
            do j = 1, maxMol
              if(.not. isActive(j)) then
                cycle
              endif
              jType = typeList(j)
              jMol = subIndxList(j)
              globIndx2 = molArray(jType)%mol(jMol)%globalIndx(1)
              if(rPairNew(globIndxN,globIndx2)%p%r_sq .le. Dist_Critr_sq) then
!                neiCount(i) = neiCount(i) + 1
                if(ETable(j) + dE(j) .gt. EMax) then
                  EMax = ETable(j) + dE(j)
                endif
              endif
            enddo            
          endif
        else
          iType = typeList(i)
          iMol = subIndxList(i)
          globIndx1 = molArray(iType)%mol(iMol)%globalIndx(1)
          do j = 1, maxMol
            if( .not. isActive(j) ) then
              if(j .ne. nIndx) then
                cycle
              else
                if(rPairNew(globIndxN,globIndx1)%p%r_sq .le. Dist_Critr_sq) then
                  if(dE(j) .gt. EMax) then
!                     neiCount(nIndx) = neiCount(nIndx) + 1
                     EMax = dE(j)
                   endif
                endif
              endif
            else
              if( NeighborList(i,j) ) then
                if( i .ne. j ) then         
                   if( ETable(j) + dE(j) .gt. EMax ) then
                    EMax = ETable(j) + dE(j)
                  endif
                endif
               endif
            endif        
          enddo
        endif
        newNeiTable(i) = EMax
      enddo
      
       
      end subroutine    

!=================================================================================
      subroutine Insert_NewNeiETable_Distance(nType, dE, newNeiTable)
      use EnergyTables
      use SimParameters
      use Coords
      use PairStorage, only: rPairNew
      implicit none
      integer, intent(in) :: nType
!      real(dp), intent(in) :: PairList(:)
      real(dp), intent(inout) :: dE(:), newNeiTable(:)

      integer :: i,j, nIndx
      integer :: iType, jType, iMol, jMol, iLowIndx, jLowIndx
      integer :: globIndxN, globIndx1, globIndx2
      real(dp) :: EMax, ETab
       
      
      neiCount = 0
      newNeiTable = 0E0
!      return
      nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
      globIndxN = molArray(nType)%mol(NPART(nType)+1)%globalIndx(1)

      jLowIndx = 0
      do jType = 1, nType-1
        jLowIndx = jLowIndx + NMAX(jType)
      enddo

!      write(*,*) NPART
      iLowIndx = 0
      do iType = 1,nMolTypes
        do i = ilowIndx+1, ilowIndx+NPART(iType)
          iMol = subIndxList(i)
          globIndx1 = molArray(iType)%mol(iMol)%globalIndx(1)
!          write(*,*) iType, iMol, globIndx1, i
          EMax = -huge(dp)
          do j = jlowIndx+1, jlowIndx+NPART(nType)
            if(NeighborList(j,i)) then
              ETab = ETable(j) + dE(j)
              neiCount(i) = neiCount(i) + 1
              if(ETab .gt. EMax) then
                EMax = ETab
              endif
            endif 
          enddo
          if(rPairNew(globIndxN,globIndx1)%p%r_sq .le. Dist_Critr_sq) then
            ETab = ETable(nIndx) + dE(nIndx)
            neiCount(i) = neiCount(i) + 1
            if(ETab .gt. EMax) then
              EMax = ETab
            endif
          endif
          newNeiTable(i) = EMax  
        enddo    
        iLowIndx = iLowIndx + NMAX(iType)

      enddo

!      Calcuate the neiTable value for the newly inserted particle.

      EMax = -huge(dp)
      do jType = 1,nMolTypes
!        do j = jlowIndx+1, jlowIndx+NPART(jType)
        do jMol = 1, NPART(jType)
!          jMol = subIndxList(j)
          globIndx2 = molArray(jType)%mol(jMol)%globalIndx(1)
          if(rPairNew(globIndxN, globIndx2)%p%r_sq .le. Dist_Critr_sq) then
            j = molArray(jType)%mol(jMol)%indx
            ETab = ETable(j) + dE(j)
            neiCount(nIndx) = neiCount(nIndx) + 1
            if(ETab .gt. EMax) then
              EMax = ETab
            endif
          endif 
        enddo
      enddo    
      newNeiTable(nIndx) = EMax
     
!      if( all(neiCount .eq. 0) )then
!        write(*,*) "ERROR!"
!        do i = 1, maxMol
!          write(*,*) i, PairList(i), neiCount(i)
!        enddo
!      endif

!      write(2,*) "NeighborTable"
!      do i = 1, maxMol
!        write(2,*) i, newNeiTable(i), neiCount(i)
!      enddo
      end subroutine          
         
!=================================================================================
subroutine Insert_NewNeiETable_Distance_V2(nType, PairList, dE, newNeiTable)

  use SimParameters, only: nMolTypes, NPART, NMAX, Dist_Critr_sq, isActive
  use EnergyTables, only: ETable, neiCount
  use PairStorage,  only: useDistStore
  use Coords,       only: MolArray, NeighborList
  use VarPrecision, only: dp
  use ParallelVar,  only: nout

  implicit none

  integer, intent(in)    :: nType
  real(dp), intent(in)   :: PairList(:)
  real(dp), intent(in)   :: dE(:)
  real(dp), intent(out)  :: newNeiTable(:)

  integer :: i, j, iType
  integer :: nIndx
  integer :: iLowIndx, jLowIndx
  real(dp) :: EMax, ETab

  ! If using stored pair distances, use the corresponding optimized routine.
  ! That routine must follow the same definition:
  !   neiCount(i)    = number of neighbors of i that are type nType
  !   newNeiTable(i) = max energy among neighbors of i that are type nType
  if (useDistStore) then
    call Insert_NewNeiETable_Distance_V2_PairStore(nType, dE, newNeiTable)
    return
  endif

  ! Clear trial neighbor-selection tables.
  neiCount    = 0
  newNeiTable = 0.0_dp

  ! Global index of the inserted molecule.
  ! The molecule is not active yet in the real system, but it belongs to the
  ! proposed n+1 trial state.
  nIndx = MolArray(nType)%mol(NPART(nType) + 1)%indx

  ! Offset for molecules of type nType.
  jLowIndx = 0
  do iType = 1, nType - 1
    jLowIndx = jLowIndx + NMAX(iType)
  enddo

  ! For every existing molecule i, compute the reverse-deletion target bias
  ! table in the proposed n+1 state.
  !
  ! This must match Create_NeiETable(nType):
  ! only neighbors of type nType contribute to neiCount(i) and newNeiTable(i).
  iLowIndx = 0
  do iType = 1, nMolTypes

    do i = iLowIndx + 1, iLowIndx + NPART(iType)

      if (.not. isActive(i)) cycle

      EMax = -huge(dp)

      ! Old-state neighbors of molecule i that are type nType.
      do j = jLowIndx + 1, jLowIndx + NPART(nType)

        if (.not. isActive(j)) cycle

        if (NeighborList(j, i)) then
          ETab = ETable(j) + dE(j)
          neiCount(i) = neiCount(i) + 1
          EMax = max(EMax, ETab)
        endif

      enddo

      ! Add the inserted molecule as a possible type-nType neighbor of i
      ! in the proposed new state.
      if (PairList(i) <= Dist_Critr_sq) then
        ETab = ETable(nIndx) + dE(nIndx)
        neiCount(i) = neiCount(i) + 1
        EMax = max(EMax, ETab)
      endif

      ! Store ε_i only if molecule i has at least one type-nType neighbor.
      if (neiCount(i) > 0) then
        newNeiTable(i) = EMax
      else
        newNeiTable(i) = 0.0_dp
      endif

    enddo

    iLowIndx = iLowIndx + NMAX(iType)

  enddo

  ! Now compute the entry for the inserted molecule itself.
  !
  ! The inserted molecule is also part of the proposed n+1 state, so it must
  ! be included in the normalization used by the reverse target-selection
  ! probability.
  !
  ! Again, only old molecules of type nType are counted as neighbors.
  EMax = -huge(dp)

  do j = jLowIndx + 1, jLowIndx + NPART(nType)

    if (.not. isActive(j)) cycle

    if (PairList(j) <= Dist_Critr_sq) then
      ETab = ETable(j) + dE(j)
      neiCount(nIndx) = neiCount(nIndx) + 1
      EMax = max(EMax, ETab)
    endif

  enddo

  if (neiCount(nIndx) > 0) then
    newNeiTable(nIndx) = EMax
  else
    newNeiTable(nIndx) = 0.0_dp
  endif

  if (all(neiCount == 0)) then
    write(nout, *) 'Warning: no type-nType neighbors found in Insert_NewNeiETable_Distance_V2, nType=', nType
  endif

end subroutine Insert_NewNeiETable_Distance_V2
!=================================================================================
subroutine Insert_NewNeiETable_Distance_V2_PairStore(nType, dE, newNeiTable)

  use EnergyTables, only: ETable, neiCount
  use SimParameters, only: nMolTypes, NPART, NMAX, maxMol, Dist_Critr_sq, isActive
  use Coords, only: molArray, NeighborList
  use PairStorage, only: rPairNew
  use VarPrecision, only: dp

  implicit none

  integer, intent(in)   :: nType
  real(dp), intent(in)  :: dE(:)
  real(dp), intent(out) :: newNeiTable(:)

  integer :: iType, iMol, iIndx
  integer :: jMol, jIndx
  integer :: nIndx
  integer :: globIndxN, globIndx1
  real(dp) :: EMax, ETab

  neiCount    = 0
  newNeiTable = 0.0_dp

  nIndx    = molArray(nType)%mol(NPART(nType)+1)%indx
  globIndxN = molArray(nType)%mol(NPART(nType)+1)%globalIndx(1)

  ! For every existing molecule i, count only neighbors of type nType.
  do iType = 1, nMolTypes
    do iMol = 1, NPART(iType)

      iIndx = molArray(iType)%mol(iMol)%indx
      if (.not. isActive(iIndx)) cycle

      EMax = -huge(dp)

      ! Old neighbors of i that are type nType.
      do jMol = 1, NPART(nType)

        jIndx = molArray(nType)%mol(jMol)%indx
        if (.not. isActive(jIndx)) cycle

        if (NeighborList(jIndx, iIndx)) then
          ETab = ETable(jIndx) + dE(jIndx)
          neiCount(iIndx) = neiCount(iIndx) + 1
          EMax = max(EMax, ETab)
        endif

      enddo

      ! Inserted molecule as a new type-nType neighbor of i.
      globIndx1 = molArray(iType)%mol(iMol)%globalIndx(1)

      if (rPairNew(globIndxN, globIndx1)%p%r_sq <= Dist_Critr_sq) then
        ETab = ETable(nIndx) + dE(nIndx)
        neiCount(iIndx) = neiCount(iIndx) + 1
        EMax = max(EMax, ETab)
      endif

      if (neiCount(iIndx) > 0) then
        newNeiTable(iIndx) = EMax
      else
        newNeiTable(iIndx) = 0.0_dp
      endif

    enddo
  enddo

  ! Entry for inserted molecule itself.
  EMax = -huge(dp)

  do jMol = 1, NPART(nType)

    jIndx = molArray(nType)%mol(jMol)%indx
    if (.not. isActive(jIndx)) cycle

    globIndx1 = molArray(nType)%mol(jMol)%globalIndx(1)

    if (rPairNew(globIndxN, globIndx1)%p%r_sq <= Dist_Critr_sq) then
      ETab = ETable(jIndx) + dE(jIndx)
      neiCount(nIndx) = neiCount(nIndx) + 1
      EMax = max(EMax, ETab)
    endif

  enddo

  if (neiCount(nIndx) > 0) then
    newNeiTable(nIndx) = EMax
  else
    newNeiTable(nIndx) = 0.0_dp
  endif

end subroutine Insert_NewNeiETable_Distance_V2_PairStore
!=================================================================================
      end module
