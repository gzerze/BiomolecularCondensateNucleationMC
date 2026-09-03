!=====================================================================
!   This module contains data arrays whose purpose is to save the distance
!   pairs during the course of the simulation. The primary storage variable
!   is allocated as a 1D array to minimize the amount of memory that needs
!   to be used to describe each distance pair.  To facilitate easy access
!   a 2D array of pointers (rPair) is used so that the programer can quickly
!   access and modify the distance storage using matrix notation.

    module PairStorage
    use VarPrecision

     !Primary distance/energy pair storage variable defintinio
    type DistArray
      logical :: storeRValue
      logical :: storeRParts
      logical :: usePair
      integer :: arrayIndx
      integer :: indx1, indx2
      real(dp) :: rx, ry, rz
      real(dp) :: r_sq, r
      real(dp) :: E_Pair
    end type

     !Variable defintion for storing trial distances
    type DistArrayNew
      integer :: indx1, indx2
      real(dp) :: rx, ry, rz
      real(dp) :: r_sq, r
      real(dp) :: E_Pair
    end type

     !Variable used to create a 2D array of pointers to make
     !accessing the distance storage easier to work with
    type DistPointer
      type(DistArray), pointer :: p 
    end type

    type DistPointerNew
      type(DistArrayNew), pointer :: p 
    end type

    logical :: useDistStore = .true.
    logical :: coordShift = .false.
    integer :: nMaxPairs, nTotalAtoms, nNewDist
    integer, allocatable, target :: oldIndxArray(:)
    type(DistArray), allocatable, target :: distStorage(:)
    type(DistPointer), allocatable :: rPair(:,:)

    type(DistArrayNew), allocatable, target :: newDist(:)
    type(DistArrayNew), target :: nullPair
    type(DistPointerNew), allocatable :: rPairNew(:,:)

    contains
!=====================================================================
      subroutine CreateDistArrays
      implicit none
      integer :: AllocationStat
      integer :: i, j, cnt

      nMaxPairs = nint( dble(nTotalAtoms * (nTotalAtoms - 1)) / 2d0 )
      allocate(distStorage(0:nMaxPairs), stat = AllocationStat)
      allocate(rPair(1:nTotalAtoms, 1:nTotalAtoms), stat = AllocationStat)
      allocate(newDist(1:nMaxPairs), stat = AllocationStat) 
      allocate(rPairNew(1:nTotalAtoms, 1:nTotalAtoms), stat = AllocationStat)
      allocate(oldIndxArray(1:nMaxPairs), stat = AllocationStat) 


!     This block assigns all the relevent indicies and flags to the distanceStorage array.
      cnt = 0
      do i = 1, nTotalAtoms-1
        do j = i+1, nTotalAtoms
          cnt = cnt + 1
          distStorage(cnt)%arrayIndx = cnt
          distStorage(cnt)%indx1 = i
          distStorage(cnt)%indx2 = j
          distStorage(cnt)%r_sq = 0E0_dp
          distStorage(cnt)%r = 0E0_dp
          distStorage(cnt)%E_Pair = 0E0_dp
          distStorage(cnt)%usePair = .true.
          distStorage(cnt)%storeRValue = .false.
          distStorage(cnt)%storeRParts = .false. 
          rPair(i,j)%p => distStorage(cnt)
          rPair(j,i)%p => distStorage(cnt)
        enddo
      enddo

!     To avoid problems with unassocaited pointers, the diagonal of the pair array pointer
!     is set to reference the 0-th value of the array. 
      distStorage(0)%arrayIndx = 0
      distStorage(0)%indx1 = 0
      distStorage(0)%indx2 = 0
      distStorage(0)%r_sq = 0d0
      distStorage(0)%E_Pair = 0d0
      distStorage(0)%usePair = .false.
      distStorage(0)%storeRValue = .false.
      distStorage(0)%storeRParts = .false. 
      do i = 1, nTotalAtoms
        rPair(i,i)%p => distStorage(0)
      enddo 

      nullPair%indx1 = 0 
      nullPair%indx2 = 0
      nullPair%rx = 0E0_dp
      nullPair%ry = 0E0_dp
      nullPair%rz = 0E0_dp
      nullPair%r_sq = 0E0_dp
      nullPair%r = 0E0_dp
      nullPair%E_Pair = 0E0_dp


     end subroutine
!=====================================================================
     subroutine SetStorageFlags(q_tab)
     use Coords
     use ForceField
     use SimParameters, only: NMAX, nMolTypes
     implicit none
     real(dp), intent(in) :: q_tab(:,:)
     integer :: i, iType,jType,iMol,jMol,iAtom,jAtom
     integer(kind=atomIntType) :: atmType1, atmType2      
     integer :: globIndx1, globIndx2 
     real(dp) :: q_ij

     do i = 1, nMaxPairs
       distStorage(i) % storeRValue = .false.
     enddo 

     do iType = 1,nMolTypes
       do jType = iType, nMolTypes
         do iMol=1,NMAX(iType)
           do jMol = 1,NMAX(jType)
             do iAtom = 1,nAtoms(iType)
               atmType1 = atomArray(iType,iAtom)
               globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
               do jAtom = 1,nAtoms(jType)       
                 globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
                 atmType2 = atomArray(jType, jAtom)
                 q_ij = q_tab(atmType1, atmType2)
                 if(q_ij .ne. 0E0) then
                   rPair(globIndx1, globIndx2) % p % storeRValue = .true.
                 else
                   rPair(globIndx1, globIndx2) % p % storeRValue = .false.
                 endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

     end subroutine
!=====================================================================
     subroutine TurnOnAllStorageFlags
     implicit none
     integer :: i

     do i = 1, nMaxPairs
       distStorage(i) % storeRValue = .true.
       distStorage(i) % storeRParts = .true.
     enddo 



     end subroutine
!=====================================================================
     subroutine CalcAllDistPairs
      use Coords
      use ForceField
      use SimParameters, only: NPART, nMolTypes
      use ParallelVar,     only: nout
      implicit none
      integer :: iType,jType,iMol,jMol,iAtom,jAtom
      integer :: iIndx, jIndx
      integer(kind=atomIntType) :: atmType1, atmType2      
      integer :: globIndx1, globIndx2 
      real(dp) :: rx, ry, rz, r_sq
      real(dp) :: rmin_ij   

      do iType = 1,nMolTypes
        do jType = iType, nMolTypes
          do iMol=1,NPART(iType)
            iIndx = molArray(iType)%mol(iMol)%indx
            do jMol = 1,NPART(jType)
              jIndx = molArray(jType)%mol(jMol)%indx
              do iAtom = 1,nAtoms(iType)
                atmType1 = atomArray(iType,iAtom)
                globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
                do jAtom = 1,nAtoms(jType)       
                  globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
                  if(globIndx1 .eq. globIndx2) then
                    cycle
                  endif
!                  if(.not. rPair(globIndx1, globIndx2)%p%usePair) then
!                    cycle
!                  endif
                  atmType2 = atomArray(jType, jAtom) 
                  rmin_ij = r_min_tab(atmType2, atmType1)
                  rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
                  ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
                  rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom) 
!                  write(*,*) rPair(globIndx1, globIndx2)% p % indx1, rPair(globIndx1, globIndx2)% p % indx2
!                  write(*,*) globIndx1, globIndx2, rx, ry, rz
                  if( rPair(globIndx1, globIndx2) % p % storeRParts ) then
                    if(globIndx1 .gt. globIndx2) then
                      rPair(globIndx1, globIndx2)% p % rx = -rx
                      rPair(globIndx1, globIndx2)% p % ry = -ry
                      rPair(globIndx1, globIndx2)% p % rz = -rz
                    else
                      rPair(globIndx1, globIndx2)% p % rx = rx
                      rPair(globIndx1, globIndx2)% p % ry = ry
                      rPair(globIndx1, globIndx2)% p % rz = rz
                    endif
                  endif
                  r_sq = rx*rx + ry*ry + rz*rz
                  if(r_sq .lt. rmin_ij) then
                    if(iIndx .ne. jIndx) then
                      write(nout,*) "Index1:",iType, iMol, iAtom
                      write(nout,*) "Index2:",jType, jMol, jAtom
                      write(nout,*) "Distance:", sqrt(r_sq)
                      write(nout,*) "R_Min:", sqrt(rmin_ij)
                      stop "ERROR! Overlaping atoms found in the current configuration!"
                    endif
                  endif 
                  rPair(globIndx1, globIndx2) % p % r_sq = r_sq
                  rPair(globIndx1, globIndx2) % p % r = sqrt(r_sq)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo



     end subroutine
!=====================================================================
     subroutine CalcNewDistPairs(disp, rejMove)
      use Coords
      use ForceField
      use SimParameters, only: NPART, nMolTypes
      implicit none
      type(Displacement), intent(in) :: disp(:)
      logical, intent(out) :: rejMove
      integer :: iType,jType,iMol,jMol,iAtom,jAtom
      integer :: iIndx, jIndx
      integer(kind=atomIntType) :: atmType1,atmType2      
      integer :: iDisp, sizeDisp, globIndx1, globIndx2
      real(dp) :: rx,ry,rz,r_sq
      real(dp) :: rmin_ij   

   
      iType = disp(1)%molType

      iMol = disp(1)%molIndx
      iIndx = molArray(iType)%mol(iMol)%indx
      sizeDisp = size(disp)
      rejMove = .false.
      nNewDist = 0




      do iDisp = 1, sizeDisp
        iAtom = disp(iDisp)%atmIndx
        atmType1 = atomArray(iType, iAtom)
        globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
        rPairNew(globIndx1, globIndx1)%p => nullPair
        do jType = 1, nMolTypes
          do jAtom = 1, nAtoms(jType)        
            atmType2 = atomArray(jType,jAtom)
            rmin_ij = r_min_tab(atmType1, atmType2)     
            do jMol = 1, NPART(jType)
              jIndx = molArray(jType)%mol(jMol)%indx
              if(iIndx .eq. jIndx) then
                cycle
              endif
              globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
              if(globIndx1 .eq. globIndx2 ) then
                cycle
              endif
              if(.not. rPair(globIndx1, globIndx2)%p%usePair) then
                cycle
              endif   

!               Distance for the New position
              rx = disp(iDisp)%x_new - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = disp(iDisp)%y_new - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = disp(iDisp)%z_new - MolArray(jType)%mol(jMol)%z(jAtom)
              r_sq = rx*rx + ry*ry + rz*rz
!             If r_new is less than r_min reject the move.              
              if(r_sq .lt. rmin_ij) then
                rejMove = .true.
                return
              endif    
              nNewDist = nNewDist + 1
              rPairNew(globIndx1, globIndx2)%p => newDist(nNewDist)
              rPairNew(globIndx2, globIndx1)%p => newDist(nNewDist)
              oldIndxArray(nNewDist) = rPair(globIndx1, globIndx2)%p%arrayIndx
              newDist(nNewDist)%indx1 = globIndx1
              newDist(nNewDist)%indx2 = globIndx2
              newDist(nNewDist)%r_sq = r_sq
!              newDist(nNewDist)%E_Pair = 0d0
              if( rPair(globIndx1, globIndx2)%p%storeRValue ) then
                newDist(nNewDist)%r = sqrt(r_sq)
              endif
              if( rPair(globIndx1, globIndx2)%p%storeRParts ) then
                newDist(nNewDist)%rx = rx
                newDist(nNewDist)%ry = ry
                newDist(nNewDist)%rz = rz
              endif
            enddo
          enddo
        enddo
      enddo


     end subroutine
!=====================================================================
     subroutine CalcSwapInDistPairs(rejMove)
      use Coords
      use ForceField
      use SimParameters, only: NPART, nMolTypes
      implicit none
      logical, intent(out) :: rejMove
      integer :: iType, jType, iMol, jMol, iAtom, jAtom
      integer(kind=atomIntType) :: atmType1, atmType2      
      integer :: globIndx1, globIndx2
      real(dp) :: rx, ry, rz, r_sq
      real(dp) :: rmin_ij   

!      iType  = 1
!      do iMol = 1, NPART(iType)
!        write(*,*) "O", MolArray(iType)%mol(iMol)%x(1), MolArray(iType)%mol(iMol)%y(1), &
!                      MolArray(iType)%mol(iMol)%z(1)
!      enddo
!      write(*,*) "O", newMol%x(1), newMol%y(1), newMol%z(1)
!      write(*,*)   


      iType = newMol%molType
      iMol = NPART(iType)+1
      rejMove = .false.
      nNewDist = 0
      do iAtom = 1, nAtoms(iType)
        atmType1 = atomArray(iType, iAtom)
        globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
        rPairNew(globIndx1, globIndx1)%p => nullPair
        do jType = 1, nMolTypes
          do jAtom = 1, nAtoms(jType)        
            atmType2 = atomArray(jType,jAtom)
            rmin_ij = r_min_tab(atmType1, atmType2)
            do jMol = 1, NPART(jType)
              globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
!              if(.not. rPair(globIndx1, globIndx2)%p%usePair) then
!                cycle
!              endif   
              rx = newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r_sq = rx*rx + ry*ry + rz*rz
              if(r_sq .lt. rmin_ij) then
                rejMove = .true.
                return
              endif
              nNewDist = nNewDist + 1
              rPairNew(globIndx1, globIndx2)%p => newDist(nNewDist)
              rPairNew(globIndx2, globIndx1)%p => newDist(nNewDist)
              oldIndxArray(nNewDist) = rPair(globIndx1, globIndx2) % p % arrayIndx
              newDist(nNewDist)%indx1 = globIndx1
              newDist(nNewDist)%indx2 = globIndx2
              newDist(nNewDist)%r_sq = r_sq
!              newDist(nNewDist)%E_Pair = 0d0

              if( rPair(globIndx1, globIndx2)%p%storeRValue ) then
                newDist(nNewDist)%r = sqrt(r_sq)
              endif
              if( rPair(globIndx1, globIndx2)%p%storeRParts ) then
                newDist(nNewDist)%rx = rx
                newDist(nNewDist)%ry = ry
                newDist(nNewDist)%rz = rz
              endif

            enddo
          enddo
        enddo
      enddo


     end subroutine
!=====================================================================
     subroutine UpdateDistArray
      implicit none
      integer :: iPair
      integer :: oldIndx
      do iPair = 1, nNewDist
        oldIndx = oldIndxArray(iPair)
        distStorage(oldIndx)%r_sq = newDist(iPair)%r_sq
        distStorage(oldIndx)%E_Pair = newDist(iPair)%E_Pair

        if( distStorage(oldIndx)%storeRValue ) then
          distStorage(oldIndx)%r = newDist(iPair)%r
        endif

        if( distStorage(oldIndx)%storeRParts ) then
          if(newDist(iPair)%indx1 .eq. distStorage(oldIndxArray(iPair))%indx1) then
            distStorage(oldIndx)%rx = newDist(iPair)%rx
            distStorage(oldIndx)%ry = newDist(iPair)%ry
            distStorage(oldIndx)%rz = newDist(iPair)%rz
          else
            distStorage(oldIndx)%rx = -newDist(iPair)%rx
            distStorage(oldIndx)%ry = -newDist(iPair)%ry
            distStorage(oldIndx)%rz = -newDist(iPair)%rz
          endif
        endif
      enddo


     end subroutine
!=====================================================================
     subroutine UpdateDistArray_SwapOut(nType, nMol)
      use Coords
      use Forcefield, only: nAtoms
      use SimParameters, only: NPART
      implicit none
      integer, intent(in) :: nType, nMol
      integer :: iAtom, nMol2
      integer :: globIndx1, globIndx2, globIndx3


      nMol2 = NPART(nType)
      if(nMol .eq. nMol2) then
        return
      endif


      do iAtom = 1, nAtoms(nType)
        globIndx1 = molArray(nType)%mol(nMol)%globalIndx(iAtom)
        globIndx2 = molArray(nType)%mol(nMol2)%globalIndx(iAtom)
        do globIndx3 = 1, nTotalAtoms
          rPair(globIndx1, globIndx3)%p%r_sq = rPair(globIndx2, globIndx3)%p%r_sq
          rPair(globIndx1, globIndx3)%p%E_Pair = rPair(globIndx2, globIndx3)%p%E_Pair
          rPair(globIndx1, globIndx3)%p%storeRValue = rPair(globIndx2, globIndx3)%p%storeRValue
          if( rPair(globIndx2, globIndx3)%p%storeRValue ) then
            rPair(globIndx1, globIndx3)%p%r = rPair(globIndx2, globIndx3)%p%r
          endif
          if( rPair(globIndx1, globIndx3)%p%storeRParts ) then
            if(globIndx1 .gt. globIndx3) then
              rPair(globIndx1, globIndx3)%p%rx = rPair(globIndx2, globIndx3)%p%rx
              rPair(globIndx1, globIndx3)%p%ry = rPair(globIndx2, globIndx3)%p%ry
              rPair(globIndx1, globIndx3)%p%rz = rPair(globIndx2, globIndx3)%p%rz
            else
              rPair(globIndx1, globIndx3)%p%rx = -rPair(globIndx2, globIndx3)%p%rx
              rPair(globIndx1, globIndx3)%p%ry = -rPair(globIndx2, globIndx3)%p%ry
              rPair(globIndx1, globIndx3)%p%rz = -rPair(globIndx2, globIndx3)%p%rz
            endif
          endif
        enddo
      enddo


     end subroutine
!=====================================================================
     subroutine PrintDistArray
      integer :: iAtom, jAtom
      write(35,*) "----------------------------------------"
      write(35,*) "Distance List"
      do iAtom = 1, nTotalAtoms-1
        do jAtom = iAtom+1, nTotalAtoms
          write(35,*) iAtom, jAtom, rPair(iAtom, jAtom)%p%rx, rPair(iAtom, jAtom)%p%ry, &
                      rPair(iAtom, jAtom)%p%rz, rPair(iAtom, jAtom)%p%r
        enddo
      enddo


     end subroutine
!=====================================================================
!     subroutine CalculateAngle(globIndx1, globIndx2, globIndx3)
!      implicit none
!      integer, intent(in) :: globIndx1, globIndx2, globIndx3
!      real(dp) :: r12, rx12, ry12, rz12
!      real(dp) :: r23, rx23, ry23, rz23
!
!
!    end subroutine
!=====================================================================
    end module
!=====================================================================


