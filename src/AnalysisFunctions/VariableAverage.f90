!======================================================================================    
!      This module contains the functions nessisary to initalize and 
!      perform single pair calculations.  This can be used
!      in the Umbrella Sampling algorithms or be used to study
!      di
      module SimpleDistPair
        use PairStorage
        use SimpleDistPair, only: pairArrayIndx
        use Q6Functions, only: q6ArrayIndx
        use Q6Functions, only: q6ArrayIndx        
        private

        integer :: nDistPair
        integer, allocatable :: variableIndx(:)

        public :: nDistPair, pairArrayIndx
        public :: Initialize_DistPair
        public :: SetPairVariables
        public :: CalcDistPairs
        public :: CalcDistPairs_New, CalcDistPairs_SwapIn, CalcDistPairs_SwapOut

        contains
     !--------------------------------------------------------------------------------
        subroutine Initialize_DistPair
           use MiscelaniousVars
           implicit none 
           integer :: AllocationStatus
           integer :: startIndx, endIndx, iPair


           if(nDistPair .eq. 0) then
             return
           endif 

           allocate(pairArrayIndx(1:nDistPair), STAT = AllocationStatus)
           allocate(pairGloIndx1(1:nDistPair), STAT = AllocationStatus)
           allocate(pairGloIndx2(1:nDistPair), STAT = AllocationStatus)
           allocate(molIndx1(1:nDistPair), STAT = AllocationStatus)
           allocate(molIndx2(1:nDistPair), STAT = AllocationStatus)
           
           call ReserveSpace_Coord(nDistPair, startIndx, endIndx)

           do iPair = 1, nDistPair
             pairArrayIndx(iPair) = startIndx + iPair - 1
           enddo


        end subroutine
     !--------------------------------------------------------------------------------
        subroutine SetPairVariables(nPair, nType1, nMol1, nAtom1, nType2, nMol2, nAtom2)
           use MiscelaniousVars
           use Coords, only: molArray
           implicit none 
           integer, intent(in) :: nPair 
           integer, intent(in) :: nType1, nMol1, nAtom1
           integer, intent(in) :: nType2, nMol2, nAtom2
           integer :: gloIndx1, gloIndx2, indx1, indx2

           gloIndx1 = molArray(nType1)%mol(nMol1)%globalIndx(nAtom1)
           gloIndx2 = molArray(nType2)%mol(nMol2)%globalIndx(nAtom2)
           indx1 = molArray(nType1)%mol(nMol1)%indx
           indx2 = molArray(nType2)%mol(nMol2)%indx

           pairGloIndx1(nPair) = gloIndx1
           pairGloIndx2(nPair) = gloIndx2
           molIndx1(nPair) = indx1
           molIndx2(nPair) = indx2
           
           rPair(gloIndx1, gloIndx2)%p%storeRValue = .true.

        end subroutine
     !--------------------------------------------------------------------------------
        subroutine CalcDistPairs
          use MiscelaniousVars
          use Coords
          implicit none 
          integer :: iDistPair
          integer :: gloIndx1, gloIndx2
!          real(dp) :: r, r_sq

          do iDistPair =1, nDistPair
            gloIndx1 = pairGloIndx1(iDistPair)
            gloIndx2 = pairGloIndx2(iDistPair)
!            r_sq = rPair(gloIndx1, gloIndx2)%p%r_sq           
!            r = dsqrt(r_sq)
!            r = rPair(gloIndx1, gloIndx2) % p % r
            miscCoord(pairArrayIndx(iDistPair)) = rPair(gloIndx1, gloIndx2) % p % r
!            write(*,*) miscCoord(pairArrayIndx(iDistPair))
           
          enddo


        end subroutine
     !--------------------------------------------------------------------------------
        subroutine CalcDistPairs_New(disp)
          use SimParameters, only: maxAtoms
          use MiscelaniousVars
          use Coords
          implicit none 
          type(Displacement), intent(in) :: disp(:)
          logical :: changed
          integer :: sizeDisp, iDistPair, iDisp
          integer :: nType, nMol, nAtom
          integer :: dispIndx, gloIndx1, gloIndx2
          integer :: gloList(1:maxAtoms)
          real(dp) :: rx, ry, rz, r

          sizeDisp = size(disp)
          do iDisp = 1, sizeDisp
            nType = disp(iDisp)%molType
            nMol = disp(iDisp)%molIndx
            nAtom = disp(iDisp)%atmIndx
            gloList(iDisp) = molArray(nType)%mol(nMol)%globalIndx(nAtom)
          enddo

          do iDistPair = 1, nDistPair
            gloIndx1 = pairGloIndx1(iDistPair)
            gloIndx2 = pairGloIndx2(iDistPair)
            changed = .false.
            do iDisp = 1, sizeDisp
              if( gloList(iDisp) .eq. gloIndx1 ) then
                nType = atomIndicies(gloIndx2)%nType
                nMol = atomIndicies(gloIndx2)%nMol
                nAtom = atomIndicies(gloIndx2)%nAtom
                dispIndx = iDisp
                changed = .true.
                exit
              elseif( gloList(iDisp) .eq. gloIndx2 ) then
                nType = atomIndicies(gloIndx1)%nType
                nMol = atomIndicies(gloIndx1)%nMol
                nAtom = atomIndicies(gloIndx1)%nAtom
                dispIndx = iDisp
                changed = .true.
                exit
              endif
            enddo
            if(.not. changed) then
              miscCoord_New(pairArrayIndx(iDistPair)) = miscCoord(pairArrayIndx(iDistPair))
              cycle
            endif
            rx = disp(dispIndx)%x_new - molArray(nType)%mol(nMol)%x(nAtom)
            ry = disp(dispIndx)%y_new - molArray(nType)%mol(nMol)%y(nAtom)
            rz = disp(dispIndx)%z_new - molArray(nType)%mol(nMol)%z(nAtom)
            r = rx*rx + ry*ry + rz*rz
            r = sqrt(r)
            miscCoord_New(pairArrayIndx(iDistPair)) = r       
          enddo


        end subroutine
     !--------------------------------------------------------------------------------
        subroutine CalcDistPairs_SwapIn
          use MiscelaniousVars
          implicit none 
          integer :: iDistPair
  
          do iDistPair = 1, nDistPair
            miscCoord_New(pairArrayIndx(iDistPair)) = miscCoord(pairArrayIndx(iDistPair)) 
          enddo


        end subroutine

     !--------------------------------------------------------------------------------
        subroutine CalcDistPairs_SwapOut(nType, nMol)
          use MiscelaniousVars
          implicit none 
          integer, intent(in) :: nType, nMol
          integer :: iDistPair
  
          do iDistPair = 1, nDistPair
            miscCoord_New(pairArrayIndx(iDistPair)) = miscCoord(pairArrayIndx(iDistPair)) 
          enddo


        end subroutine

     !--------------------------------------------------------------------------------
      end module
!====================================================================================== 


