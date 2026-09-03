!****************************************************************************************

!****************************************************************************************
      module E_Interface_Tersoff
      use CoordinateTypes
!=============================================================================
      contains
!=============================================================================      
      subroutine Detailed_ECalc_Tersoff(E_T, rejMove)
      use InterEnergy_Tersoff
      use EnergyCriteria
      use DistanceCriteria_PairStore    
      use SimParameters
      use PairStorage, only: CalcAllDistPairs, SetStorageFlags, PrintDistArray
      implicit none
      
      logical , intent(inout) :: rejMove
      real(dp), intent(inout) :: E_T
      integer :: i,j
      real(dp) :: PairList(1:maxMol, 1:maxMol)

      E_T = 0d0
      call CalcAllDistPairs
      call Detailed_ECalc_Inter(E_T, PairList)

      if(distCriteria) then
!        call Detailed_DistanceCriteria(PairList, rejMove)
        call Detailed_DistanceCriteria(rejMove)
      else
        call Detailed_EnergyCriteria(PairList, rejMove)      
      endif

      write(35,*) "Cluster Criteria Satisifed?: ", .not. rejMove
      write(35,*) "Pairlist:"
      do i = 1, maxMol
        if(isActive(i)) then
          do j = 1, maxMol
            if(isActive(j)) then
              write(35,*) i, j, PairList(i,j)            
            endif
          enddo
        endif
      enddo
      write(35,*)

      end subroutine
!=============================================================================      
!     This function contains the energy and cluster criteria functions for any move
!     where a molecule is moved within a cluster.  This function takes a set of displacement
!     vectors (the "disp" variable) and returns the change in energy for both the Intra- and
!     Inter-molecular components.  This function can be used for for moves that any number of
!     atoms in a given molecule, but it can not be used if more than one molecule changes
!     in a given move.  
      subroutine Shift_ECalc_Tersoff(E_Inter, E_Intra, disp, PairList, dETable, useIntra, rejMove, &
                                          useInter, NeighborDetailsNew)
      use SimParameters, only: distCriteria, beta, softcutoff, NTotal
      use CBMC_Variables   
      use Coords      
      use DistanceCriteria_PairStore      
      use EnergyCriteria
      use EnergyTables      
      use InterEnergy_Tersoff
      use PairStorage
      implicit none
      
      logical, intent(in), optional :: useInter
      type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)
      logical, intent(in) :: useIntra(1:4)
      type(Displacement), intent(in) :: disp(:)
      logical, intent(inout) :: rejMove
      real(dp), intent(inout) :: dETable(:)
      real(dp), intent(out) :: E_Intra, E_Inter
      real(dp), intent(InOut) :: PairList(:)
      
      logical :: interSwitch
      integer :: nIndx, nDisp


      nDisp = size(disp)
      rejMove = .false.
      E_Inter = 0d0
      E_Intra = 0d0
      dETable = 0d0

      E_Inter_Diff = 0d0
      E_NBond_Diff = 0d0
      E_Strch_Diff = 0d0
      E_Bend_Diff = 0d0
      E_Tors_Diff = 0d0



!     Begin by calculating the intermolecular potential. If any atoms overlap the move will be rejected
!     immediately.      
      if(NTotal .gt. 1) then
          call CalcNewDistPairs(disp, rejMove)
          if(rejMove) then
            return
          endif
          dETable = 0d0
          PairList = 0d0
          call Shift_ECalc_Inter(E_Inter, disp,newDist, PairList, dETable)

          if(rejMove) then
            return
          endif
          E_Inter_Diff = E_Inter

          if(E_Inter*beta .gt. softCutOff) then
            rejMove = .true.
            return
          endif
      endif

!        Using the data collected from the intermolecular function, check to see that the new position
!        satisfies the cluster criteria.      
      nIndx = MolArray( disp(1)%molType )%mol( disp(1)%molIndx )%indx
      if(distCriteria) then
!        call Shift_DistanceCriteria(PairList, nIndx, rejMove)    
        call Shift_DistanceCriteria(nIndx, rejMove)   
      else      
        call Shift_EnergyCriteria(PairList, nIndx, rejMove)
      endif        
      if(rejMove) then
        return
      endif


      
      end subroutine
!=============================================================================      
      subroutine SwapIn_ECalc_Tersoff(E_Inter, E_Intra, PairList, dETable, rejMove, useInter, NeighborDetailsNew)
      use Coords
      use CBMC_Variables      
      use EnergyCriteria
      use EnergyTables        
      use InterEnergy_Tersoff
      use PairStorage
      implicit none
      
      logical, intent(out) :: rejMove
      logical, intent(in), optional :: useInter
      type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)
      real(dp), intent(out) :: E_Inter, E_Intra
      real(dp), intent(inout) :: PairList(:), dETable(:)      

      
      rejMove = .false.      
      E_Inter = 0d0
      E_Intra = 0d0
      PairList = 0d0
      dETable = 0d0

      E_Inter_Diff = 0d0
      E_NBond_Diff = 0d0
      E_Strch_Diff = 0d0
      E_Bend_Diff = 0d0
      E_Tors_Diff = 0d0

      call CalcSwapInDistPairs(rejMove)
      if(rejMove) then
        return
      endif
      call NewMol_ECalc_Inter(E_Inter, PairList, dETable)
      if(rejMove) then
        return
      endif

      
      E_Inter_Diff = E_Inter 

      
      end subroutine
!=============================================================================
!     This function contains the energy calculations that are used when a molecule
!     has been selected for removal.  
      subroutine SwapOut_ECalc_Tersoff(E_Inter, E_Intra, nType, nMol, dETable, useInter)
      use InterEnergy_Tersoff
      use EnergyCriteria
      use EnergyTables         
      use Coords
      use CBMC_Variables
      implicit none
      
      logical, intent(in), optional :: useInter
      real(dp), intent(out) :: E_Inter, E_Intra      
      integer, intent(in) :: nType, nMol
      real(dp), intent(inout) :: dETable(:)         
      
      E_Inter = 0d0
      E_Intra = 0d0
         

      E_Inter_Diff = 0d0
      E_NBond_Diff = 0d0
      E_Strch_Diff = 0d0
      E_Bend_Diff = 0d0
      E_Tors_Diff = 0d0
      dETable = 0d0

      call Mol_ECalc_Inter(nType, nMol, dETable, E_Inter)
      E_Inter = -E_Inter
      E_Inter_Diff = E_Inter
  
      end subroutine      

!=============================================================================      
      
      end module
      
