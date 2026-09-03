!****************************************************************************************

!****************************************************************************************
      module E_Interface_Pedone
      use CoordinateTypes
!=============================================================================      
      ! interface
         ! subroutine TotalE(E_T)
           ! real(dp), intent(inout) :: E_T
         ! end subroutine

         ! subroutine ShiftE(E_Trial,disp)
           ! real(dp), intent(inout) :: E_Trial
           ! type(Displacement), intent(in) :: disp(:)
         ! end subroutine          
      ! end interface  
        
       
      ! procedure(TotalE), pointer :: TotalECalc => NULL()
      ! procedure(TotalE), pointer :: TotalECalc => NULL()      
!=============================================================================
      contains
!=============================================================================      
      subroutine Detailed_ECalc_Pedone(E_T,rejMove)
      use InterEnergy_Pedone
      use EnergyCriteria
      use DistanceCriteria      
      use SimParameters
      implicit none
      
      logical , intent(inout) :: rejMove
      real(dp), intent(inout) :: E_T
      integer :: i,j
      real(dp) :: PairList(1:maxMol,1:maxMol)
      
      E_T = 0d0
      call Detailed_ECalc_Inter(E_T,PairList)

      if(distCriteria) then
        call Detailed_DistanceCriteria(PairList,rejMove)
      else
        call Detailed_EnergyCriteria(PairList,rejMove)      
      endif
      
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
      subroutine Shift_ECalc_Pedone(E_Inter, E_Intra, disp, PairList, dETable, useIntra, rejMove, useInter, NeighborDetailsNew)
      use SimParameters, only: distCriteria, beta, softcutoff, NTotal
      use CBMC_Variables   
      use Coords      
      use DistanceCriteria            
      use EnergyCriteria
      use EnergyTables      
      use InterEnergy_Pedone
      implicit none
      
      logical, intent(in), optional :: useInter
      type(NeighVars), intent(inout), optional :: NeighborDetailsNew(:)
      logical, intent(in) :: useIntra(1:4)
      type(Displacement), intent(in) :: disp(:)
      logical, intent(inout) :: rejMove
      real(dp), intent(inout) :: dETable(:)
      real(dp), intent(out) :: E_Intra, E_Inter
      real(dp), intent(InOut) :: PairList(:)
      
      logical :: interSwitch
      integer :: nIndx, nDisp
      real(dp) :: E_NonBond, E_Stretch, E_Bend
      real(dp) :: E_Torsion, E_Improper

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

      if(present(useInter)) then
        interSwitch = useInter
      else
        interSwitch = .true.
      endif


!     Begin by calculating the intermolecular potential. If any atoms overlap the move will be rejected
!     immediately.      
      if(NTotal .gt. 1) then
          dETable = 0d0
          PairList = 0d0
          call Shift_ECalc_Inter(E_Inter,disp, PairList, dETable, rejMove)

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
        call Shift_DistanceCriteria(PairList, nIndx, rejMove)      
      else      
        call Shift_EnergyCriteria(PairList, nIndx, rejMove)
      endif        
      if(rejMove) then
        return
      endif


      
      end subroutine
!=============================================================================      
      subroutine SwapIn_ECalc_Pedone(E_Inter, E_Intra, PairList, dETable, rejMove, useInter, NeighborDetailsNew)
      use InterEnergy_Pedone
      use EnergyCriteria
      use EnergyTables        
 
      use Coords
      use CBMC_Variables      
      implicit none
      
      logical, intent(out) :: rejMove
      logical, intent(in), optional :: useInter
      type(NeighVars), intent(inout), optional :: NeighborDetailsNew(:)
      real(dp), intent(out) :: E_Inter, E_Intra
      real(dp), intent(inout) :: PairList(:), dETable(:)      

      logical :: interSwitch
      real(dp) :: E_NonBond, E_Stretch, E_Bend
      real(dp) :: E_Torsion, E_Improper

      
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

      call NewMol_ECalc_Inter(E_Inter, PairList, dETable, rejMove)
      if(rejMove) then
        return
      endif

      
      E_Inter_Diff = E_Inter 

      
      end subroutine
!=============================================================================
!     This function contains the energy calculations that are used when a molecule
!     has been selected for removal.  
      subroutine SwapOut_ECalc_Pedone(E_Inter, E_Intra, nType, nMol, dETable, useInter)
      use InterEnergy_Pedone
      use EnergyCriteria
      use EnergyTables         
      use Coords
      use CBMC_Variables
      implicit none
      
      logical, intent(in), optional :: useInter
      real(dp), intent(out) :: E_Inter, E_Intra      
      integer, intent(in) :: nType, nMol
      real(dp), intent(inout) :: dETable(:)      

      logical :: interSwitch      
      
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
      
