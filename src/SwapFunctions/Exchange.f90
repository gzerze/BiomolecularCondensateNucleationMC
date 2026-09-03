!***********************************************************************************
!     This module contains experiemntal functions used to perform an exchange move.
!     In an exchange move one molecule is swapped out of the cluster at the same time
!     another molecule is swapped into the cluster. 
!===================================================================================            
      module Exchange_Module
      contains
!===================================================================================
      subroutine Exchange(E_T, acc_x, atmp_x)
      use AcceptRates
      use AVBMC_RejectionVar
      use AVBMC_CBMC
      use CBMC_Variables
      use CBMC_Utility
      use Constants  
      use Coords
      use ForceField
      use E_Interface_LJ_Q
      use IndexingFunctions      
      use EnergyCriteria
      use DistanceCriteria
      use InterEnergy_LJ_Electro
      use EnergyTables

      use SimParameters
      use NeighborTable
      use UmbrellaFunctions
      use EnergyPressurePointers, only: SwapOut_ECalc, SwapIn_ECalc
      implicit none
      
      real(dp), intent(inout) :: E_T      
      real(dp), intent(inout) :: acc_x, atmp_x
      logical rejMove    
      logical :: isIncluded(1:maxMol) 
      integer :: NDiff(1:nMolTypes)
      integer :: i, nTargType, nTargMol, nTargIndx, nTarget
      integer :: nMol1, nNei
      integer :: nType1, nType2, nIndx, nIndx2,  bIndx
      integer :: maxIndx
      real(dp) :: grnd
      real(dp) :: genProbRatio, rosenRatio_in, rosenRatio_out
      real(dp) :: E_Inter, E_Intra, bias_Diff
      real(dp) :: biasOld, biasNew
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)
      real(dp) :: dist, rx, ry, rz

      if(NTotal .eq. 1) then
        return
      endif
      atmp_x = atmp_x + 1d0
  
      nTarget = floor(NTotal*grnd() + 1d0)
      call Get_MolIndex(nTarget, NPart, nTargType, nTargMol)
      nTargIndx = MolArray(nTargType)%mol(nTargMol)%indx  
      call Uniform_ChooseNeighbor(nTargIndx, nIndx, nNei)
      nType1 = typeList(nIndx)
      nMol1 = subIndxList(nIndx)
      rejMove = .false.

      if(.not. distCriteria) then
        rx = molArray(nTargType)%mol(nTargMol)%x(1) - molArray(nType1)%mol(nMol1)%x(1)
        ry = molArray(nTargType)%mol(nTargMol)%y(1) - molArray(nType1)%mol(nMol1)%y(1)
        rz = molArray(nTargType)%mol(nTargMol)%z(1) - molArray(nType1)%mol(nMol1)%z(1)
        dist = rx*rx + ry*ry + rz*rz
        if(dist .gt. Dist_Critr_sq) then
          return
        endif
      endif

      nType2 = nType1
      do while(nType2 .eq. nType1) 
        nType2 = floor(nMolTypes*grnd() + 1d0)
      enddo

      
      if(NPART(nType2) .eq. NMAX(nType2)) then
        return
      endif

      NDiff = 0 
      NDiff(nType1) = -1
      NDiff(nType2) = 1

    

!      Generate the configuration for the newly inserted molecule
      nIndx2 = MolArray(nType2)%mol(NPART(nType2) + 1)%indx      
      call Rosen_CreateSubset(nTargIndx, isIncluded)
      isIncluded(nIndx) = .false.

      select case(regrowType(nType2))
      case(0)
        call Ridgid_RosenConfigGen(nType2, nIndx2, nTargIndx, nTargType, isIncluded, rosenRatio_in, rejMove)
      case(1)
        call Simple_RosenConfigGen(nType2, nIndx2, nTargIndx, nTargType, isIncluded, rosenRatio_in, rejMove)   
      case(2)
        call StraightChain_RosenConfigGen(nType2, nIndx2, nTargIndx, nTargType, isIncluded, rosenRatio_in, rejMove)   
      case default
        write(*,*) "Error! EBias can not regrow a molecule of regrow type:", regrowType(nType2)
        stop
      end select        
      if(rejMove) then
        return
      endif          

!      Perform a check to see if the cluster criteria is statisfied or not.
      if(.not. distCriteria) then
        rejMove = .false.
        call QuickNei_ECalc_Inter_LJ_Q(nTargType, nTargMol, rejMove)     
        if(rejMove) then
          return
        endif  
      endif
      
      select case(regrowType(nType1))
      case(0)
         call Ridgid_RosenConfigGen_Reverse(nType1, nMol1, nTargIndx, nTargType, rosenRatio_out)
      case(1)
         call Simple_RosenConfigGen_Reverse(nType1, nMol1, nTargIndx, nTargType, rosenRatio_out)
      case(2)
         call StraightChain_RosenConfigGen_Reverse(nType1, nMol1, nTargIndx, nTargType, rosenRatio_out)
      case default
         write(*,*) "Error! EBias can not regrow a molecule of regrow type:", regrowType(nType1)
         stop
      end select 


!      Calculate the Energy Difference Associated with the move
      E_Inter = 0d0
      E_Intra = 0d0
      call Exchange_ECalc_Inter(E_Inter, nType1, nMol1, PairList, dETable, rejMove)
      if(rejMove) then
        return
      endif
!      write(2,*) E_Inter
!      do i = 1, maxMol
!        if(PairList(i) .ne. 0d0) then
!          write(35,*) i, PairList(i)
!        endif
!      enddo
!      write(35,*)  
!      flush(35)

      if(.not. distCriteria) then
        rejMove = .false.
        isIncluded = isActive
        isIncluded(nIndx) = .false.
        call MultipleSwap_EnergyCriteria(nType2, nIndx, PairList, isIncluded, rejMove)    
        if(rejMove) then
          return
        endif  
      endif
      
!     Calculate the umbrella sampling bias.
      bIndx = getBiasIndex(NPart,NMAX)
      biasOld = NBias(bIndx)
      bIndx = getNewBiasIndex(NPart,NMAX, NDiff)
      biasNew = NBias(bIndx)      
      bias_diff = biasNew - biasOld

!     Calculate acceptance probability and determine if the move is accepted or not          
      genProbRatio = (gas_dens(nType2)*rosenRatio_out) / (gas_dens(nType1)*rosenRatio_in)
      if( genProbRatio * exp(-beta*E_Inter + bias_diff) .gt. grnd() ) then
         E_T = E_T + E_Inter
         E_Inter_T = E_Inter_T + E_Inter
         ETable = ETable + dETable  
         if(distCriteria) then
!           call NeighborUpdate_Distance(PairList, nIndx2)        
           call NeighborUpdate_Distance(PairList, nIndx2) 
         else
           call NeighborUpdate(PairList, nIndx2)
         endif 
         isActive(nIndx2) = .true.

         call SwapOut_ECalc(E_Inter, E_Intra, nType1, nMol1, dETable, .false.)
         call Update_SubEnergies
         E_T = E_T + E_Intra
         molArray(nType1)%mol(nMol1)%x(1:nAtoms(nType1)) = molArray(nType1)%mol(NPART(nType1))%x(1:nAtoms(nType1))
         molArray(nType1)%mol(nMol1)%y(1:nAtoms(nType1)) = molArray(nType1)%mol(NPART(nType1))%y(1:nAtoms(nType1))
         molArray(nType1)%mol(nMol1)%z(1:nAtoms(nType1)) = molArray(nType1)%mol(NPART(nType1))%z(1:nAtoms(nType1))
!         call NeighborUpdate_Delete(nIndx)
         maxIndx = molArray(nType1)%mol(NPART(nType1))%indx
         isActive(maxIndx) = .false.         
         ETable(nIndx) = ETable(maxIndx)
         ETable(maxIndx) = 0d0

         call SwapIn_ECalc(E_Inter, E_Intra, PairList, dETable, rejMove, .false.)
         call Update_SubEnergies
         E_T = E_T + E_Intra
         do i=1,nAtoms(nType2)      
           molArray(nType2)%mol(NPART(nType2)+1)%x(i) = newMol%x(i)
           molArray(nType2)%mol(NPART(nType2)+1)%y(i) = newMol%y(i)
           molArray(nType2)%mol(NPART(nType2)+1)%z(i) = newMol%z(i)
         enddo
         acc_x = acc_x + 1d0       
         NPART(nType1) = NPART(nType1) - 1 
         NPART(nType2) = NPART(nType2) + 1 
       endif

       end subroutine
!=================================================================================	  
      subroutine Uniform_ChooseNeighbor(nTarget, nSel, nNei)
      use SimParameters  
      use Coords
      implicit none
      integer, intent(in) :: nTarget
      integer, intent(out) :: nSel, nNei
      integer :: i,j
      integer :: ListCur(1:60)
      real(dp) :: grnd
 
       
      nNei = 0
      ListCur = 0
      do j = 1, maxMol
        if( NeighborList(nTarget,j) ) then
          if( j .ne. nTarget ) then             
            nNei = nNei + 1
            ListCur(nNei) = j
          endif
        endif
      enddo

      j = floor(nNei*grnd() + 1d0)
      nSel = ListCur(j)
      
      if(nSel .eq. 0) then
        write(35,*) "ERROR"
        write(35,*) "NTARGET:", nTarget
        write(35,*) "NPART:", NPART
        do i = 1, maxMol
          write(35,*) (NeighborList(i,j),j=1,maxMol)      
        enddo
      endif
      
      end subroutine
!=================================================================================    
      end module
