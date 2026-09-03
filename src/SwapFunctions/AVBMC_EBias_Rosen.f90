!****** Energy Biased Aggregation-Volume-Bias Monte Carlo (AVBMC) algorithm *******
!   This file contains the nessisary functions to impliment the energy biased swap
!   move for cluster simulations. 
!===================================================================================            
      module AVBMC_Module
      use VarPrecision
      real(dp), allocatable :: swapProb(:)
      contains
!===================================================================================            
      subroutine AVBMC(E_T, acc_x, atmp_x)
      use SimParameters
      use Constants       
      implicit none
      real(dp), intent(inout) :: E_T, atmp_x, acc_x
      real(dp) :: grnd


      prevMoveAccepted = .false.
!      call AVBMC_EBias_Rosen_In(E_T, maxMol, acc_x, atmp_x)     
!      call AVBMC_EBias_Rosen_Out(E_T, maxMol, acc_x, atmp_x) 
!      return
              
      if(grnd() .lt. 0.5d0) then
!        write(35,*) "in"
        call AVBMC_EBias_Rosen_In(E_T, maxMol, acc_x, atmp_x)     
      else
!        write(35,*) "out"
        call AVBMC_EBias_Rosen_Out(E_T, maxMol, acc_x, atmp_x)    
      endif

      end subroutine
!===================================================================================
      subroutine AVBMC_EBias_Rosen_In(E_T, arrayMax, acc_x, atmp_x)
      use AcceptRates
      use AVBMC_RejectionVar
      use AVBMC_CBMC
      use CBMC_Utility
      use CBMC_Variables
      use Constants  
      use Coords
      use SimParameters
      use ForceField
      use EnergyPressurePointers, only: SwapIn_ECalc, SwapIn_PCalc, Quick_Nei_ECalc
      use E_Interface_LJ_Q, only: Update_SubEnergies
      use UmbrellaFunctions
      use ForceField
      use IndexingFunctions      
      use EnergyCriteria
      use DistanceCriteria
!      use InterEnergy_LJ_Electro
      use EnergyTables
      use NeighborTable
      use SwapBoundary
!      use Pressure_LJ_Electro, only: NewMol_PressCalc_Inter
      use PairStorage, only: UpdateDistArray, PrintDistArray
      use UmbrellaSamplingNew, only: GetUmbrellaBias_SwapIn
      implicit none
      
      integer, intent(in) :: arrayMax
      real(dp), intent(inout) :: E_T      
      real(dp), intent(inout) :: acc_x, atmp_x
      logical :: rejMove
      logical :: isIncluded(1:arrayMax)
      integer :: NDiff(1:nMolTypes)
      integer :: i, nTargType, nTargMol, nTargIndx, nTarget
      integer :: nType, nIndx
      real(dp) :: grnd
      real(dp) :: genProbRatio, rosenRatio
      real(dp) :: E_Inter, E_Intra, biasDiff
      real(dp) :: PairList(1:arrayMax+1)
      real(dp) :: dETable(1:arrayMax)
      real(dp) :: newNeiETable(1:arrayMax)      
      real(dp) :: ProbTarg_In, ProbTarg_Out, ProbSel_Out
      real(dp) :: Boltzterm
      real(dp) :: sumInt, ranNum

      if(NTotal .eq. maxMol) then
        boundaryRej = boundaryRej + 1d0
        totalRej = totalRej + 1d0
        return
      endif

!      Choose the type of molecule to be inserted      
      if(nMolTypes .eq. 1) then
        nType = 1
      else
!        nType = floor(nMolTypes*grnd() + 1d0)
        ranNum = grnd()
        sumInt = swapProb(1)
        nType = 1
        do while(sumInt .lt. ranNum) 
          nType = nType + 1
          sumInt = sumInt + swapProb(nType)
        enddo
      endif

      NDiff = 0
      NDiff(nType) = 1
      rejMove = boundaryFunction(NPART, NDiff)
      if(rejMove) then
         boundaryRej = boundaryRej + 1d0
         totalRej = totalRej + 1d0
         return
      endif


      atmp_x = atmp_x + 1d0
      atmpSwapIn(nType) = atmpSwapIn(nType) + 1d0
      atmpInSize(NTotal) = atmpInSize(NTotal) + 1d0

      call EBias_Insert_ChooseTarget(nType, nTarget, nTargType, nTargMol, ProbTarg_In)
!      write(35,*) "Target:", nTarget
      nTargIndx = MolArray(nTargType)%mol(nTargMol)%indx      

!      Generate the configuration for the newly inserted molecule
      nIndx = MolArray(nType)%mol(NPART(nType) + 1)%indx      
      call Rosen_CreateSubset(nTarget, isIncluded)

      select case(regrowType(nType))
      case(0)
        call Ridgid_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      case(1)
        call Simple_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)   
      case(2)
        call StraightChain_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      case(3)
        call BranchedMol_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      case default
        write(*,*) "Error! EBias can not regrow a molecule of regrow type:", regrowType(nType)
        stop
      end select        
      if(rejMove) then
        totalRej = totalRej + 1d0
        ovrlapRej = ovrlapRej + 1d0
        return
      endif      

!      call DEBUG_Output_NewConfig


!      Perform a check to see if the cluster criteria is statisfied or not.
      if(.not. distCriteria) then
        rejMove = .false.
!        call QuickNei_ECalc_Inter(nTargType, nTargMol, rejMove)     
        call Quick_Nei_ECalc(nTargType, nTargMol, rejMove)     
        if(rejMove) then
          totalRej = totalRej + 1d0
          critriaRej = critriaRej + 1d0
          clusterCritRej = clusterCritRej + 1d0
          return
        endif  
      endif

!     Calculate the umbrella sampling bias.
      NPART_new = NPART + NDiff
      NTotal_New = NTotal + 1
      call GetUmbrellaBias_SwapIn(biasDiff, rejMove)
      if(rejMove) then
        boundaryRej = boundaryRej + 1d0
        totalRej = totalRej + 1d0
        return
      endif


!      Calculate the Energy Difference Associated with the move
      E_Inter = 0d0
      E_Intra = 0d0
!      call SwapIn_EnergyCalc(E_Inter, E_Intra, PairList, dETable, rejMove) 
      call SwapIn_ECalc(E_Inter, E_Intra, PairList, dETable, rejMove) 
      if(rejMove) then
        totalRej = totalRej + 1d0
        ovrlapRej = ovrlapRej + 1d0
        return
      endif
 

!     Determine the reverse probability of this move.
      if(distCriteria) then
!        call Insert_NewNeiETable_Distance(nType, dETable, newNeiETable)  
        call Insert_NewNeiETable_Distance_V2(nType, PairList, dETable, newNeiETable)  
      else
        call Insert_NewNeiETable(nType, PairList, dETable, newNeiETable)      
      endif
      call EBias_Insert_ReverseProbTarget(nTarget, nType, newNeiETable, ProbTarg_Out)
      call EBias_Insert_ReverseProbSel(nTarget, nType, dETable, ProbSel_Out)
      


!     Calculate acceptance probability and determine if the move is accepted or not          
      genProbRatio = (ProbTarg_Out * ProbSel_Out * avbmc_vol  * gas_dens(nType)) / (ProbTarg_In * rosenRatio)
      Boltzterm = exp(-beta*E_Inter + biasDiff)
      if( genProbRatio * Boltzterm .gt. grnd() ) then
         if(calcPressure) then
           call SwapIn_PCalc(P_Diff)
           pressure = pressure + P_Diff
         endif
!         call PrintDistArray
         acptSwapIn(nType) = acptSwapIn(nType) + 1d0        
         acptInSize(NTotal) = acptInSize(NTotal) + 1d0         
         do i = 1, nAtoms(nType)      
           molArray(nType)%mol(NPART(nType)+1)%x(i) = newMol%x(i)
           molArray(nType)%mol(NPART(nType)+1)%y(i) = newMol%y(i)
           molArray(nType)%mol(NPART(nType)+1)%z(i) = newMol%z(i)
         enddo
         E_T = E_T + E_Inter + E_Intra
         acc_x = acc_x + 1d0
!         nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
         isActive(nIndx) = .true.
         if(distCriteria) then
!           call NeighborUpdate_Distance(PairList,nIndx)  
!           call NeighborUpdate_Distance(nIndx)      
           call NeighborUpdate_SwapIn_Distance(PairList, nType)        
         else
           call NeighborUpdate(PairList, nIndx)
         endif  
         call UpdateDistArray
         NTotal = NTotal + 1
         ETable = ETable + dETable         
         NPART(nType) = NPART(nType) + 1 
         call Update_SubEnergies
         prevMoveAccepted = .true.
!         call PrintDistArray
        
       else
         totalRej = totalRej + 1d0
         dbalRej = dbalRej + 1d0
       endif
       end subroutine
!===================================================================================            
      subroutine AVBMC_EBias_Rosen_Out(E_T, arrayMax, acc_x, atmp_x)
      use AVBMC_CBMC
      use AVBMC_RejectionVar
      use SimParameters
      use Constants
!      use E_Interface
      use EnergyPressurePointers, only: SwapOut_ECalc, SwapOut_PCalc
      use E_Interface_LJ_Q, only: Update_SubEnergies
      use Coords
      use UmbrellaFunctions
      use ForceField
      use IndexingFunctions
      use EnergyCriteria
      use DistanceCriteria
      use EnergyTables
      use AcceptRates
      use UmbrellaFunctions
      use CBMC_Variables
      use NeighborTable
      use SwapBoundary
      use PairStorage, only: UpdateDistArray_SwapOut
      use UmbrellaSamplingNew, only: GetUmbrellaBias_SwapOut
!      use Pressure_LJ_Electro, only: Mol_PressCalc_Inter
      implicit none
      
      integer, intent(in) :: arrayMax
      real(dp), intent(inout) :: E_T      
      real(dp), intent(inout) :: acc_x, atmp_x
      
      logical :: rejMove  
      integer :: nTarget, nIndx
      integer :: nSel,nType, nMol,nTargMol,nTargType
      integer :: NDiff(1:nMolTypes)      
      real(dp) :: grnd
      real(dp) :: genProbRatio
      real(dp) :: biasDiff           
      real(dp) :: E_Inter, E_Intra
      real(dp) :: dETable(1:arrayMax)
      real(dp) :: ProbTargOut, ProbSel, ProbTargIn
      real(dp) :: rx, ry, rz, dist, rosenRatio
      real(dp) :: ranNum, sumInt


      if(NTotal .eq. 1) then
        boundaryRej_out = boundaryRej_out + 1d0
        totalRej_out = totalRej_out + 1d0
        return
      endif
      
!     Pick a Random Particle Type
      if(nMolTypes .eq. 1) then
        nType = 1
      else
!        nType = floor(nMolTypes*grnd() + 1d0)
        ranNum = grnd()
        sumInt = swapProb(1)
        nType = 1
        do while(sumInt .lt. ranNum) 
          nType = nType + 1
          sumInt = sumInt + swapProb(nType)
        enddo
      endif

      NDiff = 0
      NDiff(nType) = -1
      rejMove = boundaryFunction(NPART, NDiff)
      if(rejMove) then
         boundaryRej_out = boundaryRej_out + 1d0
         totalRej_out = totalRej_out + 1d0
         return
      endif   

      call Create_NeiETable(nType)
      call EBias_Remove_ChooseTarget(nType, nTarget, nTargType, nTargMol, ProbTargOut)
      call EBias_Remove_ChooseNeighbor(nTarget, nType, nSel, ProbSel)
!      nType = typeList(nSel)
      nMol = subIndxList(nSel)
      atmp_x = atmp_x + 1d0
      atmpSwapOut(nType) = atmpSwapOut(nType) + 1d0     

!     Check to see that the appropriate atoms are within the insertion distance
!     in order to ensure the move is reversible. If not reject the move since
!     the reverse probility is equal to 0. 
      if(.not. distCriteria) then
        rx = molArray(nTargType)%mol(nTargMol)%x(1) - molArray(nType)%mol(nMol)%x(1)
        ry = molArray(nTargType)%mol(nTargMol)%y(1) - molArray(nType)%mol(nMol)%y(1)
        rz = molArray(nTargType)%mol(nTargMol)%z(1) - molArray(nType)%mol(nMol)%z(1)
        dist = rx*rx + ry*ry + rz*rz
        if(dist .gt. Dist_Critr_sq) then
          critriaRej_out = critriaRej_out + 1d0
          totalRej_out = totalRej_out + 1d0
          return
        endif
      endif

!     Calculate the umbrella sampling bias.
!      NDiff = 0 
!      NDiff(nType) = -1
      NPART_new = NPART + NDiff
      NTotal_New = NTotal - 1
      call GetUmbrellaBias_SwapOut(nType, nMol, biasDiff, rejMove)
      if(rejMove) then
        boundaryRej_out = boundaryRej_out + 1d0
        totalRej_out = totalRej_out + 1d0
        return
      endif
        
!     Check to see if the deletion of the particle will break the cluster
      rejMove = .false.
      call SwapOut_ClusterCriteria(nSel, rejMove)
      if(rejMove) then
        critriaRej_out = critriaRej_out + 1d0
        totalRej_out = totalRej_out + 1d0
        return
      endif


      select case(regrowType(nType))
      case(0)
        call Ridgid_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      case(1)
        call Simple_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      case(2)
        call StraightChain_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      case(3)
        call BranchedMol_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      case default
        write(*,*) "Error! EBias can not regrow a molecule of regrow type:", nType
        stop
      end select 
!      Calculate the Energy Difference Associated with the move.
      E_Inter = 0d0
      E_Intra = 0d0
!      call SwapOut_EnergyCalc(E_Inter, E_Intra, nType, nMol, dETable)
      call SwapOut_ECalc(E_Inter, E_Intra, nType, nMol, dETable)
      call EBias_Remove_ReverseProbTarget(nTarget, nSel, nType, dETable, ProbTargIn)

!      genProbRatio = (ProbTargIn * rosenRatio) / (ProbTargOut * ProbSel * dble(nMolTypes) * avbmc_vol * gas_dens(nType))
      genProbRatio = (ProbTargIn * rosenRatio) / (ProbTargOut * ProbSel * avbmc_vol * gas_dens(nType))

!      Calculate Acceptance and determine if the move is accepted or not         
      if( genProbRatio * exp(-beta*E_Inter + biasDiff) .gt. grnd() ) then
         if(calcPressure) then
           call SwapOut_PCalc(nType, nMol, P_Diff)
           pressure = pressure - P_Diff
         endif
         acptSwapOut(nType) = acptSwapOut(nType) + 1d0
         molArray(nType)%mol(nMol)%x(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%x(1:nAtoms(nType))
         molArray(nType)%mol(nMol)%y(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%y(1:nAtoms(nType))
         molArray(nType)%mol(nMol)%z(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%z(1:nAtoms(nType))
         E_T = E_T + E_Inter + E_Intra
         nIndx = molArray(nType)%mol(nMol)%indx
         call NeighborUpdate_Delete(nIndx, molArray(nType)%mol(NPART(nType))%indx )
         call UpdateDistArray_SwapOut(nType, nMol)
         isActive( molArray(nType)%mol(NPART(nType))%indx ) = .false.
         ETable = ETable - dETable
         ETable(nIndx) = ETable( molArray(nType)%mol(NPART(nType))%indx )
         ETable( molArray(nType)%mol(NPART(nType))%indx ) = 0d0
         NPART(nType) = NPART(nType) - 1
         NTotal = NTotal - 1
         acc_x = acc_x + 1d0
         call Update_SubEnergies
         gasConfig(nType)%x(1:nAtoms(nType)) = newMol%x(1:nAtoms(nType))
         gasConfig(nType)%y(1:nAtoms(nType)) = newMol%y(1:nAtoms(nType))
         gasConfig(nType)%z(1:nAtoms(nType)) = newMol%z(1:nAtoms(nType))
         prevMoveAccepted = .true.
!         call DEBUG_Output_NeighborList
       else
         dbalRej_out = dbalRej_out + 1d0
         totalRej_out = totalRej_out + 1d0
       endif
       end subroutine
!=================================================================================    
      subroutine EBias_Insert_ChooseTarget(nInsType, nTarget, nTargType, nMol, ProbSel)
      use SimParameters  
      use Coords
      use EnergyTables
      implicit none
      integer , intent(in) :: nInsType
      integer, intent(out) :: nTarget, nTargType, nMol
      real(dp), intent(out) :: ProbSel
      
      integer :: i, iType
      integer :: cnt(1:nMolTypes)
      real(dp) :: avgE(1:nMolTypes)
      real(dp) :: ProbTable(1:maxMol)
      real(dp) :: grnd, norm       
      real(dp) :: ranNum, sumInt   
        
      ProbTable = 0d0
      avgE = 0d0
      cnt = 0 
      do i = 1, maxMol
        if(isActive(i)) then
           iType = typeList(i)
           avgE(iType) = avgE(iType) + ETable(i)
           cnt(iType) = cnt(iType) + 1
        endif
      enddo

      do iType = 1, nMolTypes
        if(cnt(iType) .ne. 0) then
          avgE(iType) = avgE(iType)/dble(cnt(iType))
        endif
      enddo
      
      do i = 1, maxMol
        if(isActive(i)) then
          iType = typeList(i)        
          ProbTable(i) = exp(AlphaTargetIns(nInsType,iType) * beta*(ETable(i)-avgE(iType)))
        endif
      enddo

      norm = sum(ProbTable)
      ranNum = norm*grnd()

      sumInt = 0d0
      nTarget = 0
      do while(sumInt .lt. ranNum)
        nTarget = nTarget + 1
        sumInt = sumInt + ProbTable(nTarget)
      enddo
      
      nTargType = typeList(nTarget) 
      ProbSel = ProbTable(nTarget)/norm
      
      nMol = 0
      do i = 1, nTargType-1
        nMol = nMol + NMAX(i)
      enddo
      nMol = nTarget - nMol
      
      end subroutine
!=================================================================================    
      subroutine EBias_Insert_ReverseProbTarget(nTarget, nType, newNeiETable, ProbRev)
      use SimParameters  
      use Coords
      use EnergyTables
      implicit none
      integer, intent(in) :: nTarget,nType
      real(dp), intent(in) :: newNeiETable(:)
      real(dp), intent(out) :: ProbRev
      
      integer :: i, nIndx, iType
      real(dp) :: ProbTable(1:maxMol)
      real(dp) :: norm, EMax
  
      


      nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
      ProbTable = 0E0_dp
      EMax = -huge(dp)
      do i = 1, maxMol
        if(neiCount(i) .gt. 0) then
          if(newNeiETable(i) .gt. EMax) then
            EMax = newNeiETable(i)
          endif
        endif
      enddo
      do i = 1, maxMol
        if(neiCount(i) .gt. 0) then
!          write(2,*) beta*(newNeiETable(i)-Emax)
          iType = typeList(i)
          if(AlphaTargetDel(nType,iType) * beta*(newNeiETable(i)-Emax) .gt. -30d0) then
!            write(2,*) beta*(newNeiETable(i)-Emax)
            ProbTable(i) = exp(AlphaTargetDel(nType,iType) * beta*(newNeiETable(i)-Emax))
          endif
        endif
      enddo
!      ProbTable(nIndx) = exp(beta*(newNeiETable(nIndx)-Emax))
!      write(*,*) ProbTable
      
      norm = sum(ProbTable)
      ProbRev = ProbTable(nTarget)/norm

      
      end subroutine
!=================================================================================    
      subroutine EBias_Insert_ReverseProbSel(nTarget, nType, dE, ProbRev)
      use SimParameters  
      use Coords
      use EnergyTables
      use UmbrellaFunctions
      implicit none
      integer, intent(in) :: nTarget, nType
      real(dp), intent(in) :: dE(:)
      real(dp), intent(out) :: ProbRev
      
      integer :: iMol, iIndx, nIndx
      real(dp) :: norm    
      
      nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
      norm = 0d0
      do iMol = 1, NPART(nType)
        iIndx = molArray(nType)%mol(iMol)%indx
        if( NeighborList(iIndx,nTarget) ) then
          norm = norm + exp(beta*(ETable(iIndx)+dE(iIndx)-dE(nIndx)))
        endif
      enddo
      norm = norm + 1E0
      ProbRev = 1E0/norm
      
      end subroutine
!=================================================================================    
      subroutine EBias_Remove_ChooseTarget(nType, nTarget, nTargType, nTargMol, ProbTarget)
      use Coords
      use EnergyTables
      use SimParameters 
      implicit none
      integer, intent(in) :: nType
      integer, intent(out) :: nTarget, nTargType, nTargMol
      real(dp), intent(out) :: ProbTarget
      
      integer :: i, iType
      real(dp) :: ProbTable(1:maxMol)
      real(dp) :: grnd, norm       
      real(dp) :: ranNum, sumInt   
        
      ProbTable = 0d0
      do i = 1, maxMol
!        if(isActive(i)) then
        if(neiCount(i) .gt. 0) then
          iType = typeList(i)
          ProbTable(i) = exp(AlphaTargetDel(nType,iType) * beta * NeiETable(i))
        endif
      enddo

      norm = sum(ProbTable)
      ranNum = norm * grnd()

      sumInt = ProbTable(1)
      nTarget = 1
      do while(sumInt .lt. ranNum)
        nTarget = nTarget + 1
        sumInt = sumInt + ProbTable(nTarget)
      enddo
      
      nTargType = typeList(nTarget) 
      ProbTarget = ProbTable(nTarget)/norm
      
      nTargMol = 0
      do i = 1, nTargType-1
        nTargMol = nTargMol + NMAX(i)
      enddo
      nTargMol = nTarget - nTargMol
      
      end subroutine
!=================================================================================    
      subroutine EBias_Remove_ChooseNeighbor(nTarget, nType, nSel, ProbTarget)
      use SimParameters  
      use Coords
      use EnergyTables
      implicit none
      integer, intent(in) :: nTarget, nType
      integer, intent(out) ::  nSel
      real(dp), intent(out) :: ProbTarget
      
      integer :: iIndx, iMol
      real(dp) :: ProbTable(1:maxMol)
      real(dp) :: grnd, norm       
      real(dp) :: ranNum, sumInt   

      
      ProbTable = 0d0
!      iIndx = molArray(nType)%mol(1)%indx
      do iMol = 1, NPART(nType)
        iIndx = molArray(nType)%mol(iMol)%indx
        if(NeighborList(iIndx,nTarget)) then
!          iType = typeList(i)
          ProbTable(iIndx) = exp(beta * ETable(iIndx))
!          ProbTable(i) = exp(beta * ETable(i))
        endif
!        iIndx = iIndx + 1
      enddo

      norm = sum(ProbTable)
      ranNum = norm * grnd()
      sumInt = ProbTable(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbTable(nSel)
      enddo
      
      ProbTarget = ProbTable(nSel) / norm
      
      end subroutine 
!=================================================================================    
      subroutine EBias_Remove_ReverseProbTarget(nTarget, nSel, nType, dE, ProbSel)
      use SimParameters  
      use Coords
      use EnergyTables
      implicit none
      integer, intent(in) :: nSel, nType
      integer, intent(in) :: nTarget
      real(dp), intent(in) :: dE(:)
      real(dp), intent(out) :: ProbSel

      
      integer :: i, iType
      integer :: cnt(1:nMolTypes)

      real(dp) :: avgE(1:nMolTypes)
      real(dp) :: ProbTable(1:maxMol)
      real(dp) :: norm       
        
      if(NTotal .eq. 2) then
         ProbSel = 1d0
         return      
      endif
      
      ProbTable = 0d0
      avgE = 0d0
      cnt = 0 
      do i = 1, maxMol
        if(isActive(i)) then
          if(i .ne. nSel) then        
            iType = typeList(i)
            avgE(iType) = avgE(iType) + ETable(i) - dE(i)
            cnt(iType) = cnt(iType) + 1

          endif
        endif
      enddo

      do iType = 1, nMolTypes
        if(cnt(iType) .ne. 0) then
          avgE(iType) = avgE(iType)/dble(cnt(iType))
        endif
      enddo
      
      do i = 1, maxMol
        if(isActive(i)) then
          if(i .ne. nSel) then    
            iType = typeList(i)        
            ProbTable(i) = exp(AlphaTargetIns(nType,iType) * beta*((ETable(i)-dE(i))-avgE(iType)))
          endif
        endif
      enddo
    
      norm = sum(ProbTable)
      ProbSel = ProbTable(nTarget)/norm

      end subroutine

!=================================================================================    
      end module
