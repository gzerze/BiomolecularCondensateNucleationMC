!*********************************************************************************************************************
!     This file contains the energy functions that work for Pedone style forcefields
!     these functions are enclosed inside of the module "InterMolecularEnergy" so that
!     the energy functions can be freely exchanged from the simulation.
!     For the Pedone potential the oxygen atoms are always molecule type 1
!     The prefix naming scheme implies the following:
!           Detailed - Complete energy calculation inteded for use at the beginning and end
!                      of the simulation.  This function is not inteded for use mid-simulation.
!             Shift  - Calculates the energy difference for any move that does not result
!                      in molecules being added or removed from the cluster. This function
!                      receives any number of Displacement vectors from the parent function as input.
!              Mol   - Calculates the energy for a molecule already present in the system. For
!                      use in moves such as particle deletion moves. 
!             NewMol - Calculates the energy for a molecule that has been freshly inserted into the system.
!                      Intended for use in Rosenbluth Sampling, Swap In, etc.
!*********************************************************************************************************************
      module InterEnergy_Pedone
      use VarPrecision

      real(dp), parameter :: dieletric = 6d0

      contains
!======================================================================================      
      pure real(dp) function SolventFunction(r, q_ij, born1, born2) 
        real(dp), intent(in) :: r, q_ij, born1, born2
        real(dp) :: f

        f = sqrt(r*r + born1*born2*exp(-r*r/(4d0*born1*born2) ) ) 
        solventFunction = -0.5d0*(1d0-1d0/dieletric)*q_ij / f
      end function
!======================================================================================      
      subroutine Detailed_ECalc_Inter(E_T, PairList)
      use ParallelVar
      use ForceField
      use ForceFieldPara_Pedone
      use Coords
      use SimParameters
      use EnergyTables
      use PairStorage, only: rPair
      implicit none
      real(dp), intent(inOut) :: E_T
      real(dp), intent(inOut) :: PairList(:,:)
      integer :: iType,jType,iMol, jMol
      integer(kind=atomIntType) :: atmType1,atmType2      
      integer :: iIndx, jIndx, jMolMin
      integer :: gloIndx1, gloIndx2
      real(dp) :: r, r_sq
      real(dp) :: r_eq, repul_C, q_ij
      real(dp) :: alpha, delta
      real(dp) :: Morse, LJ, Ele, Solvent
      real(dp) :: E_Ele, E_LJ, E_Morse, E_Solvent
      real(dp) :: rmin_ij      
      real(dp) :: born1, born2

      E_LJ = 0d0
      E_Ele = 0d0
      E_Morse = 0d0
      E_Solvent = 0d0
      PairList = 0d0      
      ETable = 0d0

      Solvent = 0E0_dp

!      This first loop calculates oxide-oxide and oxide-metal interactions.
      do iType = 1, nMolTypes
        atmType1 = atomArray(iType, 1)
        do jType = iType, nMolTypes
          atmType2 = atomArray(jType, 1)
          r_eq = rEq_tab(atmType1, atmType2)
          q_ij = q_tab(atmType1, atmType2)
          alpha = alpha_Tab(atmType1, atmType2)
          delta = D_Tab(atmType1, atmType2)
          repul_C = repul_tab(atmType1, atmType2)
          rmin_ij = r_min_tab(atmType1, atmType2)
          do iMol=1,NPART(iType)
           if(iType .eq. jType) then
             jMolMin = iMol+1
           else
             jMolMin = 1
           endif
           iIndx = MolArray(iType)%mol(iMol)%indx
           gloIndx1 = MolArray(iType)%mol(iMol)%globalIndx(1)
           do jMol = jMolMin,NPART(jType)
             jIndx = MolArray(jType)%mol(jMol)%indx
             gloIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
             r_sq = rPair(gloIndx1, gloIndx2)%p%r_sq
             r = rPair(gloIndx1, gloIndx2)%p%r

             if(distCriteria) then
               PairList(iIndx, jIndx) = r_sq
               PairList(jIndx, iIndx) = PairList(iIndx,jIndx)
             endif
             LJ = 0E0_dp
             if(repul_C .ne. 0d0) then
               LJ = (1d0/r_sq)**6
               LJ = repul_C * LJ
               E_LJ = E_LJ + LJ
             endif

             Ele = q_ij/r
             if(implcSolvent) then
               born1 = bornRad(atmType1)
               born2 = bornRad(atmType2)
               Solvent = solventFunction(r, q_ij, born1, born2)
               E_Solvent = E_Solvent + Solvent
             endif
             E_Ele = E_Ele + Ele

             Morse = 0E0_dp
             if(delta .ne. 0E0_dp) then
               Morse = 1d0 - exp(-alpha*(r-r_eq))
               Morse = delta*(Morse*Morse - 1d0)
               E_Morse = E_Morse + Morse
             endif

             if(.not. distCriteria) then            
               PairList(iIndx, jIndx) = PairList(iIndx, jIndx) + Ele + LJ + Morse + Solvent
               PairList(jIndx, iIndx) = PairList(iIndx, jIndx)
             endif
             rPair(gloIndx1, gloIndx2)%p%E_Pair = + Ele + LJ + Morse + Solvent
             ETable(iIndx) = ETable(iIndx) + Ele + LJ + Morse + Solvent
             ETable(jIndx) = ETable(jIndx) + Ele + LJ + Morse + Solvent
            enddo
          enddo
        enddo
      enddo

!      if(implcSolvent) then
!        do iType = 1, nMolTypes
!          atmType1 = atomArray(iType, 1)
!          do iMol=1, NPART(iType)
!            q_ij = q_tab(atmType1, atmType1)
!            born1 = bornRad(atmType1)
!            E_Solvent = E_Solvent - 0.5d0*(1d0-1d0/dieletric)*q/(born1*born1)
!          enddo
!        enddo
!      endif
      write(nout,*) "Lennard-Jones Energy:", E_LJ
      write(nout,*) "Eletrostatic Energy:", E_Ele
      write(nout,*) "Morse Energy:", E_Morse
      if(implcSolvent) then
        write(nout,*) "Solvent Energy:", E_Solvent
      endif

!      write(35,*) "Pair List:"
!      do iMol=1,maxMol
!        write(35,*) iMol, PairList(iMol)
!      enddo
      
      E_T = E_T + E_Ele + E_LJ + E_Morse + E_Solvent
      E_Inter_T = E_Ele + E_LJ + E_Morse + E_Solvent
      
      end subroutine
!======================================================================================      
      subroutine Shift_ECalc_Inter(E_Trial, disp, newDist, PairList, dETable, rejMove)
      use ForceField
      use ForceFieldPara_Pedone
      use Coords
      use SimParameters
      use EnergyTables
      use PairStorage, only: distStorage, DistArrayNew, nNewDist, oldIndxArray
      implicit none
       
      type(Displacement), intent(in) :: disp(:)      
      type(DistArrayNew), intent(inout) :: newDist(:)
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      logical, intent(out) :: rejMove

      
      integer :: iType,jType,iMol,jMol, iPair
      integer(kind=atomIntType) :: atmType1,atmType2,iIndx,jIndx
      integer :: gloIndx1, gloIndx2
      real(dp) :: r, r_sq
      real(dp) :: r_eq, repul_C, q_ij
      real(dp) :: alpha, delta
      real(dp) :: Morse, LJ, Ele, Solvent
      real(dp) :: E_LJ, E_Ele, E_Morse, E_Solvent
      real(dp) :: E_Old, E_PairOld
      real(dp) :: born1, born2

      iType = disp(1)%molType
      iMol = disp(1)%molIndx
      atmType1 = atomArray(iType,1)
      iIndx = molArray(iType)%mol(iMol)%indx

      E_LJ = 0E0_dp
      E_Ele = 0E0_dp
      E_Morse = 0E0_dp
      E_Solvent = 0E0_dp
      E_Old = 0E0_dp
      Solvent = 0E0_dp
      rejMove = .false.

      do iPair = 1, nNewDist
        gloIndx1 = newDist(iPair)%indx1
        gloIndx2 = newDist(iPair)%indx2
        jType = atomIndicies(gloIndx2)%nType
        jMol  = atomIndicies(gloIndx2)%nMol
        atmType2 = atomArray(jType, 1)
        jIndx = MolArray(jType)%mol(jMol)%indx
        if(jIndx .ne. iIndx) then
          r_eq = rEq_tab(atmType1, atmType2)
          q_ij = q_tab(atmType1, atmType2)
          alpha = alpha_Tab(atmType1, atmType2)
          delta = D_Tab(atmType1, atmType2)
          repul_C = repul_tab(atmType1, atmType2)   
          r_sq = newDist(iPair)%r_sq
          r = newDist(iPair)%r

          if(distCriteria) then
            PairList(jIndx) = newDist(iPair)%r_sq
          endif
          LJ = 0E0_dp
          if(repul_C .ne. 0E0_dp) then
            LJ = (1E0_dp/r_sq)**6
            LJ = repul_C * LJ
            if(.not. distCriteria) then
              PairList(jIndx) = PairList(jIndx) + LJ
            endif
            E_LJ = E_LJ + LJ
          endif

          Ele = q_ij/r
          if(implcSolvent) then
            born1 = bornRad(atmType1)
            born2 = bornRad(atmType2)
            Solvent = solventFunction(r, q_ij, born1, born2)
            E_Solvent = E_Solvent + Solvent
          endif
          E_Ele = E_Ele + Ele
          if(.not. distCriteria) then
            PairList(jIndx) = PairList(jIndx) + Ele + Solvent
          endif

          Morse = 0E0_dp
          if(delta .ne. 0E0_dp) then
            Morse = 1E0_dp - exp(-alpha*(r-r_eq))
            Morse = delta*(Morse*Morse - 1E0_dp)
            if(.not. distCriteria) then            
              PairList(jIndx) = PairList(jIndx) + Morse 
            endif
            E_Morse = E_Morse + Morse
          endif
          newDist(iPair)%E_Pair = Ele + LJ + Solvent + Morse
          E_PairOld = distStorage(oldIndxArray(iPair))%E_Pair
          dETable(iIndx) = dETable(iIndx) + newDist(iPair)%E_Pair - E_PairOld
          dETable(jIndx) = dETable(jIndx) + newDist(iPair)%E_Pair - E_PairOld  
          E_Old = E_Old + E_PairOld
        endif
      enddo



      E_Trial = E_LJ + E_Ele + E_Morse - E_Old
      E_Inter_T =  E_LJ + E_Ele + E_Morse - E_Old    
      
      end subroutine

!======================================================================================      
      pure subroutine Mol_ECalc_Inter(iType, iMol, dETable, E_Trial)
      use ForceField
      use ForceFieldPara_Pedone
      use Coords
      use SimParameters
      use PairStorage
      implicit none
      integer, intent(in) :: iType, iMol     
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: dETable(:)
      
      integer :: iIndx,jType,jIndx,jMol
      integer :: gloIndx1, gloIndx2
      real(dp) :: E_Pair


      E_Trial = 0E0_dp
      dETable = 0E0_dp
      iIndx = MolArray(iType)%mol(iMol)%indx
      gloIndx1 = MolArray(iType)%mol(iMol)%globalIndx(1) 
      do jType = 1, nMolTypes
        do jMol=1, NPART(jType)
          jIndx = MolArray(jType)%mol(jMol)%indx  
          if(iIndx .ne. jIndx) then
            gloIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1) 
            E_Pair = rPair(gloIndx1, gloIndx2)%p%E_Pair
            E_Trial = E_Trial + E_Pair
            dETable(iIndx) = dETable(iIndx) + E_Pair
            dETable(jIndx) = dETable(jIndx) + E_Pair
          endif
        enddo
      enddo
      
      end subroutine

!======================================================================================      
      subroutine NewMol_ECalc_Inter(E_Trial, PairList, dETable, rejMove)
      use ForceField
      use ForceFieldPara_Pedone
      use Coords
      use SimParameters
      use PairStorage
      implicit none
      logical, intent(out) :: rejMove
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      
      integer :: iType, iMol, iIndx, jType, jIndx, jMol, iPair
      integer(kind=atomIntType) :: atmType1,atmType2
      integer :: gloIndx1, gloIndx2
      real(dp) :: r, r_sq      
      real(dp) :: r_eq, repul_C, q_ij
      real(dp) :: alpha, delta
      real(dp) :: LJ, Ele, Morse
      real(dp) :: E_Ele,E_LJ, E_Morse, E_Solvent
      real(dp) :: born1, born2

      E_LJ = 0E0_dp
      E_Ele = 0E0_dp    
      E_Morse = 0E0_dp
      E_Trial = 0E0_dp
      E_Solvent = 0E0_dp
      PairList = 0E0_dp
      dETable = 0E0_dp
      rejMove = .false.


      iType = newMol%molType
      iMol = NPART(iType)+1
      iIndx = molArray(iType)%mol(iMol)%indx
      do iPair = 1, nNewDist

        gloIndx1 = newDist(iPair)%indx1
        gloIndx2 = newDist(iPair)%indx2

        iType = atomIndicies(gloIndx1)%nType
        iMol = atomIndicies(gloIndx1)%nMol
        iIndx = molArray(iType)%mol(iMol)%indx

        jType = atomIndicies(gloIndx2)%nType
        jMol = atomIndicies(gloIndx2)%nMol
        jIndx = MolArray(jType)%mol(jMol)%indx
        if(iIndx .ne. jIndx) then
          atmType1 = atomArray(iType, 1)
          atmType2 = atomArray(jType, 1)
          r_eq = rEq_tab(atmType1, atmType2)
          q_ij = q_tab(atmType1, atmType2)
          alpha = alpha_Tab(atmType1, atmType2)
          delta = D_Tab(atmType1, atmType2)
          repul_C = repul_tab(atmType1, atmType2)   
          r_sq = newDist(iPair)%r_sq
          r = newDist(iPair)%r

          LJ = 0E0_dp
          if(repul_C .ne. 0E0_dp) then
            LJ = (1E0_dp/r_sq)**6
            LJ = repul_C * LJ
            E_LJ = E_LJ + LJ
          endif

          Ele = q_ij/r
          if(implcSolvent) then
            born1 = bornRad(atmType1)
            born2 = bornRad(atmType2)
            Ele = Ele + solventFunction(r, q_ij, born1, born2)
          endif
          E_Ele = E_Ele + Ele

          Morse = 0E0_dp
          if(delta .ne. 0E0_dp) then
            Morse = 1d0 - exp(-alpha*(r-r_eq))
            Morse = delta*(Morse*Morse - 1d0)
            E_Morse = E_Morse + Morse
          endif
          if(.not. distCriteria) then                
            PairList(jIndx) = PairList(jIndx) + LJ + Ele + Morse
          endif
          dETable(iIndx) = dETable(iIndx) + LJ + Ele + Morse
          dETable(jIndx) = dETable(jIndx) + LJ + Ele + Morse
          newDist(iPair)%E_Pair = LJ + Ele + Morse
        endif 
         
      enddo

!      if(implcSolvent) then
!        q_ij = q_tab(atmType1, atmType1)
!        born1 = bornRad(atmType1)
!        E_Solvent = -(1d0-1d0/dieletric)*q/(born1*born1)
!      endif
     
     
      E_Trial = E_LJ + E_Ele + E_Morse + E_Solvent
      
      
      end subroutine    
!======================================================================================      
      subroutine QuickNei_ECalc_Inter_Pedone(jType, jMol, rejMove)
      use ForceField
      use ForceFieldPara_Pedone
      use Coords
      use SimParameters
      implicit none
      integer, intent(in) :: jType, jMol     
      logical, intent(out) :: rejMove

      integer(kind=atomIntType)  :: atmType1,atmType2

      real(dp) :: rx,ry,rz,r
      real(dp) :: r_eq, repul_C, q_ij
      real(dp) :: alpha, delta
      real(dp) :: LJ, Ele, Morse
      real(dp) :: E_Trial, E_Ele, E_LJ, E_Morse
      real(dp) :: rmin_ij
      real(dp) :: born1, born2

      E_LJ = 0d0
      E_Ele = 0d0      
      E_Morse = 0d0
      E_Trial = 0d0
      rejMove = .false.
    
      atmType1 = atomArray(newMol%molType, 1)
      atmType2 = atomArray(jType, 1)
      r_eq = rEq_tab(atmType1, atmType2)
      q_ij = q_tab(atmType1, atmType2)
      alpha = alpha_Tab(atmType1, atmType2)
      delta = D_Tab(atmType1, atmType2)
      repul_C = repul_tab(atmType1, atmType2)          
      rmin_ij = r_min_tab(atmType1, atmType2) 

      rx = newMol%x(1) - MolArray(jType)%mol(jMol)%x(1)
      ry = newMol%y(1) - MolArray(jType)%mol(jMol)%y(1)
      rz = newMol%z(1) - MolArray(jType)%mol(jMol)%z(1)

      r = rx*rx + ry*ry + rz*rz
      if(r .lt. rmin_ij) then
         rejMove = .true.
         return
      endif          

      if(repul_C .ne. 0d0) then
        LJ = ( 1d0/r )**6
        LJ = repul_C * LJ
        E_LJ = E_LJ + LJ
      endif

      r = sqrt(r)
      Ele = q_ij/r
      if(implcSolvent) then
        born1 = bornRad(atmType1)
        born2 = bornRad(atmType2)
        Ele = Ele + solventFunction(r, q_ij, born1, born2)
      endif
      E_Ele = E_Ele + Ele

      if(delta .ne. 0d0) then
        Morse = 1d0-exp(-alpha*(r-r_eq))
        Morse = delta*(Morse*Morse - 1d0)
        E_Morse = E_Morse + Morse
      endif
     
      E_Trial = E_LJ + E_Ele + E_Morse

      if( E_Trial .gt. Eng_Critr(newMol%molType,jType) ) then
        rejMove = .true.
      endif
!      write(2,*) "E:",E_Trial , rejMove
      
      end subroutine
!======================================================================================
      end module
