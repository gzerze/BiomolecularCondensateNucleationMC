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
      implicit none
      real(dp), intent(inOut) :: E_T
      real(dp), intent(inOut) :: PairList(:,:)
      integer :: iType,jType,iMol, jMol
      integer(kind=atomIntType) :: atmType1,atmType2      
      integer :: iIndx, jIndx, jMolMin
      real(dp) :: rx,ry,rz,r
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
           do jMol = jMolMin,NPART(jType)
             jIndx = MolArray(jType)%mol(jMol)%indx  

             rx = MolArray(iType)%mol(iMol)%x(1) - MolArray(jType)%mol(jMol)%x(1)
             ry = MolArray(iType)%mol(iMol)%y(1) - MolArray(jType)%mol(jMol)%y(1)
             rz = MolArray(iType)%mol(iMol)%z(1) - MolArray(jType)%mol(jMol)%z(1) 
             r = rx**2 + ry**2 + rz**2
             if(r .lt. rmin_ij) then
               stop "ERROR: Overlaping atoms found in the configuration!"
             endif 
             if(distCriteria) then
               PairList(iIndx, jIndx) = r
               PairList(jIndx, iIndx) = PairList(iIndx,jIndx)                    
             endif
             if(repul_C .ne. 0d0) then
               LJ = (1d0/r)**6
               LJ = repul_C * LJ
               E_LJ = E_LJ + LJ
             endif

             r = sqrt(r)
             Ele = q_ij/r
             if(implcSolvent) then
               born1 = bornRad(atmType1)
               born2 = bornRad(atmType2)
               Solvent = solventFunction(r, q_ij, born1, born2)
               E_Solvent = E_Solvent + Solvent
             endif
             E_Ele = E_Ele + Ele

             if(delta .ne. 0d0) then
               Morse = 1d0 - exp(-alpha*(r-r_eq))
               Morse = delta*(Morse*Morse - 1d0)
               E_Morse = E_Morse + Morse
             endif
             if(.not. distCriteria) then            
               PairList(iIndx, jIndx) = PairList(iIndx, jIndx) + Ele + LJ + Morse + Solvent
               PairList(jIndx, iIndx) = PairList(iIndx, jIndx)
             endif
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
      write(nout,*) "Solvent Energy:", E_Solvent

!      write(35,*) "Pair List:"
!      do iMol=1,maxMol
!        write(35,*) iMol, PairList(iMol)
!      enddo
      
      E_T = E_T + E_Ele + E_LJ + E_Morse + E_Solvent
      E_Inter_T = E_Ele + E_LJ + E_Morse + E_Solvent
      
      end subroutine
!======================================================================================      
      pure subroutine Shift_ECalc_Inter(E_Trial, disp, PairList, dETable, rejMove)
      use ForceField
      use ForceFieldPara_Pedone
      use Coords
      use SimParameters
      implicit none
      
      type(Displacement), intent(in) :: disp(:)      
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      logical, intent(out) :: rejMove

      
      integer :: iType,jType,iMol,jMol,iDisp
      integer(kind=atomIntType) :: atmType1,atmType2,iIndx,jIndx
      integer :: sizeDisp 
      real(dp) :: rx,ry,rz
      real(dp) :: r_new, r_old
      real(dp) :: r_min1_sq      
      real(dp) :: r_eq, repul_C, q_ij
      real(dp) :: alpha, delta
      real(dp) :: Morse, LJ, Ele
      real(dp) :: E_Ele,E_LJ, E_Morse
      real(dp) :: rmin_ij    
      real(dp) :: born1, born2

      sizeDisp = size(disp)
      E_LJ = 0d0
      E_Ele = 0d0      
      E_Morse = 0d0
      E_Trial = 0d0
      PairList = 0d0      

      dETable = 0d0
      iType = disp(1)%molType
      iMol = disp(1)%molIndx
      iIndx = MolArray(iType)%mol(iMol)%indx
      atmType1 = atomArray(iType, 1)
      
!      !This section calculates the Intermolecular interaction between the atoms that
!      !have been modified in this trial move with the atoms that have remained stationary
       do iDisp = 1, sizeDisp
         do jType = 1, nMolTypes
           atmType2 = atomArray(jType, 1)
           r_eq = rEq_tab(atmType1, atmType2)
           q_ij = q_tab(atmType1, atmType2)
           alpha = alpha_Tab(atmType1, atmType2)
           delta = D_Tab(atmType1, atmType2)
           repul_C = repul_tab(atmType1, atmType2)          
           rmin_ij = r_min_tab(atmType1, atmType2)  
           do jMol = 1, NPART(jType)
             if(iType .eq. jType) then
               if(iMol .eq. jMol) then
                 cycle
               endif
             endif  
             
!            Distance for the New position
             rx = disp(iDisp)%x_new - MolArray(jType)%mol(jMol)%x(1)
             ry = disp(iDisp)%y_new - MolArray(jType)%mol(jMol)%y(1)
             rz = disp(iDisp)%z_new - MolArray(jType)%mol(jMol)%z(1)
             r_new = rx*rx + ry*ry + rz*rz

!            If r_new is less than r_min reject the move.              
             if(r_new .lt. rmin_ij) then
               rejMove = .true.
               return
             endif        

             jIndx = MolArray(jType)%mol(jMol)%indx   
             if(distCriteria) then   
               PairList(jIndx) = r_new
             endif  

             if(repul_C .ne. 0d0) then
               LJ = (1d0/r_new)**6
               LJ = repul_C * LJ
               dETable(iIndx) = dETable(iIndx) + LJ
               dETable(jIndx) = dETable(jIndx) + LJ
               if(.not. distCriteria) then   
                 PairList(jIndx) = PairList(jIndx) + LJ
               endif
               E_LJ = E_LJ + LJ
             endif

             r_new = sqrt(r_new)
             Ele = q_ij/r_new
             if(implcSolvent) then
               born1 = bornRad(atmType1)
               born2 = bornRad(atmType2)
               Ele = Ele + solventFunction(r_new, q_ij, born1, born2)
             endif
             dETable(iIndx) = dETable(iIndx) + Ele
             dETable(jIndx) = dETable(jIndx) + Ele
             if(.not. distCriteria) then   
               PairList(jIndx) = PairList(jIndx) + Ele
             endif
             E_Ele = E_Ele + Ele


             if(delta .ne. 0d0) then
               Morse = 1d0-exp(-alpha*(r_new-r_eq))
               Morse = delta*(Morse*Morse - 1d0)
               dETable(iIndx) = dETable(iIndx) + Morse
               dETable(jIndx) = dETable(jIndx) + Morse
               if(.not. distCriteria) then   
                 PairList(jIndx) = PairList(jIndx) + Morse
               endif
               E_Morse = E_Morse + Morse
             endif
   
!            Distance for the Old position
             rx = disp(iDisp)%x_old - MolArray(jType)%mol(jMol)%x(1)
             ry = disp(iDisp)%y_old - MolArray(jType)%mol(jMol)%y(1)
             rz = disp(iDisp)%z_old - MolArray(jType)%mol(jMol)%z(1)
             r_old = rx*rx + ry*ry + rz*rz              

             if(repul_C .ne. 0d0) then
               LJ = (1d0/r_old)**6
               LJ = repul_C * LJ
               dETable(iIndx) = dETable(iIndx) - LJ
               dETable(jIndx) = dETable(jIndx) - LJ
               E_LJ = E_LJ - LJ
             endif
 
             r_old= sqrt(r_old)
             Ele = q_ij/r_old
             if(implcSolvent) then
               born1 = bornRad(atmType1)
               born2 = bornRad(atmType2)
               Ele = Ele + solventFunction(r_old, q_ij, born1, born2)
             endif
             dETable(iIndx) = dETable(iIndx) - Ele
             dETable(jIndx) = dETable(jIndx) - Ele
             E_Ele = E_Ele - Ele

             if(delta .ne. 0d0) then
               Morse = 1d0-exp(-alpha*(r_old-r_eq))
               Morse = delta*(Morse*Morse - 1d0)
               dETable(iIndx) = dETable(iIndx) - Morse
               dETable(jIndx) = dETable(jIndx) - Morse
               E_Morse = E_Morse - Morse
             endif
          enddo
        enddo
      enddo


     
      E_Trial = E_LJ + E_Ele + E_Morse
      
      
      end subroutine

!======================================================================================      
      pure subroutine Mol_ECalc_Inter(iType, iMol, dETable, E_Trial)
      use ForceField
      use ForceFieldPara_Pedone
      use Coords
      use SimParameters
      implicit none
      integer, intent(in) :: iType, iMol     
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: dETable(:)
      
      integer :: iIndx,jType,jIndx,jMol
      integer(kind=atomIntType)  :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: r_eq, repul_C, q_ij
      real(dp) :: alpha, delta
      real(dp) :: LJ, Ele, Morse
      real(dp) :: E_Ele, E_LJ, E_Morse, E_Solvent
      real(dp) :: born1, born2
      
      E_LJ = 0d0
      E_Ele = 0d0      
      E_Morse = 0d0
      E_Trial = 0d0
      dETable = 0d0
      
       !This section calculates the Intermolecular interaction between the atoms that
       !have been modified in this trial move with the atoms that have remained stationary
      atmType1 = atomArray(iType, 1)
      iIndx = MolArray(iType) % mol(iMol) % indx
      do jType = 1, nMolTypes
        atmType2 = atomArray(jType, 1)
        r_eq = rEq_tab(atmType1, atmType2)
        q_ij = q_tab(atmType1, atmType2)
        alpha = alpha_Tab(atmType1, atmType2)
        delta = D_Tab(atmType1, atmType2)
        repul_C = repul_tab(atmType1, atmType2)          
        do jMol = 1, NPART(jType)
          if(iType .eq. jType) then
            if(iMol .eq. jMol) then
              cycle
            endif
          endif  
          rx = MolArray(iType)%mol(iMol)%x(1) - MolArray(jType)%mol(jMol)%x(1)
          ry = MolArray(iType)%mol(iMol)%y(1) - MolArray(jType)%mol(jMol)%y(1)
          rz = MolArray(iType)%mol(iMol)%z(1) - MolArray(jType)%mol(jMol)%z(1)
          r = rx*rx + ry*ry + rz*rz              

          jIndx = MolArray(jType)%mol(jMol)%indx  
          if(repul_C .ne. 0d0) then
            LJ = (1d0/r)**6
            LJ = repul_C * LJ
            dETable(iIndx) = dETable(iIndx) + LJ
            dETable(jIndx) = dETable(jIndx) + LJ
            E_LJ = E_LJ + LJ
          endif

          r = sqrt(r)
          Ele = q_ij/r
          if(implcSolvent) then
            born1 = bornRad(atmType1)
            born2 = bornRad(atmType2)
            Ele = Ele + solventFunction(r, q_ij, born1, born2)
          endif
          dETable(iIndx) = dETable(iIndx) + Ele
          dETable(jIndx) = dETable(jIndx) + Ele
          E_Ele = E_Ele + Ele

          if(delta .ne. 0d0) then
            Morse = 1d0-exp(-alpha*(r - r_eq))
            Morse = delta*(Morse*Morse - 1d0)
            dETable(iIndx) = dETable(iIndx) + Morse
            dETable(jIndx) = dETable(jIndx) + Morse
            E_Morse = E_Morse + Morse
          endif
        enddo
      enddo

!      if(implcSolvent) then
!        q_ij = q_tab(atmType1, atmType1)
!        born1 = bornRad(atmType1)
!        E_Solvent = -(1d0-1d0/dieletric)*q/(born1*born1)
!      endif
     
      E_Trial = E_LJ + E_Ele + E_Morse + E_Solvent
      
      
      end subroutine

!======================================================================================      
      pure subroutine NewMol_ECalc_Inter(E_Trial, PairList, dETable, rejMove)
      use ForceField
      use ForceFieldPara_Pedone
      use Coords
      use SimParameters
      implicit none
      logical, intent(out) :: rejMove
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      
      integer :: iAtom, iIndx, jType, jIndx, jMol, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2

      real(dp) :: rx,ry,rz
      real(dp) :: r
      real(dp) :: r_min1_sq      
      real(dp) :: r_eq, repul_C, q_ij
      real(dp) :: alpha, delta
      real(dp) :: LJ, Ele, Morse
      real(dp) :: E_Ele,E_LJ, E_Morse, E_Solvent
      real(dp) :: rmin_ij    
      real(dp) :: born1, born2

      E_LJ = 0d0
      E_Ele = 0d0      
      E_Morse = 0d0      
      E_Trial = 0d0
      PairList = 0d0      
      dETable = 0d0

      iIndx = molArray(newMol%molType)%mol(NPART(newMol%molType)+1)%indx
      atmType1 = atomArray(newMol%molType, 1)
      do jType = 1, nMolTypes
        atmType2 = atomArray(jType, 1)
        r_eq = rEq_tab(atmType1, atmType2)
        q_ij = q_tab(atmType1, atmType2)
        alpha = alpha_Tab(atmType1, atmType2)
        delta = D_Tab(atmType1, atmType2)
        repul_C = repul_tab(atmType1, atmType2)          
        rmin_ij = r_min_tab(atmType1, atmType2)  
        do jMol = 1,NPART(jType)
          rx = newMol%x(1) - MolArray(jType)%mol(jMol)%x(1)
          ry = newMol%y(1) - MolArray(jType)%mol(jMol)%y(1)
          rz = newMol%z(1) - MolArray(jType)%mol(jMol)%z(1)
          r = rx*rx + ry*ry + rz*rz
          if(r .lt. rmin_ij) then
            rejMove = .true.
            return
          endif   
          jIndx = molArray(jType)%mol(jMol)%indx
          if(distCriteria) then          
            PairList(jIndx) = r
          endif
          if(repul_C .ne. 0d0) then
            LJ = ( 1d0/r )**6
            LJ = repul_C * LJ
            dETable(iIndx) = dETable(iIndx) + LJ
            dETable(jIndx) = dETable(jIndx) + LJ
            if(.not. distCriteria) then  
              PairList(jIndx) = PairList(jIndx) + LJ
            endif
            E_LJ = E_LJ + LJ
          endif

          r = sqrt(r)
          Ele = q_ij/r
          if(implcSolvent) then
            born1 = bornRad(atmType1)
            born2 = bornRad(atmType2)
            Ele = Ele + solventFunction(r, q_ij, born1, born2)
          endif
          dETable(iIndx) = dETable(iIndx) + Ele
          dETable(jIndx) = dETable(jIndx) + Ele
          if(.not. distCriteria) then  
            PairList(jIndx) = PairList(jIndx) + Ele
          endif
          E_Ele = E_Ele + Ele

          if(delta .ne. 0d0) then
            Morse = 1d0-exp(-alpha*(r-r_eq))
            Morse = delta*(Morse*Morse - 1d0)
            dETable(iIndx) = dETable(iIndx) + Morse
            dETable(jIndx) = dETable(jIndx) + Morse
            if(.not. distCriteria) then  
              PairList(jIndx) = PairList(jIndx) + Morse
            endif
            E_Morse = E_Morse + Morse
          endif
        enddo
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
      
      integer :: iAtom,jAtom
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
