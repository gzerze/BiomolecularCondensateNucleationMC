!*********************************************************************************************************************
!     This file contains the energy functions that work for Lennard-Jones w/ Columbic style forcefields
!     these functions are enclosed inside of the module "InterMolecularEnergy" so that
!     the energy functions can be freely exchanged from the simulation.
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
      module IntraEnergy_LJ_Electro
      contains
!======================================================================================      
      subroutine Detailed_ECalc_IntraNonBonded(E_T)
      use ParallelVar
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      use EnergyTables, only: E_NBond_T
      implicit none
      real(dp), intent(inout) :: E_T
      integer :: iType,iMol,iPair,iAtom,jAtom
      integer :: atmType1,atmType2      
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      
      E_LJ = 0d0
      E_Ele = 0d0      
      
      do iType = 1,nMolTypes
        do iMol=1,NPART(iType)
          do iPair = 1,nIntraNonBond(iType)
            iAtom = nonBondArray(iType,iPair)%nonMembr(1)
            jAtom = nonBondArray(iType,iPair)%nonMembr(2)
            atmType1 = atomArray(iType,iAtom)
            atmType2 = atomArray(iType,jAtom)
            
            ep = ep_tab(atmType1,atmType2)
            sig_sq = sig_tab(atmType1,atmType2)
            q = q_tab(atmType1,atmType2)
            
            rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(iType)%mol(iMol)%x(jAtom)
            ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(iType)%mol(iMol)%y(jAtom)
            rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(iType)%mol(iMol)%z(jAtom) 
            r = rx**2 + ry**2 + rz**2

            LJ = (sig_sq/r)**3
            LJ = ep * LJ * (LJ-1d0)
            E_LJ = E_LJ + LJ
            
            r = sqrt(r)
            Ele = q/r
            E_Ele = E_Ele + Ele            
          enddo
        enddo
      enddo
      write(nout,*) "Intra Lennard-Jones Energy:", E_LJ
      write(nout,*) "Intra Eletrostatic Energy:", E_Ele
      E_T = E_T + E_Ele + E_LJ    
      E_NBond_T = E_Ele + E_LJ  
      
      end subroutine
!======================================================================================      
      pure subroutine Shift_ECalc_IntraNonBonded(E_Trial,disp)
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      implicit none
      
      real(dp), intent(inout) :: E_Trial
      type(Displacement), intent(in) :: disp(:)  

      logical :: changed(1:maxAtoms)
      integer :: dispIndx(1:maxAtoms)      
      integer :: i,nDisp
      integer :: iType,iMol,iPair,iAtom,jAtom
      integer :: atmType1,atmType2      
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: x1_New, y1_New, z1_New
      real(dp) :: x1_Old, y1_Old, z1_Old
      real(dp) :: x2_New, y2_New, z2_New
      real(dp) :: x2_Old, y2_Old, z2_Old
      

      nDisp = size(disp)
      iType = disp(1)%molType
      iMol  = disp(1)%molIndx
      changed = .false.
      E_LJ = 0d0
      E_Ele = 0d0  

      
      do i=1,nDisp
        changed(disp(i)%atmIndx) = .true.
        dispIndx(disp(i)%atmIndx) = i
      enddo
      
      do iPair = 1,nIntraNonBond(iType)
        iAtom = nonBondArray(iType,iPair)%nonMembr(1)
        jAtom = nonBondArray(iType,iPair)%nonMembr(2)
        if( all(changed([iAtom,jAtom]) .eqv. .false.) ) then
           cycle
        endif
        atmType1 = atomArray(iType,iAtom)
        atmType2 = atomArray(iType,jAtom)
        ep = ep_tab(atmType1,atmType2)
        sig_sq = sig_tab(atmType1,atmType2)
        q = q_tab(atmType1,atmType2)
        
        if(changed(iAtom)) then   
          x1_Old = disp(dispIndx(iAtom))%x_old
          y1_Old = disp(dispIndx(iAtom))%y_old
          z1_Old = disp(dispIndx(iAtom))%z_old
          x1_New = disp(dispIndx(iAtom))%x_new
          y1_New = disp(dispIndx(iAtom))%y_new
          z1_New = disp(dispIndx(iAtom))%z_new
        else
          x1_Old = MolArray(iType)%mol(iMol)%x(iAtom)
          y1_Old = MolArray(iType)%mol(iMol)%y(iAtom)
          z1_Old = MolArray(iType)%mol(iMol)%z(iAtom)      
          x1_New = x1_Old
          y1_New = y1_Old
          z1_New = z1_Old
        endif        
        if(changed(jAtom)) then   
          x2_Old = disp(dispIndx(jAtom))%x_old
          y2_Old = disp(dispIndx(jAtom))%y_old
          z2_Old = disp(dispIndx(jAtom))%z_old
          x2_New = disp(dispIndx(jAtom))%x_new
          y2_New = disp(dispIndx(jAtom))%y_new
          z2_New = disp(dispIndx(jAtom))%z_new
        else
          x2_Old = MolArray(iType)%mol(iMol)%x(jAtom)
          y2_Old = MolArray(iType)%mol(iMol)%y(jAtom)
          z2_Old = MolArray(iType)%mol(iMol)%z(jAtom)      
          x2_New = x2_Old
          y2_New = y2_Old
          z2_New = z2_Old
        endif
        
        rx = x2_New - x1_New
        ry = y2_New - y1_New
        rz = z2_New - z1_New
        r = rx*rx + ry*ry + rz*rz

        LJ = (sig_sq/r)
        LJ = LJ * LJ * LJ
        LJ = 4d0 * ep * LJ * (LJ-1d0)
        E_LJ = E_LJ + LJ
            
        r = sqrt(r)
        Ele = q/r
        E_Ele = E_Ele + Ele            
        
        rx = x2_Old - x1_Old
        ry = y2_Old - y1_Old
        rz = z2_Old - z1_Old
        r = rx*rx + ry*ry + rz*rz

        LJ = (sig_sq/r)
        LJ = LJ * LJ * LJ              
        LJ = ep * LJ * (LJ-1d0)
        E_LJ = E_LJ - LJ
            
        r = sqrt(r)
        Ele = q/r
        E_Ele = E_Ele - Ele           
      enddo

      E_Trial = E_LJ + E_Ele
      
      end subroutine
!======================================================================================
      pure subroutine Mol_ECalc_IntraNonBonded(iType, iMol, E_Trial)
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      implicit none
      integer, intent(in) :: iType, iMol     
      real(dp), intent(out) :: E_Trial
     
      integer :: iPair, iAtom, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ

      E_LJ = 0d0
      E_Ele = 0d0      
      E_Trial = 0d0
      
      do iPair = 1,nIntraNonBond(iType)
        iAtom = nonBondArray(iType,iPair)%nonMembr(1)
        jAtom = nonBondArray(iType,iPair)%nonMembr(2)

        atmType1 = atomArray(iType,iAtom)
        atmType2 = atomArray(iType,jAtom)
        ep = ep_tab(atmType1,atmType2)
        sig_sq = sig_tab(atmType1,atmType2)
        q = q_tab(atmType1,atmType2)
        
        rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(iType)%mol(iMol)%x(jAtom)
        ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(iType)%mol(iMol)%y(jAtom)
        rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(iType)%mol(iMol)%z(jAtom) 
        r = rx*rx + ry*ry + rz*rz
        if(ep .ne. 0d0) then
          LJ = (sig_sq/r)
          LJ = LJ * LJ * LJ
          LJ = ep * LJ * (LJ-1d0)
          E_LJ = E_LJ + LJ
        endif
        if(q .ne. 0d0) then            
          r = sqrt(r)
          Ele = q/r
          E_Ele = E_Ele + Ele            
        endif
      enddo
     
      E_Trial = E_LJ + E_Ele
      
      
      end subroutine
!======================================================================================      
      pure subroutine NewMol_ECalc_IntraNonBonded(E_Trial)
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      implicit none
      real(dp), intent(out) :: E_Trial
     
      integer :: iType,iPair, iAtom, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ

      E_LJ = 0d0
      E_Ele = 0d0      
      E_Trial = 0d0
      iType = newMol%molType
      do iPair = 1,nIntraNonBond(iType)
        iAtom = nonBondArray(iType,iPair)%nonMembr(1)
        jAtom = nonBondArray(iType,iPair)%nonMembr(2)

        atmType1 = atomArray(iType,iAtom)
        atmType2 = atomArray(iType,jAtom)
        ep = ep_tab(atmType1,atmType2)
        sig_sq = sig_tab(atmType1,atmType2)
        q = q_tab(atmType1,atmType2)
        
        rx = newMol%x(iAtom) - newMol%x(jAtom)
        ry = newMol%y(iAtom) - newMol%y(jAtom)
        rz = newMol%z(iAtom) - newMol%z(jAtom)
        r = rx*rx + ry*ry + rz*rz
        if(ep .ne. 0d0) then
          LJ = (sig_sq/r)
          LJ = LJ * LJ * LJ
          LJ = ep * LJ * (LJ-1d0)
          E_LJ = E_LJ + LJ
        endif
        if(q .ne. 0d0) then            
          r = sqrt(r)
          Ele = q/r
          E_Ele = E_Ele + Ele            
        endif
      enddo

      E_Trial = E_LJ + E_Ele
      
      
      end subroutine    
!======================================================================================   
      end module
      
      
