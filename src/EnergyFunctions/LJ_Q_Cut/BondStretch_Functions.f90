!*********************************************************************************************************************
!     This file contains the energy functions that calculate the energy associated with bond
!     stretching. 
!     The subroutine's prefix naming scheme implies the following:
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
      module BondStretchFunctions
      contains
!======================================================================================      
      subroutine Detailed_ECalc_BondStretch(E_T)
      use ParallelVar
      use SimParameters
      use ForceField
      use Coords
      use ForceFieldFunctions
      use EnergyTables, only: E_Stretch_T      
      implicit none
      real(dp), intent(inout) :: E_T
      integer :: iType,iMol,iBond
      integer :: bondType      
      integer :: memb1, memb2      
      real(dp) :: k_eq,r_eq
      real(dp) :: rx,ry,rz,r
      real(dp) :: E_Bond, E_Harmonic

      E_Bond = 0d0
      E_Stretch_T = 0d0
      do iType = 1, nMolTypes
        do iMol = 1,NPART(iType)
          do iBond = 1,nBonds(iType)
            bondType = bondArray(iType,iBond)%bondType
            if(bondData(bondType)%k_eq .eq. 0d0) then
              cycle
            endif
            memb1 = bondArray(iType,iBond)%bondMembr(1)
            memb2 = bondArray(iType,iBond)%bondMembr(2)
            k_eq = bondData(bondType)%k_eq
            r_eq = bondData(bondType)%r_eq
            rx = MolArray(iType)%mol(iMol)%x(memb1) - MolArray(iType)%mol(iMol)%x(memb2)
            ry = MolArray(iType)%mol(iMol)%y(memb1) - MolArray(iType)%mol(iMol)%y(memb2)
            rz = MolArray(iType)%mol(iMol)%z(memb1) - MolArray(iType)%mol(iMol)%z(memb2)            
            r = sqrt(rx*rx + ry*ry + rz*rz)
            E_Harmonic = Harmonic(r, k_eq, r_eq)
            E_Bond = E_Bond + E_Harmonic
          enddo
        enddo
      enddo
      write(nout,*) "Stretch Energy:", E_Bond
      E_T = E_T + E_Bond
      E_Stretch_T = E_Bond

      end subroutine
!======================================================================================      
      pure subroutine Shift_ECalc_BondStretch(E_Bond,disp)
      use SimParameters
      use ForceField
      use Coords
      use ForceFieldFunctions
      implicit none
      
      real(dp), intent(out) :: E_Bond
      type(Displacement), intent(in) :: disp(:)  

      logical :: changed(1:maxAtoms)
      integer :: dispIndx(1:maxAtoms)
      integer :: i      
      integer :: nDisp      
      integer :: iType,iMol,iBond
      integer :: bondType
      integer :: memb1, memb2      
      real(dp) :: k_eq,r_eq
      real(dp) :: rx,ry,rz,r
      real(dp) :: x1_new,y1_new,z1_new
      real(dp) :: x1_old,y1_old,z1_old
      real(dp) :: x2_new,y2_new,z2_new
      real(dp) :: x2_old,y2_old,z2_old

	  
      nDisp = size(disp)
      iType = disp(1)%molType
      iMol  = disp(1)%molIndx
      changed = .false.
      E_Bond = 0d0
      dispIndx = 0
      do i=1,nDisp
        changed(disp(i)%atmIndx) = .true.
        dispIndx(disp(i)%atmIndx) = i
      enddo
      do iBond = 1,nBonds(iType)
        bondType = bondArray(iType,iBond)%bondType
        if(bondData(bondType)%k_eq .eq. 0d0) then
          cycle
        endif
        memb1 = bondArray(iType,iBond)%bondMembr(1)
        memb2 = bondArray(iType,iBond)%bondMembr(2)
        if(.not. changed(memb1)) then
          if(.not. changed(memb2)) then        
             cycle
          endif
        endif
        k_eq = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq

        if(changed(memb1) ) then
          x1_new = disp(dispIndx(memb1))%x_new
          y1_new = disp(dispIndx(memb1))%y_new
          z1_new = disp(dispIndx(memb1))%z_new
		  
          x1_old = disp(dispIndx(memb1))%x_old
          y1_old = disp(dispIndx(memb1))%y_old
          z1_old = disp(dispIndx(memb1))%z_old
        else
          x1_new = MolArray(iType)%mol(iMol)%x(memb1)
          y1_new = MolArray(iType)%mol(iMol)%y(memb1)
          z1_new = MolArray(iType)%mol(iMol)%z(memb1)
		  
          x1_old = x1_new
          y1_old = y1_new
          z1_old = z1_new
        endif
        
        if(changed(memb2) ) then
          x2_new = disp(dispIndx(memb2))%x_new
          y2_new = disp(dispIndx(memb2))%y_new
          z2_new = disp(dispIndx(memb2))%z_new
		  
          x2_old = disp(dispIndx(memb2))%x_old
          y2_old = disp(dispIndx(memb2))%y_old
          z2_old = disp(dispIndx(memb2))%z_old
        else
          x2_new = MolArray(iType)%mol(iMol)%x(memb2)
          y2_new = MolArray(iType)%mol(iMol)%y(memb2)
          z2_new = MolArray(iType)%mol(iMol)%z(memb2)
		  
          x2_old = x2_new
          y2_old = y2_new
          z2_old = z2_new
        endif        
        rx = x1_new - x2_new
        ry = y1_new - y2_new
        rz = z1_new - z2_new
        r = sqrt(rx*rx + ry*ry + rz*rz)
        E_Bond = E_Bond + Harmonic(r, k_eq, r_eq)
        
        rx = x1_old - x2_old
        ry = y1_old - y2_old
        rz = z1_old - z2_old
        r = sqrt(rx*rx + ry*ry + rz*rz)
        E_Bond = E_Bond - Harmonic(r, k_eq, r_eq)
      enddo
      
      
      end subroutine
!======================================================================================      
      pure subroutine Mol_ECalc_BondStretch(iType,iMol,E_Bond)
      use SimParameters
      use ForceField
      use Coords
      use ForceFieldFunctions
      implicit none
      real(dp), intent(out) :: E_Bond
      integer, intent(in) :: iType,iMol
      integer :: bondType,iBond
      integer :: memb1, memb2      
      real(dp) :: k_eq,r_eq
      real(dp) :: rx,ry,rz,r

      E_Bond = 0d0
      do iBond = 1,nBonds(iType)
        bondType = bondArray(iType,iBond)%bondType
        if(bondData(bondType)%k_eq .eq. 0d0) then
          cycle
        endif
        memb1 = bondArray(iType,iBond)%bondMembr(1)
        memb2 = bondArray(iType,iBond)%bondMembr(2)
        k_eq = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
            
        rx = MolArray(iType)%mol(iMol)%x(memb1) - MolArray(iType)%mol(iMol)%x(memb2)
        ry = MolArray(iType)%mol(iMol)%y(memb1) - MolArray(iType)%mol(iMol)%y(memb2)
        rz = MolArray(iType)%mol(iMol)%z(memb1) - MolArray(iType)%mol(iMol)%z(memb2)
        r = sqrt(rx**2 + ry**2 + rz**2)
        E_Bond = E_Bond + Harmonic(r, k_eq, r_eq)
      enddo
      
      end subroutine
!======================================================================================      
      pure subroutine NewMol_ECalc_BondStretch(E_Bond)
      use SimParameters
      use ForceField
      use Coords
      use ForceFieldFunctions
      implicit none
      real(dp), intent(out) :: E_Bond
      integer :: iType,iBond
      integer :: bondType
      integer :: memb1, memb2      
      real(dp) :: k_eq,r_eq
      real(dp) :: rx,ry,rz,r

      E_Bond = 0d0
      iType = newMol%molType
      do iBond = 1,nBonds(iType)
        bondType = bondArray(iType,iBond)%bondType
        if(bondData(bondType)%k_eq .eq. 0d0) then
          cycle
        endif
        memb1 = bondArray(iType,iBond)%bondMembr(1)
        memb2 = bondArray(iType,iBond)%bondMembr(2)
        k_eq = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
            
        rx = newMol%x(memb1) - newMol%x(memb2)
        ry = newMol%y(memb1) - newMol%y(memb2)
        rz = newMol%z(memb1) - newMol%z(memb2)
        r = sqrt(rx**2 + ry**2 + rz**2)
        E_Bond = E_Bond + Harmonic(r, k_eq, r_eq)
      enddo
      
      end subroutine
	  
!======================================================================================
      end module
      
      
