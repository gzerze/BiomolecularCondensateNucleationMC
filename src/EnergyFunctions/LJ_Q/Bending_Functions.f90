!*********************************************************************************************************************
!     This file contains the energy functions that calculate the energy associated with bond
!     stretching. 
!     The subroutine's prefix naming scheme implies the following:
!           Detailed - Complete energy calculation intended for use at the beginning and end
!                      of the simulation.  This function is not intended for use mid-simulation.
!             Shift  - Calculates the energy difference for any move that does not result
!                      in molecules being added or removed from the cluster. This function
!                      receives any number of Displacement vectors from the parent function as input.
!              Mol   - Calculates the energy for a molecule already present in the system. For
!                      use in moves such as particle deletion moves. 
!             NewMol - Calculates the energy for a molecule that has been freshly inserted into the system.
!                      Intended for use in Rosenbluth Sampling, Swap In, etc.
!*********************************************************************************************************************
      module BendingFunctions
      contains
!======================================================================================      
      subroutine Detailed_ECalc_Bending(E_T)
      use ParallelVar
      use SimParameters 
      use ForceField
      use Coords
      use ForceFieldFunctions 
      use Constants
      use EnergyTables, only: E_Bend_T
      implicit none
      real(dp), intent(out) :: E_T
      integer :: iType,iMol,iBend
      integer :: bendType      
      integer :: memb1, memb2, memb3
      real(dp) :: k_eq,ang_eq
      real(dp) :: rx12,ry12,rz12,r12
      real(dp) :: rx23,ry23,rz23,r23
      real(dp) :: Angle      
      real(dp) :: E_Bend, E_Harmonic

      E_Bend = 0d0
      E_Bend_T = 0d0
      do iType = 1, nMolTypes
        do iMol = 1,NPART(iType)
          do iBend = 1,nAngles(iType)
            bendType = bendArray(iType,iBend)%bendType
            if(bendData(bendType)%k_eq .eq. 0d0) then
              cycle
            endif
            memb1 = bendArray(iType,iBend)%bendMembr(1)
            memb2 = bendArray(iType,iBend)%bendMembr(2)
            memb3 = bendArray(iType,iBend)%bendMembr(3)
            k_eq = bendData(bendType)%k_eq
            ang_eq = bendData(bendType)%ang_eq
            
            rx12 = MolArray(iType)%mol(iMol)%x(memb1) - MolArray(iType)%mol(iMol)%x(memb2)
            ry12 = MolArray(iType)%mol(iMol)%y(memb1) - MolArray(iType)%mol(iMol)%y(memb2)
            rz12 = MolArray(iType)%mol(iMol)%z(memb1) - MolArray(iType)%mol(iMol)%z(memb2)
            r12 = sqrt(rx12**2 + ry12**2 + rz12**2)
            
            rx23 = MolArray(iType)%mol(iMol)%x(memb3) - MolArray(iType)%mol(iMol)%x(memb2)
            ry23 = MolArray(iType)%mol(iMol)%y(memb3) - MolArray(iType)%mol(iMol)%y(memb2)
            rz23 = MolArray(iType)%mol(iMol)%z(memb3) - MolArray(iType)%mol(iMol)%z(memb2)          
            r23 = sqrt(rx23**2 + ry23**2 + rz23**2)
            
            Angle = rx12*rx23 + ry12*ry23 + rz12*rz23
            Angle = Angle/(r12*r23)
            if(abs(Angle) .gt. 1d0) then
              Angle = sign(real(1,dp), Angle)
            endif
            Angle = acos(Angle)
            E_Harmonic = Harmonic(Angle, k_eq, ang_eq)
            E_Bend = E_Bend + E_Harmonic
          enddo
        enddo
      enddo

      
      write(nout,*) "Bending Energy:", E_Bend
      E_T = E_T + E_Bend
      E_Bend_T = E_Bend
      
      end subroutine
!======================================================================================      
      pure subroutine Shift_ECalc_Bending(E_Bend,disp)
      use SimParameters
      use ForceField
      use Coords
      use ForceFieldFunctions
      implicit none
      
      real(dp), intent(out) :: E_Bend
      type(Displacement), intent(in) :: disp(:)

      logical :: changed(1:maxAtoms)
      integer :: dispIndx(1:maxAtoms)
      integer :: i      
      integer :: nDisp      
      integer :: iType,iMol,iBend
      integer :: bendType
      integer :: memb1, memb2, memb3
      real(dp) :: k_eq,ang_eq
      real(dp) :: rx12,ry12,rz12,r12
      real(dp) :: rx23,ry23,rz23,r23
      real(dp) :: x1_new,y1_new,z1_new
      real(dp) :: x1_old,y1_old,z1_old
      real(dp) :: x2_new,y2_new,z2_new
      real(dp) :: x2_old,y2_old,z2_old
      real(dp) :: x3_new,y3_new,z3_new
      real(dp) :: x3_old,y3_old,z3_old      
      real(dp) :: Angle
      
      nDisp = size(disp)
      iType = disp(1)%molType
      iMol  = disp(1)%molIndx
      changed = .false.
      E_Bend = 0d0
      dispIndx = 0
      do i=1,nDisp
        changed(disp(i)%atmIndx) = .true.
        dispIndx(disp(i)%atmIndx) = i
      enddo
      
      do iBend = 1,nAngles(iType)
        bendType = bendArray(iType,iBend)%bendType
        if(bendData(bendType)%k_eq .eq. 0d0) then
          cycle
        endif
        memb1 = bendArray(iType,iBend)%bendMembr(1)
        memb2 = bendArray(iType,iBend)%bendMembr(2)
        memb3 = bendArray(iType,iBend)%bendMembr(3)
        if(.not. changed(memb1)) then
          if(.not. changed(memb2)) then
            if(.not. changed(memb3)) then         
              cycle
            endif
          endif
        endif
        k_eq = bendData(bendType)%k_eq
        ang_eq = bendData(bendType)%ang_eq

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

        if(changed(memb3) ) then
          x3_new = disp(dispIndx(memb3))%x_new
          y3_new = disp(dispIndx(memb3))%y_new
          z3_new = disp(dispIndx(memb3))%z_new
          
          x3_old = disp(dispIndx(memb3))%x_old
          y3_old = disp(dispIndx(memb3))%y_old
          z3_old = disp(dispIndx(memb3))%z_old
        else
          x3_new = MolArray(iType)%mol(iMol)%x(memb3)
          y3_new = MolArray(iType)%mol(iMol)%y(memb3)
          z3_new = MolArray(iType)%mol(iMol)%z(memb3)
          
          x3_old = x3_new
          y3_old = y3_new
          z3_old = z3_new
        endif        

!        write(2,*) x1_new, x2_new
!        write(2,*) y1_new, y2_new
!        write(2,*) z1_new, z2_new
        
        rx12 = x1_new - x2_new
        ry12 = y1_new - y2_new
        rz12 = z1_new - z2_new
        r12 = sqrt(rx12**2 + ry12**2 + rz12**2)

        rx23 = x3_new - x2_new
        ry23 = y3_new - y2_new
        rz23 = z3_new - z2_new
        r23 = sqrt(rx23**2 + ry23**2 + rz23**2)        
       
        Angle = rx12*rx23 + ry12*ry23 + rz12*rz23
        Angle = Angle/(r12*r23)
        if(abs(Angle) .gt. 1d0) then
          Angle = sign(real(1,dp), Angle)
        endif
        Angle = acos(Angle)
        E_Bend = E_Bend + Harmonic(Angle, k_eq, ang_eq)
        
        rx12 = x1_old - x2_old
        ry12 = y1_old - y2_old
        rz12 = z1_old - z2_old
        r12 = sqrt(rx12**2 + ry12**2 + rz12**2)
        
        rx23 = x3_old - x2_old
        ry23 = y3_old - y2_old
        rz23 = z3_old - z2_old
        r23 = sqrt(rx23**2 + ry23**2 + rz23**2)        
       
        Angle = rx12*rx23 + ry12*ry23 + rz12*rz23
        Angle = Angle/(r12*r23)
        if(abs(Angle) .gt. 1d0) then
          Angle = sign(real(1,dp), Angle)
        endif
        Angle = acos(Angle)
        E_Bend = E_Bend - Harmonic(Angle, k_eq, ang_eq)
      enddo
      
      
      end subroutine
!======================================================================================      
      pure subroutine Mol_ECalc_Bending(iType,iMol,E_Bend)
      use SimParameters 
      use ForceField
      use Coords
      use ForceFieldFunctions 
      implicit none
      real(dp), intent(out) :: E_Bend
      integer, intent(in) :: iType, iMol
      integer :: iBend
      integer :: bendType
      integer :: memb1, memb2, memb3
      real(dp) :: k_eq,ang_eq
      real(dp) :: rx12,ry12,rz12,r12
      real(dp) :: rx23,ry23,rz23,r23
      real(dp) :: Angle    

      E_Bend = 0d0
      do iBend = 1,nAngles(iType)
        bendType = bendArray(iType,iBend)%bendType
        if(bendData(bendType)%k_eq .eq. 0d0) then
           cycle
        endif
        memb1 = bendArray(iType,iBend)%bendMembr(1)
        memb2 = bendArray(iType,iBend)%bendMembr(2)
        memb3 = bendArray(iType,iBend)%bendMembr(3)
        k_eq = bendData(bendType)%k_eq
        ang_eq = bendData(bendType)%ang_eq
            
        rx12 = MolArray(iType)%mol(iMol)%x(memb1) - MolArray(iType)%mol(iMol)%x(memb2)
        ry12 = MolArray(iType)%mol(iMol)%y(memb1) - MolArray(iType)%mol(iMol)%y(memb2)
        rz12 = MolArray(iType)%mol(iMol)%z(memb1) - MolArray(iType)%mol(iMol)%z(memb2)
        r12 = sqrt(rx12*rx12 + ry12*ry12 + rz12*rz12)
            
        rx23 = MolArray(iType)%mol(iMol)%x(memb3) - MolArray(iType)%mol(iMol)%x(memb2)
        ry23 = MolArray(iType)%mol(iMol)%y(memb3) - MolArray(iType)%mol(iMol)%y(memb2)
        rz23 = MolArray(iType)%mol(iMol)%z(memb3) - MolArray(iType)%mol(iMol)%z(memb2)          
        r23 = sqrt(rx23*rx23 + ry23*ry23 + rz23*rz23)
            
        Angle = rx12*rx23 + ry12*ry23 + rz12*rz23
        Angle = Angle/(r12*r23)
        if(abs(Angle) .gt. 1d0) then
           Angle = sign(real(1,dp), Angle)
        endif
        Angle = acos(Angle)
            
        E_Bend = E_Bend + Harmonic(Angle, k_eq, ang_eq)
      enddo
      
      end subroutine
    
!======================================================================================      
      pure subroutine NewMol_ECalc_Bending(E_Bend)
      use SimParameters
      use ForceField
      use Coords
      use ForceFieldFunctions 
      implicit none
      real(dp), intent(out) :: E_Bend
      integer :: iType
      integer :: iBend
      integer :: bendType
      integer :: memb1, memb2, memb3
      real(dp) :: k_eq,ang_eq
      real(dp) :: rx12,ry12,rz12,r12
      real(dp) :: rx23,ry23,rz23,r23
      real(dp) :: Angle    

      E_Bend = 0d0
      iType = newMol%molType
      do iBend = 1,nAngles(iType)
        bendType = bendArray(iType,iBend)%bendType
        if(bendData(bendType)%k_eq .eq. 0d0) then
          cycle
        endif
        memb1 = bendArray(iType,iBend)%bendMembr(1)
        memb2 = bendArray(iType,iBend)%bendMembr(2)
        memb3 = bendArray(iType,iBend)%bendMembr(3)
        k_eq = bendData(bendType)%k_eq
        ang_eq = bendData(bendType)%ang_eq
            
        rx12 = newMol%x(memb1) - newMol%x(memb2)
        ry12 = newMol%y(memb1) - newMol%y(memb2)
        rz12 = newMol%z(memb1) - newMol%z(memb2)
        r12 = sqrt(rx12*rx12 + ry12*ry12 + rz12*rz12)
        
        rx23 = newMol%x(memb3) - newMol%x(memb2)
        ry23 = newMol%y(memb3) - newMol%y(memb2)
        rz23 = newMol%z(memb3) - newMol%z(memb2)
        r23 = sqrt(rx23*rx23 + ry23*ry23 + rz23*rz23)

        Angle = rx12*rx23 + ry12*ry23 + rz12*rz23
        Angle = Angle/(r12*r23)
        if(abs(Angle) .gt. 1d0) then
           Angle = sign(real(1,dp), Angle)
        endif
        Angle = acos(Angle)
        E_Bend = E_Bend + Harmonic(Angle, k_eq, ang_eq)
      enddo

      end subroutine      
!======================================================================================
!  This function calculates the energy of a currently existing bending angle.      
      pure subroutine Single_ECalc_Bending(iType, iMol, iBend, E_Bend)
      use SimParameters 
      use ForceField
      use Coords
      use ForceFieldFunctions 
      implicit none
      real(dp), intent(out) :: E_Bend
      integer, intent(in) :: iType, iMol,iBend
      integer :: bendType
      integer :: memb1, memb2, memb3
      real(dp) :: k_eq,ang_eq
      real(dp) :: rx12,ry12,rz12,r12
      real(dp) :: rx23,ry23,rz23,r23
      real(dp) :: Angle    

      E_Bend = 0d0
      bendType = bendArray(iType,iBend)%bendType
      if(bendData(bendType)%k_eq .eq. 0d0) then
        return
      endif
      memb1 = bendArray(iType,iBend)%bendMembr(1)
      memb2 = bendArray(iType,iBend)%bendMembr(2)
      memb3 = bendArray(iType,iBend)%bendMembr(3)
      k_eq = bendData(bendType)%k_eq
      ang_eq = bendData(bendType)%ang_eq
           
      rx12 = MolArray(iType)%mol(iMol)%x(memb1) - MolArray(iType)%mol(iMol)%x(memb2)
      ry12 = MolArray(iType)%mol(iMol)%y(memb1) - MolArray(iType)%mol(iMol)%y(memb2)
      rz12 = MolArray(iType)%mol(iMol)%z(memb1) - MolArray(iType)%mol(iMol)%z(memb2)
      r12 = sqrt(rx12**2 + ry12**2 + rz12**2)
            
      rx23 = MolArray(iType)%mol(iMol)%x(memb3) - MolArray(iType)%mol(iMol)%x(memb2)
      ry23 = MolArray(iType)%mol(iMol)%y(memb3) - MolArray(iType)%mol(iMol)%y(memb2)
      rz23 = MolArray(iType)%mol(iMol)%z(memb3) - MolArray(iType)%mol(iMol)%z(memb2)          
      r23 = sqrt(rx23**2 + ry23**2 + rz23**2)
            
      Angle = rx12*rx23 + ry12*ry23 + rz12*rz23
      Angle = Angle/(r12*r23)
      if(abs(Angle) .gt. 1d0) then
         Angle = sign(real(1,dp), Angle)
      endif
      Angle = acos(Angle)
            
      E_Bend = E_Bend + Harmonic(Angle, k_eq, ang_eq)
     
      end subroutine  
      
!======================================================================================
      end module
      
      
