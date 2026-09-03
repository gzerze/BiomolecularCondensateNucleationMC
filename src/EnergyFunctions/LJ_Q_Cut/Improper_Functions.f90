!*********************************************************************************************************************
!     This file contains the energy functions that calculate the energy associated with bond
!     stretching. 
!     The subroutine's prefix naming scheme implies the following:
!           Detailed - Complete energy calculation intended for use at the beginning and end
!                      of the simulation.  This function is not intended for use mid-simulation.
!             Shift  - Calculates the energy difference for any move that does not result
!                      in molecules being added or removed from the cluster. This function
!                      receives any number of Displacement vectors from the parent function as input and
!                      calculates only the components that have changed. 
!              Mol   - Calculates the energy for a molecule already present in the system. For
!                      use in moves such as particle deletion moves. 
!             NewMol - Calculates the energy for a molecule that has been freshly inserted into the system.
!                      Intended for use in Rosenbluth Sampling, Swap In, etc.
!*********************************************************************************************************************
      module ImproperAngleFunctions
      contains
!======================================================================================      
      pure subroutine Detailed_ECalc_Improper(E_T)
      use SimParameters
      use ForceField
      use Coords
      use ForceFieldFunctions
      implicit none
      real(dp), intent(inout) :: E_T
      integer :: iType,iMol,iImprop
      integer(kind=atomIntType) :: impropType      
      integer :: memb1, memb2, memb3, memb4
      real(dp) :: x12,y12,z12
      real(dp) :: x23,y23,z23
      real(dp) :: x34,y34,z34
      real(dp) :: vx1,vy1,vz1
      real(dp) :: vx2,vy2,vz2
      real(dp) :: vx3,vy3,vz3
      real(dp) :: r1,r3,dot1,dot2
      real(dp) :: Angle      
      real(dp) :: E_Improp

      E_Improp = 0d0
      do iType = 1, nMolTypes
        do iMol = 1,NPART(iType)
          do iImprop = 1,nTorsional(iType)
            impropType = impropArray(iType,iImprop)%impropType
            
            memb1 = impropArray(iType,iImprop)%impropMembr(1)
            memb2 = impropArray(iType,iImprop)%impropMembr(2)
            memb3 = impropArray(iType,iImprop)%impropMembr(3)
            memb4 = impropArray(iType,iImprop)%impropMembr(4)
            
            x12 = MolArray(iType)%mol(iMol)%x(memb2) - MolArray(iType)%mol(iMol)%x(memb1)
            y12 = MolArray(iType)%mol(iMol)%y(memb2) - MolArray(iType)%mol(iMol)%y(memb1)
            z12 = MolArray(iType)%mol(iMol)%z(memb2) - MolArray(iType)%mol(iMol)%z(memb1)
            
            x23 = MolArray(iType)%mol(iMol)%x(memb3) - MolArray(iType)%mol(iMol)%x(memb2)
            y23 = MolArray(iType)%mol(iMol)%y(memb3) - MolArray(iType)%mol(iMol)%y(memb2)
            z23 = MolArray(iType)%mol(iMol)%z(memb3) - MolArray(iType)%mol(iMol)%z(memb2)          

            x34 = MolArray(iType)%mol(iMol)%x(memb4) - MolArray(iType)%mol(iMol)%x(memb3)
            y34 = MolArray(iType)%mol(iMol)%y(memb4) - MolArray(iType)%mol(iMol)%y(memb3)
            z34 = MolArray(iType)%mol(iMol)%z(memb4) - MolArray(iType)%mol(iMol)%z(memb3)          

!           Calculate the vector, v1 = <r12 x r23>
            vx1 =   y12*z23 - z12*y23
            vy1 = -(x12*z23 - z12*x23)
            vz1 =   x12*y23 - y12*x23
           
!           Calculate the vector, v2=<r3 x r2>
            vx2 =   y23*z34 - z23*y34
            vy2 = -(x23*z34 - z23*x34)
            vz2 =   x23*y34 - y23*x34   

!           Calculate the vector, v3=<v1 x r2>
            vx3 =   vy1*z23 - vz1*y23
            vy3 = -(vx1*z23 - vz1*x23)
            vz3 =   vx1*y23 - vy1*x23            
            
            r1 = vx1**2 + vy1**2 + vz1**2
            r3 = vx3**2 + vy3**2 + vz3**2
            r1 = sqrt(r1)
            r3 = sqrt(r3)

            dot1 = vx1*vx2 + vy1*vy2 + vz1*vz2
            dot2 = vx2*vx3 + vy2*vy3 + vz2*vz3            
            dot1 = dot1/(r1)
            dot2 = dot2/(r3)    
            angle = atan2(dot2,dot1)            
            
            E_Improp = E_Improp + TorsHarmonic(Angle, impropData(impropType)%a)
          enddo
        enddo
      enddo

      E_T = E_T + E_Improp
      
      end subroutine
!======================================================================================      
      pure subroutine Shift_ECalc_Improper(E_Improp,disp)
      use SimParameters
      use ForceField
      use Coords
      use ForceFieldFunctions
      implicit none
      
      real(dp), intent(out) :: E_Improp
      type(Displacement), intent(in) :: disp(:)


      logical :: changed(1:maxAtoms)
      integer :: dispIndx(1:maxAtoms)      
      integer :: i,nDisp,iType,iMol,iImprop
      integer :: impropType      
      integer :: memb(1:4)
      real(dp) :: x_new(1:4),y_new(1:4),z_new(1:4)  
      real(dp) :: x_old(1:4),y_old(1:4),z_old(1:4)  
      real(dp) :: x12,y12,z12
      real(dp) :: x23,y23,z23
      real(dp) :: x34,y34,z34
      real(dp) :: vx1,vy1,vz1
      real(dp) :: vx2,vy2,vz2
      real(dp) :: vx3,vy3,vz3
      real(dp) :: r1,r3,dot1,dot2
      real(dp) :: Angle      

      nDisp = size(disp)
      iType = disp(1)%molType
      iMol  = disp(1)%molIndx
      changed = .false.
      E_Improp = 0d0
      
      do i=1,nDisp
        changed(disp(i)%atmIndx) = .true.
        dispIndx(disp(i)%atmIndx) = i
      enddo
      
      do iImprop = 1,nTorsional(iType)
        impropType = impropArray(iType,iImprop)%impropType
        memb(1) = impropArray(iType,iImprop)%impropMembr(1)
        memb(2) = impropArray(iType,iImprop)%impropMembr(2)
        memb(3) = impropArray(iType,iImprop)%impropMembr(3)
        memb(4) = impropArray(iType,iImprop)%impropMembr(4)
        if(all( changed([memb(1:4)]) .eqv. .false. )) then
           cycle
        endif
        
        do i=1,4
          if(changed(memb(i)) ) then
            x_new(i) = disp(dispIndx(memb(i)))%x_new
            y_new(i) = disp(dispIndx(memb(i)))%y_new
            z_new(i) = disp(dispIndx(memb(i)))%z_new
            
            x_old(i) = disp(dispIndx(memb(i)))%x_old
            y_old(i) = disp(dispIndx(memb(i)))%y_old
            z_old(i) = disp(dispIndx(memb(i)))%z_old
          else
            x_new(i) = MolArray(iType)%mol(iMol)%x(memb(i))
            y_new(i) = MolArray(iType)%mol(iMol)%y(memb(i))
            z_new(i) = MolArray(iType)%mol(iMol)%z(memb(i))
            
            x_new(i) = x_old(i)
            y_new(i) = y_old(i)
            z_new(i) = z_old(i)
          endif
        enddo
        
        x12 = x_new(2) - x_new(1)
        y12 = y_new(2) - y_new(1)
        z12 = z_new(2) - z_new(1)
            
        x23 = x_new(3) - x_new(2)
        y23 = y_new(3) - y_new(2)
        z23 = z_new(3) - z_new(2)

        x34 = x_new(4) - x_new(3)
        y34 = y_new(4) - y_new(3)
        z34 = z_new(4) - z_new(3)

!       Calculate the vector, v1 = <r12 x r23>
        vx1 =   y12*z23 - z12*y23
        vy1 = -(x12*z23 - z12*x23)
        vz1 =   x12*y23 - y12*x23
           
!           Calculate the vector, v2=<r3 x r2>
        vx2 =   y23*z34 - z23*y34
        vy2 = -(x23*z34 - z23*x34)
        vz2 =   x23*y34 - y23*x34   

!       Calculate the vector, v3=<v1 x r2>
        vx3 =   vy1*z23 - vz1*y23
        vy3 = -(vx1*z23 - vz1*x23)
        vz3 =   vx1*y23 - vy1*x23            
          
        r1 = vx1**2 + vy1**2 + vz1**2
        r3 = vx3**2 + vy3**2 + vz3**2
        r1 = sqrt(r1)
        r3 = sqrt(r3)

        dot1 = vx1*vx2 + vy1*vy2 + vz1*vz2
        dot2 = vx2*vx3 + vy2*vy3 + vz2*vz3            
        dot1 = dot1/(r1)
        dot2 = dot2/(r3)   
        angle = atan2(dot2,dot1)            
            
        E_Improp = E_Improp + TorsHarmonic(Angle, impropData(impropType)%a)
        
        x12 = x_old(2) - x_old(1)
        y12 = y_old(2) - y_old(1)
        z12 = z_old(2) - z_old(1)
            
        x23 = x_old(3) - x_old(2)
        y23 = y_old(3) - y_old(2)
        z23 = z_old(3) - z_old(2)

        x34 = x_old(4) - x_old(3)
        y34 = y_old(4) - y_old(3)
        z34 = z_old(4) - z_old(3)

!       Calculate the vector, v1 = <r12 x r23>
        vx1 =   y12*z23 - z12*y23
        vy1 = -(x12*z23 - z12*x23)
        vz1 =   x12*y23 - y12*x23
           
!           Calculate the vector, v2=<r3 x r2>
        vx2 =   y23*z34 - z23*y34
        vy2 = -(x23*z34 - z23*x34)
        vz2 =   x23*y34 - y23*x34   

!       Calculate the vector, v3=<v1 x r2>
        vx3 =   vy1*z23 - vz1*y23
        vy3 = -(vx1*z23 - vz1*x23)
        vz3 =   vx1*y23 - vy1*x23            
          
        r1 = vx1**2 + vy1**2 + vz1**2
        r3 = vx3**2 + vy3**2 + vz3**2
        r1 = sqrt(r1)
        r3 = sqrt(r3)

        dot1 = vx1*vx2 + vy1*vy2 + vz1*vz2
        dot2 = vx2*vx3 + vy2*vy3 + vz2*vz3            
        dot1 = dot1/(r1)
        dot2 = dot2/(r3)   
        angle = atan2(dot2,dot1)            
            
        E_Improp = E_Improp - TorsHarmonic(Angle, impropData(impropType)%a)
      enddo
      
      end subroutine
!======================================================================================      
      pure subroutine Mol_ECalc_Improper(iType,iMol,E_Improp)
      use SimParameters
      use ForceField
      use Coords
      use ForceFieldFunctions
      implicit none
      real(dp), intent(out) :: E_Improp
      integer, intent(in) :: iType,iMol
      integer :: iTors
      integer(kind=atomIntType) :: impropType      
      integer :: memb1, memb2, memb3, memb4
      real(dp) :: x12,y12,z12
      real(dp) :: x23,y23,z23
      real(dp) :: x34,y34,z34
      real(dp) :: vx1,vy1,vz1
      real(dp) :: vx2,vy2,vz2
      real(dp) :: vx3,vy3,vz3
      real(dp) :: r1,r3,dot1,dot2
      real(dp) :: Angle      

      E_Improp = 0d0
      do iTors = 1,nTorsional(iType)
        impropType = impropArray(iType,iTors)%impropType
            
        memb1 = impropArray(iType,iTors)%impropMembr(1)
        memb2 = impropArray(iType,iTors)%impropMembr(2)
        memb3 = impropArray(iType,iTors)%impropMembr(3)
        memb4 = impropArray(iType,iTors)%impropMembr(4)
            
        x12 = MolArray(iType)%mol(iMol)%x(memb2) - MolArray(iType)%mol(iMol)%x(memb1)
        y12 = MolArray(iType)%mol(iMol)%y(memb2) - MolArray(iType)%mol(iMol)%y(memb1)
        z12 = MolArray(iType)%mol(iMol)%z(memb2) - MolArray(iType)%mol(iMol)%z(memb1)
            
        x23 = MolArray(iType)%mol(iMol)%x(memb3) - MolArray(iType)%mol(iMol)%x(memb2)
        y23 = MolArray(iType)%mol(iMol)%y(memb3) - MolArray(iType)%mol(iMol)%y(memb2)
        z23 = MolArray(iType)%mol(iMol)%z(memb3) - MolArray(iType)%mol(iMol)%z(memb2)          

        x34 = MolArray(iType)%mol(iMol)%x(memb4) - MolArray(iType)%mol(iMol)%x(memb3)
        y34 = MolArray(iType)%mol(iMol)%y(memb4) - MolArray(iType)%mol(iMol)%y(memb3)
        z34 = MolArray(iType)%mol(iMol)%z(memb4) - MolArray(iType)%mol(iMol)%z(memb3)          

!       Calculate the vector, v1 = <r12 x r23>
        vx1 =   y12*z23 - z12*y23
        vy1 = -(x12*z23 - z12*x23)
        vz1 =   x12*y23 - y12*x23
           
!       Calculate the vector, v2=<r3 x r2>
        vx2 =   y23*z34 - z23*y34
        vy2 = -(x23*z34 - z23*x34)
        vz2 =   x23*y34 - y23*x34   

!       Calculate the vector, v3=<v1 x r2>
        vx3 =   vy1*z23 - vz1*y23
        vy3 = -(vx1*z23 - vz1*x23)
        vz3 =   vx1*y23 - vy1*x23            
            
        r1 = vx1**2 + vy1**2 + vz1**2
        r3 = vx3**2 + vy3**2 + vz3**2
        r1 = sqrt(r1)
        r3 = sqrt(r3)
        dot1 = vx1*vx2 + vy1*vy2 + vz1*vz2
        dot2 = vx2*vx3 + vy2*vy3 + vz2*vz3            
        dot1 = dot1/(r1)
        dot2 = dot2/(r3)   
        angle = atan2(dot2,dot1)            
            
        E_Improp = E_Improp + TorsHarmonic(Angle, impropData(impropType)%a)
      enddo
     
      end subroutine
!======================================================================================      
      pure subroutine NewMol_ECalc_Improper(E_Improp)
      use SimParameters
      use ForceField
      use Coords
      use ForceFieldFunctions
      implicit none
      real(dp), intent(out) :: E_Improp
      integer :: iType,iTors
      integer(kind=atomIntType) :: impropType      
      integer :: memb1, memb2, memb3, memb4
      real(dp) :: x12,y12,z12
      real(dp) :: x23,y23,z23
      real(dp) :: x34,y34,z34
      real(dp) :: vx1,vy1,vz1
      real(dp) :: vx2,vy2,vz2
      real(dp) :: vx3,vy3,vz3
      real(dp) :: r1,r3,dot1,dot2
      real(dp) :: Angle      

      E_Improp = 0d0
      iType = newMol%molType
      do iTors = 1,nTorsional(iType)
        impropType = impropArray(iType,iTors)%impropType
            
        memb1 = impropArray(iType,iTors)%impropMembr(1)
        memb2 = impropArray(iType,iTors)%impropMembr(2)
        memb3 = impropArray(iType,iTors)%impropMembr(3)
        memb4 = impropArray(iType,iTors)%impropMembr(4)
            
        x12 = newMol%x(memb2) - newMol%x(memb1)
        y12 = newMol%y(memb2) - newMol%y(memb1)
        z12 = newMol%z(memb2) - newMol%z(memb1)
            
        x23 = newMol%x(memb3) - newMol%x(memb2)
        y23 = newMol%y(memb3) - newMol%y(memb2)
        z23 = newMol%z(memb3) - newMol%z(memb2) 

        x34 = newMol%x(memb4) - newMol%x(memb3)
        y34 = newMol%y(memb4) - newMol%y(memb3)
        z34 = newMol%z(memb4) - newMol%z(memb3)

!       Calculate the vector, v1 = <r12 x r23>
        vx1 =   y12*z23 - z12*y23
        vy1 = -(x12*z23 - z12*x23)
        vz1 =   x12*y23 - y12*x23
           
!       Calculate the vector, v2=<r3 x r2>
        vx2 =   y23*z34 - z23*y34
        vy2 = -(x23*z34 - z23*x34)
        vz2 =   x23*y34 - y23*x34   

!       Calculate the vector, v3=<v1 x r2>
        vx3 =   vy1*z23 - vz1*y23
        vy3 = -(vx1*z23 - vz1*x23)
        vz3 =   vx1*y23 - vy1*x23            
            
        r1 = vx1**2 + vy1**2 + vz1**2
        r3 = vx3**2 + vy3**2 + vz3**2
        r1 = sqrt(r1)
        r3 = sqrt(r3)
        dot1 = vx1*vx2 + vy1*vy2 + vz1*vz2
        dot2 = vx2*vx3 + vy2*vy3 + vz2*vz3            
        dot1 = dot1/(r1)
        dot2 = dot2/(r3)   
        angle = atan2(dot2,dot1)            
            
        E_Improp = E_Improp + TorsHarmonic(Angle, impropData(impropType)%a)
      enddo
     
      end subroutine
!======================================================================================      
      end module
      
      
