!*********************************************************************************************************************
!     This file contains the energy functions that calculate the energy associated with 
!     torsional bending motion.
!
!     The subroutine's prefix naming scheme implies the following:
!           Detailed - Complete energy calculation intended for use at the beginning and end
!                      of the simulation.  This function is not intended for use mid-simulation.
!             Shift  - Calculates the energy difference for any move that does not result
!                      in molecules being added or removed from the cluster. This function
!                      is only intended for moves where only one molecule is changed
!                      in a given move. 
!              Mol   - Calculates the energy for a molecule already present in the system. For
!                      use in moves such as particle deletion moves. 
!             NewMol - Calculates the energy for a molecule that has been freshly inserted into the system.
!                      Intended for use in Rosenbluth Sampling, Swap In, etc.
!*********************************************************************************************************************
      module TorsionalFunctions
      contains
!======================================================================================      
      subroutine Detailed_ECalc_Torsional(E_T)
      use ParallelVar      
      use SimParameters
      use ForceField
      use Coords
      use ForceFieldFunctions 
      use Constants
      use EnergyTables, only: E_Torsion_T      
      implicit none
      real(dp), intent(inout) :: E_T
      integer ::  iType,iMol,iTors
      integer(kind=atomIntType) :: torsType      
      integer :: memb1, memb2, memb3, memb4
      real(dp) :: x12,y12,z12
      real(dp) :: x23,y23,z23
      real(dp) :: x34,y34,z34
      real(dp) :: vx1,vy1,vz1
      real(dp) :: vx2,vy2,vz2
      real(dp) :: vx3,vy3,vz3
      real(dp) :: r1,r3,dot1,dot2
      real(dp) :: Angle      
      real(dp) :: E_Tors

      E_Tors = 0d0
      E_Torsion_T = 0d0
      do iType = 1, nMolTypes
        do iMol = 1,NPART(iType)
!          do i = 1,nAtoms(iType)
!            write(2,*)  MolArray(iType)%mol(iMol)%x(i),MolArray(iType)%mol(iMol)%y(i),MolArray(iType)%mol(iMol)%z(i)     
!          enddo
          do iTors = 1,nTorsional(iType)
            torsType = torsArray(iType,iTors)%TorsType
            
            memb1 = torsArray(iType,iTors)%torsMembr(1)
            memb2 = torsArray(iType,iTors)%torsMembr(2)
            memb3 = torsArray(iType,iTors)%torsMembr(3)
            memb4 = torsArray(iType,iTors)%torsMembr(4)
            
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

!           Calculate the vector, v3=<v1 x r2> to create an orthonormal framework
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
            E_Tors = E_Tors + Trappe_CosNx(Angle, torsData(torsType)%a)
!            write(2,*) Angle*180d0/pi,Trappe_CosNx(Angle, torsData(torsType)%a), torsData(torsType)%a
          enddo
        enddo
      enddo
!      write(2,*)
      write(nout,*) "Torsional Energy", E_Tors
      E_T = E_T + E_Tors
      E_Torsion_T = E_Tors
      
      end subroutine
!======================================================================================      
!      pure subroutine Shift_ECalc_Torsional(E_Tors,disp)
      pure subroutine Shift_ECalc_Torsional(E_Tors,disp)      
      use SimParameters
      use ForceField
      use Coords
      use ForceFieldFunctions
      implicit none
      real(dp), intent(out) :: E_Tors
      type(Displacement), intent(in) :: disp(:)

      logical :: changed(1:maxAtoms)
      integer :: dispIndx(1:maxAtoms)        
      integer :: i,nDisp,iType,iMol,iTors
      integer :: torsType      
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
      real(dp) :: Angle,eng


      iType = disp(1)%molType
      E_Tors = 0d0      
      if(nTorsional(iType) .eq. 0) then
         return
      endif      
      iMol  = disp(1)%molIndx      
      nDisp = size(disp)      
      changed = .false.
      dispIndx = 0

      do i=1,nDisp
        changed(disp(i)%atmIndx) = .true.
        dispIndx(disp(i)%atmIndx) = i
      enddo
      
      do iTors = 1,nTorsional(iType)
        torsType = torsArray(iType,iTors)%TorsType
        memb(1) = torsArray(iType,iTors)%torsMembr(1)
        memb(2) = torsArray(iType,iTors)%torsMembr(2)
        memb(3) = torsArray(iType,iTors)%torsMembr(3)
        memb(4) = torsArray(iType,iTors)%torsMembr(4)
        if(all(changed([memb(1:4)]) .eqv. .false. )) then
           cycle
        endif

        do i=1,4
          if( changed(memb(i)) ) then
            x_new(i) = disp(dispIndx(memb(i)))%x_new
            y_new(i) = disp(dispIndx(memb(i)))%y_new
            z_new(i) = disp(dispIndx(memb(i)))%z_new
            
            x_old(i) = disp(dispIndx(memb(i)))%x_old
            y_old(i) = disp(dispIndx(memb(i)))%y_old
            z_old(i) = disp(dispIndx(memb(i)))%z_old
          else
            x_old(i) = MolArray(iType)%mol(iMol)%x(memb(i))
            y_old(i) = MolArray(iType)%mol(iMol)%y(memb(i))
            z_old(i) = MolArray(iType)%mol(iMol)%z(memb(i))
            
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

!        Calculate the vector, v1 = <r12 x r23>
        vx1 =   y12*z23 - z12*y23
        vy1 = -(x12*z23 - z12*x23)
        vz1 =   x12*y23 - y12*x23
           
!        Calculate the vector, v2=<r3 x r2>
        vx2 =   y23*z34 - z23*y34
        vy2 = -(x23*z34 - z23*x34)
        vz2 =   x23*y34 - y23*x34   

!        Calculate the vector, v3=<v1 x r2> to create an orthonormal framework
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
            
        eng = Trappe_CosNx(Angle, torsData(torsType)%a)
        E_Tors = E_Tors + eng
        
        x12 = x_old(2) - x_old(1)
        y12 = y_old(2) - y_old(1)
        z12 = z_old(2) - z_old(1)
            
        x23 = x_old(3) - x_old(2)
        y23 = y_old(3) - y_old(2)
        z23 = z_old(3) - z_old(2)

        x34 = x_old(4) - x_old(3)
        y34 = y_old(4) - y_old(3)
        z34 = z_old(4) - z_old(3)


!        Calculate the vector, v1 = <r12 x r23>
        vx1 =   y12*z23 - z12*y23
        vy1 = -(x12*z23 - z12*x23)
        vz1 =   x12*y23 - y12*x23
           
!        Calculate the vector, v2=<r34 x r23>
        vx2 =   y23*z34 - z23*y34
        vy2 = -(x23*z34 - z23*x34)
        vz2 =   x23*y34 - y23*x34   

!        Calculate the vector, v3=<v1 x r23> to create an orthonormal framework
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
            
        eng = Trappe_CosNx(Angle, torsData(torsType)%a)
        E_Tors = E_Tors - eng
      enddo
      
      end subroutine
!======================================================================================      
      pure subroutine Mol_ECalc_Torsional(iType,iMol,E_Tors)
      use SimParameters
      use ForceField
      use Coords
      use ForceFieldFunctions
      implicit none
      real(dp), intent(out) :: E_Tors
      integer, intent(in) :: iType,iMol
      integer :: iTors
      integer(kind=atomIntType) :: torsType      
      integer :: memb1, memb2, memb3, memb4
      real(dp) :: x12,y12,z12
      real(dp) :: x23,y23,z23
      real(dp) :: x34,y34,z34
      real(dp) :: vx1,vy1,vz1
      real(dp) :: vx2,vy2,vz2
      real(dp) :: vx3,vy3,vz3
      real(dp) :: r1,r3,dot1,dot2
      real(dp) :: Angle      

      if(nTorsional(iType) .eq. 0) then
        return
      endif
      
      E_Tors = 0d0
      do iTors = 1,nTorsional(iType)
        torsType = torsArray(iType,iTors)%TorsType
            
        memb1 = torsArray(iType,iTors)%torsMembr(1)
        memb2 = torsArray(iType,iTors)%torsMembr(2)
        memb3 = torsArray(iType,iTors)%torsMembr(3)
        memb4 = torsArray(iType,iTors)%torsMembr(4)
            
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
            
        E_Tors = E_Tors + Trappe_CosNx(Angle, torsData(torsType)%a)
      enddo
     
      end subroutine
!======================================================================================      
      pure subroutine NewMol_ECalc_Torsional(E_Tors)
      use SimParameters
      use ForceField
      use Coords
      use ForceFieldFunctions
      implicit none
      real(dp), intent(out) :: E_Tors
      integer :: iType,iTors
      integer(kind=atomIntType) :: torsType      
      integer :: memb1, memb2, memb3, memb4
      real(dp) :: x12,y12,z12
      real(dp) :: x23,y23,z23
      real(dp) :: x34,y34,z34
      real(dp) :: vx1,vy1,vz1
      real(dp) :: vx2,vy2,vz2
      real(dp) :: vx3,vy3,vz3
      real(dp) :: r1,r3,dot1,dot2
      real(dp) :: Angle      

      E_Tors = 0d0
      iType = newMol%molType
      if(nTorsional(iType) .eq. 0) then
        return
      endif
      
      do iTors = 1,nTorsional(iType)
        torsType = torsArray(iType,iTors)%TorsType
            
        memb1 = torsArray(iType,iTors)%torsMembr(1)
        memb2 = torsArray(iType,iTors)%torsMembr(2)
        memb3 = torsArray(iType,iTors)%torsMembr(3)
        memb4 = torsArray(iType,iTors)%torsMembr(4)
            
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
            
        E_Tors = E_Tors + Trappe_CosNx(Angle, torsData(torsType)%a)
      enddo
     
      end subroutine
!======================================================================================      
      end module
      
      
