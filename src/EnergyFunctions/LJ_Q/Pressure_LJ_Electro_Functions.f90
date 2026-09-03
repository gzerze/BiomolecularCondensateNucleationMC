!*********************************************************************************************************************
!     This file contains the pressure functions that work for Lennard-Jones w/ Columbic style forcefields
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
      module Pressure_LJ_Electro
      use VarPrecision
      use InterEnergy_LJ_Electro, only: lj_cut, lj_Cut_sq, q_cut, q_cut_sq

      contains
!======================================================================================      
      pure function LJ_Func(r_sq, ep, sig) result(LJ)
      implicit none
      real(dp), intent(in) :: r_sq, ep, sig
      real(dp) :: LJ  
 
      LJ = (sig/r_sq)
      LJ = LJ * LJ * LJ
      LJ = 12d0 * ep * LJ * (LJ - 0.5E0_dp)  

      end function
!======================================================================================      
      pure function Ele_Func(r, q) result(Ele)
      implicit none
      real(dp), intent(in) :: r, q
      real(dp) :: Ele
 
!      r = sqrt(r_sq)
      Ele = q/r

      end function
!======================================================================================      
      subroutine Detailed_PCalc_LJ_Q(P_T)
!      use ParallelVar
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_LJ_Q, only: ep_tab, q_tab, sig_tab
      use Coords, only: MolArray
      use SimParameters, only: nMolTypes, NPART
      use ParallelVar, only: nout
      implicit none
      real(dp), intent(out) :: P_T

      integer :: iType,jType,iMol,jMol,iAtom,jAtom
      integer(kind=atomIntType) :: atmType1,atmType2      
      integer :: iIndx, jIndx, jMolMin
      real(dp) :: r_sq, r,  rx, ry, rz
      real(dp) :: ep, sig_sq, q
      real(dp) :: LJ, Ele
      real(dp) :: P_Ele, P_LJ    



      P_LJ = 0E0_dp
      P_Ele = 0E0_dp
      P_T = 0E0_dp
      do iType = 1,nMolTypes
        do jType = iType, nMolTypes
          do iMol=1,NPART(iType)
           if(iType .eq. jType) then
             jMolMin = iMol+1
           else
             jMolMin = 1        
           endif
           do jMol = jMolMin,NPART(jType)
             iIndx = MolArray(iType)%mol(iMol)%indx
             jIndx = MolArray(jType)%mol(jMol)%indx  
             do iAtom = 1,nAtoms(iType)
               atmType1 = atomArray(iType,iAtom)
               do jAtom = 1,nAtoms(jType)        
                 atmType2 = atomArray(jType,jAtom)
                 ep = ep_tab(atmType1,atmType2)
                 q = q_tab(atmType1,atmType2)
                 sig_sq = sig_tab(atmType1,atmType2)          
                 rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
                 ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
                 rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom) 
                 r_sq = rx**2 + ry**2 + rz**2
                 if(r_sq .lt. lj_cut_sq) then
                   LJ = LJ_Func(r_sq, ep, sig_sq)             
                   P_LJ = P_LJ + LJ
                 endif  
                 r = sqrt(r_sq)
                 Ele = q / r
                 P_Ele = P_Ele + Ele
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      
      write(nout,*) "Lennard-Jones Pressure:", P_LJ
      write(nout,*) "Eletrostatic Pressure:", P_Ele
      P_T = P_Ele + P_LJ    
      write(nout,*) "Total Pressure:", P_T
      end subroutine
!======================================================================================      
      subroutine Shift_PCalc_LJ_Q(P_Trial, disp)
      use Coords, only: Displacement,  molArray
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_LJ_Q, only: ep_tab, q_tab, sig_tab
      use SimParameters, only: nMolTypes, NPART
      use PairStorage, only: DistArrayNew
      implicit none
      
      type(Displacement), intent(in) :: disp(:)  
      real(dp), intent(out) :: P_Trial
      
      integer :: iType, jType, iMol, jMol, iAtom, jAtom
      integer(kind=atomIntType) :: atmType1, atmType2, iIndx
      integer :: sizeDisp, iDisp
!      integer, pointer :: oldIndx 
      real(dp) :: r_new, r_old, rx, ry, rz
      real(dp) :: ep, sig_sq, q
      real(dp) :: LJ, Ele, P_Old
      real(dp) :: P_Ele, P_LJ


      P_LJ = 0E0_dp
      P_Ele = 0E0_dp
      P_Trial = 0E0_dp
      P_Old = 0E0_dp
      iType = disp(1)%molType
      iMol = disp(1)%molIndx
      iIndx = MolArray(iType)%mol(iMol)%indx
      sizeDisp = size(disp)
!      !This section calculates the Intermolecular interaction between the atoms that
!      !have been modified in this trial move with the atoms that have remained stationary

!      write(*,*) iType
      do iDisp = 1, sizeDisp
        iAtom = disp(iDisp)%atmIndx
!        write(*,*) iAtom, iDisp
        atmType1 = atomArray(iType,iAtom)
!        write(*,*) iAtom, atmType1
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType,jAtom)
            ep = ep_tab(atmType2, atmType1)
            q = q_tab(atmType2, atmType1)
            sig_sq = sig_tab(atmType2,atmType1)
            do jMol=1,NPART(jType)
              if(iType .eq. jType) then
                if(iMol .eq. jMol) then
                  cycle
                endif
              endif  

!             Distance for the New position
              rx = disp(iDisp)%x_new - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = disp(iDisp)%y_new - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = disp(iDisp)%z_new - MolArray(jType)%mol(jMol)%z(jAtom)
              r_new = rx*rx + ry*ry + rz*rz

!             Distance for the Old position
              rx = disp(iDisp)%x_old - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = disp(iDisp)%y_old - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = disp(iDisp)%z_old - MolArray(jType)%mol(jMol)%z(jAtom)
              r_old = rx*rx + ry*ry + rz*rz              


!             Check to see if there is a non-zero Lennard-Jones parmaeter. If so calculate
!             the Lennard-Jones energy           
              LJ = 0E0_dp
              if(ep .ne. 0E0_dp) then
                if(r_new .lt. lj_cut_sq) then
                  LJ = LJ_Func(r_new, ep, sig_sq)                 
                endif
                
                if(r_old .lt. lj_cut_sq) then
                  LJ = LJ - LJ_Func(r_old, ep, sig_sq)                 
                endif
              endif
!             Check to see if there is a non-zero Electrostatic parmaeter. If so calculate
!             the electrostatic energy              
              Ele = 0E0_dp
              if(q .ne. 0E0_dp) then
                r_new = sqrt(r_new)
                Ele = q / r_new
                
                r_old = sqrt(r_old)
                Ele = Ele - q / r_old
              endif
              P_Trial = P_Trial + Ele + LJ
            enddo
          enddo
        enddo
      enddo

      
      end subroutine
!======================================================================================      
      pure subroutine SwapOut_PCalc_LJ_Q(iType, iMol, P_Trial)
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_LJ_Q, only: ep_tab, q_tab, sig_tab
      use Coords, only: Displacement,  molArray
      use SimParameters, only: nMolTypes, NPART
      implicit none
      integer, intent(in) :: iType, iMol     
      real(dp), intent(out) :: P_Trial

      integer :: iAtom, iIndx, jType, jIndx, jMol, jAtom
      integer  :: atmType1, atmType2
      real(dp) :: ep, sig_sq, q
      real(dp) :: Ele, LJ
      real(dp) :: r_sq, r,  rx, ry, rz

      P_Trial = 0E0_dp
      iIndx = MolArray(iType)%mol(iMol)%indx

      do iAtom = 1,nAtoms(iType)
        atmType1 = atomArray(iType, iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType,jAtom)
            sig_sq = sig_tab(atmType2, atmType1)
            ep = ep_tab(atmType2, atmType1)
            q  = q_tab(atmType2, atmType1)
            do jMol=1, NPART(jType)
              jIndx = MolArray(jType)%mol(jMol)%indx  
              if(iIndx .eq. jIndx) then
                cycle
              endif

              rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r_sq = rx*rx + ry*ry + rz*rz
              LJ = 0E0_dp
              if(ep .ne. 0E0_dp) then
                if(r_sq .lt. lj_Cut_sq) then
                  LJ = LJ_Func(r_sq, ep, sig_sq)             
                endif
              endif

              Ele = 0E0_dp
              if(q .ne. 0E0_dp) then
                r = sqrt(r_sq)
                Ele = q/r            
              endif
              P_Trial = P_Trial + Ele + LJ
            enddo
          enddo
        enddo
      enddo

  
      
      end subroutine
!======================================================================================      
      subroutine SwapIn_PCalc_LJ_Q(P_Trial)
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_LJ_Q, only: ep_tab, q_tab, sig_tab
      use Coords, only: Displacement, molArray, newmol
      use SimParameters, only: nMolTypes, NPART
      implicit none
      real(dp), intent(out) :: P_Trial
      
      integer :: iType, iMol, iAtom, iIndx, jType, jMol, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: r_sq, r,  rx, ry, rz
      real(dp) :: ep, sig_sq, q
      real(dp) :: LJ, Ele

      real(dp) :: P_Ele, P_LJ

      P_LJ = 0E0_dp
      P_Ele = 0E0_dp
      P_Trial = 0E0_dp

      iType = newMol%molType
      iMol = NPART(iType)+1
      iIndx = molArray(iType)%mol(iMol)%indx
      do iAtom = 1, nAtoms(newMol%molType)
        atmType1 = atomArray(newMol%molType,iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1, nAtoms(jType)        
            atmType2 = atomArray(jType,jAtom)
            ep = ep_tab(atmType2,atmType1)
            q = q_tab(atmType2,atmType1)
            sig_sq = sig_tab(atmType2,atmType1)
            do jMol = 1, NPART(jType)
              rx = newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r_sq = rx*rx + ry*ry + rz*rz

              if(ep .ne. 0E0) then
                if(r_sq .lt. lj_Cut_sq) then
                  LJ = LJ_Func(r_sq, ep, sig_sq)  
                  P_LJ = P_LJ + LJ
                endif
              endif
              if(q .ne. 0E0) then
                r = sqrt(r_sq)
                Ele = q / r
                P_Ele = P_Ele + Ele
              endif
            enddo
          enddo
        enddo
      enddo

      P_Trial = P_LJ + P_Ele
      
      end subroutine 
!======================================================================================
      end module
      
       
