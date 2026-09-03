!*********************************************************************************************************************
!     This file contains the pressure functions that work for Wang–Frenkel w/ Columbic style forcefields
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
      module Pressure_Mpipi
      use VarPrecision
!      use InterEnergy_Mpipi, only: lj_cut, lj_Cut_sq, q_cut, q_cut_sq

      contains
!======================================================================================      
      pure function WF_Press_Func(r_sq, epAlpha, sig_sq, Mu) result(WF)
      implicit none
      real(dp), intent(in) :: r_sq, epAlpha, sig_sq
      integer, intent(in) :: Mu
      real(dp) :: WF, x, y, z, k  

      x = sig_sq/r_sq
      if(Mu .eq. 2) then
        x = x*x
        k = 81.0_dp  ! 9^2
      elseif(Mu .eq. 3) then
        x = x*x*x
        k = 729.0_dp  ! 9^3
      else
        x = x**Mu
        k = 9.0_dp ** Mu
      endif
      y = k*x
      z = y-1E0_dp

      WF = 2E0_dp * epAlpha * Mu * (x*z*z + 2E0_dp*y*(x-1E0_dp)*z)

      end function
!======================================================================================      
      subroutine Detailed_PCalc_Mpipi(P_T)
!      use ParallelVar
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, rcutElecsq
      use Coords, only: MolArray
      use SimParameters, only: nMolTypes, NPART, kapa
      use ParallelVar, only: nout
      implicit none
      real(dp), intent(out) :: P_T

      integer :: iType,jType,iMol,jMol,iAtom,jAtom,Mu
      integer(kind=atomIntType) :: atmType1,atmType2      
      integer :: iIndx, jIndx, jMolMin
      real(dp) :: r_sq, r,  rx, ry, rz
      real(dp) :: epAlpha,sig_sq,q,rcutNBsq
      logical :: q_Nzero
      real(dp) :: WF, Ele
      real(dp) :: P_Ele, P_WF    



      P_WF = 0E0_dp
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
                 epAlpha = epAlpha_tab(atmType1,atmType2)
                 q_Nzero = q_Nonzero(atmType1,atmType2)
                 sig_sq = sigsq_tab(atmType1,atmType2)    
                 rcutNBsq = 9d0 * sig_sq
                 Mu = Mu_tab(atmType1,atmType2)      
                 rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
                 ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
                 rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom) 
                 r_sq = rx**2 + ry**2 + rz**2
                 if(r_sq .lt. rcutNBsq) then
                   WF = WF_Press_Func(r_sq, epAlpha, sig_sq, Mu)             
                   P_WF = P_WF + WF
                 endif  
                 Ele = 0E0_dp
                 if (q_Nzero) then
                  if (r_sq .lt. rcutElecsq) then
                   q = q_tab(atmType1,atmType2)
                   r = sqrt(r_sq)
                   Ele = q * exp(-kapa * r) * (kapa + (1E0_dp/r))
                  endif
                 endif
                 P_Ele = P_Ele + Ele
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      
      write(nout,*) "Wang–Frenkel Pressure:", P_WF
      write(nout,*) "Eletrostatic Pressure:", P_Ele
      P_T = P_Ele + P_WF    
      write(nout,*) "Total Pressure:", P_T
      end subroutine
!======================================================================================      
      subroutine Shift_PCalc_Mpipi(P_Trial, disp)
      use Coords, only: Displacement,  molArray
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, rcutElecsq
      use SimParameters, only: nMolTypes, NPART, kapa
      use PairStorage, only: DistArrayNew
      implicit none
      
      type(Displacement), intent(in) :: disp(:)  
      real(dp), intent(out) :: P_Trial
      
      integer :: iType, jType, iMol, jMol, iAtom, jAtom, Mu
      integer(kind=atomIntType) :: atmType1, atmType2, iIndx
      integer :: sizeDisp, iDisp
!      integer, pointer :: oldIndx 
      real(dp) :: r_new, r_old, rx, ry, rz
      real(dp) :: epAlpha,sig_sq,q,rcutNBsq
      logical :: q_Nzero
      real(dp) :: WF, Ele, P_Old
      real(dp) :: P_Ele, P_WF


      P_WF = 0E0_dp
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
            epAlpha = epAlpha_tab(atmType2, atmType1)
            q_Nzero = q_Nonzero(atmType2, atmType1)
            sig_sq = sigsq_tab(atmType2,atmType1)
            rcutNBsq = 9d0 * sig_sq
            Mu = Mu_tab(atmType2, atmType1)
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
              WF = 0E0_dp
!              if(ep .ne. 0E0_dp) then
                if(r_new .lt. rcutNBsq) then
                  WF = WF_Press_Func(r_new, epAlpha, sig_sq, Mu)                 
                endif
                
                if(r_old .lt. rcutNBsq) then
                  WF = WF - WF_Press_Func(r_old, epAlpha, sig_sq, Mu)                 
                endif
!              endif
!             Check to see if there is a non-zero Electrostatic parmaeter. If so calculate
!             the electrostatic energy              
              Ele = 0E0_dp
              if(q_Nzero) then
                q = q_tab(atmType2, atmType1)
               if (r_new .lt. rcutElecsq) then
                r_new = sqrt(r_new)
                Ele = q * exp(-kapa * r_new) * (kapa + (1E0_dp/r_new))
               endif
               if (r_old .lt. rcutElecsq) then
                r_old = sqrt(r_old)
                Ele = Ele - q * exp(-kapa * r_old) * (kapa + (1E0_dp/r_old))
               endif
              endif
              P_Trial = P_Trial + Ele + WF
            enddo
          enddo
        enddo
      enddo

      
      end subroutine
!======================================================================================      
      pure subroutine SwapOut_PCalc_Mpipi(iType, iMol, P_Trial)
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, rcutElecsq
      use Coords, only: Displacement,  molArray
      use SimParameters, only: nMolTypes, NPART, kapa
      implicit none
      integer, intent(in) :: iType, iMol     
      real(dp), intent(out) :: P_Trial

      integer :: iAtom, iIndx, jType, jIndx, jMol, jAtom, Mu
      integer  :: atmType1, atmType2
      real(dp) :: epAlpha,sig_sq,q,rcutNBsq
      logical :: q_Nzero
      real(dp) :: WF, Ele
      real(dp) :: r_sq, r,  rx, ry, rz

      P_Trial = 0E0_dp
      iIndx = MolArray(iType)%mol(iMol)%indx

      do iAtom = 1,nAtoms(iType)
        atmType1 = atomArray(iType, iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType,jAtom)
            epAlpha = epAlpha_tab(atmType2, atmType1)
            q_Nzero = q_Nonzero(atmType2, atmType1)
            sig_sq = sigsq_tab(atmType2,atmType1)
            rcutNBsq = 9d0 * sig_sq
            Mu = Mu_tab(atmType2, atmType1)
            do jMol=1, NPART(jType)
              jIndx = MolArray(jType)%mol(jMol)%indx  
              if(iIndx .eq. jIndx) then
                cycle
              endif

              rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r_sq = rx*rx + ry*ry + rz*rz
              WF = 0E0_dp
!              if(ep .ne. 0E0_dp) then
                if(r_sq .lt. rcutNBsq) then
                  WF = WF_Press_Func(r_sq, epAlpha, sig_sq, Mu)           
                endif
!              endif

              Ele = 0E0_dp
              if(q_Nzero) then
               if (r_sq .lt. rcutElecsq) then
                q = q_tab(atmType2, atmType1)
                r = sqrt(r_sq)
                Ele = q * exp(-kapa * r) * (kapa + (1E0_dp/r))            
               endif
              endif
              P_Trial = P_Trial + Ele + WF
            enddo
          enddo
        enddo
      enddo

  
      
      end subroutine
!======================================================================================      
      subroutine SwapIn_PCalc_Mpipi(P_Trial)
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, rcutElecsq
      use Coords, only: Displacement, molArray, newmol
      use SimParameters, only: nMolTypes, NPART, kapa
      implicit none
      real(dp), intent(out) :: P_Trial
      
      integer :: iType, iMol, iAtom, iIndx, jType, jMol, jAtom, Mu
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: r_sq, r,  rx, ry, rz
      real(dp) :: epAlpha,sig_sq,q,rcutNBsq
      logical :: q_Nzero
      real(dp) :: WF, Ele

      real(dp) :: P_Ele, P_WF

      P_WF = 0E0_dp
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
            epAlpha = epAlpha_tab(atmType2, atmType1)
            q_Nzero = q_Nonzero(atmType2, atmType1)
            sig_sq = sigsq_tab(atmType2,atmType1)
            rcutNBsq = 9d0 * sig_sq
            Mu = Mu_tab(atmType2, atmType1)
            do jMol = 1, NPART(jType)
              rx = newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r_sq = rx*rx + ry*ry + rz*rz

!              if(ep .ne. 0E0) then
                if(r_sq .lt. rcutNBsq) then
                  WF = WF_Press_Func(r_sq, epAlpha, sig_sq, Mu)  
                  P_WF = P_WF + WF
                endif
!              endif
              if(q_Nzero) then
               if (r_sq .lt. rcutElecsq) then
                q = q_tab(atmType2, atmType1)
                r = sqrt(r_sq)
                Ele = q * exp(-kapa * r) * (kapa + (1E0_dp/r))
                P_Ele = P_Ele + Ele
               endif
              endif
            enddo
          enddo
        enddo
      enddo

      P_Trial = P_WF + P_Ele
      
      end subroutine 
!======================================================================================
      end module
      
       
