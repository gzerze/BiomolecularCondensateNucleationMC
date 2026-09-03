!*********************************************************************************************************************
!     This file contains the pressure functions that work for Hydrophobicity scale w/ electrostatic screening style forcefields
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
      module Pressure_HPS_single
      use VarPrecision

      contains
!======================================================================================      
      subroutine Detailed_PCalc_HPS_single(P_T)
!      use ParallelVar
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_HPS_single, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
      use Coords, only: MolArray
      use SimParameters, only: nMolTypes, NPART, distCriteria, kapa
      use PairStorage, only: rPair, distStorage
      use ParallelVar, only: nout
      implicit none
      real(dp), intent(out) :: P_T

      integer :: iType,jType,iMol,jMol,iAtom,jAtom
      integer(kind=atomIntType) :: atmType1,atmType2      
      integer :: iIndx, jIndx, globIndx1, globIndx2, jMolMin
      real(dp) :: r_sq, r
      real(dp) :: eps,sig_sq,q,lambda
      logical :: q_Nzero
      real(dp) :: HPS, LJ, Ele, x,rcutNBsq
      real(dp) :: P_Ele, P_HPS    



      P_HPS = 0E0_dp
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
               globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
               do jAtom = 1,nAtoms(jType)        
                 atmType2 = atomArray(jType,jAtom)
                 eps = eps_tab(atmType1,atmType2)
                 q_Nzero = q_Nonzero(atmType1,atmType2)
                 sig_sq = sigsq_tab(atmType1,atmType2) 
                 rcutNBsq = cutoffNBsq_tab(atmType1,atmType2)
                 lambda = lambda_tab(atmType1,atmType2)       
                 globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
                 r_sq = rPair(globIndx1, globIndx2)%p%r_sq
!                 if(r_sq .gt. lj_Cut_sq) then
!                   cycle
!                 endif 
!                 if(r_sq .lt. lj_cut_sq) then
                  HPS = 0E0_dp
                  if (r_sq .lt. rcutNBsq) then
                    x = sig_sq/r_sq
                    x = x * x * x
                    LJ = 48.0_dp * eps * x * (x - 0.5_dp)    
                    HPS = lambda * LJ         
                  endif
                   P_HPS = P_HPS + HPS
!                 endif 
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
      
      write(nout,*) "Hydrophobicity scale Pressure:", P_HPS
      write(nout,*) "Eletrostatic Pressure:", P_Ele
      P_T = P_Ele + P_HPS    
      write(nout,*) "Total Pressure:", P_T
      end subroutine
!======================================================================================      
      subroutine Shift_PCalc_HPS_single(P_Trial, disp)
      use Coords, only: Displacement,  atomIndicies, molArray
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_HPS_single, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
      use SimParameters, only: distCriteria, kapa
      use PairStorage, only: distStorage, rPair, rPairNew, newDist, DistArrayNew, nNewDist, oldIndxArray
      implicit none
      
      type(Displacement), intent(in) :: disp(:)  
      real(dp), intent(out) :: P_Trial
      
      integer :: iType, jType, iMol, jMol, iAtom, jAtom, iPair
      integer(kind=atomIntType) :: atmType1, atmType2, iIndx, jIndx
      integer :: sizeDisp
!      integer, pointer :: oldIndx 
      integer :: globIndx1, globIndx2
      real(dp) :: r_new, r_new_sq, r, r_sq
      real(dp) :: eps,sig_sq,q,rcutNBsq,lambda
      logical :: q_Nzero
      real(dp) :: HPS, LJ, Ele, x, P_New, P_PairOld, P_Old
      real(dp) :: P_Ele, P_HPS


      P_HPS = 0E0_dp
      P_Ele = 0E0_dp
      P_Trial = 0E0_dp
      P_Old = 0E0_dp
      iType = disp(1)%molType
      iMol = disp(1)%molIndx
      iIndx = MolArray(iType)%mol(iMol)%indx
      
!      !This section calculates the Intermolecular interaction between the atoms that
!      !have been modified in this trial move with the atoms that have remained stationary

      do iPair = 1, nNewDist

        globIndx1 = newDist(iPair)%indx1
        globIndx2 = newDist(iPair)%indx2
        if(.not. rPair(globIndx1, globIndx2)%p%usePair) then
          cycle
        endif
        jType = atomIndicies(globIndx2)%nType
        jMol  = atomIndicies(globIndx2)%nMol
        jIndx = MolArray(jType)%mol(jMol)%indx
        jAtom = atomIndicies(globIndx2)%nAtom
        iAtom = atomIndicies(globIndx1)%nAtom
       
        atmType1 = atomArray(iType,iAtom)
        atmType2 = atomArray(jType,jAtom)

        eps = eps_tab(atmType2, atmType1)
        q_Nzero = q_Nonzero(atmType2, atmType1)
        sig_sq = sigsq_tab(atmType2,atmType1)
        rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
        lambda = lambda_tab(atmType2, atmType1)
        HPS = 0E0_dp
!        if(ep .ne. 0E0_dp) then
        r_new_sq = rPairNew(globIndx1, globIndx2)%p%r_sq
!         if(r_new_sq .lt. lj_cut_sq) then
        if (r_new_sq .lt. rcutNBsq) then
          x = sig_sq/r_new_sq
          x = x * x * x
          LJ = 48.0_dp * eps * x * (x - 0.5_dp)    
          HPS = lambda * LJ    
        endif
!          endif
        r_sq = rPair(globIndx1, globIndx2)%p%r_sq  
!          if(r_sq .lt. lj_cut_sq) then
        if (r_sq .lt. rcutNBsq) then
          x = sig_sq/r_sq
          x = x * x * x
          LJ = 48.0_dp * eps * x * (x - 0.5_dp)    
          HPS = HPS - lambda * LJ 
        endif
!          endif         
!        endif
        Ele = 0E0_dp
        if(q_Nzero) then
          q = q_tab(atmType2, atmType1)
         if (r_new_sq .lt. rcutElecsq) then
          r_new = rPairNew(globIndx1, globIndx2)%p%r
          Ele = Ele + q * exp(-kapa * r_new) * (kapa + (1E0_dp/r_new))
         endif
         if (r_sq .lt. rcutElecsq) then
          r = rPair(globIndx1, globIndx2)%p%r
          Ele = Ele - q * exp(-kapa * r) * (kapa + (1E0_dp/r))             
         endif
        endif
        P_Trial = P_Trial + Ele + HPS
        if(abs(P_Trial) > 1E6) then
          write(*,*) r_new, r, r_sq, HPS, Ele, q, epAlpha, sig_sq, rPair(globIndx1, globIndx2)%p%storeRValue
        endif
      enddo    
      
      end subroutine
!======================================================================================      
      pure subroutine SwapOut_PCalc_HPS_single(iType, iMol, P_Trial)
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_HPS_single, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
      use Coords, only: Displacement,  atomIndicies, molArray
      use SimParameters, only: distCriteria, nMolTypes, NPART, kapa
      use PairStorage, only: distStorage, rPair
      implicit none
      integer, intent(in) :: iType, iMol     
      real(dp), intent(out) :: P_Trial

      integer :: iAtom, iIndx, jType, jIndx, jMol, jAtom
      integer  :: globIndx1, globIndx2
      integer  :: atmType1, atmType2
      real(dp) :: eps,sig_sq,q,rcutNBsq,lambda
      logical :: q_Nzero
      real(dp) :: HPS, LJ, Ele, x
      real(dp) :: r_sq, r

      P_Trial = 0E0_dp
      iIndx = MolArray(iType)%mol(iMol)%indx

      do iAtom = 1,nAtoms(iType)
        atmType1 = atomArray(iType, iAtom)
        globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom) 
        do jType = 1, nMolTypes
          do jMol=1, NPART(jType)
            jIndx = MolArray(jType)%mol(jMol)%indx  
            if(iIndx .eq. jIndx) then
              cycle
            endif
            do jAtom = 1,nAtoms(jType)        
              globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
              if(.not. rPair(globIndx1, globIndx2)%p%usePair) then
                cycle
              endif 
              atmType2 = atomArray(jType,jAtom)
              eps = eps_tab(atmType2, atmType1)
              q_Nzero = q_Nonzero(atmType2, atmType1)
              sig_sq = sigsq_tab(atmType2,atmType1)
              rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
              lambda = lambda_tab(atmType2, atmType1)
              HPS = 0E0_dp
!              r_sq = rPair(globIndx1, globIndx2)%p%r_sq
!              if(ep .ne. 0E0_dp) then
              r_sq = rPair(globIndx1, globIndx2)%p%r_sq
!                if(r_sq .lt. lj_Cut_sq) then
!                  sig_sq = sig_tab(atmType2, atmType1)
              if (r_sq .lt. rcutNBsq) then
                x = sig_sq/r_sq
                x = x * x * x
                LJ = 48.0_dp * eps * x * (x - 0.5_dp)    
                HPS = lambda * LJ  
              endif
!                endif
!              endif

              Ele = 0E0_dp
              if(q_Nzero) then
               if (r_sq .lt. rcutElecsq) then
                q = q_tab(atmType2, atmType1)
                r = rPair(globIndx1, globIndx2)%p%r
                Ele = q * exp(-kapa * r) * (kapa + (1E0_dp/r))  
               endif          
              endif
              P_Trial = P_Trial + Ele + HPS
            enddo
          enddo
        enddo
      enddo

  
      
      end subroutine
!======================================================================================      
      subroutine SwapIn_PCalc_HPS_single(P_Trial)
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_HPS_single, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
      use Coords, only: Displacement,  atomIndicies, molArray, newmol
      use SimParameters, only: distCriteria, nMolTypes, NPART, kapa
      use PairStorage, only: distStorage, rPair, rPairNew, newDist, nNewDist
      implicit none
      real(dp), intent(out) :: P_Trial
      
      integer :: iPair
      integer :: iType, iMol, iAtom, iIndx, jType, jIndx, jMol, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      integer :: globIndx1, globIndx2
      real(dp) :: r_sq, r
      real(dp) :: eps,sig_sq,q,rcutNBsq,lambda
      logical :: q_Nzero
      real(dp) :: HPS, LJ, Ele, x

      real(dp) :: P_Ele, P_HPS

      P_HPS = 0E0_dp
      P_Ele = 0E0_dp
      P_Trial = 0E0_dp

      iType = newMol%molType
      iMol = NPART(iType)+1
      iIndx = molArray(iType)%mol(iMol)%indx
      do iPair = 1, nNewDist
        globIndx1 = newDist(iPair)%indx1
        globIndx2 = newDist(iPair)%indx2
        if(.not. rPair(globIndx1, globIndx2)%p%usePair) then
          cycle
        endif
        jType = atomIndicies(globIndx2)%nType
        jMol = atomIndicies(globIndx2)%nMol
        jIndx = MolArray(jType)%mol(jMol)%indx
        if(iIndx .ne. jIndx) then
          iAtom = atomIndicies(globIndx1)%nAtom
          jAtom = atomIndicies(globIndx2)%nAtom

          atmType1 = atomArray(iType, iAtom)
          atmType2 = atomArray(jType, jAtom)

          eps = eps_tab(atmType2, atmType1)
          q_Nzero = q_Nonzero(atmType2, atmType1)
          sig_sq = sigsq_tab(atmType2,atmType1)
          rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
          lambda = lambda_tab(atmType2, atmType1)
!          if(ep .ne. 0E0_dp) then
            r_sq = newDist(iPair)%r_sq
            if(r_sq .lt. rcutNBsq) then
!              sig_sq = sig_tab(atmType2,atmType1)
              x = sig_sq/r_sq
              x = x * x * x
              LJ = 48.0_dp * eps * x * (x - 0.5_dp)    
              HPS = lambda * LJ  
              P_HPS = P_HPS + HPS
            endif
!          endif

          if(q_Nzero) then
           if (r_sq .lt. rcutElecsq) then
            q = q_tab(atmType2, atmType1)
            r = newDist(iPair)%r
            Ele = q * exp(-kapa * r) * (kapa + (1E0_dp/r))
            P_Ele = P_Ele + Ele
           endif
          endif
        endif
      enddo

      P_Trial = P_HPS + P_Ele
      
      end subroutine 
!======================================================================================
      end module
      
       
