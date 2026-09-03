!   HPS_piecewise (HPS-KR, KH, FB, TSCL-M2, Urry: Piecewise HPS LJ + Debye-Hückel
!      - u_LJ(r)           = 4 * ε * [ (σ / r)**12 − (σ / r)**6 ]
!      - U_HPS(r) =
!            u_LJ(r) + ε * (1 − λ)       , if r <= 2**(1/6) * σ
!            λ * u_LJ(r)                 , if r >  2**(1/6) * σ
!      - U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!*********************************************************************************************************************
      module Rosenbluth_Functions_HPS_piecewise
      use VarPrecision
      contains
!====================================================================================== 
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which regrows an entire molecule for the given trial.
      pure subroutine Rosen_Mol_New_HPS_piecewise(nRosen, nType, included,  E_Trial, overlap)
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
      use Coords, only: rosenTrial, MolArray
      use SimParameters, only: nMolTypes, NPART, kapa
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nRosen      
      
      logical, intent(out) :: overlap
      real(dp), intent(out) :: E_Trial
      
      integer :: iAtom, jType, jIndx, jMol, jAtom
      integer(kind = atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: eps,sig_sq,q,lambda
      logical :: q_Nzero
      real(dp) :: HPS, LJ, Ele, x
      real(dp) :: E_Ele,E_HPS, rcutNBsq, rmax
      real(dp) :: rmin_ij

      E_Trial = 0E0
      overlap = .false.
      
      E_HPS = 0E0
      E_Ele = 0E0      

      do iAtom = 1,nAtoms(nType)
        atmType1 = atomArray(nType, iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType, jAtom)
            eps = eps_tab(atmType2, atmType1)
            q_Nzero = q_Nonzero(atmType2, atmType1)
            sig_sq = sigsq_tab(atmType2, atmType1)  
            rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
            lambda = lambda_tab(atmType2, atmType1)   
            rmin_ij = r_min_tab(atmType2, atmType1)
            do jMol = 1,NPART(jType)
              jIndx = molArray(jType)%mol(jMol)%indx              
              if(included(jIndx) .eqv. .false.) then
                cycle
              endif            

              rx = abs(rosenTrial(nRosen)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom))
              ry = abs(rosenTrial(nRosen)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom))
              rz = abs(rosenTrial(nRosen)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom))
              rmax = max(rx, ry, rz) 
              if (rmax .gt. rcutElec) cycle
              if (.not. q_Nzero) then
                rmax = rmax * rmax
                if (rmax .gt. rcutNBsq) cycle
              endif
              r = rx*rx + ry*ry + rz*rz
              if (r .gt. rcutElecsq) cycle
              if(r .lt. rmin_ij) then
                E_Trial = huge(dp)
                overlap = .true.
                return
              endif

              HPS = 0E0
              if (r .lt. rcutNBsq) then
                x = sig_sq/r
                x = x * x * x
                LJ = 4.0_dp * eps * x * (x - 1.0_dp)
                if (x .lt. 0.5_dp) then
                  HPS = lambda * LJ
                else
                  HPS = LJ + eps * (1.0_dp - lambda)
                endif
              endif
              E_HPS = E_HPS + HPS
              if(q_Nzero) then
!               if (r .lt. rcutElecsq) then
                q = q_tab(atmType2, atmType1)
                r = sqrt(r)
                Ele = q * exp(-kapa * r) / r
                E_Ele = E_Ele + Ele
!               endif
              endif
            enddo
          enddo
        enddo
      enddo
!      write(*,*) E_Trial
      E_Trial = E_HPS + E_Ele
      
      end subroutine 
!======================================================================================================
      pure subroutine Rosen_Mol_Old_HPS_piecewise(mol_x, mol_y, mol_z, nType, included,  E_Trial)
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
      use Coords, only: MolArray
      use SimParameters, only: NPART, nMolTypes, kapa
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType
      real(dp), intent(in) :: mol_x(:), mol_y(:), mol_z(:)
      real(dp), intent(out) :: E_Trial
      
      integer :: iAtom, jType, jIndx, jMol, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: eps,sig_sq,q,lambda
      logical :: q_Nzero
      real(dp) :: HPS, LJ, Ele, x
      real(dp) :: E_Ele,E_HPS, rcutNBsq, rmax
      real(dp) :: rmin_ij

      E_Trial = 0E0
      E_HPS = 0E0
      E_Ele = 0E0      


      do iAtom = 1,nAtoms(nType)
        atmType1 = atomArray(nType, iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType, jAtom)
            eps = eps_tab(atmType2, atmType1)
            q_Nzero = q_Nonzero(atmType2, atmType1)
            sig_sq = sigsq_tab(atmType2, atmType1)  
            rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
            lambda = lambda_tab(atmType2, atmType1)  
            rmin_ij = r_min_tab(atmType2, atmType1)

            do jMol = 1,NPART(jType)
              jIndx = molArray(jType)%mol(jMol)%indx              
              if(included(jIndx) .eqv. .false.) then
                cycle
              endif          
              rx = abs(mol_x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom))
              ry = abs(mol_y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom))
              rz = abs(mol_z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom))
              rmax = max(rx, ry, rz) 
              if (rmax .gt. rcutElec) cycle
              if (.not. q_Nzero) then
                rmax = rmax * rmax
                if (rmax .gt. rcutNBsq) cycle
              endif
              r = rx*rx + ry*ry + rz*rz
              if (r .gt. rcutElecsq) cycle
              if(r .lt. rmin_ij) then
                E_Trial = huge(dp)
                return
              endif
              HPS = 0E0
              if (r .lt. rcutNBsq) then
                x = sig_sq/r
                x = x * x * x
                LJ = 4.0_dp * eps * x * (x - 1.0_dp)
                if (x .lt. 0.5_dp) then
                  HPS = lambda * LJ
                else
                  HPS = LJ + eps * (1.0_dp - lambda)
                endif
              endif
              E_HPS = E_HPS + HPS
              if(q_Nzero) then
!               if (r .lt. rcutElecsq) then
                q = q_tab(atmType2, atmType1)
                r = sqrt(r)
                Ele = q * exp(-kapa * r) / r
                E_Ele = E_Ele + Ele
!               endif
              endif
            enddo
          enddo
        enddo
      enddo
     
      E_Trial = E_HPS + E_Ele
      
      end subroutine 
!======================================================================================================
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which each atom is regrown sequentially for the given trial.
      pure subroutine Rosen_Atom_New_HPS_piecewise(nType, nAtom, trialPos, included,  E_Trial, overlap)
      use Coords, only: MolArray, SimpleAtomCoords
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
      use SimParameters, only: nMolTypes, NPART, kapa
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      
      logical, intent(inout) :: overlap
      real(dp), intent(out) :: E_Trial
      
      integer :: jType, jIndx, jMol, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: eps,sig_sq,q,lambda
      logical :: q_Nzero
      real(dp) :: HPS, LJ, Ele, x
      real(dp) :: E_Ele,E_HPS, rcutNBsq, rmax
      real(dp) :: rmin_ij

      
      E_HPS = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      overlap = .false.

      atmType1 = atomArray(nType, nAtom)
      do jType = 1, nMolTypes
        do jMol = 1,NPART(jType)
          jIndx = molArray(jType)%mol(jMol)%indx              
          if(included(jIndx) .eqv. .false.) then
            cycle
          endif      
          do jAtom = 1,nAtoms(jType)    
            rx = abs(trialPos%x - MolArray(jType)%mol(jMol)%x(jAtom))
            ry = abs(trialPos%y - MolArray(jType)%mol(jMol)%y(jAtom))
            rz = abs(trialPos%z - MolArray(jType)%mol(jMol)%z(jAtom)) 
            rmax = max(rx, ry, rz) 
            if (rmax .gt. rcutElec) cycle
            atmType2 = atomArray(jType, jAtom)
            rmin_ij = r_min_tab(atmType2, atmType1)  
            eps = eps_tab(atmType2, atmType1)
            q_Nzero = q_Nonzero(atmType2, atmType1)
            sig_sq = sigsq_tab(atmType2, atmType1)  
            rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
            lambda = lambda_tab(atmType2, atmType1)  
            if (.not. q_Nzero) then
              rmax = rmax * rmax
              if (rmax .gt. rcutNBsq) cycle
            endif

            r = rx*rx + ry*ry + rz*rz
            if (r .gt. rcutElecsq) cycle

            if(r .lt. rmin_ij) then
              E_Trial = huge(dp)
              overlap = .true.
              return
            endif
            HPS = 0E0
            if (r .lt. rcutNBsq) then
              x = sig_sq/r
              x = x * x * x
              LJ = 4.0_dp * eps * x * (x - 1.0_dp)
              if (x .lt. 0.5_dp) then
                HPS = lambda * LJ
              else
                HPS = LJ + eps * (1.0_dp - lambda)
              endif
            endif
            E_HPS = E_HPS + HPS
            if(q_Nzero) then
!             if (r .lt. rcutElecsq) then
              q = q_tab(atmType2, atmType1)
              r = sqrt(r)
              Ele = q * exp(-kapa * r) / r
              E_Ele = E_Ele + Ele
!             endif
            endif
          enddo
        enddo
      enddo
     
      E_Trial = E_HPS + E_Ele
      
      end subroutine 
!======================================================================================================
      pure subroutine Rosen_Atom_Old_HPS_piecewise(nType, nMol, nAtom, included, E_Trial)
      use Coords, only: MolArray
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
      use SimParameters, only: nMolTypes, NPART, kapa
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nAtom, nMol
      real(dp), intent(out) :: E_Trial
      
      integer :: jType, jIndx, jMol, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: eps,sig_sq,q,lambda
      logical :: q_Nzero
      real(dp) :: HPS, LJ, Ele, x
      real(dp) :: E_Ele,E_HPS, rcutNBsq, rmax
      real(dp) :: rmin_ij

      
      E_HPS = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0

      atmType1 = atomArray(nType, nAtom)
      do jType = 1, nMolTypes
        do jMol = 1,NPART(jType)
          jIndx = molArray(jType)%mol(jMol)%indx              
          if(included(jIndx) .eqv. .false.) then
            cycle
          endif              
          do jAtom = 1,nAtoms(jType)         
            atmType2 = atomArray(jType, jAtom)
            rmin_ij = r_min_tab(atmType2, atmType1)
            rx = abs(MolArray(nType)%mol(nMol)%x(nAtom) - MolArray(jType)%mol(jMol)%x(jAtom))
            ry = abs(MolArray(nType)%mol(nMol)%y(nAtom) - MolArray(jType)%mol(jMol)%y(jAtom))
            rz = abs(MolArray(nType)%mol(nMol)%z(nAtom) - MolArray(jType)%mol(jMol)%z(jAtom))
            rmax = max(rx, ry, rz) 
            if (rmax .gt. rcutElec) cycle
            eps = eps_tab(atmType2, atmType1)
            q_Nzero = q_Nonzero(atmType2, atmType1)
            sig_sq = sigsq_tab(atmType2, atmType1) 
            rcutNBsq = cutoffNBsq_tab(atmType2, atmType1) 
            lambda = lambda_tab(atmType2, atmType1) 
            if (.not. q_Nzero) then
              rmax = rmax * rmax
              if (rmax .gt. rcutNBsq) cycle
            endif
            r = rx*rx + ry*ry + rz*rz
            if (r .gt. rcutElecsq) cycle
            if(r .lt. rmin_ij) then
              E_Trial = huge(dp)
              return
            endif
 

            HPS = 0E0
            if (r .lt. rcutNBsq) then
              x = sig_sq/r
              x = x * x * x
              LJ = 4.0_dp * eps * x * (x - 1.0_dp)
              if (x .lt. 0.5_dp) then
                HPS = lambda * LJ
              else
                HPS = LJ + eps * (1.0_dp - lambda)
              endif
            endif
            E_HPS = E_HPS + HPS
            if(q_Nzero) then
!             if (r .lt. rcutElecsq) then
              q = q_tab(atmType2, atmType1)
              r = sqrt(r)
              Ele = q * exp(-kapa * r) / r
              E_Ele = E_Ele + Ele
!             endif
            endif
          enddo
        enddo
      enddo
     
      E_Trial = E_HPS + E_Ele
      
      end subroutine 
!======================================================================================================
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which each atom is regrown sequentially for the given trial.
      pure subroutine Rosen_Atom_Intra_New_HPS_piecewise(nType, nAtom, trialPos, regrown, E_Trial, overlap)
      use Coords, only: newMol, SimpleAtomCoords
      use ForceField, only: atomArray, nonBondArray, nIntraNonBond, r_min_tab
      use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
      use SimParameters, only: kapa
      implicit none

      logical, intent(in) :: regrown(:)
      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      logical, intent(inout) :: overlap
      real(dp), intent(out) :: E_Trial

      integer :: temp, iPair, iAtom, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: eps,sig_sq,q,lambda
      logical :: q_Nzero
      real(dp) :: HPS, LJ, Ele, x
      real(dp) :: E_Ele,E_HPS, rcutNBsq, rmax
      real(dp) :: rmin_ij

      E_HPS = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      overlap = .false.
      do iPair = 1,nIntraNonBond(nType)
        if(all(nonBondArray(nType,iPair)%nonMembr .ne. nAtom)) then
          cycle
        endif
        iAtom = nonBondArray(nType,iPair)%nonMembr(1)
        jAtom = nonBondArray(nType,iPair)%nonMembr(2)
        if(iAtom .ne. nAtom) then
          temp = iAtom
          iAtom = jAtom
          jAtom = temp
        endif
        if(regrown(jAtom) .eqv. .false.) then
          cycle
        endif
        rx = abs(trialPos%x - newMol%x(jAtom))
        ry = abs(trialPos%y - newMol%y(jAtom))
        rz = abs(trialPos%z - newMol%z(jAtom))
        rmax = max(rx, ry, rz) 
        if (rmax .gt. rcutElec) cycle
        atmType1 = atomArray(nType,iAtom)
        atmType2 = atomArray(nType,jAtom)
        eps = eps_tab(atmType1,atmType2)
        q_Nzero = q_Nonzero(atmType1,atmType2)
        sig_sq = sigsq_tab(atmType1,atmType2) 
        rcutNBsq = cutoffNBsq_tab(atmType1,atmType2)
        lambda = lambda_tab(atmType1,atmType2) 
        rmin_ij = r_min_tab(atmType1,atmType2)
        if (.not. q_Nzero) then
          rmax = rmax * rmax
          if (rmax .gt. rcutNBsq) cycle
        endif
 
        r = rx*rx + ry*ry + rz*rz
        if (r .gt. rcutElecsq) cycle
        if(r .lt. rmin_ij) then
          E_Trial = huge(dp)
          overlap = .true.
          return
        endif
        HPS = 0E0
        if (r .lt. rcutNBsq) then
          x = sig_sq/r
          x = x * x * x
          LJ = 4.0_dp * eps * x * (x - 1.0_dp)
          if (x .lt. 0.5_dp) then
            HPS = lambda * LJ
          else
            HPS = LJ + eps * (1.0_dp - lambda)
          endif
        endif
        E_HPS = E_HPS + HPS
        if(q_Nzero) then   
          q = q_tab(atmType1,atmType2)     
          r = sqrt(r)
          Ele = q * exp(-kapa * r)/r
          E_Ele = E_Ele + Ele       
        endif
      enddo
      E_Trial = E_HPS + E_Ele

      end subroutine 
!======================================================================================================
      pure subroutine Rosen_Atom_Intra_Old_HPS_piecewise(nType, nMol, nAtom, trialPos, regrown, E_Trial)
      use Coords, only: MolArray, SimpleAtomCoords
      use ForceField, only: atomArray, nonBondArray, nIntraNonBond, r_min_tab
      use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
      use SimParameters, only: kapa
      implicit none

      logical, intent(in) :: regrown(:)
      integer, intent(in) :: nType, nMol, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      real(dp), intent(out) :: E_Trial

      integer :: temp, iPair, iAtom, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: eps,sig_sq,q,lambda
      logical :: q_Nzero
      real(dp) :: HPS, LJ, Ele, x
      real(dp) :: E_Ele,E_HPS, rcutNBsq, rmax
      real(dp) :: rmin_ij

      E_HPS = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
!      overlap = .false.
      do iPair = 1,nIntraNonBond(nType)
        if(all(nonBondArray(nType,iPair)%nonMembr .ne. nAtom)) then
          cycle
        endif
        iAtom = nonBondArray(nType,iPair)%nonMembr(1)
        jAtom = nonBondArray(nType,iPair)%nonMembr(2)
        if(iAtom .ne. nAtom) then
          temp = iAtom
          iAtom = jAtom
          jAtom = temp
        endif
        if(regrown(jAtom) .eqv. .false.) then
          cycle
        endif
        rx = abs(trialPos%x - MolArray(nType)%mol(nMol)%x(jAtom))
        ry = abs(trialPos%y - MolArray(nType)%mol(nMol)%y(jAtom))
        rz = abs(trialPos%z - MolArray(nType)%mol(nMol)%z(jAtom))
        rmax = max(rx, ry, rz) 
        if (rmax .gt. rcutElec) cycle
        atmType1 = atomArray(nType,iAtom)
        atmType2 = atomArray(nType,jAtom)
        eps = eps_tab(atmType1,atmType2)
        q_Nzero = q_Nonzero(atmType1,atmType2)
        sig_sq = sigsq_tab(atmType1,atmType2) 
        rcutNBsq = cutoffNBsq_tab(atmType1,atmType2)
        lambda = lambda_tab(atmType1,atmType2) 
        rmin_ij = r_min_tab(atmType1,atmType2)
        if (.not. q_Nzero) then
          rmax = rmax * rmax
          if (rmax .gt. rcutNBsq) cycle
        endif
 
        r = rx*rx + ry*ry + rz*rz
        if (r .gt. rcutElecsq) cycle
        if(r .lt. rmin_ij) then
          E_Trial = huge(dp)
          return
        endif
        HPS = 0E0
        if (r .lt. rcutNBsq) then
          x = sig_sq/r
          x = x * x * x
          LJ = 4.0_dp * eps * x * (x - 1.0_dp)
          if (x .lt. 0.5_dp) then
            HPS = lambda * LJ
          else
            HPS = LJ + eps * (1.0_dp - lambda)
          endif
        endif
        E_HPS = E_HPS + HPS
        if(q_Nzero) then   
          q = q_tab(atmType1,atmType2)     
          r = sqrt(r)
          Ele = q * exp(-kapa * r)/r
          E_Ele = E_Ele + Ele       
        endif
      enddo
      E_Trial = E_HPS + E_Ele

      end subroutine 
!======================================================================================================
      pure subroutine Rosen_Atom_Intra_Gas_Old_HPS_piecewise(nType, nAtom, trialPos, regrown, E_Trial)
      use Coords, only: gasConfig, SimpleAtomCoords
      use ForceField, only: atomArray, nonBondArray, nIntraNonBond, r_min_tab
      use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
      use SimParameters, only: kapa
      implicit none

      logical, intent(in) :: regrown(:)
      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      real(dp), intent(out) :: E_Trial

      integer :: temp, iPair, iAtom, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: eps,sig_sq,q,lambda
      logical :: q_Nzero
      real(dp) :: HPS, LJ, Ele, x
      real(dp) :: E_Ele,E_HPS, rcutNBsq, rmax
      real(dp) :: rmin_ij

      E_HPS = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
!      overlap = .false.
      do iPair = 1,nIntraNonBond(nType)
        if(all(nonBondArray(nType,iPair)%nonMembr .ne. nAtom)) then
          cycle
        endif
        iAtom = nonBondArray(nType,iPair)%nonMembr(1)
        jAtom = nonBondArray(nType,iPair)%nonMembr(2)
        if(iAtom .ne. nAtom) then
          temp = iAtom
          iAtom = jAtom
          jAtom = temp
        endif
        if(regrown(jAtom) .eqv. .false.) then
          cycle
        endif
        rx = abs(trialPos%x - gasConfig(nType)%x(jAtom))
        ry = abs(trialPos%y - gasConfig(nType)%y(jAtom))
        rz = abs(trialPos%z - gasConfig(nType)%z(jAtom))
        rmax = max(rx, ry, rz) 
        if (rmax .gt. rcutElec) cycle
        atmType1 = atomArray(nType,iAtom)
        atmType2 = atomArray(nType,jAtom)
        eps = eps_tab(atmType1,atmType2)
        q_Nzero = q_Nonzero(atmType1,atmType2)
        sig_sq = sigsq_tab(atmType1,atmType2) 
        rcutNBsq = cutoffNBsq_tab(atmType1,atmType2)
        lambda = lambda_tab(atmType1,atmType2) 
        rmin_ij = r_min_tab(atmType1,atmType2)
        if (.not. q_Nzero) then
          rmax = rmax * rmax
          if (rmax .gt. rcutNBsq) cycle
        endif
 
        r = rx*rx + ry*ry + rz*rz
        if (r .gt. rcutElecsq) cycle
        if(r .lt. rmin_ij) then
          E_Trial = huge(dp)
          return
        endif
        HPS = 0E0
        if (r .lt. rcutNBsq) then
          x = sig_sq/r
          x = x * x * x
          LJ = 4.0_dp * eps * x * (x - 1.0_dp)
          if (x .lt. 0.5_dp) then
            HPS = lambda * LJ
          else
            HPS = LJ + eps * (1.0_dp - lambda)
          endif
        endif
        E_HPS = E_HPS + HPS
        if(q_Nzero) then   
          q = q_tab(atmType1,atmType2)     
          r = sqrt(r)
          Ele = q * exp(-kapa * r)/r
          E_Ele = E_Ele + Ele       
        endif
      enddo
      E_Trial = E_HPS + E_Ele

      end subroutine 
!======================================================================================================
! Computes intermolecular non-bonded energy for a trial atom position in a CBMC/FECBMC regrowth
! using the HPS_piecewise force field in a Grand Canonical Monte Carlo nucleation simulation.
! Evaluates piecewise Lennard-Jones (HPS) with switching at (sigma/r)^6 = 0.5 and Debye-Hückel
! electrostatic interactions with specified interacting atoms (InterAtoms). Checks for overlaps
! based on r_min_tab and returns total energy (E_Trial). Optimized for long linear biomolecules
! with sparse neighbor lists.
pure subroutine RosenLL_Atom_New_HPS_piecewise(nType, nAtom, trialPos, nInterAtoms, InterAtoms, E_Trial, overlap)
  use VarPrecision, only: dp, atomIntType
  use SimParameters, only: kapa
  use Coords, only: MolArray, SimpleAtomCoords
  use ForceField, only: r_min_tab, atomArray
  use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, &
                                        lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
  implicit none

  ! Input/Output variables
  integer, intent(in) :: nType, nAtom                ! Molecule type and atom index
  type(SimpleAtomCoords), intent(in) :: trialPos     ! Trial position coordinates
  integer, intent(in) :: nInterAtoms                 ! Number of interacting atoms
  integer, intent(in) :: InterAtoms(:, :)            ! Interacting atoms (type, mol, atom)
  logical, intent(inout) :: overlap                  ! Flag for overlap detection
  real(dp), intent(out) :: E_Trial                  ! Total intermolecular energy

  ! Local variables
  integer :: iInterAtoms, jType, jMol, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_sq, r, rmax, rmin_ij
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS_total, E_Ele_total
  logical :: q_Nzero

  ! Initialize energy terms and overlap flag
  E_HPS_total = 0E0_dp
  E_Ele_total = 0E0_dp
  E_Trial = 0E0_dp
  overlap = .false.

  ! Get atom type for trial position
  atmType1 = atomArray(nType, nAtom)

  ! Calculate intermolecular energies with interacting atoms
  do iInterAtoms = 1, nInterAtoms
    jType = InterAtoms(iInterAtoms, 1)
    jMol = InterAtoms(iInterAtoms, 2)
    jAtom = InterAtoms(iInterAtoms, 3)

    ! Calculate distance
    rx = abs(trialPos%x - MolArray(jType)%mol(jMol)%x(jAtom))
    ry = abs(trialPos%y - MolArray(jType)%mol(jMol)%y(jAtom))
    rz = abs(trialPos%z - MolArray(jType)%mol(jMol)%z(jAtom))
    rmax = max(rx, ry, rz)
    if (rmax > rcutElec) cycle

    ! Retrieve force field parameters
    atmType2 = atomArray(jType, jAtom)
    rmin_ij = r_min_tab(atmType2, atmType1)
    eps = eps_tab(atmType2, atmType1)
    q_Nzero = q_Nonzero(atmType2, atmType1)
    sig_sq = sigsq_tab(atmType2, atmType1)
    rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
    lambda = lambda_tab(atmType2, atmType1)
    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif

    r_sq = rx * rx + ry * ry + rz * rz
    if (r_sq > rcutElecsq) cycle

    ! Check for overlap
    if (r_sq < rmin_ij) then
      E_Trial = huge(dp)
      overlap = .true.
      return
    endif

    ! Calculate HPS energy
    HPS = 0E0_dp
    if (r_sq < rcutNBsq) then
      x = sig_sq / r_sq
      x = x * x * x
      LJ = 4.0_dp * eps * x * (x - 1.0_dp)
      if (x < 0.5_dp) then
        HPS = lambda * LJ
      else
        HPS = LJ + eps * (1.0_dp - lambda)
      endif
    endif
    E_HPS_total = E_HPS_total + HPS

    ! Calculate electrostatic energy
    Ele = 0E0_dp
    if (q_Nzero) then
      q = q_tab(atmType2, atmType1)
      r = sqrt(r_sq)
      if (kapa * r < 50.0_dp) then  ! Prevent underflow
        Ele = q * exp(-kapa * r) / r
      endif
      E_Ele_total = E_Ele_total + Ele
    endif
  enddo

  ! Compute total energy
  E_Trial = E_HPS_total + E_Ele_total

end subroutine RosenLL_Atom_New_HPS_piecewise
!======================================================================================================
! Computes intermolecular non-bonded energy for an existing atom in a CBMC/FECBMC regrowth
! using the HPS_piecewise force field in a Grand Canonical Monte Carlo nucleation simulation.
! Evaluates piecewise Lennard-Jones (HPS) with switching at (sigma/r)^6 = 0.5 and Debye-Hückel
! electrostatic interactions with specified interacting atoms (InterAtoms). Checks for overlaps
! based on r_min_tab and returns total energy (E_Trial). Optimized for long linear biomolecules
! with sparse neighbor lists.
pure subroutine RosenLL_Atom_Old_HPS_piecewise(nType, nMol, nAtom, nInterAtoms, InterAtoms, E_Trial)
  use VarPrecision, only: dp, atomIntType
  use SimParameters, only: kapa
  use Coords, only: MolArray
  use ForceField, only: r_min_tab, atomArray
  use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, &
                                        lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
  implicit none

  ! Input/Output variables
  integer, intent(in) :: nType                    ! Molecule type index
  integer, intent(in) :: nMol                     ! Molecule index within type
  integer, intent(in) :: nAtom                    ! Atom index within molecule
  integer, intent(in) :: nInterAtoms              ! Number of interacting atoms
  integer, intent(in) :: InterAtoms(:, :)         ! Interacting atoms (type, mol, atom)
  real(dp), intent(out) :: E_Trial                ! Total non-bonded energy (HPS + Ele)

  ! Local variables
  integer :: iInterAtoms, jType, jMol, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_sq, r, rmax, rmin_ij
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS_total, E_Ele_total
  logical :: q_Nzero

  ! Initialize energy terms
  E_HPS_total = 0E0_dp
  E_Ele_total = 0E0_dp
  E_Trial = 0E0_dp

  ! Get atom type for the existing atom
  atmType1 = atomArray(nType, nAtom)

  ! Calculate intermolecular energies with interacting atoms
  do iInterAtoms = 1, nInterAtoms
    jType = InterAtoms(iInterAtoms, 1)
    jMol = InterAtoms(iInterAtoms, 2)
    jAtom = InterAtoms(iInterAtoms, 3)

    ! Calculate distance
    rx = abs(MolArray(nType)%mol(nMol)%x(nAtom) - MolArray(jType)%mol(jMol)%x(jAtom))
    ry = abs(MolArray(nType)%mol(nMol)%y(nAtom) - MolArray(jType)%mol(jMol)%y(jAtom))
    rz = abs(MolArray(nType)%mol(nMol)%z(nAtom) - MolArray(jType)%mol(jMol)%z(jAtom))
    rmax = max(rx, ry, rz)
    if (rmax > rcutElec) cycle

    ! Retrieve force field parameters
    atmType2 = atomArray(jType, jAtom)
    rmin_ij = r_min_tab(atmType2, atmType1)
    eps = eps_tab(atmType2, atmType1)
    q_Nzero = q_Nonzero(atmType2, atmType1)
    sig_sq = sigsq_tab(atmType2, atmType1)
    rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
    lambda = lambda_tab(atmType2, atmType1)
    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif

    r_sq = rx * rx + ry * ry + rz * rz
    if (r_sq > rcutElecsq) cycle

    ! Check for overlap
    if (r_sq < rmin_ij) then
      E_Trial = huge(dp)
      return
    endif

    ! Calculate HPS energy
    HPS = 0E0_dp
    if (r_sq < rcutNBsq) then
      x = sig_sq / r_sq
      x = x * x * x
      LJ = 4.0_dp * eps * x * (x - 1.0_dp)
      if (x < 0.5_dp) then
        HPS = lambda * LJ
      else
        HPS = LJ + eps * (1.0_dp - lambda)
      endif
    endif
    E_HPS_total = E_HPS_total + HPS

    ! Calculate electrostatic energy
    Ele = 0E0_dp
    if (q_Nzero) then
      q = q_tab(atmType2, atmType1)
      r = sqrt(r_sq)
      if (kapa * r < 50.0_dp) then  ! Prevent underflow
        Ele = q * exp(-kapa * r) / r
      endif
      E_Ele_total = E_Ele_total + Ele
    endif
  enddo

  ! Compute total energy
  E_Trial = E_HPS_total + E_Ele_total

end subroutine RosenLL_Atom_Old_HPS_piecewise
!======================================================================================================
! Computes intramolecular non-bonded energy for a trial atom position in a CBMC/FECBMC regrowth
! using the HPS_piecewise force field in a Grand Canonical Monte Carlo nucleation simulation.
! Evaluates piecewise Lennard-Jones (HPS) with switching at (sigma/r)^6 = 0.5 and Debye-Hückel
! electrostatic interactions with specified intramolecular atoms (IntraUnits). Checks for overlaps
! based on r_min_tab and returns total energy (E_Trial). Optimized for long linear biomolecules
! with sparse non-bonded pairs.
pure subroutine RosenLL_Atom_Intra_New_HPS_piecewise(nType, nAtom, trialPos, nIntraUnits, IntraUnits, E_Trial, overlap)
  use VarPrecision, only: dp, atomIntType
  use SimParameters, only: kapa
  use Coords, only: newMol, SimpleAtomCoords
  use ForceField, only: atomArray, r_min_tab
  use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, &
                                        lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
  implicit none

  ! Input/Output variables
  integer, intent(in) :: nType                     ! Molecule type index
  integer, intent(in) :: nAtom                     ! Atom index within molecule
  type(SimpleAtomCoords), intent(in) :: trialPos   ! Trial coordinates (x, y, z)
  integer, intent(in) :: nIntraUnits               ! Number of intramolecular interacting atoms
  integer, intent(in) :: IntraUnits(:)             ! Indices of interacting atoms
  real(dp), intent(out) :: E_Trial                 ! Total non-bonded energy (HPS + Ele)
  logical, intent(inout) :: overlap                  ! True if any pair distance < r_min_tab

  ! Local variables
  integer :: iPair, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_sq, r, rmax, rmin_ij
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS_total, E_Ele_total
  logical :: q_Nzero

  ! Initialize energy terms and overlap flag
  E_HPS_total = 0E0_dp
  E_Ele_total = 0E0_dp
  E_Trial = 0E0_dp
  overlap = .false.

  ! Get atom type for trial position
  atmType1 = atomArray(nType, nAtom)

  ! Calculate intramolecular energies with interacting atoms
  do iPair = 1, nIntraUnits
    jAtom = IntraUnits(iPair)

    ! Calculate distance
    rx = abs(trialPos%x - newMol%x(jAtom))
    ry = abs(trialPos%y - newMol%y(jAtom))
    rz = abs(trialPos%z - newMol%z(jAtom))
    rmax = max(rx, ry, rz)
    if (rmax > rcutElec) cycle

    ! Retrieve force field parameters
    atmType2 = atomArray(nType, jAtom)
    eps = eps_tab(atmType1, atmType2)
    q_Nzero = q_Nonzero(atmType1, atmType2)
    sig_sq = sigsq_tab(atmType1, atmType2)
    rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
    lambda = lambda_tab(atmType1, atmType2)
    rmin_ij = r_min_tab(atmType1, atmType2)
    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif

    r_sq = rx * rx + ry * ry + rz * rz
    if (r_sq > rcutElecsq) cycle

    ! Check for overlap
    if (r_sq < rmin_ij) then
      E_Trial = huge(dp)
      overlap = .true.
      return
    endif

    ! Calculate HPS energy
    HPS = 0E0_dp
    if (r_sq < rcutNBsq) then
      x = sig_sq / r_sq
      x = x * x * x
      LJ = 4.0_dp * eps * x * (x - 1.0_dp)
      if (x < 0.5_dp) then
        HPS = lambda * LJ
      else
        HPS = LJ + eps * (1.0_dp - lambda)
      endif
    endif
    E_HPS_total = E_HPS_total + HPS

    ! Calculate electrostatic energy
    Ele = 0E0_dp
    if (q_Nzero) then
      q = q_tab(atmType1, atmType2)
      r = sqrt(r_sq)
      if (kapa * r < 50.0_dp) then  ! Prevent underflow
        Ele = q * exp(-kapa * r) / r
      endif
      E_Ele_total = E_Ele_total + Ele
    endif
  enddo

  ! Compute total energy
  E_Trial = E_HPS_total + E_Ele_total

end subroutine RosenLL_Atom_Intra_New_HPS_piecewise
!======================================================================================================
! Computes intramolecular non-bonded energy for a trial atom position in a CBMC/FECBMC regrowth
! using the HPS_piecewise force field in a Grand Canonical Monte Carlo nucleation simulation.
! Evaluates piecewise Lennard-Jones (HPS) with switching at (sigma/r)^6 = 0.5 and Debye-Hückel
! electrostatic interactions with specified intramolecular atoms (IntraUnits). Checks for overlaps
! based on r_min_tab and returns total energy (E_Trial). Optimized for long linear biomolecules
! with sparse non-bonded pairs.
pure subroutine RosenLL_Atom_Intra_Old_HPS_piecewise(nType, nMol, nAtom, trialPos, nIntraUnits, IntraUnits, E_Trial)
  use VarPrecision, only: dp, atomIntType
  use SimParameters, only: kapa
  use Coords, only: MolArray, SimpleAtomCoords
  use ForceField, only: atomArray, r_min_tab
  use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, &
                                        lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
  implicit none

  ! Input/Output variables
  integer, intent(in) :: nType                    ! Molecule type index
  integer, intent(in) :: nMol                     ! Molecule index within type
  integer, intent(in) :: nAtom                    ! Atom index within molecule
  type(SimpleAtomCoords), intent(in) :: trialPos  ! Trial coordinates (x, y, z)
  integer, intent(in) :: nIntraUnits              ! Number of intramolecular interacting atoms
  integer, intent(in) :: IntraUnits(:)            ! Indices of interacting atoms
  real(dp), intent(out) :: E_Trial                ! Total non-bonded energy (HPS + Ele)

  ! Local variables
  integer :: iPair, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_sq, r, rmax, rmin_ij
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS_total, E_Ele_total
  logical :: q_Nzero

  ! Initialize energy terms
  E_HPS_total = 0E0_dp
  E_Ele_total = 0E0_dp
  E_Trial = 0E0_dp

  ! Get atom type for trial position
  atmType1 = atomArray(nType, nAtom)

  ! Calculate intramolecular energies with interacting atoms
  do iPair = 1, nIntraUnits
    jAtom = IntraUnits(iPair)

    ! Calculate distance
    rx = abs(trialPos%x - MolArray(nType)%mol(nMol)%x(jAtom))
    ry = abs(trialPos%y - MolArray(nType)%mol(nMol)%y(jAtom))
    rz = abs(trialPos%z - MolArray(nType)%mol(nMol)%z(jAtom))
    rmax = max(rx, ry, rz)
    if (rmax > rcutElec) cycle

    ! Retrieve force field parameters
    atmType2 = atomArray(nType, jAtom)
    eps = eps_tab(atmType1, atmType2)
    q_Nzero = q_Nonzero(atmType1, atmType2)
    sig_sq = sigsq_tab(atmType1, atmType2)
    rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
    lambda = lambda_tab(atmType1, atmType2)
    rmin_ij = r_min_tab(atmType1, atmType2)
    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif

    r_sq = rx * rx + ry * ry + rz * rz
    if (r_sq > rcutElecsq) cycle

    ! Check for overlap
    if (r_sq < rmin_ij) then
      E_Trial = huge(dp)
      return
    endif

    ! Calculate HPS energy
    HPS = 0E0_dp
    if (r_sq < rcutNBsq) then
      x = sig_sq / r_sq
      x = x * x * x
      LJ = 4.0_dp * eps * x * (x - 1.0_dp)
      if (x < 0.5_dp) then
        HPS = lambda * LJ
      else
        HPS = LJ + eps * (1.0_dp - lambda)
      endif
    endif
    E_HPS_total = E_HPS_total + HPS

    ! Calculate electrostatic energy
    Ele = 0E0_dp
    if (q_Nzero) then
      q = q_tab(atmType1, atmType2)
      r = sqrt(r_sq)
      if (kapa * r < 50.0_dp) then  ! Prevent underflow
        Ele = q * exp(-kapa * r) / r
      endif
      E_Ele_total = E_Ele_total + Ele
    endif
  enddo

  ! Compute total energy
  E_Trial = E_HPS_total + E_Ele_total

end subroutine RosenLL_Atom_Intra_Old_HPS_piecewise
!======================================================================================================
! Computes intramolecular non-bonded energy for a trial segment position in the gas phase using
! the HPS_piecewise force field in a Grand Canonical Monte Carlo nucleation simulation.
! Evaluates piecewise Lennard-Jones (HPS) with switching at (sigma/r)^6 = 0.5 and Debye-Hückel
! electrostatic interactions with specified intramolecular segments (IntraUnits) in gasConfig.
! Checks for overlaps based on r_min_tab and returns total energy (E_Trial) for the Rosenbluth
! weight in the reverse move (new to old) of LongChain_RosenConfigGen_GasPhase_Reverse.
pure subroutine RosenLL_Atom_Intra_Gas_Old_HPS_piecewise(nType, nAtom, trialPos, nIntraUnits, IntraUnits, E_Trial)
  use VarPrecision, only: dp, atomIntType
  use SimParameters, only: kapa
  use Coords, only: gasConfig, SimpleAtomCoords
  use ForceField, only: atomArray, r_min_tab
  use ForceFieldPara_HPS_piecewise, only: eps_tab, q_Nonzero, q_tab, sigsq_tab, &
                                        lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
  implicit none

  ! Input/Output variables
  integer, intent(in) :: nType                    ! Molecule type index
  integer, intent(in) :: nAtom                    ! Segment index for trial position
  type(SimpleAtomCoords), intent(in) :: trialPos  ! Trial coordinates (x, y, z)
  integer, intent(in) :: nIntraUnits              ! Number of intramolecular segments
  integer, intent(in) :: IntraUnits(:)            ! Indices of interacting segments
  real(dp), intent(out) :: E_Trial                ! Total non-bonded energy (HPS + Ele)

  ! Local variables
  integer :: iPair, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r_sq, r, rmax, rmin_ij
  real(dp) :: eps, sig_sq, q, lambda, rcutNBsq
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_HPS_total, E_Ele_total
  logical :: q_Nzero

  ! Initialize energy terms
  E_HPS_total = 0E0_dp
  E_Ele_total = 0E0_dp
  E_Trial = 0E0_dp

  ! Get atom type for trial position
  atmType1 = atomArray(nType, nAtom)

  ! Calculate intramolecular energies with interacting segments
  do iPair = 1, nIntraUnits
    jAtom = IntraUnits(iPair)

    ! Calculate distance
    rx = abs(trialPos%x - gasConfig(nType)%x(jAtom))
    ry = abs(trialPos%y - gasConfig(nType)%y(jAtom))
    rz = abs(trialPos%z - gasConfig(nType)%z(jAtom))
    rmax = max(rx, ry, rz)
    if (rmax > rcutElec) cycle

    ! Retrieve force field parameters
    atmType2 = atomArray(nType, jAtom)
    eps = eps_tab(atmType1, atmType2)
    q_Nzero = q_Nonzero(atmType1, atmType2)
    sig_sq = sigsq_tab(atmType1, atmType2)
    rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
    lambda = lambda_tab(atmType1, atmType2)
    rmin_ij = r_min_tab(atmType1, atmType2)
    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif

    r_sq = rx * rx + ry * ry + rz * rz
    if (r_sq > rcutElecsq) cycle

    ! Check for overlap
    if (r_sq < rmin_ij) then
      E_Trial = huge(dp)
      return
    endif

    ! Calculate HPS energy
    HPS = 0E0_dp
    if (r_sq < rcutNBsq) then
      x = sig_sq / r_sq
      x = x * x * x
      LJ = 4.0_dp * eps * x * (x - 1.0_dp)
      if (x < 0.5_dp) then
        HPS = lambda * LJ
      else
        HPS = LJ + eps * (1.0_dp - lambda)
      endif
    endif
    E_HPS_total = E_HPS_total + HPS

    ! Calculate electrostatic energy
    Ele = 0E0_dp
    if (q_Nzero) then
      q = q_tab(atmType1, atmType2)
      r = sqrt(r_sq)
      if (kapa * r < 50.0_dp) then  ! Prevent underflow
        Ele = q * exp(-kapa * r) / r
      endif
      E_Ele_total = E_Ele_total + Ele
    endif
  enddo

  ! Compute total energy
  E_Trial = E_HPS_total + E_Ele_total

end subroutine RosenLL_Atom_Intra_Gas_Old_HPS_piecewise
!======================================================================================================
! Identifies interacting atoms within a cutoff distance from a reference atom (RefAtm) for the
! HPS_piecewise force field in a CBMC/FECBMC move for long linear biomolecules. Computes distances
! to all atoms in other molecules, including a buffer distance (r0), and accounts for non-bonded
! (piecewise HPS) and electrostatic interactions based on atom types. Stores interacting atom indices
! (molecule type, molecule index, atom index) in InterAtoms and counts them in nInterAtoms. Excludes
! the molecule with index nIndx. Called during regrowth for accurate energy calculations.
pure subroutine Find_InterAtms_HPS_piecewise(RefAtm, r0, atmType1, nIndx, nInterAtoms, InterAtoms)
  use VarPrecision, only: dp, atomIntType
  use SimParameters, only: nMolTypes, NPART
  use Coords, only: MolArray, SimpleAtomCoords
  use ForceField, only: nAtoms, atomArray
  use ForceFieldPara_HPS_piecewise, only: cutoffNB_tab, q_Nonzero, rcutElec
  implicit none

  ! Input/Output variables
  type(SimpleAtomCoords), intent(in) :: RefAtm       ! Position of reference atom
  real(dp), intent(in) :: r0                         ! Buffer distance for cutoff
  integer(kind=atomIntType), intent(in) :: atmType1  ! Atom type of reference atom
  integer, intent(in) :: nIndx                      ! Global index of molecule to exclude
  integer, intent(inout) :: nInterAtoms              ! Number of interacting atoms
  integer, intent(out) :: InterAtoms(:, :)           ! Indices (type, mol, atom) of interacting atoms

  ! Local variables
  integer :: jType, jMol, jAtom, jIndx
  integer(kind=atomIntType) :: atmType2
  real(dp) :: rx, ry, rz, r_sq, rmax
  real(dp) :: rcutNB, rcutshort, rcutlong, rcutshortsq, rcutlongsq
  logical :: q_Nzero

  ! Initialize outputs
  nInterAtoms = 0
  InterAtoms = 0

  ! Set cutoff distances
  rcutlong = r0 + rcutElec
  rcutlongsq = rcutlong * rcutlong

  ! Loop over all molecules and atoms, excluding nIndx
  do jType = 1, nMolTypes
    do jMol = 1, NPART(jType)
      jIndx = MolArray(jType)%mol(jMol)%indx
      if (jIndx == nIndx) cycle
      do jAtom = 1, nAtoms(jType)
        ! Calculate distance
        rx = abs(RefAtm%x - MolArray(jType)%mol(jMol)%x(jAtom))
        ry = abs(RefAtm%y - MolArray(jType)%mol(jMol)%y(jAtom))
        rz = abs(RefAtm%z - MolArray(jType)%mol(jMol)%z(jAtom))
        rmax = max(rx, ry, rz)
        if (rmax > rcutlong) cycle

        r_sq = rx * rx + ry * ry + rz * rz
        if (r_sq > rcutlongsq) cycle

        ! Check interaction based on atom types
        atmType2 = atomArray(jType, jAtom)
        q_Nzero = q_Nonzero(atmType2, atmType1)
        rcutNB = cutoffNB_tab(atmType2, atmType1)
        rcutshort = r0 + rcutNB
        rcutshortsq = rcutshort * rcutshort
        if (.not. q_Nzero) then
          if (r_sq > rcutshortsq) cycle
        endif

        ! Store interacting atom
        nInterAtoms = nInterAtoms + 1
        InterAtoms(nInterAtoms, 1) = jType
        InterAtoms(nInterAtoms, 2) = jMol
        InterAtoms(nInterAtoms, 3) = jAtom
      enddo
    enddo
  enddo

end subroutine Find_InterAtms_HPS_piecewise
!======================================================================================================
  
      end module 


