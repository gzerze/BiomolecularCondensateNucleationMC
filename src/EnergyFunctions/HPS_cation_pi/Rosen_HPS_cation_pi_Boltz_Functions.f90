!   HPS_piecewise (HPS-KR, KH, FB, TSCL-M2, Urry: Piecewise HPS LJ + Debye-Hückel
!      - u_LJ(r)           = 4 * ε * [ (σ / r)**12 − (σ / r)**6 ]
!      - U_HPS(r) =
!            u_LJ(r) + ε * (1 − λ)       , if r <= 2**(1/6) * σ
!            λ * u_LJ(r)                 , if r >  2**(1/6) * σ
!      - U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!
!   HPS_cation_pi (HPS-cation-pi-i and -ii): HPS_piecewise + LJ for specific π interactions
!      - U_total(r)        = U_HPS_piecewise(r) + U_extra_LJ(r) for cation-π pairs
!      - U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!*********************************************************************************************************************
      module Rosenbluth_Functions_HPS_cation_pi
      use VarPrecision
      contains
!====================================================================================== 
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which regrows an entire molecule for the given trial.
      pure subroutine Rosen_Mol_New_HPS_cation_pi(nRosen, nType, included,  E_Trial, overlap)
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
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
      real(dp) :: eps,sig_sq,q,lambda,epsLJ
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
            epsLJ = epsLJ_tab(atmType2, atmType1)
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
                if (epsLJ .gt. 1E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
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
      pure subroutine Rosen_Mol_Old_HPS_cation_pi(mol_x, mol_y, mol_z, nType, included,  E_Trial)
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
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
      real(dp) :: eps,sig_sq,q,lambda,epsLJ
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
            epsLJ = epsLJ_tab(atmType2, atmType1)
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
                if (epsLJ .gt. 1E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
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
      pure subroutine Rosen_Atom_New_HPS_cation_pi(nType, nAtom, trialPos, included,  E_Trial, overlap)
      use Coords, only: MolArray, SimpleAtomCoords
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
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
      real(dp) :: eps,sig_sq,q,lambda,epsLJ
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
            epsLJ = epsLJ_tab(atmType2, atmType1)
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
              if (epsLJ .gt. 1E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
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
      pure subroutine Rosen_Atom_Old_HPS_cation_pi(nType, nMol, nAtom, included, E_Trial)
      use Coords, only: MolArray
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
      use SimParameters, only: nMolTypes, NPART, kapa
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nAtom, nMol
      real(dp), intent(out) :: E_Trial
      
      integer :: jType, jIndx, jMol, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: eps,sig_sq,q,lambda,epsLJ
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
            epsLJ = epsLJ_tab(atmType2, atmType1)
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
              if (epsLJ .gt. 1E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
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
      pure subroutine Rosen_Atom_Intra_New_HPS_cation_pi(nType, nAtom, trialPos, regrown, E_Trial, overlap)
      use Coords, only: newMol, SimpleAtomCoords
      use ForceField, only: atomArray, nonBondArray, nIntraNonBond, r_min_tab
      use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
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
      real(dp) :: eps,sig_sq,q,lambda,epsLJ
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
        epsLJ = epsLJ_tab(atmType1,atmType2)
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
          if (epsLJ .gt. 1E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
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
      pure subroutine Rosen_Atom_Intra_Old_HPS_cation_pi(nType, nMol, nAtom, trialPos, regrown, E_Trial)
      use Coords, only: MolArray, SimpleAtomCoords
      use ForceField, only: atomArray, nonBondArray, nIntraNonBond, r_min_tab
      use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
      use SimParameters, only: kapa
      implicit none

      logical, intent(in) :: regrown(:)
      integer, intent(in) :: nType, nMol, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      real(dp), intent(out) :: E_Trial

      integer :: temp, iPair, iAtom, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: eps,sig_sq,q,lambda,epsLJ
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
        epsLJ = epsLJ_tab(atmType1,atmType2)
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
          if (epsLJ .gt. 1E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
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
      pure subroutine Rosen_Atom_Intra_Gas_Old_HPS_cation_pi(nType, nAtom, trialPos, regrown, E_Trial)
      use Coords, only: gasConfig, SimpleAtomCoords
      use ForceField, only: atomArray, nonBondArray, nIntraNonBond, r_min_tab
      use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
      use SimParameters, only: kapa
      implicit none

      logical, intent(in) :: regrown(:)
      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      real(dp), intent(out) :: E_Trial

      integer :: temp, iPair, iAtom, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: eps,sig_sq,q,lambda,epsLJ
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
        epsLJ = epsLJ_tab(atmType1,atmType2)
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
          if (epsLJ .gt. 1E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
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
pure subroutine RosenLL_Atom_New_HPS_cation_pi(nType, nAtom, trialPos, nInterAtoms, InterAtoms, E_Trial, overlap)
  ! Computes intermolecular non-bonded energy for a trial atom position in a CBMC/FECBMC regrowth
  ! using the HPS_cation_pi force field in a grand canonical ensemble nucleation simulation.
  ! Evaluates piecewise Hydropathy Scale Lennard-Jones, cation-π interactions, and Debye-Hückel
  ! electrostatic interactions with specified interacting atoms. Checks for overlaps based on
  ! minimum distance criterion and returns total energy. Optimized for long linear biomolecules
  ! with sparse neighbor lists.
  ! HPS_cation_pi potential:
  !   u_LJ(r) = 4 * eps * [ (sigma/r)^12 - (sigma/r)^6 ]
  !   U_HPS(r) = u_LJ(r) + eps * (1 - lambda), if r <= 2^(1/6) * sigma
  !              lambda * u_LJ(r), if r > 2^(1/6) * sigma
  !   U_cation_pi(r) = 4 * epsLJ * [ (sigma/r)^12 - (sigma/r)^6 ], for specific pairs
  !   U_Electrostatic(r) = (q1 * q2 / r) * exp(-kappa * r)

  use SimParameters, only: dp, kapa               ! Precision and Debye-Hückel screening parameter
  use Coords, only: MolArray, SimpleAtomCoords    ! Molecule coordinates and trial position type
  use ForceField, only: r_min_tab, atomArray      ! Minimum distance table and atom types
  use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, &
                                          lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
  implicit none

  ! Input variables
  integer, intent(in) :: nType                    ! Molecule type index
  integer, intent(in) :: nAtom                    ! Atom index within molecule
  type(SimpleAtomCoords), intent(in) :: trialPos  ! Trial coordinates (x, y, z)
  integer, intent(in) :: nInterAtoms              ! Number of interacting atoms
  integer, intent(in) :: InterAtoms(:, :)         ! Interacting atoms (type, molecule, atom indices)

  ! Output variables
  real(dp), intent(out) :: E_Trial                ! Total non-bonded energy (HPS + cation-π + electrostatic)
  logical, intent(inout) :: overlap                 ! True if any atom pair distance < r_min_tab

  ! Local variables
  integer :: iInterAtoms, jType, jMol, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: eps, sig_sq, q, lambda, epsLJ
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS, rcutNBsq, rmax
  real(dp) :: rmin_ij

  ! Initialize energy and overlap variables
  E_HPS = 0.0_dp
  E_Ele = 0.0_dp      
  E_Trial = 0.0_dp
  overlap = .false.

  ! Get atom type for trial position
  atmType1 = atomArray(nType, nAtom)
  ! Loop over interacting atoms
  do iInterAtoms = 1, nInterAtoms
    jType = InterAtoms(iInterAtoms, 1)
    jMol = InterAtoms(iInterAtoms, 2) 
    jAtom = InterAtoms(iInterAtoms, 3)
    ! Calculate distance components
    rx = abs(trialPos%x - MolArray(jType)%mol(jMol)%x(jAtom))
    ry = abs(trialPos%y - MolArray(jType)%mol(jMol)%y(jAtom))
    rz = abs(trialPos%z - MolArray(jType)%mol(jMol)%z(jAtom)) 
    rmax = max(rx, ry, rz) 
    ! Check electrostatic cutoff
    if (rmax > rcutElec) cycle
    ! Retrieve force field parameters
    atmType2 = atomArray(jType, jAtom)
    rmin_ij = r_min_tab(atmType2, atmType1)  
    eps = eps_tab(atmType2, atmType1)
    epsLJ = epsLJ_tab(atmType2, atmType1)
    q_Nzero = q_Nonzero(atmType2, atmType1)
    sig_sq = sigsq_tab(atmType2, atmType1)  
    rcutNBsq = cutoffNBsq_tab(atmType2, atmType1) 
    lambda = lambda_tab(atmType2, atmType1)  
    ! Check non-bonded cutoff for non-charged pairs
    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif
    ! Compute squared distance
    r = rx*rx + ry*ry + rz*rz
    if (r > rcutElecsq) cycle
    ! Check for atomic overlaps
    if (r < rmin_ij) then
      E_Trial = huge(dp)
      overlap = .true.
      return
    endif
    ! Calculate HPS and LJ energies
    HPS = 0.0_dp
    if (r < rcutNBsq) then
      x = sig_sq/r
      x = x * x * x
      LJ = 4.0_dp * eps * x * (x - 1.0_dp)
      if (x < 0.5_dp) then
        HPS = lambda * LJ
      else
        HPS = LJ + eps * (1.0_dp - lambda)
      endif
      if (epsLJ > 1.0E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
    endif
    E_HPS = E_HPS + HPS
    ! Calculate electrostatic energy
    if (q_Nzero) then
      q = q_tab(atmType2, atmType1)
      r = sqrt(r)
      Ele = q * exp(-kapa * r) / r
      E_Ele = E_Ele + Ele
    endif
  enddo

  ! Compute total non-bonded energy
  E_Trial = E_HPS + E_Ele

end subroutine RosenLL_Atom_New_HPS_cation_pi
!======================================================================================================
pure subroutine RosenLL_Atom_Old_HPS_cation_pi(nType, nMol, nAtom, nInterAtoms, InterAtoms, E_Trial)
  ! Computes intermolecular non-bonded energy for an existing atom in a CBMC/FECBMC regrowth
  ! using the HPS_cation_pi force field in a grand canonical ensemble nucleation simulation.
  ! Evaluates piecewise Hydropathy Scale Lennard-Jones, cation-π interactions, and Debye-Hückel
  ! electrostatic interactions with specified interacting atoms. Checks for overlaps based on
  ! minimum distance criterion and returns total energy. Optimized for long linear biomolecules
  ! with sparse neighbor lists.
  ! HPS_cation_pi potential:
  !   u_LJ(r) = 4 * eps * [ (sigma/r)^12 - (sigma/r)^6 ]
  !   U_HPS(r) = u_LJ(r) + eps * (1 - lambda), if r <= 2^(1/6) * sigma
  !              lambda * u_LJ(r), if r > 2^(1/6) * sigma
  !   U_cation_pi(r) = 4 * epsLJ * [ (sigma/r)^12 - (sigma/r)^6 ], for specific pairs
  !   U_Electrostatic(r) = (q1 * q2 / r) * exp(-kappa * r)

  use SimParameters, only: dp, kapa               ! Precision and Debye-Hückel screening parameter
  use Coords, only: MolArray                     ! Molecule coordinates
  use ForceField, only: r_min_tab, atomArray     ! Minimum distance table and atom types
  use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, &
                                          lambda_tab, cutoffNBsq_tab, rcutElecsq, rcutElec
  implicit none

  ! Input variables
  integer, intent(in) :: nType                    ! Molecule type index
  integer, intent(in) :: nMol                    ! Molecule index within type
  integer, intent(in) :: nAtom                   ! Atom index within molecule
  integer, intent(in) :: nInterAtoms             ! Number of interacting atoms
  integer, intent(in) :: InterAtoms(:, :)        ! Interacting atoms (type, molecule, atom indices)

  ! Output variables
  real(dp), intent(out) :: E_Trial               ! Total non-bonded energy (HPS + cation-π + electrostatic)

  ! Local variables
  integer :: iInterAtoms, jType, jMol, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: eps, sig_sq, q, lambda, epsLJ
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS, rcutNBsq, rmax
  real(dp) :: rmin_ij

  ! Initialize energy variables
  E_HPS = 0.0_dp
  E_Ele = 0.0_dp      
  E_Trial = 0.0_dp

  ! Get atom type for the existing atom
  atmType1 = atomArray(nType, nAtom)
  ! Loop over interacting atoms
  do iInterAtoms = 1, nInterAtoms
    jType = InterAtoms(iInterAtoms, 1)
    jMol = InterAtoms(iInterAtoms, 2) 
    jAtom = InterAtoms(iInterAtoms, 3)
    ! Calculate distance components
    rx = abs(MolArray(nType)%mol(nMol)%x(nAtom) - MolArray(jType)%mol(jMol)%x(jAtom))
    ry = abs(MolArray(nType)%mol(nMol)%y(nAtom) - MolArray(jType)%mol(jMol)%y(jAtom))
    rz = abs(MolArray(nType)%mol(nMol)%z(nAtom) - MolArray(jType)%mol(jMol)%z(jAtom)) 
    rmax = max(rx, ry, rz) 
    ! Check electrostatic cutoff
    if (rmax > rcutElec) cycle
    ! Retrieve force field parameters
    atmType2 = atomArray(jType, jAtom)
    rmin_ij = r_min_tab(atmType2, atmType1)  
    eps = eps_tab(atmType2, atmType1)
    epsLJ = epsLJ_tab(atmType2, atmType1)
    q_Nzero = q_Nonzero(atmType2, atmType1)
    sig_sq = sigsq_tab(atmType2, atmType1)  
    rcutNBsq = cutoffNBsq_tab(atmType2, atmType1)
    lambda = lambda_tab(atmType2, atmType1)  
    ! Check non-bonded cutoff for non-charged pairs
    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif
    ! Compute squared distance
    r = rx*rx + ry*ry + rz*rz
    if (r > rcutElecsq) cycle
    ! Check for atomic overlaps
    if (r < rmin_ij) then
      E_Trial = huge(dp)
      return
    endif
    ! Calculate HPS and LJ energies
    HPS = 0.0_dp
    if (r < rcutNBsq) then
      x = sig_sq/r
      x = x * x * x
      LJ = 4.0_dp * eps * x * (x - 1.0_dp)
      if (x < 0.5_dp) then
        HPS = lambda * LJ
      else
        HPS = LJ + eps * (1.0_dp - lambda)
      endif
      if (epsLJ > 1.0E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
    endif
    E_HPS = E_HPS + HPS
    ! Calculate electrostatic energy
    if (q_Nzero) then
      q = q_tab(atmType2, atmType1)
      r = sqrt(r)
      Ele = q * exp(-kapa * r) / r
      E_Ele = E_Ele + Ele
    endif
  enddo

  ! Compute total non-bonded energy
  E_Trial = E_HPS + E_Ele

end subroutine RosenLL_Atom_Old_HPS_cation_pi
!======================================================================================================
pure subroutine RosenLL_Atom_Intra_New_HPS_cation_pi(nType, nAtom, trialPos, nIntraUnits, IntraUnits, E_Trial, overlap)
  ! Computes intramolecular non-bonded energy for a trial atom position in a CBMC/FECBMC regrowth
  ! using the HPS_cation_pi force field in a grand canonical ensemble nucleation simulation.
  ! Evaluates piecewise Hydropathy Scale Lennard-Jones, cation-π interactions, and Debye-Hückel
  ! electrostatic interactions with specified intramolecular atoms. Checks for overlaps based on
  ! minimum distance criterion and returns total energy. Optimized for long linear biomolecules
  ! with sparse non-bonded pairs.
  ! HPS_cation_pi potential:
  !   u_LJ(r) = 4 * eps * [ (sigma/r)^12 - (sigma/r)^6 ]
  !   U_HPS(r) = u_LJ(r) + eps * (1 - lambda), if r <= 2^(1/6) * sigma
  !              lambda * u_LJ(r), if r > 2^(1/6) * sigma
  !   U_cation_pi(r) = 4 * epsLJ * [ (sigma/r)^12 - (sigma/r)^6 ], for specific pairs
  !   U_Electrostatic(r) = (q1 * q2 / r) * exp(-kappa * r)

  use SimParameters, only: dp, kapa               ! Precision and Debye-Hückel screening parameter
  use Coords, only: newMol, SimpleAtomCoords      ! New molecule coordinates and trial position type
  use ForceField, only: atomArray, r_min_tab      ! Atom types and minimum distance table
  use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, &
                                          lambda_tab, cutoffNBsq_tab, rcutElecsq, rcutElec
  implicit none

  ! Input variables
  integer, intent(in) :: nType                    ! Molecule type index
  integer, intent(in) :: nAtom                   ! Atom index within molecule
  type(SimpleAtomCoords), intent(in) :: trialPos  ! Trial coordinates (x, y, z)
  integer, intent(in) :: nIntraUnits             ! Number of intramolecular interacting atoms
  integer, intent(in) :: IntraUnits(:)           ! Indices of interacting atoms in the same molecule

  ! Output variables
  real(dp), intent(out) :: E_Trial               ! Total non-bonded energy (HPS + cation-π + electrostatic)
  logical, intent(inout) :: overlap              ! True if any atom pair distance < r_min_tab

  ! Local variables
  integer :: iPair, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: eps, sig_sq, q, lambda, epsLJ
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS, rcutNBsq, rmax
  real(dp) :: rmin_ij

  ! Initialize energy and overlap variables
  E_HPS = 0.0_dp
  E_Ele = 0.0_dp      
  E_Trial = 0.0_dp
  overlap = .false.

  ! Get atom type for trial position
  atmType1 = atomArray(nType, nAtom)
  ! Loop over intramolecular interacting atoms
  do iPair = 1, nIntraUnits
    jAtom = IntraUnits(iPair)
    ! Calculate distance components
    rx = abs(trialPos%x - newMol%x(jAtom))
    ry = abs(trialPos%y - newMol%y(jAtom))
    rz = abs(trialPos%z - newMol%z(jAtom))
    rmax = max(rx, ry, rz) 
    ! Check electrostatic cutoff
    if (rmax > rcutElec) cycle
    ! Retrieve force field parameters
    atmType2 = atomArray(nType, jAtom)
    eps = eps_tab(atmType1, atmType2)
    epsLJ = epsLJ_tab(atmType1, atmType2)
    q_Nzero = q_Nonzero(atmType1, atmType2)
    sig_sq = sigsq_tab(atmType1, atmType2) 
    rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
    lambda = lambda_tab(atmType1, atmType2) 
    rmin_ij = r_min_tab(atmType1, atmType2)
    ! Check non-bonded cutoff for non-charged pairs
    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif
    ! Compute squared distance
    r = rx*rx + ry*ry + rz*rz
    if (r > rcutElecsq) cycle
    ! Check for atomic overlaps
    if (r < rmin_ij) then
      E_Trial = huge(dp)
      overlap = .true.
      return
    endif
    ! Calculate HPS and LJ energies
    HPS = 0.0_dp
    if (r < rcutNBsq) then
      x = sig_sq/r
      x = x * x * x
      LJ = 4.0_dp * eps * x * (x - 1.0_dp)
      if (x < 0.5_dp) then
        HPS = lambda * LJ
      else
        HPS = LJ + eps * (1.0_dp - lambda)
      endif
      if (epsLJ > 1.0E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
    endif
    E_HPS = E_HPS + HPS
    ! Calculate electrostatic energy
    if (q_Nzero) then   
      q = q_tab(atmType1, atmType2)     
      r = sqrt(r)
      Ele = q * exp(-kapa * r) / r
      E_Ele = E_Ele + Ele       
    endif
  enddo

  ! Compute total non-bonded energy
  E_Trial = E_HPS + E_Ele
      
end subroutine RosenLL_Atom_Intra_New_HPS_cation_pi
!======================================================================================================
pure subroutine RosenLL_Atom_Intra_Old_HPS_cation_pi(nType, nMol, nAtom, trialPos, nIntraUnits, IntraUnits, E_Trial)
  ! Computes intramolecular non-bonded energy for a trial atom position in a CBMC/FECBMC regrowth
  ! using the HPS_cation_pi force field in a grand canonical ensemble nucleation simulation.
  ! Evaluates piecewise Hydropathy Scale Lennard-Jones, cation-π interactions, and Debye-Hückel
  ! electrostatic interactions with specified intramolecular atoms. Checks for overlaps based on
  ! minimum distance criterion and returns total energy. Optimized for long linear biomolecules
  ! with sparse non-bonded pairs.
  ! HPS_cation_pi potential:
  !   u_LJ(r) = 4 * eps * [ (sigma/r)^12 - (sigma/r)^6 ]
  !   U_HPS(r) = u_LJ(r) + eps * (1 - lambda), if r <= 2^(1/6) * sigma
  !              lambda * u_LJ(r), if r > 2^(1/6) * sigma
  !   U_cation_pi(r) = 4 * epsLJ * [ (sigma/r)^12 - (sigma/r)^6 ], for specific pairs
  !   U_Electrostatic(r) = (q1 * q2 / r) * exp(-kappa * r)

  use SimParameters, only: dp, kapa               ! Precision and Debye-Hückel screening parameter
  use Coords, only: MolArray, SimpleAtomCoords   ! Molecule coordinates and trial position type
  use ForceField, only: atomArray, r_min_tab     ! Atom types and minimum distance table
  use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, &
                                          lambda_tab, cutoffNBsq_tab, rcutElecsq, rcutElec
  implicit none

  ! Input variables
  integer, intent(in) :: nType                   ! Molecule type index
  integer, intent(in) :: nMol                    ! Molecule index within type
  integer, intent(in) :: nAtom                   ! Atom index within molecule
  type(SimpleAtomCoords), intent(in) :: trialPos ! Trial coordinates (x, y, z)
  integer, intent(in) :: nIntraUnits            ! Number of intramolecular interacting atoms
  integer, intent(in) :: IntraUnits(:)          ! Indices of interacting atoms in the same molecule

  ! Output variables
  real(dp), intent(out) :: E_Trial              ! Total non-bonded energy (HPS + cation-π + electrostatic)

  ! Local variables
  integer :: iPair, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: eps, sig_sq, q, lambda, epsLJ
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS, rcutNBsq, rmax
  real(dp) :: rmin_ij

  ! Initialize energy variables
  E_HPS = 0.0_dp
  E_Ele = 0.0_dp      
  E_Trial = 0.0_dp

  ! Get atom type for trial position
  atmType1 = atomArray(nType, nAtom)
  ! Loop over intramolecular interacting atoms
  do iPair = 1, nIntraUnits
    jAtom = IntraUnits(iPair)
    ! Calculate distance components
    rx = abs(trialPos%x - MolArray(nType)%mol(nMol)%x(jAtom))
    ry = abs(trialPos%y - MolArray(nType)%mol(nMol)%y(jAtom))
    rz = abs(trialPos%z - MolArray(nType)%mol(nMol)%z(jAtom)) 
    rmax = max(rx, ry, rz) 
    ! Check electrostatic cutoff
    if (rmax > rcutElec) cycle
    ! Retrieve force field parameters
    atmType2 = atomArray(nType, jAtom)
    eps = eps_tab(atmType1, atmType2)
    epsLJ = epsLJ_tab(atmType1, atmType2)
    q_Nzero = q_Nonzero(atmType1, atmType2)
    sig_sq = sigsq_tab(atmType1, atmType2)  
    rcutNBsq = cutoffNBsq_tab(atmType1, atmType2)
    lambda = lambda_tab(atmType1, atmType2) 
    rmin_ij = r_min_tab(atmType1, atmType2)
    ! Check non-bonded cutoff for non-charged pairs
    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif
    ! Compute squared distance
    r = rx*rx + ry*ry + rz*rz
    if (r > rcutElecsq) cycle
    ! Check for atomic overlaps
    if (r < rmin_ij) then
      E_Trial = huge(dp)
      return
    endif
    ! Calculate HPS and LJ energies
    HPS = 0.0_dp
    if (r < rcutNBsq) then
      x = sig_sq/r
      x = x * x * x
      LJ = 4.0_dp * eps * x * (x - 1.0_dp)
      if (x < 0.5_dp) then
        HPS = lambda * LJ
      else
        HPS = LJ + eps * (1.0_dp - lambda)
      endif
      if (epsLJ > 1.0E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
    endif
    E_HPS = E_HPS + HPS
    ! Calculate electrostatic energy
    if (q_Nzero) then     
      q = q_tab(atmType1, atmType2)     
      r = sqrt(r)
      Ele = q * exp(-kapa * r) / r
      E_Ele = E_Ele + Ele       
    endif
  enddo

  ! Compute total non-bonded energy
  E_Trial = E_HPS + E_Ele
      
end subroutine RosenLL_Atom_Intra_Old_HPS_cation_pi
!======================================================================================================
pure subroutine RosenLL_Atom_Intra_Gas_Old_HPS_cation_pi(nType, nAtom, trialPos, nIntraUnits, IntraUnits, E_Trial)
  ! Computes intramolecular non-bonded energy for a trial segment position in the gas phase using
  ! the HPS_cation_pi force field in a grand canonical Monte Carlo nucleation simulation. Called by
  ! LongChain_RosenConfigGen_GasPhase_Reverse within LLAVBMC_EBias_Rosen_In, it calculates the
  ! energy (E_Trial) between a trial segment (nAtom, position trialPos) of molecule type nType and
  ! previously grown segments (IntraUnits) in the gas-phase configuration (gasConfig). Checks for
  ! overlaps using minimum distance (r_min_tab) and applies cutoffs for HPS/cation-pi (cutoffNBsq_tab)
  ! and electrostatics (rcutElecsq, rcutElec). Contributes to the Rosenbluth weight (W_old) for
  ! insertion move acceptance probability.
  ! HPS_cation_pi potential:
  !   u_LJ(r) = 4 * eps * [(sigma/r)^12 - (sigma/r)^6]
  !   U_HPS(r) = lambda * u_LJ(r) if (sigma/r)^6 < 0.5, else u_LJ(r) + eps * (1 - lambda)
  !   U_cation_pi(r) = 4 * epsLJ * [(sigma/r)^12 - (sigma/r)^6], for cation-aromatic pairs
  !   U_Electrostatic(r) = (q1 * q2 / r) * exp(-kappa * r)

  use Coords, only: gasConfig, SimpleAtomCoords
  use ForceField, only: atomArray, r_min_tab
  use ForceFieldPara_HPS_cation_pi, only: eps_tab, epsLJ_tab, q_Nonzero, q_tab, sigsq_tab, &
                                          lambda_tab, cutoffNBsq_tab, rcutElec, rcutElecsq
  use SimParameters, only: kapa
  use Constants, only: dp
  implicit none

  ! Input variables
  integer, intent(in) :: nType                  ! Molecule type (e.g., protein type)
  integer, intent(in) :: nAtom                 ! Segment index (e.g., residue) for trial position
  type(SimpleAtomCoords), intent(in) :: trialPos ! Trial coordinates of the segment
  integer, intent(in) :: nIntraUnits           ! Number of non-bonded intramolecular segments
  integer, intent(in) :: IntraUnits(:)         ! Indices of interacting segments

  ! Output variables
  real(dp), intent(out) :: E_Trial             ! Total intramolecular energy (HPS/cation-pi + electrostatic)

  ! Local variables
  integer :: iPair, jAtom
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: eps, sig_sq, q, lambda, epsLJ
  logical :: q_Nzero
  real(dp) :: HPS, LJ, Ele, x
  real(dp) :: E_Ele, E_HPS, rcutNBsq, rmax
  real(dp) :: rmin_ij

  ! Initialize energy variables
  E_HPS = 0.0_dp
  E_Ele = 0.0_dp      
  E_Trial = 0.0_dp

  ! Get atom type for trial position
  atmType1 = atomArray(nType, nAtom)
  ! Loop over intramolecular interacting segments
  do iPair = 1, nIntraUnits
    jAtom = IntraUnits(iPair)
    ! Calculate distance components
    rx = abs(trialPos%x - gasConfig(nType)%x(jAtom))
    ry = abs(trialPos%y - gasConfig(nType)%y(jAtom))
    rz = abs(trialPos%z - gasConfig(nType)%z(jAtom)) 
    rmax = max(rx, ry, rz) 
    ! Check electrostatic cutoff
    if (rmax > rcutElec) cycle
    ! Retrieve force field parameters
    atmType2 = atomArray(nType, jAtom)
    eps = eps_tab(atmType1, atmType2)
    epsLJ = epsLJ_tab(atmType1, atmType2)
    q_Nzero = q_Nonzero(atmType1, atmType2)
    sig_sq = sigsq_tab(atmType1, atmType2)
    rcutNBsq = cutoffNBsq_tab(atmType1, atmType2) 
    lambda = lambda_tab(atmType1, atmType2) 
    rmin_ij = r_min_tab(atmType1, atmType2)
    ! Check non-bonded cutoff for non-charged pairs
    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif
    ! Compute squared distance
    r = rx*rx + ry*ry + rz*rz
    if (r > rcutElecsq) cycle
    ! Check for atomic overlaps
    if (r < rmin_ij) then
      E_Trial = huge(dp)
      return
    endif
    ! Calculate HPS and LJ energies
    HPS = 0.0_dp
    if (r < rcutNBsq) then
      x = sig_sq/r
      x = x * x * x
      LJ = 4.0_dp * eps * x * (x - 1.0_dp)
      if (x < 0.5_dp) then
        HPS = lambda * LJ
      else
        HPS = LJ + eps * (1.0_dp - lambda)
      endif
      if (epsLJ > 1.0E-5_dp) HPS = HPS + 4.0_dp * epsLJ * x * (x - 1.0_dp)
    endif
    E_HPS = E_HPS + HPS
    ! Calculate electrostatic energy
    if (q_Nzero) then       
      q = q_tab(atmType2, atmType1)     
      r = sqrt(r)
      Ele = q * exp(-kapa * r) / r
      E_Ele = E_Ele + Ele       
    endif
  enddo

  ! Compute total intramolecular energy
  E_Trial = E_HPS + E_Ele
      
end subroutine RosenLL_Atom_Intra_Gas_Old_HPS_cation_pi
!======================================================================================================
pure subroutine Find_InterAtms_HPS_cation_pi(RefAtm, r0, atmType1, nIndx, nInterAtoms, InterAtoms)
  ! Identifies interacting atoms within a cutoff distance from a reference atom (RefAtm) for the
  ! HPS_cation_pi force field in a CBMC or FECBMC move for long linear biomolecules. Computes
  ! distances to all atoms in other molecules, including a buffer distance (r0), and accounts for
  ! non-bonded (HPS and cation-π) and electrostatic interactions based on atom types. Stores
  ! interacting atom indices (molecule type, molecule index, atom index) in InterAtoms and counts
  ! them in nInterAtoms. Excludes the molecule with index nIndx. Called during regrowth to ensure
  ! accurate energy calculations for significantly changing neighbor lists. Supports minimum
  ! distance criterion and distance storage (useDistStore).

  use SimParameters, only: dp, nMolTypes, NPART
  use Coords, only: MolArray, SimpleAtomCoords
  use ForceField, only: nAtoms, atomArray
  use ForceFieldPara_HPS_cation_pi, only: cutoffNB_tab, q_Nonzero, rcutElec
  implicit none

  ! Input variables
  type(SimpleAtomCoords), intent(in) :: RefAtm      ! Position of reference atom
  real(dp), intent(in) :: r0                       ! Buffer distance for cutoff
  integer(kind=atomIntType), intent(in) :: atmType1 ! Atom type of reference atom
  integer, intent(in) :: nIndx                     ! Global index of molecule to exclude

  ! Output variables
  integer, intent(inout) :: nInterAtoms            ! Number of interacting atoms
  integer, intent(out) :: InterAtoms(:, :)         ! Indices (type, mol, atom) of interacting atoms

  ! Local variables
  integer :: jType, jMol, jAtom, jIndx
  integer(kind=atomIntType) :: atmType2
  real(dp) :: rx, ry, rz, r
  logical :: q_Nzero
  real(dp) :: rcutNB, rmax, rcutlong, rcutshort, rcutlongsq, rcutshortsq

  ! Initialize cutoffs and output variables
  rcutlong = r0 + rcutElec
  rcutlongsq = rcutlong * rcutlong
  nInterAtoms = 0
  InterAtoms = 0

  ! Loop over all molecule types
  do jType = 1, nMolTypes
    ! Loop over molecules of current type
    do jMol = 1, NPART(jType) 
      jIndx = molArray(jType)%mol(jMol)%indx     
      ! Skip the molecule with index nIndx
      if (jIndx == nIndx) cycle
      ! Loop over atoms in the molecule
      do jAtom = 1, nAtoms(jType)    
        ! Calculate distance components
        rx = abs(RefAtm%x - MolArray(jType)%mol(jMol)%x(jAtom))
        ry = abs(RefAtm%y - MolArray(jType)%mol(jMol)%y(jAtom))
        rz = abs(RefAtm%z - MolArray(jType)%mol(jMol)%z(jAtom)) 
        rmax = max(rx, ry, rz) 
        ! Check long-range cutoff (electrostatic + buffer)
        if (rmax > rcutlong) cycle
        r = rx*rx + ry*ry + rz*rz
        if (r > rcutlongsq) cycle
        ! Retrieve force field parameters
        atmType2 = atomArray(jType, jAtom)
        q_Nzero = q_Nonzero(atmType2, atmType1)
        rcutNB = cutoffNB_tab(atmType2, atmType1)
        rcutshort = r0 + rcutNB
        rcutshortsq = rcutshort * rcutshort
        ! Check non-bonded cutoff for non-charged pairs
        if (.not. q_Nzero) then
          if (r > rcutshortsq) cycle
        endif
        ! Store interacting atom indices
        nInterAtoms = nInterAtoms + 1
        InterAtoms(nInterAtoms, 1) = jType
        InterAtoms(nInterAtoms, 2) = jMol
        InterAtoms(nInterAtoms, 3) = jAtom
      enddo
    enddo
  enddo

end subroutine Find_InterAtms_HPS_cation_pi
!======================================================================================================
  
      end module 


