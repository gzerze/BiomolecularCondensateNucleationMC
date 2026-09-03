      module Rosenbluth_Functions_Mpipi
      use VarPrecision
      contains
!======================================================================================================
pure function WF_Func(r_sq, epAlpha, sig_sq, Mu) result(WF)

  ! Computes the Wang–Frenkel potential function value given squared distance and parameters.
  ! Uses integer Mu for generalization; optimized for Mu = 2 or 3 with precomputed constants.

  implicit none
  real(dp), intent(in) :: r_sq, epAlpha, sig_sq
  integer, intent(in) :: Mu
  real(dp) :: WF, x, y, k

  x = sig_sq / r_sq
  if(Mu == 2) then
    x = x * x
    k = 81.0_dp  ! 9^2
  elseif(Mu == 3) then
    x = x * x * x
    k = 729.0_dp  ! 9^3
  else
    x = x ** Mu
    k = 9.0_dp ** Mu
  endif

  y = k * x - 1.0_dp
  WF = epAlpha * (x - 1.0_dp) * y * y

end function WF_Func
!====================================================================================== 
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which regrows an entire molecule for the given trial.
      pure subroutine Rosen_Mol_New_Mpipi(nRosen, nType, included,  E_Trial, overlap)
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, rcutElec, rcutElecsq
      use Coords, only: rosenTrial, MolArray
      use SimParameters, only: nMolTypes, NPART, kapa
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nRosen      
      
      logical, intent(out) :: overlap
      real(dp), intent(out) :: E_Trial
      
      integer :: iAtom, jType, jIndx, jMol, jAtom, Mu
      integer(kind = atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: epAlpha,sig_sq,q
      logical :: q_Nzero
      real(dp) :: WF, Ele
      real(dp) :: E_Ele,E_WF, rcutNBsq, rmax
      real(dp) :: rmin_ij

      E_Trial = 0E0
      overlap = .false.
      
      E_WF = 0E0
      E_Ele = 0E0      

      do iAtom = 1,nAtoms(nType)
        atmType1 = atomArray(nType, iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType, jAtom)
            epAlpha = epAlpha_tab(atmType2, atmType1)
            q_Nzero = q_Nonzero(atmType2, atmType1)
            sig_sq = sigsq_tab(atmType2, atmType1)  
            rcutNBsq = 9E0 * sig_sq
            Mu = Mu_tab(atmType2, atmType1)   
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

              WF = 0E0
              if (r .lt. rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
              E_WF = E_WF + WF
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
      E_Trial = E_WF + E_Ele
      
      end subroutine 
!======================================================================================================
      pure subroutine Rosen_Mol_Old_Mpipi(mol_x, mol_y, mol_z, nType, included,  E_Trial)
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, rcutElec, rcutElecsq
      use Coords, only: MolArray
      use SimParameters, only: NPART, nMolTypes, kapa
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType
      real(dp), intent(in) :: mol_x(:), mol_y(:), mol_z(:)
      real(dp), intent(out) :: E_Trial
      
      integer :: iAtom, jType, jIndx, jMol, jAtom, Mu
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: epAlpha,sig_sq,q
      logical :: q_Nzero
      real(dp) :: WF, Ele
      real(dp) :: E_Ele,E_WF, rcutNBsq, rmax
      real(dp) :: rmin_ij

      E_Trial = 0E0
      E_WF = 0E0
      E_Ele = 0E0      


      do iAtom = 1,nAtoms(nType)
        atmType1 = atomArray(nType, iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType, jAtom)
            epAlpha = epAlpha_tab(atmType2, atmType1)
            q_Nzero = q_Nonzero(atmType2, atmType1)
            sig_sq = sigsq_tab(atmType2, atmType1)  
            rcutNBsq = 9E0 * sig_sq
            Mu = Mu_tab(atmType2, atmType1)  
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
              WF = 0E0
              if (r .lt. rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
              E_WF = E_WF + WF
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
     
      E_Trial = E_WF + E_Ele
      
      end subroutine 
!======================================================================================================
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which each atom is regrown sequentially for the given trial.
      pure subroutine Rosen_Atom_New_Mpipi(nType, nAtom, trialPos, included,  E_Trial, overlap)
      use Coords, only: MolArray, SimpleAtomCoords
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, rcutElec, rcutElecsq
      use SimParameters, only: nMolTypes, NPART, kapa
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      
      logical, intent(inout) :: overlap
      real(dp), intent(out) :: E_Trial
      
      integer :: jType, jIndx, jMol, jAtom, Mu
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: epAlpha,sig_sq,q
      logical :: q_Nzero
      real(dp) :: WF, Ele
      real(dp) :: E_Ele,E_WF, rcutNBsq, rmax
      real(dp) :: rmin_ij

      
      E_WF = 0E0
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
            epAlpha = epAlpha_tab(atmType2, atmType1)
            q_Nzero = q_Nonzero(atmType2, atmType1)
            sig_sq = sigsq_tab(atmType2, atmType1)  
            rcutNBsq = 9E0 * sig_sq
            Mu = Mu_tab(atmType2, atmType1)  
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
            WF = 0E0
            if (r .lt. rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
            E_WF = E_WF + WF
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
     
      E_Trial = E_WF + E_Ele
      
      end subroutine 
!======================================================================================================
      pure subroutine Rosen_Atom_Old_Mpipi(nType, nMol, nAtom, included, E_Trial)
      use Coords, only: MolArray
      use ForceField, only: nAtoms, r_min_tab, atomArray
      use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, rcutElec, rcutElecsq
      use SimParameters, only: nMolTypes, NPART, kapa
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nAtom, nMol
      real(dp), intent(out) :: E_Trial
      
      integer :: jType, jIndx, jMol, jAtom, Mu
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: epAlpha,sig_sq,q
      logical :: q_Nzero
      real(dp) :: WF, Ele
      real(dp) :: E_Ele,E_WF, rcutNBsq, rmax
      real(dp) :: rmin_ij

      
      E_WF = 0E0
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
            epAlpha = epAlpha_tab(atmType2, atmType1)
            q_Nzero = q_Nonzero(atmType2, atmType1)
            sig_sq = sigsq_tab(atmType2, atmType1) 
            rcutNBsq = 9E0 * sig_sq 
            Mu = Mu_tab(atmType2, atmType1) 
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
 

            WF = 0E0
            if (r .lt. rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
            E_WF = E_WF + WF
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
     
      E_Trial = E_WF + E_Ele
      
      end subroutine 
!======================================================================================================
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which each atom is regrown sequentially for the given trial.
      pure subroutine Rosen_Atom_Intra_New_Mpipi(nType, nAtom, trialPos, regrown, E_Trial, overlap)
      use Coords, only: newMol, SimpleAtomCoords
      use ForceField, only: atomArray, nonBondArray, nIntraNonBond, r_min_tab
      use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, rcutElec, rcutElecsq
      use SimParameters, only: kapa
      implicit none

      logical, intent(in) :: regrown(:)
      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      logical, intent(inout) :: overlap
      real(dp), intent(out) :: E_Trial

      integer :: temp, iPair, iAtom, jAtom, Mu
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: epAlpha,sig_sq,q
      logical :: q_Nzero
      real(dp) :: WF, Ele
      real(dp) :: E_Ele,E_WF, rcutNBsq, rmax
      real(dp) :: rmin_ij

      E_WF = 0E0
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
        epAlpha = epAlpha_tab(atmType1,atmType2)
        q_Nzero = q_Nonzero(atmType1,atmType2)
        sig_sq = sigsq_tab(atmType1,atmType2) 
        rcutNBsq = 9d0 * sig_sq 
        Mu = Mu_tab(atmType1,atmType2) 
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
        WF = 0E0
        if (r .lt. rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
        E_WF = E_WF + WF
        if(q_Nzero) then   
          q = q_tab(atmType1,atmType2)     
          r = sqrt(r)
          Ele = q * exp(-kapa * r)/r
          E_Ele = E_Ele + Ele       
        endif
      enddo
      E_Trial = E_WF + E_Ele

      end subroutine 
!======================================================================================================
      pure subroutine Rosen_Atom_Intra_Old_Mpipi(nType, nMol, nAtom, trialPos, regrown, E_Trial)
      use Coords, only: MolArray, SimpleAtomCoords
      use ForceField, only: atomArray, nonBondArray, nIntraNonBond, r_min_tab
      use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, rcutElec, rcutElecsq
      use SimParameters, only: kapa
      implicit none

      logical, intent(in) :: regrown(:)
      integer, intent(in) :: nType, nMol, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      real(dp), intent(out) :: E_Trial

      integer :: temp, iPair, iAtom, jAtom, Mu
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: epAlpha,sig_sq,q
      logical :: q_Nzero
      real(dp) :: WF, Ele
      real(dp) :: E_Ele,E_WF, rcutNBsq, rmax
      real(dp) :: rmin_ij

      E_WF = 0E0
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
        epAlpha = epAlpha_tab(atmType1,atmType2)
        q_Nzero = q_Nonzero(atmType1,atmType2)
        sig_sq = sigsq_tab(atmType1,atmType2) 
        rcutNBsq = 9d0 * sig_sq 
        Mu = Mu_tab(atmType1,atmType2) 
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
        WF = 0E0
        if (r .lt. rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
        E_WF = E_WF + WF
        if(q_Nzero) then   
          q = q_tab(atmType1,atmType2)     
          r = sqrt(r)
          Ele = q * exp(-kapa * r)/r
          E_Ele = E_Ele + Ele       
        endif
      enddo
      E_Trial = E_WF + E_Ele

      end subroutine 
!======================================================================================================
      pure subroutine Rosen_Atom_Intra_Gas_Old_Mpipi(nType, nAtom, trialPos, regrown, E_Trial)
      use Coords, only: gasConfig, SimpleAtomCoords
      use ForceField, only: atomArray, nonBondArray, nIntraNonBond, r_min_tab
      use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, rcutElec, rcutElecsq
      use SimParameters, only: kapa
      implicit none

      logical, intent(in) :: regrown(:)
      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      real(dp), intent(out) :: E_Trial

      integer :: temp, iPair, iAtom, jAtom, Mu
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: epAlpha,sig_sq,q
      logical :: q_Nzero
      real(dp) :: WF, Ele
      real(dp) :: E_Ele,E_WF, rcutNBsq, rmax
      real(dp) :: rmin_ij

      E_WF = 0E0
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
        epAlpha = epAlpha_tab(atmType1,atmType2)
        q_Nzero = q_Nonzero(atmType1,atmType2)
        sig_sq = sigsq_tab(atmType1,atmType2) 
        rcutNBsq = 9d0 * sig_sq 
        Mu = Mu_tab(atmType1,atmType2) 
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
        WF = 0E0
        if (r .lt. rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
        E_WF = E_WF + WF
        if(q_Nzero) then   
          q = q_tab(atmType1,atmType2)     
          r = sqrt(r)
          Ele = q * exp(-kapa * r)/r
          E_Ele = E_Ele + Ele       
        endif
      enddo
      E_Trial = E_WF + E_Ele

      end subroutine 
!======================================================================================================
! Computes intermolecular non-bonded energy for a trial atom position in a CBMC/FECBMC regrowth
! using the Mpipi force field in a grand canonical ensemble nucleation simulation. Evaluates
! Wang–Frankel and Debye-Hückel electrostatic interactions with specified interacting atoms,
! checks for overlaps based on minimum distance criterion, and returns total energy. Optimized
! for long linear biomolecules with sparse neighbor lists.
pure subroutine RosenLL_Atom_New_Mpipi(nType, nAtom, trialPos, nInterAtoms, InterAtoms, E_Trial, overlap)

  ! Computes intermolecular non-bonded energy for a trial atom position in a CBMC/FECBMC regrowth
  ! using the Mpipi force field in a grand canonical ensemble nucleation simulation. Evaluates
  ! Wang–Frankel and Debye-Hückel electrostatic interactions with specified interacting atoms,
  ! checks for overlaps based on minimum distance criterion, and returns total energy. Optimized
  ! for long linear biomolecules with sparse neighbor lists.

  use VarPrecision, only: dp, atomIntType
  use SimParameters, only: kapa
  use Coords, only: MolArray, SimpleAtomCoords
  use ForceField, only: r_min_tab, atomArray
  use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, &
                                  rcutElec, rcutElecsq

  implicit none

  integer, intent(in) :: nType, nAtom
  type(SimpleAtomCoords), intent(in) :: trialPos
  integer, intent(in) :: nInterAtoms
  integer, intent(in) :: InterAtoms(:,:)
  logical, intent(inout) :: overlap
  real(dp), intent(out) :: E_Trial

  integer :: iInterAtoms, jType, jMol, jAtom, Mu
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: epAlpha, sig_sq, q
  logical :: q_Nzero
  real(dp) :: WF, Ele
  real(dp) :: E_Ele, E_WF, rcutNBsq, rmax, rmin_ij

  E_WF = 0.0_dp
  E_Ele = 0.0_dp
  E_Trial = 0.0_dp
  overlap = .false.

  atmType1 = atomArray(nType, nAtom)

  do iInterAtoms = 1, nInterAtoms
    jType = InterAtoms(iInterAtoms, 1)
    jMol  = InterAtoms(iInterAtoms, 2)
    jAtom = InterAtoms(iInterAtoms, 3)

    rx = abs(trialPos%x - MolArray(jType)%mol(jMol)%x(jAtom))
    ry = abs(trialPos%y - MolArray(jType)%mol(jMol)%y(jAtom))
    rz = abs(trialPos%z - MolArray(jType)%mol(jMol)%z(jAtom))
    rmax = max(rx, ry, rz)
    if (rmax > rcutElec) cycle

    atmType2 = atomArray(jType, jAtom)
    rmin_ij  = r_min_tab(atmType2, atmType1)
    epAlpha  = epAlpha_tab(atmType2, atmType1)
    q_Nzero  = q_Nonzero(atmType2, atmType1)
    sig_sq   = sigsq_tab(atmType2, atmType1)
    rcutNBsq = 9.0_dp * sig_sq
    Mu       = Mu_tab(atmType2, atmType1)

    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif

    r = rx*rx + ry*ry + rz*rz
    if (r > rcutElecsq) cycle

    if (r < rmin_ij) then
      E_Trial = huge(dp)
      overlap = .true.
      return
    endif

    WF = 0.0_dp
    if (r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
    E_WF = E_WF + WF

    if (q_Nzero) then
      q = q_tab(atmType2, atmType1)
      r = sqrt(r)
      Ele = q * exp(-kapa * r) / r
      E_Ele = E_Ele + Ele
    endif
  enddo

  E_Trial = E_WF + E_Ele

end subroutine RosenLL_Atom_New_Mpipi
!======================================================================================================
! Computes intermolecular non-bonded energy for an existing atom in a CBMC/FECBMC regrowth
! using the Mpipi force field in a grand canonical ensemble nucleation simulation. Evaluates
! Wang–Frankel and Debye-Hückel electrostatic interactions with specified interacting atoms
! and returns total energy. Optimized for long linear biomolecules with sparse neighbor lists.
pure subroutine RosenLL_Atom_Old_Mpipi(nType, nMol, nAtom, nInterAtoms, InterAtoms, E_Trial)

  ! Computes intermolecular non-bonded energy for an existing atom in a CBMC/FECBMC regrowth
  ! using the Mpipi force field in a grand canonical ensemble nucleation simulation. Evaluates
  ! Wang–Frankel and Debye-Hückel electrostatic interactions with specified interacting atoms
  ! and returns total energy. Optimized for long linear biomolecules with sparse neighbor lists.

  use SimParameters, only: kapa
  use Coords, only: MolArray
  use ForceField, only: r_min_tab, atomArray
  use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, &
                                  rcutElec, rcutElecsq
  use VarPrecision, only: dp, atomIntType

  implicit none

  integer, intent(in) :: nType, nMol, nAtom
  integer, intent(in) :: nInterAtoms
  integer, intent(in) :: InterAtoms(:,:)
  real(dp), intent(out) :: E_Trial

  integer :: iInterAtoms, jType, jMol, jAtom, Mu
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: epAlpha, sig_sq, q
  logical :: q_Nzero
  real(dp) :: WF, Ele
  real(dp) :: E_Ele, E_WF, rcutNBsq, rmax, rmin_ij

  E_WF = 0.0_dp
  E_Ele = 0.0_dp
  E_Trial = 0.0_dp

  atmType1 = atomArray(nType, nAtom)

  do iInterAtoms = 1, nInterAtoms
    jType = InterAtoms(iInterAtoms, 1)
    jMol  = InterAtoms(iInterAtoms, 2)
    jAtom = InterAtoms(iInterAtoms, 3)

    rx = abs(MolArray(nType)%mol(nMol)%x(nAtom) - MolArray(jType)%mol(jMol)%x(jAtom))
    ry = abs(MolArray(nType)%mol(nMol)%y(nAtom) - MolArray(jType)%mol(jMol)%y(jAtom))
    rz = abs(MolArray(nType)%mol(nMol)%z(nAtom) - MolArray(jType)%mol(jMol)%z(jAtom))
    rmax = max(rx, ry, rz)
    if (rmax > rcutElec) cycle

    atmType2 = atomArray(jType, jAtom)
    rmin_ij  = r_min_tab(atmType2, atmType1)
    epAlpha  = epAlpha_tab(atmType2, atmType1)
    q_Nzero  = q_Nonzero(atmType2, atmType1)
    sig_sq   = sigsq_tab(atmType2, atmType1)
    rcutNBsq = 9.0_dp * sig_sq
    Mu       = Mu_tab(atmType2, atmType1)

    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif

    r = rx*rx + ry*ry + rz*rz
    if (r > rcutElecsq) cycle

    if (r < rmin_ij) then
      E_Trial = huge(dp)
      return
    endif

    WF = 0.0_dp
    if (r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
    E_WF = E_WF + WF

    if (q_Nzero) then
      q = q_tab(atmType2, atmType1)
      r = sqrt(r)
      Ele = q * exp(-kapa * r) / r
      E_Ele = E_Ele + Ele
    endif
  enddo

  E_Trial = E_WF + E_Ele

end subroutine RosenLL_Atom_Old_Mpipi
!======================================================================================================
! Computes intramolecular non-bonded energy for a trial atom position in a CBMC/FECBMC regrowth
! using the Mpipi force field in a grand canonical ensemble nucleation simulation. Evaluates
! Wang–Frankel and Debye-Hückel electrostatic interactions with specified intramolecular atoms,
! checks for overlaps based on minimum distance criterion, and returns total energy. Optimized
! for long linear biomolecules with sparse non-bonded pairs.
pure subroutine RosenLL_Atom_Intra_New_Mpipi(nType, nAtom, trialPos, nIntraUnits, IntraUnits, E_Trial, overlap)

  ! Computes intramolecular non-bonded energy for a trial atom position in a CBMC/FECBMC regrowth
  ! using the Mpipi force field in a grand canonical ensemble nucleation simulation. Evaluates
  ! Wang–Frankel and Debye-Hückel electrostatic interactions with specified intramolecular atoms,
  ! checks for overlaps based on minimum distance criterion, and returns total energy. Optimized
  ! for long linear biomolecules with sparse non-bonded pairs.

  use Coords, only: newMol, SimpleAtomCoords
  use ForceField, only: atomArray, r_min_tab
  use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, rcutElec, rcutElecsq
  use SimParameters, only: kapa
  use VarPrecision, only: dp, atomIntType

  implicit none

  integer, intent(in) :: nType, nAtom, nIntraUnits
  integer, intent(in) :: IntraUnits(:)
  type(SimpleAtomCoords), intent(in) :: trialPos
  real(dp), intent(out) :: E_Trial
  logical, intent(inout) :: overlap

  integer :: iPair, jAtom, Mu
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: epAlpha, sig_sq, q
  logical :: q_Nzero
  real(dp) :: WF, Ele
  real(dp) :: E_Ele, E_WF, rcutNBsq, rmax, rmin_ij

  E_WF = 0.0_dp
  E_Ele = 0.0_dp
  E_Trial = 0.0_dp
  overlap = .false.

  atmType1 = atomArray(nType, nAtom)

  do iPair = 1, nIntraUnits
    jAtom = IntraUnits(iPair)

    rx = abs(trialPos%x - newMol%x(jAtom))
    ry = abs(trialPos%y - newMol%y(jAtom))
    rz = abs(trialPos%z - newMol%z(jAtom))
    rmax = max(rx, ry, rz)
    if (rmax > rcutElec) cycle

    atmType2 = atomArray(nType, jAtom)
    epAlpha  = epAlpha_tab(atmType1, atmType2)
    sig_sq   = sigsq_tab(atmType1, atmType2)
    q_Nzero  = q_Nonzero(atmType1, atmType2)
    Mu       = Mu_tab(atmType1, atmType2)
    rmin_ij  = r_min_tab(atmType1, atmType2)
    rcutNBsq = 9.0_dp * sig_sq

    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif

    r = rx*rx + ry*ry + rz*rz
    if (r > rcutElecsq) cycle

    if (r < rmin_ij) then
      E_Trial = huge(dp)
      overlap = .true.
      return
    endif

    WF = 0.0_dp
    if (r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
    E_WF = E_WF + WF

    if (q_Nzero) then
      q = q_tab(atmType1, atmType2)
      r = sqrt(r)
      Ele = q * exp(-kapa * r) / r
      E_Ele = E_Ele + Ele
    endif

  enddo

  E_Trial = E_WF + E_Ele

end subroutine RosenLL_Atom_Intra_New_Mpipi
!======================================================================================================
! Computes intramolecular non-bonded energy for a trial atom position in a CBMC/FECBMC regrowth
! using the Mpipi force field in a grand canonical ensemble nucleation simulation. Evaluates
! Wang–Frankel and Debye-Hückel electrostatic interactions with specified intramolecular atoms
! and returns total energy. Optimized for long linear biomolecules with sparse non-bonded pairs.
pure subroutine RosenLL_Atom_Intra_Old_Mpipi(nType, nMol, nAtom, trialPos, nIntraUnits, IntraUnits, E_Trial)

  ! Computes intramolecular non-bonded energy for a trial atom position in a CBMC/FECBMC regrowth
  ! using the Mpipi force field in a grand canonical ensemble nucleation simulation. Evaluates
  ! Wang–Frankel and Debye-Hückel electrostatic interactions with specified intramolecular atoms
  ! and returns total energy. Optimized for long linear biomolecules with sparse non-bonded pairs.

  use Coords, only: MolArray, SimpleAtomCoords
  use ForceField, only: atomArray, r_min_tab
  use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, rcutElec, rcutElecsq
  use SimParameters, only: kapa
  use VarPrecision, only: dp, atomIntType

  implicit none

  integer, intent(in) :: nType, nMol, nAtom, nIntraUnits
  integer, intent(in) :: IntraUnits(:)
  type(SimpleAtomCoords), intent(in) :: trialPos
  real(dp), intent(out) :: E_Trial

  integer :: iPair, jAtom, Mu
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: epAlpha, sig_sq, q
  logical :: q_Nzero
  real(dp) :: WF, Ele
  real(dp) :: E_Ele, E_WF, rcutNBsq, rmax, rmin_ij

  E_WF    = 0.0_dp
  E_Ele   = 0.0_dp
  E_Trial = 0.0_dp

  atmType1 = atomArray(nType, nAtom)

  do iPair = 1, nIntraUnits
    jAtom = IntraUnits(iPair)

    rx = abs(trialPos%x - MolArray(nType)%mol(nMol)%x(jAtom))
    ry = abs(trialPos%y - MolArray(nType)%mol(nMol)%y(jAtom))
    rz = abs(trialPos%z - MolArray(nType)%mol(nMol)%z(jAtom))
    rmax = max(rx, ry, rz)
    if (rmax > rcutElec) cycle

    atmType2 = atomArray(nType, jAtom)
    epAlpha  = epAlpha_tab(atmType1, atmType2)
    sig_sq   = sigsq_tab(atmType1, atmType2)
    q_Nzero  = q_Nonzero(atmType1, atmType2)
    Mu       = Mu_tab(atmType1, atmType2)
    rmin_ij  = r_min_tab(atmType1, atmType2)
    rcutNBsq = 9.0_dp * sig_sq

    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif

    r = rx*rx + ry*ry + rz*rz
    if (r > rcutElecsq) cycle

    if (r < rmin_ij) then
      E_Trial = huge(dp)
      return
    endif

    WF = 0.0_dp
    if (r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
    E_WF = E_WF + WF

    if (q_Nzero) then
      q = q_tab(atmType1, atmType2)
      r = sqrt(r)
      Ele = q * exp(-kapa * r) / r
      E_Ele = E_Ele + Ele
    endif

  enddo

  E_Trial = E_WF + E_Ele

end subroutine RosenLL_Atom_Intra_Old_Mpipi
!======================================================================================================
! Computes the non-bonded intramolecular energy for a trial segment position in the gas phase using 
! the Mpipi force field, as part of the reverse move (new to old state) in a grand canonical Monte 
! Carlo nucleation simulation. Called by LongChain_RosenConfigGen_GasPhase_Reverse within 
! LLAVBMC_EBias_Rosen_In, it calculates the energy (E_Trial) between a trial segment (nAtom, 
! position trialPos) of molecule type nType and previously grown segments (IntraUnits) in the 
! gas-phase configuration (gasConfig). The Mpipi force field includes Wang–Frankel (WF) 
! interactions (via WF_Func, with parameters epAlpha, sig_sq, Mu) and screened Coulombic 
! interactions (exp(-kapa * r)/r). Energies contribute to the Rosenbluth weight (W_old) for the 
! insertion move’s acceptance probability. Checks for overlaps using minimum distance (r_min_tab) 
! and applies cutoffs for WF (rcutNBsq = 9 * sig_sq) and electrostatic (rcutElecsq) interactions. 
! Compatible with the minimum-distance cluster criterion (minDistCriteria) and the Mpipi force 
! field (extensible to LJ_Q, HPS_single, HPS_piecewise, HPS_cation_pi).
pure subroutine RosenLL_Atom_Intra_Gas_Old_Mpipi(nType, nAtom, trialPos, nIntraUnits, IntraUnits, E_Trial)

  ! Computes the non-bonded intramolecular energy for a trial segment position in the gas phase using 
  ! the Mpipi force field, as part of the reverse move (new to old state) in a grand canonical Monte 
  ! Carlo nucleation simulation. Called by LongChain_RosenConfigGen_GasPhase_Reverse within 
  ! LLAVBMC_EBias_Rosen_In, it calculates the energy (E_Trial) between a trial segment (nAtom, 
  ! position trialPos) of molecule type nType and previously grown segments (IntraUnits) in the 
  ! gas-phase configuration (gasConfig). The Mpipi force field includes Wang–Frankel (WF) 
  ! interactions (via WF_Func, with parameters epAlpha, sig_sq, Mu) and screened Coulombic 
  ! interactions (exp(-kapa * r)/r). Energies contribute to the Rosenbluth weight (W_old) for the 
  ! insertion move’s acceptance probability. Checks for overlaps using minimum distance (r_min_tab) 
  ! and applies cutoffs for WF (rcutNBsq = 9 * sig_sq) and electrostatic (rcutElecsq) interactions. 
  ! Compatible with the minimum-distance cluster criterion and the Mpipi force field.

  use Coords, only: gasConfig, SimpleAtomCoords
  use ForceField, only: atomArray, r_min_tab
  use ForceFieldPara_Mpipi, only: q_tab, epAlpha_tab, sigsq_tab, q_Nonzero, Mu_tab, rcutElec, rcutElecsq
  use SimParameters, only: kapa
  use VarPrecision, only: dp, atomIntType

  implicit none

  integer, intent(in) :: nType, nAtom, nIntraUnits
  integer, intent(in) :: IntraUnits(:)
  type(SimpleAtomCoords), intent(in) :: trialPos
  real(dp), intent(out) :: E_Trial

  integer :: iPair, jAtom, Mu
  integer(kind=atomIntType) :: atmType1, atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: epAlpha, sig_sq, q
  logical :: q_Nzero
  real(dp) :: WF, Ele
  real(dp) :: E_Ele, E_WF, rcutNBsq, rmax, rmin_ij

  E_WF    = 0.0_dp
  E_Ele   = 0.0_dp
  E_Trial = 0.0_dp

  atmType1 = atomArray(nType, nAtom)

  do iPair = 1, nIntraUnits
    jAtom = IntraUnits(iPair)

    rx = abs(trialPos%x - gasConfig(nType)%x(jAtom))
    ry = abs(trialPos%y - gasConfig(nType)%y(jAtom))
    rz = abs(trialPos%z - gasConfig(nType)%z(jAtom))
    rmax = max(rx, ry, rz)
    if (rmax > rcutElec) cycle

    atmType2 = atomArray(nType, jAtom)
    epAlpha  = epAlpha_tab(atmType1, atmType2)
    sig_sq   = sigsq_tab(atmType1, atmType2)
    q_Nzero  = q_Nonzero(atmType1, atmType2)
    Mu       = Mu_tab(atmType1, atmType2)
    rmin_ij  = r_min_tab(atmType1, atmType2)
    rcutNBsq = 9.0_dp * sig_sq

    if (.not. q_Nzero) then
      rmax = rmax * rmax
      if (rmax > rcutNBsq) cycle
    endif

    r = rx*rx + ry*ry + rz*rz
    if (r > rcutElecsq) cycle

    if (r < rmin_ij) then
      E_Trial = huge(dp)
      return
    endif

    WF = 0.0_dp
    if (r < rcutNBsq) WF = WF_Func(r, epAlpha, sig_sq, Mu)
    E_WF = E_WF + WF

    if (q_Nzero) then
      q = q_tab(atmType2, atmType1)
      r = sqrt(r)
      Ele = q * exp(-kapa * r) / r
      E_Ele = E_Ele + Ele
    endif
  enddo

  E_Trial = E_WF + E_Ele

end subroutine RosenLL_Atom_Intra_Gas_Old_Mpipi
!======================================================================================================
! Identifies interacting atoms within a cutoff distance from a reference atom (RefAtm) for the Mpipi
! force field in a Configurational Bias Monte Carlo (CBMC) or Fixed Endpoint CBMC (FECBMC) move for
! long linear biomolecules. Computes distances to all atoms in other molecules, including a buffer
! distance (r0), and accounts for non-bonded (Lennard-Jones) and electrostatic interactions based on
! atom types. Stores interacting atom indices (molecule type, molecule index, atom index) in
! InterAtoms and counts them in nInterAtoms. Excludes the molecule with index nIndx. Called during
! regrowth to ensure accurate energy calculations for significantly changing neighbor lists.
! Supports minimum distance criterion and distance storage (useDistStore).
subroutine Find_InterAtms_Mpipi(RefAtm, r0, atmType1, nIndx, nInterAtoms, InterAtoms)

  ! Identifies interacting atoms within a cutoff distance from a reference atom (RefAtm) for the Mpipi
  ! force field in a Configurational Bias Monte Carlo (CBMC) or Fixed Endpoint CBMC (FECBMC) move for
  ! long linear biomolecules. Computes distances to all atoms in other molecules, including a buffer
  ! distance (r0), and accounts for non-bonded (Lennard-Jones) and electrostatic interactions based on
  ! atom types. Stores interacting atom indices (molecule type, molecule index, atom index) in
  ! InterAtoms and counts them in nInterAtoms. Excludes the molecule with index nIndx. Called during
  ! regrowth to ensure accurate energy calculations for significantly changing neighbor lists.
  ! Supports minimum distance criterion and distance storage (useDistStore).

  use Coords, only: MolArray, SimpleAtomCoords
  use ForceField, only: nAtoms, atomArray
  use ForceFieldPara_Mpipi, only: sig_tab, q_Nonzero, rcutElec
  use SimParameters, only: nMolTypes, NPART
  use VarPrecision, only: dp, atomIntType

  implicit none

  type(SimpleAtomCoords), intent(in) :: RefAtm
  real(dp), intent(in) :: r0
  integer(kind=atomIntType), intent(in) :: atmType1
  integer, intent(in) :: nIndx
  integer, intent(inout) :: nInterAtoms
  integer, intent(out) :: InterAtoms(:,:)

  integer :: jType, jMol, jAtom, jIndx
  integer(kind=atomIntType) :: atmType2
  real(dp) :: rx, ry, rz, r
  real(dp) :: rcutNB
  real(dp) :: rmax, rcutlong, rcutshort, rcutlongsq, rcutshortsq
  logical :: q_Nzero

  rcutlong   = r0 + rcutElec
  rcutlongsq = rcutlong * rcutlong

  nInterAtoms = 0
  InterAtoms = 0

  do jType = 1, nMolTypes
    do jMol = 1, NPART(jType)
      jIndx = MolArray(jType)%mol(jMol)%indx
      if (jIndx == nIndx) cycle

      do jAtom = 1, nAtoms(jType)

        rx = abs(RefAtm%x - MolArray(jType)%mol(jMol)%x(jAtom))
        ry = abs(RefAtm%y - MolArray(jType)%mol(jMol)%y(jAtom))
        rz = abs(RefAtm%z - MolArray(jType)%mol(jMol)%z(jAtom))
        rmax = max(rx, ry, rz)
        if (rmax > rcutlong) cycle

        r = rx*rx + ry*ry + rz*rz
        if (r > rcutlongsq) cycle

        atmType2 = atomArray(jType, jAtom)
        q_Nzero  = q_Nonzero(atmType2, atmType1)
        rcutNB   = 3.0_dp * sig_tab(atmType2, atmType1)
        rcutshort = r0 + rcutNB
        rcutshortsq = rcutshort * rcutshort

        if (.not. q_Nzero) then
          if (r > rcutshortsq) cycle
        endif

        nInterAtoms = nInterAtoms + 1
        InterAtoms(nInterAtoms, 1) = jType
        InterAtoms(nInterAtoms, 2) = jMol
        InterAtoms(nInterAtoms, 3) = jAtom

      enddo
    enddo
  enddo

end subroutine Find_InterAtms_Mpipi
!======================================================================================================
  
      end module 


