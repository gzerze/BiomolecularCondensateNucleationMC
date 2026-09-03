      module Rosenbluth_Functions_Tersoff
      use InterEnergy_Tersoff, only: ep, sig, LJ_Func
      contains

!======================================================================================================
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which regrows an entire molecule for the given trial.
      pure subroutine Rosen_Mol_New_Tersoff(nRosen, nType, included,  E_Trial, overlap)
      use ForceField
      use ForceFieldPara_Tersoff
      use Coords
      use SimParameters
      use InterEnergy_Tersoff, only: Fc_Func
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nRosen      
      
      logical, intent(out) :: overlap
      real(dp), intent(out) :: E_Trial
      
      integer :: jType
      integer :: jMol
      integer(kind=atomIntType) :: atmType1, atmType2      
      integer :: jIndx

      real(dp) :: A, B, c, d, R_eq, D2 
      real(dp) :: lam1, lam2
      real(dp) :: BetaPar, n, h
      real(dp) :: LJ
      real(dp) :: rx, ry, rz, r, rmin_ij

      E_Trial = 0E0
      overlap = .false.

      atmType1 = atomArray(nType, 1)
      A = tersoffData(atmType1)%A
      B = tersoffData(atmType1)%B
      c = tersoffData(atmType1)%c
      d = tersoffData(atmType1)%d
      R_eq = tersoffData(atmType1)%R
      D2 = tersoffData(atmType1)%D2
      BetaPar = tersoffData(atmType1)%beta
      n = tersoffData(atmType1)%n
      h = tersoffData(atmType1)%h
      lam1 = tersoffData(atmType1)%lam1
      lam2 = tersoffData(atmType1)%lam2

      jType = 1
      atmType2 = atomArray(jType, 1)
      rmin_ij = r_min_tab(atmType2, atmType1)
      do jMol = 1,NPART(jType)
        jIndx = molArray(jType)%mol(jMol)%indx              
        if(.not. included(jIndx)) then
          cycle
        endif              

        rx = rosenTrial(nRosen)%x(1) - MolArray(jType)%mol(jMol)%x(1)
        ry = rosenTrial(nRosen)%y(1) - MolArray(jType)%mol(jMol)%y(1)
        rz = rosenTrial(nRosen)%z(1) - MolArray(jType)%mol(jMol)%z(1)
        r = rx*rx + ry*ry + rz*rz
        LJ = LJ_Func(r, ep, sig)
        if(r .lt. rmin_ij) then
          E_Trial = huge(dp)
          overlap = .true.
          return
        endif   
        r = sqrt(r)    
        
        E_Trial = E_Trial + Fc_Func(r, R_eq, D2) * (A*exp(-lam1*r) - 0.95E0_dp*B*exp(-lam2*r))
        E_Trial = E_Trial + LJ
      enddo

      E_Trial = E_Trial * 0.5E0_dp
!      write(*,*) E_Trial
      
      end subroutine 
!======================================================================================================
      pure subroutine Rosen_Mol_Old_Tersoff(mol_x, mol_y, mol_z, nType, included,  E_Trial)
      use ForceField
      use ForceFieldPara_Tersoff
      use Coords
      use SimParameters
      use InterEnergy_Tersoff, only: Fc_Func
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType
      real(dp), intent(in) :: mol_x(:), mol_y(:), mol_z(:)
      real(dp), intent(out) :: E_Trial
      
      integer :: jType
      integer :: jMol
!      integer(kind=atomIntType) :: atmType1, atmType2      
      integer :: atmType1, atmType2      
      integer :: jIndx

      real(dp) :: A, B, c, d, R_eq, D2 
      real(dp) :: lam1, lam2
      real(dp) :: BetaPar, n, h
      real(dp) :: LJ
      real(dp) :: rx, ry, rz, r, rmin_ij


      E_Trial = 0E0
      jType = 1
      atmType1 = atomArray(nType, 1)
      A = tersoffData(atmType1)%A
      B = tersoffData(atmType1)%B
      c = tersoffData(atmType1)%c
      d = tersoffData(atmType1)%d
      R_eq = tersoffData(atmType1)%R
      D2 = tersoffData(atmType1)%D2
      BetaPar = tersoffData(atmType1)%beta
      n = tersoffData(atmType1)%n
      h = tersoffData(atmType1)%h
      lam1 = tersoffData(atmType1)%lam1
      lam2 = tersoffData(atmType1)%lam2

      atmType2 = atomArray(jType, 1)
      rmin_ij = r_min_tab(atmType2, atmType1)
      do jMol = 1,NPART(jType)
        jIndx = molArray(jType)%mol(jMol)%indx              
        if(.not. included(jIndx)) then
          cycle
        endif         
        rx = mol_x(1) - MolArray(jType)%mol(jMol)%x(1)
        ry = mol_y(1) - MolArray(jType)%mol(jMol)%y(1)
        rz = mol_z(1) - MolArray(jType)%mol(jMol)%z(1)
        r = rx*rx + ry*ry + rz*rz
        LJ = LJ_Func(r, ep, sig)
        if(r .lt. rmin_ij) then
          E_Trial = huge(dp)
          return
        endif
        r = sqrt(r)
        E_Trial = E_Trial + Fc_Func(r, R_eq, D2) * (A*exp(-lam1*r) - 0.95E0_dp*B*exp(-lam2*r))
        E_Trial = E_Trial + LJ
      enddo

      E_Trial = E_Trial * 0.5E0_dp
!      write(*,*) E_Trial      
      end subroutine 
!=====================================================================================    
      end module 


