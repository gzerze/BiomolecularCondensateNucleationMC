!======================================================================================================
      module Rosenbluth_Functions_Pedone
      contains
!======================================================================================================
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which regrows an entire molecule for the given trial.
      pure subroutine Rosen_Mol_New_Pedone(nRosen, nType, included,  E_Trial, overlap)
      use ForceField
      use ForceFieldPara_Pedone
      use InterEnergy_Pedone, only: SolventFunction
      use Coords
      use SimParameters
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nRosen      
      
      logical, intent(out) :: overlap
      real(dp), intent(out) :: E_Trial
      
      integer :: jType, jIndx, jMol
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: q, alpha, delta, repul_C, r_eq
      real(dp) :: LJ, Ele, Morse
      real(dp) :: E_Ele, E_LJ, E_Morse
      real(dp) :: rmin_ij
      real(dp) :: born1, born2

      E_Trial = 0E0_dp
      overlap = .false.
      
      E_LJ = 0E0_dp
      E_Ele = 0E0_dp
      E_Morse = 0E0_dp

      atmType1 = atomArray(nType, 1)
      do jType = 1, nMolTypes
        atmType2 = atomArray(jType, 1)
        r_eq = rEq_tab(atmType2, atmType1)
        q = q_tab(atmType2, atmType1)
        alpha = alpha_Tab(atmType2, atmType1)
        delta = D_Tab(atmType2, atmType1)
        repul_C = repul_tab(atmType2, atmType1)
        rmin_ij = r_min_tab(atmType2, atmType1)
        do jMol = 1,NPART(jType)
          jIndx = molArray(jType)%mol(jMol)%indx              
          if(included(jIndx) .eqv. .false.) then
            cycle
          endif              
          rx = rosenTrial(nRosen)%x(1) - MolArray(jType)%mol(jMol)%x(1)
          ry = rosenTrial(nRosen)%y(1) - MolArray(jType)%mol(jMol)%y(1)
          rz = rosenTrial(nRosen)%z(1) - MolArray(jType)%mol(jMol)%z(1)
          r = rx*rx + ry*ry + rz*rz
          if(r .lt. rmin_ij) then
            E_Trial = huge(dp)
            overlap = .true.
            return
          endif           

          if(repul_C .ne. 0E0) then
            LJ = ( 1E0/r )**6
            LJ = repul_C * LJ
            E_LJ = E_LJ + LJ
          endif

          r = sqrt(r)
          Ele = q/r
          if(implcSolvent) then
            born1 = bornRad(atmType1)
            born2 = bornRad(atmType2)
            Ele = Ele + solventFunction(r, q, born1, born2)
          endif
          E_Ele = E_Ele + Ele

          if(delta .ne. 0E0) then
            Morse = 1E0 - exp( -alpha*(r-r_eq) )
            Morse = delta*(Morse*Morse - 1E0)
            E_Morse = E_Morse + Morse
          endif
        enddo
      enddo

      E_Trial = E_LJ + E_Ele + E_Morse
      
      end subroutine 
!======================================================================================================
      pure subroutine Rosen_Mol_Old_Pedone(mol_x, mol_y, mol_z, nType, included,  E_Trial)
      use ForceField
      use ForceFieldPara_Pedone
      use InterEnergy_Pedone, only: SolventFunction
      use Coords
      use SimParameters
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType
      real(dp), intent(in) :: mol_x(:), mol_y(:), mol_z(:)
      
      real(dp), intent(out) :: E_Trial
      
      integer :: jType, jIndx, jMol
      integer(kind=atomIntType) :: atmType1, atmType2
      real(dp) :: rx, ry, rz, r
      real(dp) :: q, alpha, delta, repul_C, r_eq
      real(dp) :: LJ, Ele, Morse
      real(dp) :: E_Ele, E_LJ, E_Morse
      real(dp) :: rmin_ij
      real(dp) :: born1, born2

      E_Trial = 0E0_dp
      E_LJ = 0E0_dp
      E_Ele = 0E0_dp 
      E_Morse = 0E0_dp
      atmType1 = atomArray(nType, 1)
      do jType = 1, nMolTypes
        atmType2 = atomArray(jType, 1)
        r_eq = rEq_tab(atmType1, atmType2)
        q = q_tab(atmType1, atmType2)
        alpha = alpha_Tab(atmType1, atmType2)
        delta = D_Tab(atmType1, atmType2)
        repul_C = repul_tab(atmType1, atmType2)          
        rmin_ij = r_min_tab(atmType1, atmType2) 
        do jMol = 1,NPART(jType)
          jIndx = molArray(jType)%mol(jMol)%indx              
          if(included(jIndx) .eqv. .false.) then
            cycle
          endif         
          rx = mol_x(1) - MolArray(jType)%mol(jMol)%x(1)
          ry = mol_y(1) - MolArray(jType)%mol(jMol)%y(1)
          rz = mol_z(1) - MolArray(jType)%mol(jMol)%z(1)
          r = rx*rx + ry*ry + rz*rz
          if(r .lt. rmin_ij) then
            E_Trial = huge(dp)
            return
          endif
          if(repul_C .ne. 0E0) then
            LJ = ( 1E0/r )**6
            LJ = repul_C * LJ
            E_LJ = E_LJ + LJ
          endif

          r = sqrt(r)
          Ele = q/r
          if(implcSolvent) then
            born1 = bornRad(atmType1)
            born2 = bornRad(atmType2)
            Ele = Ele + solventFunction(r, q, born1, born2)
          endif
          E_Ele = E_Ele + Ele

          if(delta .ne. 0E0) then
            Morse = 1d0-exp(-alpha*(r-r_eq))
            Morse = delta*(Morse*Morse - 1d0)
            E_Morse = E_Morse + Morse
          endif
        enddo
      enddo

     
      E_Trial = E_LJ + E_Ele + E_Morse
      
      end subroutine 
!=====================================================================================    
      end module 


