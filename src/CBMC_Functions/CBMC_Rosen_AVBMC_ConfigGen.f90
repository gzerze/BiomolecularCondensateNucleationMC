      module AVBMC_CBMC
      contains
!=======================================================================
!     The purpose of this subroutine is generate a set of trial configurations that will
!     used for the Aggregation-Volume-Bias insertion move. This subroutine is intended to generate
!     configurations for perfectly ridgid molecules. Once a set of configurations is generated
!     one of these configurations will be randomly chosen based on the Rosenbluth weight for each
!     configuraiton. 
      subroutine Ridgid_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions
      use EnergyPressurePointers, only: Rosen_Mol_New
      use CBMC_Variables
      use CBMC_Utility
      use CoordinateFunctions
      implicit none

      integer, intent(in) :: nType, nTarget, nTargType, nIndx
      logical, intent(out) :: rejMove      
      real(dp), intent(out) :: rosenRatio
      logical, intent(in) :: isIncluded(:)

      logical :: overlap(1:maxRosenTrial)
      integer :: atmType1, atmType2
      integer :: i, iRosen, nSel, nTargetMol
      real(dp) :: E_Trial(1:maxRosenTrial)      
      real(dp) :: grnd,rotang
      real(dp) :: c_term,s_term
      real(dp) :: dx, dy, dz
      real(dp) :: r, x1, y1, z1      
      real(dp) :: x_shift,y_shift,z_shift
      real(dp) :: x_rid_cm, y_rid_cm,z_rid_cm
      real(dp) :: E_Max, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: ranNum, sumInt, rmin_ij
      real(dp) :: rx, ry, rz
      
      newMol%molType = nType      
!      call Rosen_CreateSubset(nTarget, isIncluded)
      E_Trial = 0d0
      rejMove = .false.   
      nTargetMol = subIndxList(nTarget)

      atmType1 = atomArray(nType, 1)
      atmType2 = atomArray(nTargType, 1)
      rmin_ij = r_min_tab(atmType2, atmType1)
      rmin_ij = sqrt(rmin_ij)

      do iRosen = 1, nRosenTrials(nType)
        rosenTrial(iRosen)%x(1:nAtoms(nType)) = gasConfig(nType)%x(1:nAtoms(nType))
        rosenTrial(iRosen)%y(1:nAtoms(nType)) = gasConfig(nType)%y(1:nAtoms(nType))
        rosenTrial(iRosen)%z(1:nAtoms(nType)) = gasConfig(nType)%z(1:nAtoms(nType))
        x1 = rosenTrial(iRosen)%x(1)
        y1 = rosenTrial(iRosen)%y(1)
        z1 = rosenTrial(iRosen)%z(1)
        do i = 1, nAtoms(nType)
          rosenTrial(iRosen)%x(i) = rosenTrial(iRosen)%x(i) - x1
          rosenTrial(iRosen)%y(i) = rosenTrial(iRosen)%y(i) - y1
          rosenTrial(iRosen)%z(i) = rosenTrial(iRosen)%z(i) - z1
        enddo
        if(nAtoms(nType) .ne. 1) then
          x_rid_cm = gasConfig(nType)%x(1)
          y_rid_cm = gasConfig(nType)%y(1)
          z_rid_cm = gasConfig(nType)%z(1)    
        
!          Rotate xz axis
          rotang = two_pi*grnd()
          c_term = cos(rotang)
          s_term = sin(rotang)
          do i = 1,nAtoms(nType)
            x_shift = rosenTrial(iRosen)%x(i) - x_rid_cm
            z_shift = rosenTrial(iRosen)%z(i) - z_rid_cm        
            rosenTrial(iRosen)%x(i) = c_term*x_shift - s_term*z_shift + x_rid_cm
            rosenTrial(iRosen)%z(i) = s_term*x_shift + c_term*z_shift + z_rid_cm
          enddo   
        
!          Rotate yz axis
          rotang = two_pi*grnd()
          c_term = cos(rotang)
          s_term = sin(rotang)
          do i = 1,nAtoms(nType)
            y_shift = rosenTrial(iRosen)%y(i) - y_rid_cm
            z_shift = rosenTrial(iRosen)%z(i) - z_rid_cm        
            rosenTrial(iRosen)%y(i) = c_term*y_shift - s_term*z_shift + y_rid_cm
            rosenTrial(iRosen)%z(i) = s_term*y_shift + c_term*z_shift + z_rid_cm
          enddo   
  
!          Rotate xy axis
          rotang = two_pi*grnd()
          c_term = cos(rotang)
          s_term = sin(rotang)
          do i = 1,nAtoms(nType)
            x_shift = rosenTrial(iRosen)%x(i) - x_rid_cm
            y_shift = rosenTrial(iRosen)%y(i) - y_rid_cm        
            rosenTrial(iRosen)%x(i) = c_term*x_shift - s_term*y_shift + x_rid_cm
            rosenTrial(iRosen)%y(i) = s_term*x_shift + c_term*y_shift + y_rid_cm
          enddo    
        endif

        r = Dist_Critr * grnd()**(1d0/3d0)
!        if(r .lt. rmin_ij) then
!          overlap(iRosen) = .true.
!        endif

        call Generate_UnitSphere(dx, dy, dz)
        dx = r * dx
        dy = r * dy
        dz = r * dz      

        
        
        x1 = molArray(nTargType)%mol(nTargetMol)%x(1) + dx - rosenTrial(iRosen)%x(1)
        y1 = molArray(nTargType)%mol(nTargetMol)%y(1) + dy - rosenTrial(iRosen)%y(1)
        z1 = molArray(nTargType)%mol(nTargetMol)%z(1) + dz - rosenTrial(iRosen)%z(1)

        do i=1,nAtoms(nType)
          rosenTrial(iRosen)%x(i) = rosenTrial(iRosen)%x(i) + x1
          rosenTrial(iRosen)%y(i) = rosenTrial(iRosen)%y(i) + y1
          rosenTrial(iRosen)%z(i) = rosenTrial(iRosen)%z(i) + z1
        enddo 

!        call Rosen_BoltzWeight_Molecule_New(iRosen, nType, isIncluded, E_Trial(iRosen), overlap(iRosen))
        call Rosen_Mol_New(iRosen, nType, isIncluded, E_Trial(iRosen), overlap(iRosen))
      enddo

!      To prevent overflow errors the energies are scaled with respect to the most negative energy
!      or in other words the energy with the largest boltzmann weight. This rescaling naturally cancels
!      out when the probability is normalized. 
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen) - E_Max))
      enddo
!      if(all(ProbRosen .eq. 0d0)) then
!        rejMove = .true.
!        return
!      endif
      rosenNorm = sum(ProbRosen)
      
      ranNum = grnd()*rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
         nSel = nSel+1
         sumInt = sumInt + ProbRosen(nSel)
      enddo
      
      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif

!      Update the coordinates      
      do i = 1, nAtoms(nType)
        newMol%x(i) = rosenTrial(nSel)%x(i)
        newMol%y(i) = rosenTrial(nSel)%y(i)
        newMol%z(i) = rosenTrial(nSel)%z(i)
      enddo  
      rx = newMol%x(1) - molArray(nTargType)%mol(nTargetMol)%x(1)
      ry = newMol%y(1) - molArray(nTargType)%mol(nTargetMol)%y(1)
      rz = newMol%z(1) - molArray(nTargType)%mol(nTargetMol)%z(1)
      r = rx*rx + ry*ry + rz*rz
      if(r  .gt. Dist_Critr_sq) then
        write(*,*) "SCREWED UP!", rx, ry, rz, r, Dist_Critr, Dist_Critr_sq
      endif
      rosenRatio = (ProbRosen(nSel)*dble(nRosenTrials(nType)))/rosenNorm
      end subroutine
!=======================================================================
!     This subroutine is used for
      subroutine Ridgid_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions
      use EnergyPressurePointers, only: Rosen_Mol_Old
      use CBMC_Variables
      use CBMC_Utility
      use CoordinateFunctions
      implicit none

      integer, intent(in) :: nType, nMol, nTarget, nTargType
      real(dp), intent(out) :: rosenRatio

      logical :: isIncluded(1:maxMol)
      integer :: i, iRosen, nIndx, nTargetMol
      real(dp) :: E_Trial(1:maxRosenTrial)      
      real(dp) :: grnd,rotang
      real(dp) :: c_term,s_term
      real(dp) :: dx, dy, dz
      real(dp) :: x_shift,y_shift,z_shift
      real(dp) :: x_rid_cm, y_rid_cm,z_rid_cm
      real(dp) :: E_Max, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: r, x1, y1, z1       

      ProbRosen = 0d0      
      E_Trial = 0d0      
      nIndx = molArray(nType)%mol(nMol)%indx
!      write(2,*) NPART
      call Rosen_CreateSubset_Reverse(nTarget, nIndx, isIncluded)
!      call Rosen_BoltzWeight_Molecule_Old(molArray(nType)%mol(nMol)%x(:), molArray(nType)%mol(nMol)%y(:), &
!                                 molArray(nType)%mol(nMol)%z(:), nType, isIncluded, E_Trial(1))
      call Rosen_Mol_Old(molArray(nType)%mol(nMol)%x(:), molArray(nType)%mol(nMol)%y(:), &
                                 molArray(nType)%mol(nMol)%z(:), nType, isIncluded, E_Trial(1))
      nTargetMol = subIndxList(nTarget)
      do iRosen = 2, nRosenTrials(nType)
        newMol%x(1:nAtoms(nType)) = gasConfig(nType)%x(1:nAtoms(nType))
        newMol%y(1:nAtoms(nType)) = gasConfig(nType)%y(1:nAtoms(nType))
        newMol%z(1:nAtoms(nType)) = gasConfig(nType)%z(1:nAtoms(nType))
      
        if(nAtoms(nType) .ne. 1) then
          x_rid_cm = gasConfig(nType)%x(1)
          y_rid_cm = gasConfig(nType)%y(1)
          z_rid_cm = gasConfig(nType)%z(1)    
        
!          Rotate xz axis
          rotang = two_pi * grnd()
          c_term = cos(rotang)
          s_term = sin(rotang)
          do i = 1,nAtoms(nType)
            x_shift = newMol%x(i) - x_rid_cm
            z_shift = newMol%z(i) - z_rid_cm        
            newMol%x(i) = c_term*x_shift - s_term*z_shift + x_rid_cm
            newMol%z(i) = s_term*x_shift + c_term*z_shift + z_rid_cm
          enddo   
        
!          Rotate yz axis
          rotang = two_pi*grnd()
          c_term = cos(rotang)
          s_term = sin(rotang)
          do i = 1,nAtoms(nType)
            y_shift = newMol%y(i) - y_rid_cm
            z_shift = newMol%z(i) - z_rid_cm        
            newMol%y(i) = c_term*y_shift - s_term*z_shift + y_rid_cm
            newMol%z(i) = s_term*y_shift + c_term*z_shift + z_rid_cm
          enddo   
  
!          Rotate xy axis
          rotang = two_pi*grnd()
          c_term = cos(rotang)
          s_term = sin(rotang)
          do i = 1,nAtoms(nType)
            x_shift = newMol%x(i) - x_rid_cm
            y_shift = newMol%y(i) - y_rid_cm        
            newMol%x(i) = c_term*x_shift - s_term*y_shift + x_rid_cm
            newMol%y(i) = s_term*x_shift + c_term*y_shift + y_rid_cm
          enddo    
        endif

        r = Dist_Critr * grnd()**(1d0/3d0)
        call Generate_UnitSphere(dx, dy, dz)
        dx = r * dx
        dy = r * dy
        dz = r * dz      
        
        x1 = molArray(nTargType)%mol(nTargetMol)%x(1) + dx - newMol%x(1)
        y1 = molArray(nTargType)%mol(nTargetMol)%y(1) + dy - newMol%y(1)
        z1 = molArray(nTargType)%mol(nTargetMol)%z(1) + dz - newMol%z(1)

        do i=1,nAtoms(nType)
          newMol%x(i) = newMol%x(i) + x1
          newMol%y(i) = newMol%y(i) + y1
          newMol%z(i) = newMol%z(i) + z1
        enddo 

!        call Rosen_BoltzWeight_Molecule_Old( newMol%x(:), newMol%y(:), newMol%z(:), &
!                                   nType, isIncluded, E_Trial(iRosen) )
        call Rosen_Mol_Old( newMol%x(:), newMol%y(:), newMol%z(:), &
                                   nType, isIncluded, E_Trial(iRosen) )
      enddo

      E_Max = minval(E_Trial)
!      write(2,*) NPART
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen) - E_Max))
!        write(2,*) iRosen, E_Trial(iRosen), ProbRosen(iRosen)
      enddo
!      write(2,*)
!      flush(2)
      rosenNorm = sum(ProbRosen)
      rosenRatio = (ProbRosen(1)*dble(nRosenTrials(nType)))/rosenNorm
      
      end subroutine
!=======================================================================
      subroutine Simple_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions
      use EnergyPressurePointers, only: Rosen_Mol_New
      use CBMC_Variables
      use CBMC_Utility
      use CoordinateFunctions
      implicit none

      integer, intent(in) :: nType, nTarget, nTargType, nIndx
      logical, intent(out) :: rejMove      
      real(dp), intent(out):: rosenRatio
      logical, intent(in) :: isIncluded(:)
      
!      logical :: isIncluded(1:maxMol)
      logical :: overlap(1:maxRosenTrial)
      integer :: atmType1, atmType2
      integer :: i, iRosen, nTargetMol
      integer :: Atm1, Atm2, Atm3, Atm4, Atm5
      integer :: bondType, nSel  
      integer :: bendType, bendType2, bendType3, bendType4, bendType5, bendType6
      real(dp) :: E_Trial(1:maxRosenTrial)      
      real(dp) :: grnd
      real(dp) :: E_Max, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: ranNum, sumInt
      real(dp) :: k_bond, r_eq, r, Prob, r_sigma, rmax_sq
      real(dp) :: r1, r2, r3, r4
      real(dp) :: ang
      real(dp) :: ang1, ang2, ang3, dihed, dihed2
      real(dp) :: x1, y1, z1, dx, dy, dz 
      real(dp) :: rmin_ij
      type(SimpleAtomCoords) :: v1, v2, v3, v4, v5
      real(dp) :: wBending(1:maxRosenTrial)
      
      newMol%molType = nType      
!      call Rosen_CreateSubset(nTarget, isIncluded)
      E_Trial = 0d0
      rejMove = .false.      
      nTargetMol = subIndxList(nTarget)

      atmType1 = atomArray(nType, 1)
      atmType2 = atomArray(nTargType, 1)
      rmin_ij = r_min_tab(atmType2, atmType1)
      rmin_ij = sqrt(rmin_ij)

      do iRosen = 1, nRosenTrials(nType)
      
!        Initialize the first atom coordinates to 0      
        rosenTrial(iRosen)%x = 0d0
        rosenTrial(iRosen)%y = 0d0
        rosenTrial(iRosen)%z = 0d0

!        Depending the molecule type, choose an appropriate algorithm
!        and regrow the atoms starting from atom 1. 
        select case(nAtoms(nType))
        case(1)
          wBending(iRosen) = 1d0
        case(2)
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)        
          bondType = bondArray(nType,1)%bondType
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          rosenTrial(iRosen)%x(Atm2) = r * dx
          rosenTrial(iRosen)%y(Atm2) = r * dy
          rosenTrial(iRosen)%z(Atm2) = r * dz
          wBending(iRosen) = 1d0
        case(3)
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)
          Atm3 = pathArray(nType)%path(1, 3)
          call FindBond(nType, Atm1, Atm2, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          v1%x = -r*dx
          v1%y = -r*dy
          v1%z = -r*dz
          rosenTrial(iRosen)%x(atm2) = r * dx
          rosenTrial(iRosen)%y(atm2) = r * dy
          rosenTrial(iRosen)%z(atm2) = r * dz
          call FindBond(nType, Atm2, Atm3, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          bendType = bendArray(nType,1)%bendType
!          k_bend = bendData(bendType)%k_eq
!          ang_eq = bendData(bendType)%ang_eq
!          call GenerateBendAngle(ang, k_bend, ang_eq, Prob)
          call GenerateBendAngle(ang, bendType, Prob)
          call Generate_UnitCone(v1, r, ang, v2)
          rosenTrial(iRosen)%x(atm3) = rosenTrial(iRosen)%x(atm2) + v2%x 
          rosenTrial(iRosen)%y(atm3) = rosenTrial(iRosen)%y(atm2) + v2%y
          rosenTrial(iRosen)%z(atm3) = rosenTrial(iRosen)%z(atm2) + v2%z
          wBending(iRosen) = 1d0
        case(4)
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)
          Atm3 = pathArray(nType)%path(2, 1)
          Atm4 = pathArray(nType)%path(3, 1)
          call FindBond(nType, Atm1, Atm2, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r1, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          v1%x = -r1*dx
          v1%y = -r1*dy
          v1%z = -r1*dz
          rosenTrial(iRosen)%x(atm2) = r1 * dx
          rosenTrial(iRosen)%y(atm2) = r1 * dy
          rosenTrial(iRosen)%z(atm2) = r1 * dz
          call FindBond(nType, Atm2, Atm3, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r2, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindBond(nType, Atm2, Atm4, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r3, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
          call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
          call FindAngle(nType, Atm3, Atm2, Atm4, bendType3)
          wBending(iRosen) = 0d0
          call GenerateTwoBranches(ang1, ang2, dihed, bendType, bendType2, bendType3, wBending(iRosen))
          call Generate_UnitPyramid(v1, r2, r3, ang1, ang2, dihed, v2, v3)
          rosenTrial(iRosen)%x(atm3) = rosenTrial(iRosen)%x(atm2) + v2%x 
          rosenTrial(iRosen)%y(atm3) = rosenTrial(iRosen)%y(atm2) + v2%y
          rosenTrial(iRosen)%z(atm3) = rosenTrial(iRosen)%z(atm2) + v2%z
          rosenTrial(iRosen)%x(atm4) = rosenTrial(iRosen)%x(atm2) + v3%x 
          rosenTrial(iRosen)%y(atm4) = rosenTrial(iRosen)%y(atm2) + v3%y
          rosenTrial(iRosen)%z(atm4) = rosenTrial(iRosen)%z(atm2) + v3%z
        case(5)
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)
          Atm3 = pathArray(nType)%path(2, 1)
          Atm4 = pathArray(nType)%path(3, 1)
          Atm5 = pathArray(nType)%path(4, 1)
          call FindBond(nType, Atm1, Atm2, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r1, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          v1%x = -r1*dx
          v1%y = -r1*dy
          v1%z = -r1*dz
          rosenTrial(iRosen)%x(atm2) = r1 * dx
          rosenTrial(iRosen)%y(atm2) = r1 * dy
          rosenTrial(iRosen)%z(atm2) = r1 * dz
          call FindBond(nType, Atm2, Atm3, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r2, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindBond(nType, Atm2, Atm4, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r3, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindBond(nType, Atm2, Atm5, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r4, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
          call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
          call FindAngle(nType, Atm1, Atm2, Atm5, bendType3)
          call FindAngle(nType, Atm3, Atm2, Atm4, bendType4)
          call FindAngle(nType, Atm4, Atm2, Atm5, bendType5)
          call FindAngle(nType, Atm5, Atm2, Atm3, bendType6)
          wBending(iRosen) = 0d0
          call GenerateThreeBranches(ang1, ang2, ang3, dihed, dihed2, &
               bendType, bendType2, bendType3, bendType4, bendType5, bendType6, wBending(iRosen))
          call Generate_UnitTetrahedral(v1, r2, r3, r4, ang1, ang2, ang3, dihed, dihed2, v2, v3, v4)
          rosenTrial(iRosen)%x(atm3) = rosenTrial(iRosen)%x(atm2) + v2%x 
          rosenTrial(iRosen)%y(atm3) = rosenTrial(iRosen)%y(atm2) + v2%y
          rosenTrial(iRosen)%z(atm3) = rosenTrial(iRosen)%z(atm2) + v2%z
          rosenTrial(iRosen)%x(atm4) = rosenTrial(iRosen)%x(atm2) + v3%x 
          rosenTrial(iRosen)%y(atm4) = rosenTrial(iRosen)%y(atm2) + v3%y
          rosenTrial(iRosen)%z(atm4) = rosenTrial(iRosen)%z(atm2) + v3%z
          rosenTrial(iRosen)%x(atm5) = rosenTrial(iRosen)%x(atm2) + v4%x 
          rosenTrial(iRosen)%y(atm5) = rosenTrial(iRosen)%y(atm2) + v4%y
          rosenTrial(iRosen)%z(atm5) = rosenTrial(iRosen)%z(atm2) + v4%z
        case default
          write(*,*) "Error! Molecule has too many atoms for a simple regrowth"
          stop
        end select

        r = Dist_Critr * grnd()**(1d0/3d0)
        if(r .lt. rmin_ij) then
          overlap(iRosen) = .true.
        endif        
        call Generate_UnitSphere(dx, dy, dz)
        dx = r * dx
        dy = r * dy
        dz = r * dz   
        
        x1 = molArray(nTargType)%mol(nTargetMol)%x(1) + dx - rosenTrial(iRosen)%x(1)
        y1 = molArray(nTargType)%mol(nTargetMol)%y(1) + dy - rosenTrial(iRosen)%y(1)
        z1 = molArray(nTargType)%mol(nTargetMol)%z(1) + dz - rosenTrial(iRosen)%z(1)
        
        do i=1,nAtoms(nType)
          rosenTrial(iRosen)%x(i) = rosenTrial(iRosen)%x(i) + x1 
          rosenTrial(iRosen)%y(i) = rosenTrial(iRosen)%y(i) + y1 
          rosenTrial(iRosen)%z(i) = rosenTrial(iRosen)%z(i) + z1 
        enddo 

!        call Rosen_BoltzWeight_Molecule_New(iRosen, nType, isIncluded, E_Trial(iRosen), overlap(iRosen))
        call Rosen_Mol_New(iRosen, nType, isIncluded, E_Trial(iRosen), overlap(iRosen))
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = wBending(iRosen) * exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo

!      if(all(ProbRosen .le. 0d0)) then
!        rejMove = .true.
!        return
!      endif
      rosenNorm = sum(ProbRosen)
      
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
         nSel = nSel+1
         sumInt = sumInt + ProbRosen(nSel)
      enddo

      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif

!      Update the coordinates      
      newMol%x(1:nAtoms(nType)) = rosenTrial(nSel)%x(1:nAtoms(nType))
      newMol%y(1:nAtoms(nType)) = rosenTrial(nSel)%y(1:nAtoms(nType))
      newMol%z(1:nAtoms(nType)) = rosenTrial(nSel)%z(1:nAtoms(nType))

      wBending(1) = 1d0 
      if (nAtoms(nType) .eq. 4) then
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)
          Atm3 = pathArray(nType)%path(2, 1)
          Atm4 = pathArray(nType)%path(3, 1)
          call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
          call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
          call FindAngle(nType, Atm3, Atm2, Atm4, bendType3)
          v1%x = gasConfig(nType)%x(Atm1) - gasConfig(nType)%x(Atm2)
          v1%y = gasConfig(nType)%y(Atm1) - gasConfig(nType)%y(Atm2)
          v1%z = gasConfig(nType)%z(Atm1) - gasConfig(nType)%z(Atm2)
          v3%x = gasConfig(nType)%x(Atm3) - gasConfig(nType)%x(Atm2)
          v3%y = gasConfig(nType)%y(Atm3) - gasConfig(nType)%y(Atm2)
          v3%z = gasConfig(nType)%z(Atm3) - gasConfig(nType)%z(Atm2)
          v4%x = gasConfig(nType)%x(Atm4) - gasConfig(nType)%x(Atm2)
          v4%y = gasConfig(nType)%y(Atm4) - gasConfig(nType)%y(Atm2)
          v4%z = gasConfig(nType)%z(Atm4) - gasConfig(nType)%z(Atm2)
          wBending(1) = 0d0
          call GenerateTwoBranches_Reverse(v1, v3, v4, bendType, bendType2, bendType3, wBending(1))
      elseif (nAtoms(nType) .eq. 5) then
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)
          Atm3 = pathArray(nType)%path(2, 1)
          Atm4 = pathArray(nType)%path(3, 1)
          Atm5 = pathArray(nType)%path(4, 1)
          call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
          call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
          call FindAngle(nType, Atm1, Atm2, Atm5, bendType3)
          call FindAngle(nType, Atm3, Atm2, Atm4, bendType4)
          call FindAngle(nType, Atm4, Atm2, Atm5, bendType5)
          call FindAngle(nType, Atm5, Atm2, Atm3, bendType6)
          v1%x = gasConfig(nType)%x(Atm1) - gasConfig(nType)%x(Atm2)
          v1%y = gasConfig(nType)%y(Atm1) - gasConfig(nType)%y(Atm2)
          v1%z = gasConfig(nType)%z(Atm1) - gasConfig(nType)%z(Atm2)
          v3%x = gasConfig(nType)%x(Atm3) - gasConfig(nType)%x(Atm2)
          v3%y = gasConfig(nType)%y(Atm3) - gasConfig(nType)%y(Atm2)
          v3%z = gasConfig(nType)%z(Atm3) - gasConfig(nType)%z(Atm2)
          v4%x = gasConfig(nType)%x(Atm4) - gasConfig(nType)%x(Atm2)
          v4%y = gasConfig(nType)%y(Atm4) - gasConfig(nType)%y(Atm2)
          v4%z = gasConfig(nType)%z(Atm4) - gasConfig(nType)%z(Atm2)
          v5%x = gasConfig(nType)%x(Atm5) - gasConfig(nType)%x(Atm2)
          v5%y = gasConfig(nType)%y(Atm5) - gasConfig(nType)%y(Atm2)
          v5%z = gasConfig(nType)%z(Atm5) - gasConfig(nType)%z(Atm2)
          wBending(1) = 0d0
          call GenerateThreeBranches_Reverse(v1,v3,v4,v5,bendType,bendType2,bendType3,bendType4,bendType5,bendType6,wBending(1))
      endif
      do iRosen = 2, nRosenTrials(nType)
     
!        Initialize the first atom coordinates to 0      
        rosenTrial(iRosen)%x = 0d0
        rosenTrial(iRosen)%y = 0d0
        rosenTrial(iRosen)%z = 0d0

!        Depending the molecule type, choose an appropriate algorithm
!        and regrow the atoms starting from atom 1. 
        select case(nAtoms(nType))
        case(1)
          wBending(iRosen) = 1d0
        case(2)
          wBending(iRosen) = 1d0
        case(3)
          wBending(iRosen) = 1d0
        case(4)
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)
          Atm3 = pathArray(nType)%path(2, 1)
          Atm4 = pathArray(nType)%path(3, 1)
          call FindBond(nType, Atm1, Atm2, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r1, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          v1%x = -r1*dx
          v1%y = -r1*dy
          v1%z = -r1*dz
          rosenTrial(iRosen)%x(atm2) = r1 * dx
          rosenTrial(iRosen)%y(atm2) = r1 * dy
          rosenTrial(iRosen)%z(atm2) = r1 * dz
          call FindBond(nType, Atm2, Atm3, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r2, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindBond(nType, Atm2, Atm4, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r3, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
          call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
          call FindAngle(nType, Atm3, Atm2, Atm4, bendType3)
          wBending(iRosen) = 0d0
          call GenerateTwoBranches(ang1, ang2, dihed, bendType, bendType2, bendType3, wBending(iRosen))
          call Generate_UnitPyramid(v1, r2, r3, ang1, ang2, dihed, v2, v3)
          rosenTrial(iRosen)%x(atm3) = rosenTrial(iRosen)%x(atm2) + v2%x 
          rosenTrial(iRosen)%y(atm3) = rosenTrial(iRosen)%y(atm2) + v2%y
          rosenTrial(iRosen)%z(atm3) = rosenTrial(iRosen)%z(atm2) + v2%z
          rosenTrial(iRosen)%x(atm4) = rosenTrial(iRosen)%x(atm2) + v3%x 
          rosenTrial(iRosen)%y(atm4) = rosenTrial(iRosen)%y(atm2) + v3%y
          rosenTrial(iRosen)%z(atm4) = rosenTrial(iRosen)%z(atm2) + v3%z
        case(5)
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)
          Atm3 = pathArray(nType)%path(2, 1)
          Atm4 = pathArray(nType)%path(3, 1)
          Atm5 = pathArray(nType)%path(4, 1)
          call FindBond(nType, Atm1, Atm2, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r1, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          v1%x = -r1*dx
          v1%y = -r1*dy
          v1%z = -r1*dz
          rosenTrial(iRosen)%x(atm2) = r1 * dx
          rosenTrial(iRosen)%y(atm2) = r1 * dy
          rosenTrial(iRosen)%z(atm2) = r1 * dz
          call FindBond(nType, Atm2, Atm3, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r2, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindBond(nType, Atm2, Atm4, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r3, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindBond(nType, Atm2, Atm5, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r4, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
          call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
          call FindAngle(nType, Atm1, Atm2, Atm5, bendType3)
          call FindAngle(nType, Atm3, Atm2, Atm4, bendType4)
          call FindAngle(nType, Atm4, Atm2, Atm5, bendType5)
          call FindAngle(nType, Atm5, Atm2, Atm3, bendType6)
          wBending(iRosen) = 0d0
          call GenerateThreeBranches(ang1, ang2, ang3, dihed, dihed2, &
               bendType, bendType2, bendType3, bendType4, bendType5, bendType6, wBending(iRosen))
          call Generate_UnitTetrahedral(v1, r2, r3, r4, ang1, ang2, ang3, dihed, dihed2, v2, v3, v4)
          rosenTrial(iRosen)%x(atm3) = rosenTrial(iRosen)%x(atm2) + v2%x 
          rosenTrial(iRosen)%y(atm3) = rosenTrial(iRosen)%y(atm2) + v2%y
          rosenTrial(iRosen)%z(atm3) = rosenTrial(iRosen)%z(atm2) + v2%z
          rosenTrial(iRosen)%x(atm4) = rosenTrial(iRosen)%x(atm2) + v3%x 
          rosenTrial(iRosen)%y(atm4) = rosenTrial(iRosen)%y(atm2) + v3%y
          rosenTrial(iRosen)%z(atm4) = rosenTrial(iRosen)%z(atm2) + v3%z
          rosenTrial(iRosen)%x(atm5) = rosenTrial(iRosen)%x(atm2) + v4%x 
          rosenTrial(iRosen)%y(atm5) = rosenTrial(iRosen)%y(atm2) + v4%y
          rosenTrial(iRosen)%z(atm5) = rosenTrial(iRosen)%z(atm2) + v4%z
        case default
         stop "Error! Molecule has too many atoms for a simple regrowth"
        end select
      enddo
      
      rosenRatio = (ProbRosen(nSel)*sum(wBending))/rosenNorm
      
      end subroutine
!=======================================================================
      subroutine Simple_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions
      use EnergyPressurePointers, only: Rosen_Mol_Old
      use CBMC_Variables
      use CBMC_Utility
      use CoordinateFunctions
      implicit none

      integer, intent(in) :: nType, nTarget, nTargType, nMol
      real(dp), intent(out):: rosenRatio
      
      logical :: isIncluded(1:maxMol)

      integer :: i, iRosen, nIndx, nTargetMol
      integer :: Atm1, Atm2, Atm3, Atm4, Atm5
      integer :: bondType
      integer :: bendType, bendType2, bendType3, bendType4, bendType5, bendType6
      real(dp) :: E_Trial(1:maxRosenTrial)      
      real(dp) :: grnd
      real(dp) :: E_Max, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: k_bond, r_eq, r, Prob, r_sigma, rmax_sq
      real(dp) :: r1, r2, r3, r4
      real(dp) :: ang
      real(dp) :: ang1, ang2, ang3, dihed, dihed2
      real(dp) :: x1, y1, z1, dx, dy, dz 
      type(SimpleAtomCoords) :: v1, v2, v3, v4
      real(dp) :: wBending(1:maxRosenTrial)
      
      ProbRosen = 0d0      
      E_Trial = 0d0      
      nTargetMol = subIndxList(nTarget)
      nIndx = molArray(nType)%mol(nMol)%indx
      call Rosen_CreateSubset_Reverse(nTarget, nIndx, isIncluded)
!      call Rosen_BoltzWeight_Molecule_Old(molArray(nType)%mol(nMol)%x(:), molArray(nType)%mol(nMol)%y(:), &
!                                 molArray(nType)%mol(nMol)%z(:), nType, isIncluded, E_Trial(1))      

      call Rosen_Mol_Old(molArray(nType)%mol(nMol)%x(:), molArray(nType)%mol(nMol)%y(:), &
                                 molArray(nType)%mol(nMol)%z(:), nType, isIncluded, E_Trial(1))     

      do iRosen = 2, nRosenTrials(nType)
      
        newMol%x = 0d0
        newMol%y = 0d0
        newMol%z = 0d0

!        Depending the molecule type, choose an appropriate algorithm
!        and regrow the atoms starting from atom 1. 
        select case(nAtoms(nType))
        case(1)
          wBending(iRosen) = 1d0
        case(2)
          bondType = bondArray(nType,1)%bondType
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          newMol%x(2) = r*dx
          newMol%y(2) = r*dy
          newMol%z(2) = r*dz
          wBending(iRosen) = 1d0
        case(3)
          Atm1 = pathArray(nType)%path(1,1)
          Atm2 = pathArray(nType)%path(1,2)
          Atm3 = pathArray(nType)%path(1,3)
          call FindBond(nType,Atm1, Atm2, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          v1%x = -r*dx
          v1%y = -r*dy
          v1%z = -r*dz
          newMol%x(atm2) = r*dx
          newMol%y(atm2) = r*dy
          newMol%z(atm2) = r*dz
          call FindBond(nType, Atm2, Atm3, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          bendType = bendArray(nType,1)%bendType
!          k_bend = bendData(bendType)%k_eq
!          ang_eq = bendData(bendType)%ang_eq
!          call GenerateBendAngle(ang, k_bend, ang_eq, Prob)
          call GenerateBendAngle(ang, bendType, Prob)
          call Generate_UnitCone(v1, r, ang, v2)
          newMol%x(atm3) = newMol%x(atm2) + v2%x 
          newMol%y(atm3) = newMol%y(atm2) + v2%y
          newMol%z(atm3) = newMol%z(atm2) + v2%z
          wBending(iRosen) = 1d0
        case(4)
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)
          Atm3 = pathArray(nType)%path(2, 1)
          Atm4 = pathArray(nType)%path(3, 1)
          call FindBond(nType, Atm1, Atm2, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r1, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          v1%x = -r1*dx
          v1%y = -r1*dy
          v1%z = -r1*dz
          newMol%x(atm2) = r1 * dx
          newMol%y(atm2) = r1 * dy
          newMol%z(atm2) = r1 * dz
          call FindBond(nType, Atm2, Atm3, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r2, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindBond(nType, Atm2, Atm4, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r3, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
          call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
          call FindAngle(nType, Atm3, Atm2, Atm4, bendType3)
          wBending(iRosen) = 0d0
          call GenerateTwoBranches(ang1, ang2, dihed, bendType, bendType2, bendType3, wBending(iRosen))
          call Generate_UnitPyramid(v1, r2, r3, ang1, ang2, dihed, v2, v3)
          newMol%x(atm3) = newMol%x(atm2) + v2%x 
          newMol%y(atm3) = newMol%y(atm2) + v2%y
          newMol%z(atm3) = newMol%z(atm2) + v2%z
          newMol%x(atm4) = newMol%x(atm2) + v3%x 
          newMol%y(atm4) = newMol%y(atm2) + v3%y
          newMol%z(atm4) = newMol%z(atm2) + v3%z
        case(5)
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)
          Atm3 = pathArray(nType)%path(2, 1)
          Atm4 = pathArray(nType)%path(3, 1)
          Atm5 = pathArray(nType)%path(4, 1)
          call FindBond(nType, Atm1, Atm2, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r1, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          v1%x = -r1*dx
          v1%y = -r1*dy
          v1%z = -r1*dz
          rosenTrial(iRosen)%x(atm2) = r1 * dx
          rosenTrial(iRosen)%y(atm2) = r1 * dy
          rosenTrial(iRosen)%z(atm2) = r1 * dz
          call FindBond(nType, Atm2, Atm3, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r2, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindBond(nType, Atm2, Atm4, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r3, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindBond(nType, Atm2, Atm5, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r4, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
          call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
          call FindAngle(nType, Atm1, Atm2, Atm5, bendType3)
          call FindAngle(nType, Atm3, Atm2, Atm4, bendType4)
          call FindAngle(nType, Atm4, Atm2, Atm5, bendType5)
          call FindAngle(nType, Atm5, Atm2, Atm3, bendType6)
          wBending(iRosen) = 0d0
          call GenerateThreeBranches(ang1, ang2, ang3, dihed, dihed2, &
               bendType, bendType2, bendType3, bendType4, bendType5, bendType6, wBending(iRosen))
          call Generate_UnitTetrahedral(v1, r2, r3, r4, ang1, ang2, ang3, dihed, dihed2, v2, v3, v4)
          newMol%x(atm3) = newMol%x(atm2) + v2%x 
          newMol%y(atm3) = newMol%y(atm2) + v2%y
          newMol%z(atm3) = newMol%z(atm2) + v2%z
          newMol%x(atm4) = newMol%x(atm2) + v3%x 
          newMol%y(atm4) = newMol%y(atm2) + v3%y
          newMol%z(atm4) = newMol%z(atm2) + v3%z
          newMol%x(atm5) = newMol%x(atm2) + v4%x 
          newMol%y(atm5) = newMol%y(atm2) + v4%y
          newMol%z(atm5) = newMol%z(atm2) + v4%z
        case default
          write(*,*) "Error! Molecule has too many atoms for a simple regrowth"
          stop
        end select
        
        r = Dist_Critr + 1d0
        do while(r*r .gt. Dist_Critr_sq )
          r = Dist_Critr * grnd()**(1d0/3d0)
        enddo
        call Generate_UnitSphere(dx, dy, dz)
        dx = r * dx
        dy = r * dy
        dz = r * dz   
        
        x1 = molArray(nTargType)%mol(nTargetMol)%x(1) + dx - newMol%x(1)
        y1 = molArray(nTargType)%mol(nTargetMol)%y(1) + dy - newMol%y(1)
        z1 = molArray(nTargType)%mol(nTargetMol)%z(1) + dz - newMol%z(1)
        
        do i=1,nAtoms(nType)
          newMol%x(i) = newMol%x(i) + x1
          newMol%y(i) = newMol%y(i) + y1
          newMol%z(i) = newMol%z(i) + z1
        enddo 

!        call Rosen_BoltzWeight_Molecule_Old( newMol%x(:), newMol%y(:), newMol%z(:), nType, isIncluded, E_Trial(iRosen) )
        call Rosen_Mol_Old( newMol%x(:), newMol%y(:), newMol%z(:), nType, isIncluded, E_Trial(iRosen) )

      enddo

      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = wBending(iRosen) * exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      rosenNorm = sum(ProbRosen)

      do iRosen = 1, nRosenTrials(nType)
     
!        Initialize the first atom coordinates to 0      
        rosenTrial(iRosen)%x = 0d0
        rosenTrial(iRosen)%y = 0d0
        rosenTrial(iRosen)%z = 0d0

!        Depending the molecule type, choose an appropriate algorithm
!        and regrow the atoms starting from atom 1. 
        select case(nAtoms(nType))
        case(1)
          wBending(iRosen) = 1d0
        case(2)
          wBending(iRosen) = 1d0
        case(3)
          wBending(iRosen) = 1d0
        case(4)
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)
          Atm3 = pathArray(nType)%path(2, 1)
          Atm4 = pathArray(nType)%path(3, 1)
          call FindBond(nType, Atm1, Atm2, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r1, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          v1%x = -r1*dx
          v1%y = -r1*dy
          v1%z = -r1*dz
          rosenTrial(iRosen)%x(atm2) = r1 * dx
          rosenTrial(iRosen)%y(atm2) = r1 * dy
          rosenTrial(iRosen)%z(atm2) = r1 * dz
          call FindBond(nType, Atm2, Atm3, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r2, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindBond(nType, Atm2, Atm4, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r3, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
          call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
          call FindAngle(nType, Atm3, Atm2, Atm4, bendType3)
          wBending(iRosen) = 0d0
          call GenerateTwoBranches(ang1, ang2, dihed, bendType, bendType2, bendType3, wBending(iRosen))
          call Generate_UnitPyramid(v1, r2, r3, ang1, ang2, dihed, v2, v3)
          rosenTrial(iRosen)%x(atm3) = rosenTrial(iRosen)%x(atm2) + v2%x 
          rosenTrial(iRosen)%y(atm3) = rosenTrial(iRosen)%y(atm2) + v2%y
          rosenTrial(iRosen)%z(atm3) = rosenTrial(iRosen)%z(atm2) + v2%z
          rosenTrial(iRosen)%x(atm4) = rosenTrial(iRosen)%x(atm2) + v3%x 
          rosenTrial(iRosen)%y(atm4) = rosenTrial(iRosen)%y(atm2) + v3%y
          rosenTrial(iRosen)%z(atm4) = rosenTrial(iRosen)%z(atm2) + v3%z
        case(5)
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)
          Atm3 = pathArray(nType)%path(2, 1)
          Atm4 = pathArray(nType)%path(3, 1)
          Atm5 = pathArray(nType)%path(4, 1)
          call FindBond(nType, Atm1, Atm2, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r1, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          v1%x = -r1*dx
          v1%y = -r1*dy
          v1%z = -r1*dz
          rosenTrial(iRosen)%x(atm2) = r1 * dx
          rosenTrial(iRosen)%y(atm2) = r1 * dy
          rosenTrial(iRosen)%z(atm2) = r1 * dz
          call FindBond(nType, Atm2, Atm3, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r2, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindBond(nType, Atm2, Atm4, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r3, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindBond(nType, Atm2, Atm5, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          r_sigma = bondData(bondType)%r_sigma
          rmax_sq = bondData(bondType)%rmax_sq
          call GenerateGaussianBondLength(r4, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
          call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
          call FindAngle(nType, Atm1, Atm2, Atm5, bendType3)
          call FindAngle(nType, Atm3, Atm2, Atm4, bendType4)
          call FindAngle(nType, Atm4, Atm2, Atm5, bendType5)
          call FindAngle(nType, Atm5, Atm2, Atm3, bendType6)
          wBending(iRosen) = 0d0
          call GenerateThreeBranches(ang1, ang2, ang3, dihed, dihed2, &
               bendType, bendType2, bendType3, bendType4, bendType5, bendType6, wBending(iRosen))
          call Generate_UnitTetrahedral(v1, r2, r3, r4, ang1, ang2, ang3, dihed, dihed2, v2, v3, v4)
          rosenTrial(iRosen)%x(atm3) = rosenTrial(iRosen)%x(atm2) + v2%x 
          rosenTrial(iRosen)%y(atm3) = rosenTrial(iRosen)%y(atm2) + v2%y
          rosenTrial(iRosen)%z(atm3) = rosenTrial(iRosen)%z(atm2) + v2%z
          rosenTrial(iRosen)%x(atm4) = rosenTrial(iRosen)%x(atm2) + v3%x 
          rosenTrial(iRosen)%y(atm4) = rosenTrial(iRosen)%y(atm2) + v3%y
          rosenTrial(iRosen)%z(atm4) = rosenTrial(iRosen)%z(atm2) + v3%z
          rosenTrial(iRosen)%x(atm5) = rosenTrial(iRosen)%x(atm2) + v4%x 
          rosenTrial(iRosen)%y(atm5) = rosenTrial(iRosen)%y(atm2) + v4%y
          rosenTrial(iRosen)%z(atm5) = rosenTrial(iRosen)%z(atm2) + v4%z
        case default
         stop "Error! Molecule has too many atoms for a simple regrowth"
        end select
      enddo
      
      rosenRatio = (ProbRosen(1)*sum(wBending))/rosenNorm
      
      end subroutine
!=======================================================================
      subroutine StraightChain_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions_LJ_Q
      use EnergyPressurePointers, only: Rosen_Atom_New, Rosen_Atom_Intra_New
      use CBMC_Variables
      use CBMC_Utility
      use CoordinateFunctions
      implicit none

      integer, intent(in) :: nType, nTarget, nTargType
      real(dp), intent(out):: rosenRatio
      logical, intent(out) :: rejMove
      logical, intent(in) :: isIncluded(:)

!      logical :: isIncluded(1:maxMol)
      integer :: iRosen, nSel, nIndx, nTargetMol
      integer :: nRegrown
      integer :: Atm1, Atm2, Atm3, Atm4
      integer :: cnt
      integer :: bondType, bendType, torsType
      logical :: overlap(1:maxRosenTrial)
      logical :: regrown(1:maxAtoms)
      real(dp) :: E_Trial(1:maxRosenTrial), E_Complete, E_Intra
      real(dp) :: grnd
      real(dp) :: E_Max, ProbRosen(1:maxRosenTrial), rosenNorm, wReverse
      real(dp) :: ranNum, sumInt
      real(dp) :: k_bond, r_eq, r, Prob, r_sigma, rmax_sq
      real(dp) :: bend_angle, tors_angle
      real(dp) :: dx, dy, dz 
      type(SimpleAtomCoords) :: trialPos(1:maxRosenTrial)
      type(SimpleAtomCoords) :: v1, v2, v3
      
   
  
!      call Rosen_CreateSubset(nTarget, isIncluded)
      nTargetMol = subIndxList(nTarget)
      newMol%molType = nType 
      newMol%x = 0d0
      newMol%y = 0d0
      newMol%z = 0d0
      regrown = .false.
      rejMove = .false.
      nRegrown = 0

!      Begin the regrowth process by choosing an insertion site for the first atom in the chain
      
      E_Trial = 0d0
      E_Complete = 0d0
      overlap = .false.
      do iRosen = 1, nRosenTrials(nType)
        r = Dist_Critr * grnd()**(1d0/3d0)
        call Generate_UnitSphere(dx, dy, dz)
        trialPos(iRosen)%x = r * dx + molArray(nTargType)%mol(nTargetMol)%x(1)
        trialPos(iRosen)%y = r * dy + molArray(nTargType)%mol(nTargetMol)%y(1)
        trialPos(iRosen)%z = r * dz + molArray(nTargType)%mol(nTargetMol)%z(1)
        call Rosen_Atom_New(nType, 1, trialPos(iRosen), isIncluded,  E_Trial(iRosen), overlap(iRosen))
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      if(all(ProbRosen .le. 0d0)) then
        rejMove = .true.
        return
      endif
      rosenNorm = sum(ProbRosen)
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbRosen(nSel)
      enddo

      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif
      rosenRatio = ProbRosen(nSel)/rosenNorm
      regrown(1) = .true.
      nRegrown = nRegrown + 1
      newMol%x(1) = trialPos(nSel)%x 
      newMol%y(1) = trialPos(nSel)%y 
      newMol%z(1) = trialPos(nSel)%z

!      Having inserted the first atom, when we go to insert the second atom we must account for the fact that
!      the second atom must be inserted at the correct bond distance from the first atom.
      Atm1 = 1
      Atm2 = regrowOrder(nType, 2) 
      call FindBond(nType,Atm1, Atm2, bondType)
      k_bond = bondData(bondType)%k_eq
      r_eq = bondData(bondType)%r_eq
      r_sigma = bondData(bondType)%r_sigma
      rmax_sq = bondData(bondType)%rmax_sq
      overlap = .false.
      E_Trial = 0d0
      do iRosen = 1, nRosenTrials(nType)
        call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
        call Generate_UnitSphere(dx, dy, dz)
        trialPos(iRosen)%x = r*dx + newMol%x(Atm1) 
        trialPos(iRosen)%y = r*dy + newMol%y(Atm1)
        trialPos(iRosen)%z = r*dz + newMol%z(Atm1)
        call Rosen_Atom_New(nType, Atm2, trialPos(iRosen), isIncluded,  E_Trial(iRosen), overlap(iRosen))
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      if(all(ProbRosen .le. 0d0)) then
        rejMove = .true.
        return
      endif
      rosenNorm = sum(ProbRosen)
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbRosen(nSel)
      enddo
      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif
      rosenRatio = rosenRatio*ProbRosen(nSel)/rosenNorm
      regrown(Atm2) = .true.
      nRegrown = nRegrown + 1
      newMol%x(Atm2) = trialPos(nSel)%x 
      newMol%y(Atm2) = trialPos(nSel)%y 
      newMol%z(Atm2) = trialPos(nSel)%z

!       For the third atom we must begin to include the bending angle in choosing its trial positions. The first thing
!       we must consider is if the second regrown atom was a terminal atom or a linker atom. If it was a terminal atom 
!       then the bending angle must use the 1st atom as the central atom for the angle generation.  However if it was a link in the chain then the 2nd 
!       atom regrown should be used as the central atom.
      if(topolArray(nType)%atom(Atm2) .eq. 2) then
        Atm1 = regrowOrder(nType, 1) 
        Atm2 = regrowOrder(nType, 2) 
        Atm3 = regrowOrder(nType, 3) 
      else
        Atm1 = regrowOrder(nType, 2) 
        Atm2 = regrowOrder(nType, 1) 
        Atm3 = regrowOrder(nType, 3) 
      endif
      v1%x = newMol%x(Atm1) - newMol%x(Atm2)
      v1%y = newMol%y(Atm1) - newMol%y(Atm2)
      v1%z = newMol%z(Atm1) - newMol%z(Atm2)
      call FindBond(nType, Atm2, Atm3, bondType)
      k_bond = bondData(bondType)%k_eq
      r_eq = bondData(bondType)%r_eq
      r_sigma = bondData(bondType)%r_sigma
      rmax_sq = bondData(bondType)%rmax_sq
      call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
!      k_bend = bendData(bendType)%k_eq
!      ang_eq = bendData(bendType)%ang_eq
      overlap = .false.
      do iRosen = 1, nRosenTrials(nType)
        call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
!        call GenerateBendAngle(bend_angle, k_bend, ang_eq, Prob)
        call GenerateBendAngle(bend_angle, bendType, Prob)
        call Generate_UnitCone(v1, r, bend_angle, v2)
        trialPos(iRosen)%x = v2%x + newMol%x(Atm2) 
        trialPos(iRosen)%y = v2%y + newMol%y(Atm2)
        trialPos(iRosen)%z = v2%z + newMol%z(Atm2)
        call Rosen_Atom_New(nType, Atm3, trialPos(iRosen), isIncluded,  E_Trial(iRosen), overlap(iRosen))
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
!      if(all(ProbRosen .le. 0d0)) then
!        rejMove = .true.
!        return
!      endif
      rosenNorm = sum(ProbRosen)
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbRosen(nSel)
      enddo
      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif
      rosenRatio = rosenRatio*ProbRosen(nSel)/rosenNorm
      regrown(Atm3) = .true.
      nRegrown = nRegrown + 1
      newMol%x(Atm3) = trialPos(nSel)%x 

      newMol%y(Atm3) = trialPos(nSel)%y 
      newMol%z(Atm3) = trialPos(nSel)%z


!      Now that three atoms have been regrown, all remaining atoms must take the torsional angles into account.
      cnt = 3
      do while(any(regrown .eqv. .false.))
        cnt = cnt + 1
        atm4 = regrowOrder(nType, cnt) 
        call FindAtomsFromPath(nType, regrown, 1, Atm4, Atm1, Atm2, Atm3)
        call FindBond(nType, Atm3, Atm4, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        r_sigma = bondData(bondType)%r_sigma
        rmax_sq = bondData(bondType)%rmax_sq
        call FindAngle(nType, Atm2, Atm3, Atm4, bendType)
!        k_bend = bendData(bendType)%k_eq
!        ang_eq = bendData(bendType)%ang_eq
        call FindTorsion(nType, Atm1, Atm2, Atm3, Atm4, torsType)
        v1%x = newMol%x(Atm1) - newMol%x(Atm3)
        v1%y = newMol%y(Atm1) - newMol%y(Atm3)
        v1%z = newMol%z(Atm1) - newMol%z(Atm3)

        v2%x = newMol%x(Atm2) - newMol%x(Atm3)
        v2%y = newMol%y(Atm2) - newMol%y(Atm3)
        v2%z = newMol%z(Atm2) - newMol%z(Atm3)
        overlap = .false.
        E_Trial = 0d0
        do iRosen = 1, nRosenTrials(nType)
          call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
!          call GenerateBendAngle(bend_angle, k_bend, ang_eq, Prob)
          call GenerateBendAngle(bend_angle, bendType, Prob)
          call GenerateTorsAngle(tors_angle, torsType, Prob)
          call Generate_UnitTorsion(v1, v2, r, bend_angle, tors_angle, v3)
          trialPos(iRosen)%x = v3%x + newMol%x(Atm3) 
          trialPos(iRosen)%y = v3%y + newMol%y(Atm3)
          trialPos(iRosen)%z = v3%z + newMol%z(Atm3)
          call Rosen_Atom_New(nType, Atm4, trialPos(iRosen), isIncluded,  E_Trial(iRosen), overlap(iRosen))
          if (.not. overlap(iRosen)) then
            call Rosen_Atom_Intra_New(nType, Atm4, trialPos(iRosen), regrown, E_Intra, overlap(iRosen))
          endif
          E_Trial(iRosen) = E_Trial(iRosen) + E_Intra
        enddo
        E_Max = minval(E_Trial)
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
        enddo
!        if(all(ProbRosen .le. 0d0)) then
!          rejMove = .true.
!          return
!        endif
        rosenNorm = sum(ProbRosen)
        ranNum = grnd() * rosenNorm
        sumInt = ProbRosen(1)
        nSel = 1
        do while(sumInt .lt. ranNum)
          nSel = nSel + 1
          sumInt = sumInt + ProbRosen(nSel)
        enddo
        if(overlap(nSel) .eqv. .true.) then
          rejMove = .true.
          return
        endif
        rosenRatio = rosenRatio*ProbRosen(nSel)/rosenNorm
        regrown(Atm4) = .true.
        nRegrown = nRegrown + 1
        newMol%x(Atm4) = trialPos(nSel)%x 
        newMol%y(Atm4) = trialPos(nSel)%y 
        newMol%z(Atm4) = trialPos(nSel)%z
      enddo
      call StraightChain_RosenConfigGen_GasPhase_Reverse(nType, wReverse)
      rosenRatio = rosenRatio*wReverse


      end subroutine
!=======================================================================
      subroutine StraightChain_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions_LJ_Q
      use EnergyPressurePointers, only: Rosen_Atom_New, Rosen_Atom_Old, Rosen_Atom_Intra_Old
      use CBMC_Variables
      use CBMC_Utility
      use CoordinateFunctions
      implicit none

      integer, intent(in) :: nType, nMol, nTarget, nTargType
      real(dp), intent(out):: rosenRatio

      integer :: iRosen, nIndx, nTargetMol
      integer :: Atm1, Atm2, Atm3, Atm4
      integer :: bondType, bendType, torsType, cnt
      logical :: overlap, rejMove
      logical :: isIncluded(1:maxMol)
      logical :: regrown(1:maxAtoms)
      real(dp) :: E_Trial(1:maxRosenTrial), E_Intra
      real(dp) :: grnd
      real(dp) :: E_Max, ProbRosen(1:maxRosenTrial), rosenNorm, wForward
      real(dp) :: k_bond, r_eq, r, Prob, r_sigma, rmax_sq
      real(dp) :: bend_angle, tors_angle
      real(dp) :: dx, dy, dz 
      type(SimpleAtomCoords) :: trialPos
      type(SimpleAtomCoords) :: v1, v2, v3
      
      ProbRosen = 0d0      
      E_Trial = 0d0      
      nTargetMol = subIndxList(nTarget)
      nIndx = molArray(nType)%mol(nMol)%indx
      call Rosen_CreateSubset_Reverse(nTarget, nIndx, isIncluded)
      regrown = .false.

!      Begin the regrowth process by choosing an insertion site for the first atom in the chain
      E_Trial = 0d0
!      E_Complete = 0d0
      rosenRatio = 1d0
      call Rosen_Atom_Old(nType, nMol, 1, isIncluded,  E_Trial(1))
      do iRosen = 2, nRosenTrials(nType)
        r = Dist_Critr * grnd()**(1d0/3d0)
        call Generate_UnitSphere(dx, dy, dz)
        trialPos%x = r * dx + molArray(nTargType)%mol(nTargetMol)%x(1)
        trialPos%y = r * dy + molArray(nTargType)%mol(nTargetMol)%y(1)
        trialPos%z = r * dz + molArray(nTargType)%mol(nTargetMol)%z(1)
        call Rosen_Atom_New(nType, 1, trialPos, isIncluded,  E_Trial(iRosen), overlap)
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      rosenNorm = sum(ProbRosen)
      rosenRatio = rosenRatio*ProbRosen(1)/rosenNorm
      regrown(1) = .true.


!      Having inserted the first atom, when we go to insert the second atom we must account for the fact that
!      the second atom must be inserted at the correct bond distance from the first atom.
      Atm1 = 1
      Atm2 = regrowOrder(nType, 2) 
      call FindBond(nType,Atm1, Atm2, bondType)
      k_bond = bondData(bondType)%k_eq
      r_eq = bondData(bondType)%r_eq
      r_sigma = bondData(bondType)%r_sigma
      rmax_sq = bondData(bondType)%rmax_sq
      call Rosen_Atom_Old(nType, nMol, Atm2, isIncluded,  E_Trial(1))
      do iRosen = 2, nRosenTrials(nType)
        call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
        call Generate_UnitSphere(dx, dy, dz)
        trialPos%x = r*dx + molArray(nType)%mol(nMol)%x(Atm1)
        trialPos%y = r*dy + molArray(nType)%mol(nMol)%y(Atm1)
        trialPos%z = r*dz + molArray(nType)%mol(nMol)%z(Atm1)
        call Rosen_Atom_New(nType, Atm2, trialPos, isIncluded,  E_Trial(iRosen), overlap)
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      rosenNorm = sum(ProbRosen)
      rosenRatio = rosenRatio*ProbRosen(1)/rosenNorm
      regrown(Atm2) = .true.

!       For the third atom we must begin to include the bending angle in choosing its trial positions. The first thing
!       we must consider is if the second regrown atom was a terminal atom or a linker atom. If it was a terminal atom 
!       then the bending angle must use the 1st atom as the central atom for the angle generation.  However if it was a link in the chain then the 2nd 
!       atom regrown should be used as the central atom.
      if(topolArray(nType)%atom(Atm2) .eq. 2) then
        Atm1 = regrowOrder(nType, 1) 
        Atm2 = regrowOrder(nType, 2) 
        Atm3 = regrowOrder(nType, 3) 
      else
        Atm1 = regrowOrder(nType, 2) 
        Atm2 = regrowOrder(nType, 1) 
        Atm3 = regrowOrder(nType, 3) 
      endif

      v1%x = molArray(nType)%mol(nMol)%x(Atm1) - molArray(nType)%mol(nMol)%x(Atm2)
      v1%y = molArray(nType)%mol(nMol)%y(Atm1) - molArray(nType)%mol(nMol)%y(Atm2)
      v1%z = molArray(nType)%mol(nMol)%z(Atm1) - molArray(nType)%mol(nMol)%z(Atm2)
      call FindBond(nType, Atm2, Atm3, bondType)
      k_bond = bondData(bondType)%k_eq
      r_eq = bondData(bondType)%r_eq
      r_sigma = bondData(bondType)%r_sigma
      rmax_sq = bondData(bondType)%rmax_sq
      call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
!      k_bend = bendData(bendType)%k_eq
!      ang_eq = bendData(bendType)%ang_eq

      call Rosen_Atom_Old(nType, nMol, Atm3, isIncluded,  E_Trial(1))
      trialPos%x = molArray(nType)%mol(nMol)%x(Atm3) 
      trialPos%y = molArray(nType)%mol(nMol)%y(Atm3)
      trialPos%z = molArray(nType)%mol(nMol)%z(Atm3)
      call Rosen_Atom_Intra_Old(nType, nMol, Atm3, trialPos, regrown, E_Intra)
      E_Trial(1) = E_Trial(1) + E_Intra
      do iRosen = 2, nRosenTrials(nType)
        call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
!        call GenerateBendAngle(bend_angle, k_bend, ang_eq, Prob)
        call GenerateBendAngle(bend_angle, bendType, Prob)
        call Generate_UnitCone(v1, r, bend_angle, v2)
        trialPos%x = v2%x + molArray(nType)%mol(nMol)%x(Atm2)
        trialPos%y = v2%y + molArray(nType)%mol(nMol)%y(Atm2)
        trialPos%z = v2%z + molArray(nType)%mol(nMol)%z(Atm2)
        call Rosen_Atom_New(nType, Atm3, trialPos, isIncluded,  E_Trial(iRosen), overlap)
        call Rosen_Atom_Intra_Old(nType, nMol, Atm3, trialPos, regrown, E_Intra)
        E_Trial(iRosen) = E_Trial(iRosen) + E_Intra
      enddo

      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      rosenNorm = sum(ProbRosen)
      rosenRatio = rosenRatio*ProbRosen(1)/rosenNorm
      regrown(Atm3) = .true.

!      Now that three atoms have been regrown, all remaining atoms must take the torsional angles into account.
      cnt = 3
      do while(any(regrown .eqv. .false.))
        cnt = cnt + 1
        atm4 = regrowOrder(nType, cnt) 
        call FindAtomsFromPath(nType, regrown, 1, Atm4, Atm1, Atm2, Atm3)
        call FindBond(nType, Atm3, Atm4, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        r_sigma = bondData(bondType)%r_sigma
        rmax_sq = bondData(bondType)%rmax_sq
        call FindAngle(nType, Atm2, Atm3, Atm4, bendType)
!        k_bend = bendData(bendType)%k_eq
!        ang_eq = bendData(bendType)%ang_eq
        call FindTorsion(nType, Atm1, Atm2, Atm3, Atm4, torsType)
        v1%x = molArray(nType)%mol(nMol)%x(Atm1) - molArray(nType)%mol(nMol)%x(Atm3)
        v1%y = molArray(nType)%mol(nMol)%y(Atm1) - molArray(nType)%mol(nMol)%y(Atm3)
        v1%z = molArray(nType)%mol(nMol)%z(Atm1) - molArray(nType)%mol(nMol)%z(Atm3)

        v2%x = molArray(nType)%mol(nMol)%x(Atm2) - molArray(nType)%mol(nMol)%x(Atm3)
        v2%y = molArray(nType)%mol(nMol)%y(Atm2) - molArray(nType)%mol(nMol)%y(Atm3)
        v2%z = molArray(nType)%mol(nMol)%z(Atm2) - molArray(nType)%mol(nMol)%z(Atm3)
        overlap = .false.
        call Rosen_Atom_Old(nType, nMol, Atm4, isIncluded,  E_Trial(1))
        trialPos%x = molArray(nType)%mol(nMol)%x(Atm4) 
        trialPos%y = molArray(nType)%mol(nMol)%y(Atm4)
        trialPos%z = molArray(nType)%mol(nMol)%z(Atm4)
        call Rosen_Atom_Intra_Old(nType, nMol, Atm4, trialPos, regrown, E_Intra)
        E_Trial(1) = E_Trial(1) + E_Intra
        do iRosen = 2, nRosenTrials(nType)
          call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
!          call GenerateBendAngle(bend_angle, k_bend, ang_eq, Prob)
          call GenerateBendAngle(bend_angle, bendType, Prob)
          call GenerateTorsAngle(tors_angle, torsType, Prob)
          call Generate_UnitTorsion(v1, v2, r, bend_angle, tors_angle, v3)
          trialPos%x = v3%x + molArray(nType)%mol(nMol)%x(Atm3)
          trialPos%y = v3%y + molArray(nType)%mol(nMol)%y(Atm3)
          trialPos%z = v3%z + molArray(nType)%mol(nMol)%z(Atm3)
          call Rosen_Atom_New(nType, Atm4, trialPos, isIncluded,  E_Trial(iRosen), overlap)
          call Rosen_Atom_Intra_Old(nType, nMol, Atm4, trialPos, regrown, E_Intra)
          E_Trial(iRosen) = E_Trial(iRosen) + E_Intra
        enddo
        E_Max = minval(E_Trial)
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
        enddo
        rosenNorm = sum(ProbRosen)
        rosenRatio = rosenRatio*ProbRosen(1)/rosenNorm
        regrown(Atm4) = .true.
      enddo

      call StraightChain_RosenConfigGen_GasPhase(nType, wForward, rejMove) 
      if (rejMove) then
        rosenRatio = 0d0
        return
      else
        rosenRatio = rosenRatio * wForward
      endif


      end subroutine

!=======================================================================

      subroutine StraightChain_RosenConfigGen_GasPhase(nType, wForward, rejMove)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions_LJ_Q
      use EnergyPressurePointers, only: Rosen_Atom_Intra_New
      use CBMC_Variables
      use CBMC_Utility
      use CoordinateFunctions
      implicit none

      integer, intent(in) :: nType
      real(dp) :: rosenRatio, wForward
      logical, intent(out) :: rejMove

!      logical :: isIncluded(1:maxMol)
      integer :: iRosen, nSel
      integer :: nRegrown
      integer :: Atm1, Atm2, Atm3, Atm4
      integer :: cnt
      integer :: bondType, bendType, torsType
      logical :: overlap(1:maxRosenTrial)
      logical :: regrown(1:maxAtoms)
      real(dp) :: E_Trial(1:maxRosenTrial), E_Complete
      real(dp) :: grnd
      real(dp) :: E_Max, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: ranNum, sumInt
      real(dp) :: k_bond, r_eq, r, Prob, r_sigma, rmax_sq
      real(dp) :: bend_angle, tors_angle
      real(dp) :: dx, dy, dz 
      type(SimpleAtomCoords) :: trialPos(1:maxRosenTrial)
      type(SimpleAtomCoords) :: v1, v2, v3
      
   
  
!      call Rosen_CreateSubset(nTarget, isIncluded)

      newMol%molType = nType 
      newMol%x = 0d0
      newMol%y = 0d0
      newMol%z = 0d0
      regrown = .false.
      rejMove = .false.
      nRegrown = 0

!      Since Growth is in gas phase, there is no need for choosing first atom position, so the first atom is located at the origin
      
      E_Trial = 0d0
      E_Complete = 0d0
      overlap = .false.

      rosenRatio = 1d0
      wForward = dble(nRosenTrials(nType))
      regrown(1) = .true.
      nRegrown = nRegrown + 1
      newMol%x(1) = 0d0
      newMol%y(1) = 0d0
      newMol%z(1) = 0d0

!      Having inserted the first atom, when we go to insert the second atom we must account for the fact that
!      the second atom must be inserted at the correct bond distance from the first atom.
      Atm1 = 1
      Atm2 = regrowOrder(nType, 2) 
      call FindBond(nType,Atm1, Atm2, bondType)
      k_bond = bondData(bondType)%k_eq
      r_eq = bondData(bondType)%r_eq
      r_sigma = bondData(bondType)%r_sigma
      rmax_sq = bondData(bondType)%rmax_sq
      overlap = .false.
      E_Trial = 0d0
      do iRosen = 1, nRosenTrials(nType)
        call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
        call Generate_UnitSphere(dx, dy, dz)
        trialPos(iRosen)%x = r*dx + newMol%x(Atm1) 
        trialPos(iRosen)%y = r*dy + newMol%y(Atm1)
        trialPos(iRosen)%z = r*dz + newMol%z(Atm1)
        E_Trial(iRosen) = 0d0
      enddo
!      E_Max = minval(E_Trial)
      E_Max = 0d0
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      if(all(ProbRosen .le. 0d0)) then
        rejMove = .true.
        return
      endif
      rosenNorm = sum(ProbRosen)
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbRosen(nSel)
      enddo
      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif
      rosenRatio = rosenRatio*ProbRosen(nSel)/rosenNorm
      wForward = wForward * rosenNorm
      regrown(Atm2) = .true.
      nRegrown = nRegrown + 1
      newMol%x(Atm2) = trialPos(nSel)%x 
      newMol%y(Atm2) = trialPos(nSel)%y 
      newMol%z(Atm2) = trialPos(nSel)%z

!       For the third atom we must begin to include the bending angle in choosing its trial positions. The first thing
!       we must consider is if the second regrown atom was a terminal atom or a linker atom. If it was a terminal atom 
!       then the bending angle must use the 1st atom as the central atom for the angle generation.  However if it was a link in the chain then the 2nd 
!       atom regrown should be used as the central atom.
      if(topolArray(nType)%atom(Atm2) .eq. 2) then
        Atm1 = regrowOrder(nType, 1) 
        Atm2 = regrowOrder(nType, 2) 
        Atm3 = regrowOrder(nType, 3) 
      else
        Atm1 = regrowOrder(nType, 2) 
        Atm2 = regrowOrder(nType, 1) 
        Atm3 = regrowOrder(nType, 3) 
      endif
      v1%x = newMol%x(Atm1) - newMol%x(Atm2)
      v1%y = newMol%y(Atm1) - newMol%y(Atm2)
      v1%z = newMol%z(Atm1) - newMol%z(Atm2)
      call FindBond(nType, Atm2, Atm3, bondType)
      k_bond = bondData(bondType)%k_eq
      r_eq = bondData(bondType)%r_eq
      r_sigma = bondData(bondType)%r_sigma
      rmax_sq = bondData(bondType)%rmax_sq
      call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
!      k_bend = bendData(bendType)%k_eq
!      ang_eq = bendData(bendType)%ang_eq
      overlap = .false.
      do iRosen = 1, nRosenTrials(nType)
        call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
!        call GenerateBendAngle(bend_angle, k_bend, ang_eq, Prob)
        call GenerateBendAngle(bend_angle, bendType, Prob)
        call Generate_UnitCone(v1, r, bend_angle, v2)
        trialPos(iRosen)%x = v2%x + newMol%x(Atm2) 
        trialPos(iRosen)%y = v2%y + newMol%y(Atm2)
        trialPos(iRosen)%z = v2%z + newMol%z(Atm2)
        E_Trial(iRosen) = 0d0
        call Rosen_Atom_Intra_New(nType, Atm3, trialPos(iRosen), regrown, E_Trial(iRosen), overlap(iRosen))
      enddo
!      E_Max = minval(E_Trial)
      E_Max = 0d0
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
!      if(all(ProbRosen .le. 0d0)) then
!        rejMove = .true.
!        return
!      endif
      rosenNorm = sum(ProbRosen)
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbRosen(nSel)
      enddo
      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif
!      rosenRatio = rosenRatio*ProbRosen(nSel)*dble(nRosenTrials(nType))/rosenNorm
      rosenRatio = rosenRatio*ProbRosen(nSel)/rosenNorm
      wForward = wForward * rosenNorm
      regrown(Atm3) = .true.
      nRegrown = nRegrown + 1
      newMol%x(Atm3) = trialPos(nSel)%x 
      newMol%y(Atm3) = trialPos(nSel)%y 
      newMol%z(Atm3) = trialPos(nSel)%z


!      Now that three atoms have been regrown, all remaining atoms must take the torsional angles into account.
      cnt = 3
      do while(any(regrown .eqv. .false.))
        cnt = cnt + 1
        atm4 = regrowOrder(nType, cnt) 
        call FindAtomsFromPath(nType, regrown, 1, Atm4, Atm1, Atm2, Atm3)
        call FindBond(nType, Atm3, Atm4, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        r_sigma = bondData(bondType)%r_sigma
        rmax_sq = bondData(bondType)%rmax_sq
        call FindAngle(nType, Atm2, Atm3, Atm4, bendType)
!        k_bend = bendData(bendType)%k_eq
!        ang_eq = bendData(bendType)%ang_eq
        call FindTorsion(nType, Atm1, Atm2, Atm3, Atm4, torsType)
        v1%x = newMol%x(Atm1) - newMol%x(Atm3)
        v1%y = newMol%y(Atm1) - newMol%y(Atm3)
        v1%z = newMol%z(Atm1) - newMol%z(Atm3)

        v2%x = newMol%x(Atm2) - newMol%x(Atm3)
        v2%y = newMol%y(Atm2) - newMol%y(Atm3)
        v2%z = newMol%z(Atm2) - newMol%z(Atm3)
        overlap = .false.
        E_Trial = 0d0
        do iRosen = 1, nRosenTrials(nType)
          call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
!          call GenerateBendAngle(bend_angle, k_bend, ang_eq, Prob)
          call GenerateBendAngle(bend_angle, bendType, Prob)
          call GenerateTorsAngle(tors_angle, torsType, Prob)
          call Generate_UnitTorsion(v1, v2, r, bend_angle, tors_angle, v3)
          trialPos(iRosen)%x = v3%x + newMol%x(Atm3) 
          trialPos(iRosen)%y = v3%y + newMol%y(Atm3)
          trialPos(iRosen)%z = v3%z + newMol%z(Atm3)
          E_Trial(iRosen) = 0d0
!          call Rosen_IntraNonBond_Atom_New(newMol2%x(:), newMol2%y(:), newMol2%z(:), nType, Atm4, trialPos(iRosen), &
!                                             regrown, E_Trial(iRosen), overlap(iRosen))
          call Rosen_Atom_Intra_New(nType, Atm4, trialPos(iRosen), regrown, &
                                                E_Trial(iRosen), overlap(iRosen))
        enddo
!        E_Max = minval(E_Trial)
        E_Max = 0d0
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
        enddo
!        if(all(ProbRosen .le. 0d0)) then
!          rejMove = .true.
!          return
!        endif
        rosenNorm = sum(ProbRosen)
        ranNum = grnd() * rosenNorm
        sumInt = ProbRosen(1)
        nSel = 1
        do while(sumInt .lt. ranNum)
          nSel = nSel + 1
          sumInt = sumInt + ProbRosen(nSel)
        enddo
        if(overlap(nSel) .eqv. .true.) then
          rejMove = .true.
          return
        endif
        rosenRatio = rosenRatio*ProbRosen(nSel)/rosenNorm
        wForward = wForward * rosenNorm
        regrown(Atm4) = .true.
        nRegrown = nRegrown + 1
        newMol%x(Atm4) = trialPos(nSel)%x 
        newMol%y(Atm4) = trialPos(nSel)%y 
        newMol%z(Atm4) = trialPos(nSel)%z
      enddo
	  

      end subroutine
!=======================================================================

      subroutine StraightChain_RosenConfigGen_GasPhase_Reverse(nType, wReverse)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions_LJ_Q
      use EnergyPressurePointers, only: Rosen_Atom_Intra_Gas_Old
      use CBMC_Variables
      use CBMC_Utility
      use CoordinateFunctions
      implicit none

      integer, intent(in) :: nType
      real(dp) :: wReverse

      integer :: iRosen
      integer :: Atm1, Atm2, Atm3, Atm4
      integer :: bondType, bendType, torsType, cnt
      logical :: overlap
      logical :: regrown(1:maxAtoms)
      real(dp) :: E_Trial(1:maxRosenTrial)
      real(dp) :: E_Max, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: k_bond, r_eq, r, Prob, r_sigma, rmax_sq
      real(dp) :: bend_angle, tors_angle
      type(SimpleAtomCoords) :: trialPos
      type(SimpleAtomCoords) :: v1, v2, v3
      
      ProbRosen = 0d0      
      E_Trial = 0d0      
      regrown = .false.
	  
!      Since the growth of first three atoms does not include nonBonded intera-molecular energies, we can write

      wReverse = dble(nRosenTrials(nType)) * dble(nRosenTrials(nType)) * dble(nRosenTrials(nType))


!       For the third atom we must begin to include the bending angle in choosing its trial positions. The first thing
!       we must consider is if the second regrown atom was a terminal atom or a linker atom. If it was a terminal atom 
!       then the bending angle must use the 1st atom as the central atom for the angle generation.  However if it was a link in the chain then the 2nd 
!       atom regrown should be used as the central atom.
      if(topolArray(nType)%atom(Atm2) .eq. 2) then
        Atm1 = regrowOrder(nType, 1) 
        Atm2 = regrowOrder(nType, 2) 
        Atm3 = regrowOrder(nType, 3) 
      else
        Atm1 = regrowOrder(nType, 2) 
        Atm2 = regrowOrder(nType, 1) 
        Atm3 = regrowOrder(nType, 3) 
      endif
      regrown(Atm1) = .true.
      regrown(Atm2) = .true.
      regrown(Atm3) = .true.

!      Now that three atoms have been regrown, all remaining atoms must take the torsional angles into account.
      cnt = 3
      do while(any(regrown .eqv. .false.))
        cnt = cnt + 1
        atm4 = regrowOrder(nType, cnt) 
        call FindAtomsFromPath(nType, regrown, 1, Atm4, Atm1, Atm2, Atm3)
        call FindBond(nType, Atm3, Atm4, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        r_sigma = bondData(bondType)%r_sigma
        rmax_sq = bondData(bondType)%rmax_sq
        call FindAngle(nType, Atm2, Atm3, Atm4, bendType)
!        k_bend = bendData(bendType)%k_eq
!        ang_eq = bendData(bendType)%ang_eq
        call FindTorsion(nType, Atm1, Atm2, Atm3, Atm4, torsType)
        v1%x = gasConfig(nType)%x(Atm1) - gasConfig(nType)%x(Atm3)
        v1%y = gasConfig(nType)%y(Atm1) - gasConfig(nType)%y(Atm3)
        v1%z = gasConfig(nType)%z(Atm1) - gasConfig(nType)%z(Atm3)

        v2%x = gasConfig(nType)%x(Atm2) - gasConfig(nType)%x(Atm3)
        v2%y = gasConfig(nType)%y(Atm2) - gasConfig(nType)%y(Atm3)
        v2%z = gasConfig(nType)%z(Atm2) - gasConfig(nType)%z(Atm3)
        overlap = .false.
        E_Trial(1) = 0d0
        trialPos%x = gasConfig(nType)%x(Atm4)
        trialPos%y = gasConfig(nType)%y(Atm4)
        trialPos%z = gasConfig(nType)%z(Atm4)
        call Rosen_Atom_Intra_Gas_Old(nType, Atm4, trialPos, regrown, E_Trial(1))
!        call Rosen_IntraNonBond_Atom_Old(gasConfig(nType)%x(:), gasConfig(nType)%y(:), &
!                     gasConfig(nType)%z(:), nType, Atm4, regrown, E_Trial(1))
        do iRosen = 2, nRosenTrials(nType)
          call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
!          call GenerateBendAngle(bend_angle, k_bend, ang_eq, Prob)
          call GenerateBendAngle(bend_angle, bendType, Prob)
          call GenerateTorsAngle(tors_angle, torsType, Prob)
          call Generate_UnitTorsion(v1, v2, r, bend_angle, tors_angle, v3)
          trialPos%x = v3%x + gasConfig(nType)%x(Atm3)
          trialPos%y = v3%y + gasConfig(nType)%y(Atm3)
          trialPos%z = v3%z + gasConfig(nType)%z(Atm3)
          E_Trial(iRosen) = 0d0
          call Rosen_Atom_Intra_Gas_Old(nType, Atm4, trialPos, regrown, E_Trial(iRosen))
!          call Rosen_IntraNonBond_Atom_New(gasConfig(nType)%x(:), gasConfig(nType)%y(:), &
!                       gasConfig(nType)%z(:), nType, Atm4, trialPos, regrown, E_Trial(iRosen), overlap)
        enddo
!        E_Max = minval(E_Trial)
        E_Max = 0d0
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
        enddo
        rosenNorm = sum(ProbRosen)
        wReverse = wReverse*rosenNorm
        regrown(Atm4) = .true.
      enddo


      end subroutine
!=======================================================================
      subroutine BranchedMol_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions_LJ_Q
      use EnergyPressurePointers, only: Rosen_Atom_New, Rosen_Atom_Intra_New
      use CBMC_Variables
      use CBMC_Utility
      use CoordinateFunctions
      implicit none
	  
      integer, intent(in) :: nType, nTarget, nTargType
      real(dp), intent(out):: rosenRatio
      logical, intent(out) :: rejMove
      logical, intent(in) :: isIncluded(:)
	  
      integer :: iAtom, jAtom, iRosen, nIndx, iGrow, nTargetMol, nSel, totalRegrown
      integer :: atmCur, atmPrev, atmGrow(1:3), atmTor(1:3), nToGrow, nTor
      integer :: bondType, bendTypes(1:6), torsTypes(1:3,1:3)
      logical :: regrown(1:maxAtoms)
      logical :: overlap(1:maxRosenTrial)
	  
      real(dp) :: E_Trial(1:maxRosenTrial), E_Trial_i, wBending(1:maxRosenTrial), wTorsion(1:maxRosenTrial), E_Max
      real(dp) :: grnd, ranNum, sumInt
      real(dp) :: k_bonds(1:3), r_eqs(1:3), r, rGrow(1:3), Prob, ProbRosen(1:maxRosenTrial), rosenNorm, wReverse
      real(dp) :: bendingAngles(1:3), dihedralAngles(1:2), r_sigmas(1:3), rmax_sqs(1:3)
      real(dp) :: dx, dy, dz 
      type(SimpleAtomCoords) :: trialPos(1:maxRosenTrial), trialPos_Branched(1:maxRosenTrial,1:3)
      type(SimpleAtomCoords) :: vPrev, vTor(1:3), vGrow(1:3)
	  
      nTargetMol = subIndxList(nTarget)
      newMol%molType = nType 
      newMol%x = 0d0
      newMol%y = 0d0
      newMol%z = 0d0
      regrown = .false.
      rejMove = .false.
      totalRegrown = 0

!      Begin the regrowth process by choosing an insertion site for the first atom in the chain
      
      E_Trial = 0d0
      overlap = .false.
      do iRosen = 1, nRosenTrials(nType)
        r = Dist_Critr * grnd()**(1d0/3d0)
        call Generate_UnitSphere(dx, dy, dz)
        trialPos(iRosen)%x = r * dx + molArray(nTargType)%mol(nTargetMol)%x(1)
        trialPos(iRosen)%y = r * dy + molArray(nTargType)%mol(nTargetMol)%y(1)
        trialPos(iRosen)%z = r * dz + molArray(nTargType)%mol(nTargetMol)%z(1)
        call Rosen_Atom_New(nType, 1, trialPos(iRosen), isIncluded,  E_Trial(iRosen), overlap(iRosen))
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      if(all(ProbRosen .le. 0d0)) then
        rejMove = .true.
        return
      endif
      rosenNorm = sum(ProbRosen)
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbRosen(nSel)
      enddo
      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif
      rosenRatio = ProbRosen(nSel)/rosenNorm
      regrown(1) = .true.
      totalRegrown = totalRegrown + 1
      newMol%x(1) = trialPos(nSel)%x 
      newMol%y(1) = trialPos(nSel)%y 
      newMol%z(1) = trialPos(nSel)%z
	  
!      Having inserted the first atom, when we go to insert the second atom we must account for the fact that
!      the second atom must be inserted at the correct bond distance from the first atom.
      atmPrev = 1
      atmCur = 2
      call FindBond(nType, atmPrev, atmCur, bondType)
      k_bonds(1) = bondData(bondType)%k_eq
      r_eqs(1) = bondData(bondType)%r_eq
      r_sigmas(1) = bondData(bondType)%r_sigma
      rmax_sqs(1) = bondData(bondType)%rmax_sq
      overlap = .false.
      E_Trial = 0d0
      do iRosen = 1, nRosenTrials(nType)
        call GenerateGaussianBondLength(r, k_bonds(1), r_eqs(1), r_sigmas(1), rmax_sqs(1), Prob)
        call Generate_UnitSphere(dx, dy, dz)
        trialPos(iRosen)%x = r*dx + newMol%x(atmPrev) 
        trialPos(iRosen)%y = r*dy + newMol%y(atmPrev)
        trialPos(iRosen)%z = r*dz + newMol%z(atmPrev)
        call Rosen_Atom_New(nType, atmCur, trialPos(iRosen), isIncluded,  E_Trial(iRosen), overlap(iRosen))
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      if(all(ProbRosen .le. 0d0)) then
        rejMove = .true.
        return
      endif
      rosenNorm = sum(ProbRosen)
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbRosen(nSel)
      enddo
      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif
      rosenRatio = rosenRatio*ProbRosen(nSel)/rosenNorm
      regrown(atmCur) = .true.
      totalRegrown = totalRegrown + 1
      newMol%x(atmCur) = trialPos(nSel)%x 
      newMol%y(atmCur) = trialPos(nSel)%y 
      newMol%z(atmCur) = trialPos(nSel)%z
	  
	  
      do iGrow = 1, SwapGrowOrder(nType)%GrowthSteps
        atmGrow = 0
        atmTor = 0
        nToGrow = SwapGrowOrder(nType)%GrowNum(iGrow)
        nTor = SwapGrowOrder(nType)%TorNum(iGrow)
        atmPrev = SwapGrowOrder(nType)%GrowPrev(iGrow)
        atmCur = SwapGrowOrder(nType)%GrowFrom(iGrow)
        vPrev%x = newMol%x(atmPrev) - newMol%x(atmCur)
        vPrev%y = newMol%y(atmPrev) - newMol%y(atmCur)
        vPrev%z = newMol%z(atmPrev) - newMol%z(atmCur)
        do iAtom = 1, nToGrow
          atmGrow(iAtom) = SwapGrowOrder(nType)%GrowList(iGrow,iAtom)
          call FindBond(nType, atmCur, atmGrow(iAtom), bondType)
          k_bonds(iAtom) = bondData(bondType)%k_eq
          r_eqs(iAtom) = bondData(bondType)%r_eq
          r_sigmas(iAtom) = bondData(bondType)%r_sigma
          rmax_sqs(iAtom) = bondData(bondType)%rmax_sq
          call FindAngle(nType, atmPrev, atmCur, atmGrow(iAtom), bendTypes(iAtom))
          do jAtom = 1, nTor
            if (iAtom .eq. 1) then
              atmTor(jAtom) = SwapGrowOrder(nType)%TorList(iGrow,jAtom)
              vTor(jAtom)%x = newMol%x(atmTor(jAtom)) - newMol%x(atmCur)
              vTor(jAtom)%y = newMol%y(atmTor(jAtom)) - newMol%y(atmCur)
              vTor(jAtom)%z = newMol%z(atmTor(jAtom)) - newMol%z(atmCur)
            endif
            call FindTorsion(nType, atmTor(jAtom), atmPrev, atmCur, atmGrow(iAtom), torsTypes(iAtom,jAtom))
          enddo
        enddo
        if (nToGrow .eq. 2) then
          call FindAngle(nType, atmGrow(1), atmCur, atmGrow(2), bendTypes(3))
        elseif (nToGrow .eq. 3) then
          call FindAngle(nType, atmGrow(1), atmCur, atmGrow(2), bendTypes(4))
          call FindAngle(nType, atmGrow(2), atmCur, atmGrow(3), bendTypes(5))
          call FindAngle(nType, atmGrow(3), atmCur, atmGrow(1), bendTypes(6))
        endif
        overlap = .false.
        do iRosen = 1, nRosenTrials(nType)
          bendingAngles = 0d0
          dihedralAngles = 0d0
          do iAtom = 1, nToGrow
            call GenerateGaussianBondLength(rGrow(iAtom), k_bonds(iAtom), r_eqs(iAtom), &
                                    r_sigmas(iAtom), rmax_sqs(iAtom), Prob)
          enddo
          select case(nToGrow)
          case(1)
            call GenerateBendAngle(bendingAngles(1), bendTypes(1), Prob)
            wBending(iRosen) = 1d0
          case(2)
            wBending(iRosen) = 0d0
            call GenerateTwoBranches(bendingAngles(1), bendingAngles(2), dihedralAngles(1), & 
                                     bendTypes(1), bendTypes(2), bendTypes(3), wBending(iRosen))
          case(3)
            wBending(iRosen) = 0d0
            call GenerateThreeBranches(bendingAngles(1), bendingAngles(1), bendingAngles(3), &
                                       dihedralAngles(1), dihedralAngles(2), bendTypes(1), bendTypes(2), bendTypes(3), &
                                       bendTypes(4), bendTypes(5), bendTypes(6), wBending(iRosen))
          case default
            write(*,*) "Error! Atom", atmCur, "in molecule type", nType, "has too many branches for a simple regrowth"
            stop
          end select
          call GenerateRotation(nTor, nToGrow, vTor, vPrev, rGrow, bendingAngles, dihedralAngles, &
                                torsTypes, vGrow, wTorsion(iRosen))
          do iAtom = 1, nToGrow
            trialPos_Branched(iRosen,iAtom)%x = vGrow(iAtom)%x + newMol%x(atmCur)
            trialPos_Branched(iRosen,iAtom)%y = vGrow(iAtom)%y + newMol%y(atmCur)
            trialPos_Branched(iRosen,iAtom)%z = vGrow(iAtom)%z + newMol%z(atmCur)
          enddo
          E_Trial(iRosen) = 0d0
          do iAtom = 1, nToGrow
            call Rosen_Atom_New(nType, atmGrow(iAtom), trialPos_Branched(iRosen,iAtom), isIncluded,  E_Trial_i, &
                                            overlap(iRosen))
            E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_i
!            if (nIntraNonBond(nType) .gt. 0) then
!              call Rosen_IntraNonBond_Atom_New(newMol%x(:), newMol%y(:), newMol%z(:), nType, atmGrow(iAtom), &
!                                               trialPos_Branched(iRosen,iAtom), regrown, E_Trial(iRosen), overlap(iRosen))
!            endif
            call Rosen_Atom_Intra_New(nType, atmGrow(iAtom), trialPos_Branched(iRosen,iAtom), &
                                                  regrown, E_Trial_i, overlap(iRosen))
            E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_i
          enddo
        enddo
        E_Max = minval(E_Trial)
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = wBending(iRosen) * wTorsion(iRosen) * exp(-beta*(E_Trial(iRosen)-E_Max))         
        enddo
        rosenNorm = sum(ProbRosen)
        ranNum = grnd() * rosenNorm
        sumInt = ProbRosen(1)
        nSel = 1
        do while(sumInt .lt. ranNum)
          nSel = nSel + 1
          sumInt = sumInt + ProbRosen(nSel)
        enddo
        if(overlap(nSel) .eqv. .true.) then
          rejMove = .true.
          return
        endif
        rosenRatio = rosenRatio*ProbRosen(nSel)/rosenNorm
        do iAtom = 1, nToGrow
          regrown(atmGrow(iAtom)) = .true.
          totalRegrown = totalRegrown + 1
          newMol%x(atmGrow(iAtom)) = trialPos_Branched(nSel,iAtom)%x 
          newMol%y(atmGrow(iAtom)) = trialPos_Branched(nSel,iAtom)%y 
          newMol%z(atmGrow(iAtom)) = trialPos_Branched(nSel,iAtom)%z
        enddo
      enddo
	  
      call BranchedMol_RosenConfigGen_GasPhase_Reverse(nType, wReverse)
      rosenRatio = rosenRatio * wReverse
	  
	  
      end subroutine
!=======================================================================
      subroutine BranchedMol_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions_LJ_Q
      use EnergyPressurePointers, only: Rosen_Atom_New, Rosen_Atom_Old, Rosen_Atom_Intra_Old
      use CBMC_Variables
      use CBMC_Utility
      use CoordinateFunctions
      implicit none
	  
      integer, intent(in) :: nType, nMol, nTarget, nTargType
      real(dp), intent(out):: rosenRatio
	  
      integer :: iAtom, jAtom, iRosen, nIndx, iGrow, nTargetMol
      integer :: atmCur, atmPrev, atmGrow(1:3), atmTor(1:3), nToGrow, nTor, atm2
      integer :: bondType, bendTypes(1:6), torsTypes(1:3,1:3)
      logical :: regrown(1:maxAtoms)
      logical :: overlap, rejMove, isIncluded(1:maxMol)
	  
      real(dp) :: E_Trial(1:maxRosenTrial), E_Trial_i, wBending(1:maxRosenTrial), wTorsion(1:maxRosenTrial), E_Max
      real(dp) :: grnd
      real(dp) :: k_bonds(1:3), r_eqs(1:3), r, rGrow(1:3), Prob, ProbRosen(1:maxRosenTrial), rosenNorm, wForward
      real(dp) :: bendingAngles(1:3), dihedralAngles(1:2), r_sigmas(1:3), rmax_sqs(1:3)
      real(dp) :: dx, dy, dz 
      type(SimpleAtomCoords) :: trialPos, trialPos_Branched(1:3)
      type(SimpleAtomCoords) :: vPrev, vTor(1:3), vGrow(1:3)
	  
      ProbRosen = 0d0      
      E_Trial = 0d0      
      nTargetMol = subIndxList(nTarget)
      nIndx = molArray(nType)%mol(nMol)%indx
      call Rosen_CreateSubset_Reverse(nTarget, nIndx, isIncluded)
      regrown = .false.
	  
!      Begin the regrowth process by choosing an insertion site for the first atom in the chain
      E_Trial = 0d0
      rosenRatio = 1d0
      call Rosen_Atom_Old(nType, nMol, 1, isIncluded,  E_Trial(1))
      do iRosen = 2, nRosenTrials(nType)
        r = Dist_Critr * grnd()**(1d0/3d0)
        call Generate_UnitSphere(dx, dy, dz)
        trialPos%x = r * dx + molArray(nTargType)%mol(nTargetMol)%x(1)
        trialPos%y = r * dy + molArray(nTargType)%mol(nTargetMol)%y(1)
        trialPos%z = r * dz + molArray(nTargType)%mol(nTargetMol)%z(1)
        call Rosen_Atom_New(nType, 1, trialPos, isIncluded,  E_Trial(iRosen), overlap)
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      rosenNorm = sum(ProbRosen)
      rosenRatio = rosenRatio*ProbRosen(1)/rosenNorm
      regrown(1) = .true.
	  
!      Having inserted the first atom, when we go to insert the second atom we must account for the fact that
!      the second atom must be inserted at the correct bond distance from the first atom.
      atmPrev = 1
      atmCur = 2
      call FindBond(nType, atmPrev, atmCur, bondType)
      k_bonds(1) = bondData(bondType)%k_eq
      r_eqs(1) = bondData(bondType)%r_eq
      r_sigmas(1) = bondData(bondType)%r_sigma
      rmax_sqs(1) = bondData(bondType)%rmax_sq
      E_Trial = 0d0
      call Rosen_Atom_Old(nType, nMol, Atm2, isIncluded,  E_Trial(1))
      do iRosen = 2, nRosenTrials(nType)
        call GenerateGaussianBondLength(r, k_bonds(1), r_eqs(1), r_sigmas(1), rmax_sqs(1), Prob)
        call Generate_UnitSphere(dx, dy, dz)
        trialPos%x = r*dx + newMol%x(atmPrev) 
        trialPos%y = r*dy + newMol%y(atmPrev)
        trialPos%z = r*dz + newMol%z(atmPrev)
        call Rosen_Atom_New(nType, atmCur, trialPos, isIncluded,  E_Trial(iRosen), overlap)
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      rosenNorm = sum(ProbRosen)
      rosenRatio = rosenRatio*ProbRosen(1)/rosenNorm
      regrown(atmCur) = .true.
	  
	  
      do iGrow = 1, SwapGrowOrder(nType)%GrowthSteps
        atmGrow = 0
        atmTor = 0
        nToGrow = SwapGrowOrder(nType)%GrowNum(iGrow)
        nTor = SwapGrowOrder(nType)%TorNum(iGrow)
        atmPrev = SwapGrowOrder(nType)%GrowPrev(iGrow)
        atmCur = SwapGrowOrder(nType)%GrowFrom(iGrow)
        vPrev%x = molArray(nType)%mol(nMol)%x(atmPrev) - molArray(nType)%mol(nMol)%x(atmCur)
        vPrev%y = molArray(nType)%mol(nMol)%y(atmPrev) - molArray(nType)%mol(nMol)%y(atmCur)
        vPrev%z = molArray(nType)%mol(nMol)%z(atmPrev) - molArray(nType)%mol(nMol)%z(atmCur)
        do iAtom = 1, nToGrow
          atmGrow(iAtom) = SwapGrowOrder(nType)%GrowList(iGrow,iAtom)
          vGrow(iAtom)%x = molArray(nType)%mol(nMol)%x(atmGrow(iAtom)) - molArray(nType)%mol(nMol)%x(atmCur)
          vGrow(iAtom)%y = molArray(nType)%mol(nMol)%y(atmGrow(iAtom)) - molArray(nType)%mol(nMol)%y(atmCur)
          vGrow(iAtom)%z = molArray(nType)%mol(nMol)%z(atmGrow(iAtom)) - molArray(nType)%mol(nMol)%z(atmCur)
          call FindBond(nType, atmCur, atmGrow(iAtom), bondType)
          k_bonds(iAtom) = bondData(bondType)%k_eq
          r_eqs(iAtom) = bondData(bondType)%r_eq
          r_sigmas(iAtom) = bondData(bondType)%r_sigma
          rmax_sqs(iAtom) = bondData(bondType)%rmax_sq
          call FindAngle(nType, atmPrev, atmCur, atmGrow(iAtom), bendTypes(iAtom))
          do jAtom = 1, nTor
            if (iAtom .eq. 1) then
              atmTor(jAtom) = SwapGrowOrder(nType)%TorList(iGrow,jAtom)
              vTor(jAtom)%x = molArray(nType)%mol(nMol)%x(atmTor(jAtom)) - molArray(nType)%mol(nMol)%x(atmCur)
              vTor(jAtom)%y = molArray(nType)%mol(nMol)%y(atmTor(jAtom)) - molArray(nType)%mol(nMol)%y(atmCur)
              vTor(jAtom)%z = molArray(nType)%mol(nMol)%z(atmTor(jAtom)) - molArray(nType)%mol(nMol)%z(atmCur)
            endif
            call FindTorsion(nType, atmTor(jAtom), atmPrev, atmCur, atmGrow(iAtom), torsTypes(iAtom,jAtom))
          enddo
        enddo
        E_Trial(1) = 0d0
        do iAtom = 1, nToGrow
          call Rosen_Atom_Old(nType, nMol, atmGrow(iAtom), isIncluded,  E_Trial_i)
          E_Trial(1) = E_Trial(1) + E_Trial_i
          trialPos_Branched(iAtom)%x = molArray(nType)%mol(nMol)%x(atmGrow(iAtom))
          trialPos_Branched(iAtom)%y = molArray(nType)%mol(nMol)%y(atmGrow(iAtom))
          trialPos_Branched(iAtom)%z = molArray(nType)%mol(nMol)%z(atmGrow(iAtom))
          call Rosen_Atom_Intra_Old(nType, nMol, atmGrow(iAtom), trialPos_Branched(iAtom), &
                                                regrown, E_Trial_i)
          E_Trial(1) = E_Trial(1) + E_Trial_i
!          if (nIntraNonBond(nType) .gt. 0) then
!            call Rosen_IntraNonBond_Atom_Old(molArray(nType)%mol(nMol)%x(:), molArray(nType)%mol(nMol)%y(:), &
!                       molArray(nType)%mol(nMol)%z(:), nType, atmGrow(iAtom), regrown, E_Trial(1))
!          endif
        enddo
        select case(nToGrow)
        case(1)
          wBending(1) = 1d0
        case(2)
          call FindAngle(nType, atmGrow(1), atmCur, atmGrow(2), bendTypes(3))
          wBending(1) = 0d0
          call GenerateTwoBranches_Reverse(vPrev, vGrow(1), vGrow(2), bendTypes(1), bendTypes(2), bendTypes(3), wBending(1))
        case(3)
          call FindAngle(nType, atmGrow(1), atmCur, atmGrow(2), bendTypes(4))
          call FindAngle(nType, atmGrow(2), atmCur, atmGrow(3), bendTypes(5))
          call FindAngle(nType, atmGrow(3), atmCur, atmGrow(1), bendTypes(6))
          wBending(1) = 0d0
          call GenerateThreeBranches_Reverse(vPrev, vGrow(1) ,vGrow(2) ,vGrow(3) , bendTypes(1), bendTypes(2), bendTypes(3), &
                                             bendTypes(4), bendTypes(5), bendTypes(6), wBending(1))
        case default
          write(*,*) "Error! Atom", atmCur, "in molecule type", nType, "has too many branches for a simple regrowth"
          stop
        end select
        call GenerateRotation_Reverse(nTor, nToGrow, vTor, vPrev, torsTypes, vGrow, wTorsion(1))
        do iRosen = 2, nRosenTrials(nType)
          bendingAngles = 0d0
          dihedralAngles = 0d0
          do iAtom = 1, nToGrow
            call GenerateGaussianBondLength(rGrow(iAtom), k_bonds(iAtom), r_eqs(iAtom), &
                                    r_sigmas(iAtom), rmax_sqs(iAtom), Prob)
          enddo
          select case(nToGrow)
          case(1)
            call GenerateBendAngle(bendingAngles(1), bendTypes(1), Prob)
            wBending(iRosen) = 1d0
          case(2)
            wBending(iRosen) = 0d0
            call GenerateTwoBranches(bendingAngles(1), bendingAngles(2), dihedralAngles(1), & 
                                     bendTypes(1), bendTypes(2), bendTypes(3), wBending(iRosen))
          case(3)
            wBending(iRosen) = 0d0
            call GenerateThreeBranches(bendingAngles(1), bendingAngles(1), bendingAngles(3), &
                                       dihedralAngles(1), dihedralAngles(2), bendTypes(1), bendTypes(2), bendTypes(3), &
                                       bendTypes(4), bendTypes(5), bendTypes(6), wBending(iRosen))
          case default
            write(*,*) "Error! Atom", atmCur, "in molecule type", nType, "has too many branches for a simple regrowth"
            stop
          end select
          call GenerateRotation(nTor, nToGrow, vTor, vPrev, rGrow, bendingAngles, dihedralAngles, &
                                torsTypes, vGrow, wTorsion(iRosen))
          do iAtom = 1, nToGrow
            trialPos_Branched(iAtom)%x = vGrow(iAtom)%x + molArray(nType)%mol(nMol)%x(atmCur)
            trialPos_Branched(iAtom)%y = vGrow(iAtom)%y + molArray(nType)%mol(nMol)%y(atmCur)
            trialPos_Branched(iAtom)%z = vGrow(iAtom)%z + molArray(nType)%mol(nMol)%z(atmCur)
          enddo
          E_Trial(iRosen) = 0d0
          do iAtom = 1, nToGrow
            call Rosen_Atom_New(nType, atmGrow(iAtom), trialPos_Branched(iAtom), isIncluded,  E_Trial_i, overlap)
            E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_i
            call Rosen_Atom_Intra_Old(nType, nMol, atmGrow(iAtom), trialPos_Branched(iAtom), &
                                                regrown, E_Trial_i)
            E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_i
!            if (nIntraNonBond(nType) .gt. 0) then
!              call Rosen_IntraNonBond_Atom_New(molArray(nType)%mol(nMol)%x(:), molArray(nType)%mol(nMol)%y(:), &
!              molArray(nType)%mol(nMol)%z(:), nType, atmGrow(iAtom), trialPos_Branched(iAtom), regrown, E_Trial(iRosen), overlap)
!            endif
          enddo
        enddo
        E_Max = minval(E_Trial)
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = wBending(iRosen) * wTorsion(iRosen) * exp(-beta*(E_Trial(iRosen)-E_Max))         
        enddo
        rosenNorm = sum(ProbRosen)
        rosenRatio = rosenRatio*ProbRosen(1)/rosenNorm
        do iAtom = 1, nToGrow
          regrown(atmGrow(iAtom)) = .true.
        enddo
      enddo
      call BranchedMol_RosenConfigGen_GasPhase(nType, wForward, rejMove)
      if (rejMove .eqv. .true.) then
        rosenRatio = 0d0
        return
      endif
      rosenRatio = rosenRatio * wForward

	  
      end subroutine
	  
!=======================================================================
      subroutine BranchedMol_RosenConfigGen_GasPhase(nType, wForward, rejMove)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions_LJ_Q
      use EnergyPressurePointers, only: Rosen_Atom_Intra_New
      use CBMC_Variables
      use CBMC_Utility
      use CoordinateFunctions
      implicit none
	  
      integer, intent(in) :: nType
      real(dp) :: rosenRatio, wForward
      logical, intent(out) :: rejMove
	  
      integer :: iAtom, jAtom, iRosen, iGrow, nSel, totalRegrown
      integer :: atmCur, atmPrev, atmGrow(1:3), atmTor(1:3), nToGrow, nTor
      integer :: bondType, bendTypes(1:6), torsTypes(1:3,1:3)
      logical :: regrown(1:maxAtoms)
      logical :: overlap(1:maxRosenTrial)
	  
      real(dp) :: E_Trial(1:maxRosenTrial), E_Trial_i, wBending(1:maxRosenTrial), wTorsion(1:maxRosenTrial), E_Max
      real(dp) :: grnd, ranNum, sumInt
      real(dp) :: k_bonds(1:3), r_eqs(1:3), r, rGrow(1:3), Prob, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: bendingAngles(1:3), dihedralAngles(1:2), r_sigmas(1:3), rmax_sqs(1:3)
      real(dp) :: dx, dy, dz 
      type(SimpleAtomCoords) :: trialPos(1:maxRosenTrial), trialPos_Branched(1:maxRosenTrial,1:3)
      type(SimpleAtomCoords) :: vPrev, vTor(1:3), vGrow(1:3)
	  
      newMol%molType = nType 
      newMol%x = 0d0
      newMol%y = 0d0
      newMol%z = 0d0
      regrown = .false.
      rejMove = .false.
      totalRegrown = 0
	  	  
!      Since Growth is in gas phase, there is no need for choosing first atom position, so the first atom is located at the origin
      
      E_Trial = 0d0
!      E_Complete = 0d0
      overlap = .false.

      rosenRatio = 1d0
      wForward = dble(nRosenTrials(nType))
      regrown(1) = .true.
      totalRegrown = totalRegrown + 1
      newMol%x(1) = 0d0
      newMol%y(1) = 0d0
      newMol%z(1) = 0d0
	  
!      Having inserted the first atom, when we go to insert the second atom we must account for the fact that
!      the second atom must be inserted at the correct bond distance from the first atom.
      atmPrev = 1
      atmCur = 2
      call FindBond(nType, atmPrev, atmCur, bondType)
      k_bonds(1) = bondData(bondType)%k_eq
      r_eqs(1) = bondData(bondType)%r_eq
      r_sigmas(1) = bondData(bondType)%r_sigma
      rmax_sqs(1) = bondData(bondType)%rmax_sq
      overlap = .false.
      E_Trial = 0d0
      do iRosen = 1, nRosenTrials(nType)
        call GenerateGaussianBondLength(r, k_bonds(1), r_eqs(1), r_sigmas(1), rmax_sqs(1), Prob)
        call Generate_UnitSphere(dx, dy, dz)
        trialPos(iRosen)%x = r*dx + newMol%x(atmPrev) 
        trialPos(iRosen)%y = r*dy + newMol%y(atmPrev)
        trialPos(iRosen)%z = r*dz + newMol%z(atmPrev)
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      if(all(ProbRosen .le. 0d0)) then
        rejMove = .true.
        return
      endif
      rosenNorm = sum(ProbRosen)
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbRosen(nSel)
      enddo
      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif
      rosenRatio = rosenRatio*ProbRosen(nSel)/rosenNorm
      wForward = wForward * rosenNorm
      regrown(atmCur) = .true.
      totalRegrown = totalRegrown + 1
      newMol%x(atmCur) = trialPos(nSel)%x 
      newMol%y(atmCur) = trialPos(nSel)%y 
      newMol%z(atmCur) = trialPos(nSel)%z

      do iGrow = 1, SwapGrowOrder(nType)%GrowthSteps
        atmGrow = 0
        atmTor = 0
        nToGrow = SwapGrowOrder(nType)%GrowNum(iGrow)
        nTor = SwapGrowOrder(nType)%TorNum(iGrow)
        atmPrev = SwapGrowOrder(nType)%GrowPrev(iGrow)
        atmCur = SwapGrowOrder(nType)%GrowFrom(iGrow)
        vPrev%x = newMol%x(atmPrev) - newMol%x(atmCur)
        vPrev%y = newMol%y(atmPrev) - newMol%y(atmCur)
        vPrev%z = newMol%z(atmPrev) - newMol%z(atmCur)
        do iAtom = 1, nToGrow
          atmGrow(iAtom) = SwapGrowOrder(nType)%GrowList(iGrow,iAtom)
          call FindBond(nType, atmCur, atmGrow(iAtom), bondType)
          k_bonds(iAtom) = bondData(bondType)%k_eq
          r_eqs(iAtom) = bondData(bondType)%r_eq
          r_sigmas(iAtom) = bondData(bondType)%r_sigma
          rmax_sqs(iAtom) = bondData(bondType)%rmax_sq
          call FindAngle(nType, atmPrev, atmCur, atmGrow(iAtom), bendTypes(iAtom))
          do jAtom = 1, nTor
            if (iAtom .eq. 1) then
              atmTor(jAtom) = SwapGrowOrder(nType)%TorList(iGrow,jAtom)
              vTor(jAtom)%x = newMol%x(atmTor(jAtom)) - newMol%x(atmCur)
              vTor(jAtom)%y = newMol%y(atmTor(jAtom)) - newMol%y(atmCur)
              vTor(jAtom)%z = newMol%z(atmTor(jAtom)) - newMol%z(atmCur)
            endif
            call FindTorsion(nType, atmTor(jAtom), atmPrev, atmCur, atmGrow(iAtom), torsTypes(iAtom,jAtom))
          enddo
        enddo
        if (nToGrow .eq. 2) then
          call FindAngle(nType, atmGrow(1), atmCur, atmGrow(2), bendTypes(3))
        elseif (nToGrow .eq. 3) then
          call FindAngle(nType, atmGrow(1), atmCur, atmGrow(2), bendTypes(4))
          call FindAngle(nType, atmGrow(2), atmCur, atmGrow(3), bendTypes(5))
          call FindAngle(nType, atmGrow(3), atmCur, atmGrow(1), bendTypes(6))
        endif
        overlap = .false.
        do iRosen = 1, nRosenTrials(nType)
          bendingAngles = 0d0
          dihedralAngles = 0d0
          do iAtom = 1, nToGrow
            call GenerateGaussianBondLength(rGrow(iAtom), k_bonds(iAtom), r_eqs(iAtom), &
                                    r_sigmas(iAtom), rmax_sqs(iAtom), Prob)
          enddo
          select case(nToGrow)
          case(1)
            call GenerateBendAngle(bendingAngles(1), bendTypes(1), Prob)
            wBending(iRosen) = 1d0
          case(2)
            wBending(iRosen) = 0d0
            call GenerateTwoBranches(bendingAngles(1), bendingAngles(2), dihedralAngles(1), & 
                                     bendTypes(1), bendTypes(2), bendTypes(3), wBending(iRosen))
          case(3)
            wBending(iRosen) = 0d0
            call GenerateThreeBranches(bendingAngles(1), bendingAngles(1), bendingAngles(3), &
                                       dihedralAngles(1), dihedralAngles(2), bendTypes(1), bendTypes(2), bendTypes(3), &
                                       bendTypes(4), bendTypes(5), bendTypes(6), wBending(iRosen))
          case default
            write(*,*) "Error! Atom", atmCur, "in molecule type", nType, "has too many branches for a simple regrowth"
            stop
          end select
          call GenerateRotation(nTor, nToGrow, vTor, vPrev, rGrow, bendingAngles, dihedralAngles, &
                                torsTypes, vGrow, wTorsion(iRosen))
          do iAtom = 1, nToGrow
            trialPos_Branched(iRosen,iAtom)%x = vGrow(iAtom)%x + newMol%x(atmCur)
            trialPos_Branched(iRosen,iAtom)%y = vGrow(iAtom)%y + newMol%y(atmCur)
            trialPos_Branched(iRosen,iAtom)%z = vGrow(iAtom)%z + newMol%z(atmCur)
          enddo
          E_Trial(iRosen) = 0d0
          do iAtom = 1, nToGrow
!            if (nIntraNonBond(nType) .gt. 0) then
!              call Rosen_IntraNonBond_Atom_New(newMol%x(:), newMol%y(:), newMol%z(:), nType, atmGrow(iAtom), &
!                                               trialPos_Branched(iRosen,iAtom), regrown, E_Trial(iRosen), overlap(iRosen))
!            endif
            call Rosen_Atom_Intra_New(nType, atmGrow(iAtom), trialPos_Branched(iRosen,iAtom), &
                                                  regrown,  E_Trial_i, overlap(iRosen))
            E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_i
          enddo
        enddo
        E_Max = 0d0
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = wBending(iRosen) * wTorsion(iRosen) * exp(-beta*(E_Trial(iRosen)-E_Max))         
        enddo
        rosenNorm = sum(ProbRosen)
        ranNum = grnd() * rosenNorm
        sumInt = ProbRosen(1)
        nSel = 1
        do while(sumInt .lt. ranNum)
          nSel = nSel + 1
          sumInt = sumInt + ProbRosen(nSel)
        enddo
        if(overlap(nSel) .eqv. .true.) then
          rejMove = .true.
          return
        endif
        rosenRatio = rosenRatio*ProbRosen(nSel)/rosenNorm
        wForward = wForward * rosenNorm
        do iAtom = 1, nToGrow
          regrown(atmGrow(iAtom)) = .true.
          totalRegrown = totalRegrown + 1
          newMol%x(atmGrow(iAtom)) = trialPos_Branched(nSel,iAtom)%x 
          newMol%y(atmGrow(iAtom)) = trialPos_Branched(nSel,iAtom)%y 
          newMol%z(atmGrow(iAtom)) = trialPos_Branched(nSel,iAtom)%z
        enddo
      enddo
	  
      end subroutine
!=======================================================================
      subroutine BranchedMol_RosenConfigGen_GasPhase_Reverse(nType, wReverse)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions_LJ_Q
      use EnergyPressurePointers, only: Rosen_Atom_Intra_Gas_Old
      use CBMC_Variables
      use CBMC_Utility
      use CoordinateFunctions
      implicit none
	  
      integer, intent(in) :: nType
      real(dp), intent(out):: wReverse
	  
      integer :: iAtom, jAtom, iRosen, iGrow
      integer :: atmCur, atmPrev, atmGrow(1:3), atmTor(1:3), nToGrow, nTor
      integer :: bondType, bendTypes(1:6), torsTypes(1:3,1:3)
      logical :: regrown(1:maxAtoms)
	  
      real(dp) :: E_Trial(1:maxRosenTrial), E_Trial_i, wBending(1:maxRosenTrial), wTorsion(1:maxRosenTrial), E_Max
      real(dp) :: k_bonds(1:3), r_eqs(1:3), rGrow(1:3), Prob, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: bendingAngles(1:3), dihedralAngles(1:2), r_sigmas(1:3), rmax_sqs(1:3)
      type(SimpleAtomCoords) :: trialPos_Branched(1:3)
      type(SimpleAtomCoords) :: vPrev, vTor(1:3), vGrow(1:3)
	  
      ProbRosen = 0d0      
      E_Trial = 0d0      
      regrown = .false.
	  
!      Since the growth of first two atoms does not include nonBonded intera-molecular energies, we can write

      wReverse = dble(nRosenTrials(nType)) * dble(nRosenTrials(nType))
      regrown(1) = .true.
      regrown(2) = .true.
	  
      do iGrow = 1, SwapGrowOrder(nType)%GrowthSteps
        atmGrow = 0
        atmTor = 0
        nToGrow = SwapGrowOrder(nType)%GrowNum(iGrow)
        nTor = SwapGrowOrder(nType)%TorNum(iGrow)
        atmPrev = SwapGrowOrder(nType)%GrowPrev(iGrow)
        atmCur = SwapGrowOrder(nType)%GrowFrom(iGrow)
        vPrev%x = gasConfig(nType)%x(atmPrev) - gasConfig(nType)%x(atmCur)
        vPrev%y = gasConfig(nType)%y(atmPrev) - gasConfig(nType)%y(atmCur)
        vPrev%z = gasConfig(nType)%z(atmPrev) - gasConfig(nType)%z(atmCur)
        do iAtom = 1, nToGrow
          atmGrow(iAtom) = SwapGrowOrder(nType)%GrowList(iGrow,iAtom)
          vGrow(iAtom)%x = gasConfig(nType)%x(atmGrow(iAtom)) - gasConfig(nType)%x(atmCur)
          vGrow(iAtom)%y = gasConfig(nType)%y(atmGrow(iAtom)) - gasConfig(nType)%y(atmCur)
          vGrow(iAtom)%z = gasConfig(nType)%z(atmGrow(iAtom)) - gasConfig(nType)%z(atmCur)
          call FindBond(nType, atmCur, atmGrow(iAtom), bondType)
          k_bonds(iAtom) = bondData(bondType)%k_eq
          r_eqs(iAtom) = bondData(bondType)%r_eq
          r_sigmas(iAtom) = bondData(bondType)%r_sigma
          rmax_sqs(iAtom) = bondData(bondType)%rmax_sq
          call FindAngle(nType, atmPrev, atmCur, atmGrow(iAtom), bendTypes(iAtom))
          do jAtom = 1, nTor
            if (iAtom .eq. 1) then
              atmTor(jAtom) = SwapGrowOrder(nType)%TorList(iGrow,jAtom)
              vTor(jAtom)%x = gasConfig(nType)%x(atmTor(jAtom)) - gasConfig(nType)%x(atmCur)
              vTor(jAtom)%y = gasConfig(nType)%y(atmTor(jAtom)) - gasConfig(nType)%y(atmCur)
              vTor(jAtom)%z = gasConfig(nType)%z(atmTor(jAtom)) - gasConfig(nType)%z(atmCur)
            endif
            call FindTorsion(nType, atmTor(jAtom), atmPrev, atmCur, atmGrow(iAtom), torsTypes(iAtom,jAtom))
          enddo
        enddo
        E_Trial(1) = 0d0
        do iAtom = 1, nToGrow
          trialPos_Branched(iAtom)%x = gasConfig(nType)%x(atmGrow(iAtom))
          trialPos_Branched(iAtom)%y = gasConfig(nType)%y(atmGrow(iAtom))
          trialPos_Branched(iAtom)%z = gasConfig(nType)%z(atmGrow(iAtom))
!          if (nIntraNonBond(nType) .gt. 0) then
!            call Rosen_IntraNonBond_Atom_Old(gasConfig(nType)%x(:), gasConfig(nType)%y(:), &
!                       gasConfig(nType)%z(:), nType, atmGrow(iAtom), regrown, E_Trial(1))
!          endif
          call Rosen_Atom_Intra_Gas_Old(nType, atmGrow(iAtom), trialPos_Branched(iAtom), regrown, E_Trial_i)
          E_Trial(1) = E_Trial(1) + E_Trial_i
        enddo
        select case(nToGrow)
        case(1)
          wBending(1) = 1d0
        case(2)
          call FindAngle(nType, atmGrow(1), atmCur, atmGrow(2), bendTypes(3))
          wBending(1) = 0d0
          call GenerateTwoBranches_Reverse(vPrev, vGrow(1), vGrow(2), bendTypes(1), bendTypes(2), bendTypes(3), wBending(1))
        case(3)
          call FindAngle(nType, atmGrow(1), atmCur, atmGrow(2), bendTypes(4))
          call FindAngle(nType, atmGrow(2), atmCur, atmGrow(3), bendTypes(5))
          call FindAngle(nType, atmGrow(3), atmCur, atmGrow(1), bendTypes(6))
          wBending(1) = 0d0
          call GenerateThreeBranches_Reverse(vPrev, vGrow(1) ,vGrow(2) ,vGrow(3) , bendTypes(1), bendTypes(2), bendTypes(3), &
                                             bendTypes(4), bendTypes(5), bendTypes(6), wBending(1))
        case default
          write(*,*) "Error! Atom", atmCur, "in molecule type", nType, "has too many branches for a simple regrowth"
          stop
        end select
        call GenerateRotation_Reverse(nTor, nToGrow, vTor, vPrev, torsTypes, vGrow, wTorsion(1))
        do iRosen = 2, nRosenTrials(nType)
          bendingAngles = 0d0
          dihedralAngles = 0d0
          do iAtom = 1, nToGrow
            call GenerateGaussianBondLength(rGrow(iAtom), k_bonds(iAtom), r_eqs(iAtom), &
                                    r_sigmas(iAtom), rmax_sqs(iAtom), Prob)
          enddo
          select case(nToGrow)
          case(1)
            call GenerateBendAngle(bendingAngles(1), bendTypes(1), Prob)
            wBending(iRosen) = 1d0
          case(2)
            wBending(iRosen) = 0d0
            call GenerateTwoBranches(bendingAngles(1), bendingAngles(2), dihedralAngles(1), & 
                                     bendTypes(1), bendTypes(2), bendTypes(3), wBending(iRosen))
          case(3)
            wBending(iRosen) = 0d0
            call GenerateThreeBranches(bendingAngles(1), bendingAngles(1), bendingAngles(3), &
                                       dihedralAngles(1), dihedralAngles(2), bendTypes(1), bendTypes(2), bendTypes(3), &
                                       bendTypes(4), bendTypes(5), bendTypes(6), wBending(iRosen))
          case default
            write(*,*) "Error! Atom", atmCur, "in molecule type", nType, "has too many branches for a simple regrowth"
            stop
          end select
          call GenerateRotation(nTor, nToGrow, vTor, vPrev, rGrow, bendingAngles, dihedralAngles, &
                                torsTypes, vGrow, wTorsion(iRosen))
          E_Trial(iRosen) = 0d0
          do iAtom = 1, nToGrow
            trialPos_Branched(iAtom)%x = vGrow(iAtom)%x + gasConfig(nType)%x(atmCur)
            trialPos_Branched(iAtom)%y = vGrow(iAtom)%y + gasConfig(nType)%y(atmCur)
            trialPos_Branched(iAtom)%z = vGrow(iAtom)%z + gasConfig(nType)%z(atmCur)
            call Rosen_Atom_Intra_Gas_Old(nType, atmGrow(iAtom), trialPos_Branched(iAtom), &
                regrown, E_Trial_i)
            E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_i
          enddo

!          do iAtom = 1, nToGrow
!            if (nIntraNonBond(nType) .gt. 0) then
!              call Rosen_IntraNonBond_Atom_New(gasConfig(nType)%x(:), gasConfig(nType)%y(:), &
!                         gasConfig(nType)%z(:), nType, atmGrow(iAtom), trialPos_Branched(iAtom), regrown, E_Trial(iRosen), overlap)
!            endif
!          enddo
        enddo
        E_Max = 0d0
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = wBending(iRosen) * wTorsion(iRosen) * exp(-beta*(E_Trial(iRosen)-E_Max))         
        enddo
        rosenNorm = sum(ProbRosen)
        wReverse = wReverse*rosenNorm
        do iAtom = 1, nToGrow
          regrown(atmGrow(iAtom)) = .true.
        enddo
      enddo
	  
      end subroutine
!=======================================================================
      end module
