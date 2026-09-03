!=======================================================================
      subroutine Simple_Partial_ConfigGen(nType, nIndx, nMol, rosenRatio_New, rosenRatio_Old, rejMove)
      use SimParameters
      use Coords
      use CoordinateFunctions
      use ForceField
      use Constants
!      use Rosenbluth_Functions_LJ_Q
      use EnergyPressurePointers, only: Rosen_Mol_New, Rosen_Mol_Old
      use CBMC_Variables
      use CBMC_Utility
      implicit none

      integer, intent(in) :: nType, nMol, nIndx
      logical, intent(out) :: rejMove      
      real(dp), intent(out):: rosenRatio_New, rosenRatio_Old
      
      logical :: isIncluded(1:maxMol)
      logical :: overlap(1:maxRosenTrial)
      integer :: i, iRosen, nSel
      integer :: Atm1, Atm2, Atm3, Atm4, Atm5
      integer :: bondType
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
      type(SimpleAtomCoords) :: v1, v2, v3, v4, v5
      real(dp) :: wBending(1:maxRosenTrial)
      
      newMol%molType = nType      
      call Rosen_CreateSubset(nIndx, isIncluded)
      isIncluded(nIndx) = .false.
      E_Trial = 0d0
      rejMove = .false.      
      do iRosen = 1, nRosenTrials(nType)
     
         !Initialize the first atom coordinates to 0      
        rosenTrial(iRosen)%x = 0d0
        rosenTrial(iRosen)%y = 0d0
        rosenTrial(iRosen)%z = 0d0

         !Depending the number of atoms in the molecule, choose an appropriate algorithm
         !and regrow the atoms starting from atom 1. 
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
          rosenTrial(iRosen)%x(Atm2) = r1 * dx
          rosenTrial(iRosen)%y(Atm2) = r1 * dy
          rosenTrial(iRosen)%z(Atm2) = r1 * dz
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
        
        x1 = molArray(nType)%mol(nMol)%x(1) - rosenTrial(iRosen)%x(1)
        y1 = molArray(nType)%mol(nMol)%y(1) - rosenTrial(iRosen)%y(1)
        z1 = molArray(nType)%mol(nMol)%z(1) - rosenTrial(iRosen)%z(1)
        
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

      if(all(ProbRosen .le. 0d0)) then
        rejMove = .true.
        return
      endif
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
      rosenRatio_New = ProbRosen(nSel)/rosenNorm

!      Update the coordinates      
      newMol%x(1:nAtoms(nType)) = rosenTrial(nSel)%x(1:nAtoms(nType))
      newMol%y(1:nAtoms(nType)) = rosenTrial(nSel)%y(1:nAtoms(nType))
      newMol%z(1:nAtoms(nType)) = rosenTrial(nSel)%z(1:nAtoms(nType))
      


!      To simulate the reverse move, calculate the 
!      call Rosen_BoltzWeight_Molecule_Old(molArray(nType)%mol(nMol)%x(:), molArray(nType)%mol(nMol)%y(:), &
!                                 molArray(nType)%mol(nMol)%z(:), nType, isIncluded, E_Trial(1))   
      wBending = 1d0
      call Rosen_Mol_Old(molArray(nType)%mol(nMol)%x(:), molArray(nType)%mol(nMol)%y(:), &
                                 molArray(nType)%mol(nMol)%z(:), nType, isIncluded, E_Trial(1))  
      if (nAtoms(nType) .eq. 4) then
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)
          Atm3 = pathArray(nType)%path(2, 1)
          Atm4 = pathArray(nType)%path(3, 1)
          call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
          call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
          call FindAngle(nType, Atm3, Atm2, Atm4, bendType3)
          v1%x = MolArray(nType)%mol(nMol)%x(Atm1) - MolArray(nType)%mol(nMol)%x(Atm2)
          v1%y = MolArray(nType)%mol(nMol)%y(Atm1) - MolArray(nType)%mol(nMol)%y(Atm2)
          v1%z = MolArray(nType)%mol(nMol)%z(Atm1) - MolArray(nType)%mol(nMol)%z(Atm2)
          v3%x = MolArray(nType)%mol(nMol)%x(Atm3) - MolArray(nType)%mol(nMol)%x(Atm2)
          v3%y = MolArray(nType)%mol(nMol)%y(Atm3) - MolArray(nType)%mol(nMol)%y(Atm2)
          v3%z = MolArray(nType)%mol(nMol)%z(Atm3) - MolArray(nType)%mol(nMol)%z(Atm2)
          v4%x = MolArray(nType)%mol(nMol)%x(Atm4) - MolArray(nType)%mol(nMol)%x(Atm2)
          v4%y = MolArray(nType)%mol(nMol)%y(Atm4) - MolArray(nType)%mol(nMol)%y(Atm2)
          v4%z = MolArray(nType)%mol(nMol)%z(Atm4) - MolArray(nType)%mol(nMol)%z(Atm2)
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
          v1%x = MolArray(nType)%mol(nMol)%x(Atm1) - MolArray(nType)%mol(nMol)%x(Atm2)
          v1%y = MolArray(nType)%mol(nMol)%y(Atm1) - MolArray(nType)%mol(nMol)%y(Atm2)
          v1%z = MolArray(nType)%mol(nMol)%z(Atm1) - MolArray(nType)%mol(nMol)%z(Atm2)
          v3%x = MolArray(nType)%mol(nMol)%x(Atm3) - MolArray(nType)%mol(nMol)%x(Atm2)
          v3%y = MolArray(nType)%mol(nMol)%y(Atm3) - MolArray(nType)%mol(nMol)%y(Atm2)
          v3%z = MolArray(nType)%mol(nMol)%z(Atm3) - MolArray(nType)%mol(nMol)%z(Atm2)
          v4%x = MolArray(nType)%mol(nMol)%x(Atm4) - MolArray(nType)%mol(nMol)%x(Atm2)
          v4%y = MolArray(nType)%mol(nMol)%y(Atm4) - MolArray(nType)%mol(nMol)%y(Atm2)
          v4%z = MolArray(nType)%mol(nMol)%z(Atm4) - MolArray(nType)%mol(nMol)%z(Atm2)
          v5%x = MolArray(nType)%mol(nMol)%x(Atm5) - MolArray(nType)%mol(nMol)%x(Atm2)
          v5%y = MolArray(nType)%mol(nMol)%y(Atm5) - MolArray(nType)%mol(nMol)%y(Atm2)
          v5%z = MolArray(nType)%mol(nMol)%z(Atm5) - MolArray(nType)%mol(nMol)%z(Atm2)
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
         stop "Error! Molecule has too many atoms for a simple regrowth"
        end select
        
        x1 = molArray(nType)%mol(nMol)%x(1) - rosenTrial(iRosen)%x(1)
        y1 = molArray(nType)%mol(nMol)%y(1) - rosenTrial(iRosen)%y(1)
        z1 = molArray(nType)%mol(nMol)%z(1) - rosenTrial(iRosen)%z(1)
        
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
      rosenNorm = sum(ProbRosen)
      rosenRatio_Old = ProbRosen(1)/rosenNorm

      
      end subroutine

!=======================================================================
      subroutine StraightChain_Partial_ConfigGen(nType, nMol, regrownIn, regrowDirection, startAtom, rosenRatio, rejMove)
      use SimParameters
      use Coords
      use CoordinateFunctions
      use ForceField
      use Constants
!      use Rosenbluth_Functions_LJ_Q
      use EnergyPressurePointers, only: Rosen_Atom_New, Rosen_Atom_Intra_New
      use CBMC_Variables
      use CBMC_Utility
      implicit none

      integer, intent(in) :: nType, nMol
      logical, intent(in) :: regrowDirection
      integer, intent(in) :: startAtom
      real(dp), intent(out):: rosenRatio
      logical, intent(out) :: rejMove
      logical, intent(in) :: regrownIn(1:maxAtoms)

      integer :: i, iRosen, nIndx
      integer :: Atm1, Atm2, Atm3, Atm4
      integer :: Incrmt
      integer :: bondType, bendType, torsType
      integer :: totalRegrown, curPos
      logical :: regrown(1:maxAtoms)
      logical :: isIncluded(1:maxMol)
      logical :: overlap(1:maxRosenTrial)

      real(dp) :: E_Trial(1:maxRosenTrial), E_Intra
      real(dp) :: k_bond, r_eq, r, Prob, r_sigma, rmax_sq
      real(dp) :: bend_angle, tors_angle
      real(dp) :: dx, dy, dz 
      type(SimpleAtomCoords) :: trialPos(1:maxRosenTrial)
      type(SimpleAtomCoords) :: v1, v2, v3
      
      do i = 1, maxAtoms
        regrown(i) = regrownIn(i)
      enddo

      if(regrowDirection) then
        Incrmt =  1
      else
        Incrmt = -1
      endif  
 
      nIndx = molArray(nType)%mol(nMol)%indx
      call Rosen_CreateSubset(nIndx, isIncluded)
      isIncluded(nIndx) = .false. 
      newMol%molType = nType 
      newMol%x = 0d0
      newMol%y = 0d0
      newMol%z = 0d0
      rejMove = .false.

!      For any atoms that are not chosen for regrowth, assign 
      totalRegrown = 0
      do i = 1, nAtoms(nType)
        if(regrown(i)) then
         newMol%x(i) = molArray(nType)%mol(nMol)%x(i)
         newMol%y(i) = molArray(nType)%mol(nMol)%y(i)
         newMol%z(i) = molArray(nType)%mol(nMol)%z(i)
         totalRegrown = totalRegrown + 1
        endif
      enddo
       
!      Having inserted the first atom, when we go to insert the second atom we must account for the fact that
!      the second atom must be inserted at the correct bond distance from the first atom.
      do i = 1, patharray(nType)%pathmax(1)
        if(patharray(nType)%path(1,i) .eq. startAtom) then
          curPos = i
        endif
      enddo

      rosenRatio = 1d0      
      if(totalRegrown .eq. 1) then
        if(curPos .eq. 1) then
          curPos = curPos + 1
        endif
        if(curPos .eq. patharray(nType)%pathmax(1)) then
          curPos = curPos - 1
        endif
        Atm2 = patharray(nType)%path(1,curPos)
        Atm1 = patharray(nType)%path(1,curPos-Incrmt)
!        call FindNextBondFromPath(nType, regrown, 1, Atm2, Atm1)
        call FindBond(nType,Atm1, Atm2, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        r_sigma = bondData(bondType)%r_sigma
        rmax_sq = bondData(bondType)%rmax_sq
        overlap = .false.
        do iRosen = 1, nRosenTrials(nType)
          call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          trialPos(iRosen)%x = r*dx + newMol%x(Atm1) 
          trialPos(iRosen)%y = r*dy + newMol%y(Atm1)
          trialPos(iRosen)%z = r*dz + newMol%z(Atm1)
          call Rosen_Atom_New(nType, Atm2, trialPos(iRosen), isIncluded,  E_Trial(iRosen), overlap(iRosen))
        enddo
        call ChooseRosenTrial(Atm2, nType, trialPos, regrown, E_Trial, overlap, rosenRatio, rejMove)  
        if(rejMove) then
          return
        endif
        regrown(Atm2) = .true.
        curPos = curPos + Incrmt
        totalRegrown = totalRegrown + 1
      endif

!       For the third atom we must begin to include the bending angle in choosing its trial positions. The first thing
!       we must consider is if the second regrown atom was a terminal atom or a linker atom. If it was a terminal atom 
!       then the bending angle must use the 1st atom as the central atom for the angle generation.  However if it was a link in the chain then the 2nd 
!       atom regrown should be used as the central atom.
      if(totalRegrown .eq. 2) then
        Atm3 = patharray(nType)%path(1,curPos)
        Atm2 = patharray(nType)%path(1,curPos - 1*Incrmt)
        Atm1 = patharray(nType)%path(1,curPos - 2*Incrmt)
        v1%x = newMol%x(Atm1) - newMol%x(Atm2)
        v1%y = newMol%y(Atm1) - newMol%y(Atm2)
        v1%z = newMol%z(Atm1) - newMol%z(Atm2)
        call FindBond(nType, Atm2, Atm3, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        r_sigma = bondData(bondType)%r_sigma
        rmax_sq = bondData(bondType)%rmax_sq
        call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
!        k_bend = bendData(bendType)%k_eq
!        ang_eq = bendData(bendType)%ang_eq
        overlap = .false.
        do iRosen = 1, nRosenTrials(nType)
          call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
!          call GenerateBendAngle(bend_angle, k_bend, ang_eq, Prob)
          call GenerateBendAngle(bend_angle, bendType, Prob)
          call Generate_UnitCone(v1, r, bend_angle, v2)
          trialPos(iRosen)%x = v2%x + newMol%x(Atm2) 
          trialPos(iRosen)%y = v2%y + newMol%y(Atm2)
          trialPos(iRosen)%z = v2%z + newMol%z(Atm2)
          call Rosen_Atom_New(nType, Atm3, trialPos(iRosen), isIncluded, E_Trial(iRosen), overlap(iRosen))
          if (.not. overlap(iRosen)) then
            call Rosen_Atom_Intra_New(nType, Atm3, trialPos(iRosen), regrown, E_Intra, overlap(iRosen))
            E_Trial(iRosen) = E_Trial(iRosen) + E_Intra
          endif
        enddo
        call ChooseRosenTrial(Atm3, nType, trialPos, regrown, E_Trial, overlap, rosenRatio, rejMove)
        if(rejMove) then
          return
        endif
        regrown(Atm3) = .true.
        curPos = curPos + Incrmt   
        totalRegrown = totalRegrown + 1
      endif


!      Now that three atoms have been regrown, all remaining atoms must take the torsional angles into account.
      do while(any(regrown .eqv. .false.))
        Atm4 = patharray(nType)%path(1,curPos           )
        Atm3 = patharray(nType)%path(1,curPos - 1*Incrmt)
        Atm2 = patharray(nType)%path(1,curPos - 2*Incrmt)
        Atm1 = patharray(nType)%path(1,curPos - 3*Incrmt)
!        call FindAtomsFromPath(nType, regrown, 1, Atm4, Atm1, Atm2, Atm3)
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
            E_Trial(iRosen) = E_Trial(iRosen) + E_Intra
          endif
        enddo
        call ChooseRosenTrial(Atm4, nType, trialPos, regrown, E_Trial, overlap, rosenRatio, rejMove)
        if(rejMove) then
          return
        endif
        regrown(Atm4) = .true.
        curPos = curPos + Incrmt   
        totalRegrown = totalRegrown + 1
      enddo

      end subroutine
!=======================================================================
      subroutine StraightChain_Partial_ConfigGen_Reverse(nType, nMol, regrownIn, regrowDirection, startAtom, rosenRatio)
      use SimParameters
      use Coords
      use CoordinateFunctions
      use ForceField
      use Constants
!      use Rosenbluth_Functions_LJ_Q
      use EnergyPressurePointers, only: Rosen_Atom_Old, Rosen_Atom_New, Rosen_Atom_Intra_Old
      use CBMC_Variables
      use CBMC_Utility
      implicit none

      integer, intent(in) :: nType, nMol
      logical, intent(in) :: regrowDirection
      integer, intent(in) :: startAtom
      real(dp), intent(out):: rosenRatio
      logical, intent(in) :: regrownIn(1:maxAtoms)


      integer :: i, iRosen, nIndx
      integer :: Atm1, Atm2, Atm3, Atm4
      integer :: Incrmt
      integer :: bondType, bendType, torsType
      integer :: totalRegrown, curPos
      logical :: regrown(1:maxAtoms)
      logical :: isIncluded(1:maxMol)
      logical :: overlap

      real(dp) :: E_Trial(1:maxRosenTrial), E_Intra
      real(dp) :: E_Min, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: k_bond, r_eq, r, Prob, r_sigma, rmax_sq
      real(dp) :: bend_angle, tors_angle
      real(dp) :: dx, dy, dz 
      type(SimpleAtomCoords) :: trialPos
      type(SimpleAtomCoords) :: v1, v2, v3
      
   
      regrown(1:maxAtoms) = regrownIn(1:maxAtoms)
      if(regrowDirection) then
        Incrmt =  1
      else
        Incrmt = -1
      endif  
 
      nIndx = molArray(nType)%mol(nMol)%indx
      call Rosen_CreateSubset(nIndx, isIncluded)
      isIncluded(nIndx) = .false. 
!      For any atoms that are not chosen for regrowth, assign 
      totalRegrown = 0
      do i = 1, nAtoms(nType)
        if(regrown(i)) then
          totalRegrown = totalRegrown + 1
        endif
      enddo
       
!      Having inserted the first atom, when we go to insert the second atom we must account for the fact that
!      the second atom must be inserted at the correct bond distance from the first atom.
      do i = 1, patharray(nType)%pathmax(1)
        if(patharray(nType)%path(1,i) .eq. startAtom) then
          curPos = i
        endif
      enddo

      rosenRatio = 1d0
      
      if(totalRegrown .eq. 1) then
        if(curPos .eq. 1) then
          curPos = curPos + 1
        endif
        if(curPos .eq. patharray(nType)%pathmax(1)) then
          curPos = curPos - 1
        endif
        Atm2 = patharray(nType)%path(1,curPos)
        Atm1 = patharray(nType)%path(1,curPos-Incrmt)
!        call FindNextBondFromPath(nType, regrown, 1, Atm2, Atm1)
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
        E_Min = minval(E_Trial)
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Min))         
        enddo
        rosenNorm = sum(ProbRosen)
        rosenRatio = rosenRatio*ProbRosen(1)/rosenNorm
        regrown(Atm2) = .true.
        curPos = curPos + Incrmt
        totalRegrown = totalRegrown + 1
      endif

!       For the third atom we must begin to include the bending angle in choosing its trial positions. The first thing
!       we must consider is if the second regrown atom was a terminal atom or a linker atom. If it was a terminal atom 
!       then the bending angle must use the 1st atom as the central atom for the angle generation.  However if it was a link in the chain then the 2nd 
!       atom regrown should be used as the central atom.
      if(totalRegrown .eq. 2) then
        Atm3 = patharray(nType)%path(1,curPos)
        Atm2 = patharray(nType)%path(1,curPos - 1*Incrmt)
        Atm1 = patharray(nType)%path(1,curPos - 2*Incrmt)
        v1%x = molArray(nType)%mol(nMol)%x(Atm1) - molArray(nType)%mol(nMol)%x(Atm2)
        v1%y = molArray(nType)%mol(nMol)%y(Atm1) - molArray(nType)%mol(nMol)%y(Atm2)
        v1%z = molArray(nType)%mol(nMol)%z(Atm1) - molArray(nType)%mol(nMol)%z(Atm2)
        call FindBond(nType, Atm2, Atm3, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        r_sigma = bondData(bondType)%r_sigma
        rmax_sq = bondData(bondType)%rmax_sq
        call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
!        k_bend = bendData(bendType)%k_eq
!        ang_eq = bendData(bendType)%ang_eq
        call Rosen_Atom_Old(nType, nMol, Atm3, isIncluded, E_Trial(1))
        trialPos%x = molArray(nType)%mol(nMol)%x(Atm3) 
        trialPos%y = molArray(nType)%mol(nMol)%y(Atm3)
        trialPos%z = molArray(nType)%mol(nMol)%z(Atm3)
        call Rosen_Atom_Intra_Old(nType, nMol, Atm3, trialPos, regrown, E_Intra)
        E_Trial(1) = E_Trial(1) + E_Intra
        do iRosen = 2, nRosenTrials(nType)
          call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
!          call GenerateBendAngle(bend_angle, k_bend, ang_eq, Prob)
          call GenerateBendAngle(bend_angle, bendType, Prob)
          call Generate_UnitCone(v1, r, bend_angle, v2)
          trialPos%x = v2%x + molArray(nType)%mol(nMol)%x(Atm2) 
          trialPos%y = v2%y + molArray(nType)%mol(nMol)%y(Atm2)
          trialPos%z = v2%z + molArray(nType)%mol(nMol)%z(Atm2)
          call Rosen_Atom_New(nType, Atm3, trialPos, isIncluded,  E_Trial(iRosen), overlap)
          call Rosen_Atom_Intra_Old(nType, nMol, Atm3, trialPos, regrown, E_Intra)
          E_Trial(iRosen) = E_Trial(iRosen) + E_Intra
        enddo
        E_Min = minval(E_Trial)
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Min))         
        enddo
        rosenNorm = sum(ProbRosen)
        rosenRatio = rosenRatio*ProbRosen(1)/rosenNorm
        regrown(Atm3) = .true.
        curPos = curPos + Incrmt   
        totalRegrown = totalRegrown + 1
      endif


!      Now that three atoms have been regrown, all remaining atoms must take the torsional angles into account.
      do while(any(regrown .eqv. .false.))
        Atm4 = patharray(nType)%path(1,curPos           )
        Atm3 = patharray(nType)%path(1,curPos - 1*Incrmt)
        Atm2 = patharray(nType)%path(1,curPos - 2*Incrmt)
        Atm1 = patharray(nType)%path(1,curPos - 3*Incrmt)
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
        E_Min = minval(E_Trial)
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Min))         
        enddo
        rosenNorm = sum(ProbRosen)
        rosenRatio = rosenRatio*ProbRosen(1)/rosenNorm
        regrown(Atm4) = .true.
        curPos = curPos + Incrmt   
        totalRegrown = totalRegrown + 1
      enddo

      end subroutine
!=======================================================================
      subroutine BranchedMol_Partial_ConfigGen(nType, nMol, regrownIn, nGrow,  &
                                               GrowFrom, GrowPrev, GrowNum, GrowList, TorNum, TorList, rosenRatio, rejMove)
      use SimParameters
      use Coords
      use CoordinateFunctions
      use ForceField
      use Constants
!      use Rosenbluth_Functions_LJ_Q
      use EnergyPressurePointers, only: Rosen_Atom_New, Rosen_Atom_Intra_New
      use CBMC_Variables
      use CBMC_Utility
      implicit none
	  
      integer, intent(in) :: nType, nMol, nGrow
      real(dp), intent(out):: rosenRatio
      logical, intent(out) :: rejMove
      logical, intent(in) :: regrownIn(1:maxAtoms)
      integer, intent(in) :: GrowFrom(1:maxAtoms), GrowPrev(1:maxAtoms), GrowNum(1:maxAtoms), TorNum(1:maxAtoms)
      integer, intent(in) :: GrowList(1:maxAtoms,1:maxBranches), TorList(1:maxAtoms,1:maxBranches)
	  
      integer :: iAtom, jAtom, iRosen, nIndx, iGrow, jGrow
      integer :: atmCur, atmPrev, atmGrow(1:3), atmTor(1:3)
      integer :: bondType, bendTypes(1:6), torsTypes(1:3,1:3)
      integer :: totalRegrown
      logical :: regrown(1:maxAtoms)
      logical :: isIncluded(1:maxMol)
      logical :: overlap(1:maxRosenTrial)
	  
      real(dp) :: E_Trial(1:maxRosenTrial), E_Trial_i, wBending(1:maxRosenTrial), wTorsion(1:maxRosenTrial)
      real(dp) :: k_bonds(1:3), r_eqs(1:3), r, rGrow(1:3), Prob, r_sigmas(1:3), rmax_sqs(1:3)
      real(dp) :: bendingAngles(1:3), dihedralAngles(1:2)
      real(dp) :: dx, dy, dz 
      type(SimpleAtomCoords) :: trialPos(1:maxRosenTrial), trialPos_Branched(1:maxRosenTrial,1:3)
      type(SimpleAtomCoords) :: vPrev, vTor(1:3), vGrow(1:3)
	  
      do iAtom = 1, maxAtoms
        regrown(iAtom) = regrownIn(iAtom)
      enddo
	  
      nIndx = molArray(nType)%mol(nMol)%indx
      call Rosen_CreateSubset(nIndx, isIncluded)
      isIncluded(nIndx) = .false.
      newMol%molType = nType 
      newMol%x = 0d0
      newMol%y = 0d0
      newMol%z = 0d0
      rejMove = .false.

!      For any atoms that are not chosen for regrowth, assign 
      totalRegrown = 0
      do iAtom = 1, nAtoms(nType)
        if(regrown(iAtom)) then
         newMol%x(iAtom) = molArray(nType)%mol(nMol)%x(iAtom)
         newMol%y(iAtom) = molArray(nType)%mol(nMol)%y(iAtom)
         newMol%z(iAtom) = molArray(nType)%mol(nMol)%z(iAtom)
         totalRegrown = totalRegrown + 1
        endif
      enddo
	  
      rosenRatio = 1d0 
      jGrow = 1
! In the case of regrowing the whole molecule, only one atom (one terminal) presents and the first 
! growing atom is bonded to this terminal
      if(totalRegrown .eq. 1) then
        atmCur = GrowFrom(jGrow)
        atmGrow(1) = GrowList(jGrow,1)
        call FindBond(nType, atmCur, atmGrow(1), bondType)
        k_bonds(1) = bondData(bondType)%k_eq
        r_eqs(1) = bondData(bondType)%r_eq
        r_sigmas(1) = bondData(bondType)%r_sigma
        rmax_sqs(1) = bondData(bondType)%rmax_sq
        overlap = .false.
        do iRosen = 1, nRosenTrials(nType)
          call GenerateGaussianBondLength(r, k_bonds(1), r_eqs(1), r_sigmas(1), rmax_sqs(1), Prob)
          call Generate_UnitSphere(dx, dy, dz)
          trialPos(iRosen)%x = r*dx + newMol%x(atmCur) 
          trialPos(iRosen)%y = r*dy + newMol%y(atmCur)
          trialPos(iRosen)%z = r*dz + newMol%z(atmCur)
          call Rosen_Atom_New(nType, atmGrow(1), trialPos(iRosen), isIncluded,  E_Trial(iRosen), overlap(iRosen))
        enddo
        call ChooseRosenTrial(atmGrow(1), nType, trialPos, regrown, E_Trial, overlap, rosenRatio, rejMove)  
        if(rejMove) then
          return
        endif
        regrown(atmGrow(1)) = .true.
        totalRegrown = totalRegrown + 1
        jGrow = jGrow + 1
      endif
	  
      do iGrow = jGrow, nGrow
        atmGrow = 0
        atmTor = 0
        atmPrev = GrowPrev(iGrow)
        atmCur = GrowFrom(iGrow)
        vPrev%x = newMol%x(atmPrev) - newMol%x(atmCur)
        vPrev%y = newMol%y(atmPrev) - newMol%y(atmCur)
        vPrev%z = newMol%z(atmPrev) - newMol%z(atmCur)
        do iAtom = 1, GrowNum(iGrow)
          atmGrow(iAtom) = GrowList(iGrow,iAtom)
          call FindBond(nType, atmCur, atmGrow(iAtom), bondType)
          k_bonds(iAtom) = bondData(bondType)%k_eq
          r_eqs(iAtom) = bondData(bondType)%r_eq
          r_sigmas(iAtom) = bondData(bondType)%r_sigma
          rmax_sqs(iAtom) = bondData(bondType)%rmax_sq
          call FindAngle(nType, atmPrev, atmCur, atmGrow(iAtom), bendTypes(iAtom))
          do jAtom = 1, TorNum(iGrow)
            if (iAtom .eq. 1) then
              atmTor(jAtom) = TorList(iGrow,jAtom)
              vTor(jAtom)%x = newMol%x(atmTor(jAtom)) - newMol%x(atmCur)
              vTor(jAtom)%y = newMol%y(atmTor(jAtom)) - newMol%y(atmCur)
              vTor(jAtom)%z = newMol%z(atmTor(jAtom)) - newMol%z(atmCur)
            endif
            call FindTorsion(nType, atmTor(jAtom), atmPrev, atmCur, atmGrow(iAtom), torsTypes(iAtom,jAtom))
          enddo
        enddo
        if (GrowNum(iGrow) .eq. 2) then
          call FindAngle(nType, atmGrow(1), atmCur, atmGrow(2), bendTypes(3))
        elseif (GrowNum(iGrow) .eq. 3) then
          call FindAngle(nType, atmGrow(1), atmCur, atmGrow(2), bendTypes(4))
          call FindAngle(nType, atmGrow(2), atmCur, atmGrow(3), bendTypes(5))
          call FindAngle(nType, atmGrow(3), atmCur, atmGrow(1), bendTypes(6))
        endif
        overlap = .false.
        do iRosen = 1, nRosenTrials(nType)
          bendingAngles = 0d0
          dihedralAngles = 0d0
          do iAtom = 1, GrowNum(iGrow)
            call GenerateGaussianBondLength(rGrow(iAtom), k_bonds(iAtom), r_eqs(iAtom), &
                                            r_sigmas(iAtom), rmax_sqs(iAtom), Prob)
          enddo
          select case(GrowNum(iGrow))
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
          call GenerateRotation(TorNum(iGrow), GrowNum(iGrow), vTor, vPrev, rGrow, bendingAngles, dihedralAngles, &
                                torsTypes, vGrow, wTorsion(iRosen))
          do iAtom = 1, GrowNum(iGrow)
            trialPos_Branched(iRosen,iAtom)%x = vGrow(iAtom)%x + newMol%x(atmCur)
            trialPos_Branched(iRosen,iAtom)%y = vGrow(iAtom)%y + newMol%y(atmCur)
            trialPos_Branched(iRosen,iAtom)%z = vGrow(iAtom)%z + newMol%z(atmCur)
          enddo
          E_Trial(iRosen) = 0d0
          do iAtom = 1, GrowNum(iGrow)
            call Rosen_Atom_New(nType, atmGrow(iAtom), trialPos_Branched(iRosen,iAtom), isIncluded,  E_Trial_i, &
                                            overlap(iRosen))
            E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_i
            call Rosen_Atom_Intra_New(nType, atmGrow(iAtom), trialPos_Branched(iRosen,iAtom), regrown, &
                                                  E_Trial_i, overlap(iRosen))
            E_Trial(iRosen) = E_Trial(iRosen) + E_Trial_i
!            if (nIntraNonBond(nType) .gt. 0) then
!              call Rosen_IntraNonBond_Atom_New(newMol%x(:), newMol%y(:), newMol%z(:), nType, atmGrow(iAtom), &
!                                               trialPos_Branched(iRosen,iAtom), regrown, E_Trial(iRosen), overlap(iRosen))
!            endif
          enddo
        enddo
        call ChooseRosenTrial_Branched(GrowNum(iGrow), atmGrow, nType, trialPos_Branched, regrown, E_Trial, wBending, &
                                       wTorsion, overlap, rosenRatio, rejMove)
        if(rejMove) then
          return
        endif
        do iAtom = 1, GrowNum(iGrow)
          totalRegrown = totalRegrown + 1
        enddo
      enddo
	  	  
      end subroutine
!=======================================================================
      subroutine BranchedMol_Partial_ConfigGen_Reverse(nType, nMol, regrownIn, nGrow, &
                                                       GrowFrom, GrowPrev, GrowNum, GrowList, TorNum, TorList, rosenRatio)
      use SimParameters
      use Coords
      use CoordinateFunctions
      use ForceField
      use Constants
!      use Rosenbluth_Functions_LJ_Q
      use EnergyPressurePointers, only: Rosen_Atom_New, Rosen_Atom_Old, Rosen_Atom_Intra_Old
      use CBMC_Variables
      use CBMC_Utility
      implicit none

      integer, intent(in) :: nType, nMol, nGrow
      real(dp), intent(out):: rosenRatio
      logical, intent(in) :: regrownIn(1:maxAtoms)
      integer, intent(in) :: GrowFrom(1:maxAtoms), GrowPrev(1:maxAtoms), GrowNum(1:maxAtoms), TorNum(1:maxAtoms)
      integer, intent(in) :: GrowList(1:maxAtoms,1:maxBranches), TorList(1:maxAtoms,1:maxBranches)
	  
      integer :: iAtom, jAtom, iRosen, nIndx, iGrow, jGrow
      integer :: atmCur, atmPrev, atmGrow(1:3), atmTor(1:3)
      integer :: bondType, bendTypes(1:6), torsTypes(1:3,1:3)
      integer :: totalRegrown
      logical :: regrown(1:maxAtoms)
      logical :: isIncluded(1:maxMol)
      logical :: overlap
	  
      real(dp) :: E_Trial(1:maxRosenTrial), E_Trial_i, wBending(1:maxRosenTrial), wTorsion(1:maxRosenTrial)
      real(dp) :: ProbRosen(1:maxRosenTrial), E_Min, rosenNorm
      real(dp) :: k_bonds(1:3), r_eqs(1:3), r, rGrow(1:3), Prob, r_sigmas(1:3), rmax_sqs(1:3)
      real(dp) :: bendingAngles(1:3), dihedralAngles(1:2)
      real(dp) :: dx, dy, dz 
      type(SimpleAtomCoords) :: trialPos, trialPos_Branched(1:3)
      type(SimpleAtomCoords) :: vPrev, vTor(1:3), vGrow(1:3)
	  
      do iAtom = 1, maxAtoms
        regrown(iAtom) = regrownIn(iAtom)
      enddo
	  
      nIndx = molArray(nType)%mol(nMol)%indx
      call Rosen_CreateSubset(nIndx, isIncluded)
      isIncluded(nIndx) = .false.

!      For any atoms that are not chosen for regrowth, assign 
      totalRegrown = 0
      do iAtom = 1, nAtoms(nType)
        if(regrown(iAtom)) then
         totalRegrown = totalRegrown + 1
        endif
      enddo
	  
      rosenRatio = 1d0 
      jGrow = 1
! In the case of regrowing the whole molecule, only one atom (one terminal) presents and the first 
! growing atom is bonded to this terminal
      if(totalRegrown .eq. 1) then
        atmCur = GrowFrom(jGrow)
        atmGrow(1) = GrowList(jGrow,1)
        call FindBond(nType, atmCur, atmGrow(1), bondType)
        k_bonds(1) = bondData(bondType)%k_eq
        r_eqs(1) = bondData(bondType)%r_eq
        r_sigmas(1) = bondData(bondType)%r_sigma
        rmax_sqs(1) = bondData(bondType)%rmax_sq
        call Rosen_Atom_Old(nType, nMol, atmGrow(1), isIncluded,  E_Trial(1))
        do iRosen = 2, nRosenTrials(nType)
          call GenerateGaussianBondLength(r, k_bonds(1), r_eqs(1), r_sigmas(1), rmax_sqs(1), Prob)
          call Generate_UnitSphere(dx, dy, dz)
          trialPos%x = r*dx + molArray(nType)%mol(nMol)%x(atmCur) 
          trialPos%y = r*dy + molArray(nType)%mol(nMol)%y(atmCur)
          trialPos%z = r*dz + molArray(nType)%mol(nMol)%z(atmCur)
          call Rosen_Atom_New(nType, atmGrow(1), trialPos, isIncluded,  E_Trial(iRosen), overlap)
        enddo
        E_Min = minval(E_Trial)
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Min))         
        enddo
        rosenNorm = sum(ProbRosen)
        rosenRatio = rosenRatio*ProbRosen(1)/rosenNorm
        regrown(atmGrow(1)) = .true.
        totalRegrown = totalRegrown + 1
        jGrow = jGrow + 1
      endif
	  
      do iGrow = jGrow, nGrow
        atmGrow = 0
        atmTor = 0
        atmPrev = GrowPrev(iGrow)
        atmCur = GrowFrom(iGrow)
        vPrev%x = molArray(nType)%mol(nMol)%x(atmPrev) - molArray(nType)%mol(nMol)%x(atmCur)
        vPrev%y = molArray(nType)%mol(nMol)%y(atmPrev) - molArray(nType)%mol(nMol)%y(atmCur)
        vPrev%z = molArray(nType)%mol(nMol)%z(atmPrev) - molArray(nType)%mol(nMol)%z(atmCur)
        do iAtom = 1, GrowNum(iGrow)
          atmGrow(iAtom) = GrowList(iGrow,iAtom)
          vGrow(iAtom)%x = molArray(nType)%mol(nMol)%x(atmGrow(iAtom)) - molArray(nType)%mol(nMol)%x(atmCur)
          vGrow(iAtom)%y = molArray(nType)%mol(nMol)%y(atmGrow(iAtom)) - molArray(nType)%mol(nMol)%y(atmCur)
          vGrow(iAtom)%z = molArray(nType)%mol(nMol)%z(atmGrow(iAtom)) - molArray(nType)%mol(nMol)%z(atmCur)
          call FindBond(nType, atmCur, atmGrow(iAtom), bondType)
          k_bonds(iAtom) = bondData(bondType)%k_eq
          r_eqs(iAtom) = bondData(bondType)%r_eq
          r_sigmas(iAtom) = bondData(bondType)%r_sigma
          rmax_sqs(iAtom) = bondData(bondType)%rmax_sq
          call FindAngle(nType, atmPrev, atmCur, atmGrow(iAtom), bendTypes(iAtom))
          do jAtom = 1, TorNum(iGrow)
            if (iAtom .eq. 1) then
              atmTor(jAtom) = TorList(iGrow,jAtom)
              vTor(jAtom)%x = molArray(nType)%mol(nMol)%x(atmTor(jAtom)) - molArray(nType)%mol(nMol)%x(atmCur)
              vTor(jAtom)%y = molArray(nType)%mol(nMol)%y(atmTor(jAtom)) - molArray(nType)%mol(nMol)%y(atmCur)
              vTor(jAtom)%z = molArray(nType)%mol(nMol)%z(atmTor(jAtom)) - molArray(nType)%mol(nMol)%z(atmCur)
            endif
            call FindTorsion(nType, atmTor(jAtom), atmPrev, atmCur, atmGrow(iAtom), torsTypes(iAtom,jAtom))
          enddo
        enddo
        E_Trial(1) = 0d0
        do iAtom = 1, GrowNum(iGrow)
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
        select case(GrowNum(iGrow))
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
        call GenerateRotation_Reverse(TorNum(iGrow), GrowNum(iGrow), vTor, vPrev, torsTypes, vGrow, wTorsion(1))
        do iRosen = 2, nRosenTrials(nType)
          bendingAngles = 0d0
          dihedralAngles = 0d0
          do iAtom = 1, GrowNum(iGrow)
            call GenerateGaussianBondLength(rGrow(iAtom), k_bonds(iAtom), r_eqs(iAtom), &
                                    r_sigmas(iAtom), rmax_sqs(iAtom), Prob)
          enddo
          select case(GrowNum(iGrow))
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
          call GenerateRotation(TorNum(iGrow), GrowNum(iGrow), vTor, vPrev, rGrow, bendingAngles, dihedralAngles, &
                                torsTypes, vGrow, wTorsion(iRosen))
          do iAtom = 1, GrowNum(iGrow)
            trialPos_Branched(iAtom)%x = vGrow(iAtom)%x + molArray(nType)%mol(nMol)%x(atmCur)
            trialPos_Branched(iAtom)%y = vGrow(iAtom)%y + molArray(nType)%mol(nMol)%y(atmCur)
            trialPos_Branched(iAtom)%z = vGrow(iAtom)%z + molArray(nType)%mol(nMol)%z(atmCur)
          enddo
          E_Trial(iRosen) = 0d0
          do iAtom = 1, GrowNum(iGrow)
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
        E_Min = minval(E_Trial)
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = wBending(iRosen) * wTorsion(iRosen) * exp(-beta*(E_Trial(iRosen)-E_Min))         
        enddo
        rosenNorm = sum(ProbRosen)
        rosenRatio = rosenRatio*ProbRosen(1)/rosenNorm
        do iAtom = 1, GrowNum(iGrow)
          regrown(atmGrow(iAtom)) = .true.
          totalRegrown = totalRegrown + 1
        enddo
      enddo
	  
      end subroutine
