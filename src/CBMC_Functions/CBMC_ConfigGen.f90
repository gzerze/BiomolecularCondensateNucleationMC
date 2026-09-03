!=======================================================================
      subroutine Ridgid_ConfigGen(nType)
      use Coords
      use CoordinateFunctions
      use ForceField
      use Constants
      implicit none
      integer, intent(in) :: nType
      integer :: i

      real(dp) :: grnd,rotang
      real(dp) :: c_term,s_term
      real(dp) :: x_shift,y_shift,z_shift
      real(dp) :: x_rid_cm, y_rid_cm,z_rid_cm
      
      newMol%molType = nType
      newMol%x(1:nAtoms(nType)) = gasConfig(nType)%x(1:nAtoms(nType))
      newMol%y(1:nAtoms(nType)) = gasConfig(nType)%y(1:nAtoms(nType))
      newMol%z(1:nAtoms(nType)) = gasConfig(nType)%z(1:nAtoms(nType))

      x_rid_cm = gasConfig(nType)%x(1)
      y_rid_cm = gasConfig(nType)%y(1)
      z_rid_cm = gasConfig(nType)%z(1)      
      
      if(nAtoms(nType) .eq. 1) then
         return
      endif
      
!     Rotate xz axis
      rotang = two_pi*grnd()
      c_term = cos(rotang)
      s_term = sin(rotang)
      do i = 1,nAtoms(nType)
        x_shift = newMol%x(i) - x_rid_cm
        z_shift = newMol%z(i) - z_rid_cm        
        newMol%x(i) = c_term*x_shift - s_term*z_shift + x_rid_cm
        newMol%z(i) = s_term*x_shift + c_term*z_shift + z_rid_cm
      enddo   

        
!     Rotate yz axis
      rotang = two_pi*grnd()
      c_term = cos(rotang)
      s_term = sin(rotang)
      do i = 1,nAtoms(nType)
        y_shift = newMol%y(i) - y_rid_cm
        z_shift = newMol%z(i) - z_rid_cm        
        newMol%y(i) = c_term*y_shift - s_term*z_shift + y_rid_cm
        newMol%z(i) = s_term*y_shift + c_term*z_shift + z_rid_cm
      enddo   
  
!     Rotate xy axis
      rotang = two_pi*grnd()
      c_term = cos(rotang)
      s_term = sin(rotang)
      do i = 1,nAtoms(nType)
        x_shift = newMol%x(i) - x_rid_cm
        y_shift = newMol%y(i) - y_rid_cm        
        newMol%x(i) = c_term*x_shift - s_term*y_shift + x_rid_cm
        newMol%y(i) = s_term*x_shift + c_term*y_shift + y_rid_cm
      enddo    
      
      end subroutine
!=======================================================================
      subroutine Simple_ConfigGen(nType)
      use Coords
      use CoordinateFunctions
      use ForceField
      use Constants
      use CBMC_Variables
      implicit none
      integer, intent(in) :: nType
      integer :: bondType, bendType
      integer :: Atm1, Atm2, Atm3
!      real(dp) :: grnd
      real(dp) :: dx,dy,dz
      real(dp) :: k_bond, r_eq, r, Prob, r_sigma, rmax_sq
      real(dp) :: k_bend, ang_eq, ang
      type(SimpleAtomCoords) :: v1, v2
      
      newMol%molType = nType
      newMol%x = 0d0
      newMol%y = 0d0
      newMol%z = 0d0
      
      select case(nAtoms(nType))
      case(1)
        return
      case(2)
        bondType = bondArray(nType,1)%bondType
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        r_sigma = bondData(bondType)%r_sigma
        rmax_sq = bondData(bondType)%rmax_sq
        call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
        call Generate_UnitSphere(dx,dy,dz)
        newMol%x(2) = r*dx
        newMol%y(2) = r*dy
        newMol%z(2) = r*dz
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
        call Generate_UnitSphere(dx,dy,dz)
        v1%x = -r*dx
        v1%y = -r*dy
        v1%z = -r*dz
        newMol%x(atm2) = r*dx
        newMol%y(atm2) = r*dy
        newMol%z(atm2) = r*dz
        call FindBond(nType,Atm2, Atm3, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        r_sigma = bondData(bondType)%r_sigma
        rmax_sq = bondData(bondType)%rmax_sq
        call GenerateGaussianBondLength(r, k_bond, r_eq, r_sigma, rmax_sq, Prob)
        bendType = bendArray(nType,1)%bendType
        k_bend = bendData(bendType)%k_eq
        ang_eq = bendData(bendType)%ang_eq
!        call GenerateBendAngle(ang, k_bend, ang_eq, Prob)
        call Generate_UnitCone(v1,r,ang,v2)
        newMol%x(atm3) = newMol%x(atm2) + v2%x 
        newMol%y(atm3) = newMol%y(atm2) + v2%y
        newMol%z(atm3) = newMol%z(atm2) + v2%z
      case default
        stop "Error! Molecule has too many atoms for a simple regrowth"
      end select
      
       
      
      end subroutine

