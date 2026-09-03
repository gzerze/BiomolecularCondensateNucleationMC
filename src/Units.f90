!===================================================================
! Module: Constants
! Purpose:
!   Defines global mathematical constants used throughout the code.
!
!   Constants are declared using double precision (`dp`) to ensure
!   numerical accuracy across the entire simulation.
!
! Dependencies:
!   - Uses VarPrecision for precision consistency
!
! Author:
!   Aliasghar Sepehri
!===================================================================
module Constants
  use VarPrecision
  implicit none

  !------------------------------------------------------
  ! Mathematical constants (double precision)
  real(dp), parameter :: pi      = 4.0d0 * datan(1.0d0)   ! π ≈ 3.141592653589793
  real(dp), parameter :: two_pi  = 8.0d0 * datan(1.0d0)   ! 2π ≈ 6.283185307179586

end module Constants

!===================================================================      
!==============================================================
! Module: Units
! Purpose: Provides unit conversion functions for energy, length, and angle
!==============================================================

module Units
  use VarPrecision
  contains

  !----------------------------------------------------------
  ! Converts various energy units to kB*T units (default internal unit)
  real(dp) function FindEngUnit(unitName)
    implicit none
    character(len=*), intent(in) :: unitName

    select case(trim(adjustl(unitName)))
      case("j-mol")
        FindEngUnit = 1d0 / 8.3144621d0
      case("kj-mol")
        FindEngUnit = 1d0 / 8.3144621d-3
      case("cal-mol")
        FindEngUnit = 1d0 / 1.9872041d0
      case("kcal-mol")
        FindEngUnit = 1d0 / 1.9872041d-3
      case("eV", "ev")
        FindEngUnit = 1d0 / 8.6173303d-5
      case("kB", "kb")
        FindEngUnit = 1d0
      case default
        write(*,*) "Error! Invalid Energy Unit Type!"
        write(*,*) unitName
        stop
    end select
  end function FindEngUnit

  !----------------------------------------------------------
  ! Converts length units to Angstroms (default internal unit)
  real(dp) function FindLengthUnit(unitName)
    implicit none
    character(len=*), intent(in) :: unitName

    select case(trim(adjustl(unitName)))
      case("nm")
        FindLengthUnit = 1d1
      case("a", "ang", "sigma")
        FindLengthUnit = 1d0
      case default
        write(*,*) "Error! Invalid Length Unit Type!"
        write(*,*) unitName
        stop
    end select
  end function FindLengthUnit

  !----------------------------------------------------------
  ! Converts angle units to radians (default internal unit)
  real(dp) function FindAngularUnit(unitName)
    use Constants
    implicit none
    character(len=*), intent(in) :: unitName

    select case(trim(adjustl(unitName)))
      case("deg", "degree", "degrees")
        FindAngularUnit = pi / 180d0
      case("rad", "radian", "radians")
        FindAngularUnit = 1d0
      case default
        write(*,*) "Error! Invalid Angular Unit Type!"
        write(*,*) unitName
        stop
    end select
  end function FindAngularUnit

end module Units
  
!===================================================================
