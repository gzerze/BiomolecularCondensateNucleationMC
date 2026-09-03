!======================================================
! Module: VarPrecision
! Purpose:
!   Defines precision types used consistently across the codebase.
!
! Contents:
!   - dp           : kind parameter for double-precision real numbers
!   - atomIntType  : kind parameter for default integer type used in atom indexing
!
! Range and Precision:
!   For dp (real(kind=dp)):
!     - Smallest positive number     = tiny(1.0_dp)       ≈ 2.225074e-308
!     - Largest representable number = huge(1.0_dp)       ≈ 1.797693e+308
!     - Decimal precision            = precision(1.0_dp)  ≈ 15–17 digits
!
!   For atomIntType (integer(kind=atomIntType)):
!     - Largest representable value  = huge(1_atomIntType) ≈ 2,147,483,647
!     - Smallest representable value = -huge(1_atomIntType) - 1 ≈ -2,147,483,648
!
! Usage:
!   Use `use VarPrecision` in any module or subroutine where precise
!   real or integer types are required.
!
! Author:
!   Aliasghar Sepehri
!======================================================
module VarPrecision
  implicit none

  ! Double-precision real kind parameter (typically 8 bytes)
  integer, parameter :: dp = kind(0.0d0)

  ! Integer kind parameter for atom indexing (typically 4 bytes)
  integer, parameter :: atomIntType = kind(0)

end module VarPrecision
!======================================================
