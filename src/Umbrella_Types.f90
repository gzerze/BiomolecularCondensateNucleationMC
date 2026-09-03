!=====================================================================================
module UmbrellaTypes
  use VarPrecision

  !---------------------------------------------------------------------
  ! Module: UmbrellaTypes
  ! Contains type definitions and interfaces for umbrella sampling bias
  ! variable handling and umbrella update routines.
  !---------------------------------------------------------------------

  ! Bias variable can be either an integer or a real variable,
  ! depending on the umbrella sampling scheme
  type :: BiasVariablePointer
    integer :: varType                          ! Type identifier for the variable (e.g., 1 for int, 2 for real)
    integer, pointer :: intVar                  ! Pointer to an integer bias variable
    real(dp), pointer :: realVar                ! Pointer to a real bias variable
  end type BiasVariablePointer


  !---------------------------------------------------------------------
  ! Abstract interface for umbrella displacement update
  ! This subroutine operates on a displacement array
  !---------------------------------------------------------------------
  abstract interface
    subroutine UDispFunc(disp)
      use CoordinateTypes
      implicit none
      type(Displacement), intent(in) :: disp(:)
    end subroutine UDispFunc
  end interface

  !---------------------------------------------------------------------
  ! Abstract interface for umbrella swap-out operation
  ! Called when a molecule is removed from the system
  !---------------------------------------------------------------------
  abstract interface
    subroutine USwapOutFunc(nType, nMol)
      use CoordinateTypes
      implicit none
      integer, intent(in) :: nType   ! Molecule type
      integer, intent(in) :: nMol    ! Molecule index
    end subroutine USwapOutFunc
  end interface


  !---------------------------------------------------------------------
  ! Arrays of procedure pointers for displacement, swap-in, and swap-out
  ! umbrella handlers. These are dynamically assigned in the main input
  ! processing subroutines.
  !---------------------------------------------------------------------

  type :: DispUmbrellaArray
    procedure(UDispFunc), pointer, nopass :: func  ! Pointer to displacement umbrella function
  end type DispUmbrellaArray

  type :: SwapInUmbrellaArray
    procedure(), pointer, nopass :: func           ! Generic pointer to swap-in function (no args)
  end type SwapInUmbrellaArray

  type :: SwapOutUmbrellaArray
    procedure(USwapOutFunc), pointer, nopass :: func  ! Pointer to swap-out umbrella function
  end type SwapOutUmbrellaArray

end module UmbrellaTypes

!=====================================================================================
!===========================================================
!  Analysis and Umbrella Sampling – Code Documentation
!===========================================================

! ---------------
! DISPLACEMENT STRUCTURE
! ---------------
! The following type is used to represent the displacement of a specific atom 
! during Monte Carlo (MC) trial moves.
!
! molType   : Molecule type
! atmIndx   : Atom index within a molecule
! molIndx   : Molecule index
! Displaced : Whether this atom was displaced
! x_new, y_new, z_new : New trial coordinates
! x_old, y_old, z_old : Pointers to original coordinates

! type :: Displacement
!   integer(atomIntType) :: molType
!   integer(atomIntType) :: atmIndx
!   integer(atomIntType) :: molIndx
!   logical :: Displaced
!   real(dp) :: x_new, y_new, z_new
!   real(dp), pointer :: x_old, y_old, z_old
! end type Displacement

! ---------------
! BIAS VARIABLE STRUCTURE (FOR UMBRELLA SAMPLING)
! ---------------
! Each umbrella variable is either:
!   - Integer (varType = 1)
!   - Real/float (varType = 2)
! This structure holds the bias variable used in umbrella sampling.

! type :: BiasVariablePointer
!   integer :: varType
!   integer, pointer :: intVar
!   real(dp), pointer :: realVar
! end type BiasVariablePointer

! ---------------
! COLLECTIVE VARIABLE (CV) HANDLERS FOR UMBRELLA SAMPLING
! ---------------
! These abstract interfaces are used to register CVs that respond to:
!   - displacement moves,
!   - swap-in moves,
!   - swap-out moves.

! abstract interface
!   subroutine UDispFunc(disp)
!     use CoordinateTypes
!     implicit none
!     type(Displacement), intent(in) :: disp(:)
!   end subroutine
! end interface

! abstract interface
!   subroutine USwapOutFunc(nType, nMol)
!     use CoordinateTypes
!     implicit none
!     integer, intent(in) :: nType, nMol
!   end subroutine
! end interface

! type :: DispUmbrellaArray
!   procedure(UDispFunc), pointer, nopass :: func
! end type

! type :: SwapInUmbrellaArray
!   procedure(), pointer, nopass :: func
! end type

! type :: SwapOutUmbrellaArray
!   procedure(USwapOutFunc), pointer, nopass :: func
! end type

! ---------------
! WHERE THESE CVS ARE USED
! ---------------
! - CVs such as pairdist, q3, q4, q6 may be used in umbrella sampling.
! - Others such as radialdistribution and radialdensity are purely histogram-based 
!   and cannot act as umbrella biasing variables.

! ---------------
! ANALYSIS FRAMEWORK
! ---------------
! Analysis functions can be added to the simulation in a modular way. 
! This is managed through an interface and function arrays.

! interface
!   subroutine UmbrellaLoader(iUmbrella, varIndx, biasVar, biasVarNew, outputFormat, &
!                             iDisp, DispUmbrella, iSwapIn, SwapInUmbrella, iSwapOut, SwapOutUmbrella)
!     use MiscelaniousVars
!     use UmbrellaTypes
!     implicit none
!     integer, intent(in) :: iUmbrella, varIndx
!     integer, intent(inout) :: iDisp, iSwapIn, iSwapOut
!     type(BiasVariablePointer), intent(inout) :: biasVar(:), biasVarNew(:)
!     type(DispUmbrellaArray), intent(inout) :: DispUmbrella(:)
!     type(SwapInUmbrellaArray), intent(inout) :: SwapInUmbrella(:)
!     type(SwapOutUmbrellaArray), intent(inout) :: SwapOutUmbrella(:)
!     character(len=10), intent(inout) :: outputFormat(:)
!   end subroutine
! end interface

! type :: UmbrellaLoaderArray
!   procedure(UmbrellaLoader), pointer, nopass :: func
! end type

! type(UmbrellaLoaderArray), allocatable :: loadUmbArray(:)

! type :: AnalysisFunctionArray
!   procedure(), pointer, nopass :: func
! end type

! postMoveArray(:) is called after each Monte Carlo move (e.g., update CVs)
! type(AnalysisFunctionArray), allocatable :: postMoveArray(:)
!           if(useAnalysis) then !Property Calculation Variables
!             call PostMoveAnalysis
!           endif

! outputArray(:) is called when writing output files
! type(AnalysisFunctionArray), allocatable :: outputArray(:)
!      if(myid .eq. 0) then
!        if(useAnalysis) then
!          call OutputAnalysis
!        endif
!      endif

! Example syntax for analysis block:
!   radialdistribution: type1, atom1, type2, atom2, binSize, nBins, filename
!   pairdist: type1, mol1, atom1, type2, mol2, atom2
!   radialdensity: type1, binSize, nBins, filename
!   q3 / q4 / q6: (requires only distance cutoff)

! IMPORTANT: If you use analysis, set useDistStore = .true.

! ---------------
! UMBRELLA SAMPLING OVERVIEW
! ---------------
! Supported umbrella variables:
!   - clustersize
!   - totalclustersize
!   - temperature
!   - pairdist
!   - q6
!   - analysisvar

! Most common use:
!   clustersize → number of molecules of a given type in the cluster
!
!   Example:
!     biasvar(iUmbrella)%varType    = 1
!     biasvar(iUmbrella)%intVar     => NPart(indxVar)
!     biasvarnew(iUmbrella)%varType = 1
!     biasvarnew(iUmbrella)%intVar  => NPart_New(indxVar)
!     binMax(iUmbrella)             = NMAX(indxVar)
!     binMin(iUmbrella)             = NMIN(indxVar)
!     UBinSize(iUmbrella)           = 1.0_dp
!     outputFormat(iUmbrella)       = "2x,F5.1,"

! ---------------
! MULTIVARIATE → 1D INDEX TRANSFORMATION
! ---------------
! If you bias on multiple umbrella variables (e.g., two cluster sizes), the code
! flattens the multidimensional bin grid into a 1D array:
!
! For each umbrella variable `i`, we define:
!   - binMin(i), binMax(i): min/max integer bin values
!   - indexCoeff(i): stride for this variable's contribution to 1D index
!
! The 1D index is calculated as:
!   biasIndx = 1 + Σ_i [ indexCoeff(i) * (binIndx(i) - binMin(i)) ]
!
! Example (2D biasing):
!   binMin(1) = 0, binMax(1) = 3     → values: 0,1,2,3
!   binMin(2) = 1, binMax(2) = 3     → values: 1,2,3
!   UBinSize(:) = 1.0_dp
!
!   From these:
!     indexCoeff(1) = 1
!     indexCoeff(2) = 4   ! (max(1)-min(1)+1) = 4
!
!   Now if varArray = [2.0, 3.0] → binIndx = [2,3]
!   biasIndx = 1 + 1*(2-0) + 4*(3-1) = 1 + 2 + 8 = 11
!
! Reverse mapping (to get bin values from UIndx) is handled by findVarValues()

!===========================================================
