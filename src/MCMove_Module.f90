!======================================================
module MoveTypeModule
  use VarPrecision

  ! This module links each Monte Carlo (MC) move routine to the master move array.
  ! To integrate new MC move subroutines, import them here and ensure they conform
  ! to the MCMoveSub interface.

  use AVBMC_Module, only: AVBMC
  use CBMC_Module, only: CBMC
  use Exchange_Module, only: Exchange
  use LLAVBMC_Module, only: LLAVBMC
  use LLCBMC_Module, only: LLCBMC
  use SimpleMCMoves_Module, only: Translation, Rotation, SingleAtom_Translation, &
                                   TemperatureMove, LLTranslation, LLRotation

  integer, parameter :: maxLineLen = 500

  ! Abstract interface for Monte Carlo move subroutines.
  abstract interface
    subroutine MCMoveSub(E_T, acc_x, atmp_x)
      use VarPrecision
      implicit none
      real(dp), intent(inout) :: E_T      ! Total energy of the cluster
      real(dp), intent(inout) :: acc_x    ! Counter for accepted moves
      real(dp), intent(inout) :: atmp_x   ! Counter for attempted moves
    end subroutine
  end interface

  ! Type that holds a pointer to a move subroutine
  type MoveArray
    procedure(MCMoveSub), pointer, nopass :: moveFunction => NULL()
  end type

  ! Flags to check if AVBMC or CBMC are used in the simulation
  logical :: avbmcUsed, cbmcUsed

  ! Number of Monte Carlo move types defined
  integer :: nMoveTypes

  ! Arrays to manage MC moves and their statistics
  type(MoveArray), allocatable :: mcMoveArray(:)      ! Array of MC move subroutine pointers
  real(dp), allocatable       :: moveProbability(:)   ! Probability of each move
  real(dp), allocatable       :: movesAccepted(:)     ! Acceptance count per move
  real(dp), allocatable       :: movesAttempt(:)      ! Attempt count per move
  real(dp), allocatable, target :: accptRate(:)       ! Acceptance rate per move
  character(len=35), allocatable :: moveName(:)       ! Name of each move type

  contains
!======================================================
      subroutine CalcAcceptanceRates
      implicit none
      integer :: iMoves

      do iMoves = 1, nMoveTypes
        accptRate(iMoves) = movesAccepted(iMoves)/movesAttempt(iMoves)
      enddo

      end subroutine CalcAcceptanceRates

!======================================================
subroutine ScriptInput_MCMove(lines)
  !---------------------------------------------------------
  ! This subroutine parses the user-defined MC move types
  ! and initializes:
  !   - Function pointers to the move routines
  !   - Move probabilities
  !   - Acceptance and attempt counters
  !   - Move labels
  !---------------------------------------------------------

  use AcceptRates, only: acptInSize, atmpInSize
  use AVBMC_Module, only: swapProb
  use SimParameters, only: nMolTypes, NMIN, NMAX, maxMol
  use VarPrecision, only: dp
  use ParallelVar,     only: nout
  implicit none

  character(len=maxLineLen), intent(in) :: lines(:)  ! Input block of lines
  integer :: nLines                                 ! Total number of input lines
  integer :: i, iMoves                              ! Loop indices
  integer :: AllocateStatus                         ! Status for memory allocation

  real(dp) :: norm                                  ! Normalization factor for probabilities
  character(len=30) :: moveName_temp                ! Temporary storage for move name

  !---------------------------------------------------------
  ! Determine number of defined move types (skip headers)
  !---------------------------------------------------------
  nLines = size(lines)
  nMoveTypes = nLines - 2  ! Skip header and footer lines

  ! Log number of moves
  write(34,*) "Number of Monte Carlo moves specified:", nMoveTypes

  !---------------------------------------------------------
  ! Basic input validation
  !---------------------------------------------------------
  if(nMoveTypes .le. 0) then
    write(nout,*) "ERROR! The user has specified an invalid number of Monte Carlo moves"
    write(nout,*) "Please specify at least one valid Monte Carlo move to continue"
    do i = 1, nLines
      write(*,*) lines(i)
    enddo
    write(34,*) "ERROR: Invalid number of MC moves, aborting."
    stop
  endif
  !---------------------------------------------------------
  ! Allocate arrays for:
  !   - Pointers to move routines
  !   - Move probabilities
  !   - Acceptance and attempt counters
  !   - Acceptance rates
  !   - Move names
  !---------------------------------------------------------
  allocate(mcMoveArray(1:nMoveTypes),      STAT = AllocateStatus)
  allocate(moveProbability(1:nMoveTypes),  STAT = AllocateStatus)
  allocate(movesAccepted(1:nMoveTypes),    STAT = AllocateStatus)
  allocate(movesAttempt(1:nMoveTypes),     STAT = AllocateStatus)
  allocate(accptRate(1:nMoveTypes),        STAT = AllocateStatus)
  allocate(moveName(1:nMoveTypes),         STAT = AllocateStatus)

  if (AllocateStatus /= 0) then
    write(nout,*) "ERROR: Memory allocation failed for move type arrays."
    write(34,*) "ERROR: Memory allocation failed in ScriptInput_MCMove."
    stop
  endif

  write(34,*) "MC move arrays allocated successfully for", nMoveTypes, "moves."

  !------------------------------------------------------------
  ! Loop through each user-specified move
  !   - Parse move name and probability
  !   - Normalize later
  !   - Assign corresponding function pointer
  !------------------------------------------------------------
  norm = 0.0_dp
  avbmcUsed = .false.
  cbmcUsed = .false.

  do iMoves = 1, nMoveTypes
    read(lines(iMoves+1), *) moveName_temp, moveProbability(iMoves)
    norm = norm + moveProbability(iMoves)
    call LowerCaseLine(moveName_temp)

    select case( trim(adjustl(moveName_temp)) )

    case("translation")
      mcMoveArray(iMoves) % moveFunction => Translation
      moveName(iMoves) = "Translation"

    case("lltranslation")
      mcMoveArray(iMoves) % moveFunction => LLTranslation
      moveName(iMoves) = "LLTranslation"

    case("rotation")
      mcMoveArray(iMoves) % moveFunction => Rotation
      moveName(iMoves) = "Rotation"

    case("llrotation")
      mcMoveArray(iMoves) % moveFunction => LLRotation
      moveName(iMoves) = "LLRotation"

    case("avbmc")
      mcMoveArray(iMoves) % moveFunction => AVBMC
      moveName(iMoves) = "AVBMC"
      allocate(acptInSize(1:maxMol), STAT = AllocateStatus)
      allocate(atmpInSize(1:maxMol), STAT = AllocateStatus)
      avbmcUsed = .true.

    case("llavbmc")
      mcMoveArray(iMoves) % moveFunction => LLAVBMC
      moveName(iMoves) = "LLAVBMC"
      allocate(acptInSize(1:maxMol), STAT = AllocateStatus)
      allocate(atmpInSize(1:maxMol), STAT = AllocateStatus)
      avbmcUsed = .true.

    case("cbmc")
      mcMoveArray(iMoves) % moveFunction => CBMC
      moveName(iMoves) = "CBMC"
      cbmcUsed = .true.

    case("llcbmc")
      mcMoveArray(iMoves) % moveFunction => LLCBMC
      moveName(iMoves) = "LLCBMC"
      cbmcUsed = .true.

    case("exchange")
      mcMoveArray(iMoves) % moveFunction => Exchange
      moveName(iMoves) = "Exchange"

    case("temperature")
      mcMoveArray(iMoves) % moveFunction => TemperatureMove
      moveName(iMoves) = "Temperature"

    case("singleatom_translation")
      mcMoveArray(iMoves) % moveFunction => SingleAtom_Translation
      moveName(iMoves) = "Single Atom Translation"

    case default
      write(nout,*) "ERROR! Invalid move type specified in input file."
      write(nout,*) trim(adjustl(moveName_temp)), moveProbability(iMoves)
      write(34,*) "ERROR: Invalid MC move type encountered:", trim(adjustl(moveName_temp))
      stop
    end select

    write(34,*) "Move", iMoves, ":", trim(moveName(iMoves)), "with probability", moveProbability(iMoves)

  enddo


  !------------------------------------------------------------
  ! Normalize move probabilities so their sum equals 1
  ! Then convert them to cumulative form for selection
  ! Example: input = [0.4, 0.3, 0.3] → cumulative = [0.4, 0.7, 1.0]
  !------------------------------------------------------------

  do iMoves = 1, nMoveTypes
    moveProbability(iMoves) = moveProbability(iMoves) / norm
  enddo

  if(nMoveTypes > 1) then
    do iMoves = 2, nMoveTypes
      moveProbability(iMoves) = moveProbability(iMoves) + moveProbability(iMoves - 1)
    enddo
  endif

  ! Write cumulative move probabilities to output unit 34
  write(34,*) "Cumulative Move Probabilities:"
  do iMoves = 1, nMoveTypes
    write(34,*) trim(moveName(iMoves)), ":", moveProbability(iMoves)
  enddo
  !------------------------------------------------------------
  ! Setup swap probabilities for AVBMC move if used
  ! Only molecule types with NMIN ≠ NMAX are eligible for swaps
  !------------------------------------------------------------

  if(avbmcUsed .eqv. .true.) then
    allocate(swapProb(1:nMolTypes), stat = AllocateStatus)
    swapProb = 0.0_dp

    ! Assign preliminary swap probability = 1.0 if type is eligible
    do i = 1, nMolTypes
      if(NMIN(i) .ne. NMAX(i)) then
        swapProb(i) = 1.0_dp
      endif
    enddo

    ! Normalize swap probabilities
    norm = sum(swapProb)

    if(norm == 0.0_dp) then 
      write(nout,*) "ERROR! AVBMC has been used, but no swap moves can occur"
      write(nout,*) "since NMIN is equal to NMAX for all molecule types"
      write(nout,*) "NMIN:", NMIN
      write(nout,*) "NMAX:", NMAX
      stop
    endif

    write(35,*) "Probability of swapping molecule type i:"
    do i = 1, nMolTypes
      swapProb(i) = swapProb(i) / norm
      write(35,*) i, swapProb(i)
    enddo
  endif



end subroutine ScriptInput_MCMove
!======================================================
end module MoveTypeModule
!======================================================
