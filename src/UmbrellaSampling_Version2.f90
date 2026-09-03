!==========================================================================
    module UmbrellaSamplingNew
    use SimParameters, only: maxMol, maxAtoms
    use VarPrecision
    use SimpleDistPair
    use UmbrellaTypes
    implicit none
    private

    integer, parameter :: maxLineLen = 500  

    logical :: useUmbrella = .false.
    logical :: energyAnalytics = .true.
    integer :: nBiasVariables = 0
    integer :: curUIndx, umbrellaLimit
    integer, allocatable :: binIndx(:)
    integer, allocatable :: binMax(:), binMin(:)
    integer, allocatable :: indexCoeff(:)
    integer, allocatable :: UArray(:)
    real(dp), allocatable :: UBias(:)
    real(dp), allocatable :: UHist(:)
    real(dp), allocatable :: UHistTotal(:)
    real(dp), allocatable :: UBinSize(:)
    real(dp), allocatable :: varValues(:)
    character(len=50):: inputFile
    character(len=10), allocatable :: outputFormat(:)
    character(len=100) :: screenFormat
    type(BiasVariablePointer), allocatable :: biasvar(:)
    type(BiasVariablePointer), allocatable :: biasvarnew(:)

    integer :: nDispFunc, nSwapInFunc, nSwapOutFunc
    type(DispUmbrellaArray), allocatable :: DispUmbrella(:)
    type(SwapInUmbrellaArray), allocatable :: SwapInUmbrella(:)
    type(SwapOutUmbrellaArray), allocatable :: SwapOutUmbrella(:)

    integer, parameter :: E_Bins = 1000
    real(dp), parameter :: dE = -100000/E_Bins
    real(dp), allocatable :: U_EAvg(:)
    real(dp), allocatable :: U_EHist(:, :)

!    public :: ReadInput_Umbrella
    public :: AllocateUmbrellaVariables,  AllocateUmbrellaArray, UmbrellaHistAdd
    public :: useUmbrella, OutputUmbrellaHist, GetUmbrellaBias_Disp, findVarValues, getBiasIndex
    public :: nBiasVariables, umbrellaLimit, UBias, UHist, UBinSize, outputFormat, curUIndx
    public :: DispUmbrella, SwapInUmbrella, SwapOutUmbrella, biasVar, biasVarnew
    public :: GetUmbrellaBias_SwapIn, GetUmbrellaBias_SwapOut, ScreenOutputUmbrella, screenFormat
    public :: CheckInitialValues, energyAnalytics, OutputUmbrellaAnalytics
    public :: ScriptInput_Umbrella, GetUmbrellaBias_Temperature
    public :: inputFile
!==========================================================================================
    contains
!==========================================================================================
subroutine AllocateUmbrellaVariables
  implicit none
  integer :: AllocateStatus
  integer :: i

  ! Allocate arrays that define binning parameters for umbrella sampling variables
  allocate( binMin(1:nBiasVariables),     STAT = AllocateStatus )
  allocate( binMax(1:nBiasVariables),     STAT = AllocateStatus )
  allocate( UBinSize(1:nBiasVariables),   STAT = AllocateStatus )

  ! Allocate and initialize pointers to bias variables (old state)
  allocate( biasvar(1:nBiasVariables),    STAT = AllocateStatus )
  do i = 1, nBiasVariables
    biasvar(i)%varType  = 0
    biasvar(i)%intVar   => null()
    biasvar(i)%realVar  => null()
  enddo

  ! Allocate and initialize pointers to bias variables (new state)
  allocate( biasvarnew(1:nBiasVariables), STAT = AllocateStatus )
  do i = 1, nBiasVariables
    biasvarnew(i)%varType = 0
    biasvarnew(i)%intVar  => null()
    biasvarnew(i)%realVar => null()
  enddo

  ! Allocate arrays for indexing and formatting umbrella sampling output
  allocate( binIndx(1:nBiasVariables),    STAT = AllocateStatus )
  allocate( indexCoeff(1:nBiasVariables), STAT = AllocateStatus )
  allocate( outputFormat(1:nBiasVariables), STAT = AllocateStatus )
  do i = 1, nBiasVariables
    outputFormat(i) = " "
  enddo

  ! Allocate and initialize function pointer arrays for displacement functions
  allocate( DispUmbrella(1:nDispFunc),    STAT = AllocateStatus )
  do i = 1, nDispFunc
    DispUmbrella(i)%func => null()
  enddo

  ! Allocate and initialize function pointer arrays for swap-in functions
  allocate( SwapInUmbrella(1:nSwapInFunc), STAT = AllocateStatus )
  do i = 1, nSwapInFunc
    SwapInUmbrella(i)%func => null()
  enddo

  ! Allocate and initialize function pointer arrays for swap-out functions
  allocate( SwapOutUmbrella(1:nSwapOutFunc), STAT = AllocateStatus )
  do i = 1, nSwapOutFunc
    SwapOutUmbrella(i)%func => null()
  enddo

  ! Stop the program if any allocation failed
  if (AllocateStatus /= 0) stop "*** Not enough memory ***"

  ! Log that umbrella allocation was successful
  write(34,*) "Umbrella variable arrays allocated successfully."

end subroutine AllocateUmbrellaVariables
!==========================================================================================
subroutine ScriptInput_Umbrella(inputLines)
!-----------------------------------------------------------------------
!  This subroutine parses the "umbrella" input block and sets up
!  bias variables, reference values, and umbrella potential handlers
!  for umbrella sampling.
!-----------------------------------------------------------------------

  use AnalysisMain, only: loadUmbArray, nAnalysisVar
  use MiscelaniousVars
  use SimpleDistPair, only: nDistPair, pairArrayIndx, CalcDistPairs_New, UmbrellaVar_DistPair
  use SimParameters, only: NMAX, NMIN, NPART, NPart_New, nMolTypes, maxMol, echoInput, &
                           NTotal, NTotal_New, temperature, TempNew
  use Q6Functions, only: q6ArrayIndx, CalcQ6_Disp, CalcQ6_SwapIn, CalcQ6_SwapOut, UmbrellaVar_Q6
  use ParallelVar, only: nout
  use WHAM_Module, only: refBin, refSizeNumbers

  implicit none

!-----------------------------------------------------------------------
!  Input arguments
!-----------------------------------------------------------------------
  character(len=maxLineLen), intent(in) :: inputLines(:)

!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
  integer :: nLines
  integer :: iUmbrella, AllocateStatus
  integer :: indxVar, stat
  integer :: iDisp, iSwapIn, iSwapOut
  real(dp) :: binSize, valMax, valMin
  real(dp), allocatable :: refVals(:)
  character(len=30) :: labelField
  character(len=30) :: umbrellaName

!-----------------------------------------------------------------------
!  Count number of umbrella bias variables from input lines
!-----------------------------------------------------------------------
  nLines = size(inputLines)
  nBiasVariables = nLines - 3

  if (nBiasVariables < 0) then
    write(nout,*) "ERROR! The user has specified an invalid number of Umbrella Sampling Variables"
    write(nout,*) labelField, nBiasVariables
    stop
  endif

!-----------------------------------------------------------------------
!  Umbrella sampling activation and reporting
!-----------------------------------------------------------------------
  if (nBiasVariables == 0) then
    useUmbrella = .false.
    write(nout,*) "Umbrella Sampling Used? = ", useUmbrella
    write(34,*)   "Umbrella Sampling Used? = ", useUmbrella
    return
  else
    useUmbrella = .true.
    write(nout,*) "Umbrella Sampling Used? = ", useUmbrella
    write(nout,*) "Number of Umbrella Variables:", nBiasVariables
    write(34,*)   "Umbrella Sampling Used? = ", useUmbrella
    write(34,*)   "Number of Umbrella Variables:", nBiasVariables
  endif

!-----------------------------------------------------------------------
!  Allocate reference values for bias variables
!-----------------------------------------------------------------------
  allocate(refVals(1:nBiasVariables), stat=AllocateStatus)
  if (AllocateStatus /= 0) then
    write(nout,*) "ERROR! Memory allocation failed for refVals"
    stop
  endif

!-----------------------------------------------------------------------
!  Read the file name for the initial Umbrella Sampling Bias
!-----------------------------------------------------------------------
  read(inputLines(1), *) labelField, inputFile

!-----------------------------------------------------------------------
!  Initial scan through the umbrella input lines to count how many 
!  displacement, swap-in, and swap-out functions are needed. 
!  This will be used to allocate corresponding arrays.
!-----------------------------------------------------------------------
  nDispFunc    = 0
  nSwapInFunc  = 0
  nSwapOutFunc = 0

  do iUmbrella = 1, nBiasVariables
    read(inputLines(iUmbrella + 2), *) umbrellaName

    select case( trim(adjustl(umbrellaName)) )
    
    case("clustersize")
      ! Cluster size does not require handler functions
      ! Cluster size is the mostly used variable
      continue

    case("totalclustersize")
      ! Total cluster size also does not require handler functions
      continue

    case("temperature")
      ! Temperature bias variable also handled separately
      continue

    case("pairdist")
      ! Pairwise distance variable uses all three handlers
      nDispFunc    = nDispFunc + 1
      nSwapInFunc  = nSwapInFunc + 1
      nSwapOutFunc = nSwapOutFunc + 1

    case("q6")
      ! Q6 order parameter requires all three umbrella handlers
      nDispFunc    = nDispFunc + 1
      nSwapInFunc  = nSwapInFunc + 1
      nSwapOutFunc = nSwapOutFunc + 1

    case("analysisvar")
      ! Generic analysis variable with all three handlers
      nDispFunc    = nDispFunc + 1
      nSwapInFunc  = nSwapInFunc + 1
      nSwapOutFunc = nSwapOutFunc + 1

    case default
      ! Invalid umbrella variable type
      write(nout,*) "ERROR! Invalid variable type specified in input file"
      write(nout,*) umbrellaName
      stop

    end select
  enddo

! Assign umbrella sampling variable pointers and setup sampling parameters
  call AllocateUmbrellaVariables

! Loop over each umbrella variable and process its setup
  iDisp     = 0
  iSwapIn   = 0
  iSwapOut  = 0

  do iUmbrella = 1, nBiasVariables
    read(inputLines(iUmbrella+2), *) umbrellaName

    select case( trim(adjustl(umbrellaName)) )

  ! --- CLUSTERSIZE: Sample cluster size for all molecule types ---
      case("clustersize")
        read(inputLines(iUmbrella+2), *) labelField, indxVar
        if(indxVar <= 0 .or. indxVar > nMolTypes) then
          write(*,*) "Error! Invalid molecule type for clustersize:"
          write(*,*) "  Defined Mol Types:", nMolTypes
          write(*,*) "  Chosen Mol Type:", indxVar
          stop
        endif
        biasvar(iUmbrella)%varType    = 1
        biasvar(iUmbrella)%intVar     => NPart(indxVar)
        biasvarnew(iUmbrella)%varType = 1
        biasvarnew(iUmbrella)%intVar  => NPart_New(indxVar)
        binMax(iUmbrella)             = NMAX(indxVar)
        binMin(iUmbrella)             = NMIN(indxVar)
        UBinSize(iUmbrella)           = 1.0_dp
        outputFormat(iUmbrella)       = "2x,F5.1,"

  ! --- TOTALCLUSTERSIZE: Sample the total cluster size over all types ---
      case("totalclustersize")
        biasvar(iUmbrella)%varType    = 1
        biasvar(iUmbrella)%intVar     => NTotal
        biasvarnew(iUmbrella)%varType = 1
        biasvarnew(iUmbrella)%intVar  => NTotal_New
        binMax(iUmbrella)             = maxMol
        binMin(iUmbrella)             = sum(NMin)
        UBinSize(iUmbrella)           = 1.0_dp
        outputFormat(iUmbrella)       = "2x,F5.1,"

  ! --- TEMPERATURE: Sample temperature range ---
      case("temperature")
        read(inputLines(iUmbrella+2), *) labelField, valMin, valMax, binSize
        biasvar(iUmbrella)%varType    = 2
        biasvar(iUmbrella)%realVar    => temperature
        biasvarnew(iUmbrella)%varType = 2
        biasvarnew(iUmbrella)%realVar => TempNew
        UBinSize(iUmbrella)           = binSize
        binMin(iUmbrella)             = nint(valMin / binSize)
        binMax(iUmbrella)             = nint(valMax / binSize)
        outputFormat(iUmbrella)       = "2x,F12.6,"

  ! --- PAIRDIST: Sample distance between specific atom pairs ---
      case("pairdist")
        indxVar = 0
        read(inputLines(iUmbrella+2), *) labelField, indxVar, valMin, valMax, binSize
        if(indxVar <= 0 .or. indxVar > nDistPair) then
          write(*,*) "Error! Invalid pairdist index:", indxVar
          stop
        endif
        biasvar(iUmbrella)%varType    = 2
        biasvar(iUmbrella)%realVar    => miscCoord(pairArrayIndx(indxVar))
        biasvarnew(iUmbrella)%varType = 2
        biasvarnew(iUmbrella)%realVar => miscCoord_New(pairArrayIndx(indxVar))
        UBinSize(iUmbrella)           = binSize
        binMin(iUmbrella)             = nint(valMin / binSize)
        binMax(iUmbrella)             = nint(valMax / binSize)
        outputFormat(iUmbrella)       = "2x,F12.6,"

        iDisp     = iDisp + 1
        iSwapIn   = iSwapIn + 1
        iSwapOut  = iSwapOut + 1
        DispUmbrella(iDisp)%func       => CalcDistPairs_New
        SwapInUmbrella(iSwapIn)%func   => CalcDistPairs_SwapIn
        SwapOutUmbrella(iSwapOut)%func => CalcDistPairs_SwapOut

  ! --- Q6: Sample the bond order parameter Q6 ---
      case("q6")
        read(inputLines(iUmbrella+2), *) labelField, valMin, valMax, binSize
        biasvar(iUmbrella)%varType    = 2
        biasvar(iUmbrella)%realVar    => miscCoord(q6ArrayIndx)
        biasvarnew(iUmbrella)%varType = 2
        biasvarnew(iUmbrella)%realVar => miscCoord_New(q6ArrayIndx)
        UBinSize(iUmbrella)           = binSize
        binMin(iUmbrella)             = nint(valMin / binSize)
        binMax(iUmbrella)             = nint(valMax / binSize)
    ! outputFormat(iUmbrella)    = "2x,F12.6," (optional)

        iDisp     = iDisp + 1
        iSwapIn   = iSwapIn + 1
        iSwapOut  = iSwapOut + 1
        DispUmbrella(iDisp)%func       => CalcQ6_Disp
        SwapInUmbrella(iSwapIn)%func   => CalcQ6_SwapIn
        SwapOutUmbrella(iSwapOut)%func => CalcQ6_SwapOut

  ! --- ANALYSISVAR: Custom user-defined analysis variable ---
      case("analysisvar")
        read(inputLines(iUmbrella+2), *) labelField, indxVar, valMin, valMax, binSize
        if(indxVar > nAnalysisVar) then
          write(*,*) "ERROR! Analysis variable index out of bounds:", indxVar
          write(*,*) inputLines(iUmbrella)
          stop
        endif
        call loadUmbArray(indxVar)%func(iUmbrella, indxVar, biasVar, biasVarNew, outputFormat, &
                                    iDisp, DispUmbrella, iSwapIn, SwapInUmbrella, iSwapOut, SwapOutUmbrella)
        UBinSize(iUmbrella) = binSize
        binMin(iUmbrella)   = nint(valMin / binSize)
        binMax(iUmbrella)   = nint(valMax / binSize)
        write(34,*) valMin, valMax, binSize
        write(34,*) binMin(iUmbrella), binMax(iUmbrella)

  ! --- Invalid Entry ---
      case default
        write(*,*) "ERROR! Invalid umbrella variable type:"
        write(*,*) umbrellaName
        stop
    end select
  enddo

  ! Construct the format string for screen or file output based on umbrella variable formats
  write(screenFormat,*) "(", (trim(adjustl(outputFormat(iUmbrella))), iUmbrella = 1, nBiasVariables), "1x)"

  ! Allocate arrays to store current umbrella variable values and umbrella bias potential
  allocate(varValues(1:nBiasVariables))
  allocate(UArray(1:nBiasVariables))

  ! Allocate the main umbrella sampling bias array (e.g., WHAM histogram)
  call AllocateUmbrellaArray

  ! Read the initial bias values from file (specified earlier in inputLines)
  call ReadInitialBias

! Use the input to specify the reference bin for umbrella sampling
  allocate(refSizeNumbers(1:nBiasVariables), STAT = AllocateStatus)

! Read the label and reference values from the second line of input
  read(inputLines(2), *) labelField, (refVals(iUmbrella), iUmbrella = 1, nBiasVariables)

! Convert the reference values to bin indices using the defined bin structure
  call getUIndexArray(refVals, refBin, stat)

! Check for errors in reference bin conversion
  if (stat .ne. 0) then
    if (stat .eq. 1) then
      write(34,*) "ERROR! Reference bin is above the maximum allowed value for one or more umbrella variables."
      write(34,*) trim(adjustl(inputLines(2)))
      stop
    elseif (stat .eq. -1) then
      write(34,*) "ERROR! Reference bin is below the minimum allowed value for one or more umbrella variables."
      write(34,*) trim(adjustl(inputLines(2)))
      stop
    else
      write(34,*) "ERROR! Invalid reference bin value provided for umbrella sampling variables."
      write(34,*) trim(adjustl(inputLines(2)))
      stop
    endif
  endif

! Store the raw reference values into refSizeNumbers for later use
  do iUmbrella = 1, nBiasVariables
    refSizeNumbers(iUmbrella) = refVals(iUmbrella)
  enddo

! Deallocate temporary storage
  deallocate(refVals)


end subroutine ScriptInput_Umbrella
!==========================================================================================
subroutine AllocateUmbrellaArray
  use ParallelVar
  implicit none

  integer :: i, j
  integer :: AllocateStatus

  !-------------------------------------------
  ! STEP 1: Set index coefficients
  ! These coefficients are used to convert multi-dimensional umbrella variables 
  ! into a single 1D array index (flattening the multidimensional space).
  !-------------------------------------------

  indexCoeff(1) = 1
  do i = 2, nBiasVariables
    indexCoeff(i) = 1
    do j = 1, i-1
      indexCoeff(i) = indexCoeff(i) + indexCoeff(j) * (binMax(j) - binMin(j))
    enddo
  enddo

  !---------------------------------------------------
  ! Example (2 umbrella variables):
  !
  !   binMin(1) = 0, binMax(1) = 3   → values: 0,1,2,3   → 4 bins
  !   binMin(2) = 1, binMax(2) = 3   → values: 1,2,3     → 3 bins
  !
  ! Then:
  !   indexCoeff(1) = 1
  !   indexCoeff(2) = 1 + indexCoeff(1) * (3 - 0) = 1 + 1*3 = 4
  !
  ! Mapping to 1D index:
  !   index(v1, v2) = indexCoeff(1)*(v1 - binMin(1)) 
  !                + indexCoeff(2)*(v2 - binMin(2))
  !
  ! Example:
  !   (v1=2, v2=3) → 1 + (2 - 0)*1 + (3 - 1)*4 = 1 + 2 + 8 = 11
  !---------------------------------------------------

  !-------------------------------------------
  ! STEP 2: Calculate total size of umbrella array
  !-------------------------------------------
  umbrellaLimit = 1
  do i = 1, nBiasVariables
    umbrellaLimit = umbrellaLimit + indexCoeff(i) * (int(binMax(i),4) - int(binMin(i),4))
  enddo

  write(nout,*) "Number of Umbrella Bins:", umbrellaLimit

  !-------------------------------------------
  ! STEP 3: Allocate umbrella bias and histogram arrays
  !-------------------------------------------
  allocate(UBias(1:umbrellaLimit+1), STAT = AllocateStatus)
  allocate(UHist(1:umbrellaLimit+1), STAT = AllocateStatus)
  UBias = 0E0_dp
  UHist = 0E0_dp

  !-------------------------------------------
  ! STEP 4: If energy analysis is enabled, allocate extra arrays
  !-------------------------------------------
  if (energyAnalytics) then
    allocate(U_EAvg(1:umbrellaLimit+1), STAT = AllocateStatus)
    allocate(U_EHist(1:umbrellaLimit+1, 0:E_Bins), STAT = AllocateStatus)
    allocate(UHistTotal(1:umbrellaLimit+1), STAT = AllocateStatus)
    U_EAvg = 0E0_dp
    U_EHist = 0E0_dp
    UHistTotal = 0E0_dp
  endif

  !-------------------------------------------
  ! Error check on allocation
  !-------------------------------------------
  if (AllocateStatus /= 0) stop "*** Not enough memory ***"

end subroutine AllocateUmbrellaArray
!==========================================================================================
subroutine ReadInitialBias
  implicit none

  integer :: AllocateStatus
  integer :: j, iInput, inStat, biasIndx
  real(dp), allocatable :: varValue(:)
  real(dp) :: curBias

  ! Open the file containing initial umbrella bias values
  open(unit=80, file=trim(adjustl(inputFile)))

  ! Allocate array to temporarily store variable values
  allocate(varValue(1:nBiasVariables), STAT = AllocateStatus)
  if (AllocateStatus /= 0) stop "*** Not enough memory ***"

  ! Initialize umbrella bias array
  UBias = 0.0_dp

  ! Loop over each line of the file (upper bound is an arbitrarily large number)
  do iInput = 1, nint(1d7)
    
    ! Read nBiasVariables real values + 1 bias value from the line
    read(80, *, IOSTAT=inStat) (varValue(j), j=1,nBiasVariables), curBias

    ! Exit loop if EOF or error encountered
    if (inStat < 0) exit

    ! Convert varValue → biasIndx in 1D array using index transformation
    call getUIndexArray(varValue, biasIndx, inStat)

    ! If variable values are out of valid bounds, skip this entry
    if (inStat == 1) cycle

    ! Assign the bias value to the correct index
    UBias(biasIndx) = curBias

  enddo

  ! Clean up
  deallocate(varValue)
  close(80)

end subroutine ReadInitialBias
!==========================================================================================
     subroutine UmbrellaHistAdd(E_T)
     implicit none
     real(dp), intent(in) :: E_T 

     integer :: bin

     curUIndx = getBiasIndex()
     UHist(curUIndx) = UHist(curUIndx) + 1E0_dp
!     write(*,*) curUIndx, UHist(curUIndx)

     if(energyAnalytics) then
       U_EAvg(curUIndx) = U_EAvg(curUIndx) + E_T
       UHistTotal(curUIndx) = UHistTotal(curUIndx) + 1E0_dp
       bin = floor(E_T / dE)

       if( (bin .ge. 0) .and. (bin .lt. E_Bins) ) then
         U_EHist(curUIndx, bin) = U_EHist(curUIndx, bin) + 1d0
       else
         U_EHist(curUIndx, E_Bins) = U_EHist(curUIndx, E_Bins) + 1d0
       endif

     endif

     end subroutine
!==========================================================================================
! Calculates the biasing potential difference for a displacement move in a Grand Canonical Monte
! Carlo nucleation simulation when umbrella sampling is enabled (useUmbrella=.true.). Returns
! biasDiff=0 if useUmbrella=.false. or if cluster size is the biasing variable (no change in NPART).
! Updates NPART_New, NTotal_New, and TempNew to reflect the trial configuration. Calls
! DispUmbrella functions to process the displacement and computes the new bias index (newUIndx).
! Sets rejMove=.true. if the new configuration is invalid. Compatible with all force fields
! (Mpipi, LJ_Q, HPS_single, HPS_piecewise, HPS_cation_pi) and cluster criteria.
subroutine GetUmbrellaBias_Disp(disp, biasDiff, rejMove)
  use CoordinateTypes, only: Displacement
  use SimParameters, only: NPART, NPART_new, NTotal, NTotal_New, temperature, TempNew
  use VarPrecision, only: dp
  implicit none

  ! Input/Output variables
  type(Displacement), intent(in) :: disp(:)  ! Displacement array for moved atoms
  real(dp), intent(out) :: biasDiff         ! Difference in biasing potential (new - old)
  logical, intent(out) :: rejMove           ! Flag to reject move if invalid

  ! Local variables
  integer :: iDispFunc, newUIndx, sizeDisp  ! Loop index, new bias index, displacement size
  real(dp) :: biasOld, biasNew              ! Old and new biasing potentials

  ! Step 1: Initialize outputs and check umbrella sampling
  rejMove = .false.
  biasDiff = 0.0_dp
  if (.not. useUmbrella) return  ! No bias if umbrella sampling is disabled

  ! Step 2: Update trial configuration variables
  NPART_New = NPART
  NTotal_New = NTotal
  TempNew = temperature
  if (size(NPART_New) /= size(NPART) .or. size(NPART) > maxMol) then
    write(35, *) "Error: NPART_New or NPART size mismatch in GetUmbrellaBias_Disp"
    rejMove = .true.
    return
  endif

  ! Step 3: Get old bias
  curUIndx = getBiasIndex()
  if (curUIndx < 1 .or. curUIndx > size(UBias)) then
    write(35, *) "Error: Invalid curUIndx=", curUIndx, " in GetUmbrellaBias_Disp"
    rejMove = .true.
    return
  endif
  biasOld = UBias(curUIndx)

  ! Step 4: Process displacement with umbrella functions
  if (nDispFunc > 0) then
    sizeDisp = size(disp)
    if (sizeDisp < 1 .or. sizeDisp > maxAtoms) then
      write(35, *) "Error: Invalid disp size=", sizeDisp, " in GetUmbrellaBias_Disp"
      rejMove = .true.
      return
    endif
    do iDispFunc = 1, nDispFunc
      call DispUmbrella(iDispFunc)%func(disp(1:sizeDisp))
    enddo
  endif

  ! Step 5: Get new bias index and check validity
  call getNewBiasIndex(newUIndx, rejMove)
  if (rejMove) return
  biasNew = UBias(newUIndx)

  ! Step 6: Calculate bias difference
  biasDiff = biasNew - biasOld
end subroutine GetUmbrellaBias_Disp
!==========================================================================================
  ! Computes the umbrella sampling bias difference for an insertion move in a grand
  ! canonical ensemble nucleation simulation. Calculates biasDiff = f_bias(n+1) - f_bias(n),
  ! where f_bias is the biasing potential based on cluster size (NTotal_New) or molecule
  ! counts (NPART_new). Sets rejMove to .true. if the new configuration violates umbrella
  ! constraints (e.g., cluster size limits). Supports all force fields and molecule types.
subroutine GetUmbrellaBias_SwapIn(biasDiff, rejMove)
  use SimParameters, only:temperature, TempNew
  use VarPrecision, only: dp
  implicit none

  ! Output variables
  logical, intent(out) :: rejMove         ! Rejection flag for umbrella constraints
  real(dp), intent(out) :: biasDiff       ! Bias difference: f_bias(n+1) - f_bias(n)

  ! Local variables
  integer :: iSwapFunc                    ! Loop index for swap-in functions
  integer :: newUIndx                     ! New bias index for proposed state
  real(dp) :: biasOld, biasNew            ! Old and new bias values

  ! Skip if umbrella sampling is disabled
  if (.not. useUmbrella) then
    biasDiff = 0.0_dp                    ! No bias applied
    return
  endif

  ! Initialize temperature and rejection flag
  TempNew = temperature                  ! Set proposed temperature
  rejMove = .false.                      ! Assume valid configuration

  ! Get old bias for current state
  curUIndx = getBiasIndex()              ! Current bias index (NTotal, NPART)
  biasOld = UBias(curUIndx)              ! Current bias value

  ! Apply swap-in umbrella functions (if any)
  if (nSwapInFunc /= 0) then
    do iSwapFunc = 1, nSwapInFunc
      call SwapInUmbrella(iSwapFunc) % func  ! Execute swap-in function
    enddo
  endif

  ! Get new bias for proposed state
  call getNewBiasIndex(newUIndx, rejMove)  ! New bias index (NTotal_New, NPART_new)
  if (rejMove) then
    return                              ! Reject if constraints violated
  endif
  biasNew = UBias(newUIndx)             ! New bias value

  ! Compute bias difference
  biasDiff = biasNew - biasOld          ! f_bias(n+1) - f_bias(n)

end subroutine GetUmbrellaBias_SwapIn
!==========================================================================================
! Calculates the umbrella sampling bias difference for a deletion move of molecule
! nMol of type nType, computing biasDiff = -β * [f_bias(n-1) - f_bias(n)] for use in
! AVBMC_EBias_Rosen_Out’s deletion acceptance probability. Calls swap-out functions
! and updates bias indices. Correctness is ensured by:
! - Skipping calculation if useUmbrella = .false.
! - Applying multiple swap-out umbrella functions if nSwapOutFunc > 0.
! - Validating inputs and bias indices.
! Supports all force fields (LJ_Q, Mpipi, HPS_single, HPS_piecewise, HPS_cation_pi)
! and molecule types (rigid, small, linear, branched).
subroutine GetUmbrellaBias_SwapOut(nType, nMol, biasDiff, rejMove)
  use VarPrecision, only: dp
  implicit none

  integer, intent(in) :: nType, nMol
  real(dp), intent(out) :: biasDiff
  logical, intent(out) :: rejMove

  integer :: iSwapFunc, curUIndx, newUIndx
  real(dp) :: biasOld, biasNew

  ! Skip if umbrella sampling is disabled
  ! Correctness: No bias contribution
  if (.not. useUmbrella) then
    biasDiff = 0.0_dp
    rejMove = .false.
    return
  endif

  ! Initialize rejection flag
  rejMove = .false.

  ! Get current bias index and old bias
  ! Correctness: Retrieves bias for current state
  curUIndx = getBiasIndex()
  biasOld = UBias(curUIndx)

  ! Apply swap-out umbrella functions if needed
  ! Correctness: Updates state for deletion
  if (nSwapOutFunc > 0) then
    do iSwapFunc = 1, nSwapOutFunc
      call SwapOutUmbrella(iSwapFunc)%func(nType, nMol)
    enddo
  endif

  ! Get new bias index for post-deletion state
  ! Correctness: Checks validity of new state
  call getNewBiasIndex(newUIndx, rejMove)
  if (rejMove) return

  ! Compute bias difference
  ! Correctness: biasDiff = biasNew - biasOld = -β * [f_bias(n-1) - f_bias(n)]
  biasNew = UBias(newUIndx)
  biasDiff = biasNew - biasOld
end subroutine GetUmbrellaBias_SwapOut
!==========================================================================================
     subroutine GetUmbrellaBias_Temperature(biasDiff, rejMove)
     use SimParameters, only: NPART, NPART_new, NTotal, NTotal_New
     implicit none
     logical, intent(out) :: rejMove
     real(dp), intent(out) :: biasDiff
     integer :: newUIndx
     real(dp) :: biasNew, biasOld

     if(.not. useUmbrella) then
       biasDiff = 0E0_dp
       return
     endif

     rejMove = .false.
     curUIndx = getBiasIndex()
     biasOld = UBias(curUIndx)

     NPART_New = NPART
     NTotal_New = NTotal

!     do iUmbrella = 1, nBiasVariables
!       if(biasvarnew(iUmbrella) % varType .eq. 1) then
!         biasvarnew(iUmbrella) % intVar = biasvar(iUmbrella) % intVar 
!       elseif(biasvarnew(iUmbrella) % varType .eq. 2) then
!         biasvarnew(iUmbrella) % realVar = biasvar(iUmbrella) % realVar
!       endif
!     enddo

     call getNewBiasIndex(newUIndx, rejMove)
     if(rejMove) then
       return
     endif
     biasNew = UBias(newUIndx)
     biasDiff = biasNew - biasOld

 
     end subroutine
!==========================================================================================
    subroutine ScreenOutputUmbrella
    use ParallelVar
    implicit none
    integer :: iBias, j

    do iBias = 1, nBiasVariables
      if(biasvar(iBias)%varType .eq. 1) then
        varValues(iBias) = real(biasvar(iBias)%intVar,dp)
      elseif(biasvar(iBias)%varType .eq. 2) then
        varValues(iBias) = biasvar(iBias)%realVar
      endif
    enddo
    
    write(nout,screenFormat) (varValues(j), j=1,nBiasVariables)

    end subroutine
!==========================================================================================
    subroutine OutputUmbrellaHist
    implicit none
    integer :: iUmbrella, iBias
!    integer, allocatable :: UArray(:)
    character(len = 100) :: outputString


!    allocate(UArray(1:nBiasVariables))
!    allocate(varValues(1:nBiasVariables))

    write(outputString, *) "(", (trim(outputFormat(iBias)), iBias =1,nBiasVariables), "2x, F18.1)"
    open(unit=60, file="TemporaryHist.txt")
    do iUmbrella = 1, umbrellaLimit
      call findVarValues(iUmbrella, UArray)
      do iBias = 1, nBiasVariables
        varValues(iBias) = real( UArray(iBias), dp) * UBinSize(iBias)
      enddo
      if(UHist(iUmbrella) .ne. 0E0_dp) then
        write(60,outputString) (varValues(iBias), iBias =1,nBiasVariables), UHist(iUmbrella)
      endif
    enddo 
    flush(60)
    close(60)

!    deallocate(UArray)
!    deallocate(varValues)

    end subroutine

!==========================================================================================
    subroutine OutputUmbrellaAnalytics
    use ParallelVar
    use MPI
    implicit none
    integer :: iUmbrella, iBias, iBin
    integer :: arraySize
    real(dp), allocatable :: TempHist(:), TempAvg(:), Temp2D(:,:)

    if(.not. useUmbrella) then
      return
    endif


!    if(myid .eq. 0) then
      allocate( TempHist(1:umbrellaLimit+1) ) 
      allocate( TempAvg(1:umbrellaLimit+1) ) 
      allocate( Temp2D(1:umbrellaLimit+1, 0:E_Bins) )
      TempHist = 0E0_dp
      TempAvg = 0E0_dp
      Temp2D = 0E0_dp
!    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierror) 
    arraySize = size(U_EAvg)   
!    write(*,*) arraySize, size(TempHist)
    call MPI_REDUCE(U_EAvg, TempAvg, arraySize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 

    call MPI_REDUCE(UHistTotal, TempHist, arraySize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 
!    write(*,*) myid, ierror
!    do iUmbrella = 1, umbrellaLimit
!      write(*,*) iUmbrella, U_EAvg(iUmbrella), TempHist(iUmbrella)
!    enddo

    if(myid .eq. 0) then
      do iUmbrella = 1, umbrellaLimit
        U_EAvg(iUmbrella) = TempAvg(iUmbrella)
        UHistTotal(iUmbrella) = TempHist(iUmbrella)
      enddo
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)    
    arraySize = size(U_EHist)   
    call MPI_REDUCE(U_EHist, Temp2D, arraySize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)    
!    write(*,*) myid, ierror   

    if(myid .eq. 0) then
      do iUmbrella = 1, umbrellaLimit
        do iBin = 0, E_Bins
          U_EHist(iUmbrella, iBin) = Temp2D(iUmbrella, iBin)
        enddo
      enddo
    endif




    if(myid .eq. 0) then
      open(unit=60, file = "Umbrella_AvgE.txt")
      do iUmbrella = 1, umbrellaLimit
        if(UHistTotal(iUmbrella) .ne. 0E0_dp) then
          call findVarValues(iUmbrella, UArray)
          do iBias = 1, nBiasVariables
            varValues(iBias) = real( UArray(iBias), dp) * UBinSize(iBias)
          enddo
          write(60, *) (varValues(iBias), iBias =1,nBiasVariables), U_EAvg(iUmbrella)/UHistTotal(iUmbrella)
        endif
      enddo
      close(60)

      open(unit=60, file = "Umbrella_DensityStates.txt")
      do iUmbrella = 1, umbrellaLimit
        if(UHistTotal(iUmbrella) .ne. 0E0_dp) then
          call findVarValues(iUmbrella, UArray)
          do iBias = 1, nBiasVariables
            varValues(iBias) = real( UArray(iBias), dp) * UBinSize(iBias)
          enddo
          write(60,*) (varValues(iBias), iBias =1,nBiasVariables)
          do iBin = 0, E_Bins
            if(U_EHist(iUmbrella, iBin) .ne. 0d0) then
              write(60, *) iBin*dE, U_EHist(iUmbrella, iBin)
            endif
          enddo
        endif
      enddo
      close(60)
    endif


    deallocate(TempHist) 
    deallocate(TempAvg)
    deallocate(Temp2D)

    end subroutine
!==========================================================================
function getBiasIndex() result(biasIndx)
  !---------------------------------------------------------------------------
  ! Computes the umbrella sampling bias index for the current configuration.
  !
  ! This function uses the current values pointed to by `biasvar()`:
  !   - If the umbrella variable is an integer (e.g., cluster size), use it directly
  !   - If the umbrella variable is real-valued, bin it using its bin size
  !
  ! OUTPUT:
  !   biasIndx - 1D index into umbrella sampling arrays (UBias, UHist, etc.)
  !
  ! NOTE:
  !   - The type of variable is identified using `varType`
  !     varType = 1 → integer variable
  !     varType = 2 → real variable
  !
  ! EXAMPLE:
  !   For cluster size sampling:
  !     → biasvar(iBias)%intVar points to NPart(indxVar), which stores the cluster size
  !     → The corresponding biasIndx is computed via indexCoeff and binMin
  !---------------------------------------------------------------------------

  integer :: biasIndx
  integer :: iBias

  !---------------------------------------------
  ! STEP 1: Fill binIndx array based on type
  !         of each umbrella sampling variable
  !---------------------------------------------
  do iBias = 1, nBiasVariables
    if (biasvar(iBias)%varType == 1) then
      ! Integer variable (e.g., cluster size, total cluster size)
      binIndx(iBias) = biasvar(iBias)%intVar

    elseif (biasvar(iBias)%varType == 2) then
      ! Real-valued variable (e.g., temperature, distance)
      binIndx(iBias) = floor(biasvar(iBias)%realVar / UBinSize(iBias))

    else
      write(*,*) "ERROR! Unknown umbrella variable type:", biasvar(iBias)%varType
      stop
    endif
  enddo

  !---------------------------------------------
  ! STEP 2: Compute unique 1D index
  !         using precomputed index coefficients
  !---------------------------------------------
  biasIndx = 1
  do iBias = 1, nBiasVariables
    biasIndx = biasIndx + indexCoeff(iBias) * (binIndx(iBias) - binMin(iBias))
  enddo

end function getBiasIndex
!==========================================================================
subroutine getNewBiasIndex(biasIndx, rejMove)
  implicit none

  !> Output: rejection flag
  logical, intent(out) :: rejMove

  !> Output: index into 1D umbrella array
  integer, intent(out) :: biasIndx

  integer :: iBias

  ! Assume move is acceptable initially
  rejMove = .false.

  ! Loop through all umbrella sampling variables
  do iBias = 1, nBiasVariables

    ! Case 1: Integer-type umbrella variable (e.g., cluster size)
    if (biasvarnew(iBias) % varType == 1) then
      binIndx(iBias) = biasvarnew(iBias) % intVar

    ! Case 2: Real-type umbrella variable (e.g., distance, Q6)
    elseif (biasvarnew(iBias) % varType == 2) then
      binIndx(iBias) = floor( biasvarnew(iBias) % realVar / UBinSize(iBias) )

    endif

    ! Reject if new value falls below min bin
    if (binIndx(iBias) < binMin(iBias)) then
      rejMove = .true.
      return
    endif

    ! Reject if new value exceeds max bin
    if (binIndx(iBias) > binMax(iBias)) then
      rejMove = .true.
      return
    endif

  enddo

  ! Compute the linearized index from the bin indices
  biasIndx = 1
  do iBias = 1, nBiasVariables
    biasIndx = biasIndx + indexCoeff(iBias) * (binIndx(iBias) - binMin(iBias))
  enddo

  ! Reject if index exceeds allocated limit
  if (biasIndx > umbrellaLimit) then
    rejMove = .true.
    return
  endif

end subroutine getNewBiasIndex
!==========================================================================
subroutine getUIndexArray(varArray, biasIndx, stat)
  !---------------------------------------------------------------------------
  ! Maps a set of umbrella variable values to a unique 1D bias index
  !
  ! INPUT:
  !   varArray  - real array of umbrella variable values
  !
  ! OUTPUT:
  !   biasIndx  - the corresponding 1D index for the umbrella bias/histogram
  !   stat      - status flag:
  !                = 0 : OK
  !                = 1 : At least one variable is above the allowed max
  !                = -1: At least one variable is below the allowed min
  !---------------------------------------------------------------------------

  implicit none
  real(dp), intent(in)  :: varArray(:)
  integer,  intent(out) :: biasIndx, stat

  integer :: iBias

  stat = 0

  !---------------------------------------------
  ! STEP 1: Compute bin index for each variable
  !---------------------------------------------
  do iBias = 1, nBiasVariables
    binIndx(iBias) = nint(varArray(iBias) / UBinSize(iBias))

    if (binIndx(iBias) > binMax(iBias)) then
      stat = 1
      return
    elseif (binIndx(iBias) < binMin(iBias)) then
      stat = -1
      return
    endif
  enddo

  !---------------------------------------------
  ! STEP 2: Map bin indices to 1D bias index
  !         using precomputed index coefficients
  !---------------------------------------------
  biasIndx = 1
  do iBias = 1, nBiasVariables
    biasIndx = biasIndx + indexCoeff(iBias) * (binIndx(iBias) - binMin(iBias))
  enddo

  !---------------------------------------------
  ! EXAMPLE:
  !   Suppose:
  !     binMin(1) = 0, binMax(1) = 3 → bin values: 0,1,2,3
  !     binMin(2) = 1, binMax(2) = 3 → bin values: 1,2,3
  !     indexCoeff(1) = 1
  !     indexCoeff(2) = 4
  !
  !   Given varArray = [1.0, 2.0]
  !     → binIndx(1) = 1
  !     → binIndx(2) = 2
  !     → biasIndx = 1 + 1*(1 - 0) + 4*(2 - 1) = 1 + 1 + 4 = 6
  !   → biasIndx = 6
  !---------------------------------------------

end subroutine getUIndexArray
!===========================================================================
subroutine findVarValues(UIndx, UArray)
  !---------------------------------------------------------------------------
  ! Converts a 1D umbrella index (UIndx) into the multi-dimensional bin
  ! indices (UArray) for each umbrella sampling variable.
  !
  ! This is the reverse of the getBiasIndex() function.
  !
  ! INPUT:
  !   UIndx  - 1D index used in UBias, UHist, etc.
  !
  ! OUTPUT:
  !   UArray - Array containing the bin value for each umbrella variable
  !
  ! EXAMPLE:
  !   Suppose we have 2 variables:
  !     binMin(1) = 0, binMax(1) = 3  → values: 0, 1, 2, 3
  !     binMin(2) = 1, binMax(2) = 3  → values: 1, 2, 3
  !
  !   Let’s assume indexCoeff(1) = 1
  !                      indexCoeff(2) = 4
  !
  !   For UIndx = 6:
  !     remainder = UIndx - 1 = 5
  !
  !     i = 1 → iBias = 2 → curVal = 5 / 4 = 1 → value = 1 + 1 = 2
  !     remainder = 5 - 1*4 = 1
  !
  !     i = 2 → iBias = 1 → curVal = 1 / 1 = 1 → value = 1 + 0 = 1
  !     remainder = 1 - 1*1 = 0
  !
  !   Output: UArray = [1, 2] → values of var1 and var2
  !---------------------------------------------------------------------------

  implicit none
  integer, intent(in)    :: UIndx         ! 1D umbrella index
  integer, intent(inout) :: UArray(:)     ! Output variable values (bins)
  integer :: i, iBias
  integer :: remainder, curVal

  ! Adjust for Fortran 1-based indexing
  remainder = UIndx - 1

  ! Reverse mapping from 1D index to individual umbrella variable bins
  do i = 1, nBiasVariables
    iBias = nBiasVariables - i + 1

    ! Get how many times indexCoeff(iBias) fits into the current remainder
    curVal = int(real(remainder, dp) / real(indexCoeff(iBias), dp))

    ! Reconstruct the umbrella variable value by adding binMin
    UArray(iBias) = curVal + binMin(iBias)

    ! Update remainder to reflect how much has been accounted for
    remainder = remainder - curVal * indexCoeff(iBias)
  enddo

end subroutine findVarValues
!==========================================================================================
     subroutine CheckInitialValues
     implicit none
     integer :: iBias

     do iBias = 1, nBiasVariables
       if(biasvar(iBias) % varType .eq. 1) then
         binIndx(iBias) = nint( biasvar(iBias) % intVar / UBinSize(iBias) )
       elseif(biasvar(iBias) % varType .eq. 2) then
         binIndx(iBias) = floor( biasvar(iBias) % realVar / UBinSize(iBias) )
       endif
       
       if(binIndx(iBias) .gt. binMax(iBias)) then
         write(*,*) "The initital system state is above the upper bounds"         
         write(*,*) "specified by the umbrella sampling input"   
         write(*,*) "Umbrella Variable:", iBias
         write(*,*) "Bin Index:", binIndx(iBias)  
         write(*,*) "Largest allowed bin:", binMax(iBias)          
         stop
       endif
       if(binIndx(iBias) .lt. binMin(iBias)) then
         write(*,*) "The initital system state is below the lower bounds"         
         write(*,*) "specified by the umbrella sampling input"
         write(*,*) "Umbrella Variable:", iBias
         write(*,*) "Bin Index:", binIndx(iBias)  
         write(*,*) "Smallest allowed bin:", binMin(iBias)              
         stop
       endif
     enddo
 
     end subroutine
!==========================================================================
    end module
!==========================================================================


