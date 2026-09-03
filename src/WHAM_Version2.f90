!=========================================================================
      module WHAM_Functions
      use WHAM_Module
      implicit none

!=========================================================================
      contains
!=========================================================================
!---------------------------------------------------------------
! Subroutine: WHAM_Initialize
! Purpose   : Allocate and initialize arrays needed for
!             Weighted Histogram Analysis Method (WHAM)
!---------------------------------------------------------------
subroutine WHAM_Initialize
  use ParallelVar, only: myid, nout
  use UmbrellaSamplingNew, only: umbrellaLimit
  implicit none

  integer :: AllocateStatus

  ! Only the root processor handles global WHAM arrays and output
  if (myid .eq. 0) then

    ! Allocate global WHAM arrays (on root only)
    allocate(WHAM_Numerator(1:umbrellaLimit),               STAT = AllocateStatus)
    allocate(WHAM_Denominator(1:umbrellaLimit, 1:nWhamItter+1), STAT = AllocateStatus)
    allocate(HistStorage(1:umbrellaLimit),                  STAT = AllocateStatus)
    allocate(BiasStorage(1:umbrellaLimit, 1:nWhamItter+1),  STAT = AllocateStatus)
    allocate(FreeEnergyEst(1:umbrellaLimit),                STAT = AllocateStatus)
    allocate(ProbArray(1:umbrellaLimit),                    STAT = AllocateStatus)
    allocate(LogProbArray(1:umbrellaLimit),                 STAT = AllocateStatus)

    write(nout,*) "Allocated WHAM Variables"

    ! Initialize arrays to zero
    WHAM_Numerator  = 0E0_dp
    WHAM_Denominator = 0E0_dp
    HistStorage     = 0E0_dp
    BiasStorage     = 0E0_dp
    FreeEnergyEst   = 0E0_dp
    ProbArray       = 0E0_dp

    ! Initialize WHAM iteration counter
    nCurWhamItter = 1

    ! Open files to store WHAM intermediate results
    open(unit = 95, file = "WHAM_InstantHists.txt")
    open(unit = 96, file = "WHAM_TempHist.incomp")
    open(unit = 97, file = "WHAM_Potential.incomp")
    open(unit = 98, file = "WHAM_Mid_DG.incomp")
  endif

  ! Allocate local arrays for all processes
  allocate(TempHist(1:umbrellaLimit), STAT = AllocateStatus)
  allocate(NewBias(1:umbrellaLimit),  STAT = AllocateStatus)

  ! Initialize local arrays
  TempHist = 0E0_dp
  NewBias  = 0E0_dp

end subroutine WHAM_Initialize
!=========================================================================
! Adjusts the umbrella sampling bias using the Weighted Histogram Analysis Method (WHAM) in a
! parallel Grand Canonical Monte Carlo nucleation simulation with MPI. Collects histogram data
! (UHist) from all processes, updates the biasing potential f_bias(n) to achieve a flat histogram,
! and computes the unbiased probability P(n) and free energy F(n) = -k_B * T * ln[P(n)/P(n_ref)].
! For biomolecules, WHAM combines histograms across iterations to address sampling challenges due
! to low acceptance rates for insertion/deletion moves. Updates:
!   f_bias_new(n) = f_bias_old(n) - k_B * T * ln[N(n)/N(n_ref)]
!   P(n) = [ Σ_t N_t(n)^2 / N̅_t ] / [ Σ_t N_t(n) * exp(f_bias^t(n) - f_t) ]
!   f_t = ln[ Σ_n P(n) * exp(-f_bias^t(n)) ]
! Called periodically when useUmbrella=.true. to refine the bias for cluster size sampling.
subroutine WHAM_AdjustHist
  use ParallelVar, only: myid, nout, ierror
  use UmbrellaSamplingNew, only: useUmbrella, UHist, UBias, umbrellaLimit
  use VarPrecision, only: dp
  use MPI
  implicit none

  ! Local variables
  integer :: i, j, cnt                           ! Loop indices and iteration counter
  integer :: maxbin, maxbin2                     ! Indices of maximum histogram bins
  integer :: arraySize                           ! Size of arrays for MPI communication
  real(dp) :: norm                               ! Normalization factor for histograms
  real(dp) :: denomSum, fSum                     ! (kept) sums for WHAM, now used as work vars
  real(dp) :: tol                                ! Convergence tolerance for self-consistent loop
  real(dp) :: refBias                            ! Reference bias for rescaling
  real(dp), allocatable :: F_Estimate(:)         ! Current estimate of f_t
  real(dp), allocatable :: F_Old(:)              ! Previous estimate of f_t
  real(dp) :: amax, maxLogP, val, logDen, logfSum

  ! Step 1: Exit if umbrella sampling is disabled
  if (.not. useUmbrella) return

  ! Step 2: Synchronize processes and collect histogram data
  write(nout, *) "Halting for WHAM"
  call MPI_BARRIER(MPI_COMM_WORLD, ierror)

  ! Reduce UHist from all processes to TempHist on root (myid = 0)
  arraySize = size(UHist)
  if (myid == 0) TempHist = 0.0_dp
  call MPI_REDUCE(UHist, TempHist, arraySize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

  ! Step 3: Process histogram data on root process
  if (myid == 0) then
    ! Initialize WHAM terms
    norm = sum(TempHist)
    do i = 1, umbrellaLimit
      BiasStorage(i, nCurWhamItter) = UBias(i)
      if (TempHist(i) > 0.0_dp) then
        WHAM_Numerator(i) = WHAM_Numerator(i) + TempHist(i) * TempHist(i) / norm

        ! ------------------------------------------------------------
        ! Store WHAM denominator in LOG form:
        !   ln[ N_t(n) * exp(f_bias^t(n)) ] = ln N_t(n) + f_bias^t(n)
        ! This prevents overflow when f_bias becomes very large.
        ! ------------------------------------------------------------
        WHAM_Denominator(i, nCurWhamItter) = log(TempHist(i)) + UBias(i)

        HistStorage(i) = HistStorage(i) + TempHist(i)
      else
        WHAM_Denominator(i, nCurWhamItter) = -huge(1.0_dp)
      endif
    enddo

    ! Allocate arrays for self-consistent WHAM solution
    allocate(F_Estimate(nCurWhamItter), F_Old(nCurWhamItter))
    F_Estimate = 0.0_dp
    F_Old      = 0.0_dp
    ProbArray  = 0.0_dp
    LogProbArray = -huge(1.0_dp) ! −∞ = log(0)
    tol = tolLimit + 1.0_dp
    cnt = 0

    ! Step 3.1: Solve WHAM equations self-consistently
    do while (tol > tolLimit)
      cnt = cnt + 1
      if (cnt > maxSelfConsist) then
        write(nout, *) "WHAM: Self-consistent limit reached"
        exit
      endif

      ! Store previous F_Estimate
      F_Old = F_Estimate

      ! Compute unbiased probabilities P(n) in log space to avoid underflow
      ! LogProbArray(i) ∝ ln(WHAM_Numerator(i)) - ln( Σ_j WHAM_Denominator(i,j)*exp(-F_j) )
    !---------------------------------------------------------------------------
    ! Numerical stability note:
    !
    ! In WHAM we must evaluate sums of exponentials of the form
    !       Σ_j  exp( log(WHAM_Denominator(i,j)) - F_Estimate(j) ).
    !
    ! The argument of the exponential can be very large or very small:
    !   - Overflow occurs if exp(x) with x > ~709  → Inf
    !   - Underflow occurs if exp(x) with x < ~ -745 → 0.0
    !
    ! For biomolecular nucleation, BiasStorage(i,j) - F_Estimate(j) can easily
    ! exceed these limits (free energies > 700 k_BT, or probabilities < 1e-300),
    ! making direct evaluation of exp(...) numerically unsafe.
    !
    ! To prevent overflow/underflow, we use the standard "log-sum-exp" technique:
    !   1. Find amax = max_j [ log(WHAM_Denominator(i,j)) - F_Estimate(j) ].
    !   2. Rewrite the sum as:
    !          Σ_j exp(a_j) = exp(amax) * Σ_j exp(a_j - amax),
    !      where (a_j - amax) ≤ 0 for all j, so exp(a_j - amax) ≤ 1.
    !
    ! This guarantees:
    !   - No overflow (argument to exp never exceeds 0)
    !   - No catastrophic underflow (small terms naturally drop out)
    !   - The logarithm of the sum can be computed safely as:
    !          logDen = amax + log(Σ_j exp(a_j - amax))
    !
    ! The variable amax therefore acts as a numerical stabilizer that keeps all
    ! exponentials in a safe floating-point range, avoiding Inf/NaN in WHAM.
    !---------------------------------------------------------------------------

      do i = 1, umbrellaLimit
        if (WHAM_Numerator(i) > 0.0_dp) then
          amax = -huge(1.0_dp)
          ! First pass: find max exponent for log-sum-exp
          do j = 1, nCurWhamItter
            if (WHAM_Denominator(i, j) > -0.5_dp*huge(1.0_dp)) then
              val = WHAM_Denominator(i, j) - F_Estimate(j)
              if (val > amax) amax = val
            endif
          enddo
          if (amax > -0.5_dp*huge(1.0_dp)) then
            denomSum = 0.0_dp
            do j = 1, nCurWhamItter
              if (WHAM_Denominator(i, j) > -0.5_dp*huge(1.0_dp)) then
                val = WHAM_Denominator(i, j) - F_Estimate(j)
                denomSum = denomSum + exp(val - amax)
              endif
            enddo
            logDen = amax + log(denomSum)        ! ln denominator
            LogProbArray(i) = log(WHAM_Numerator(i)) - logDen
          else
            LogProbArray(i) = -huge(1.0_dp)
          endif
        else
          LogProbArray(i) = -huge(1.0_dp)
        endif
      enddo

      ! Convert LogProbArray → ProbArray and normalize:
      !   ProbArray(i) = exp(LogProbArray(i) - maxLogP) / Σ_k exp(LogProbArray(k) - maxLogP)
! computes the maximum value of the array LogProbArray, but only over valid elements (not equal to −huge)
      maxLogP = maxval(LogProbArray, mask = (LogProbArray > -0.5_dp*huge(1.0_dp)))
      ProbArray = 0.0_dp
      if (maxLogP > -0.5_dp*huge(1.0_dp)) then
        do i = 1, umbrellaLimit
          if (LogProbArray(i) > -0.5_dp*huge(1.0_dp)) then
            ProbArray(i) = exp(LogProbArray(i) - maxLogP)
          else
            ProbArray(i) = 0.0_dp
          endif
        enddo
        norm = sum(ProbArray)
        if (norm > 0.0_dp) ProbArray = ProbArray / norm
      endif

      ! Update F_Estimate (f_t) using log-sum-exp over i
      do j = 1, nCurWhamItter
        amax = -huge(1.0_dp)
        ! First pass: max exponent
        do i = 1, umbrellaLimit
          if (ProbArray(i) > 0.0_dp) then
            val = log(ProbArray(i)) + BiasStorage(i, j)
            if (val > amax) amax = val
          endif
        enddo
        if (amax > -0.5_dp*huge(1.0_dp)) then
          denomSum = 0.0_dp
          do i = 1, umbrellaLimit
            if (ProbArray(i) > 0.0_dp) then
              val = log(ProbArray(i)) + BiasStorage(i, j)
              denomSum = denomSum + exp(val - amax)
            endif
          enddo
          logfSum = amax + log(denomSum)
          F_Estimate(j) = logfSum
          ! Dampen updates for stability
          F_Estimate(j) = 0.5_dp * (F_Estimate(j) + F_Old(j))
        endif
      enddo

      ! Check convergence
      tol = sum(abs(F_Estimate - F_Old))
    enddo

    ! Step 3.2: Compute free energy and update bias
    NewBias = 0.0_dp
    maxbin = maxloc(HistStorage, dim=1)  ! Reference bin for cumulative histogram
    maxbin2 = maxloc(TempHist, dim=1)    ! Reference bin for current histogram

    ! Update free energy and bias
    do i = 1, umbrellaLimit
      ! Free energy from LogProbArray: F(i) = -ln[P(i)/P(maxbin)] = -(LogProbArray(i)-LogProbArray(maxbin))
      if (LogProbArray(i) > -0.5_dp*huge(1.0_dp) .and. LogProbArray(maxbin) > -0.5_dp*huge(1.0_dp)) then
        FreeEnergyEst(i) = -(LogProbArray(i) - LogProbArray(maxbin))
      else
        FreeEnergyEst(i) = huge(1.0_dp)
      endif

      if (TempHist(i) >= 1.0_dp) then
        NewBias(i) = UBias(i) - UBias(maxbin2) - log(TempHist(i) / TempHist(maxbin2))
      else
        NewBias(i) = UBias(i) - UBias(maxbin2) + log(TempHist(maxbin2))
      endif
    enddo

    ! Rescale bias and free energy to set reference bin to zero
    refBias = NewBias(refBin)
    NewBias = NewBias - refBias
    refBias = FreeEnergyEst(refBin)
    FreeEnergyEst = FreeEnergyEst - refBias

    ! Deallocate temporary arrays
    deallocate(F_Estimate, F_Old)
  endif

  ! Step 4: Synchronize and distribute new bias
  call MPI_BARRIER(MPI_COMM_WORLD, ierror)
  arraySize = size(NewBias)
  call MPI_BCAST(NewBias, arraySize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

  ! Step 5: Update local bias and reset histogram
  UBias = NewBias
  UHist = 0.0_dp
  if (myid == 0) call WHAM_MidSimOutput
  nCurWhamItter = nCurWhamItter + 1
end subroutine WHAM_AdjustHist
!==================================================================================
! Outputs intermediate WHAM results to files in a Grand Canonical Monte Carlo nucleation simulation.
! Writes four files: instantaneous histogram (WHAM_InstantHistogram.txt), total histogram
! (WHAM_TotalHistogram.txt), current umbrella bias (WHAM_Bias.txt), and free energy estimate
! (WHAM_FreeEnergy.txt). For each file, outputs bias variable values (scaled by UBinSize) and
! corresponding data (TempHist, HistStorage, UBias, FreeEnergyEst) using a formatted string.
! Skips zero entries for total histogram and free energy where appropriate. Called by
! WHAM_AdjustHist to log WHAM progress, typically on the root process in parallel (MPI) mode.
! Supports multiple bias variables and is compatible with all force fields and simulation moves.
subroutine WHAM_MidSimOutput
  use VarPrecision, only: dp
  use UmbrellaSamplingNew
  use SwapBoundary
  implicit none

  ! Local variables
  integer :: i, j                            ! Loop indices
  integer :: UArray(nBiasVariables)          ! Array of bias variable indices
  character(len=100) :: formatString         ! Dynamic format string for output
  integer, parameter :: histFile = 95        ! File unit for instantaneous histogram
  integer, parameter :: totalHistFile = 96   ! File unit for total histogram
  integer, parameter :: biasFile = 97        ! File unit for umbrella bias
  integer, parameter :: freeEnergyFile = 98  ! File unit for free energy

  ! Step 1: Build format string for output
  write(formatString, *) "(", (trim(outputFormat(j)), j =1,nBiasVariables), "2x, F18.1)"

  ! Step 2: Write instantaneous histogram to WHAM_InstantHistogram.txt
  do i = 1, umbrellaLimit
    call findVarValues(i, UArray)
    write(histFile, formatString) (real(UArray(j), dp)*UBinSize(j), j=1,nBiasVariables), TempHist(i)
  enddo
  write(histFile, *)  ! Add separators
  write(histFile, *)
  write(histFile, *)
  flush(histFile)

  ! Step 3: Write total histogram to WHAM_TotalHistogram.txt
  rewind(totalHistFile)
  do i = 1, umbrellaLimit
    if (HistStorage(i) > 0.0_dp) then
      call findVarValues(i, UArray)
      write(totalHistFile, formatString) (real(UArray(j), dp)*UBinSize(j), j=1,nBiasVariables), HistStorage(i)
    endif
  enddo
  flush(totalHistFile)

  write(formatString, *) "(", (trim(outputFormat(j)), j =1,nBiasVariables), "2x, F18.10)"

  ! Step 4: Write current umbrella bias to WHAM_Bias.txt
  rewind(biasFile)
  do i = 1, umbrellaLimit
    call findVarValues(i, UArray)
    write(biasFile, formatString) (real(UArray(j), dp)*UBinSize(j), j=1,nBiasVariables), UBias(i)
  enddo
  flush(biasFile)

  ! Step 5: Write free energy estimate to WHAM_FreeEnergy.txt
  rewind(freeEnergyFile)
  do i = 1, umbrellaLimit
    if (LogProbArray(i) > -0.5_dp*huge(1.0_dp)) then
      call findVarValues(i, UArray)
      write(freeEnergyFile, formatString) (real(UArray(j), dp)*UBinSize(j), j=1,nBiasVariables), FreeEnergyEst(i)
    endif
  enddo
  flush(freeEnergyFile)
end subroutine WHAM_MidSimOutput
!==================================================================================
! Finalizes WHAM by outputting results to files in a Grand Canonical Monte Carlo nucleation
! simulation. Writes free energy (WHAM_FreeEnergy.txt), final bias potential (WHAM_FinalBias.txt),
! normalized probabilities (WHAM_Probabilities.txt), overall histogram (WHAM_OverallHistogram.txt),
! and, if WHAM_ExtensiveOutput=.true., histogram numerator and denominator
! (WHAM_Hist_Numerator.txt, WHAM_Hist_Denominator.txt). Outputs bias variable values (scaled by
! UBinSize) and corresponding data, skipping zero entries where appropriate. Executes on the root
! process (myid == 0) in parallel (MPI) mode. Called at simulation end when useWham=.true.
! Supports multiple bias variables and all force fields and simulation moves.
subroutine WHAM_Finalize
  use VarPrecision, only: dp
  use ParallelVar, only: myid
  use UmbrellaSamplingNew
  use SwapBoundary
  implicit none

  ! Local variables
  integer :: i, j                            ! Loop indices
  integer :: UArray(nBiasVariables)          ! Array of bias variable indices
  real(dp) :: probNorm                       ! Normalization for probabilities
  real(dp) :: refBias                        ! Reference bias for rescaling
  character(len=100) :: formatPrecise        ! Format string for precise output (F18.10)
  character(len=100) :: formatHist           ! Format string for histogram output (F18.1)
  integer, parameter :: outputFile = 92      ! File unit for main outputs
  integer, parameter :: histFile = 36        ! File unit for histogram outputs

  ! Step 1: Exit if WHAM is disabled
  if (.not. useWham) return

  ! Step 2: Build format strings for output
  write(formatPrecise, *) "(", (trim(outputFormat(j)), j =1,nBiasVariables), "2x, F18.10)"
  write(formatHist, *) "(", (trim(outputFormat(j)), j =1,nBiasVariables), "2x, F18.1)"

  ! Step 3: Output results on root process
  if (myid == 0) then
    ! Step 3.1: Write free energy to WHAM_FreeEnergy.txt
    open(unit=outputFile, file="WHAM_FreeEnergy.txt")
    refBias = FreeEnergyEst(refBin)
    do i = 1, umbrellaLimit
      if (LogProbArray(i) > -0.5_dp*huge(1.0_dp)) then
        call findVarValues(i, UArray)
        write(outputFile, formatPrecise) (real(UArray(j), dp)*UBinSize(j), j=1,nBiasVariables), &
                                         FreeEnergyEst(i) - refBias
      endif
    enddo
    close(outputFile)

    ! Step 3.2: Write final bias potential to WHAM_FinalBias.txt
    open(unit=outputFile, file="WHAM_FinalBias.txt")
    refBias = UBias(refBin)
    do i = 1, umbrellaLimit
      call findVarValues(i, UArray)
      write(outputFile, formatPrecise) (real(UArray(j), dp)*UBinSize(j), j=1,nBiasVariables), &
                                       UBias(i) - refBias
    enddo
    close(outputFile)

    ! Step 3.3: Write normalized probabilities to WHAM_Probabilities.txt
    open(unit=outputFile, file="WHAM_Probabilities.txt")
    probNorm = sum(ProbArray)
    do i = 1, umbrellaLimit
      if (ProbArray(i) > 0.0_dp) then
        call findVarValues(i, UArray)
        write(outputFile, formatPrecise) (real(UArray(j), dp)*UBinSize(j), j=1,nBiasVariables), &
                                         ProbArray(i) / probNorm
      endif
    enddo
    close(outputFile)

    ! Step 3.4: Write overall histogram to WHAM_OverallHistogram.txt
    open(unit=histFile, file="WHAM_OverallHistogram.txt")
    do i = 1, umbrellaLimit
      if (HistStorage(i) > 0.0_dp) then
        call findVarValues(i, UArray)
        write(histFile, formatHist) (real(UArray(j), dp)*UBinSize(j), j=1,nBiasVariables), HistStorage(i)
      endif
    enddo
    close(histFile)

    ! Step 3.5: Write extensive histogram data if enabled
    if (WHAM_ExtensiveOutput) then
      ! Write numerator to WHAM_Hist_Numerator.txt
      open(unit=histFile, file="WHAM_Hist_Numerator.txt")
      do i = 1, umbrellaLimit
        if (HistStorage(i) > 0.0_dp) then
          call findVarValues(i, UArray)
          write(histFile, formatPrecise) (real(UArray(j), dp)*UBinSize(j), j=1,nBiasVariables), &
                                         WHAM_Numerator(i)
        endif
      enddo
      close(histFile)

      ! Write denominator to WHAM_Hist_Denominator.txt
      open(unit=histFile, file="WHAM_Hist_Denominator.txt")
      do i = 1, umbrellaLimit
        if (HistStorage(i) > 0.0_dp) then
          call findVarValues(i, UArray)
          write(histFile, formatPrecise) (real(UArray(j), dp)*UBinSize(j), j=1,nBiasVariables), &
                                         (WHAM_Denominator(i, j), j=1,nCurWhamItter)
        endif
      enddo
      close(histFile)
    endif
  endif
end subroutine WHAM_Finalize
!=========================================================================     
      end module
