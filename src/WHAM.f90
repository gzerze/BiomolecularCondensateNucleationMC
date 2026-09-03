!=========================================================================
      subroutine WHAM_Initialize
      use ParallelVar
      use SimParameters
      use WHAM_Module
      implicit none
      integer :: AllocateStatus


      if(myid .eq. 0) then
        allocate(WHAM_Numerator(1:umbrellaLimit), STAT = AllocateStatus)
        allocate(WHAM_Denominator(1:umbrellaLimit,1:nWhamItter+1), STAT = AllocateStatus)
        allocate(HistStorage(1:umbrellaLimit), STAT = AllocateStatus)
        allocate(BiasStorage(1:umbrellaLimit,1:nWhamItter+1), STAT = AllocateStatus)
        allocate(FreeEnergyEst(1:umbrellaLimit), STAT = AllocateStatus)
        allocate(ProbArray(1:umbrellaLimit), STAT = AllocateStatus)

!        allocate(NewBias(1:umbrellaLimit), STAT = AllocateStatus)
      
        write(nout,*) "Allocated WHAM Variables"
      
  
        WHAM_Numerator = 0E0
        WHAM_Denominator = 0E0
        HistStorage = 0E0
        BiasStorage = 0E0
        nCurWhamItter = 1
!        tolLimit = 1d-2
        open(unit = 97, file="WHAM_Potential.incomp")
        open(unit = 98, file="WHAM_Mid_DG.incomp")
      endif

      allocate(TempHist(1:umbrellaLimit), STAT = AllocateStatus)      
      allocate(NewBias(1:umbrellaLimit), STAT = AllocateStatus)
      NewBias = 0E0
      TempHist = 0E0
      end subroutine

!=========================================================================
!     This subroutine periodically adjusts the Umbrella Sampling Bias
!     by collecting histogram data from across 
      subroutine WHAM_AdjustHist
      use ParallelVar     
      use SimParameters
      use WHAM_Module
      implicit none
      include 'mpif.h' 

      interface
       subroutine WHAM_CurveSmoothing(NewBias,HistStorage)
         use SimParameters
         use UmbrellaFunctions
         implicit none
         real(dp), intent(inout) :: NewBias(:)
         real(dp), intent(in) :: HistStorage(:)
       end subroutine
      end interface

      interface
       subroutine WHAM_CurveSmoothing2(NewBias)
         use SimParameters
         use UmbrellaFunctions
         implicit none
         real(dp), intent(inout) :: NewBias(:)
       end subroutine
      end interface

      integer :: arraySize, i, j, cnt, maxbin, maxbin2
      integer :: NArray(1:nMolTypes)
      real(dp) :: norm, ratio, maxBias, denomSum
      real(dp) :: F_Estimate(1:nWhamItter),F_Old(1:nWhamItter), fSum     
      real(dp) :: tol, refBias, perError, probNorm


      write(nout,*) "Halting for WHAM"
      call MPI_BARRIER(MPI_COMM_WORLD, ierror) 
      
!      This block condences the histogram data from all the different processors
!      into one collective array on the root (myid = 0) processor.        
      arraySize = size(NHist)     
      if(myid .eq. 0) then
        TempHist = 0E0
      endif
      call MPI_REDUCE(NHist, TempHist, arraySize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)       

    
      if(myid .eq. 0) then
!        This block calculates the terms needed to 
        norm = sum(TempHist)
        do i = 1, umbrellaLimit
          BiasStorage(i,nCurWhamItter) = NBias(i) 
          if(TempHist(i) .ne. 0E0) then
            WHAM_Numerator(i) = WHAM_Numerator(i) + TempHist(i)*TempHist(i)/norm
            WHAM_Denominator(i, nCurWhamItter) = TempHist(i)*exp(NBias(i))
!            WHAM_Numerator(i) = WHAM_Numerator(i) + TempHist(i)
!            WHAM_Denominator(i, nCurWhamItter) = norm*exp(NBias(i))
            HistStorage(i) = HistStorage(i) + TempHist(i)
          endif
        enddo 

!        This block solves for the free energy terms required by WHAM.  This is done
!        self-consistently.
        ProbArray = 0E0
        F_Estimate = 0E0
        tol = tolLimit + 1E0
        cnt = 0
        do while(tol .gt. tolLimit)
          cnt = cnt + 1
!          Infinite Loop protection
          if(cnt .gt. maxSelfConsist) then 
            write(35,*) "Self Consistent Limit Hit"
            exit
          endif

          ProbArray = 0E0
          do j = 1, nCurWhamItter
            F_Old(j) = F_Estimate(j)
          enddo
!            If bin #i has been sampled at any point in the simulation, estimate the unbiased probability
!            based on the current guess value for F
          do i = 1, umbrellaLimit
            if(WHAM_Numerator(i) .ne. 0E0) then
              denomSum = 0E0
!              maxBias = minval(F_Estimate)
              do j = 1, nCurWhamItter
                if(WHAM_Denominator(i,j) .gt. 0E0) then
                  denomSum = denomSum + WHAM_Denominator(i,j)*exp(-F_Estimate(j))
                endif
              enddo
              if(denomSum .ne. 0E0) then
                ProbArray(i) = WHAM_Numerator(i)/denomSum
              endif
            else
              ProbArray(i) = 0E0
            endif
          enddo 

          norm = sum(ProbArray)
          do i = 1, umbrellaLimit
            ProbArray(i) = ProbArray(i)/norm
          enddo 
!          Once all the unbiased probabilities have been estimated, use these unbiased estimates
!          to calculate a new estimate for F
          do j = 1, nCurWhamItter
            fSum = 0E0
!            maxBias = maxval(BiasStorage(:,j))
            do i = 1, umbrellaLimit
              if(ProbArray(i) .ne. 0E0) then
                fSum = fSum + ProbArray(i)*exp(BiasStorage(i,j))
!                fSum = fSum + ProbArray(i)*exp(BiasStorage(i,j) - maxBias)
              endif
            enddo
!            F_Estimate(j) = log(fSum)
            F_Estimate(j) = log(fSum)
            F_Estimate(j) = (F_Estimate(j) + F_Old(j))*0.5E0
         enddo 
!         Calculate the average change in F from the previous estimate and determine 
!         if there has been a significant change to the F values.
         tol = 0E0
         do j = 1, nCurWhamItter
           tol = tol + abs(F_Estimate(j) - F_Old(j))
!           tol = max(tol, abs(F_Estimate(j) - F_Old(j)) )
         enddo
!         tol = tol/dble(nCurWhamItter)
       enddo

!        Using the new estimates for the unbiased probability, calculate the free energy of nucleation
!        and modify the umbrella sampling bias to
        NewBias = 0E0
        maxbin = maxloc(HistStorage,1)
        if(mod(nCurWhamItter,whamEstInterval) .eq. 0) then
!          maxbin = maxloc(HistStorage,1)
          do i = 1, umbrellaLimit
            if(ProbArray(i) .gt. 0E0) then
              FreeEnergyEst(i) = -log(ProbArray(i)/ProbArray(maxbin))
              NewBias(i) = FreeEnergyEst(i)
            endif
          enddo
!          maxBias = maxval(NewBias)
          maxBias = -huge(dp)
          do i = 1, umbrellaLimit
            if(ProbArray(i) .gt. 0E0) then
              if(maxBias .lt. FreeEnergyEst(i)) then
                maxBias = FreeEnergyEst(i)
              endif
            endif
          enddo
          do i = 1, umbrellaLimit
            if(ProbArray(i) .le. 0E0) then
              NewBias(i) = maxBias + 1E0
              FreeEnergyEst(i) = maxBias + 1E0
            endif
          enddo
!          write(*,*) "WHAM Potenial Used"
!          call WHAM_CurveSmoothing(NewBias, HistStorage)
        else
          maxbin2 = maxloc(TempHist,1)
          do i = 1, umbrellaLimit
            if(ProbArray(i) .gt. 0E0) then
              FreeEnergyEst(i) = -log(ProbArray(i)/ProbArray(maxbin))
            endif
            if(TempHist(i) .gt. 0E0) then
              NewBias(i) = NBias(i) - NBias(maxbin2) - log(TempHist(i)/TempHist(maxbin2))
            endif
          enddo
!          maxBias = NBias(maxbin2)
          do i = 1, umbrellaLimit
            if(TempHist(i) .le. 0E0) then
!              NBias(i) = NBias(i) + log(3E0)
!              if(ProbArray(i) .gt. 0E0) then
!                NewBias(i) = NBias(i) - NBias(maxbin2) + log(10E0)
!              else
!                NewBias(i) = NBias(i) - NBias(maxbin2) + nCurWhamItter*log(10E0)
                NewBias(i) = NBias(i) - NBias(maxbin2) + log(TempHist(maxbin2))
!                write(*,*) i, Newbias(i)
!              endif
            endif
          enddo
!          call WHAM_CurveSmoothing(NewBias, TempHist)
!          call WHAM_CurveSmoothing2(NewBias)
        endif
!        Rescale the pontential such that the reference free energy is set to 0
        refBias = NewBias(refBin)
        do i = 1, umbrellaLimit
          NewBias(i) = NewBias(i) - refBias
        enddo
        refBias = FreeEnergyEst(refBin)
        do i = 1, umbrellaLimit
          FreeEnergyEst(i) = FreeEnergyEst(i) - refBias
        enddo
        
        call WHAM_MidSimOutput
      endif
!     End of processor 0 only block

      call MPI_BARRIER(MPI_COMM_WORLD, ierror) 
!      Distribute the new free energy estimate to all threads so that they can continue the simulation
!      with the new free energy. 
      arraySize = size(NewBias)      
      call MPI_BCast(NewBias, arraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror) 

      do i = 1, umbrellaLimit
        NBias(i) = NewBias(i)
        NHist(i) = 0E0
      enddo 

      nCurWhamItter = nCurWhamItter + 1
       
      end subroutine
!==================================================================================
      subroutine WHAM_MidSimOutput
      use SimParameters
      use ParallelVar
      use WHAM_Module
      use SwapBoundary
      implicit none
      logical :: goToNext
      integer :: arraySize
      integer :: i,j
      integer :: NArray(1:nMolTypes), NDiff(1:nMolTypes)
      real(dp) :: refBias

!        This block exports the calculated free energy to a file
      rewind(97)
      NArray = 0
      NDiff = 0
      do i = 1, umbrellaLimit
        if(nMolTypes .gt. 1) then
          do j = 1,nMolTypes-1
            if(NArray(nMolTypes - j + 1) .gt. NMAX(nMolTypes - j + 1))then
             NArray(nMolTypes - j + 1) = 0
             NArray(nMolTypes - j) = NArray(nMolTypes - j) + 1          
            endif
          enddo
        endif
        goToNext = BoundaryFunction(NArray, NDiff)
        if(goToNext) then
          NArray(nMolTypes) = NArray(nMolTypes) + 1
          cycle
        endif
        if(i .ne. 1) then
          write(97, *) (NArray(j),j=1,nMolTypes), NBias(i)
        endif
        NArray(nMolTypes) = NArray(nMolTypes) + 1
      enddo
      flush(97)
 
      rewind(98)
      NArray = 0
      do i = 1, umbrellaLimit
        if(nMolTypes .gt. 1) then
          do j = 1,nMolTypes-1
            if(NArray(nMolTypes - j + 1) .gt. NMAX(nMolTypes - j + 1))then
             NArray(nMolTypes - j + 1) = 0
             NArray(nMolTypes - j) = NArray(nMolTypes - j) + 1
            endif
          enddo
        endif
        goToNext = BoundaryFunction(NArray, NDiff)
        if(goToNext) then
          NArray(nMolTypes) = NArray(nMolTypes) + 1
          cycle
        endif
        if(i .ne. 1) then
          if(ProbArray(i) .ne. 0E0 ) then
            write(98, *) (NArray(j),j=1,nMolTypes), FreeEnergyEst(i)
          endif
        endif
        NArray(nMolTypes) = NArray(nMolTypes) + 1
      enddo
      flush(98)

           
      end subroutine
!==================================================================================
      subroutine WHAM_Finalize
      use SimParameters
      use ParallelVar
      use WHAM_Module
      use SwapBoundary
      implicit none
      include 'mpif.h' 
      logical :: goToNext
      integer :: arraySize
      integer :: i,j
      integer :: NArray(1:nMolTypes)
      integer :: nDiff(1:nMolTypes)
      real(dp) :: norm, ratio, maxBias, refBias, probNorm

!      call WHAM_AdjustHist

      if(myid .eq. 0) then
!        This block exports the calculated free energy to a file
        open(unit = 92, file="WHAM_DG_Output.txt")
        NArray = 0
        nDiff = 0
        refBias = FreeEnergyEst(refBin)
        do i = 1, umbrellaLimit
          if(nMolTypes .gt. 1) then
            do j = 1,nMolTypes-1
              if(NArray(nMolTypes - j + 1) .gt. NMAX(nMolTypes - j + 1))then
               NArray(nMolTypes - j + 1) = 0
               NArray(nMolTypes - j) = NArray(nMolTypes - j) + 1          
              endif
            enddo
          endif
          goToNext = BoundaryFunction(NArray, NDiff)
          if(goToNext) then
            NArray(nMolTypes) = NArray(nMolTypes) + 1
            cycle
          endif
          if(ProbArray(i) .gt. 0E0) then
            write(92, *) (NArray(j),j=1,nMolTypes), FreeEnergyEst(i)
!          else
!            if(i .ne. 1) then
!              write(92, *) (NArray(j),j=1,nMolTypes), NBias(i)
!            endif
          endif
          NArray(nMolTypes) = NArray(nMolTypes) + 1
        enddo
        close(92)

!        This block exports the final bias potential to a text file. 
        open(unit = 92, file="WHAM_FinalBias.txt")
        NArray = 0
        refBias = NBias(refBin)
        do i = 1, umbrellaLimit
          if(nMolTypes .gt. 1) then
            do j = 1,nMolTypes-1
              if(NArray(nMolTypes - j + 1) .gt. NMAX(nMolTypes - j + 1))then
               NArray(nMolTypes - j + 1) = 0
               NArray(nMolTypes - j) = NArray(nMolTypes - j) + 1          
              endif
            enddo
          endif
          goToNext = BoundaryFunction(NArray, NDiff)
          if(goToNext) then
            NArray(nMolTypes) = NArray(nMolTypes) + 1
            cycle
          endif
          if(i .ne. 1) then
            write(92, *) (NArray(j),j=1,nMolTypes), NBias(i)
          endif
          NArray(nMolTypes) = NArray(nMolTypes) + 1
        enddo
        close(92)


!        This block exports the calculated probabilities to a file
        open(unit = 92, file="WHAM_Probabilities.txt")
        NArray = 0
        probNorm = sum(ProbArray)
        do i = 1, umbrellaLimit
          if(nMolTypes .gt. 1) then
            do j = 1,nMolTypes-1
              if(NArray(nMolTypes - j + 1) .gt. NMAX(nMolTypes - j + 1))then
                NArray(nMolTypes - j + 1) = 0
                NArray(nMolTypes - j) = NArray(nMolTypes - j) + 1
              endif
            enddo
          endif
          goToNext = BoundaryFunction(NArray, NDiff)
          if(goToNext) then
            NArray(nMolTypes) = NArray(nMolTypes) + 1
            cycle
          endif
          if(ProbArray(i) .gt. 0E0) then
            write(92, *) (NArray(j),j=1,nMolTypes), ProbArray(i)/probNorm
          endif
          NArray(nMolTypes) = NArray(nMolTypes) + 1
        enddo
        close(92)



        open(unit=36, file="WHAM_OverallHist.txt")
        NArray = 0
        do i = 1, umbrellaLimit
          if(nMolTypes .gt. 1) then
            do j = 1,nMolTypes-1
              if(NArray(nMolTypes - j + 1) .gt. NMAX(nMolTypes - j + 1))then
                NArray(nMolTypes - j + 1) = 0
                NArray(nMolTypes - j) = NArray(nMolTypes - j) + 1          
              endif
            enddo
          endif
          goToNext = BoundaryFunction(NArray, NDiff)
          if(goToNext) then
            NArray(nMolTypes) = NArray(nMolTypes) + 1
            cycle
          endif
          if(HistStorage(i) .gt. 0E0) then
            write(36, *) (NArray(j),j=1,nMolTypes), HistStorage(i)
          endif
          NArray(nMolTypes) = NArray(nMolTypes) + 1
        enddo
        close(36)

        if(WHAM_ExtensiveOutput) then
          open(unit=36, file="WHAM_Hist_Numerator.txt")
          NArray = 0
          do i = 1, umbrellaLimit
            if(nMolTypes .gt. 1) then
              do j = 1,nMolTypes-1
                if(NArray(nMolTypes - j + 1) .gt. NMAX(nMolTypes - j + 1))then
                  NArray(nMolTypes - j + 1) = 0
                  NArray(nMolTypes - j) = NArray(nMolTypes - j) + 1          
                endif
              enddo
            endif
            if(HistStorage(i) .gt. 0E0) then
              write(36, *) (NArray(j),j=1,nMolTypes), WHAM_Numerator(i)
            endif
            NArray(nMolTypes) = NArray(nMolTypes) + 1
          enddo
          close(36)

          open(unit=36, file="WHAM_Hist_Denominator.txt")
          NArray = 0
          do i = 1, umbrellaLimit
            if(nMolTypes .gt. 1) then
              do j = 1,nMolTypes-1
                if(NArray(nMolTypes - j + 1) .gt. NMAX(nMolTypes - j + 1))then
                  NArray(nMolTypes - j + 1) = 0
                  NArray(nMolTypes - j) = NArray(nMolTypes - j) + 1          
                endif
              enddo
            endif
            if(HistStorage(i) .gt. 0E0) then
              write(36, *) (NArray(j),j=1,nMolTypes), (WHAM_Denominator(i,j), j=1,nCurWhamItter)
            endif
            NArray(nMolTypes) = NArray(nMolTypes) + 1
          enddo
          close(36)
        endif
      endif

             
      end subroutine
!==================================================================================
      subroutine WHAM_CurveSmoothing(NewBias, HistStorage)
      use SimParameters
      use UmbrellaFunctions
      implicit none
      include 'mpif.h' 
      real(dp), intent(inout) :: NewBias(:)
      real(dp), intent(in) :: HistStorage(:)
      real(dp), allocatable :: TempBias(:), weightArray(:)
      real(dp), parameter :: weightLimit = 0.5E0
      real(dp), parameter :: gaussPara = log(1d3)
      logical :: arrayCycle
      integer :: arraySize, cnt
      integer :: i,j,k
      integer :: NArray(1:nMolTypes), NDiff(1:nMolTypes)
      integer :: nNeighCells, bIndx
      real(dp) :: curBias, curWeight, sumWeight, maxHist
      real(dp) :: sumHist

      nNeighCells = 3**nMolTypes 
      arraySize = size(NewBias)
!      gaussPara = 1E0/dble(ncycle*ncycle2)
      allocate( TempBias(1:arraySize) )
      allocate( weightArray(1:arraySize) )
  
      sumHist = sum(HistStorage)
      maxHist = maxval(HistStorage)
      do i = 1, arraySize
!        if(exp(-1d-2*(1E0-HistStorage(i)/maxHist)) .gt. 0.01E0) then
!          weightArray(i) = weightLimit*abs(1E0-HistStorage(i)/maxHist)
          weightArray(i) = weightLimit*exp(-gaussPara*((HistStorage(i)-maxHist)/sumHist)**2)
!        else
!          weightArray(i) = weightLimit*(exp(-1d-2)-0.01E0)
!        endif
!        write(*,*) i, weightArray(i), HistStorage(i),(HistStorage(i)-maxHist)/sumHist
      enddo

      NArray = 0
      do i = 2, arraySize
        sumWeight = 1E0
        curBias = NewBias(i)
        NArray(nMolTypes) = NArray(nMolTypes) + 1
        if(nMolTypes .gt. 1) then
          do j = 1, nMolTypes - 1
            if(NArray(nMolTypes - j + 1) .gt. NMAX(nMolTypes - j + 1))then
              NArray(nMolTypes - j + 1) = 0
              NArray(nMolTypes - j) = NArray(nMolTypes - j) + 1          
            endif
          enddo
        endif

        NDiff = -1
        NDiff(nMolTypes) = -2
        cnt = 0
        do k = 1, nNeighCells
          NDiff(nMolTypes) = NDiff(nMolTypes) + 1
          if(nMolTypes .gt. 1) then
            do j = 1, nMolTypes - 1
              if(NDiff(nMolTypes - j + 1) .gt. 1)then
                NDiff(nMolTypes - j + 1) = -1
                NDiff(nMolTypes - j) = NDiff(nMolTypes - j) + 1          
              endif
            enddo
          endif
          arrayCycle = .false.
          if( all(NDiff .eq. 0) ) then
            arrayCycle = .true.
          endif
          do j = 1, nMolTypes
            if(NArray(j) + NDiff(j) .lt. NMIN(j)) then
              arrayCycle = .true.
            endif
            if(NArray(j) + NDiff(j) .gt. NMAX(j))  then
              arrayCycle = .true.
            endif
          enddo
          if(arrayCycle) then
            cycle
          endif
          cnt = cnt + 1
          bIndx = getNewBiasIndex(NArray, NMAX, NDiff)
          curWeight = (weightLimit - weightArray(i)) * weightArray(bIndx)
          curWeight = sqrt(curWeight)
          sumWeight = sumWeight + curWeight
          curBias = curBias + curWeight * NewBias(bIndx)
        enddo
!        if(sumWeight .ne. 0E0) then
        curBias = curBias/(sumWeight)
        TempBias(i) = curBias
!        write(*,*) i, NArray, NewBias(i), TempBias(i)
      enddo
    
      NewBias = TempBias
      deallocate(TempBias)             

      end subroutine
!=========================================================================     
