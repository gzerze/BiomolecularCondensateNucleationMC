!======================================================
     module AnalysisMain

!     To add new analysis functions to this code, insert the module for the
!     new function set here.
     use RadialDensity
     use RadialDistribution
     use SimpleDistPair
     use MiscelaniousVars

     integer, parameter :: maxLineLen = 500  

     interface 
       subroutine TrialFunction(disp)
        use CoordinateTypes
        type(Displacement), intent(in) :: disp(:)
       end subroutine
     end interface

     interface 
       subroutine UmbrellaLoader(iUmbrella,varIndx, biasVar, biasVarNew, outputFormat, &
                                 iDisp, DispUmbrella, iSwapIn, SwapInUmbrella, iSwapOut, SwapOutUmbrella)
         use MiscelaniousVars
         use UmbrellaTypes
         implicit none 
         integer, intent(in) :: iUmbrella, varIndx
         integer, intent(inout) :: iDisp, iSwapIn, iSwapOut
         
         type(BiasVariablePointer), intent(inout) :: biasvar(:)
         type(BiasVariablePointer), intent(inout) :: biasvarnew(:)
         type(DispUmbrellaArray), intent(inout)  :: DispUmbrella(:)
         type(SwapInUmbrellaArray), intent(inout)  :: SwapInUmbrella(:)
         type(SwapOutUmbrellaArray), intent(inout)  :: SwapOutUmbrella(:)
         character(len=10), intent(inout) :: outputFormat(:)
       end subroutine
     end interface

     type TrialFunctionArray
       procedure(TrialFunction), pointer, nopass :: func
     end type

     type AnalysisFunctionArray
       procedure(), pointer, nopass :: func
     end type

     type UmbrellaLoaderArray
       procedure(UmbrellaLoader), pointer, nopass :: func
     end type



     private
     logical :: useAnalysis = .false.
     integer :: nAnalysisVar = 0
     integer :: nPostMove, nOutput
!     type(TrialFunctionArray), allocatable :: TrialArray(:)
     type(AnalysisFunctionArray), allocatable :: postMoveArray(:)
     type(AnalysisFunctionArray), allocatable :: outputArray(:)

     integer, allocatable :: internalIndx(:)
     type(UmbrellaLoaderArray), allocatable :: loadUmbArray(:)
!     type(TrialFunctionArray), allocatable :: UpdateArray(:)

     public :: useAnalysis, PostMoveAnalysis, OutputAnalysis
     public :: ScriptAnalysisInput, internalIndx, loadUmbArray, nAnalysisVar

!======================================================
     contains
!======================================================
     subroutine ScriptAnalysisInput(inputLines)
      use MiscelaniousVars
      use SimpleDistPair, only: nDistPair, UmbrellaVar_DistPair
      use Q3Functions, only: CalcQ3, Initialize_q3, useQ3, q3Dist, q3DistSq, &
                             UmbrellaVar_Q3
      use Q4Functions, only: CalcQ4, Initialize_q4, useQ4, q4Dist, q4DistSq, &
                             UmbrellaVar_Q4
      use Q6Functions, only: CalcQ6, Initialize_q6, useQ6, q6Dist, q6DistSq, &
                             UmbrellaVar_Q6
      implicit none
      character(len=maxLineLen), intent(in) :: inputLines(:)
      integer :: nLines
      integer :: iAnalysis
      integer :: nBins
      integer :: iRadial, iDistPair, iRadDens, iPostMove, iOutput
      integer :: type1, type2, mol1, mol2, atom1, atom2
      real(dp) :: binSize, realVar
      character(len=30) :: labelField 
      character(len=30) :: analysisName, fileName

      nLines = size(inputLines)
      nAnalysisVar = nLines - 2
      if(nAnalysisVar .lt. 0) then
        write(*,*) "ERROR! The user has specified an invalid number of Analysis Variables"
        write(*,*) labelField, nAnalysisVar
        stop
      endif
      useAnalysis = .true.

!      Begin by counting how many entries are required for each function array.
!      These will be used in the next step to allocate the arrays. 
      nPostMove = 0
      nOutput = 0
      nDistPair = 0 
      nRadialDist = 0
      nRadialDens = 0

      allocate( internalIndx(1:nAnalysisVar) ) 
      allocate( loadUmbArray(1:nAnalysisVar) )

      do iAnalysis = 1, nAnalysisVar
        loadUmbArray(iAnalysis)%func => Null()
      enddo


      internalIndx = 0
      do iAnalysis = 1, nAnalysisVar
        read(inputLines(iAnalysis+1), *) analysisName

        select case( trim(adjustl(analysisName)) )
        case("radialdistribution")
          nRadialDist = nRadialDist + 1
          internalIndx(iAnalysis) = nRadialDist
          if(nRadialDist .eq. 1) then
            nPostMove = nPostMove + 1
            nOutput = nOutput + 1
          endif
        case("pairdist")
          nDistPair = nDistPair + 1
          internalIndx(iAnalysis) = nDistPair
          if(nDistPair .eq. 1) then
            nPostMove = nPostMove + 1
          endif

        case("q3")
          if(useQ3 .eqv. .false.) then
            useQ3 = .true.
            nPostMove = nPostMove + 1
            internalIndx(iAnalysis) = 1
          else
            stop "ERROR! The Q4 analysis function has been defined more than once in the input script."
          endif

        case("q4")
          if(useQ4 .eqv. .false.) then
            useQ4 = .true.
            nPostMove = nPostMove + 1
            internalIndx(iAnalysis) = 1
          else
            stop "ERROR! The Q4 analysis function has been defined more than once in the input script."
          endif

        case("q6")
          if(useQ6 .eqv. .false.) then
            useQ6 = .true.
            nPostMove = nPostMove + 1
            internalIndx(iAnalysis) = 1
          else
            stop "ERROR! The Q6 analysis function has been defined more than once in the input script."
          endif

        case("radialdensity")
          nRadialDens = nRadialDens + 1
          internalIndx(iAnalysis) = nRadialDens
          if(nRadialDens .eq. 1) then
            nPostMove = nPostMove + 1
            nOutput = nOutput + 1
          endif
        case default
          write(*,*) "ERROR! Invalid variable type specified in input file"
          write(*,*) analysisName
          stop
        end select
      enddo

!      Now that we know how much memory is required, allocate all relevent arrays
      if(nPostMove .ne. 0) then
        allocate( postMoveArray(1:nPostMove) )
        do iPostMove = 1, nPostMove
          postMoveArray(iPostMove)%func => null()
        enddo
      endif
      if(nOutput .ne. 0) then
        allocate( outputArray(1:nOutput) )
        do iOutPut = 1, nOutput
          outputArray(iOutPut)%func => null()
        enddo
      endif 



      call Initialize_RadialDist
      call Initialize_DistPair
      call Initialize_RadialDens
      call Initialize_Q3
      call Initialize_Q4
      call Initialize_Q6
      call AllocateMiscArrays

      iOutPut = 0
      iPostMove = 0
      iRadial = 0
      iDistPair = 0
      iRadDens = 0
      do iAnalysis = 1, nAnalysisVar
        read(inputLines(iAnalysis+1), *) analysisName

        select case( trim(adjustl(analysisName)) )
        case("radialdistribution")
          iRadial = iRadial + 1
          read(inputLines(iAnalysis+1), *)  analysisName, type1, atom1, type2, atom2, binSize, nBins, fileName

          call SetRadialParameters(iRadial, type1, type2, atom1, atom2)
          call SetRadialHist(iRadial, binSize, nBins, fileName)

          if(iRadial .eq. 1) then
            iPostMove = iPostMove + 1
            postMoveArray(iPostMove)%func => Calc_RadialDist
            iOutPut = iOutPut + 1
            outputArray(iOutPut)%func => Output_RadialDist
          endif 
        case("pairdist")
          iDistPair = iDistPair + 1
          read(inputLines(iAnalysis+1), *)  analysisName, type1, mol1, atom1, type2, mol2, atom2
          call SetPairVariables(iDistPair, Type1, Mol1, Atom1, Type2, Mol2, Atom2, iAnalysis)
          if(iDistPair .eq. 1) then
            iPostMove = iPostMove + 1
            postMoveArray(iPostMove)%func => CalcDistPairs
          endif 
          loadUmbArray(iAnalysis)%func => UmbrellaVar_DistPair
        case("radialdensity")
          read(inputLines(iAnalysis+1), *)  analysisName, type1, binSize, nBins, fileName
          iRadDens = iRadDens + 1
          call SetDensityParameters(iRadDens, type1)
          call SetDensityHist(iRadDens, binSize, nBins, fileName)
          if(iRadDens .eq. 1) then
            iPostMove = iPostMove + 1
            iOutPut = iOutPut + 1
            postMoveArray(iPostMove)%func => Calc_RadialDensity
            outputArray(iOutPut)%func => Output_RadialDensity
          endif
        case("q3")
          read(inputLines(iAnalysis+1), *)  analysisName, realVar
          iPostMove = iPostMove + 1
          q3Dist = realVar
          q3DistSq = q3Dist*q3Dist
          postMoveArray(iPostMove)%func => CalcQ3
          loadUmbArray(iAnalysis)%func => UmbrellaVar_Q3
        case("q4")
          read(inputLines(iAnalysis+1), *)  analysisName, realVar
          iPostMove = iPostMove + 1
          q4Dist = realVar
          q4DistSq = q4Dist*q4Dist
          postMoveArray(iPostMove)%func => CalcQ4
          loadUmbArray(iAnalysis)%func => UmbrellaVar_Q4
        case("q6")
          read(inputLines(iAnalysis+1), *)  analysisName, realVar
          iPostMove = iPostMove + 1
          q6Dist = realVar
          q6DistSq = q6Dist*q6Dist
          postMoveArray(iPostMove)%func => CalcQ6
          loadUmbArray(iAnalysis)%func => UmbrellaVar_Q6

        case default
          write(*,*) "ERROR! Invalid variable type specified in input file"
          write(*,*) analysisName
          stop
        end select
      enddo

      call AllocateHistBins

     end subroutine
!======================================================
!     subroutine TrialPositionAnalysis(disp)
!     use CoordinateTypes
!     implicit none
!     type(Displacement), intent(in) :: disp(:)
!     integer :: iTrialVar
! 
!     do iTrialVar = 1, nTrialVar
!       call TrialArray(iTrialVar)%func(disp)
!     enddo
! 
!     end subroutine
!======================================================
     subroutine PostMoveAnalysis
     implicit none
     integer :: iPostMove

     do iPostMove = 1, nPostMove
       call postMoveArray(iPostMove)%func
     enddo

     end subroutine
!======================================================
     subroutine OutputAnalysis
     use ParallelVar
     implicit none
     integer :: iOutput

     do iOutput = 1, nOutput
       call outputArray(iOutput)%func
     enddo

     end subroutine
!======================================================
     end module
!======================================================
