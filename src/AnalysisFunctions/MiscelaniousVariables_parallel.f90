!===========================================================
!   The purpose of this module is to provide a set of general purpose
!   memory variables that can be used in the creation of analysis funcitons
!   such as radial distribution functions.  
    module MiscelaniousVars
    use VarPrecision

    type Histograms
      character(len=10) :: histName
      character(len=50) :: fileName
      integer :: nBins
      real(dp) :: binSize, sizeInv
      real(dp), allocatable :: binCount(:)
    end type

    integer :: nIntegers = 0    
    integer :: nGeometry = 0
    integer :: nHistArrays = 0

    integer, allocatable, target  :: miscInt(:)
    real(dp), allocatable, target  :: miscCoord(:)
    real(dp), allocatable, target  :: miscCoord_New(:)
    type(Histograms), allocatable :: miscHist(:)

    contains
!===========================================================
    subroutine ReserveSpace_Integers(nReserve, startIndx, endIndx)
    implicit none
    integer, intent(in) :: nReserve
    integer, intent(out) :: startIndx, endIndx

    startIndx = nIntegers + 1
    nIntegers = nIntegers + nReserve
    endIndx = nIntegers

    end subroutine
!===========================================================
    subroutine ReserveSpace_Histograms(nReserve, startIndx, endIndx)
    implicit none
    integer, intent(in) :: nReserve
    integer, intent(out) :: startIndx, endIndx

    startIndx = nHistArrays + 1
    nHistArrays = nHistArrays + nReserve
    endIndx = nHistArrays

!    write(*,*) "Spaced Reserved:", nHistArrays, startIndx, endIndx

    end subroutine

!===========================================================
    subroutine ReserveSpace_Coord(nReserve, startIndx, endIndx)
    implicit none
    integer, intent(in) :: nReserve
    integer, intent(out) :: startIndx, endIndx

    startIndx = nGeometry + 1
    nGeometry = nGeometry + nReserve
    endIndx = nGeometry

    end subroutine
!===========================================================
    subroutine AllocateMiscArrays
    implicit none
    integer :: AllocationStat

    if(nGeometry .ne. 0) then
      allocate(miscCoord(1:nGeometry), Stat = AllocationStat)
      allocate(miscCoord_New(1:nGeometry), Stat = AllocationStat)
    endif

    if(nHistArrays .ne. 0) then
      allocate(miscHist(1:nHistArrays), Stat = AllocationStat)
    endif

    end subroutine
!===========================================================
    subroutine AllocateHistBins
    implicit none
    integer :: iHist, iBin, AllocationStat, nBins

    do iHist = 1, nHistArrays
      nBins = miscHist(iHist)%nBins
      allocate(miscHist(iHist)%binCount(0:nBins+1), Stat = AllocationStat)
      do iBin = 0, nBins+1
        miscHist(iHist)%binCount(iBin) = 0E0
      enddo
    enddo

    end subroutine
!===========================================================
    subroutine CollectHistograms
    use ParallelVar
    use MPI
    implicit none
    integer :: iHist, iBin, nBins 
    real(dp), allocatable :: TempHist(:)

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
!    write(*,*) "Collecting General Histogram Data"

    do iHist = 1, nHistArrays
      nBins = miscHist(iHist)%nBins
      if(myid .eq. 0) then
        allocate( TempHist(0:nBins+1) )
        TempHist = 0E0
      endif
      call MPI_REDUCE(miscHist(iHist)%binCount, TempHist, nBins+2, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)  
      if(myid .eq. 0) then
        do iBin= 0, nBins+1
          miscHist(iHist)%binCount(iBin) = TempHist(iBin)
        enddo
        deallocate(TempHist)
      endif
      call MPI_BARRIER(MPI_COMM_WORLD, ierror) 
    enddo

    end subroutine
!===========================================================
    end module
!===========================================================
