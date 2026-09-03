!================================================================
      subroutine BlankUmbrellaBias
      use SimParameters
      use WHAM_Module
      implicit none
      integer :: AllocateStatus
      integer :: i 
      
      umbrellaLimit = 1
      do i=1,nMolTypes
        umbrellaLimit = umbrellaLimit*(NMAX(i) + 1)
      enddo
        
      allocate(NBias(1:umbrellaLimit), STAT = AllocateStatus)
      allocate(NHist(1:umbrellaLimit), STAT = AllocateStatus)
      allocate(E_Avg(1:umbrellaLimit), STAT = AllocateStatus)  
     
      NBias = 0d0
      NHist = 0d0



      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      end subroutine
!==========================================================================================
      subroutine AllocateUmbrellaBias(fileName)
      use UmbrellaFunctions
      use SimParameters
      use WHAM_Module
      implicit none
      integer :: i,j
      integer :: AllocateStatus
      integer,allocatable :: arrayIndx(:)   
      integer :: curIndx
      character(len=30) :: fileName   
      real(dp) :: curValue
        
      umbrellaLimit = 1
      do i=1,nMolTypes
        umbrellaLimit = umbrellaLimit*(NMAX(i) + 1)
      enddo
        
      allocate(NBias(1:umbrellaLimit), STAT = AllocateStatus)
      allocate(NHist(1:umbrellaLimit), STAT = AllocateStatus)
      allocate(E_Avg(1:umbrellaLimit), STAT = AllocateStatus)     
      allocate(arrayIndx(1:nMolTypes), STAT = AllocateStatus)
       
      NHist = 0d0
      NBias = 0d0
      open( unit = 25, file = trim( adjustl(fileName) ), status='OLD')

      write(35,*) ""
      do i = 1,umbrellaLimit
!      do i = 1,nint(1d7)
         read(25,*,end=30) (arrayIndx(j), j=1,nMolTypes), curValue
!         write(35,*) (arrayIndx(j), j=1,nMolTypes), curValue
         curIndx = getBiasIndex(arrayIndx,NMAX)
         if(curIndx .le. umbrellaLimit) then
           if(curIndx .gt. 0) then
             NBias(curIndx) = curValue
           endif
         endif
      enddo
30    close(25)

!      if(useWHAM) then
!        call WHAM_Initialize
!      endif
       
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      end subroutine
!==========================================================================================
      subroutine NHistAdd(E_T)
      use SimParameters
      implicit none
      integer :: i, sizeN
      integer :: curIndx,maxIndx
      real(dp),intent(inout) :: E_T
      
      sizeN = size(NPART)
!       curIndx = 0
      curIndx = NPart(sizeN)
      maxIndx = NMAX(sizeN) + 1 
      do i=1,sizeN-1
        curIndx = curIndx + maxIndx*NPart(sizeN-i)
        maxIndx = maxIndx * (NMAX(sizeN-i) + 1)
      enddo
      NHist(curIndx+1) = NHist(curIndx+1) + 1d0
      E_Avg(curIndx+1) = E_Avg(curIndx+1) + E_T

!      NHist(curIndx) = NHist(curIndx) + 1d0
!      E_Avg(curIndx) = E_Avg(curIndx) + E_T
      
      end subroutine
  !========================================================    
!     This subrotuine collects the cluster size Histogram from across each processor
!     and collapses them into a single output file for data analysis.  
      subroutine CollapseN
      use SimParameters
      use ParallelVar 
      use MPI 
      implicit none
      integer :: i, arraySize, lowerBound, upperBound
      real(dp), allocatable :: Hist_temp(:)


      lowerBound = lbound(NHist,1)
      upperBound = ubound(NHist,1)
      arraySize = size(NHist)
!      write(35,*) "Array Bounds:", lowerBound, upperBound
!      flush(35)
      allocate( Hist_Temp(lowerBound:upperBound) )

      Hist_temp = 0d0
      
      call Single_N_HistOutput
      write(nout,*) "Compressing Histogram Data..."           
      call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      call MPI_REDUCE(NHist, Hist_temp, arraySize, &
                 MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, &
                 ierror)  

      write(nout,*) "Finished Compressing Histogram Data..."      
      if(myid .eq. 0) then
       do i = lowerBound, upperBound
         NHist(i) = Hist_temp(i)
       enddo
      endif

 
      end subroutine
!==================================================================================
      subroutine N_HistOutput
      use SimParameters
      use ParallelVar
      implicit none
      integer :: i,j
      integer :: NArray(1:nMolTypes)
            
      open(unit = 92, file="Histogram_ClusterSize.txt")

      NArray = 0
!      NArray(nMolTypes) = NArray(nMolTypes) + 1      
      do i = 1, umbrellaLimit
       if(nMolTypes .gt. 1) then
         do j = 1,nMolTypes-1
          if(NArray(nMolTypes - j + 1) .gt. NMAX(nMolTypes - j + 1))then
           NArray(nMolTypes - j + 1) = 0
           NArray(nMolTypes - j) = NArray(nMolTypes - j) + 1          
          endif
         enddo
       endif
!       write(nout,*) NArray
       if(NHist(i) .ne. 0d0) then
         write(92, *) (NArray(j),j=1,nMolTypes), NHist(i)    
       endif       

       NArray(nMolTypes) = NArray(nMolTypes) + 1
      enddo
      close(92)
       
      end subroutine
!==================================================================================
      subroutine Single_N_HistOutput
      use SimParameters
      use ParallelVar
      implicit none
      integer :: i,j
      integer :: NArray(1:nMolTypes)
     
      NArray = 0
      write(35,*)"---------------------------------"
      write(35,*)"Cluster Size Histogram:"
!      NArray(nMolTypes) = NArray(nMolTypes) + 1      
      do i = 1, umbrellaLimit
       do j = 1,nMolTypes-1
         if(NArray(nMolTypes - j + 1) .gt. NMAX(nMolTypes - j + 1)) then
          NArray(nMolTypes - j + 1) = 0
          NArray(nMolTypes - j) = NArray(nMolTypes - j) + 1          
         endif
       enddo
!       write(nout,*) NArray
       if(NHist(i) .ne. 0d0) then
         write(35,*) (NArray(j),j=1,nMolTypes), NHist(i)    
       endif       
       NArray(nMolTypes) = NArray(nMolTypes) + 1       
      enddo   
      write(35,*)"---------------------------------"       
!      flush(35)
      end subroutine      
!==========================================================================================
      subroutine getNFromIndex(nIndx, NArray)
      use SimParameters
      implicit none
      integer,intent(in) :: nIndx
      integer,intent(out) :: NArray(1:nMolTypes)
      integer :: i
      integer :: curN, curScalar,curRemain, curMod

      curN = nIndx
      curMod = 1
      do i = 1, nMolTypes - 1
        curMod = curMod*(NMAX(i)+1)
      enddo

      do i = 2, nMolTypes
        curRemain = mod(curN,curMod)
        curScalar = nint(dble((curN - curRemain)/curMod))
        NArray(nMolTypes-i+1) = curScalar
        curN = curN - curScalar*curMod
        curMod = nint(dble(curMod/(NMAX(i-1)+1)))
      enddo
      
      NArray(1) = curN
       
      end subroutine        
!==========================================================================================
