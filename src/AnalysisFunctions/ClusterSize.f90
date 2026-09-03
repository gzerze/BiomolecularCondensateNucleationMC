!==========================================================================
     module ClusterSizeUmbrella
        private
        integer :: nClusterVar
        integer, allocatable :: clusterArrayIndx(:)
        integer, allocatable :: typeIndxArray(:)

        public :: nClusterVar
        public :: Initialize_DistPair
        public :: SetPairVariables
        public :: CalcDistPairs
        public :: CalcDistPairs_New

        contains
     !--------------------------------------------------------------------------------
        subroutine Initialize_Cluster
           use MiscelaniousVars
           implicit none 
           integer :: AllocationStatus
           integer :: startIndx, endIndx, iCluster

           allocate(clusterArrayIndx(1:nClusterVar), STAT = AllocationStatus)
           allocate(typeIndxArray(1:nClusterVar), STAT = AllocationStatus)
          
           call ReserveSpace_Coord(nClusterVar, startIndx, endIndx)

           do iCluster = 1, nClusterVar
             pairArrayIndx(iCluster) = startIndx + iCluster - 1
           enddo
        end subroutine
     !--------------------------------------------------------------------------------
        subroutine CalcClusterSize
          use MiscelaniousVars
          use Coords
          use SimParameters, only: NPART
          implicit none 
          integer :: iCluster, iType, iIndx

          do iCluster =1, nClusterVar
            iType = typeIndxArray(iCluster)
            iIndx = clusterArrayIndx(iCluster)
            miscVar(iIndx) = real(NPART(iType), dp)
          enddo
        end subroutine
     !--------------------------------------------------------------------------------
        subroutine CalcClusterSize
          use MiscelaniousVars
          use Coords
          use SimParameters, only: NPART
          implicit none 
          integer :: iCluster, iType, iIndx

          do iCluster =1, nClusterVar
            iType = typeIndxArray(iCluster)
            iIndx = clusterArrayIndx(iCluster)
            miscVar(iIndx) = real(NPART(iType), dp)
          enddo
        end subroutine
     !--------------------------------------------------------------------------------
     end module
!==========================================================================
