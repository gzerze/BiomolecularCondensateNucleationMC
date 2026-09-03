!======================================================================================    
!      This module contains the functions nessisary to initalize and 
!      perform radial density distribution studies
      module RadialDensity
      use PairStorage
      use MiscelaniousVars
      
      private
      integer :: nRadialDens
      integer, allocatable :: radHistIndx(:)
      integer, allocatable :: radType(:)

      public :: nRadialDens
      public :: SetDensityHist
      public :: Initialize_RadialDens
      public :: SetDensityParameters
      public :: Calc_RadialDensity
      public :: Output_RadialDensity

      contains

!======================================================================================    
      subroutine Initialize_RadialDens
      implicit none 
      integer :: iRadial
      integer :: AllocationStatus
      integer :: startIndx, endIndx

      if(nRadialDens .eq. 0) then
        return
      endif

      allocate( radHistIndx(1:nRadialDens), stat = AllocationStatus )
      allocate( radType(1:nRadialDens), stat = AllocationStatus )

      call ReserveSpace_Histograms(nRadialDens, startIndx, endIndx)

      do iRadial = 1, nRadialDens
        radHistIndx(iRadial) = startIndx + iRadial - 1
      enddo

      end subroutine
!======================================================================================    
      subroutine SetDensityParameters(iRadial, type1)
      implicit none 
      integer, intent(in) :: iRadial, type1

      radType(iRadial) = type1

      end subroutine
!======================================================================================    
      subroutine SetDensityHist(iRadial, binSize, nBins, fileName)
      implicit none 
      integer, intent(in) :: iRadial, nBins
      real(dp), intent(in) :: binSize
      character(len=30), intent(in) :: fileName
      integer :: binIndx

      binIndx = radHistIndx(iRadial)
      miscHist(binIndx)%binSize = binSize
      miscHist(binIndx)%sizeInv = 1E0_dp/binSize
      miscHist(binIndx)%nBins = nBins
      miscHist(binIndx)%fileName = fileName

      end subroutine
!======================================================================================    
      subroutine Calc_RadialDensity
        use SimParameters, only: NPART, nMolTypes
        use Coords
        implicit none
        integer :: iRadial
        integer :: bin, nBins
        integer :: iType, nType
        integer :: iMol
        integer :: radialIndx
        integer :: cnt
        real(dp) :: r_sq, r, rx, ry, rz
        real(dp) :: xcm, ycm, zcm

        xcm = 0d0 
        cnt = 0
        ycm = 0d0 
        zcm = 0d0 
        do iType = 1, nMolTypes
          do iMol = 1, NPART(iType)
!            do iAtom = 1, nAtoms(iType)
              xcm = xcm + molArray(iType)%mol(iMol)%x(1)
              ycm = ycm + molArray(iType)%mol(iMol)%y(1)
              zcm = zcm + molArray(iType)%mol(iMol)%z(1) 
              cnt = cnt + 1
!            enddo
          enddo
        enddo

        xcm = xcm/real(cnt, dp)
        ycm = ycm/real(cnt, dp)
        zcm = zcm/real(cnt, dp)

        do iRadial = 1, nRadialDens
          nType = radType(iRadial)
          radialIndx = radHistIndx(iRadial)
          do iMol = 1, NPART(nType)
            rx = molArray(nType)%mol(iMol)%x(1) - xcm
            ry = molArray(nType)%mol(iMol)%y(1) - ycm
            rz = molArray(nType)%mol(iMol)%z(1) - zcm
            r_sq = rx*rx + ry*ry + rz*rz
            r = sqrt(r_sq)
            bin = floor(r * miscHist(radialIndx)%sizeInv)
            if(bin .le. miscHist(radialIndx)%nBins) then
              miscHist(radialIndx)%binCount(bin) = miscHist(radialIndx)%binCount(bin) + 1d0
            else 
              nBins = miscHist(radialIndx)%nBins
              miscHist(radialIndx)%binCount(nBins+1) = miscHist(radialIndx)%binCount(nBins+1) + 1d0
            endif            
          enddo
        enddo
      
      end subroutine  
!======================================================================================    
      subroutine Output_RadialDensity
        use Constants
        implicit none
        integer :: iRadial, iBin
        real(dp) :: d_bin
        real(dp) :: r, norm, rFactor


        open(unit=81, file="RadialDensity_Info.txt")
        do iRadial = 1, nRadialDens
          open(unit = 80, file = miscHist(iRadial)%fileName)
          d_bin = miscHist(iRadial)%binSize
          norm = 0E0
          do iBin = 0, miscHist(iRadial)%nBins
            norm = norm + miscHist(iRadial)%binCount(iBin)
          enddo
!          write(*,*) norm
          write(81,*) "Number of Bins:", miscHist(iRadial)%nBins
          write(81,*) "Bin Size:", d_bin
          write(81,*) "Total Count:", norm
          write(81,*)
          do iBin = 0, miscHist(iRadial)%nBins
            r = iBin * d_bin
            rFactor = norm * 4E0/3E0 * pi * ( (r+d_bin)**3 - r**3 )
            write(80,*) r, miscHist(iRadial)%binCount(iBin)/rFactor
          enddo
          close(80)
        enddo
        close(81)
      
      end subroutine 
!======================================================================================    
      end module
!======================================================================================    
