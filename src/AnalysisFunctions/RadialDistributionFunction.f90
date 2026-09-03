

!======================================================================================    
!      This module contains the functions nessisary to initalize and 
!      perform radial distribution studies
      module RadialDistribution
      use PairStorage
      use MiscelaniousVars
      
      private
      integer :: nRadialDist
      integer, allocatable :: radHistIndx(:)
      integer, allocatable :: radType1(:), radAtom1(:)
      integer, allocatable :: radType2(:), radAtom2(:)

      public :: nRadialDist
      public :: Initialize_RadialDist
      public :: SetRadialHist
      public :: Calc_RadialDist
      public :: Output_RadialDist
      public :: SetRadialParameters

      contains

!======================================================================================    
      subroutine Initialize_RadialDist
      implicit none 
      integer :: iRadial
      integer :: AllocationStatus
      integer :: startIndx, endIndx

      if(nRadialDist .eq. 0) then
        return
      endif

      allocate( radType1(1:nRadialDist), stat = AllocationStatus )
      allocate( radType2(1:nRadialDist), stat = AllocationStatus )
      allocate( radAtom1(1:nRadialDist), stat = AllocationStatus )
      allocate( radAtom2(1:nRadialDist), stat = AllocationStatus )
      allocate( radHistIndx(1:nRadialDist), stat = AllocationStatus )

      call ReserveSpace_Histograms(nRadialDist, startIndx, endIndx)

      do iRadial = 1, nRadialDist
        radHistIndx(iRadial) = startIndx + iRadial - 1
!        write(*,*) iRadial, radHistIndx(iRadial)
      enddo

      end subroutine
!======================================================================================    
      subroutine SetRadialParameters(iRadial, type1, type2, atom1, atom2)
      use PairStorage, only: rPair
      use SimParameters, only: NMAX
      use Coords
      implicit none 
      integer, intent(in) :: iRadial, type1, type2, atom1, atom2
      integer :: iMol, jMol
      integer :: gloIndx1, gloIndx2

      radType1(iRadial) = type1
      radType2(iRadial) = type2
      radAtom1(iRadial) = atom1
      radAtom2(iRadial) = atom2

      if(type1 .eq. type2) then
        do iMol = 1, NMAX(type1) - 1
          gloIndx1 = molArray(type1)%mol(iMol)%globalIndx(atom1)
          do jMol = iMol+1, NMAX(type1)
            gloIndx2 = molArray(type2)%mol(jMol)%globalIndx(atom2)            
            rPair(gloIndx1, gloIndx2)%p%storeRValue = .true.
          enddo
        enddo
      else
        do iMol = 1, NMAX(type1)
          gloIndx1 = molArray(type1)%mol(iMol)%globalIndx(atom1)
          do jMol = 1, NMAX(type2)
            gloIndx2 = molArray(type2)%mol(jMol)%globalIndx(atom2)
            rPair(gloIndx1, gloIndx2)%p%storeRValue = .true.
          enddo
        enddo
      endif

      end subroutine
!======================================================================================    
      subroutine SetRadialHist(iRadial, binSize, nBins, fileName)
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
      subroutine Calc_RadialDist
        use SimParameters, only: NPART
        use Coords
        implicit none
        integer :: iRadial
        integer :: bin, nBins
        integer :: nType1, nType2, nAtom1, nAtom2
        integer :: jMol, iMol
        integer :: gloIndx1, gloIndx2, radialIndx

        do iRadial = 1, nRadialDist
          nType1 = radType1(iRadial)
          nType2 = radType2(iRadial)
          nAtom1 = radAtom1(iRadial) 
          nAtom2 = radAtom2(iRadial) 
          radialIndx = radHistIndx(iRadial)
!          write(*,*) iRadial, radialIndx
          if(nType1 .eq. nType2) then
            do iMol = 1, NPART(nType1) - 1
              gloIndx1 = molArray(nType1)%mol(iMol)%globalIndx(nAtom1)
              do jMol = iMol+1, NPART(nType1)
                gloIndx2 = molArray(nType2)%mol(jMol)%globalIndx(nAtom2)
                bin = floor(rPair(gloIndx1, gloIndx2)%p%r * miscHist(radialIndx)%sizeInv)
                if(bin .le. miscHist(radialIndx)%nBins) then
                  miscHist(radialIndx)%binCount(bin) = miscHist(radialIndx)%binCount(bin) + 2d0
                else 
                  nBins = miscHist(radialIndx)%nBins
                  miscHist(radialIndx)%binCount(nBins+1) = miscHist(radialIndx)%binCount(nBins+1) + 2d0
                endif
              enddo
            enddo
          else
            do iMol = 1, NPART(nType1)
              gloIndx1 = molArray(nType1)%mol(iMol)%globalIndx(nAtom1)
              do jMol = 1, NPART(nType2)
                gloIndx2 = molArray(nType2)%mol(jMol)%globalIndx(nAtom2)
                bin = floor(rPair(gloIndx1, gloIndx2)%p%r * miscHist(radialIndx)%sizeInv)
                if(bin .le. miscHist(radialIndx)%nBins) then
                  miscHist(radialIndx)%binCount(bin) = miscHist(radialIndx)%binCount(bin) + 1d0
                else 
                  nBins = miscHist(radialIndx)%nBins
                  miscHist(radialIndx)%binCount(nBins+1) = miscHist(radialIndx)%binCount(nBins+1) + 1d0
                endif
              enddo
            enddo
          endif
        enddo
      
      end subroutine  
!======================================================================================    
      subroutine Output_RadialDist
        use Constants
        implicit none
        integer :: iRadial, iBin
        real(dp) :: d_bin
        real(dp) :: r, norm, rFactor

        open(unit = 81, file = "RadialDist_Info.txt")
        do iRadial = 1, nRadialDist
          open(unit = 80, file = miscHist(iRadial)%fileName)
          d_bin = miscHist(iRadial)%binSize
          norm = 0E0
          do iBin = 0, miscHist(iRadial)%nBins
            norm = norm + miscHist(iRadial)%binCount(iBin)
          enddo
!          write(*,*) norm
          write(81,*) "File Name:", miscHist(iRadial)%fileName
          write(81,*) "Number of Bins:", miscHist(iRadial)%nBins
          write(81,*) "Bin Size:", d_bin
          write(81,*) "Total Counts:", norm
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
