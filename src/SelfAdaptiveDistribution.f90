!===================================================================
     module SelfAdaptive
     use VarPrecision
     contains
!===================================================================
     subroutine UpdateDihedralProbabilites
     use Constants
     use CBMC_Variables
     implicit none
     integer :: iDihed, iBin
     real(dp) :: histNorm, probNorm
     real(dp) :: avg

     do iDihed = 1, totalDihed
       histNorm = sum(dihedData(iDihed)%Hist)
!       Compute the new distribution by averaging the new and old probability estimates
       do iBin = 0, nDihBins
         avg = dihedData(iDihed)%Prob(iBin) + dihedData(iDihed)%Hist(iBin)/histNorm
         avg = avg * 0.5E0
         dihedData(iDihed)%Prob(iBin) = avg
       enddo
!       Normalize the new distribution
       probNorm = sum(dihedData(iDihed)%Prob)
       do iBin = 0, nDihBins
         dihedData(iDihed)%Prob(iBin) = dihedData(iDihed)%Prob(iBin)/probNorm
       enddo



!       Tabulate the integral table
       dihedData(iDihed)%Integral(0) = dihedData(iDihed)%Prob(iBin)
       do iBin = 1, nDihBins
         dihedData(iDihed)%Integral(iBin) = dihedData(iDihed)%Integral(iBin-1) + dihedData(iDihed)%Prob(iBin)
       enddo
       dihedData(iDihed)%Integral(nDihBins) = 1E0

!       Determine the new acceptance constant
       dihedData(iDihed)%accConst = diBinSize

 

!       write(2,*) iDihed, dihedData(iDihed)%accConst 
!       do iBin = 0, nDihBins
!         write(2,*) iBin*diBinSize, dihedData(iDihed)%Integral(iBin)
!       enddo
!       write(2,*)
!       flush(2)
     enddo

     end subroutine
!===================================================================
     subroutine OutputDihedral
     use Constants
     use CBMC_Variables
     implicit none
     integer :: iDihed, iBin
     real(dp) :: histNorm

     do iDihed = 1, totalDihed
       histNorm = sum(dihedData(iDihed)%Hist)
       do iBin = 0, nDihBins
         write(2,*) iBin*diBinSize, dihedData(iDihed)%Hist(iBin)/histNorm
       enddo
     enddo

     end subroutine

!=================================================================== 
     end module
