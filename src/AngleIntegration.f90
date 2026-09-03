!=========================================================
       subroutine IntegrateBendAngleProb
       use SimParameters
       use ForceField
       use ForceFieldFunctions
       implicit none
 
       integer :: iAngle, iBin
       integer :: maxBin
       real(dp) :: ang_eq, k_eq
       real(dp) :: angle1, angle2, angle3
       real(dp) :: eng1, eng2, eng3
       real(dp) :: prob1, prob2, prob3, probMax
       real(dp) :: totalProb, norm
       real(dp) :: maxima, accConstant

!       write(35,*) "---------------------------"
       do iAngle = 1, nAngleTypes
         ang_eq = bendData(iAngle)%ang_eq
         k_eq = bendData(iAngle)%k_eq
!         write(35,*) k_eq,ang_eq
         bendData(iAngle)%Prob = 0d0
         call FindMaxima(ang_eq, k_eq, beta, maxima)
         maxBin = ceiling(maxima/bendBinWidth)
         eng1 = Harmonic(maxima, k_eq, ang_eq)
         probmax = sin(maxima)*exp(-beta*eng1)

         do iBin = 1, nBendHistBins
           angle1 = (iBin-1)*bendBinWidth
           angle2 = (iBin - 0.5d0)*bendBinWidth
           angle3 = iBin*bendBinWidth
           eng1 = Harmonic(angle1, k_eq, ang_eq)
           eng2 = Harmonic(angle2, k_eq, ang_eq)
           eng3 = Harmonic(angle3, k_eq, ang_eq)
           prob1 = sin(angle1)*exp(-beta*eng1)
           prob2 = sin(angle2)*exp(-beta*eng2)
           prob3 = sin(angle3)*exp(-beta*eng3)
           if(max(prob1,prob3)/probMax .lt. 1d-1) then
             totalProb = max(prob1,prob3)*bendBinWidth
           else
!             totalProb = 0.5d0*bendBinWidth*(prob1 + prob2)
             totalProb = (bendBinWidth/6d0)*(prob1 + 4d0*prob2 + prob3)
           endif
           bendData(iAngle)%Prob(iBin) = totalProb
         enddo
         norm = sum(bendData(iAngle)%Prob)
         do iBin = 1, nBendHistBins
           bendData(iAngle)%Prob(iBin) = bendData(iAngle)%Prob(iBin)/norm
!           write(35,*) (iBin-0.50d0)*bendBinWidth, bendData(iAngle)%Prob(iBin)
         enddo
!         write(35,*)

         
         accConstant = 0d0
         do iBin = 1, nBendHistBins
           angle1 = (iBin-1)*bendBinWidth
           angle2 = iBin*bendBinWidth
           eng1 = Harmonic(angle1, k_eq, ang_eq)
           eng2 = Harmonic(angle2, k_eq, ang_eq)
           prob1 = sin(angle1)*exp(-beta*eng1)
           prob2 = sin(angle2)*exp(-beta*eng2)

           prob1 = max(prob1,prob2)
           totalProb = bendBinWidth*prob1/bendData(iAngle)%Prob(iBin)

!           write(*,*) accConstant, totalProb
           if(totalProb .ge. accConstant) then
             accConstant = totalProb
           endif

         enddo

!         call FindMaxima(ang_eq, k_eq, beta, maxima)
!         maxBin = ceiling(maxima/bendBinWidth)
         eng1 = Harmonic(maxima, k_eq, ang_eq)
         prob1 = sin(maxima)*exp(-beta*eng1)
         totalProb = bendBinWidth*prob1/bendData(iAngle)%Prob(maxBin)

         if(totalProb .ge. accConstant) then
           accConstant = totalProb
         endif
!         write(*,*) accConstant, totalProb
         bendData(iAngle)%accptConstant = bendBinWidth/accConstant
!         write(35,*) iBin*bendBinWidth, bendData(iAngle)%Prob(1)
         bendData(iAngle)%startBin = 0
         do iBin = 2, nBendHistBins
           bendData(iAngle)%Prob(iBin) = bendData(iAngle)%Prob(iBin) + bendData(iAngle)%Prob(iBin-1)
           if(bendData(iAngle)%startBin .eq. 0) then     
             if(bendData(iAngle)%Prob(iBin) .ge. startProb) then
               bendData(iAngle)%startBin = iBin - 1
             endif
           endif
         enddo
       enddo
 

       end subroutine
!=========================================================
       subroutine FindMaxima(ang_eq, k_eq,beta, maxima)
       use VarPrecision
       use ForceField
       use Constants
       implicit none
 
       real(dp),intent(in) :: ang_eq, k_eq, beta
       real(dp),intent(out) :: maxima
       real(dp) :: x_new, x_old, tol
       real(dp) :: slopeTerm

       tol = 1d0 
       x_old = ang_eq
       do while(tol .gt. 1d-8)
         slopeTerm = (beta*k_eq*(x_old-ang_eq)*sin(x_old)-cos(x_old))
         slopeTerm = slopeTerm*sin(x_old)/(beta*k_Eq*sin(x_old)+1d0)
         x_new = x_old - slopeTerm
         tol = abs(slopeTerm)
         x_old = x_new
       enddo
       maxima = x_old
!        write(*,*) "max:",x_old

       end subroutine
!=========================================================
!       subroutine IntegrateTorsionalAngleProb
!
!       end subroutine
!=========================================================
