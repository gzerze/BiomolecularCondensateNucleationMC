!*********************************************************************************************************************
!     This file contains the energy functions that work the Tersoff Model
!     The prefix naming scheme implies the following:
!           Detailed - Complete energy calculation inteded for use at the beginning and end
!                      of the simulation.  This function is not inteded for use mid-simulation.
!             Shift  - Calculates the energy difference for any move that does not result
!                      in molecules being added or removed from the cluster. This function
!                      receives any number of Displacement vectors from the parent function as input.
!              Mol   - Calculates the energy for a molecule already present in the system. For
!                      use in moves such as particle deletion moves. 
!             NewMol - Calculates the energy for a molecule that has been freshly inserted into the system.
!                      Intended for use in Rosenbluth Sampling, Swap In, etc.
!           Exchange - Combines the Mol and New Mol routines for moves that simultaniously add and remove a particle at the same time.
!*********************************************************************************************************************
      module InterEnergy_Tersoff
      use VarPrecision



!      real(dp), parameter :: ep = 4d0*634.62810518d0
!      real(dp), parameter :: ep = 4d0*78d0
      real(dp), parameter :: ep = 0d0
      
      real(dp), parameter :: sig = 3.174116d0**2
!====================================================================================== 
      contains
!======================================================================================      
      pure function LJ_Func(r_sq, epx, sigx) result(LJ)
      implicit none
      real(dp), intent(in) :: r_sq, epx, sigx
      real(dp) :: LJ  
 
      LJ = 0E0_dp
      return

      LJ = (sigx/r_sq)
      LJ = LJ * LJ * LJ
      LJ = epx * LJ * (LJ-1E0_dp)  

      end function
!====================================================================================== 
      pure function Fc_Func(r, R_eq, D) result(val)
      use Constants, only: pi
      implicit none
      real(dp), intent(in) :: r, R_eq, D
      real(dp) :: val  
 
      if( r .lt. (R_eq-D) ) then
        val = 1E0_dp
      elseif( r .lt. (R_eq+D) ) then
        val = 0.5E0_dp * (1E0_dp - sin(pi*(r-R_eq)/(2E0_dp*D)) )
      else
        val = 0E0_dp
      endif

      end function
!======================================================================================      
      pure function gik_Func(theta, c, d, h) result(val)
      implicit none
      real(dp), intent(in) :: theta, c, d, h
      real(dp) :: c_sq, d_sq
      real(dp) :: val  
 
      c_sq = c * c
      d_sq = d * d
      val = 1E0_dp + c_sq/d_sq - c_sq/(d_sq + (cos(theta) - h)**2) 

      end function

!======================================================================================
      pure function angleCalc(rx12, ry12, rz12, r12, rx23, ry23, rz23, r23) result(Angle)
      implicit none
      real(dp), intent(in) :: rx12, ry12, rz12, r12, rx23, ry23, rz23, r23
      real(dp) :: Angle  
   
            
      Angle = rx12*rx23 + ry12*ry23 + rz12*rz23
      Angle = Angle/(r12*r23)
      if(abs(Angle) .gt. 1E0_dp) then
        Angle = sign(1E0_dp, Angle)
      endif
      Angle = acos(Angle)


      end function
!======================================================================================      
      subroutine Detailed_ECalc_Inter(E_T, PairList)
      use ParallelVar
      use ForceField
      use ForceFieldPara_Tersoff
      use Coords
      use Constants, only: pi
      use SimParameters
      use EnergyTables
      use PairStorage, only: rPair
      implicit none
      real(dp), intent(inOut) :: E_T
      real(dp), intent(inOut) :: PairList(:,:)
      integer :: iType, jType, kType
      integer :: iMol, jMol, kMol
      integer(kind=atomIntType) :: atmType1     
      integer :: iIndx, jIndx, globIndx1, globIndx2, globIndx3
      real(dp) :: rMax, rMax_sq

      real(dp) :: A, B, c, d, R_eq, D2 
      real(dp) :: E_Short
      real(dp) :: lam1, lam2
      real(dp) :: Zeta
      real(dp) :: BetaPar, n, h
      real(dp) :: b1, V1
      real(dp) :: angijk

      real(dp) :: rxij, ryij, rzij, rij
      real(dp) :: rxik, ryik, rzik, rik

      real(dp) :: E_LJ, LJ

      E_T = 0E0_dp
      E_Short = 0E0_dp
      E_LJ = 0E0_dp
      PairList = 0E0_dp
      ETable = 0E0_dp
      iType = 1
      jType = 1
      kType = 1
      atmType1 = atomArray(iType, 1) 

      A = tersoffData(atmType1)%A
      B = tersoffData(atmType1)%B
      c = tersoffData(atmType1)%c
      d = tersoffData(atmType1)%d
      R_eq = tersoffData(atmType1)%R
      D2 = tersoffData(atmType1)%D2
      BetaPar = tersoffData(atmType1)%beta
      n = tersoffData(atmType1)%n
      h = tersoffData(atmType1)%h
      lam1 = tersoffData(atmType1)%lam1
      lam2 = tersoffData(atmType1)%lam2

!      write(*,*) A, B, c, d, R_eq, D2, BetaPar, n, h, lam1, lam2

      rMax = R_eq + D2
      rMax_sq = rMax * rMax
!      write(*,*) rMax, R_eq - D2

      do iMol = 1, NPART(iType)
        globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(1)
        iIndx = MolArray(iType)%mol(iMol)%indx
        do jMol = 1, NPART(jType)
          if(iMol .eq. jMol) then
            cycle
          endif
          globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
          rij  = rPair(globIndx1, globIndx2)%p%r
!          write(*,*) rij, Fc_Func(rij, R_eq, D2)
          jIndx = MolArray(jType)%mol(jMol)%indx
          if(distCriteria) then
            PairList(iIndx, jIndx) = rPair(globIndx1, globIndx2)%p%r_sq
            PairList(jIndx, iIndx) = PairList(iIndx,jIndx)
          endif


          LJ = 0.5d0*LJ_Func(rPair(globIndx1, globIndx2)%p%r_sq, ep, sig) 
          E_LJ = E_LJ + LJ
!          write(*,*) "LJ:", rPair(globIndx1, globIndx2)%p%r_sq, LJ, ep, sig
          PairList(iIndx, jIndx) = PairList(iIndx, jIndx) + LJ
          PairList(jIndx, iIndx) = PairList(jIndx, iIndx) + LJ
          ETable(iIndx) = ETable(iIndx) + LJ
          ETable(jIndx) = ETable(jIndx) + LJ


          if(rij .gt. rMax) then
            cycle
          endif

          Zeta = 0E0_dp

          rxij = rPair(globIndx1, globIndx2)%p%rx
          ryij = rPair(globIndx1, globIndx2)%p%ry
          rzij = rPair(globIndx1, globIndx2)%p%rz
          if(globIndx2 .gt. globIndx1) then
            rxij = -rxij
            ryij = -ryij
            rzij = -rzij
          endif
          do kMol = 1, nPart(kType)
            if((kMol .eq. iMol) .or. (kMol .eq. jMol)) then
              cycle
            endif
            globIndx3 = MolArray(kType)%mol(kMol)%globalIndx(1)
            rik  = rPair(globIndx1, globIndx3)%p%r
            if(rik .lt. rMax) then
              rxik  = rPair(globIndx1, globIndx3)%p%rx
              ryik  = rPair(globIndx1, globIndx3)%p%ry
              rzik  = rPair(globIndx1, globIndx3)%p%rz

              if(globIndx3 .gt. globIndx1) then
                rxik = -rxik
                ryik = -ryik
                rzik = -rzik
              endif
              angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
              Zeta = Zeta + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2)
            endif
          enddo
          if(Zeta .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b1 = 1E0_dp
          endif      
    
          V1 = 0E0_dp
!          write(*,*) "rij:", rij
          V1 = 0.5E0_dp * Fc_Func(rij, R_eq, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij))
          E_Short = E_Short + V1
!          write(*,*) "V1:", b1, V1
          PairList(iIndx, jIndx) = PairList(iIndx, jIndx) + V1
          PairList(jIndx, iIndx) = PairList(jIndx, iIndx) + V1
          ETable(iIndx) = ETable(iIndx) + V1
          ETable(jIndx) = ETable(jIndx) + V1

!          write(*,*) 
        enddo
      enddo
!      E_Short = 0.5E0_dp*E_Short 
      write(nout,*) "Lennard-Jones Energy:", E_LJ/outputEConv
      write(nout,*) "ShortRange Energy:", E_Short/outputEConv

      E_T = E_T + E_Short + E_LJ
      E_Inter_T = E_Short + E_LJ
      
      end subroutine
!======================================================================================      
      subroutine Shift_ECalc_Inter(E_Trial, disp, newDist, PairList, dETable)
      use ForceField
      use ForceFieldPara_Tersoff
      use Coords
      use SimParameters
      use Constants, only: pi
      use PairStorage, only: rPair, DistArrayNew, nNewDist, rPairNew 
      implicit none
      
      type(Displacement), intent(in) :: disp(:)  
      type(DistArrayNew), intent(inout) :: newDist(:)
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
!      logical, intent(out) :: rejMove

      
      integer :: iType, jType, kType, nType, iPair
      integer :: iMol, jMol, kMol, nMol
      integer :: iNei, jNei
      integer(kind=atomIntType) :: atmType1    
      integer :: iIndx, jIndx, nIndx 
      integer :: globIndx1, globIndx2, globIndx3, globIndxN
      integer :: neiList(1:60), nNei
!      integer :: pairIndxNew(1:6), nPair
      real(dp) :: r_new,rMax, rMax_sq

      real(dp) :: A, B, c, d, R_eq, D2 
      real(dp) :: E_Short, Short, LJ, E_LJ
      real(dp) :: lam1, lam2
      real(dp) :: Zeta, Zeta2
      real(dp) :: BetaPar, n, h
      real(dp) :: b1, b2, V1
      real(dp) :: angijk

      real(dp) :: rxij, ryij, rzij, rij
      real(dp) :: rxjk, ryjk, rzjk, rjk
      real(dp) :: rxik, ryik, rzik, rik

      E_Trial = 0E0_dp
      E_Short = 0E0_dp
      PairList = 0E0_dp
      dETable = 0E0_dp
      nType = 1
      iType = 1
      jType = 1
      kType = 1
      atmType1 = atomArray(iType, 1) 

      A = tersoffData(atmType1)%A
      B = tersoffData(atmType1)%B
      c = tersoffData(atmType1)%c
      d = tersoffData(atmType1)%d
      R_eq = tersoffData(atmType1)%R
      D2 = tersoffData(atmType1)%D2
      BetaPar = tersoffData(atmType1)%beta
      n = tersoffData(atmType1)%n
      h = tersoffData(atmType1)%h
      lam1 = tersoffData(atmType1)%lam1
      lam2 = tersoffData(atmType1)%lam2

      rMax = R_eq + D2
      rMax_sq = rMax * rMax

      nType = disp(1)%molType
      nMol = disp(1)%molIndx
      nIndx = MolArray(nType)%mol(nMol)%indx
      globIndx1 = MolArray(nType)%mol(nMol)%globalIndx(1)
      globIndxN = globIndx1
      nNei = 0
      neiList = 0
!      pairIndxNew = 0

      E_LJ = 0E0_dp
      do jMol = 1, NPART(jType)
        if(jMol .eq. nMol) then
          cycle
        endif
        jIndx = MolArray(jType)%mol(jMol)%indx
        globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
        LJ = LJ_Func(rPairNew(globIndxN,globIndx2)%p%r_sq, ep, sig)
        LJ = LJ - LJ_Func(rPair(globIndxN, globIndx2)%p%r_sq, ep, sig)
        dETable(nIndx) = dETable(nIndx) + LJ
        dETable(jIndx) = dETable(jIndx) + LJ
        E_LJ = E_LJ + LJ
      enddo

!       In the Tersoff model, the strength of a given molecular bond is dependent on both the distance and the local environment around the bond.  As a result
!       when one shifts a particle's location one must not only calculate the bonds that changed during this time frame, but also calculate how that
!       shift impacts the bonds of it's neighbors.  Thus it is nessisary to compute a list of particles who may have been impacted. 


!      This block creates a list of neighbors that are located near the particle's new position.
      do iPair = 1, nNewDist
        r_new = newDist(iPair)%r
        globIndx2 = newDist(iPair)%indx2
        jIndx = globIndx2
        if(r_new .lt. rMax) then
          nNei = nNei + 1
          neiList(nNei) = newDist(iPair)%indx2
!          pairIndxNew(nNei) = iPair
        endif
      enddo

!      This block adds to the list of neighbors particles that are located near the mobile particle's old position.

      do jMol = 1, NPART(jType)
        if(jMol .eq. nMol) then
          cycle
        endif
        jIndx = MolArray(jType)%mol(jMol)%indx

        globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
        if(distCriteria) then
          PairList(jIndx) = rPairNew(globIndxN,globIndx2)%p%r_sq
        endif
        if(globIndx2 .eq. globIndxN) then
          cycle
        endif
        if(rPair(globIndxN, globIndx2)%p%r .lt. rMax ) then
          if( all(neiList(1:nNei) .ne. globIndx2) ) then
            nNei = nNei + 1
            neiList(nNei) = globIndx2
          endif
        endif
      enddo

!      This portion of the code calculates the bonds that originate from the particle that just moved.
      do jNei = 1, nNei
        jMol = neiList(jNei)
        globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
        jIndx = MolArray(jType)%mol(jMol)%indx
        Short = 0E0_dp
!        Compute the bonds at the new position
        rij  = rPairNew(globIndxN,globIndx2)%p%r
        if(rij .lt. rMax) then
          Zeta = 0E0_dp
          Zeta2 = 0E0_dp
          rxij = -rPairNew(globIndxN,globIndx2)%p%rx
          ryij = -rPairNew(globIndxN,globIndx2)%p%ry
          rzij = -rPairNew(globIndxN,globIndx2)%p%rz 
          do kMol = 1, NPART(kType)
            if((kMol .eq. nMol) .or. (kMol .eq. jMol)) then
              cycle
            endif
            globIndx3 = MolArray(kType)%mol(kMol)%globalIndx(1)
            rik  = rPairNew(globIndxN,globIndx3)%p%r

            if(rik .lt. rMax) then
              rxik  = -rPairNew(globIndxN,globIndx3)%p%rx
              ryik  = -rPairNew(globIndxN,globIndx3)%p%ry
              rzik  = -rPairNew(globIndxN,globIndx3)%p%rz
              angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
              Zeta = Zeta + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2)
            endif     

            rjk  = rPair(globIndx2, globIndx3)%p%r
            if(rjk .lt. rMax) then
              rxjk  = rPair(globIndx2, globIndx3)%p%rx
              ryjk  = rPair(globIndx2, globIndx3)%p%ry
              rzjk  = rPair(globIndx2, globIndx3)%p%rz
              if(globIndx3 .gt. globIndx2) then
                rxjk = -rxjk
                ryjk = -ryjk
                rzjk = -rzjk
              endif
              angijk = angleCalc(-rxij, -ryij, -rzij, rij, rxjk, ryjk, rzjk, rjk)
              Zeta2 = Zeta2 + gik_Func(angijk, c, d, h) *  Fc_Func(rjk, R_eq, D2)
            endif  
          enddo
          if(Zeta .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b1 = 1E0_dp
          endif
          if(Zeta2 .ne. 0E0_dp) then
            b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b2 = 1E0_dp
          endif
!          write(2,*) b1
!          write(2,*) b2

          V1 = Fc_Func(rij, R_eq, D2) * (2d0*A*exp(-lam1*rij) - (b1+b2)*B*exp(-lam2*rij)) 
          if(.not. distCriteria) then                
            PairList(jIndx) = PairList(jIndx) + 0.5E0_dp * V1
          endif

          Short = Short + V1
        endif

   

!        Compute the bonds at the old position
        rij  = rPair(globIndx1, globIndx2)%p%r
        if(rij .lt. rMax) then
          Zeta = 0E0_dp
          Zeta2 = 0E0_dp
          rxij = rPair(globIndx1, globIndx2)%p%rx
          ryij = rPair(globIndx1, globIndx2)%p%ry
          rzij = rPair(globIndx1, globIndx2)%p%rz
          if(globIndx2 .gt. globIndx1) then
            rxij = -rxij
            ryij = -ryij
            rzij = -rzij
          endif
          do kMol = 1, NPART(kType)
            if((kMol .eq. nMol) .or. (kMol .eq. jMol)) then
              cycle
            endif
            globIndx3 = MolArray(kType)%mol(kMol)%globalIndx(1)
            rik = rPair(globIndx1, globIndx3)%p%r
            if(rik .lt. rMax) then
              rxik  = rPair(globIndx1, globIndx3)%p%rx
              ryik  = rPair(globIndx1, globIndx3)%p%ry
              rzik  = rPair(globIndx1, globIndx3)%p%rz
              if(globIndx3 .gt. globIndx1) then
                rxik = -rxik
                ryik = -ryik
                rzik = -rzik
              endif
              angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
              Zeta = Zeta + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2)
            endif

            rjk = rPair(globIndx2, globIndx3)%p%r
            if(rjk .lt. rMax) then
              rxjk  = rPair(globIndx2, globIndx3)%p%rx
              ryjk  = rPair(globIndx2, globIndx3)%p%ry
              rzjk  = rPair(globIndx2, globIndx3)%p%rz
              if(globIndx3 .gt. globIndx2) then
                rxjk = -rxjk
                ryjk = -ryjk
                rzjk = -rzjk
              endif
              angijk = angleCalc(-rxij, -ryij, -rzij, rij, rxjk, ryjk, rzjk, rjk)
              Zeta2 = Zeta2 + gik_Func(angijk, c, d, h) *  Fc_Func(rjk, R_eq, D2)
            endif  
          enddo
          if(Zeta .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b1 = 1E0_dp
          endif
          if(Zeta2 .ne. 0E0_dp) then
            b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b2 = 1E0_dp
          endif
          V1 = Fc_Func(rij, R_eq, D2) * (2E0_dp*A*exp(-lam1*rij) - (b1+b2)*B*exp(-lam2*rij))
          Short = Short - V1
        endif

        dETable(nIndx) = dETable(nIndx) + 0.5E0_dp*Short
        dETable(jIndx) = dETable(jIndx) + 0.5E0_dp*Short
        E_Short = E_Short + Short
      enddo

      if(NTotal .eq. 2) then
        E_Trial = 0.5E0_dp * E_Short + E_LJ
        return
      endif



!      This portion of the code computes the energy of the atoms which neighbored
!      the particle that moved. 
      do iNei = 1, nNei
        iMol = neiList(iNei)
        globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(1)
        iIndx = MolArray(iType)%mol(iMol)%indx

        do jMol = 1, NPART(jType)
          if(iMol .eq. jMol) then
            cycle
          endif
          if(nMol .eq. jMol) then
            cycle
          endif

          globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
          rij  = rPair(globIndx1, globIndx2)%p%r
          if(rij .gt. rMax) then
            cycle
          endif

          Zeta = 0E0_dp        !Zeta for the New Configuration
          Zeta2 = 0E0_dp       !Zeta for the Old Configuration
          jIndx = MolArray(jType)%mol(jMol)%indx
 
          rxij = rPair(globIndx1, globIndx2)%p%rx
          ryij = rPair(globIndx1, globIndx2)%p%ry
          rzij = rPair(globIndx1, globIndx2)%p%rz
          if(globIndx2 .gt. globIndx1) then
            rxij = -rxij
            ryij = -ryij
            rzij = -rzij 
          endif

          do kMol = 1, nPart(kType)
            if((kMol .eq. iMol) .or. (kMol .eq. jMol)) then
              cycle
            endif
            globIndx3 = MolArray(kType)%mol(kMol)%globalIndx(1)
            if(kMol .eq. nMol) then
              rik  = rPairNew(globIndxN,globIndx1)%p%r
              if(rik .lt. rMax) then
                rxik  = rPairNew(globIndxN,globIndx1)%p%rx
                ryik  = rPairNew(globIndxN,globIndx1)%p%ry
                rzik  = rPairNew(globIndxN,globIndx1)%p%rz
                angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
                Zeta = Zeta + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2) 
              endif

              rik  = rPair(globIndx1, globIndx3)%p%r
              if(rik .lt. rMax) then
                rxik  = rPair(globIndx1, globIndx3)%p%rx
                ryik  = rPair(globIndx1, globIndx3)%p%ry
                rzik  = rPair(globIndx1, globIndx3)%p%rz
                if(globIndx3 .gt. globIndx1) then
                  rxik = -rxik
                  ryik = -ryik
                  rzik = -rzik
                endif
                angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
                Zeta2 = Zeta2 + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2) 
              endif
              
            else
              rik  = rPair(globIndx1, globIndx3)%p%r
              if(rik .lt. rMax) then
                rxik  = rPair(globIndx1, globIndx3)%p%rx
                ryik  = rPair(globIndx1, globIndx3)%p%ry
                rzik  = rPair(globIndx1, globIndx3)%p%rz
                rik   = rPair(globIndx1, globIndx3)%p%r
                if(globIndx3 .gt. globIndx1) then
                  rxik = -rxik
                  ryik = -ryik
                  rzik = -rzik
                endif
                angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
                V1 = gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2) 
                Zeta = Zeta + V1
                Zeta2 = Zeta2 + V1 
              endif
            endif
          enddo
          if(Zeta .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b1 = 1E0_dp
          endif
          if(Zeta2 .ne. 0E0_dp) then
            b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b2 = 1E0_dp
          endif

          V1 = Fc_Func(rij, R_eq, D2) * (B*exp(-lam2*rij))*(b2 - b1)
          dETable(iIndx) = dETable(iIndx) + 0.5E0_dp*V1
          dETable(jIndx) = dETable(jIndx) + 0.5E0_dp*V1
          E_Short = E_Short + V1
        enddo
      enddo



      E_Trial = 0.5E0_dp*E_Short + E_LJ
!      write(*,*) E_Trial, E_Short, E_LJ

      end subroutine 
!======================================================================================      
      subroutine Mol_ECalc_Inter(nType, nMol, dETable, E_Trial)
      use ForceField
      use ForceFieldPara_Tersoff
      use Coords
      use SimParameters
      use PairStorage
      implicit none
      integer, intent(in) :: nType, nMol     
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: dETable(:)
      
      
      integer :: iType, jType, kType
      integer :: iMol, jMol, kMol
      integer :: iNei, jNei
      integer(kind=atomIntType) :: atmType1      
      integer :: iIndx, jIndx, nIndx 
      integer :: globIndx1, globIndx2, globIndx3
      integer :: neiList(1:60), nNei
!      integer :: pairIndxNew(1:6), nPair
      real(dp) :: rMax, rMax_sq

      real(dp) :: A, B, c, d, R_eq, D2 
      real(dp) :: Short, LJ, E_LJ
      real(dp) :: lam1, lam2
      real(dp) :: Zeta, Zeta2
      real(dp) :: BetaPar, n, h
      real(dp) :: b1, b2, V1
      real(dp) :: angijk

      real(dp) :: rxij, ryij, rzij, rij
      real(dp) :: rxjk, ryjk, rzjk, rjk
      real(dp) :: rxik, ryik, rzik, rik

      E_Trial = 0E0_dp
      dETable = 0E0_dp
      iType = 1
      jType = 1
      kType = 1
      atmType1 = atomArray(nType, 1) 

      A = tersoffData(atmType1)%A
      B = tersoffData(atmType1)%B
      c = tersoffData(atmType1)%c
      d = tersoffData(atmType1)%d
      R_eq = tersoffData(atmType1)%R
      D2 = tersoffData(atmType1)%D2
      BetaPar = tersoffData(atmType1)%beta
      n = tersoffData(atmType1)%n
      h = tersoffData(atmType1)%h
      lam1 = tersoffData(atmType1)%lam1
      lam2 = tersoffData(atmType1)%lam2

      rMax = R_eq + D2
      rMax_sq = rMax * rMax
      nIndx = MolArray(nType)%mol(nMol)%indx

!       In the Tersoff model, the strength of a given molecular bond is dependent on both the distance and the local environment around the bond.  As a result
!       when one shifts a particle's location one must not only calculate the bonds that changed during this time frame, but also calculate how that
!       shift impacts the bonds of it's neighbors.  Thus it is nessisary to compute a list of particles who may have been implacted. 
      globIndx1 = molArray(nType)%mol(nMol)%globalIndx(1)
      nNei = 0
      neiList = 0
      E_LJ = 0E0_dp
      do jMol = 1, NPART(jType)
        if(jMol .eq. nMol) then
          cycle
        endif
        jIndx = MolArray(jType)%mol(jMol)%indx

        globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
        if(rPair(globIndx1, globIndx2)%p%r .lt. rMax ) then
          if( all(neiList(1:nNei) .ne. globIndx2) ) then
            nNei = nNei + 1
            neiList(nNei) = globIndx2
          endif
        endif
        LJ = LJ_Func(rPair(globIndx1, globIndx2)%p%r_sq, ep, sig) 
        dETable(nIndx) = dETable(nIndx) + LJ
        dETable(jIndx) = dETable(jIndx) + LJ
        E_LJ = E_LJ + LJ
      enddo

!      write(*,*) 1
!      This block calculates the energy penalty for removing a molecule from the cluster. 
      do jNei = 1, nNei
        jMol = neiList(jNei)
        globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
        jIndx = MolArray(jType)%mol(jMol)%indx
        Short = 0E0_dp
        rij  = rPair(globIndx1, globIndx2)%p%r
        Zeta = 0E0_dp
        Zeta2 = 0E0_dp
        rxij = rPair(globIndx1, globIndx2)%p%rx
        ryij = rPair(globIndx1, globIndx2)%p%ry
        rzij = rPair(globIndx1, globIndx2)%p%rz
        if(globIndx2 .gt. globIndx1) then
          rxij = -rxij
          ryij = -ryij
          rzij = -rzij
        endif
        do kMol = 1, NPART(kType)
          if((kMol .eq. nMol) .or. (kMol .eq. jMol)) then
            cycle
          endif
          globIndx3 = MolArray(kType)%mol(kMol)%globalIndx(1)
          rik = rPair(globIndx1, globIndx3)%p%r
!          This block calculates the energy related to the 1,2,3 angle
          if(rik .lt. rMax) then
            rxik  = rPair(globIndx1, globIndx3)%p%rx
            ryik  = rPair(globIndx1, globIndx3)%p%ry
            rzik  = rPair(globIndx1, globIndx3)%p%rz
            if(globIndx3 .gt. globIndx1) then
              rxik = -rxik
              ryik = -ryik
              rzik = -rzik
            endif
            angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
            Zeta = Zeta + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2)
          endif

!          This block calculates the energy related to the 2,1,3 angle
          rjk = rPair(globIndx2, globIndx3)%p%r
          if(rjk .lt. rMax) then
            rxjk  = rPair(globIndx2, globIndx3)%p%rx
            ryjk  = rPair(globIndx2, globIndx3)%p%ry
            rzjk  = rPair(globIndx2, globIndx3)%p%rz
            if(globIndx3 .gt. globIndx2) then
              rxjk = -rxjk
              ryjk = -ryjk
              rzjk = -rzjk
            endif
            angijk = angleCalc(-rxij, -ryij, -rzij, rij, rxjk, ryjk, rzjk, rjk)
            Zeta2 = Zeta2 + gik_Func(angijk, c, d, h) *  Fc_Func(rjk, R_eq, D2)
          endif  
        enddo
        if(Zeta .ne. 0E0_dp) then
          b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
        else
          b1 = 1E0_dp
        endif
        if(Zeta2 .ne. 0E0_dp) then
          b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
        else
          b2 = 1E0_dp
        endif
        V1 = 0.5E0_dp*Fc_Func(rij, R_eq, D2) * (2d0*A*exp(-lam1*rij) - (b1+b2)*B*exp(-lam2*rij))
        Short = Short + V1
!        write(*,*) 1.5
        dETable(nIndx) = dETable(nIndx) + Short
        dETable(jIndx) = dETable(jIndx) + Short
        E_Trial = E_Trial + Short
      enddo


!      This portion of the code computes the energy of the atoms which neighbored
!      the particle that is being removed. 
      do iNei = 1, nNei
        iMol = neiList(iNei)
        globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(1)
        iIndx = MolArray(iType)%mol(iMol)%indx

        do jMol = 1, NPART(jType)
          if(iMol .eq. jMol) then
            cycle
          endif
          if(nMol .eq. jMol) then
            cycle
          endif

          globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
          rij  = rPair(globIndx1, globIndx2)%p%r
          if(rij .gt. rMax) then
            cycle
          endif

          Zeta = 0E0_dp        !Zeta for the New Configuration
          Zeta2 = 0E0_dp       !Zeta for the Old Configuration
          jIndx = MolArray(jType)%mol(jMol)%indx
 
          rxij = rPair(globIndx1, globIndx2)%p%rx
          ryij = rPair(globIndx1, globIndx2)%p%ry
          rzij = rPair(globIndx1, globIndx2)%p%rz
          if(globIndx2 .gt. globIndx1) then
            rxij = -rxij
            ryij = -ryij
            rzij = -rzij 
          endif

          do kMol = 1, nPart(kType)
            if((kMol .eq. iMol) .or. (kMol .eq. jMol)) then
              cycle
            endif
            globIndx3 = MolArray(kType)%mol(kMol)%globalIndx(1)
            if(kMol .eq. nMol) then
              rik  = rPair(globIndx1, globIndx3)%p%r
              if(rik .lt. rMax) then
                rxik  = rPair(globIndx1, globIndx3)%p%rx
                ryik  = rPair(globIndx1, globIndx3)%p%ry
                rzik  = rPair(globIndx1, globIndx3)%p%rz
                if(globIndx3 .gt. globIndx1) then
                  rxik = -rxik
                  ryik = -ryik
                  rzik = -rzik
                endif
                angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
                Zeta2 = Zeta2 + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2) 
              endif              
            else
              rik  = rPair(globIndx1, globIndx3)%p%r
              if(rik .lt. rMax) then
                rxik  = rPair(globIndx1, globIndx3)%p%rx
                ryik  = rPair(globIndx1, globIndx3)%p%ry
                rzik  = rPair(globIndx1, globIndx3)%p%rz
                rik   = rPair(globIndx1, globIndx3)%p%r
                if(globIndx3 .gt. globIndx1) then
                  rxik = -rxik
                  ryik = -ryik
                  rzik = -rzik
                endif
                angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
                V1 = gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2) 
                Zeta = Zeta + V1
                Zeta2 = Zeta2 + V1 
              endif
            endif
          enddo
          if(Zeta .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b1 = 1E0_dp
          endif
          if(Zeta2 .ne. 0E0_dp) then
            b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b2 = 1E0_dp
          endif

          V1 = 0.5E0_dp * Fc_Func(rij, R_eq, D2) * (B*exp(-lam2*rij))*(b1 - b2)
          dETable(iIndx) = dETable(iIndx) + V1
          dETable(jIndx) = dETable(jIndx) + V1
          E_Trial = E_Trial + V1
        enddo
      enddo

      E_Trial = E_Trial + E_LJ

      end subroutine
!======================================================================================      
      subroutine NewMol_ECalc_Inter(E_Trial, PairList, dETable)
      use ForceField
      use ForceFieldPara_Tersoff
      use Coords
      use Constants, only: pi
      use SimParameters
      use PairStorage
      implicit none
!      logical, intent(out) :: rejMove
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      
      integer :: iType, jType, kType, nType, iPair
      integer :: iMol, jMol, kMol, nMol
      integer :: iNei, jNei
      integer(kind=atomIntType) :: atmType1      
      integer :: iIndx, jIndx, nIndx 
      integer :: globIndx1, globIndx2, globIndx3, globIndxN
      integer :: neiList(1:60), nNei
!      integer :: pairIndxNew(1:6), nPair
      real(dp) :: r, rMax, rMax_sq

      real(dp) :: A, B, c, d, R_eq, D2 
      real(dp) :: E_LJ, LJ
      real(dp) :: lam1, lam2
      real(dp) :: Zeta, Zeta2
      real(dp) :: BetaPar, n, h
      real(dp) :: b1, b2, V1
      real(dp) :: angijk

      real(dp) :: rxij, ryij, rzij, rij
      real(dp) :: rxjk, ryjk, rzjk, rjk
      real(dp) :: rxik, ryik, rzik, rik

      E_Trial = 0E0_dp
      PairList = 0E0_dp
      dETable = 0E0_dp
      iType = 1
      jType = 1
      kType = 1
      nType = newMol%molType
      nMol = NPART(iType)+1
      nIndx = MolArray(nType)%mol(nMol)%indx
      atmType1 = atomArray(1, 1) 

      A = tersoffData(atmType1)%A
      B = tersoffData(atmType1)%B
      c = tersoffData(atmType1)%c
      d = tersoffData(atmType1)%d
      R_eq = tersoffData(atmType1)%R
      D2 = tersoffData(atmType1)%D2
      BetaPar = tersoffData(atmType1)%beta
      n = tersoffData(atmType1)%n
      h = tersoffData(atmType1)%h
      lam1 = tersoffData(atmType1)%lam1
      lam2 = tersoffData(atmType1)%lam2

      rMax = R_eq + D2
      rMax_sq = rMax * rMax

!       In the Tersoff model, the strength of a given molecular bond is dependent on both the distance and the local environment around the bond.  As a result
!       when one shifts a particle's location one must not only calculate the bonds that changed during this time frame, but also calculate how that
!       shift impacts the bonds of it's neighbors.  Thus it is nessisary to compute a list of particles who may have been implacted. 


!      This block creates a list of neighbors that are located near the particle's new position.
      globIndx1 = MolArray(nType)%mol(nMol)%globalIndx(1)
      globIndxN  = globIndx1
      nNei = 0
      neiList = 0
      E_LJ = 0E0_dp
      do iPair = 1, nNewDist
        globIndx2 = newDist(iPair)%indx2
        jMol = globIndx2
        jIndx = MolArray(jType)%mol(jMol)%indx
        r = newDist(iPair)%r
        if(distCriteria) then                
          PairList(jIndx) = newDist(iPair)%r_sq
        endif
        if(r .lt. rMax) then
          nNei = nNei + 1
          neiList(nNei) = newDist(iPair)%indx2
!          write(*,*) nNei, neiList(nNei)
        endif
        LJ = LJ_Func(newDist(iPair)%r_sq, ep, sig) 
        dETable(nIndx) = dETable(nIndx) + LJ
        dETable(jIndx) = dETable(jIndx) + LJ
        E_LJ = E_LJ + LJ
      enddo


!      This portion of the code calculates the bonds that are created when the new molecule is inserted into the cluster. 
      do jNei = 1, nNei
        jMol = neiList(jNei)
        globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
        jIndx = MolArray(jType)%mol(jMol)%indx
        V1 = 0E0_dp
!        Compute the bonds at the new position
        rij  = rPairNew(globIndxN,globIndx2)%p%r
        if(rij .lt. rMax) then
          Zeta = 0E0_dp
          Zeta2 = 0E0_dp
          rxij = -rPairNew(globIndxN,globIndx2)%p%rx
          ryij = -rPairNew(globIndxN,globIndx2)%p%ry
          rzij = -rPairNew(globIndxN,globIndx2)%p%rz 
          do kMol = 1, NPART(kType)
            if((kMol .eq. nMol) .or. (kMol .eq. jMol)) then
              cycle
            endif
            globIndx3 = MolArray(kType)%mol(kMol)%globalIndx(1)
            rik  = rPairNew(globIndxN,globIndx3)%p%r
            if(rik .lt. rMax) then
              rxik  = -rPairNew(globIndxN,globIndx3)%p%rx
              ryik  = -rPairNew(globIndxN,globIndx3)%p%ry
              rzik  = -rPairNew(globIndxN,globIndx3)%p%rz
              angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
!              write(*,*) 1,angijk * 180d0/pi
              Zeta = Zeta + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2)
            endif

            rjk  = rPair(globIndx2, globIndx3)%p%r
            if(rjk .lt. rMax) then
              rxjk  = rPair(globIndx2, globIndx3)%p%rx
              ryjk  = rPair(globIndx2, globIndx3)%p%ry
              rzjk  = rPair(globIndx2, globIndx3)%p%rz
              if(globIndx3 .gt. globIndx2) then
                rxjk = -rxjk
                ryjk = -ryjk
                rzjk = -rzjk
              endif
              angijk = angleCalc(-rxij, -ryij, -rzij, rij, rxjk, ryjk, rzjk, rjk)
!              write(*,*) 2, angijk * 180d0/pi
              Zeta2 = Zeta2 + gik_Func(angijk, c, d, h) *  Fc_Func(rjk, R_eq, D2)
            endif
          enddo
          if(Zeta .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b1 = 1E0_dp
          endif
          if(Zeta2 .ne. 0E0_dp) then
            b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b2 = 1E0_dp
          endif
!          write(*,*) b1, b2
          V1 = 0.5E0_dp * Fc_Func(rij, R_eq, D2) * (2E0_dp*A*exp(-lam1*rij) - (b1+b2)*B*exp(-lam2*rij))
!          write(*,*) 1, V1
          if(.not. distCriteria) then                
            PairList(jIndx) = PairList(jIndx) + V1
          endif
        endif

!        Compute the bonds at the old position
        dETable(nIndx) = dETable(nIndx) + V1
        dETable(jIndx) = dETable(jIndx) + V1
        E_Trial = E_Trial + V1
      enddo

!      This portion of the code computes the energy of the atoms which neighbored
!      the particle that is being inserted into the system. 
      do iNei = 1, nNei
        iMol = neiList(iNei)
        globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(1)
        iIndx = MolArray(iType)%mol(iMol)%indx

        do jMol = 1, NPART(jType)
          if(iMol .eq. jMol) then
            cycle
          endif

          globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
          rij  = rPair(globIndx1, globIndx2)%p%r
          if(rij .gt. rMax) then
            cycle
          endif

          Zeta = 0E0_dp        !Zeta for the New Configuration
          Zeta2 = 0E0_dp       !Zeta for the Old Configuration
          jIndx = MolArray(jType)%mol(jMol)%indx
 
          rxij = rPair(globIndx1, globIndx2)%p%rx
          ryij = rPair(globIndx1, globIndx2)%p%ry
          rzij = rPair(globIndx1, globIndx2)%p%rz
          if(globIndx2 .gt. globIndx1) then
            rxij = -rxij
            ryij = -ryij
            rzij = -rzij 
          endif

          do kMol = 1, nPart(kType)
            if((kMol .eq. iMol) .or. (kMol .eq. jMol)) then
              cycle
            endif
            globIndx3 = MolArray(kType)%mol(kMol)%globalIndx(1)
            rik  = rPair(globIndx1, globIndx3)%p%r
            if(rik .lt. rMax) then
              rxik  = rPair(globIndx1, globIndx3)%p%rx
              ryik  = rPair(globIndx1, globIndx3)%p%ry
              rzik  = rPair(globIndx1, globIndx3)%p%rz
              rik   = rPair(globIndx1, globIndx3)%p%r
              if(globIndx3 .gt. globIndx1) then
                rxik = -rxik
                ryik = -ryik
                rzik = -rzik
              endif
              angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
!              write(*,*) 3, angijk * 180d0/pi
              V1 = gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2) 
              Zeta = Zeta + V1
              Zeta2 = Zeta2 + V1 
            endif
          enddo

          rik  = rPairNew(globIndxN,globIndx1)%p%r
          if(rik .lt. rMax) then
            rxik  = rPairNew(globIndxN,globIndx1)%p%rx
            ryik  = rPairNew(globIndxN,globIndx1)%p%ry
            rzik  = rPairNew(globIndxN,globIndx1)%p%rz
            angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
!            write(*,*) 4, angijk * 180d0/pi
            Zeta = Zeta + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2) 
          endif

          if(Zeta .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b1 = 1E0_dp
          endif
          if(Zeta2 .ne. 0E0_dp) then
            b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b2 = 1E0_dp
          endif

          V1 = -0.5E0_dp * Fc_Func(rij, R_eq, D2) * (B*exp(-lam2*rij))*(b1 - b2)
!          write(*,*) 2, b1, b2, V1
          dETable(iIndx) = dETable(iIndx) + V1
          dETable(jIndx) = dETable(jIndx) + V1
          E_Trial = E_Trial + V1
        enddo
      enddo


      E_Trial = E_Trial + E_LJ
!      E_Trial = 0.5E0_dp*E_Short
!      write(*,*) "E_Trial", E_Trial
      
      end subroutine    
!======================================================================================      
!      pure subroutine Exchange_ECalc_Inter(E_Trial, nType, nMol, PairList, dETable, rejMove)
!      use ForceField
!      use ForceFieldPara_Tersoff
!      use Coords
!      use SimParameters
!      implicit none
!      logical, intent(out) :: rejMove
!      integer, intent(in) :: nType, nMol
!      real(dp), intent(out) :: E_Trial
!      real(dp), intent(inout) :: PairList(:), dETable(:)
!
!
!      
!      end subroutine    
!======================================================================================      
      subroutine QuickNei_ECalc_Inter_Tersoff(jType, jMol, rejMove)
      implicit none
      integer, intent(in) :: jType, jMol     
      logical, intent(out) :: rejMove
      integer :: dummy

      dummy = jType + jMol
      rejMove = .false.

      
      end subroutine
!======================================================================================
      end module
      
       
