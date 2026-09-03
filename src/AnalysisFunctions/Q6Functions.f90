!======================================================================================    
!      This module contains the functions nessisary to initalize and 
!      perform single pair calculations.  This can be used
!      in the Umbrella Sampling algorithms or be used to study
!      di
     module Q6Functions
       use PairStorage
       use Constants
       private


       logical :: useQ6 = .false.
       integer :: q6ArrayIndx
!       integer :: analysisIndx

       logical :: initialized = .false.
       logical :: newData = .false.
       integer :: q6Neigh, dNeigh
       real(dp) :: Q6real(0:12) , Q6Img(0:12)
       real(dp) :: dQ6real(0:12) , dQ6Img(0:12)
!       real(dp) :: Q6r_t(0:12) , Q6i_t(0:12)
       real(dp) :: q6Dist, q6DistSq
       real(dp), parameter :: q6Constant = 4d0*pi/13d0



       public :: Initialize_q6, CalcQ6, useQ6, q6Dist, q6DistSq, q6ArrayIndx
       public :: CalcQ6_Disp, CalcQ6_SwapIn, CalcQ6_SwapOut, q6Neigh
       public :: UmbrellaVar_Q6
       contains
     !--------------------------------------------------------------------------------
       subroutine Initialize_q6
         use MiscelaniousVars
         use Coords
         use SimParameters, only: nMolTypes, NMAX
         implicit none 
         integer :: startIndx, endIndx
         integer :: gloIndx1, gloIndx2
         integer :: iType, iMol, jType, jMol

         if(.not. useQ6) then
           return
         endif
         call ReserveSpace_Coord(1, startIndx, endIndx)
         q6ArrayIndx = startIndx

         do iType = 1, nMolTypes
           do iMol = 1, NMAX(iType)
             gloIndx1 = molArray(iType)%mol(iMol)%globalIndx(1)
             do jType = 1, nMolTypes
               do jMol = 1, NMAX(jType)
                 gloIndx2 = molArray(jType)%mol(jMol)%globalIndx(1)
                 if(gloIndx1 .ne. gloIndx2) then
                   rPair(gloIndx1,gloIndx2) % p % storeRParts = .true.
                   rPair(gloIndx1,gloIndx2) % p % storeRValue = .true.
                 endif
               enddo
             enddo
           enddo
         enddo 
!         q6DistSq = q6Dist*q6Dist

       end subroutine
     !--------------------------------------------------------------------------------
       subroutine UmbrellaVar_Q6(iUmbrella,varIndx, biasVar, biasVarNew, outputFormat, &
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


         biasvar(iUmbrella) % varType = 2
         biasvar(iUmbrella) % realVar => miscCoord(q6ArrayIndx)
         biasvarnew(iUmbrella) % varType = 2
         biasvarnew(iUmbrella) % realVar => miscCoord_New(q6ArrayIndx)
         outputFormat(iUmbrella) = "2x,F12.6,"
!         write(*,*) outputFormat

         iDisp = iDisp + 1
         DispUmbrella(iDisp) % func => CalcQ6_Disp
         iSwapIn = iSwapIn + 1
         SwapInUmbrella(iSwapIn) % func => CalcQ6_SwapIn
         iSwapOut = iSwapOut + 1
         SwapOutUmbrella(iSwapOut) % func => CalcQ6_SwapOut
       end subroutine
     !--------------------------------------------------------------------------------
       subroutine CalcQ6
         use MiscelaniousVars
         use Coords
         use SimParameters, only: nMolTypes, NPART, NTotal, prevMoveAccepted
         implicit none 
         integer :: m
         integer :: iType, iMol, jType, jMol, jMolMin
         integer :: gloIndx1, gloIndx2
         real(dp) :: rx,ry,rz,r
         real(dp) :: theta,phi, cTerm
         real(dp) :: q6Par
         real(dp) :: rPart,iPart

  
         if(NTotal .eq. 1) then
           miscCoord(q6ArrayIndx) = 1E0_dp
           Q6real = 0E0_dp
           Q6img = 0E0_dp
           q6neigh = 0
           newData = .false.
           return  
         endif

         if(prevMoveAccepted) then
           if( newData ) then
             do m = 0,12
               Q6real(m) = Q6real(m) + dQ6Real(m)
               Q6img(m) = Q6img(m) + dQ6img(m)
             enddo
             q6Neigh = q6Neigh + dNeigh
             miscCoord(q6ArrayIndx) = miscCoord_New(q6ArrayIndx)
             newData = .false.
             return
           endif
         endif

         if(initialized) then
           return
         endif

         q6Neigh = 0
         Q6real = 0E0_dp
         Q6img = 0E0_dp

         do iType = 1, nMolTypes
           do iMol = 1, NPART(iType)
             gloIndx1 = molArray(iType)%mol(iMol)%globalIndx(1)
             do jType = 1, nMolTypes
               if(iType .eq. jType) then
                 jMolMin = iMol+1
               else
                 jMolMin = 1        
               endif
               do jMol = jMolMin, NPART(jType)
                 gloIndx2 = molArray(jType)%mol(jMol)%globalIndx(1)
                 if(gloIndx1 .eq. gloIndx2) then
                   cycle
                 endif
                 r = rPair(gloIndx2, gloIndx1)%p%r

                 if(r .le. q6Dist) then
                   q6Neigh = q6Neigh + 1
                   rx = rPair(gloIndx2, gloIndx1)%p%rx
                   ry = rPair(gloIndx2, gloIndx1)%p%ry
                   rz = rPair(gloIndx2, gloIndx1)%p%rz
                   phi = atan2(ry,rx)
                   theta = acos(rz/r)   
                   do m = 0, 12
                     call Harmonics(theta, phi, m-6, rPart, iPart)
                     Q6real(m) = Q6real(m) + rPart
                     Q6img(m) = Q6img(m) + iPart
                   enddo
                 endif
               enddo
             enddo
           enddo
         enddo

         q6par = 0E0_dp
         do m = 0, 12
           q6Par = q6par + Q6real(m)*Q6real(m) + Q6img(m)*Q6img(m) 
         enddo
         q6par = q6par * q6Constant
         q6par = sqrt(q6par)/q6Neigh
         miscCoord(q6ArrayIndx) = q6par
!         if(prevMoveAccepted) then
!           if(abs(miscCoord(q6ArrayIndx) - miscCoord_New(q6ArrayIndx)) .gt. 1d-7) then
!             write(2,*) "DISCRPANCY:", miscCoord(q6ArrayIndx), miscCoord_New(q6ArrayIndx)
!           endif
!           do m = 0, 12
!             if(abs(Q6real(m) - Q6r_t(m)) .gt. 1d-7) then
!               write(2,*) "DISCRPANCY:",m, Q6real(m) , Q6r_t(m)
!             endif
!           enddo
!           if(q6Neigh .ne. q6Neigh_t) then
!             write(2,*) "DISCRPANCY:", q6Neigh, q6Neigh_t
!           endif
!         endif
         initialized = .true.
      end subroutine
     !--------------------------------------------------------------------------------
       subroutine CalcQ6_Disp(disp)
         use MiscelaniousVars
         use Coords
         use SimParameters, only: nMolTypes, NPART, NTotal
         implicit none 
         type(Displacement), intent(in) :: disp(:)

         logical :: changed
         integer :: m
         integer :: iDisp, iType, iMol, jType, jMol
         integer :: gloIndx1, gloIndx2
         integer :: sizeDisp
         integer :: dispAtom1, newNeigh
         real(dp) :: q6r_new, q6i_new
         real(dp) :: rx,ry,rz,r
         real(dp) :: theta,phi,cTerm
         real(dp) :: q6Par
         real(dp) :: rPart,iPart

         dNeigh = 0
         dQ6real = 0E0_dp
         dQ6img = 0E0_dp
         if(NTotal .eq. 1) then
           miscCoord_New(q6ArrayIndx) = 1E0_dp
           newData = .false.
           return
         endif

         sizeDisp = size(disp)
         changed = .false.
         do iDisp = 1, sizeDisp
           if(disp(iDisp)%atmIndx .eq. 1) then
             changed = .true.
             dispAtom1 = iDisp
             exit
           endif
         enddo

         if(.not. changed) then
           miscCoord_New(q6ArrayIndx) = miscCoord(q6ArrayIndx)
           newData = .false.
           return
         endif

         iType = disp(1)%molType
         iMol = disp(1)%molIndx

         gloIndx1 =  molArray(iType)%mol(iMol)%globalIndx(1)
         do jType = 1, nMolTypes
           do jMol = 1, NPART(jType)
             gloIndx2 =  molArray(jType)%mol(jMol)%globalIndx(1)
             if(gloIndx2 .eq. gloIndx1) then
               cycle
             endif
             rx = molArray(jType)%mol(jMol)%x(1) - disp(dispAtom1)%x_new 
             ry = molArray(jType)%mol(jMol)%y(1) - disp(dispAtom1)%y_new 
             rz = molArray(jType)%mol(jMol)%z(1) - disp(dispAtom1)%z_new
             r = rx*rx + ry*ry + rz*rz
             if(r .le. q6DistSq) then
               r = sqrt(r)
               dNeigh = dNeigh + 1
               phi = atan2( ry, rx )
               theta = acos( rz / r )  
               do m = 0, 12
                 call Harmonics(theta, phi, m-6, rPart, iPart)
                 dQ6real(m) = dQ6real(m) + rPart
                 dQ6img(m) = dQ6img(m) + iPart
               enddo               
             endif
             

             if(rPair(gloIndx2, gloIndx1)%p%r_sq .le. q6DistSq) then
                dNeigh = dNeigh - 1
!                rx = rPair(gloIndx2, gloIndx1)%p%rx
!                ry = rPair(gloIndx2, gloIndx1)%p%ry
!                rz = rPair(gloIndx2, gloIndx1)%p%rz
                phi = atan2(rPair(gloIndx2, gloIndx1)%p%ry, rPair(gloIndx2, gloIndx1)%p%rx)
                theta = acos(rPair(gloIndx2, gloIndx1)%p%rz / rPair(gloIndx2, gloIndx1)%p%r )   
                do m = 0, 12
                  call Harmonics(theta, phi, m-6, rPart, iPart)
                  dQ6real(m) = dQ6real(m) - rPart
                  dQ6img(m) = dQ6img(m) - iPart
                enddo
             endif   
           enddo
         enddo
          
         q6par = 0E0_dp
         do m = 0, 12
           q6r_new = Q6real(m) + dQ6real(m)
           q6i_new = Q6img(m) + dQ6img(m)
           q6Par = q6par + q6r_new*q6r_new + q6i_new*q6i_new
         enddo
         q6par = q6par * q6Constant
         newNeigh = q6Neigh + dNeigh
         q6par = sqrt(q6par)/newNeigh
         miscCoord_New(q6ArrayIndx) = q6par
         newData = .true.
      end subroutine
     !--------------------------------------------------------------------------------
       subroutine CalcQ6_SwapIn
         use MiscelaniousVars
         use Coords
         use SimParameters, only: nMolTypes, NPART, NTotal
         implicit none 
         integer :: m
         integer :: iType, iMol, jType, jMol
         integer :: gloIndx1, gloIndx2
         integer :: jIndx, iIndx, sizeDisp
         integer :: dispAtom1
         real(dp) :: q6r_new, q6i_new
         real(dp) :: rx,ry,rz,r
         real(dp) :: theta,phi,junk,sTerm,cTerm
         real(dp) :: q6Par
         real(dp) :: rPart,iPart

 
         iType = newMol%molType
         iMol = NPART(iType) + 1
         dNeigh = 0
!         gloIndx1 =  molArray(iType)%mol(iMol)%globalIndx(1)
         dQ6real = 0E0_dp
         dQ6img = 0E0_dp
         do jType = 1, nMolTypes
           do jMol = 1, NPART(jType)
             rx = newMol%x(1) -  molArray(jType)%mol(jMol)%x(1)
             ry = newMol%y(1) -  molArray(jType)%mol(jMol)%y(1)
             rz = newMol%z(1) -  molArray(jType)%mol(jMol)%z(1)
             r = rx*rx + ry*ry + rz*rz
             if(r .le. q6DistSq) then
               r = sqrt(r)
               dNeigh = dNeigh + 1
               phi = atan2(ry,rx)
               theta = acos(rz/r)  
               do m = 0, 12
                 call Harmonics(theta, phi, m-6, rPart, iPart)
                 dQ6real(m) = dQ6real(m) + rPart
                 dQ6img(m) = dQ6img(m) + iPart
               enddo               
             endif
           enddo
         enddo
          
         q6par = 0E0_dp
         do m = 0, 12
           q6r_new = Q6real(m) + dQ6real(m)
           q6i_new = Q6img(m) + dQ6img(m)
           q6Par = q6par + q6r_new*q6r_new + q6i_new*q6i_new
         enddo
         q6par = q6par * q6Constant
         q6par = sqrt(q6par)/real((q6Neigh + dNeigh), dp)
         miscCoord_New(q6ArrayIndx) = q6par
         newData = .true.
      end subroutine
     !--------------------------------------------------------------------------------
       subroutine CalcQ6_SwapOut(nType, nMol)
         use MiscelaniousVars
         use Coords
         use SimParameters
         implicit none 
         integer, intent(in) :: nType, nMol
         integer :: m
         integer :: jType, jMol
         integer :: gloIndx1, gloIndx2
         integer :: jIndx, sizeDisp
         real(dp) :: q6r_new, q6i_new
         real(dp) :: rx,ry,rz,r
         real(dp) :: theta,phi,junk
         real(dp) :: q6Par
         real(dp) :: rPart,iPart

         if(NTotal .eq. 2) then
!           dNeigh = -q6Neigh
!           dQ6real = -Q6real
!           dQ6real = -Q6img
            miscCoord_New(q6ArrayIndx) = 1E0
           return
         endif

 
         dNeigh = 0
         dQ6real = 0E0
         dQ6img = 0E0
         gloIndx1 =  molArray(nType)%mol(nMol)%globalIndx(1)
!         iIndx =  molArray(nType)%mol(nMol)%indx
         do jType = 1, nMolTypes
           do jMol = 1, NPART(jType)
!             jIndx =  molArray(jType)%mol(jMol)%indx
             gloIndx2 =  molArray(jType)%mol(jMol)%globalIndx(1)
             if(gloIndx1 .eq. gloIndx2) then
               cycle
             endif
             if(rPair(gloIndx2, gloIndx1)%p%r_sq .le. q6DistSq) then
               r = rPair(gloIndx2, gloIndx1)%p%r
               dNeigh = dNeigh - 1
               rx = rPair(gloIndx2, gloIndx1)%p%rx
               ry = rPair(gloIndx2, gloIndx1)%p%ry
               rz = rPair(gloIndx2, gloIndx1)%p%rz
               phi = atan2(ry,rx)
               theta = acos(rz/r)  
               do m = 0, 12
                 call Harmonics(theta, phi, m-6, rPart, iPart)
                 dQ6real(m) = dQ6real(m) - rPart
                 dQ6img(m) = dQ6img(m) - iPart
               enddo               
             endif
           enddo
         enddo
          
         q6par = 0E0
         do m = 0, 12
           q6r_new = Q6real(m) + dQ6real(m)
           q6i_new = Q6img(m) + dQ6img(m)
           q6Par = q6par + q6r_new*q6r_new + q6i_new*q6i_new
         enddo
         q6par = q6par * q6Constant
         q6par = sqrt(q6par)/real((q6Neigh + dNeigh), dp)
         miscCoord_New(q6ArrayIndx) = q6par
         newData = .true.
      end subroutine
     !--------------------------------------------------------------------------------
!       subroutine Q6_Update
!         use MiscelaniousVars
!         use Coords
!         use SimParameters
!         implicit none 


!      end subroutine
     !--------------------------------------------------------------------------------
      subroutine Harmonics(Theta, Phi, m, RealVal, ImgVal)
      use Constants
      use VarPrecision  
      implicit none
      integer, intent(in) :: m
      real(dp), intent(in) :: Theta,Phi
      real(dp), intent(out) :: RealVal, ImgVal

      real(dp) :: cTerm, sTerm, cTermsq, sTermsq
      real(dp), parameter :: const0 = sqrt(13d0/pi)/32d0
      real(dp), parameter :: const1 = sqrt(273d0/(2d0*pi))/16d0
      real(dp), parameter :: const2 = sqrt(1365d0/pi)/64d0
      real(dp), parameter :: const3 = sqrt(1365d0/pi)/32d0
      real(dp), parameter :: const4 = sqrt(91d0/(2d0*pi))*3d0/32d0
      real(dp), parameter :: const5 = sqrt(1001d0/pi)*3d0/32d0
      real(dp), parameter :: const6 = sqrt(3003d0/pi)/64d0

      RealVal = 0E0_dp
      ImgVal = 0E0_dp
      cterm = cos(Theta)
      sterm = sin(Theta)

      select case (abs(m))
      case(0)
        ctermsq = cterm * cterm
        RealVal = 231d0*ctermsq*ctermsq*ctermsq - 315d0*ctermsq*ctermsq
        RealVal = RealVal + 105d0*ctermsq - 5d0
        RealVal = RealVal * const0
      case(1) 
        ctermsq = cterm * cterm
        RealVal = 33d0*ctermsq*ctermsq*cterm
        RealVal = RealVal - 30d0*cterm*ctermsq + 5d0*cterm 
        RealVal = sterm * RealVal 
        RealVal = RealVal * const1
        if(m .gt. 0) then
          RealVal = -RealVal
        endif
      case(2)
        ctermsq = cterm * cterm
        RealVal = sterm * sterm
        RealVal = RealVal * (33d0*ctermsq*ctermsq - 18d0*ctermsq + 1d0)
        RealVal = RealVal * const2
      case(3)
        RealVal = sTerm * sTerm * sTerm
        RealVal = RealVal * cterm * (11d0*cterm*cterm-3d0)
        RealVal = RealVal * const3
        if(m .gt. 0) then
          RealVal = -RealVal
        endif
      case(4)
        sTermsq = sterm * sterm
        RealVal = sTermsq * sTermsq * (11d0*cterm*cterm - 1d0)
        RealVal = RealVal * const4
      case(5)
        sTermsq = sTerm * sTerm
        RealVal = sTerm * sTermsq * sTermsq * cterm
        RealVal = RealVal * const5
        if(m .gt. 0) then
          RealVal = -RealVal
        endif
      case(6)
        sTermsq = sTerm * sTerm
        RealVal = sTermsq * sTermsq *sTermsq
        RealVal = RealVal * const6
      end select

      ImgVal = RealVal * sin(m * Phi)
      RealVal = RealVal * cos(m * Phi)


      If(Abs(ImgVal) .lt. 1d-13) then
        ImgVal = 0E0_dp
      endif

      If(Abs(RealVal) .lt. 1d-13) then
        RealVal = 0E0_dp
      endif

      end subroutine
     !--------------------------------------------------------------------------------
      end module
!====================================================================================== 



