!======================================================================================    

     module q3Functions
       use PairStorage
       use Constants
       private


       logical :: useq3 = .false.
       integer :: q3ArrayIndx

       logical :: initialized = .false.
       logical :: newData = .false.
       integer :: q3Neigh, dNeigh
       integer, parameter :: mOrder = 3
       real(dp) :: q3real(0:2*mOrder) , q3Img(0:2*mOrder)
       real(dp) :: dq3real(0:2*mOrder) , dq3Img(0:2*mOrder)
       real(dp) :: q3Dist, q3DistSq
       real(dp), parameter :: q3Constant = 4d0*pi/(2d0*mOrder+1d0)

       public :: Initialize_q3, Calcq3, useq3, q3Dist, q3DistSq, q3ArrayIndx
       public :: Calcq3_Disp, Calcq3_SwapIn, Calcq3_SwapOut, q3Neigh
       public :: UmbrellaVar_q3

       contains
     !--------------------------------------------------------------------------------
       subroutine Initialize_q3
         use MiscelaniousVars
         use Coords
         use SimParameters, only: nMolTypes, NMAX
         implicit none 
         integer :: startIndx, endIndx
         integer :: gloIndx1, gloIndx2
         integer :: iType, iMol, jType, jMol

         if(.not. useq3) then
           return
         endif
         call ReserveSpace_Coord(1, startIndx, endIndx)
         q3ArrayIndx = startIndx

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
!         q3DistSq = q3Dist*q3Dist
       end subroutine
     !--------------------------------------------------------------------------------
       subroutine UmbrellaVar_q3(iUmbrella,varIndx, biasVar, biasVarNew, outputFormat, &
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
         biasvar(iUmbrella) % realVar => miscCoord(q3ArrayIndx)
         biasvarnew(iUmbrella) % varType = 2
         biasvarnew(iUmbrella) % realVar => miscCoord_New(q3ArrayIndx)
         outputFormat(iUmbrella) = "2x,F12.6,"
!         write(*,*) outputFormat

         iDisp = iDisp + 1
         DispUmbrella(iDisp) % func => Calcq3_Disp
         iSwapIn = iSwapIn + 1
         SwapInUmbrella(iSwapIn) % func => Calcq3_SwapIn
         iSwapOut = iSwapOut + 1
         SwapOutUmbrella(iSwapOut) % func => Calcq3_SwapOut
       end subroutine
     !--------------------------------------------------------------------------------
       subroutine Calcq3
         use MiscelaniousVars
         use Coords
         use SimParameters, only: nMolTypes, NPART, NTotal, prevMoveAccepted
         implicit none 
         integer :: m
         integer :: iType, iMol, jType, jMol, jMolMin
         integer :: gloIndx1, gloIndx2
         integer :: jIndx, iIndx
         real(dp) :: rx,ry,rz,r
         real(dp) :: theta,phi,junk,sTerm,cTerm
         real(dp) :: q3Par
         real(dp) :: rPart,iPart


  
         if(NTotal .eq. 1) then
           miscCoord(q3ArrayIndx) = 1E0_dp
           q3real = 0E0_dp
           q3img = 0E0_dp
           q3neigh = 0
           newData = .false.

           return  
         endif

         if(prevMoveAccepted) then
           if( newData ) then
             do m = 0, 2*mOrder
               q3real(m) = q3real(m) + dq3Real(m)
               q3img(m) = q3img(m) + dq3img(m)
             enddo
             q3Neigh = q3Neigh + dNeigh
             miscCoord(q3ArrayIndx) = miscCoord_New(q3ArrayIndx)
             newData = .false.
             return
           endif
         endif

         if(initialized) then
           return
         endif

         q3Neigh = 0
         q3real = 0E0_dp
         q3img = 0E0_dp

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

                 if(r .le. q3Dist) then	
                   q3Neigh = q3Neigh + 1
                   rx = rPair(gloIndx2, gloIndx1)%p%rx
                   ry = rPair(gloIndx2, gloIndx1)%p%ry
                   rz = rPair(gloIndx2, gloIndx1)%p%rz
!                   if(gloIndx2 .lt. gloIndx1) then
!                     rPart = -rPart
!                     iPart = -iPart
!                   endif
                   phi = atan2(ry,rx)
                   theta = acos(rz/r)   
                   do m = 0, 2*mOrder
                     call Harmonics(theta, phi, m-mOrder, rPart, iPart)	
                     q3real(m) = q3real(m) + rPart
                     q3img(m) = q3img(m) + iPart	
                   enddo
                 endif
               enddo
             enddo
           enddo
         enddo

         q3par = 0E0_dp
         do m = 0, 2*mOrder
!           write(2,*) m-mOrder, q3real(m), q3img(m)
           q3Par = q3par + q3real(m)*q3real(m) + q3img(m)*q3img(m)
         enddo
!         flush(2)
         q3par = q3par * q3Constant
         q3par = sqrt(q3par)/q3Neigh
 	       miscCoord(q3ArrayIndx) = q3par
!         if(prevMoveAccepted) then
!           if(abs(miscCoord(q3ArrayIndx) - miscCoord_New(q3ArrayIndx)) .gt. 1d-7) then
!             write(2,*) "DISCRPANCY:", miscCoord(q3ArrayIndx), miscCoord_New(q3ArrayIndx)
!           endif
!           do m = 0, 12
!             if(abs(q3real(m) - q3r_t(m)) .gt. 1d-7) then
!               write(2,*) "DISCRPANCY:",m, q3real(m) , q3r_t(m)
!             endif
!           enddo
!           if(q3Neigh .ne. q3Neigh_t) then
!             write(2,*) "DISCRPANCY:", q3Neigh, q3Neigh_t
!           endif
!         endif
         initialized = .true.
      end subroutine
     !--------------------------------------------------------------------------------
       subroutine Calcq3_Disp(disp)
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
         real(dp) :: q3r_new, q3i_new
         real(dp) :: rx,ry,rz,r
         real(dp) :: theta,phi,junk,sTerm,cTerm
         real(dp) :: q3Par
         real(dp) :: rPart,iPart

         dNeigh = 0
         dq3real = 0E0_dp
         dq3img = 0E0_dp
         if(NTotal .eq. 1) then
 	         miscCoord_New(q3ArrayIndx) = 1E0_dp
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
           miscCoord_New(q3ArrayIndx) = miscCoord(q3ArrayIndx)
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
!             do m = 0, 2*mOrder
!               write(2,*) "total",m-mOrder, q3real(m), q3img(m)
!             enddo    
             rx = disp(dispAtom1)%x_new  - molArray(jType)%mol(jMol)%x(1) 
             ry = disp(dispAtom1)%y_new  - molArray(jType)%mol(jMol)%y(1) 
             rz = disp(dispAtom1)%z_new  - molArray(jType)%mol(jMol)%z(1)
             r = rx*rx + ry*ry + rz*rz
             
             if(r .le. q3DistSq) then	
               if(gloIndx2 .lt. gloIndx1) then
                 rx = -rx
                 ry = -ry
                 rz = -rz
               endif
               r = sqrt(r)
               dNeigh = dNeigh + 1
               phi = atan2( ry, rx )
               theta = acos( rz / r )  
!               write(2,*) rx, ry, rz, r
!               write(2,*) phi, theta
               do m = 0, 2*mOrder
                 call Harmonics(theta, phi, m-mOrder, rPart, iPart)	
!                 if(gloIndx2 .lt. gloIndx1) then
!                   rPart = -rPart
!                   iPart = -iPart
!                 endif
!                 write(2,*) "new",m-mOrder, rPart, iPart
                 dq3real(m) = dq3real(m) + rPart
                 dq3img(m) = dq3img(m) + iPart	
               enddo               
             endif
             

             if(rPair(gloIndx2, gloIndx1)%p%r_sq .le. q3DistSq) then	
                dNeigh = dNeigh - 1
!                rx = rPair(gloIndx2, gloIndx1)%p%rx
!                ry = rPair(gloIndx2, gloIndx1)%p%ry
!                rz = rPair(gloIndx2, gloIndx1)%p%rz
   
                phi = atan2(rPair(gloIndx2, gloIndx1)%p%ry, rPair(gloIndx2, gloIndx1)%p%rx)
                theta = acos(rPair(gloIndx2, gloIndx1)%p%rz / rPair(gloIndx2, gloIndx1)%p%r )   
!                write(2,*) phi, theta
                do m = 0, 2*mOrder
                  call Harmonics(theta, phi, m-mOrder, rPart, iPart)	
!                  if(gloIndx2 .lt. gloIndx1) then
!                    rPart = -rPart
!                    iPart = -iPart
!                  endif
!                  write(2,*) "old",m-mOrder, rPart, iPart
                  dq3real(m) = dq3real(m) - rPart
                  dq3img(m) = dq3img(m) - iPart	
                enddo
             endif   
           enddo
         enddo
          
!         write(2,*) q3Neigh, dNeigh
         q3par = 0E0_dp
         do m = 0, 2*mOrder
!           write(2,*) m-mOrder, q3real(m) + dq3real(m), q3img(m) + dq3img(m)
           q3r_new = q3real(m) + dq3real(m)
           q3i_new = q3img(m) + dq3img(m)
           q3Par = q3par + q3r_new*q3r_new + q3i_new*q3i_new   	   
         enddo
!         write(2,*)
!         flush(2)
         q3par = q3par * q3Constant
         newNeigh = q3Neigh + dNeigh
         q3par = sqrt(q3par)/newNeigh
!         write(2,*) "q3", q3par
 	       miscCoord_New(q3ArrayIndx) = q3par
         newData = .true.
      end subroutine
     !--------------------------------------------------------------------------------
       subroutine Calcq3_SwapIn
         use MiscelaniousVars
         use Coords
         use SimParameters, only: nMolTypes, NPART, NTotal
         implicit none 
         integer :: m
         integer :: iType, iMol, jType, jMol
         integer :: gloIndx1, gloIndx2
         integer :: jIndx, iIndx, sizeDisp
         integer :: dispAtom1
         real(dp) :: q3r_new, q3i_new
         real(dp) :: rx,ry,rz,r
         real(dp) :: theta,phi,junk,sTerm,cTerm
         real(dp) :: q3Par
         real(dp) :: rPart,iPart

 
         iType = newMol%molType
         iMol = NPART(iType) + 1
         gloIndx1 = molArray(iType)%mol(iMol)%globalIndx(1)
         dNeigh = 0
!         gloIndx1 =  molArray(iType)%mol(iMol)%globalIndx(1)
         dq3real = 0E0_dp
         dq3img = 0E0_dp
         do jType = 1, nMolTypes
           do jMol = 1, NPART(jType)
             gloIndx2 = molArray(iType)%mol(iMol)%globalIndx(1)
             rx = newMol%x(1) -  molArray(jType)%mol(jMol)%x(1)
             ry = newMol%y(1) -  molArray(jType)%mol(jMol)%y(1)
             rz = newMol%z(1) -  molArray(jType)%mol(jMol)%z(1)
             r = rx*rx + ry*ry + rz*rz
             if(r .le. q3DistSq) then	
               if(gloIndx2 .lt. gloIndx1) then
                 rx = -rx
                 ry = -ry
                 rz = -rz
               endif
               r = sqrt(r)
               dNeigh = dNeigh + 1
               phi = atan2(ry,rx)
               theta = acos(rz/r)  
               do m = 0, 2*mOrder
                 call Harmonics(theta, phi, m-mOrder, rPart, iPart)	

                 dq3real(m) = dq3real(m) + rPart
                 dq3img(m) = dq3img(m) + iPart	
               enddo               
             endif
           enddo
         enddo
          
         q3par = 0E0_dp
         do m = 0, 2*mOrder
           q3r_new = q3real(m) + dq3real(m)
           q3i_new = q3img(m) + dq3img(m)
           q3Par = q3par + q3r_new*q3r_new + q3i_new*q3i_new   	   
         enddo
         q3par = q3par * q3Constant
         q3par = sqrt(q3par)/real((q3Neigh + dNeigh), dp)
 	       miscCoord_New(q3ArrayIndx) = q3par
         newData = .true.
      end subroutine
     !--------------------------------------------------------------------------------
       subroutine Calcq3_SwapOut(nType, nMol)
         use MiscelaniousVars
         use Coords
         use SimParameters
         implicit none 
         integer, intent(in) :: nType, nMol
         integer :: m
         integer :: i, iType, iMol, jType, jMol
         integer :: gloIndx1, gloIndx2
         integer :: jIndx, iIndx, sizeDisp
         integer :: dispAtom1
         real(dp) :: q3r_new, q3i_new
         real(dp) :: rx,ry,rz,r
         real(dp) :: theta,phi,junk,sTerm,cTerm
         real(dp) :: q3Par
         real(dp) :: rPart,iPart

         if(NTotal .eq. 2) then
!           dNeigh = -q3Neigh
!           dq3real = -q3real
!           dq3real = -q3img
           miscCoord_New(q3ArrayIndx) = 1E0
           return
         endif

 
         dNeigh = 0
         dq3real = 0E0
         dq3img = 0E0
         gloIndx1 =  molArray(nType)%mol(nMol)%globalIndx(1)
!         iIndx =  molArray(nType)%mol(nMol)%indx
         do jType = 1, nMolTypes
           do jMol = 1, NPART(jType)
!             jIndx =  molArray(jType)%mol(jMol)%indx
             gloIndx2 =  molArray(jType)%mol(jMol)%globalIndx(1)
             if(gloIndx1 .eq. gloIndx2) then
               cycle
             endif
             if(rPair(gloIndx2, gloIndx1)%p%r_sq .le. q3DistSq) then
               r = rPair(gloIndx2, gloIndx1)%p%r
               dNeigh = dNeigh - 1
               rx = rPair(gloIndx2, gloIndx1)%p%rx
               ry = rPair(gloIndx2, gloIndx1)%p%ry
               rz = rPair(gloIndx2, gloIndx1)%p%rz
               phi = atan2(ry,rx)
               theta = acos(rz/r)  
               do m = 0, 2*mOrder
                 call Harmonics(theta, phi, m-mOrder, rPart, iPart)
                 if(gloIndx2 .lt. gloIndx1) then
                   rPart = -rPart
                   iPart = -iPart
                 endif
                 dq3real(m) = dq3real(m) - rPart
                 dq3img(m) = dq3img(m) - iPart
               enddo               
             endif
           enddo
         enddo
          
         q3par = 0E0
         do m = 0, 2*mOrder
           q3r_new = q3real(m) + dq3real(m)
           q3i_new = q3img(m) + dq3img(m)
           q3Par = q3par + q3r_new*q3r_new + q3i_new*q3i_new   
         enddo
         q3par = q3par * q3Constant
         q3par = sqrt(q3par)/real((q3Neigh + dNeigh), dp)
         miscCoord_New(q3ArrayIndx) = q3par
         newData = .true.
      end subroutine
     !--------------------------------------------------------------------------------
      subroutine Harmonics(Theta, Phi, m, RealVal, ImgVal)
      use Constants	
      use VarPrecision  
      implicit none	
      integer, intent(in) :: m
      real(dp), intent(in) :: Theta,Phi
      real(dp), intent(out) :: RealVal, ImgVal

      real(dp) :: cTerm, sTerm, cTermsq, sTermsq
      real(dp), parameter :: const0 = (1d0/4d0)*sqrt(7d0/pi)
      real(dp), parameter :: const1 = (1d0/8d0)*sqrt(21d0/pi)
      real(dp), parameter :: const2 = (1d0/4d0)*sqrt(105d0/(2d0*pi))
      real(dp), parameter :: const3 = (1d0/8d0)*sqrt(35d0/pi)

      RealVal = 0E0_dp
      ImgVal = 0E0_dp
      cterm = cos(Theta)
      sterm = sin(Theta)
  
      select case (abs(m))
      case(0)
        ctermsq = cterm * cterm
        RealVal = 5d0*ctermsq*cterm - 3d0*cterm
        RealVal = RealVal * const0
      case(1) 
        ctermsq = cterm * cterm
        RealVal = sterm * (5d0*ctermsq - 1d0)
        if(m .gt. 0) then
          RealVal = -RealVal
        endif
        RealVal = RealVal * const1
      case(2)
        RealVal = sterm * sterm * cterm
        RealVal = RealVal * const2
      case(3)
        RealVal = sTerm * sTerm * sTerm
        if(m .gt. 0) then
          RealVal = -RealVal
        endif
        RealVal = RealVal * const3
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



