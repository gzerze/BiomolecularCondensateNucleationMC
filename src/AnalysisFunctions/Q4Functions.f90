!======================================================================================    

     module q4Functions
       use PairStorage
       use Constants
       private


       logical :: useq4 = .false.
       integer :: q4ArrayIndx

       logical :: initialized = .false.
       logical :: newData = .false.
       integer :: q4Neigh, dNeigh
       real(dp) :: q4real(0:8) , q4Img(0:8)
       real(dp) :: dq4real(0:8) , dq4Img(0:8)
       real(dp) :: q4Dist, q4DistSq
       real(dp), parameter :: q4Constant = 4d0*pi/(2d0*4d0+1d0)

       public :: Initialize_q4, Calcq4, useq4, q4Dist, q4DistSq, q4ArrayIndx
       public :: Calcq4_Disp, Calcq4_SwapIn, Calcq4_SwapOut, q4Neigh
       public :: UmbrellaVar_Q4

       contains
     !--------------------------------------------------------------------------------
       subroutine Initialize_q4
         use MiscelaniousVars
         use Coords
         use SimParameters, only: nMolTypes, NMAX
         implicit none 
         integer :: startIndx, endIndx
         integer :: gloIndx1, gloIndx2
         integer :: iType, iMol, jType, jMol

         if(.not. useQ4) then
           return
         endif
         call ReserveSpace_Coord(1, startIndx, endIndx)
         Q4ArrayIndx = startIndx

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
!         Q4DistSq = Q4Dist*Q4Dist
       end subroutine
     !--------------------------------------------------------------------------------
       subroutine UmbrellaVar_Q4(iUmbrella,varIndx, biasVar, biasVarNew, outputFormat, &
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
         biasvar(iUmbrella) % realVar => miscCoord(q4ArrayIndx)
         biasvarnew(iUmbrella) % varType = 2
         biasvarnew(iUmbrella) % realVar => miscCoord_New(q4ArrayIndx)
         outputFormat(iUmbrella) = "2x,F12.6,"
!         write(*,*) outputFormat

         iDisp = iDisp + 1
         DispUmbrella(iDisp) % func => CalcQ4_Disp
         iSwapIn = iSwapIn + 1
         SwapInUmbrella(iSwapIn) % func => CalcQ4_SwapIn
         iSwapOut = iSwapOut + 1
         SwapOutUmbrella(iSwapOut) % func => CalcQ4_SwapOut
       end subroutine
     !--------------------------------------------------------------------------------
       subroutine Calcq4
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
         real(dp) :: q4Par
         real(dp) :: rPart,iPart


  
         if(NTotal .eq. 1) then
 	       miscCoord(q4ArrayIndx) = 1E0_dp
           q4real = 0E0_dp
           q4img = 0E0_dp
           q4neigh = 0
           newData = .false.
           return  
         endif

         if(prevMoveAccepted) then
           if( newData ) then
             do m = 0,8
               q4real(m) = q4real(m) + dq4Real(m)
               q4img(m) = q4img(m) + dq4img(m)
             enddo
             q4Neigh = q4Neigh + dNeigh
             miscCoord(q4ArrayIndx) = miscCoord_New(q4ArrayIndx)
             newData = .false.
             return
           endif
         endif

         if(initialized) then
           return
         endif

         q4Neigh = 0
         q4real = 0E0_dp
         q4img = 0E0_dp

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

                 if(r .le. q4Dist) then	
                   q4Neigh = q4Neigh + 1
                   rx = rPair(gloIndx2, gloIndx1)%p%rx
                   ry = rPair(gloIndx2, gloIndx1)%p%ry
                   rz = rPair(gloIndx2, gloIndx1)%p%rz
                   phi = atan2(ry,rx)
                   theta = acos(rz/r)   
                   do m = 0, 8
                     call Harmonics(theta, phi, m-4, rPart, iPart)	
                     q4real(m) = q4real(m) + rPart
                     q4img(m) = q4img(m) + iPart	
                   enddo
                 endif
               enddo
             enddo
           enddo
         enddo

         q4par = 0E0_dp
         do m = 0, 8
           q4Par = q4par + q4real(m)*q4real(m) + q4img(m)*q4img(m)	  	   
         enddo
         q4par = q4par * q4Constant
         q4par = sqrt(q4par)/q4Neigh
 	     miscCoord(q4ArrayIndx) = q4par
!         if(prevMoveAccepted) then
!           if(abs(miscCoord(q4ArrayIndx) - miscCoord_New(q4ArrayIndx)) .gt. 1d-7) then
!             write(2,*) "DISCRPANCY:", miscCoord(q4ArrayIndx), miscCoord_New(q4ArrayIndx)
!           endif
!           do m = 0, 12
!             if(abs(q4real(m) - q4r_t(m)) .gt. 1d-7) then
!               write(2,*) "DISCRPANCY:",m, q4real(m) , q4r_t(m)
!             endif
!           enddo
!           if(q4Neigh .ne. q4Neigh_t) then
!             write(2,*) "DISCRPANCY:", q4Neigh, q4Neigh_t
!           endif
!         endif
         initialized = .true.
      end subroutine
     !--------------------------------------------------------------------------------
       subroutine Calcq4_Disp(disp)
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
         real(dp) :: q4r_new, q4i_new
         real(dp) :: rx,ry,rz,r
         real(dp) :: theta,phi,junk,sTerm,cTerm
         real(dp) :: q4Par
         real(dp) :: rPart,iPart

         dNeigh = 0
         dq4real = 0E0_dp
         dq4img = 0E0_dp
         if(NTotal .eq. 1) then
 	       miscCoord_New(q4ArrayIndx) = 1E0_dp
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
           miscCoord_New(q4ArrayIndx) = miscCoord(q4ArrayIndx)
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
             if(r .le. q4DistSq) then	
               r = sqrt(r)
               dNeigh = dNeigh + 1
               phi = atan2( ry, rx )
               theta = acos( rz / r )  
               do m = 0, 8
                 call Harmonics(theta, phi, m-4, rPart, iPart)	
                 dq4real(m) = dq4real(m) + rPart
                 dq4img(m) = dq4img(m) + iPart	
               enddo               
             endif
             

             if(rPair(gloIndx2, gloIndx1)%p%r_sq .le. q4DistSq) then	
                dNeigh = dNeigh - 1
!                rx = rPair(gloIndx2, gloIndx1)%p%rx
!                ry = rPair(gloIndx2, gloIndx1)%p%ry
!                rz = rPair(gloIndx2, gloIndx1)%p%rz
                phi = atan2(rPair(gloIndx2, gloIndx1)%p%ry, rPair(gloIndx2, gloIndx1)%p%rx)
                theta = acos(rPair(gloIndx2, gloIndx1)%p%rz / rPair(gloIndx2, gloIndx1)%p%r )   
                do m = 0, 8
                  call Harmonics(theta, phi, m-4, rPart, iPart)	
                  dq4real(m) = dq4real(m) - rPart
                  dq4img(m) = dq4img(m) - iPart	
                enddo
             endif   
           enddo
         enddo
          
         q4par = 0E0_dp
         do m = 0, 8
           q4r_new = q4real(m) + dq4real(m)
           q4i_new = q4img(m) + dq4img(m)
           q4Par = q4par + q4r_new*q4r_new + q4i_new*q4i_new   	   
         enddo
         q4par = q4par * q4Constant
         newNeigh = q4Neigh + dNeigh
         q4par = sqrt(q4par)/newNeigh
 	     miscCoord_New(q4ArrayIndx) = q4par
         newData = .true.
      end subroutine
     !--------------------------------------------------------------------------------
       subroutine Calcq4_SwapIn
         use MiscelaniousVars
         use Coords
         use SimParameters, only: nMolTypes, NPART, NTotal
         implicit none 
         integer :: m
         integer :: iType, iMol, jType, jMol
         integer :: gloIndx1, gloIndx2
         integer :: jIndx, iIndx, sizeDisp
         integer :: dispAtom1
         real(dp) :: q4r_new, q4i_new
         real(dp) :: rx,ry,rz,r
         real(dp) :: theta,phi,junk,sTerm,cTerm
         real(dp) :: q4Par
         real(dp) :: rPart,iPart

 
         iType = newMol%molType
         iMol = NPART(iType) + 1
         dNeigh = 0
!         gloIndx1 =  molArray(iType)%mol(iMol)%globalIndx(1)
         dq4real = 0E0_dp
         dq4img = 0E0_dp
         do jType = 1, nMolTypes
           do jMol = 1, NPART(jType)
             rx = newMol%x(1) -  molArray(jType)%mol(jMol)%x(1)
             ry = newMol%y(1) -  molArray(jType)%mol(jMol)%y(1)
             rz = newMol%z(1) -  molArray(jType)%mol(jMol)%z(1)
             r = rx*rx + ry*ry + rz*rz
             if(r .le. q4DistSq) then	
               r = sqrt(r)
               dNeigh = dNeigh + 1
               phi = atan2(ry,rx)
               theta = acos(rz/r)  
               do m = 0, 8
                 call Harmonics(theta, phi, m-4, rPart, iPart)	
                 dq4real(m) = dq4real(m) + rPart
                 dq4img(m) = dq4img(m) + iPart	
               enddo               
             endif
           enddo
         enddo
          
         q4par = 0E0_dp
         do m = 0, 8
           q4r_new = q4real(m) + dq4real(m)
           q4i_new = q4img(m) + dq4img(m)
           q4Par = q4par + q4r_new*q4r_new + q4i_new*q4i_new   	   
         enddo
         q4par = q4par * q4Constant
         q4par = sqrt(q4par)/real((q4Neigh + dNeigh), dp)
 	     miscCoord_New(q4ArrayIndx) = q4par
         newData = .true.
      end subroutine
     !--------------------------------------------------------------------------------
       subroutine Calcq4_SwapOut(nType, nMol)
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
         real(dp) :: q4r_new, q4i_new
         real(dp) :: rx,ry,rz,r
         real(dp) :: theta,phi,junk,sTerm,cTerm
         real(dp) :: q4Par
         real(dp) :: rPart,iPart

         if(NTotal .eq. 2) then
!           dNeigh = -q4Neigh
!           dq4real = -q4real
!           dq4real = -q4img
 	       miscCoord_New(q4ArrayIndx) = 1E0
           return	   
         endif

 
         dNeigh = 0
         dq4real = 0E0
         dq4img = 0E0
         gloIndx1 =  molArray(nType)%mol(nMol)%globalIndx(1)
!         iIndx =  molArray(nType)%mol(nMol)%indx
         do jType = 1, nMolTypes
           do jMol = 1, NPART(jType)
!             jIndx =  molArray(jType)%mol(jMol)%indx
             gloIndx2 =  molArray(jType)%mol(jMol)%globalIndx(1)
             if(gloIndx1 .eq. gloIndx2) then
               cycle
             endif
             if(rPair(gloIndx2, gloIndx1)%p%r_sq .le. q4DistSq) then	
               r = rPair(gloIndx2, gloIndx1)%p%r
               dNeigh = dNeigh - 1
               rx = rPair(gloIndx2, gloIndx1)%p%rx
               ry = rPair(gloIndx2, gloIndx1)%p%ry
               rz = rPair(gloIndx2, gloIndx1)%p%rz
               phi = atan2(ry,rx)
               theta = acos(rz/r)  
               do m = 0, 8
                 call Harmonics(theta, phi, m-4, rPart, iPart)	
                 dq4real(m) = dq4real(m) - rPart
                 dq4img(m) = dq4img(m) - iPart	
               enddo               
             endif
           enddo
         enddo
          
         q4par = 0E0
         do m = 0, 8
           q4r_new = q4real(m) + dq4real(m)
           q4i_new = q4img(m) + dq4img(m)
           q4Par = q4par + q4r_new*q4r_new + q4i_new*q4i_new   	   
         enddo
         q4par = q4par * q4Constant
         q4par = sqrt(q4par)/real((q4Neigh + dNeigh), dp)
 	     miscCoord_New(q4ArrayIndx) = q4par
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
      real(dp), parameter :: const0 = (3d0/16d0)*sqrt(1d0/pi)
      real(dp), parameter :: const1 = (3d0/8d0)*sqrt(5d0/pi)
      real(dp), parameter :: const2 = (3d0/8d0)*sqrt(5d0/(2d0*pi))
      real(dp), parameter :: const3 = (3d0/8d0)*sqrt(35d0/pi)
      real(dp), parameter :: const4 = (3d0/16d0)*sqrt(35d0/(2d0*pi))

      RealVal = 0E0_dp
      ImgVal = 0E0_dp
      cterm = cos(Theta)
      sterm = sin(Theta)
	  
      select case (abs(m))
      case(0)
        ctermsq = cterm * cterm
        RealVal = 35d0*ctermsq*ctermsq - 30d0*ctermsq + 3d0
        RealVal = RealVal * const0
      case(1) 
        ctermsq = cterm * cterm
        RealVal = sterm * (7d0*ctermsq*cterm - 3d0 * cterm)
        if(m .gt. 0) then
	       RealVal = -RealVal
        endif
        RealVal = RealVal * const1
      case(2)
        ctermsq = cterm * cterm
        stermsq = sterm * sterm
        RealVal = stermsq * (7d0*ctermsq - 1d0)
        RealVal = RealVal * const2
      case(3)	   
        RealVal = sTerm * sTerm * sTerm
        RealVal = RealVal * cterm 
        if(m .gt. 0) then
	      RealVal = -RealVal
        endif
        RealVal = RealVal * const3
      case(4)
        sTermsq = sterm * sterm
        RealVal = sTermsq * sTermsq
        RealVal = RealVal * const4
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



