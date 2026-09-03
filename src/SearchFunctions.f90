!=====================================================
      integer function getBond(atm1, atm2, iType)
      use ForceField
      implicit none
      integer, intent(in) :: atm1, atm2, iType
      integer :: iBond
      
      do iBond = 1, nBonds(iType)
        if(any(bondArray(iType)%bondMembr .eq. atm1) then
          if(any(bondArray(iType)%bondMembr .eq. atm2) then        
            getBond = iBond
            return
          endif
        endif
      enddo
      write(6,*) "Unable to find Bond! [Error in getBond function]"
      stop
            
      end function
!=====================================================
      integer function getBending(atm1, atm2, atm3, iType)
      use ForceField
      implicit none
      integer, intent(in) :: atm1, atm2,atm3, iType
      integer :: iBend
      
      do iBend = 1, nAngles(iType)
        if(any(bendArray(iType)%bendMembr .eq. atm1) then
          if(any(bendArray(iType)%bendMembr .eq. atm2) then        
           if(any(bendArray(iType)%bendMembr .eq. atm3) then
             getTorsion = iBend
             return
           endif            
          endif
        endif
      enddo
      write(6,*) "Unable to find Torsional Angle! [Error in getTorsion function]"
      stop
            
      end function      
!=====================================================
      integer function getTorsion(atm1, atm2, atm3, atm4, iType)
      use ForceField
      implicit none
      integer, intent(in) :: atm1, atm2,atm3,atm4, iType
      integer :: iTors
      
      do iTors = 1, nTorsional(iType)
        if(any(torsArray(iType)%torsMembr .eq. atm1) then
          if(any(torsArray(iType)%torsMembr .eq. atm2) then        
           if(any(torsArray(iType)%torsMembr .eq. atm3) then
             if(any(torsArray(iType)%torsMembr .eq. atm4) then           
               getTorsion = iTors
               return
             endif
           endif            
          endif
        endif
      enddo
      write(6,*) "Unable to find Torsional Angle! [Error in getTorsion function]"
      stop
            
      end function
      
!=====================================================