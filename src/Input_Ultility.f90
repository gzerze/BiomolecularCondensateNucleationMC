      subroutine GCD_Calc(a,b,divi_last)
      implicit none     
      integer :: a,b,divi
      integer :: divi_last,temp

      if(a .gt. b) then   
        divi_last = a
        divi = b
      else
        divi_last = b
        divi = a
      endif   
        
        
      do while (divi .ne. 0) 
        temp=divi
        divi=mod(divi_last,divi)
        divi_last=temp        
      enddo
        
      end