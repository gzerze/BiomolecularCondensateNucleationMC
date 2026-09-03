      subroutine Output_VMD_Final
      use Coords
      use Forcefield
      use SimParameters
      implicit none
      integer :: iType,iMol,iAtom
      integer :: atmType
      integer :: nTotalAtoms
      
      nTotalAtoms = 0
      
      do iType=1,nMolTypes
        nTotalAtoms = nTotalAtoms + NPART(iType)*nAtoms(iType)
      enddo
      
      open(unit=15, file="VMD_Output.xyz")      
      write(15,*) nTotalAtoms
      write(15,*) NPART
      do iType=1, nMolTypes
        do iMol = 1,NPART(iType)
          do iAtom = 1,nAtoms(iType)
            atmType = atomArray(iType,iAtom)
            write(15,*) atomData(atmType)%Symb &
                        ,MolArray(iType)%Mol(iMol)%x(iAtom) &
                        ,MolArray(iType)%Mol(iMol)%y(iAtom) &
                        ,MolArray(iType)%Mol(iMol)%z(iAtom)     
          enddo
        enddo
      enddo
      close(15)
      
      end subroutine
