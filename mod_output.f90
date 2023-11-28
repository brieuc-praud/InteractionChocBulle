module mod_output

   implicit none

contains

   subroutine output(U, x, y, nb_fich)

      use mod_parameters

      !--- Intent In
      Integer, Intent(In) :: nb_fich
      Real(pr), Dimension(0:), Intent(In) :: x, y
      Real(pr), Dimension(:,:,:), Intent(In) :: U

      !--- Locals
      Integer :: imax, jmax, i, j
      Character(len=30) :: fich

      imax = size(x) - 1
      jmax = size(y) - 1

      Write(fich,*) nb_fich
      Open(unit=20, file='T_'//trim(adjustl(fich))//'.vtk')

      Write(20,'(a)') '# vtk DataFile Version 2.0'
      Write(20,'(a)') 'Euler Equations'
      Write(20,'(a)') 'ASCII'
      Write(20,'(a)') 'DATASET RECTILINEAR_GRID'
      Write(20,'(1a10,3i4)') 'DIMENSIONS ',imax+1,jmax+1,1
      Write(20,fmt='(1a13,1i10,1a10)') 'X_COORDINATES',imax+1,' float'
      Do i=0,imax
         Write(20,*) x(i)
      End Do
      Write(20,fmt='(1a13,1i10,1a10)') 'Y_COORDINATES',jmax+1,' float'
      Do j=0,jmax
         Write(20,*) y(j)
      End Do
      Write(20,fmt='(1a13,1i10,1a10)') 'Z_COORDINATES',1,' float'
      Write(20,*) 0.0_PR

      Write(20,fmt='(1a9,1i10)') 'CELL_DATA', imax*jmax
      Write(20,'(a)') 'SCALARS density double'
      Write(20,'(a)') 'LOOKUP_TABLE default'
      Do j=1,jmax
         Do i=1,imax
            Write(20,*) U(1,i,j)
         End Do
      End Do
      Write(20,fmt='(1a9,1i10)') 'CELL_DATA', imax*jmax
      Write(20,'(a)') 'VECTORS density double'
      Write(20,'(a)') 'LOOKUP_TABLE default'
      Do j=1,jmax
         Do i=1,imax
            Write(20,*) U(2,i,j)/U(1,i,j), U(3,i,j)/U(1,i,j), 0.
         End do
      End do
      Write(20,fmt='(1a9,1i10)') 'CELL_DATA', imax*jmax
      Write(20,'(a)') 'SCALARS energy double'
      Write(20,'(a)') 'LOOKUP_TABLE default'
      Do j=1,jmax
         Do i=1,imax
            Write(20,*) U(4,i,j)
         End Do
      End Do

      Close(20)

   end subroutine output

end module mod_output
