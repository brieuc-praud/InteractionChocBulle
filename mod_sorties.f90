module mod_sorties

  implicit none
  
contains
  
  subroutine sorties_fichier(T, x, y, nb_fich)

    use mod_params
    
    !--- entrees
    integer, intent(in) :: nb_fich
    real(pr), dimension(0:), intent(in) :: x, y
    real(pr), dimension(:,:), intent(in) :: T
    
    !--- locales
    integer :: imax, jmax, i, j
    character(len=30) :: fich

    
    imax = size(x) - 1
    jmax = size(y) - 1

    write(fich,*) nb_fich
    open(unit=20, file='T_'//trim(adjustl(fich))//'.vtk')

    write(20,'(a)') '# vtk DataFile Version 2.0'
    write(20,'(a)') 'Temperature'
    write(20,'(a)') 'ASCII'
    write(20,'(a)') 'DATASET RECTILINEAR_GRID'
    write(20,'(1a10,3i4)') 'DIMENSIONS ',imax+1,jmax+1,1
    write(20,fmt='(1a13,1i10,1a10)') 'X_COORDINATES',imax+1,' float'
    do i=0,imax
       write(20,*) x(i)
    end do
    write(20,fmt='(1a13,1i10,1a10)') 'Y_COORDINATES',jmax+1,' float'
    do j=0,jmax
       write(20,*) y(j)
    end do
    write(20,fmt='(1a13,1i10,1a10)') 'Z_COORDINATES',1,' float'
    write(20,*) 0.0_pr

    write(20,fmt='(1a9,1i10)') 'CELL_DATA',imax*jmax
    write(20,'(a)') 'SCALARS T double'
    write(20,'(a)') 'LOOKUP_TABLE default'
    do j=1,jmax
       do i=1,imax
          write(20,*) T(i,j)
       end do
    end do

    close(20)

  end subroutine sorties_fichier

end module mod_sorties
