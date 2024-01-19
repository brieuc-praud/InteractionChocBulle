Module mod_output

    Use mod_parameters

    Implicit None

Contains
    Subroutine output(U, nb_filename, filename_prefix)
        !--- InOut
        Real(PR), Dimension(4,-1:imax+2,-1:jmax+2) :: U
        Integer, Intent(In) :: nb_filename
        Character(len=*), Intent(In) :: filename_prefix
        !--- Locals
        Integer :: i, j
        Character(len=30) :: filename
        Real(PR), Dimension(imax,jmax) :: density, velocity_u, velocity_v, energy, q, pressure

        density = U(1,1:imax,1:jmax)
        velocity_u = U(2,1:imax,1:jmax) / density
        velocity_v = U(3,1:imax,1:jmax) / density
        energy = U(4,1:imax,1:jmax)

        q = .5_PR * ( velocity_u**2 + velocity_v**2 )
        pressure = (gammagp - 1._PR)*(energy - density*q)

        Write(filename,*) nb_filename
        Open(Unit=111, File="output/"//TRIM(ADJUSTL(filename_prefix))//"_"//TRIM(ADJUSTL(filename))//".vtk")

        Write(111,"(a)") "# vtk DataFile Version 2.0"
        Write(111,"(a)") "Euler Equations"
        Write(111,"(a)") "ASCII"
        Write(111,"(a)") "DATASET RECTILINEAR_GRID"
        Write(111,"(1a10,3i5)") "DIMENSIONS ",imax+1,jmax+1,1
        Write(111,Fmt="(1a13,1i10,1a10)") "X_COORDINATES",imax+1," float"
        Do i=0,imax
            Write(111,*) x(i)
        End Do
        Write(111,Fmt="(1a13,1i10,1a10)") "Y_COORDINATES",jmax+1," float"
        Do j=0,jmax
            Write(111,*) y(j)
        End Do
        Write(111,Fmt="(1a13,1i10,1a10)") "Z_COORDINATES",1," float"
        Write(111,*) 0.0_PR

        Write(111,Fmt="(1a9,1i10)") "CELL_DATA", imax*jmax
        Write(111,"(a)") "SCALARS density double"
        Write(111,"(a)") "LOOKUP_TABLE default"
        Do j=1,jmax
            Do i=1,imax
                Write(111,*) density(i,j)
            End Do
        End Do
        Write(111,"(a)") "SCALARS pressure double"
        Write(111,"(a)") "LOOKUP_TABLE default"
        Do j=1,jmax
            Do i=1,imax
                Write(111,*) pressure(i,j)
            End Do
        End Do
        Write(111,"(a)") "VECTORS velocity double"
        Do j=1,jmax
            Do i=1,imax
                Write(111,*) velocity_u(i,j), velocity_v(i,j), 0._PR
            End do
        End do

        Close(111)

    End Subroutine output

End Module mod_output
