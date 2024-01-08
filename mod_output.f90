Module mod_output

    Implicit None

Contains
    Subroutine output(U, gammagp, x, y, nb_filename, filename_prefix)
        Use mod_parameters

        !--- Intent In
        Integer, Intent(In) :: nb_filename
        Character(len=*) :: filename_prefix
        Real(PR), Dimension(0:), Intent(In) :: x, y
        Real(PR), Dimension(:,:,:), Intent(In) :: U
        Real(PR), Intent(In) :: gammagp

        !--- Locals
        Integer :: imax, jmax, i, j
        Character(len=30) :: filename
        Real(PR), Dimension(Size(x)-1,Size(y)-1) :: density, velocity_u, velocity_v, energy, q, pressure

        imax = Size(x) - 1
        jmax = Size(y) - 1

        density = U(1,:,:)
        velocity_u = U(2,:,:) / density
        velocity_v = U(3,:,:) / density
        energy = U(4,:,:)

        q = .5_PR * ( velocity_u**2 + velocity_v**2 )
        pressure = (gammagp - 1._PR)*(energy - density*q)

        Write(filename,*) nb_filename
        Open(Unit=111, File="output/"//TRIM(ADJUSTL(filename_prefix))//"_"//TRIM(ADJUSTL(filename))//".vtk")

        Write(111,"(a)") "# vtk DataFile Version 2.0"
        Write(111,"(a)") "Euler Equations"
        Write(111,"(a)") "ASCII"
        Write(111,"(a)") "DATASET RECTILINEAR_GRID"
        Write(111,"(1a10,3i4)") "DIMENSIONS ",imax+1,jmax+1,1
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
