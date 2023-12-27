Program euler
    Use mod_functions
    Use mod_schemes
    Use mod_output
    Use mod_test

    Implicit None

    Character(14), Parameter :: parameters = "parameters.dat"

    ! Parameters
    Real(PR) xmin, xmax, ymin, ymax, time_max, cfl, gammagp
    Integer :: imax, jmax, output_modulo, case_number
    Character(len=10) :: numflux_name
    ! Input related variables
    Character(len=100) :: buffer, label
    Integer :: pos
    Integer :: ios = 0
    Integer :: line_number = 0
    ! Arrays
    Real(PR), Dimension(:,:,:), Allocatable :: Uvect, Uvect_e, fluxF, fluxG
    Real(PR), Dimension(:), Allocatable :: x, y, xm, ym
    ! Loop indices
    Integer :: i, j, nb_iterations
    ! Other
    Real(PR) :: deltax, deltay, deltat, time


    ! Read parameters
    Open(111, File=parameters)
    Do While (ios == 0)
        Read(111, '(A)', IOstat=ios) buffer
        If (ios == 0) Then
            line_number = line_number + 1

            pos = SCAN(buffer, ' ')
            label = buffer(1:pos)
            buffer = buffer(pos+1:)

            Select Case (label)
            Case ('xmin')
                Read(buffer, *, iostat=ios) xmin
            Case ('xmax')
                Read(buffer, *, iostat=ios) xmax
            Case ('ymin')
                Read(buffer, *, iostat=ios) ymin
            Case ('ymax')
                Read(buffer, *, iostat=ios) ymax
            Case ('Nx')
                Read(buffer, *, iostat=ios) imax
            Case ('Ny')
                Read(buffer, *, iostat=ios) jmax
            Case ('tmax')
                Read(buffer, *, iostat=ios) time_max
            Case ('CFL')
                Read(buffer, *, iostat=ios) cfl
            Case ('gamma')
                Read(buffer, *, iostat=ios) gammagp
            Case ('output_modulo')
                Read(buffer, *, iostat=ios) output_modulo
            Case ('case')
                Read(buffer, *, iostat=ios) case_number
            Case ('flux')
                Read(buffer, *, iostat=ios) numflux_name
            Case ('')
                ! Do nothing if it is an empty line
            Case Default
                Write(STDOUT,*) "Invalid label", label, " at line", line_number, "(skipping)"
            End Select
        End If
    End Do
    Close(111)
    ! Allocate
    Allocate(x(0:imax), y(0:jmax), xm(imax), ym(jmax))
    Allocate(Uvect(4,imax,jmax), Uvect_e(4,imax,jmax), fluxF(4,0:imax, 0:jmax), fluxG(4,0:imax, 0:jmax))
    ! Compute the grid
    deltax = (xmax - xmin) / imax
    deltay = (ymax - ymin) / jmax
    x = (/ Real(PR) :: (xmin + i*deltax, i=0, imax) /)
    y = (/ Real(PR) :: (ymin + j*deltay, j=0, jmax) /)
    xm = (/ Real(PR) :: (xmin + .5_PR*deltax + i*deltax, i=0, imax-1) /)
    ym = (/ Real(PR) :: (ymin + .5_PR*deltay + j*deltay, j=0, jmax-1) /)
    ! Initialise U
    Do i=1, imax
        Do j=1, jmax
            Uvect(:,i,j) = Uinit(case_number, xm(i), ym(j), gammagp)
        End Do
    End Do

    !Call PerformTestsAndExit(gammagp)

    ! Output initial state
    Call output(Uvect, gammagp, x, y, 0, 'sol')
    Call output(Uvect, gammagp, x, y, 0, 'exact')

    ! Time loop
    nb_iterations = 0
    time = 0._PR
    Do While (time < time_max)
        Call compute_CFL(Uvect, deltax, deltay, deltat, cfl)
        time = MIN( time + deltat, time_max )

        Do j=1, jmax
            Do i=1, imax-1
                Select Case (TRIM(ADJUSTL(numflux_name)))
                Case ('Rusanov')
                    fluxF(:,i,j) = Rusanov('x', Uvect(:,i,j), Uvect(:,i+1,j), gammagp)
                Case ('HLLC')
                    fluxF(:,i,j) = HLLC('x', Uvect(:,i,j), Uvect(:,i+1,j), gammagp)
                Case Default ! Case ('HLL')
                    fluxF(:,i,j) = HLL('x', Uvect(:,i,j), Uvect(:,i+1,j), gammagp)
                End Select
            End Do
            ! Boundary
            Select Case (case_number)
            Case (2) ! Periodic
                ! Periodic
                Select Case (TRIM(ADJUSTL(numflux_name)))
                Case ('Rusanov')
                    fluxF(:,0,j) = Rusanov('x', Uvect(:,imax,j), Uvect(:,1,j), gammagp)
                Case ('HLLC')
                    fluxF(:,0,j) = HLLC('x', Uvect(:,imax,j), Uvect(:,1,j), gammagp)
                Case Default ! Case ('HLL')
                    fluxF(:,0,j) = HLL('x', Uvect(:,imax,j), Uvect(:,1,j), gammagp)
                End Select
                fluxF(:,imax,j) = fluxF(:,0,j)
            Case Default ! Absorbing
                fluxF(:,0,j) = fluxFuncF( Uvect(:,1,j), gammagp )
                fluxF(:,imax,j) = fluxFuncF( Uvect(:,imax,j), gammagp )
            End Select
        End Do
        Do i=1, imax
            Do j=1, jmax-1
                Select Case (TRIM(ADJUSTL(numflux_name)))
                Case ('Rusanov')
                    fluxG(:,i,j) = Rusanov('y', Uvect(:,i,j), Uvect(:,i,j+1), gammagp)
                Case ('HLLC')
                    fluxG(:,i,j) = HLLC('y', Uvect(:,i,j), Uvect(:,i,j+1), gammagp)
                Case Default ! Case ('HLL')
                    fluxG(:,i,j) = HLL('y', Uvect(:,i,j), Uvect(:,i,j+1), gammagp)
                End Select
            End Do
            ! Boundary
            Select Case (case_number)
            Case (2) ! Periodic
                Select Case (TRIM(ADJUSTL(numflux_name)))
                Case ('Rusanov')
                    fluxG(:,i,0) = Rusanov('y', Uvect(:,i,jmax), Uvect(:,i,1), gammagp)
                Case ('HLLC')
                    fluxG(:,i,0) = HLLC('y', Uvect(:,i,jmax), Uvect(:,i,1), gammagp)
                Case Default ! Case ('HLL')
                    fluxG(:,i,0) = HLL('y', Uvect(:,i,jmax), Uvect(:,i,1), gammagp)
                End Select
                fluxG(:,i,jmax) = fluxG(:,i,0)
            Case Default ! Absorbing
                fluxG(:,i,0) = fluxFuncG( Uvect(:,i,0), gammagp )
                fluxG(:,i,jmax) = fluxFuncG( Uvect(:,i,jmax), gammagp )
            End Select
        End Do

        Do i=1, imax
            Do j=1, jmax
                Uvect(:,i,j) = Uvect(:,i,j) &
                    & - deltat/deltax * (fluxF(:,i,j) - fluxF(:,i-1,j)) &
                    & - deltat/deltay * (fluxG(:,i,j) - fluxG(:,i,j-1))
            End Do
        End Do

        If ( output_modulo > 0 .AND. Modulo(nb_iterations, output_modulo) == 0 ) Then
            Write(STDOUT, *) time, time_max
            Do i=1, imax
                Do j=1, jmax
                    Uvect_e(:,i,j) = Uexact(case_number, xm(i), ym(j), time, gammagp)
                End Do
            End Do
            Call output(Uvect, gammagp, x, y, nb_iterations / output_modulo + 1, 'sol')
            Call output(Uvect_e, gammagp, x, y, nb_iterations / output_modulo + 1, 'exact')
        End If
        nb_iterations = nb_iterations + 1
    End Do

    Write(STDOUT, *) "Error:", error('Linfty', case_number, Uvect, time_max, gammagp)

    Deallocate(x, y, xm, ym)
    Deallocate(Uvect, Uvect_e, fluxF, fluxG)

Contains
    Subroutine compute_CFL(U, dx, dy, dt, cfl)
        ! --- InOut ---
        Real(PR), Dimension(:,:,:), Intent(In) :: U
        Real(PR), Intent(In) :: dx, dy,  cfl
        Real(PR), Intent(Out) :: dt
        ! --- Locals ---
        Real(PR) :: rho, velocity_u, velocity_v, e, q, p, a, bx, by, bx_max, by_max, l1, l3

        bx_max = 0._PR
        by_max = 0._PR
        Do i=1, imax
            Do j=1 , jmax
                rho = U(1,i,j)
                velocity_u = U(2,i,j) / rho
                velocity_v = U(3,i,j) / rho
                e = U(4,i,j)

                q = .5_PR * ( velocity_u**2 + velocity_v**2 )
                p = (gammagp - 1._PR)*(e - rho*q)
                a = SQRT(gammagp*p/rho)

                l1 = ABS(velocity_u - a)
                l3 = ABS(velocity_u + a)
                bx = MAX(l1, l3)
                bx_max = MAX( bx, bx_max)
                l1 = ABS(velocity_v - a)
                l3 = ABS(velocity_v + a)
                by = MAX(l1, l3)
                by_max = MAX( by, by_max)
            End Do
        End Do

        dt = cfl * 0.5_PR * MIN(dx/bx_max, dy/by_max)
    End Subroutine compute_CFL

    Function error(norm, case_number, U, time, gammagp)
        ! --- InOut
        Real(PR), Dimension(4,imax,jmax), Intent(In) :: U
        Real(PR), Intent(In) :: gammagp, time
        Integer, Intent(In) :: case_number
        Character(len=*), Intent(In) :: norm
        Real(PR), Dimension(4) :: error
        ! --- Locals
        Real(PR), Dimension(4) :: exact_value

        error = 0._PR
        Select Case (norm)
        Case ('L1') ! L1 norm
            Do i=1, imax
                Do j=1 , jmax
                    exact_value = Uexact(case_number, xm(i), ym(j), time, gammagp)
                    error = error + ABS( U(:,i,j) - exact_value )
                End Do
            End Do
            error = error / ( imax * jmax )
        Case ('L2') ! L2 norm
            Do i=1, imax
                Do j=1 , jmax
                    exact_value = Uexact(case_number, xm(i), ym(j), time, gammagp)
                    error = error + ( U(:,i,j) - exact_value )**2
                End Do
            End Do
            error = SQRT( error / ( imax * jmax ) )
        Case Default ! L_infinity norm
            Do i=1, imax
                Do j=1 , jmax
                    exact_value = Uexact(case_number, xm(i), ym(j), time, gammagp)
                    error = MAX( error, ABS( U(:,i,j) - exact_value ) )
                End Do
            End Do
        End Select
    End Function error

End Program euler
