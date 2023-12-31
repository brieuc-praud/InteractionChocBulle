Program euler
    Use mod_fluxes
    Use mod_output
    Use mod_cases
    Use mod_test

    Implicit None

    Character(14), Parameter :: parameters = "parameters.dat"

    ! Parameters
    Real(PR) xmin, xmax, ymin, ymax, time_max, cfl, gammagp
    Integer :: imax, jmax, output_modulo
    Character(len=10) :: numflux_name, norm_str, time_scheme_name
    Character(len=20) :: case_name
    ! Input related variables
    Character(len=100) :: buffer
    Character(len=20) :: label
    Integer :: pos
    Integer :: ios = 0
    Integer :: line_number = 0
    ! Arrays
    Real(PR), Dimension(:,:,:), Allocatable :: Uvect, Uvect_e, K1vect, K2vect, K3vect
    Real(PR), Dimension(:,:,:), Allocatable :: fluxF, fluxG, fluxK1F, fluxK1G, fluxK2F, fluxK2G, fluxK3F, fluxK3G
    Real(PR), Dimension(:), Allocatable :: x, y, xm, ym
    ! Loop indices
    Integer :: i, j, nb_iterations
    ! Other
    Real(PR) :: deltax, deltay, deltat, time
    Type(space_scheme) :: space_scheme_specs
    Logical :: exact_solution_available


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
            Case ('flux')
                Read(buffer, *, iostat=ios) numflux_name
            Case ('time_scheme')
                Read(buffer, *, iostat=ios) time_scheme_name
            Case ('space_scheme')
                Read(buffer, *, iostat=ios) space_scheme_specs%order, space_scheme_specs%slope_str, &
                    & space_scheme_specs%limiter_str, space_scheme_specs%generalised_minmod_parameter
            Case ('case')
                Read(buffer, *, iostat=ios) case_name
                exact_solution_available = exactSolutionAvailable(case_name)
            Case ('error_norm')
                Read(buffer, *, iostat=ios) norm_str
            Case ('','#')
                ! Do nothing if it is an empty line or a comment
            Case Default
                If ( label(1:1) /= '#' ) Then ! Special case where there is no space after '#'
                    Write(STDOUT,*) "Invalid label ", label, " at line ", line_number, " (skipping)"
                End If
            End Select
        End If
    End Do
    Close(111)
    ! Allocate
    Allocate(x(0:imax), y(0:jmax), xm(imax), ym(jmax))
    Allocate(Uvect(4,imax,jmax), Uvect_e(4,imax,jmax), fluxF(4,0:imax, 0:jmax), fluxG(4,0:imax, 0:jmax))
    Select Case (TRIM(ADJUSTL(time_scheme_name)))
    Case ('SSP-RK2') ! Strong-Stability preserving Runge-Kutta
        Allocate(K1vect(4,imax,jmax), fluxK1F(4,0:imax, 0:jmax), fluxK1G(4,0:imax, 0:jmax))
        Allocate(K2vect(4,imax,jmax), fluxK2F(4,0:imax, 0:jmax), fluxK2G(4,0:imax, 0:jmax))
    Case ('SSP-RK3') ! Strong-Stability preserving Runge-Kutta
        Allocate(K1vect(4,imax,jmax), fluxK1F(4,0:imax, 0:jmax), fluxK1G(4,0:imax, 0:jmax))
        Allocate(K2vect(4,imax,jmax), fluxK2F(4,0:imax, 0:jmax), fluxK2G(4,0:imax, 0:jmax))
        Allocate(K3vect(4,imax,jmax), fluxK3F(4,0:imax, 0:jmax), fluxK3G(4,0:imax, 0:jmax))
    End Select
    ! Compute the grid
    Call getGridDimensions(case_name, xmin, xmax, ymin, ymax)
    deltax = (xmax - xmin) / imax
    deltay = (ymax - ymin) / jmax
    x = (/ Real(PR) :: (xmin + i*deltax, i=0, imax) /)
    y = (/ Real(PR) :: (ymin + j*deltay, j=0, jmax) /)
    xm = (/ Real(PR) :: (xmin + .5_PR*deltax + i*deltax, i=0, imax-1) /)
    ym = (/ Real(PR) :: (ymin + .5_PR*deltay + j*deltay, j=0, jmax-1) /)
    ! Initialise U
    Do i=1, imax
        Do j=1, jmax
            Uvect(:,i,j) = Uinit(case_name, xm(i), ym(j), gammagp)
        End Do
    End Do

    !Call PerformTestsAndExit(gammagp)

    ! Output initial state
    If ( output_modulo > 0 ) Then
        Call output(Uvect, gammagp, x, y, 0, 'sol')
        Call output(Uvect, gammagp, x, y, 0, 'exact')
    End If

    ! Time loop
    nb_iterations = 0
    time = 0._PR
    Do While (time < time_max)
        Call compute_CFL(Uvect, time, deltax, deltay, deltat, cfl)
        time = MIN( time + deltat, time_max )

        Call step(time, time_scheme_name)

        If ( output_modulo > 0 .AND. Modulo(nb_iterations, output_modulo) == 0 ) Then
            Write(STDOUT, *) time, time_max
            If ( exact_solution_available ) Then
                Do i=1, imax
                    Do j=1, jmax
                        Uvect_e(:,i,j) = Uexact(case_name, xm(i), ym(j), time, gammagp)
                    End Do
                End Do
                Call output(Uvect_e, gammagp, x, y, nb_iterations / output_modulo + 1, 'exact')
            End If
            Call output(Uvect, gammagp, x, y, nb_iterations / output_modulo + 1, 'sol')
        End If
        nb_iterations = nb_iterations + 1
    End Do

    If (exact_solution_available) Then
        Write(STDOUT, *) "Error:", error(norm_str, case_name, Uvect, time_max, gammagp, .False.)
    End If

    Deallocate(x, y, xm, ym)
    Deallocate(Uvect, Uvect_e, fluxF, fluxG)
    Select Case (TRIM(ADJUSTL(time_scheme_name)))
    Case ('SSP-RK2') ! Strong-Stability preserving Runge-Kutta
        Deallocate(K1vect, fluxK1F, fluxK1G)
        Deallocate(K2vect, fluxK2F, fluxK2G)
    Case ('SSP-RK3') ! Strong-Stability preserving Runge-Kutta
        Deallocate(K1vect, fluxK1F, fluxK1G)
        Deallocate(K2vect, fluxK2F, fluxK2G)
        Deallocate(K3vect, fluxK3F, fluxK3G)
    End Select
Contains

    Subroutine step(time, time_scheme_name)
        ! --- InOut
        Character(len=*), Intent(In) :: time_scheme_name
        Real(PR), Intent(In) :: time

        Select Case (TRIM(ADJUSTL(time_scheme_name)))
        Case ('EE','ExplicitEuler') ! Explicit Euler
            Call compute_Fluxes('x', time, Uvect, fluxF)
            Call compute_Fluxes('y', time, Uvect, fluxG)

            Call ExplicitEuler(Uvect, Uvect, fluxF, fluxG, &
                & deltax, deltay, deltat, imax, jmax)
        Case ('SSP-RK2') ! Strong-Stability preserving Runge-Kutta
            ! First stage
            Call compute_Fluxes('x', time, Uvect, fluxF)
            Call compute_Fluxes('y', time, Uvect, fluxG)

            Call ExplicitEuler(K1vect, Uvect, fluxF, fluxG, &
                & deltax, deltay, deltat, imax, jmax)

            ! Second stage
            Call compute_Fluxes('x', time, K1vect, fluxK1F)
            Call compute_Fluxes('y', time, K1vect, fluxK1G)

            Call ExplicitEuler(K2vect, K1vect, fluxK1F, fluxK1G, &
                & deltax, deltay, deltat, imax, jmax)

            Uvect = .5_PR * Uvect + .5_PR * K2vect
        Case ('SSP-RK3') ! Strong-Stability preserving Runge-Kutta
            ! First stage
            Call compute_Fluxes('x', time, Uvect, fluxF)
            Call compute_Fluxes('y', time, Uvect, fluxG)

            Call ExplicitEuler(K1vect, Uvect, fluxF, fluxG, &
                & deltax, deltay, deltat, imax, jmax)

            ! Second stage
            Call compute_Fluxes('x', time, K1vect, fluxK1F)
            Call compute_Fluxes('y', time, K1vect, fluxK1G)

            Call ExplicitEuler(K2vect, K1vect, fluxK1F, fluxK1G, &
                & deltax, deltay, deltat, imax, jmax)

            K2vect = .75_PR * Uvect + .25_PR * K2vect

            ! Third stage
            Call compute_Fluxes('x', time, K2vect, fluxK2F)
            Call compute_Fluxes('y', time, K2vect, fluxK2G)

            Call ExplicitEuler(K3vect, K2vect, fluxK2F, fluxK2G, &
                & deltax, deltay, deltat, imax, jmax)

            Uvect = 1._PR/3._PR * Uvect + 2._PR/3._PR * K3vect
        Case Default
            Write(STDERR, *) "Unknown time scheme ", TRIM(ADJUSTL(time_scheme_name))
            Call Exit(1)
        End Select
    End Subroutine step

    Subroutine ExplicitEuler(Ures, Uvect, fluxF, fluxG, deltax, deltay, deltat, imax, jmax)
        ! --- InOut
        Real(PR), Dimension(:,:,:), Intent(In) :: Uvect
        Real(PR), Intent(In) :: deltax, deltay, deltat
        Integer, Intent(In) :: imax, jmax
        Real(PR), Dimension(:,0:,0:), Intent(In) :: fluxF, fluxG
        Real(PR), Dimension(:,:,:), Intent(InOut) :: Ures
        ! --- Locals
        Integer :: i, j

        Do i=1, imax
            Do j=1, jmax
                Ures(:,i,j) = Uvect(:,i,j) &
                    & - deltat/deltax * (fluxF(:,i,j) - fluxF(:,i-1,j)) &
                    & - deltat/deltay * (fluxG(:,i,j) - fluxG(:,i,j-1))
            End Do
        End Do
    End Subroutine ExplicitEuler

    Subroutine compute_Fluxes(axis, t, Uvect, flux)
        ! --- InOut
        Character, Intent(In) :: axis
        Real(PR), Intent(In) :: t
        Real(PR), Dimension(:,:,:), Intent(In) :: Uvect
        Real(PR), Dimension(:,0:,0:), Intent(InOut) :: flux
        ! --- Locals
        Real(PR), Dimension(4) :: ULL, UL, UR, URR

        Select Case (axis)
        Case ('x')
            Do j=1, jmax
                Do i=2, imax-2
                    flux(:,i,j) = numericalFlux(numflux_name, 'x', space_scheme_specs, &
                        & Uvect(:,i-1,j), Uvect(:,i,j), Uvect(:,i+1,j), Uvect(:,i+2,j), &
                        & gammagp)
                End Do
                ! Boundary
                Call valuesAtBoundary(case_name, 'x', 0, j, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
                flux(:,0,j) = numericalFlux(numflux_name, 'x', space_scheme_specs, &
                    & ULL, UL, UR, URR, &
                    & gammagp)
                Call valuesAtBoundary(case_name, 'x', 1, j, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
                flux(:,1,j) = numericalFlux(numflux_name, 'x', space_scheme_specs, &
                    & ULL, UL, UR, URR, &
                    & gammagp)
                Call valuesAtBoundary(case_name, 'x', imax-1, j, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
                flux(:,imax-1,j) = numericalFlux(numflux_name, 'x', space_scheme_specs, &
                    & ULL, UL, UR, URR, &
                    & gammagp)
                Call valuesAtBoundary(case_name, 'x', imax, j, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
                flux(:,imax,j) = numericalFlux(numflux_name, 'x', space_scheme_specs, &
                    & ULL, UL, UR, URR, &
                    & gammagp)
            End Do
        Case ('y')
            Do i=1, imax
                Do j=2, jmax-2
                    flux(:,i,j) = numericalFlux(numflux_name, 'y', space_scheme_specs, &
                        & Uvect(:,i,j-1) ,Uvect(:,i,j), Uvect(:,i,j+1), Uvect(:,i,j+2), &
                        & gammagp)
                End Do
                ! Boundary
                Call valuesAtBoundary(case_name, 'y', i, 0, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
                flux(:,i,0) = numericalFlux(numflux_name, 'y', space_scheme_specs, &
                    & ULL, UL, UR, URR, &
                    & gammagp)
                Call valuesAtBoundary(case_name, 'y', i, 1, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
                flux(:,i,1) = numericalFlux(numflux_name, 'y', space_scheme_specs, &
                    & ULL, UL, UR, URR, &
                    & gammagp)
                Call valuesAtBoundary(case_name, 'y', i, jmax-1, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
                flux(:,i,jmax-1) = numericalFlux(numflux_name, 'y', space_scheme_specs, &
                    & ULL, UL, UR, URR, &
                    & gammagp)
                Call valuesAtBoundary(case_name, 'y', i, jmax, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
                flux(:,i,jmax) = numericalFlux(numflux_name, 'y', space_scheme_specs, &
                    & ULL, UL, UR, URR, &
                    & gammagp)
            End Do
        Case Default
            Write(STDERR, *) "Unknown axis ", axis
            Call Exit(1)
        End Select
    End Subroutine

    Subroutine compute_CFL(U, t, dx, dy, dt, cfl)
        ! --- InOut ---
        Real(PR), Dimension(:,:,:), Intent(In) :: U
        Real(PR), Intent(In) :: t
        Real(PR), Intent(In) :: dx, dy,  cfl
        Real(PR), Intent(Out) :: dt
        ! --- Locals ---
        Real(PR), Dimension(4) :: ULi, URi, ULL, UL, UR, URR
        Real(PR) :: bx_max, by_max
        Integer :: i, j

        bx_max = 0._PR
        by_max = 0._PR

        ! --- x
        Do j=1, jmax
            Do i=2, imax-2
                Call reconstructAtInterface('x', ULi, URi, space_scheme_specs, &
                    & U(:,i-1,j), U(:,i,j), U(:,i+1,j), U(:,i+2,j) &
                    & )
                Call compute_bmax( 'x', ULi, gammagp, bx_max )
                Call compute_bmax( 'x', URi, gammagp, bx_max )
            End Do
            ! Boundary
            Call valuesAtBoundary(case_name, 'x', 0, j, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
            Call reconstructAtInterface('x', ULi, URi, space_scheme_specs, &
                & ULL, UL, UR, URR &
                & )
            Call compute_bmax( 'x', ULi, gammagp, bx_max )
            Call compute_bmax( 'x', URi, gammagp, bx_max )

            Call valuesAtBoundary(case_name, 'x', 1, j, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
            Call reconstructAtInterface('x', ULi, URi, space_scheme_specs, &
                & ULL, UL, UR, URR &
                & )
            Call compute_bmax( 'x', ULi, gammagp, bx_max )
            Call compute_bmax( 'x', URi, gammagp, bx_max )

            Call valuesAtBoundary(case_name, 'x', imax-1, j, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
            Call reconstructAtInterface('x', ULi, URi, space_scheme_specs, &
                & ULL, UL, UR, URR &
                & )
            Call compute_bmax( 'x', ULi, gammagp, bx_max )
            Call compute_bmax( 'x', URi, gammagp, bx_max )

            Call valuesAtBoundary(case_name, 'x', imax, j, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
            Call reconstructAtInterface('x', ULi, URi, space_scheme_specs, &
                & ULL, UL, UR, URR &
                & )
            Call compute_bmax( 'x', ULi, gammagp, bx_max )
            Call compute_bmax( 'x', URi, gammagp, bx_max )
        End Do

        ! --- y
        Do i=1, imax
            Do j=2, jmax-2
                Call reconstructAtInterface('y', ULi, URi, space_scheme_specs, &
                    & U(:,i,j-1) ,U(:,i,j), U(:,i,j+1), U(:,i,j+2) &
                    )
                Call compute_bmax( 'y', ULi, gammagp, by_max )
                Call compute_bmax( 'y', URi, gammagp, by_max )
            End Do
            ! Boundary
            Call valuesAtBoundary(case_name, 'y', i, 0, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
            Call reconstructAtInterface('y', ULi, URi, space_scheme_specs, &
                & ULL, UL, UR, URR &
                & )
            Call compute_bmax( 'y', ULi, gammagp, by_max )
            Call compute_bmax( 'y', URi, gammagp, by_max )

            Call valuesAtBoundary(case_name, 'y', i, 1, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
            Call reconstructAtInterface('y', ULi, URi, space_scheme_specs, &
                & ULL, UL, UR, URR &
                & )
            Call compute_bmax( 'y', ULi, gammagp, by_max )
            Call compute_bmax( 'y', URi, gammagp, by_max )

            Call valuesAtBoundary(case_name, 'y', i, jmax-1, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
            Call reconstructAtInterface('y', ULi, URi, space_scheme_specs, &
                & ULL, UL, UR, URR &
                & )
            Call compute_bmax( 'y', ULi, gammagp, by_max )
            Call compute_bmax( 'y', URi, gammagp, by_max )

            Call valuesAtBoundary(case_name, 'y', i, jmax, Uvect, deltax, deltay, t, ULL, UL, UR, URR, gammagp);
            Call reconstructAtInterface('y', ULi, URi, space_scheme_specs, &
                & ULL, UL, UR, URR &
                & )
            Call compute_bmax( 'y', ULi, gammagp, by_max )
            Call compute_bmax( 'y', URi, gammagp, by_max )
        End Do

        dt = cfl * MIN(dx/bx_max, dy/by_max)

    End Subroutine compute_CFL

    Subroutine compute_bmax( axis, Ui, gammagp, bmax )
        ! --- InOut
        Character, Intent(In) :: axis
        Real(PR), Dimension(4), Intent(In) :: Ui
        Real(PR), Intent(In) :: gammagp
        Real(PR), Intent(InOut) :: bmax
        ! --- Locals
        Real(PR) :: rho, u, v, e, q, p, a, b, l1, l3, velocity

        rho = Ui(1)
        u = Ui(2) / rho
        v = Ui(3) / rho
        e = Ui(4)

        q = .5_PR * ( u**2 + v**2 )
        p = (gammagp - 1._PR)*(e - rho*q)
        a = SQRT(gammagp*p/rho)

        Select Case (axis)
        Case ('x')
            velocity = u
        Case ('y')
            velocity = v
        Case Default
            Write(STDERR, *) "Unknown axis ", axis
            Call Exit(1)
        End Select

        l1 = ABS(velocity - a)
        l3 = ABS(velocity + a)
        b = MAX(l1, l3)
        bmax = MAX( bmax, b )
    End Subroutine compute_bmax

    Function error(norm, case_name, U, time, gammagp, compute_on_boundary)
        ! --- InOut
        Real(PR), Dimension(4,imax,jmax), Intent(In) :: U
        Real(PR), Intent(In) :: gammagp, time
        Character(len=*), Intent(In) :: case_name
        Character(len=*), Intent(In) :: norm
        Logical, Intent(In) :: compute_on_boundary
        Real(PR), Dimension(4) :: error
        ! --- Locals
        Real(PR), Dimension(4) :: exact_value
        Integer :: i_range_min, i_range_max, j_range_min, j_range_max

        If (compute_on_boundary) Then
            i_range_min = 1
            i_range_max = imax
            j_range_min = 1
            j_range_max = jmax
        Else
            i_range_min = 2
            i_range_max = imax-2
            j_range_min = 2
            j_range_max = jmax-2
        End If

        error = 0._PR
        Do i=i_range_min, i_range_max
            Do j=j_range_min, j_range_max
                exact_value = Uexact(case_name, xm(i), ym(j), time, gammagp)
                Select Case (TRIM(ADJUSTL(norm)))
                Case ('Linfty') ! L_infinity norm
                    error = MAX( error, ABS( U(:,i,j) - exact_value ) )
                Case ('L1') ! L1 norm
                    error = error + ABS( U(:,i,j) - exact_value )
                Case ('L2') ! L2 norm
                    error = error + ( U(:,i,j) - exact_value )**2
                Case Default
                    Write(STDERR,*) "Unknown norm ", norm
                    Call Exit(1)
                End Select
            End Do
        End Do

        Select Case (TRIM(ADJUSTL(norm)))
        Case ('Linfty') ! L_infinity norm
        Case ('L1') ! L1 norm
            error = error / ( imax * jmax )
        Case ('L2') ! L2 norm
            error = SQRT( error / ( imax * jmax ) )
        Case Default
            Write(STDERR,*) "Unknown norm ", norm
            Call Exit(1)
        End Select
    End Function error

End Program euler
