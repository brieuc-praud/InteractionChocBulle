Program main
    Use mod_input
    Use mod_output

    Use mod_cases
    Use mod_fluxes

    Implicit None

    ! --- Declare the variables
    Character(128) :: parameters_file ! Name of the parameters file
    Integer :: i, j ! Loop indices
    Integer :: nb_iterations
    

    ! --- Read the parameters
    Call GETARG(1, parameters_file)
    Call read_parameters(parameters_file)
    ! Get if an exact solution is available from the name of the case read in the parameters
    exact_solution_available = exactSolutionAvailable()

    ! --- Allocate
    Allocate(x(0:imax), y(0:jmax), xm(imax), ym(jmax))
    Allocate(Uvect(4,-1:imax+2,-1:jmax+2), fluxF(4,0:imax, 0:jmax), fluxG(4,0:imax, 0:jmax))
    ! Those arrays can be quite big so they are allocated only if necessary
    If (num_scheme%time_scheme_order >= 2) Then
        ! Strong-Stability preserving Runge-Kutta 2
        Allocate(K1vect(4,-1:imax+2,-1:jmax+2), fluxK1F(4,0:imax, 0:jmax), fluxK1G(4,0:imax, 0:jmax))
        Allocate(K2vect(4,-1:imax+2,-1:jmax+2), fluxK2F(4,0:imax, 0:jmax), fluxK2G(4,0:imax, 0:jmax))
    End If
    If (num_scheme%time_scheme_order >= 3) Then
        ! Strong-Stability preserving Runge-Kutta 3
        Allocate(K3vect(4,-1:imax+2,-1:jmax+2), fluxK3F(4,0:imax, 0:jmax), fluxK3G(4,0:imax, 0:jmax))
    End If
    If (exact_solution_available) Then
        Allocate(Uvect_e(4,-1:imax+2,-1:jmax+2))
    End If

    ! --- Grid construction
    Call getGridDimensions(xmin, xmax, ymin, ymax)
    deltax = (xmax - xmin) / imax
    deltay = (ymax - ymin) / jmax
    x = (/ Real(PR) :: (xmin + i*deltax, i=0, imax) /)
    y = (/ Real(PR) :: (ymin + j*deltay, j=0, jmax) /)
    xm = (/ Real(PR) :: (xmin + .5_PR*deltax + i*deltax, i=0, imax-1) /)
    ym = (/ Real(PR) :: (ymin + .5_PR*deltay + j*deltay, j=0, jmax-1) /)
    ! Initialise Uvect
    Call getInitState(quadrature_points_number)

    ! Output initial state
    Call output(Uvect, 0, 'sol')
    Call output(Uvect, 0, 'exact')

    ! Time loop
    nb_iterations = 0
    time = 0._PR
    Do While (time < time_max)
        Call step()

        ! Output a file if the output modulo corresponds OR if it is the last iteration
        If ( (output_modulo > 0 .AND. Modulo(nb_iterations, output_modulo) == 0) &
            & .OR. time > (time_max - .5_PR*deltat) ) Then

            Write(STDOUT, *) time, time_max

            If ( exact_solution_available ) Then
                Call getExactSolution(quadrature_points_number)
                Call output(Uvect_e, nb_iterations / ABS(output_modulo) + 1, 'exact')
            End If
            Call output(Uvect, nb_iterations / ABS(output_modulo) + 1, 'sol')
        End If
        nb_iterations = nb_iterations + 1
    End Do

    If (exact_solution_available) Then
        Write(STDOUT, *) "dx: ", deltax, ", Error:", error()
    End If

    Deallocate(x, y, xm, ym)
    Deallocate(Uvect, fluxF, fluxG)
    If (num_scheme%time_scheme_order >= 2) Then
        Deallocate(K1vect)
        Deallocate(K2vect)
    End If
    If (num_scheme%time_scheme_order >= 3) Then
        Deallocate(K3vect)
    End If
    If (exact_solution_available) Then
        Deallocate(Uvect_e)
    End If

Contains

    Subroutine step()
        Call fillGhosts(Uvect)
        Call compute_CFL()
        deltat = MIN( deltat, time_max - time ) ! Adjust the time step to end at time_max
        time = time + deltat

        Select Case (num_scheme%time_scheme_order)
        Case (1) ! Explicit Euler
            Call ExplicitEuler(Uvect, Uvect)
        Case (2) ! Strong-Stability Preserving Runge-Kutta 2
            ! First stage
            Call ExplicitEuler(K1vect, Uvect)
            ! Second stage
            Call fillGhosts(K1vect)
            Call ExplicitEuler(K2vect, K1vect)
            Uvect = .5_PR * Uvect + .5_PR * K2vect
        Case (3) ! Strong-Stability Preserving Runge-Kutta 3
            ! First stage
            Call ExplicitEuler(K1vect, Uvect)
            ! Second stage
            Call fillGhosts(K1vect)
            Call ExplicitEuler(K2vect, K1vect)
            K2vect = .75_PR * Uvect + .25_PR * K2vect
            ! Third stage
            Call fillGhosts(K2vect)
            Call ExplicitEuler(K3vect, K2vect)
            Uvect = 1._PR/3._PR * Uvect + 2._PR/3._PR * K3vect
        Case Default
            Write(STDERR, *) "Unsupported time scheme order ", num_scheme%time_scheme_order
            Call Exit(1)
        End Select
    End Subroutine step

    Subroutine ExplicitEuler(Up, Um)
        ! --- InOut
        Real(PR), Dimension(4,-1:imax+2,-1:jmax+2), Intent(In) :: Um
        Real(PR), Dimension(4,-1:imax+2,-1:jmax+2), Intent(InOut) :: Up
        ! --- Locals
        Integer :: i, j

        Call compute_Fluxes(Um, fluxF, fluxG)

        Do j=1, jmax
            Do i=1, imax
                Up(:,i,j) = Um(:,i,j) &
                    & - deltat/deltax * (fluxF(:,i,j) - fluxF(:,i-1,j)) &
                    & - deltat/deltay * (fluxG(:,i,j) - fluxG(:,i,j-1))! &
            End Do
        End Do
    End Subroutine ExplicitEuler

    Subroutine compute_Fluxes(U, flux_x, flux_y)
        ! --- InOut
        Real(PR), Dimension(4,-1:imax+2,-1:jmax+2), Intent(In) :: U
        Real(PR), Dimension(4,0:imax,0:jmax), Intent(InOut) :: flux_x, flux_y

        Do j=1, jmax
            Do i=1, imax
                flux_x(:,i,j) = numericalFlux('x', U(:,i-1,j), U(:,i,j), U(:,i+1,j), U(:,i+2,j))
                flux_y(:,i,j) = numericalFlux('y', U(:,i,j-1), U(:,i,j), U(:,i,j+1), U(:,i,j+2))
            End Do
        End Do
        ! Boundary
        Do j=1, jmax
            flux_x(:,0,j) = numericalFlux('x', U(:,-1,j), U(:,0,j), U(:,1,j), U(:,2,j))
        End Do
        Do i=1, imax
            flux_y(:,i,0) = numericalFlux('y', U(:,i,-1), U(:,i,0), U(:,i,1), U(:,i,2))
        End Do
    End Subroutine

    Subroutine compute_CFL()
        ! --- Locals ---
        Real(PR), Dimension(4) :: ULi, URi
        Real(PR) :: bx_max, by_max
        Integer :: i, j

        bx_max = 0._PR
        by_max = 0._PR

        ! --- x
        Do j=1, jmax
            Do i=1, imax
                Call reconstructAtInterface('x', ULi, URi, &
                    & Uvect(:,i-1,j), Uvect(:,i,j), Uvect(:,i+1,j), Uvect(:,i+2,j) )
                Call compute_bmax( 'x', ULi, bx_max )
                Call compute_bmax( 'x', URi, bx_max )
            End Do
        End Do

        ! --- y
        Do j=1, jmax
            Do i=1, imax
                Call reconstructAtInterface('y', ULi, URi, &
                    & Uvect(:,i,j-1) ,Uvect(:,i,j), Uvect(:,i,j+1), Uvect(:,i,j+2) )
                Call compute_bmax( 'y', ULi, by_max )
                Call compute_bmax( 'y', URi, by_max )
            End Do
        End Do

        deltat = cfl * MIN(deltax/bx_max, deltay/by_max)
    End Subroutine compute_CFL

    Subroutine compute_bmax(axis, Ui, bmax)
        ! --- InOut
        Character, Intent(In) :: axis
        Real(PR), Dimension(4), Intent(In) :: Ui
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

    Subroutine getExactSolution(nb_quadrature_points)
        ! --- InOut
        Integer, Intent(In) :: nb_quadrature_points
        ! --- Locals
        Integer :: i, j

        Do i=1, imax
            Do j=1, jmax
                Uvect_e(:,i,j) = Uexact_avg(xm(i), ym(j), time, nb_quadrature_points)
            End Do
        End Do
    End Subroutine getExactSolution

    Subroutine getInitState(nb_quadrature_points)
        ! --- InOut
        Integer, Intent(In) :: nb_quadrature_points
        ! --- Locals
        Integer :: i, j

        Do i=1, imax
            Do j=1, jmax
                Uvect(:,i,j) = Uinit_avg(xm(i), ym(j), nb_quadrature_points)
                Uvect(:,i,j) = Uinit(xm(i), ym(j))
            End Do
        End Do
    End Subroutine getInitState

    Function error()
        ! --- InOut
        Real(PR), Dimension(4) :: error
        ! --- Locals
        Integer :: nb_cells

        Call getExactSolution(quadrature_points_number)

        nb_cells = imax*jmax
        error = 0._PR
        Do i=1, imax
            Do j=1, jmax
                Select Case (TRIM(ADJUSTL(norm_str)))
                Case ('Linfty') ! L_infinity norm
                    error = MAX( error, ABS( Uvect(:,i,j) - Uvect_e(:,i,j) ) )
                Case ('L1') ! L1 norm
                    error = error + ABS( Uvect(:,i,j) - Uvect_e(:,i,j) )
                Case ('L2') ! L2 norm
                    error = error + ( Uvect(:,i,j) - Uvect_e(:,i,j) )**2
                Case Default
                    Write(STDERR,*) "Unknown norm ", norm_str
                    Call Exit(1)
                End Select
            End Do
        End Do

        Select Case (TRIM(ADJUSTL(norm_str)))
        Case ('Linfty') ! L_infinity norm
        Case ('L1') ! L1 norm
            error = error / nb_cells
        Case ('L2') ! L2 norm
            error = SQRT( error / nb_cells )
        Case Default
            Write(STDERR,*) "Unknown norm ", norm_str
            Call Exit(1)
        End Select
    End Function error

End Program main
