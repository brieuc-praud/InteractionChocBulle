Program euler
    Use mod_functions
    Use mod_schemes
    Use mod_output
    Use mod_test

    Implicit None

    Character(14), Parameter :: parameters = "parameters.dat"

    ! Parameters
    Real(PR) xmin, xmax, ymin, ymax, time_max, cfl, gammagp
    Integer :: imax, jmax
    Integer :: noutput
    ! Arrays
    Real(PR), Dimension(:,:,:), Allocatable :: Uvect, Uvect_e, fluxF, fluxG
    Real(PR), Dimension(:), Allocatable :: x, y, xm, ym
    ! Loop indices
    Integer :: i, j, nb_iterations
    ! Other
    Real(PR) :: deltax, deltay, deltat, time

    ! Read parameters
    Open(111, File=parameters)
    Read(111, *) xmin, xmax, ymin, ymax, time_max, cfl, imax, jmax, gammagp, noutput
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
        Do j=1 , jmax
            Uvect(:,i,j) = Uinit(2, xm(i), ym(j), gammagp)
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
                fluxF(:,i,j) = Rusanov('x', Uvect(:,i,j), Uvect(:,i+1,j), gammagp)
                !fluxF(:,i,j) = HLL('x', Uvect(:,i,j), Uvect(:,i+1,j), gammagp)
            End Do
            ! DIRICHLET
            !F(0,j) = -2._PR*D(xmin, ym(j)) * (T(1, j) - Touest(tn, ym(j))) / deltax
            !F(imax,j) = 2._PR*D(xmax, ym(j)) * (T(imax,j) - Test(tn,ym(j))) / deltax
            ! NEUMANN
            If (.True.) Then
                ! Periodic
                fluxF(:,0,j) = Rusanov('x', Uvect(:,imax,j), Uvect(:,1,j), gammagp)
                fluxF(:,imax,j) = fluxF(:,0,j)
            Else If (.False.) Then
                ! Absorbing
                fluxF(:,0,j) = fluxFuncF( Uvect(:,1,j), gammagp )
                fluxF(:,imax,j) = fluxFuncF( Uvect(:,imax,j), gammagp )
            Else
                ! Homogeneous Neumann
                fluxF(:,0,j) = 0._PR
                fluxF(:,imax,j) = 0._PR
            End If
        End Do
        Do i=1, imax
            Do j=1, jmax-1
                fluxG(:,i,j) = Rusanov('y', Uvect(:,i,j), Uvect(:,i,j+1), gammagp)
                !fluxG(:,i,j) = HLL('y', Uvect(:,i,j), Uvect(:,i,j+1), gammagp)
            End Do
            ! DIRICHLET
            !G(i,0) = -2._PR*D(xm(i), ymin) * ( T(i,1) - Tnord(tn,xm(i)) ) / deltay
            !G(i,jmax) = 2._PR*D(xm(i), ymax) * (T(i,jmax) - Tsud(tn,xm(i))) / deltay
            ! NEUMANN
            If (.True.) Then
                ! Periodic
                fluxG(:,i,0) = Rusanov('y', Uvect(:,i,jmax), Uvect(:,i,1), gammagp)
                fluxG(:,i,jmax) = fluxG(:,i,0)
            Else If (.False.) Then
                ! Absorbing
                fluxG(:,i,0) = fluxFuncG( Uvect(:,i,0), gammagp )
                fluxG(:,i,jmax) = fluxFuncG( Uvect(:,i,jmax), gammagp )
            Else
                ! Homogeneous Neumann
                fluxG(:,i,0) = 0._PR
                fluxG(:,i,jmax) = 0._PR
            End If
        End Do

        !Do i=1, imax
        !    Do j=1, jmax
        !            Write(*,*) i, j, "F", -1._PR/deltax * (fluxF(:,i,j) - fluxF(:,i-1,j))
        !            Write(*,*) i, j, "F1", fluxF(:,i,j)
        !            Write(*,*) i, j, "F2", fluxF(:,i-1,j)
        !            Write(*,*) Uvect(:,i,j)
        !            Write(*,*) Uvect(:,i+1,j)
        !            Write(*,*) i, j, "G", -1._PR/deltay * (fluxG(:,i,j) - fluxG(:,i,j-1))
        !    End Do
        !End Do
        !Call Exit()
        Do i=1, imax
            Do j=1, jmax
                Uvect(:,i,j) = Uvect(:,i,j) &
                    & - deltat/deltax * (fluxF(:,i,j) - fluxF(:,i-1,j)) &
                    & - deltat/deltay * (fluxG(:,i,j) - fluxG(:,i,j-1))
            End Do
        End Do

        If ( Modulo(nb_iterations, noutput) == 0) Then
            Write(STDOUT, *) time, time_max
            Do i=1, imax
                Do j=1 , jmax
                    Uvect_e(:,i,j) = Uexact(2, xm(i), ym(j), time, gammagp)
                End Do
            End Do
            Call output(Uvect, gammagp, x, y, nb_iterations / noutput + 1, 'sol')
            Call output(Uvect_e, gammagp, x, y, nb_iterations / noutput + 1, 'exact')
        End If
        nb_iterations = nb_iterations + 1
    End Do

    Write(STDOUT, *) "Error:", error(2, Uvect, gammagp)

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

    Function error(case_number, U, gammagp)
        ! --- InOut
        Real(PR), Dimension(4,imax,jmax), Intent(In) :: U
        Real(PR), Intent(In) :: gammagp
        Integer, Intent(In) :: case_number
        Real(PR), Dimension(4) :: error
        ! --- Locals
        Real(PR), Dimension(4) :: maximum, exact_value
        
        maximum = Huge( U(:,1,1) )
        Do i=1, imax
            Do j=1 , jmax
                exact_value = Uexact(case_number, xm(i), ym(j), time_max, gammagp)
                error = MAX( error, ABS( Uvect(:,i,j) - exact_value ) )
            End Do
        End Do
    End Function error

End Program euler
