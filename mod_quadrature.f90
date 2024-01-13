Module mod_quadrature

    Use mod_cases
    Implicit None

Contains

    Function Uexact_avg(x, y, t, nb_quadrature_points)
        ! --- InOut
        Real(PR), Intent(In) :: x, y, t
        Integer, Intent(In) :: nb_quadrature_points
        Real(PR), Dimension(4) :: Uexact_avg
        ! --- Locals
        Real(PR) :: p, px, py

        Select Case (nb_quadrature_points)
        Case (1)
            Uexact_avg = Uexact(x, y, t)
        Case (4)
            p = SQRT(3._PR)/3._PR
            px = p*deltax
            py = p*deltay
            Uexact_avg = Uexact(x-px , y-py, t)
            Uexact_avg = Uexact_avg + Uexact(x-px , y+py, t)
            Uexact_avg = Uexact_avg + Uexact(x+px , y-py, t)
            Uexact_avg = Uexact_avg + Uexact(x+px , y+py, t)
        Case Default
            Write(STDERR, *) "No quadrature formula available with ", nb_quadrature_points, " points"
            Call Exit(1)
        End Select

        Uexact_avg = Uexact_avg/Real(nb_quadrature_points, PR)

    End Function Uexact_avg

    Function Uinit_avg(x, y, nb_quadrature_points)
        ! --- InOut
        Real(PR), Intent(In) :: x, y
        Integer, Intent(In) :: nb_quadrature_points
        Real(PR), Dimension(4) :: Uinit_avg
        ! --- Locals
        Real(PR) :: p, px, py

        Select Case (nb_quadrature_points)
        Case (1)
            Uinit_avg = Uinit(x, y)
        Case (4)
            p = SQRT(3._PR)/3._PR
            px = p*deltax
            py = p*deltay
            Uinit_avg = Uinit(x-px , y-py)
            Uinit_avg = Uinit_avg + Uinit(x-px , y+py)
            Uinit_avg = Uinit_avg + Uinit(x+px , y-py)
            Uinit_avg = Uinit_avg + Uinit(x+px , y+py)
        Case Default
            Write(STDERR, *) "No quadrature formula available with ", nb_quadrature_points, " points"
            Call Exit(1)
        End Select

        Uinit_avg = Uinit_avg/Real(nb_quadrature_points, PR)

    End Function Uinit_avg


End Module mod_quadrature
