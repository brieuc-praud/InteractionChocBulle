Module mod_quadrature

    Use mod_cases
    Implicit None

Contains

    Function Uexact_avg(case_name, x, y, t, gammagp, deltax, deltay, nb_quadrature_points)
        ! --- InOut
        Character(len=*), Intent(In) :: case_name
        Real(PR), Intent(In) :: x, y, t, gammagp, deltax, deltay
        Integer, Intent(In) :: nb_quadrature_points
        Real(PR), Dimension(4) :: Uexact_avg
        ! --- Locals
        Real(PR) :: p, px, py

        Select Case (nb_quadrature_points)
        Case (1)
            Uexact_avg = Uexact(case_name, x, y, t, gammagp)
        Case (4)
            p = SQRT(3._PR)/3._PR
            px = p*deltax
            py = p*deltay
            Uexact_avg = Uexact(case_name, x-px , y-py, t, gammagp)
            Uexact_avg = Uexact_avg + Uexact(case_name, x-px , y+py, t, gammagp)
            Uexact_avg = Uexact_avg + Uexact(case_name, x+px , y-py, t, gammagp)
            Uexact_avg = Uexact_avg + Uexact(case_name, x+px , y+py, t, gammagp)
        Case Default
            Write(STDERR, *) "No quadrature formula available with ", nb_quadrature_points, " points"
            Call Exit(1)
        End Select

        Uexact_avg = Uexact_avg/Real(nb_quadrature_points, PR)

    End Function Uexact_avg

    Function Uinit_avg(case_name, x, y, gammagp, deltax, deltay, nb_quadrature_points)
        ! --- InOut
        Character(len=*), Intent(In) :: case_name
        Real(PR), Intent(In) :: x, y, gammagp, deltax, deltay
        Integer, Intent(In) :: nb_quadrature_points
        Real(PR), Dimension(4) :: Uinit_avg
        ! --- Locals
        Real(PR) :: p, px, py

        Select Case (nb_quadrature_points)
        Case (1)
            Uinit_avg = Uinit(case_name, x, y, gammagp)
        Case (4)
            p = SQRT(3._PR)/3._PR
            px = p*deltax
            py = p*deltay
            Uinit_avg = Uinit(case_name, x-px , y-py, gammagp)
            Uinit_avg = Uinit_avg + Uinit(case_name, x-px , y+py, gammagp)
            Uinit_avg = Uinit_avg + Uinit(case_name, x+px , y-py, gammagp)
            Uinit_avg = Uinit_avg + Uinit(case_name, x+px , y+py, gammagp)
        Case Default
            Write(STDERR, *) "No quadrature formula available with ", nb_quadrature_points, " points"
            Call Exit(1)
        End Select

        Uinit_avg = Uinit_avg/Real(nb_quadrature_points, PR)

    End Function Uinit_avg


End Module mod_quadrature
