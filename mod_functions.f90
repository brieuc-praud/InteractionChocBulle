Module mod_functions

    Use mod_parameters
    Implicit None

Contains
    Function Uinit(case_number, x, y, gammagp)
        ! --- InOut
        Real(PR), Intent(In) :: x, y, gammagp
        Integer, Intent(In) :: case_number
        Real(PR), Dimension(4) :: Uinit
        ! --- Locals
        Real(PR) :: r, u, v, q, p, e


        Select Case (case_number)
        Case (1) ! Shock tube
            If (x < .5_PR) Then
                ! Left
                r = 1._PR
                u = 0._PR
                v = 0._PR
                p = 1._PR
            Else
                ! Right
                r = .125_PR
                u = 0._PR
                v = 0._PR
                p = .1_PR
            End If
        Case (2) ! Isentropic vortex
            Call IsentropicVortex(x, y, 0._PR, &
                & 45._PR, &
                & SQRT(2._PR/gammagp), &
                & 1._PR, 1._PR, &
                & 1._PR, 5._PR, 1._PR, &
                & SQRT(EXP(1._PR)/gammagp)*5._PR/(2._PR*PI), &
                & gammagp, &
                & r, u, v, p)
        Case Default
            r = 1._PR
            u = 0._PR
            v = 0._PR
            p = 1._PR
        End Select

        q = .5_PR * ( u**2 + v**2 )
        e = p / (gammagp - 1._PR) + r*q

        Uinit = (/ Real(PR) :: r, r*u, r*v, e /)
    End Function Uinit

    Function Uexact(case_number, x, y, t, gammagp)
        ! --- InOut
        Real(PR), Intent(In) :: x, y, t, gammagp
        Integer, Intent(In) :: case_number
        Real(PR), Dimension(4) :: Uexact
        ! --- Locals
        Real(PR) :: r, u, v, q, p, e

        Select Case (case_number)
        Case (2) ! Isentropic vortex
            Call IsentropicVortex(x, y, t, &
                & 45._PR, &
                & SQRT(2._PR/gammagp), &
                & 1._PR, 1._PR, &
                & 1._PR, 5._PR, 1._PR, &
                & SQRT(EXP(1._PR)/gammagp)*5._PR/(2._PR*PI), &
                & gammagp, &
                & r, u, v, p)
        Case Default
            Write(STDERR, *) "No solution available for this case (case ", case_number, ")"
            Call Exit(1)
        End Select

        q = .5_PR * ( u**2 + v**2 )
        e = p / (gammagp - 1._PR) + r*q

        Uexact = (/ Real(PR) :: r, r*u, r*v, e /)
    End Function Uexact

    Subroutine IsentropicVortex(x, y, time, alpha_in_degrees, Mach_infty, density_infty, pressure_infty, &
            & vortexRadius, domain_half_length, std, perturbation_strength, gammagp, &
            & r, u, v, p)
        ! --- InOut
        Real(PR), Intent(In) :: x, y, time, alpha_in_degrees, Mach_infty, density_infty, pressure_infty, &
            & vortexRadius, domain_half_length, std, perturbation_strength, gammagp
        Real(PR), Intent(Out) :: r, u, v, p
        ! --- Locals
        Real(PR) :: f, disturbance, sound_of_speed_infty
        Real(PR) :: u_infty, v_infty
        Real(PR) :: x_old, y_old
        Real(PR) :: alpha

        sound_of_speed_infty = SQRT( gammagp * pressure_infty/density_infty )
        alpha = alpha_in_degrees * PI / 180._PR
        u_infty = sound_of_speed_infty * Mach_infty * COS( alpha )
        v_infty = sound_of_speed_infty * Mach_infty * SIN( alpha )

        x_old = MODULO( x - u_infty*time + domain_half_length, 2*domain_half_length ) - domain_half_length
        y_old = MODULO( y - v_infty*time + domain_half_length, 2*domain_half_length ) - domain_half_length
        f = ( x_old/vortexRadius )**2 + ( y_old/vortexRadius )**2
        f = - f / ( 2._PR * std**2 )
        disturbance = perturbation_strength * EXP( f )

        r = density_infty * ( 1._PR  -  .5_PR*(gammagp - 1._PR)*disturbance**2 )**(1._PR / (gammagp - 1._PR)) 
        u = u_infty - y/vortexRadius*disturbance
        v = v_infty + x/vortexRadius*disturbance
        p = pressure_infty / gammagp * ( 1._PR -.5_PR * ( gammagp - 1._PR ) * disturbance**2 )**(gammagp / (gammagp - 1._PR))
    End Subroutine IsentropicVortex
End Module mod_functions
