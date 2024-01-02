Module mod_cases

    Use mod_parameters
    Implicit None

    Type :: vortex_parameters
        Real(PR) :: gammagp
        Real(PR) :: alpha_in_degrees
        Real(PR) :: Mach_infty
        Real(PR) :: density_infty
        Real(PR) :: pressure_infty
        Real(PR) :: radius
        Real(PR) :: domain_half_length
        Real(PR) :: std
        Real(PR) :: perturbation_strength
    End Type vortex_parameters
Contains

    Function Uinit(case_name, x, y, gammagp)
        ! --- InOut
        Character(len=*), Intent(In) :: case_name
        Real(PR), Intent(In) :: x, y, gammagp
        Real(PR), Dimension(4) :: Uinit
        ! --- Locals
        Real(PR) :: r, u, v, q, p, e
        Type(vortex_parameters) :: vtx

        Select Case (TRIM(ADJUSTL(case_name)))
        Case ('ShockTube') ! 
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
        Case ('IsentropicVortex')
            Call set_vortex_parameters(vtx, gammagp)
            Call IsentropicVortex(x, y, 0._PR, &
                & vtx, &
                & r, u, v, p)
        Case Default
            Write(STDERR, *) "Unknown case ", TRIM(ADJUSTL(case_name))
            Call Exit(1)
        End Select

        q = .5_PR * ( u**2 + v**2 )
        e = p / (gammagp - 1._PR) + r*q

        Uinit = (/ Real(PR) :: r, r*u, r*v, e /)
    End Function Uinit

    Function Uexact(case_name, x, y, t, gammagp)
        ! --- InOut
        Character(len=*), Intent(In) :: case_name
        Real(PR), Intent(In) :: x, y, t, gammagp
        Real(PR), Dimension(4) :: Uexact
        ! --- Locals
        Real(PR) :: r, u, v, q, p, e
        Type(vortex_parameters) :: vtx

        Select Case (TRIM(ADJUSTL(case_name)))
        Case ('IsentropicVortex')
            Call set_vortex_parameters(vtx, gammagp)
            Call IsentropicVortex(x, y, t, &
                & vtx, &
                & r, u, v, p)
        Case Default
            Write(STDERR, *) "No solution available for case ", TRIM(ADJUSTL(case_name))
            Call Exit(1)
        End Select

        q = .5_PR * ( u**2 + v**2 )
        e = p / (gammagp - 1._PR) + r*q

        Uexact = (/ Real(PR) :: r, r*u, r*v, e /)
    End Function Uexact

    Subroutine getGridDimensions(case_name, xmin, xmax, ymin, ymax)
        ! --- InOut
        Character(len=*), Intent(In) :: case_name
        Real(PR), Intent(InOut) :: xmin, xmax, ymin, ymax
        ! --- Locals
        Type(vortex_parameters) :: vtx

        Select Case (TRIM(ADJUSTL(case_name)))
        Case ('IsentropicVortex')
            Call set_vortex_parameters(vtx, 1.4_PR)! The parameter 1.4 doesn't matter here
            xmin = -vtx%domain_half_length
            xmax = vtx%domain_half_length
            ymin = -vtx%domain_half_length
            ymax = vtx%domain_half_length
        Case Default
            Write(STDERR, *) "Unknown case ", TRIM(ADJUSTL(case_name))
            Call Exit(1)
        End Select
    End Subroutine getGridDimensions

    ! ===== Vortex =====
    Subroutine IsentropicVortex(x, y, time, &
            vtx, &
            & r, u, v, p)
        ! --- InOut
        Real(PR), Intent(In) :: x, y, time
        Type(vortex_parameters), Intent(In) :: vtx
        Real(PR), Intent(Out) :: r, u, v, p
        ! --- Locals
        Real(PR) :: f, disturbance, sound_of_speed_infty
        Real(PR) :: u_infty, v_infty
        Real(PR) :: x_old, y_old
        Real(PR) :: alpha

        sound_of_speed_infty = SQRT( vtx%gammagp * vtx%pressure_infty/vtx%density_infty )
        alpha = vtx%alpha_in_degrees * PI / 180._PR
        u_infty = sound_of_speed_infty * vtx%Mach_infty * COS( alpha )
        v_infty = sound_of_speed_infty * vtx%Mach_infty * SIN( alpha )

        x_old = MODULO( x - u_infty*time + vtx%domain_half_length, 2*vtx%domain_half_length ) - vtx%domain_half_length
        y_old = MODULO( y - v_infty*time + vtx%domain_half_length, 2*vtx%domain_half_length ) - vtx%domain_half_length
        f = ( x_old/vtx%radius )**2 + ( y_old/vtx%radius )**2
        f = - f / ( 2._PR * vtx%std**2 )
        disturbance = vtx%perturbation_strength * EXP( f )

        r = vtx%density_infty * ( 1._PR  -  .5_PR*(vtx%gammagp - 1._PR)*disturbance**2 )**(1._PR / (vtx%gammagp - 1._PR)) 
        u = u_infty - y/vtx%radius*disturbance
        v = v_infty + x/vtx%radius*disturbance
        p = vtx%pressure_infty / vtx%gammagp * &
            & ( 1._PR -.5_PR * ( vtx%gammagp - 1._PR ) * disturbance**2 )**(vtx%gammagp / (vtx%gammagp - 1._PR))
    End Subroutine IsentropicVortex

    Subroutine set_vortex_parameters(vtx, gammagp)
        Real(PR), Intent(In) :: gammagp
        Type(vortex_parameters), Intent(InOut) :: vtx

        vtx%gammagp = 1.4_PR
        vtx%alpha_in_degrees = 45._PR
        vtx%Mach_infty = SQRT(2._PR/gammagp) 
        vtx%density_infty = 1._PR
        vtx%pressure_infty = 1._PR
        vtx%radius = 1._PR
        vtx%domain_half_length = 5._PR
        vtx%std = 1._PR
        vtx%perturbation_strength = SQRT(EXP(1._PR)/gammagp)*5._PR/(2._PR*PI)
    End Subroutine
End Module mod_cases
