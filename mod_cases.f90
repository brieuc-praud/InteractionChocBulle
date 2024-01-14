Module mod_cases
    ! This modules contains the following cases:
    ! - IsentropicVortex: Shu, C.-W., "Essentially Non-oscillatory and Weighted Essentially Non-oscillatory Schemes for Hyperbolic Conservation Laws", 1998
    ! - ShockTube: Sod, G. A. (1978). "A Survey of Several Finite Difference Methods for Systems of Nonlinear Hyperbolic Conservation Laws"

    Use mod_fluxes

    Implicit None

    Type :: shock_tube_parameters
        Real(PR) :: gammagp
        Real(PR) :: densityL, densityR
        Real(PR) :: pressureL, pressureR
        Real(PR) :: interface_abscissa
        Real(PR) :: length, width
    End Type shock_tube_parameters

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
    Function exactSolutionAvailable()
        ! --- InOut
        Logical :: exactSolutionAvailable

        Select Case (TRIM(ADJUSTL(case_name)))
        Case ('ShockBubble')
            exactSolutionAvailable = .False.
        Case ('Manufactured', 'IsentropicVortex', 'ShockTube')
            exactSolutionAvailable = .True.
        Case Default
            Write(STDERR, *) "Unknown case ", TRIM(ADJUSTL(case_name))
            Call Exit(1)
        End Select
    End Function exactSolutionAvailable

    Function Uinit(x, y)
        ! --- InOut
        Real(PR), Intent(In) :: x, y
        Real(PR), Dimension(4) :: Uinit
        ! --- Locals
        Real(PR) :: r, u, v, p
        Real(PR) :: rho, ru, rv, e
        Type(vortex_parameters) :: vtx
        Type(shock_tube_parameters) :: st

        Select Case (TRIM(ADJUSTL(case_name)))
        Case ('Manufactured')
            Call manufactured(x, y, 0._PR, r, u, v, p)
        Case ('ShockBubble')
            If (x < 0._PR) Then
                ! Left
                r = 3.81_PR
                u = 2.85_PR
                v = 0._PR
                p = 10._PR
            Else
                ! Right
                If ( (x-.3_PR)**2 + y**2 < .2_PR**2 ) Then
                    r = .1_PR
                    u = 0._PR
                    v = 0._PR
                    p = 1._PR
                Else
                    r = 1._PR
                    u = 0._PR
                    v = 0._PR
                    p = 1._PR
                End If
            End If
        Case ('ShockTube')
            Call set_shock_tube_parameters(st, gammagp)
            If (x < st%interface_abscissa) Then
                ! Left
                r = st%densityL
                u = 0._PR
                v = 0._PR
                p = st%pressureL
            Else
                ! Right
                r = st%densityR
                u = 0._PR
                v = 0._PR
                p = st%pressureR
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

        Call primitiveToConservative( (/Real(PR) :: r, u, v, p/), rho, ru, rv, e )

        Uinit = (/ Real(PR) :: rho, ru, rv, e /)
    End Function Uinit

    Function Uexact(x, y, t)
        ! --- InOut
        Real(PR), Intent(In) :: x, y, t
        Real(PR), Dimension(4) :: Uexact
        ! --- Locals
        Real(PR) :: r, u, v, p
        Real(PR) :: rho, ru, rv, e
        Type(vortex_parameters) :: vtx
        Type(shock_tube_parameters) :: st

        Select Case (TRIM(ADJUSTL(case_name)))
        Case ('Manufactured')
            Call manufactured(x, y, t, r, u, v, p)
        Case ('ShockTube')
            Call set_shock_tube_parameters(st, gammagp)
            Call ShockTube(x, t, &
                & st, &
                & r, u, v, p)
        Case ('IsentropicVortex')
            Call set_vortex_parameters(vtx, gammagp)
            Call IsentropicVortex(x, y, t, &
                & vtx, &
                & r, u, v, p)
        Case Default
            Write(STDERR, *) "No solution available for case ", TRIM(ADJUSTL(case_name))
            Call Exit(1)
        End Select

        Call primitiveToConservative( (/Real(PR) :: r, u, v, p/), rho, ru, rv, e )

        Uexact = (/ Real(PR) :: rho, ru, rv, e /)
    End Function Uexact

    Function sourceTerm(x, y, t)
        ! --- InOut
        Real(PR), Intent(In) :: x, y, t
        Real(PR), Dimension(4) :: sourceTerm

        Select Case (TRIM(ADJUSTL(case_name)))
        Case ('Manufactured')
            sourceTerm = manufacturedSourceTerm(x, y, t)
        Case Default
            sourceTerm = 0._PR
        End Select
    End Function sourceTerm

    Subroutine boundaryFlux(tag, flux, boolean_flag)
        ! --- InOut
        Character(len=2), Intent(In) :: tag
        Real(PR), Dimension(4), Intent(InOut) :: flux
        Logical, Intent(InOut) :: boolean_flag

        Select Case (TRIM(ADJUSTL(case_name)))
        Case ('Manufactured', 'ShockBubble', 'ShockTube', 'IsentropicVortex')
            Select Case (tag)
            Case ('x')
                flux = 0._PR
            Case ('y')
                flux = 0._PR
            Case Default
                Write(STDERR, *) "Invalid axis ", tag
                Call Exit(1)
            End Select
            boolean_flag = .False. ! Use this flux ? No !
        Case Default
            Write(STDERR, *) "Unknown case ", TRIM(ADJUSTL(case_name))
            Call Exit(1)
        End Select
    End Subroutine boundaryFlux

    Subroutine fillGhosts(U)
        ! --- InOut
        Real(PR), Dimension(4,-1:imax+2,-1:jmax+2), Intent(InOut) :: U
        ! --- Locals
        Integer :: i, j

        Select Case (TRIM(ADJUSTL(case_name)))
        Case ('Manufactured')
            ! Enforce the solution in the ghosts cells
            Do j=1,jmax
                U(:,-1,j) =     Uexact_avg(xm(1)-2._PR*deltax, ym(j), time, quadrature_points_number)
                U(:,0,j) =      Uexact_avg(xm(1)-deltax, ym(j), time, quadrature_points_number)
                U(:,imax+1,j) = Uexact_avg(xm(imax)+deltax, ym(j), time, quadrature_points_number)
                U(:,imax+2,j) = Uexact_avg(xm(imax)+2._PR*deltax, ym(j), time, quadrature_points_number)
            End Do
            Do i=1,imax
                U(:,i,-1) =     Uexact_avg(xm(i), ym(1)-2._PR*deltay, time, quadrature_points_number)
                U(:,i,0) =      Uexact_avg(xm(i), ym(1)-deltay, time, quadrature_points_number)
                U(:,i,jmax+1) = Uexact_avg(xm(i), ym(jmax)+deltay, time, quadrature_points_number)
                U(:,i,jmax+2) = Uexact_avg(xm(i), ym(jmax)+2._PR*deltay, time, quadrature_points_number)
            End Do
        Case ('ShockBubble')
            Do j=1,jmax
                ! Inflow
                U(:,-1,j) = Uinit(-.1_PR, 0._PR)
                U(:,0,j) = Uinit(-.1_PR, 0._PR)
                ! Linear extrapolation
                U(:,imax+1,j) = 2._PR*U(:,imax,j) - U(:,imax-1,j)
                U(:,imax+2,j) = 2._PR*U(:,imax+1,j) - U(:,imax,j)
            End Do
            Do i=1,imax
                ! Reflecting / Symmetry
                U(:,i,-1) = U(:,i,2)
                U(3,i,-1) = -U(3,i,-1)
                U(:,i,0) = U(:,i,1)
                U(3,i,0) = -U(3,i,0)
                U(:,i,jmax+1) = U(:,i,jmax)
                U(3,i,jmax+1) = -U(3,i,jmax+1)
                U(:,i,jmax+2) = U(:,i,jmax-1)
                U(3,i,jmax+2) = -U(3,i,jmax+2)
            End Do
        Case ('ShockTube')
            ! Linear extrapolation
            Do j=1,jmax
                U(:,0,j) = 2._PR*U(:,1,j) - U(:,2,j)
                U(:,-1,j) = 2._PR*U(:,0,j) - U(:,1,j)
                U(:,imax+1,j) = 2._PR*U(:,imax,j) - U(:,imax-1,j)
                U(:,imax+2,j) = 2._PR*U(:,imax+1,j) - U(:,imax,j)
            End Do
            Do i=1,imax
                U(:,i,0) = 2._PR*U(:,i,1) - U(:,i,2)
                U(:,i,-1) = 2._PR*U(:,i,0) - U(:,i,1)
                U(:,i,jmax+1) = 2._PR*U(:,i,jmax) - U(:,i,jmax-1)
                U(:,i,jmax+2) = 2._PR*U(:,i,jmax+1) - U(:,i,jmax)
            End Do
        Case ('IsentropicVortex')
            ! Periodic
            Do j=1,jmax
                U(:,-1,j) = U(:,imax-1,j)
                U(:,0,j) = U(:,imax,j)
                U(:,imax+1,j) = U(:,1,j)
                U(:,imax+2,j) = U(:,2,j)
            End Do
            Do i=1,imax
                U(:,i,-1) = U(:,i,jmax-1)
                U(:,i,0) = U(:,i,jmax)
                U(:,i,jmax+1) = U(:,i,1)
                U(:,i,jmax+2) = U(:,i,2)
            End Do
        Case Default
            Write(STDERR, *) "Unknown case ", TRIM(ADJUSTL(case_name))
            Call Exit(1)
        End Select
    End Subroutine fillGhosts

    Subroutine getGridDimensions(xmin, xmax, ymin, ymax)
        ! --- InOut
        Real(PR), Intent(InOut) :: xmin, xmax, ymin, ymax
        ! --- Locals
        Type(vortex_parameters) :: vtx
        Type(shock_tube_parameters) :: st

        Select Case (TRIM(ADJUSTL(case_name)))
        Case ('Manufactured')
            xmin = -5._PR
            xmax = 5._PR
            ymin = -5._PR
            ymax = 5._PR
        Case ('ShockBubble')
            xmin = -.1_PR
            xmax = 1.6_PR
            ymin = 0._PR
            ymax = .5_PR
        Case ('ShockTube')
            Call set_shock_tube_parameters(st, gammagp)! The parameter 'gammagp' doesn't matter here
            xmin = 0._PR
            xmax = st%length
            ymin = -.5_PR*st%width
            ymax = .5_PR*st%width
        Case ('IsentropicVortex')
            Call set_vortex_parameters(vtx, gammagp)! The parameter 'gammagp' doesn't matter here
            xmin = -vtx%domain_half_length
            xmax = vtx%domain_half_length
            ymin = -vtx%domain_half_length
            ymax = vtx%domain_half_length
        Case Default
            Write(STDERR, *) "Unknown case ", TRIM(ADJUSTL(case_name))
            Call Exit(1)
        End Select
    End Subroutine getGridDimensions

    ! ===== Shock Tube =====
    Subroutine set_shock_tube_parameters(st, gammagp)
        Real(PR), Intent(In) :: gammagp
        Type(shock_tube_parameters), Intent(InOut) :: st

        st%gammagp = gammagp
        st%densityL = 1._PR
        st%densityR = .125_PR
        st%pressureL = 10._PR
        st%pressureR = 1._PR
        st%interface_abscissa = 3._PR

        st%length = 10._PR
        st%width = 1._PR
    End Subroutine set_shock_tube_parameters

    Subroutine ShockTube(x, time, &
            st, &
            & r, u, v, p)
        ! --- InOut
        Real(PR), Intent(In) :: x, time
        Type(shock_tube_parameters), Intent(In) :: st
        Real(PR), Intent(Out) :: r, u, v, p
        ! --- Parameters
        Real(PR), Parameter :: newton_eps = 1.e-12
        Integer, Parameter :: newton_max_iter = 20
        ! --- Locals
        Real(PR) :: cs1, cs5, c2, c3, c4
        Real(PR) :: alpha, beta
        Real(PR) :: u3, du3, u4, du4, p3
        Real(PR) :: density3
        Real(PR) :: x1, x2, x3, x4
        Real(PR) :: k
        Real(PR) :: iteration_count

        Select Case (time < 1.e3*EPSILON(Real(PR)))
        Case (.True.)
            If ( x < st%interface_abscissa ) Then
                r = st%densityL
                u = 0._PR
                p = st%pressureL
            Else
                r = st%densityR
                u = 0._PR
                p = st%pressureR
            End If
        Case Default
            ! see "Sod shock tube", Wikipedia
            cs1 = SQRT(st%gammagp*st%pressureL/st%densityL)
            cs5 = SQRT(st%gammagp*st%pressureR/st%densityR)

            alpha = ( st%gammagp - 1._PR ) / ( st%gammagp + 1._PR )
            beta = .5_PR * ( st%gammagp - 1._PR ) / st%gammagp

            ! Newton method to find the pressure of state 3
            p3 = .5_PR * (st%pressureL + st%pressureR)
            iteration_count = 0
            Do
                u4 = ( p3 - st%pressureR )*SQRT( (1._PR-alpha) / &
                    ( st%densityR * ( p3 + alpha*st%pressureR ) ) )
                u3 = ( st%pressureL**beta - p3**beta )*SQRT( &
                    & ( (1._PR-alpha**2)*st%pressureL**(1._PR/st%gammagp) &
                    & / ( alpha**2 * st%densityL ) ) )

                du4 = ( 1._PR + .5_PR*(st%pressureR - p3)/(p3 + alpha*st%pressureR)) &
                    & * SQRT( (1._PR-alpha)/( st%densityR * ( p3 + alpha*st%pressureR ) ) )
                du3 = -beta * p3**(beta-1._PR) * SQRT( &
                    & ( (1._PR-alpha**2)*st%pressureL**(1._PR/st%gammagp) &
                    & / ( alpha**2 * st%densityL ) ) )

                p3 = p3 - (u3 - u4)/(du3 - du4)

                iteration_count = iteration_count + 1
                If ( ABS(u3 - u4) < newton_eps &
                    .OR. iteration_count > newton_max_iter ) Exit
            End Do
            u4 = ( p3 - st%pressureR )*SQRT( (1._PR-alpha) / &
                ( st%densityR * ( p3 + alpha*st%pressureR ) ) )
            u3 = ( st%pressureL**beta - p3**beta )*SQRT( &
                & ( (1._PR-alpha**2)*st%pressureL**(1._PR/st%gammagp) &
                & / ( alpha**2 * st%densityL ) ) )
            c3 = 2._PR*cs1/(st%gammagp-1._PR) &
                & * (1._PR-(p3/st%pressureL)**(.5_PR*(st%gammagp-1._PR)/st%gammagp) )
            density3 = st%densityL*(p3/st%pressureL)**(1._PR/st%gammagp)
            c2 = c3 - SQRT( st%gammagp * p3/density3 )
            c4 = cs5*SQRT( .5_PR * (st%gammagp - 1._PR + p3/st%pressureR*(st%gammagp+1._PR) ) / st%gammagp )

            x1 = st%interface_abscissa - cs1 * time
            x2 = st%interface_abscissa + c2 * time
            x3 = st%interface_abscissa + c3 * time
            x4 = st%interface_abscissa + c4 * time
            If ( x < x1 ) Then
                r = st%densityL
                u = 0._PR
                p = st%pressureL
            Else If ( x < x2 ) Then
                u = 2._PR/(st%gammagp+1._PR)*( cs1 + ( x - st%interface_abscissa )/time )
                k = 1._PR - .5_PR*(st%gammagp-1._PR)*u/cs1
                r = st%densityL * k**(2._PR/(st%gammagp-1))
                p = st%pressureL * k**(2._PR*st%gammagp/(st%gammagp-1))
            Else If ( x < x3 ) Then
                r = density3
                u = u3
                p = p3
            Else If ( x < x4 ) Then
                r = st%densityR*( (p3 + alpha*st%pressureR) &
                    & / ( st%pressureR + alpha* p3 ) )
                u = u4
                p = p3
            Else
                r = st%densityR
                u = 0._PR
                p = st%pressureR
            End If
        End Select
        v = 0._PR
    End Subroutine ShockTube

    ! ===== Vortex =====
    Subroutine set_vortex_parameters(vtx, gammagp)
        Real(PR), Intent(In) :: gammagp
        Type(vortex_parameters), Intent(InOut) :: vtx

        vtx%gammagp = gammagp
        vtx%alpha_in_degrees = 45._PR
        vtx%Mach_infty = SQRT(2._PR/gammagp) 
        vtx%density_infty = 1._PR
        vtx%pressure_infty = 1._PR
        vtx%radius = 1._PR
        vtx%domain_half_length = 5._PR
        vtx%std = 1._PR
        vtx%perturbation_strength = SQRT(EXP(1._PR)/gammagp)*5._PR/(2._PR*PI)
    End Subroutine set_vortex_parameters

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

    Subroutine manufactured(x, y, time, r, u, v, p)
        ! --- InOut
        Real(PR), Intent(In) :: x, y, time
        Real(PR), Intent(Out) :: r, u, v, p
        ! --- Locals
        Real(PR) :: Lx, Ly, f

        Lx = xmax-xmin
        Ly = ymax-ymin
        f = COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + COS(time/time_max)

        r = 2._PR + f
        u = f
        v = f
        p = 2._PR + f
    End Subroutine manufactured

    Function manufacturedSourceTerm(x, y, t)
        ! --- InOut
        Real(PR), Intent(In) :: x, y, t
        Real(PR), Dimension(4) :: manufacturedSourceTerm
        ! --- Locals
        Real(PR) :: Lx, Ly

        Lx = xmax-xmin
        Ly = ymax-ymin

        ! Expressions obtained from SageMath 9.3
        manufacturedSourceTerm(1) = -(SIN(x/Lx)/Lx + y*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) &
            & + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR) &
            & - (SIN(y/Ly)/Ly + x*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) &
            & + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR) &
            & - (SIN(x/Lx)/Lx + y*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) &
            & + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly))) - &
            & (SIN(y/Ly)/Ly + x*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) &
            & + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly))) - (COS(t/time_max) &
            & + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR)*SIN(t/time_max)/time_max &
            & - (COS(t/time_max) + COS(x/Lx) + COS(y/Ly) &
            & + COS(x*y/(Lx*Ly)))*SIN(t/time_max)/time_max
        manufacturedSourceTerm(2) = -2._PR*(SIN(x/Lx)/Lx + y*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) &
            & + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR)*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) &
            & + COS(x*y/(Lx*Ly))) - 2._PR*(SIN(y/Ly)/Ly + x*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) &
            & + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR)*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) &
            & + COS(x*y/(Lx*Ly))) - (SIN(x/Lx)/Lx + y*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) + COS(x/Lx) &
            & + COS(y/Ly) + COS(x*y/(Lx*Ly)))**2 - (SIN(y/Ly)/Ly + x*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) &
            & + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)))**2 - 2._PR*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) &
            & + COS(x*y/(Lx*Ly)) + 2._PR)*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) &
            & + COS(x*y/(Lx*Ly)))*SIN(t/time_max)/time_max - (COS(t/time_max) + COS(x/Lx) + COS(y/Ly) &
            & + COS(x*y/(Lx*Ly)))**2*SIN(t/time_max)/time_max - SIN(t/time_max)/time_max - SIN(x/Lx)/Lx - SIN(y/Ly)/Ly &
            & - x*SIN(x*y/(Lx*Ly))/(Lx*Ly) - y*SIN(x*y/(Lx*Ly))/(Lx*Ly)
        manufacturedSourceTerm(3) = -2._PR*(SIN(x/Lx)/Lx + y*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) + COS(x/Lx) &
            & + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR)*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly))) &
            & - 2._PR*(SIN(y/Ly)/Ly + x*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) &
            & + COS(x*y/(Lx*Ly)) + 2._PR)*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly))) - (SIN(x/Lx)/Lx &
            & + y*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)))**2 - (SIN(y/Ly)/Ly &
            & + x*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)))**2 &
            & - 2._PR*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR)*(COS(t/time_max) &
            & + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)))*SIN(t/time_max)/time_max - (COS(t/time_max) + COS(x/Lx) &
            & + COS(y/Ly) + COS(x*y/(Lx*Ly)))**2*SIN(t/time_max)/time_max
        manufacturedSourceTerm(4) = -((COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR)*(COS(t/time_max) &
            & + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)))**2 + (COS(t/time_max) + COS(x/Lx) + COS(y/Ly) &
            & + COS(x*y/(Lx*Ly)) + 2._PR)/(gammagp - 1._PR) + COS(t/time_max) + COS(x/Lx) + COS(y/Ly) &
            & + COS(x*y/(Lx*Ly)) + 2._PR)*(SIN(x/Lx)/Lx + y*SIN(x*y/(Lx*Ly))/(Lx*Ly)) - ((COS(t/time_max) &
            & + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR)*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) &
            & + COS(x*y/(Lx*Ly)))**2 + (COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR)/(gammagp - 1._PR) &
            & + COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR)*(SIN(y/Ly)/Ly &
            & + x*SIN(x*y/(Lx*Ly))/(Lx*Ly)) - (2._PR*(SIN(x/Lx)/Lx + y*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) + COS(x/Lx) &
            & + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR)*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly))) + (SIN(x/Lx)/Lx &
            & + y*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)))**2 &
            & + (SIN(x/Lx)/Lx + y*SIN(x*y/(Lx*Ly))/(Lx*Ly))/(gammagp - 1._PR) + SIN(x/Lx)/Lx &
            & + y*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly))) - (2._PR*(SIN(y/Ly)/Ly &
            & + x*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR)*(COS(t/time_max) &
            & + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly))) + (SIN(y/Ly)/Ly + x*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) &
            & + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)))**2 + (SIN(y/Ly)/Ly + x*SIN(x*y/(Lx*Ly))/(Lx*Ly))/(gammagp - 1._PR) &
            & + SIN(y/Ly)/Ly + x*SIN(x*y/(Lx*Ly))/(Lx*Ly))*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly))) &
            & - (2._PR*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR)*(COS(t/time_max) + COS(x/Lx) &
            & + COS(y/Ly) + COS(x*y/(Lx*Ly)))*SIN(t/time_max)/time_max + (COS(t/time_max) + COS(x/Lx) + COS(y/Ly) &
            & + COS(x*y/(Lx*Ly)))**2*SIN(t/time_max)/time_max + SIN(t/time_max)/time_max &
            & + SIN(t/time_max)/((gammagp - 1._PR)*time_max))*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly))) &
            & - ((COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR)*(COS(t/time_max) + COS(x/Lx) + COS(y/Ly) &
            &  + COS(x*y/(Lx*Ly)))**2 + (COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR)/(gammagp - 1._PR) &
            & + COS(t/time_max) + COS(x/Lx) + COS(y/Ly) + COS(x*y/(Lx*Ly)) + 2._PR)*SIN(t/time_max)/time_max
    End Function manufacturedSourceTerm

    ! ======= Quadratures ==========
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

    Function sourceTerm_avg(x, y, t, nb_quadrature_points)
        ! --- InOut
        Real(PR), Intent(In) :: x, y, t
        Integer, Intent(In) :: nb_quadrature_points
        Real(PR), Dimension(4) :: sourceTerm_avg
        ! --- Locals
        Real(PR) :: p, px, py

        Select Case (nb_quadrature_points)
        Case (1)
            sourceTerm_avg = sourceTerm(x, y, t)
        Case (4)
            p = SQRT(3._PR)/3._PR
            px = p*deltax
            py = p*deltay
            sourceTerm_avg = sourceTerm(x-px , y-py, t)
            sourceTerm_avg = sourceTerm_avg + sourceTerm(x-px , y+py, t)
            sourceTerm_avg = sourceTerm_avg + sourceTerm(x+px , y-py, t)
            sourceTerm_avg = sourceTerm_avg + sourceTerm(x+px , y+py, t)
        Case Default
            Write(STDERR, *) "No quadrature formula available with ", nb_quadrature_points, " points"
            Call Exit(1)
        End Select

        sourceTerm_avg = sourceTerm_avg/Real(nb_quadrature_points, PR)

    End Function sourceTerm_avg

End Module mod_cases
