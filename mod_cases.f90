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

    Function Uinit(case_name, x, y, gammagp)
        ! --- InOut
        Character(len=*), Intent(In) :: case_name
        Real(PR), Intent(In) :: x, y, gammagp
        Real(PR), Dimension(4) :: Uinit
        ! --- Locals
        Real(PR) :: r, u, v, q, p, e
        Type(vortex_parameters) :: vtx
        Type(shock_tube_parameters) :: st

        Select Case (TRIM(ADJUSTL(case_name)))
        Case ('ShockTube') ! 
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
        Type(shock_tube_parameters) :: st

        Select Case (TRIM(ADJUSTL(case_name)))
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

        q = .5_PR * ( u**2 + v**2 )
        e = p / (gammagp - 1._PR) + r*q

        Uexact = (/ Real(PR) :: r, r*u, r*v, e /)
    End Function Uexact

    Function numericalFluxAtBoundary(case_name, numflux_name, axis, &
            & i, j, space_scheme_specs, Uvect, gammagp)
        ! --- InOut
        Character(len=*), Intent(In) :: case_name, numflux_name
        Character, Intent(In) :: axis
        Integer, Intent(In) :: i,j
        Real(PR), Dimension(:,:,:), Intent(In) :: Uvect
        Real(PR), Intent(In) :: gammagp
        Type(space_scheme), Intent(In) :: space_scheme_specs
        Real(PR), Dimension(4) :: numericalFluxAtBoundary
        ! --- Locals
        Integer :: imax, jmax
        Integer, Dimension(3) :: shape_Uvect
        Real(PR), Dimension(4) :: ULL, UL, UR, URR

        shape_Uvect = SHAPE(Uvect)
        imax = shape_Uvect(2)
        jmax = shape_Uvect(3)

        Select Case (TRIM(ADJUSTL(case_name)))
        Case ('ShockTube')
            ! Absorbing
            Select Case (axis)
            Case ('x')
                If (i == 0) Then
                    ULL = Uvect(:,1,j)
                    UL = Uvect(:,1,j)
                    UR = Uvect(:,1,j)
                    URR = Uvect(:,2,j)
                Else If (i == 1) Then
                    ULL = Uvect(:,1,j)
                    UL = Uvect(:,1,j)
                    UR = Uvect(:,2,j)
                    URR = Uvect(:,3,j)
                Else If (i == imax-1) Then
                    ULL = Uvect(:,imax-2,j)
                    UL = Uvect(:,imax-1,j)
                    UR = Uvect(:,imax,j)
                    URR = Uvect(:,imax,j)
                Else If (i == imax) Then
                    ULL = Uvect(:,imax-1,j)
                    UL = Uvect(:,imax,j)
                    UR = Uvect(:,imax,j)
                    URR = Uvect(:,imax,j)
                Else
                    Write(STDERR, *) "Invalid index for boundary condition ", i
                    Call Exit(1)
                End If
                numericalFluxAtBoundary = numericalFlux(numflux_name, 'x', space_scheme_specs, &
                    & ULL, UL, UR, URR, &
                    & gammagp)
            Case ('y')
                If (j == 0) Then
                    ULL = Uvect(:,i,1)
                    UL = Uvect(:,i,1)
                    UR = Uvect(:,i,1)
                    URR = Uvect(:,i,2)
                Else If (j == 1) Then
                    ULL = Uvect(:,i,1)
                    UL = Uvect(:,i,1)
                    UR = Uvect(:,i,2)
                    URR = Uvect(:,i,3)
                Else If (j == jmax-1) Then
                    ULL = Uvect(:,i,jmax-2)
                    UL = Uvect(:,i,jmax-1)
                    UR = Uvect(:,i,jmax)
                    URR = Uvect(:,i,jmax)
                Else If (j == jmax) Then
                    ULL = Uvect(:,i,jmax-1)
                    UL = Uvect(:,i,jmax)
                    UR = Uvect(:,i,jmax)
                    URR = Uvect(:,i,jmax)
                Else
                    Write(STDERR, *) "Invalid index for boundary condition ", j
                    Call Exit(1)
                End If
                numericalFluxAtBoundary = numericalFlux(numflux_name, 'y', space_scheme_specs, &
                    & ULL, UL, UR, URR, &
                    & gammagp)
            Case Default
                Write(STDERR, *) "Unknown axis ", axis
                Call Exit(1)
            End Select
        Case ('IsentropicVortex')
            ! Periodic
            Select Case (axis)
            Case ('x')
                If (i == 0 .OR. i == imax) Then
                    ULL = Uvect(:,imax-1,j)
                    UL = Uvect(:,imax,j)
                    UR = Uvect(:,1,j)
                    URR = Uvect(:,2,j)
                Else If (i == 1) Then
                    ULL = Uvect(:,imax,j)
                    UL = Uvect(:,1,j)
                    UR = Uvect(:,2,j)
                    URR = Uvect(:,3,j)
                Else If (i == imax-1) Then
                    ULL = Uvect(:,imax-2,j)
                    UL = Uvect(:,imax-1,j)
                    UR = Uvect(:,imax,j)
                    URR = Uvect(:,1,j)
                Else
                    Write(STDERR, *) "Invalid index for boundary condition ", i
                    Call Exit(1)
                End If
                numericalFluxAtBoundary = numericalFlux(numflux_name, 'x', space_scheme_specs, &
                    & ULL, UL, UR, URR, &
                    & gammagp)
            Case ('y')
                If (j == 0 .OR. j == jmax) Then
                    ULL = Uvect(:,i,jmax-1)
                    UL = Uvect(:,i,jmax)
                    UR = Uvect(:,i,1)
                    URR = Uvect(:,i,2)
                Else If (j == 1) Then
                    ULL = Uvect(:,i,jmax)
                    UL = Uvect(:,i,1)
                    UR = Uvect(:,i,2)
                    URR = Uvect(:,i,3)
                Else If (j == jmax-1) Then
                    ULL = Uvect(:,i,jmax-2)
                    UL = Uvect(:,i,jmax-1)
                    UR = Uvect(:,i,jmax)
                    URR = Uvect(:,i,1)
                Else
                    Write(STDERR, *) "Invalid index for boundary condition ", j
                    Call Exit(1)
                End If
                numericalFluxAtBoundary = numericalFlux(numflux_name, 'y', space_scheme_specs, &
                    & ULL, UL, UR, URR, &
                    & gammagp)
            Case Default
                Write(STDERR, *) "Unknown axis ", axis
                Call Exit(1)
            End Select
        Case Default
            Write(STDERR, *) "Unknown case ", TRIM(ADJUSTL(case_name))
            Call Exit(1)
        End Select
    End Function numericalFluxAtBoundary

    Subroutine getGridDimensions(case_name, xmin, xmax, ymin, ymax)
        ! --- InOut
        Character(len=*), Intent(In) :: case_name
        Real(PR), Intent(InOut) :: xmin, xmax, ymin, ymax
        ! --- Locals
        Type(vortex_parameters) :: vtx
        Type(shock_tube_parameters) :: st

        Select Case (TRIM(ADJUSTL(case_name)))
        Case ('ShockTube')
            Call set_shock_tube_parameters(st, 1.4_PR)! The parameter 1.4 doesn't matter here
            xmin = 0._PR
            xmax = st%length
            ymin = -.5_PR*st%width
            ymax = .5_PR*st%width
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

End Module mod_cases
