Module mod_limiters

    Use mod_parameters
    Implicit None

    Type :: space_scheme
        Integer :: order
        Real(PR) :: generalised_minmod_parameter
        Character(len=20) :: slope_str
        Character(len=20) :: limiter_str
    End Type space_scheme

Contains

    Subroutine reconstructAtInterface(axis, ULi, URi, space_scheme_specs, ULL, UL, UR, URR)
        ! --- InOut
        Character, Intent(In) :: axis
        Type(space_scheme), Intent(In) :: space_scheme_specs
        Real(PR), Dimension(4), Intent(In) :: ULL, UL, UR, URR
        Real(PR), Dimension(4), Intent(InOut) :: ULi, URi ! approximations at the interface
        ! --- Parameters
        Real(PR), Parameter :: eps = 1.e0_PR*EPSILON(Real(PR)) !to avoid division by 0. when computing the smoothness
        ! --- Locals
        Real(PR), Dimension(4) :: dUL, dUR 
        Real(PR), Dimension(4) :: smoothnessL, limiterL
        Real(PR), Dimension(4) :: smoothnessR, limiterR
        Real(PR) :: momentum

        Select Case (space_scheme_specs%order)
        Case (1)
            ULi = UL
            URi = UR
        Case (2) ! second order TVD-MUSCL
            ! local smoothness measure
            smoothnessL = ( UL - ULL + SIGN(eps, UL-ULL)) / ( UR - UL + SIGN(eps, UR-UL))
            smoothnessR = ( UR - UL + SIGN(eps, UR-UL)) / ( URR - UR + SIGN(eps, URR-UR))

            Select Case (axis)
            Case ('x')
                momentum = UL(2)
            Case ('y')
                momentum = UL(3)
            Case Default
                Write(STDERR, *) "Unknown axis ", axis
                Call Exit(1)
            End Select

            Select Case (space_scheme_specs%slope_str)
            Case ('left')
                dUL = .5_PR * ( UL - ULL )
                dUR = .5_PR * ( UR - UL )
            Case ('right')
                dUL = .5_PR * ( UR - UL )
                dUR = .5_PR * ( URR - UR )
            Case ('upwind')
                If ( momentum >= 0._PR ) Then
                    dUL = .5_PR * ( UL - ULL )
                    dUR = .5_PR * ( UR - UL )
                Else
                    dUL = .5_PR * ( UR - UL )
                    dUR = .5_PR * ( URR - UR )
                End If
            Case ('downwind')
                If ( momentum >= 0._PR ) Then
                    dUL = .5_PR * ( UR - UL )
                    dUR = .5_PR * ( URR - UR )
                Else
                    dUL = .5_PR * ( UL - ULL )
                    dUR = .5_PR * ( UR - UL )
                End If
            Case ('centered')
                dUL = .25_PR * ( UR - ULL )
                dUR = .25_PR * ( URR - UL )
            Case Default
                Write(STDERR,*) "Unknown slope ", space_scheme_specs%slope_str
                Call Exit(1)
            End Select

            Select Case (space_scheme_specs%limiter_str)
            Case ('None')
                limiterL = 1._PR
                limiterR = 1._PR
            Case ('minmod')
                limiterL = MAX( 0._PR, MIN(1._PR, smoothnessL) )
                limiterR = MAX( 0._PR, MIN(1._PR, smoothnessR) )
            Case ('vanLeer')
                limiterL = ( smoothnessL + ABS(smoothnessL) )/( 1._PR + ABS(smoothnessL) )
                limiterR = ( smoothnessR + ABS(smoothnessR) )/( 1._PR + ABS(smoothnessR) )
            Case ('superbee')
                limiterL = MAX( 0._PR, MIN(1._PR, 2._PR*smoothnessL), MIN(2._PR, smoothnessL) )
                limiterR = MAX( 0._PR, MIN(1._PR, 2._PR*smoothnessR), MIN(2._PR, smoothnessR) )
            Case ('UMIST')
                limiterL = MAX( 0._PR, MIN(2._PR*smoothnessL, .25_PR*smoothnessL + .75_PR, &
                    & .75_PR*smoothnessL + .25_PR, 2._PR ) )
                limiterR = MAX( 0._PR, MIN(2._PR*smoothnessR, .25_PR*smoothnessR + .75_PR, &
                    & .75_PR*smoothnessR + .25_PR, 2._PR ) )
            Case ('genminmod') ! Generalised minmod limiter
                limiterL = MAX( 0._PR, &
                    & MIN( space_scheme_specs%generalised_minmod_parameter * smoothnessL, &
                    & .5_PR * ( smoothnessL + 1._PR ), &
                    & space_scheme_specs%generalised_minmod_parameter ) &
                    & )
                limiterR = MAX( 0._PR, &
                    & MIN( space_scheme_specs%generalised_minmod_parameter * smoothnessR, &
                    & .5_PR * ( smoothnessR + 1._PR ), &
                    & space_scheme_specs%generalised_minmod_parameter ) &
                    & )
            Case Default
                Write(STDERR,*) "Unknown limiter ", space_scheme_specs%limiter_str
                Call Exit(1)
            End Select

            ULi = UL + limiterL*dUL
            URi = UR - limiterR*dUR

        Case Default
            Write(STDERR, *) "Order greater than 2 not available (Order", space_scheme_specs%order, "specified)"
            Call Exit(1)
        End Select
    End Subroutine reconstructAtInterface

End Module mod_limiters
