Module mod_limiters

    Use mod_parameters
    Implicit None

Contains

    Subroutine reconstructAtInterface(axis, ULi, URi, ULL, UL, UR, URR)
        ! --- InOut
        Character, Intent(In) :: axis
        Real(PR), Dimension(4), Intent(In) :: ULL, UL, UR, URR
        Real(PR), Dimension(4), Intent(InOut) :: ULi, URi ! approximations at the interface
        ! --- Parameters
        Real(PR), Parameter :: eps = 1.e0_PR*EPSILON(Real(PR)) !to avoid division by 0. when computing the smoothness
        ! --- Locals
        Real(PR), Dimension(4) :: dUL, dUR 
        Real(PR), Dimension(4) :: smoothnessL, limiterL
        Real(PR), Dimension(4) :: smoothnessR, limiterR
        Real(PR) :: momentum

        Select Case (num_scheme%space_scheme_order)
        Case (1) ! No reconstruction to do for the first order
            ULi = UL
            URi = UR
        Case (2) ! Second order TVD-MUSCL : linear reconstruction
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

            !Select Case (space_scheme_specs%slope_str)
            !Case ('left')
            !    dUL = .5_PR * ( UL - ULL )
            !    dUR = .5_PR * ( UR - UL )
            !Case ('right')
            !    dUL = .5_PR * ( UR - UL )
            !    dUR = .5_PR * ( URR - UR )
            !Case ('upwind')
            !    If ( momentum >= 0._PR ) Then
            !        dUL = .5_PR * ( UL - ULL )
            !        dUR = .5_PR * ( UR - UL )
            !    Else
            !        dUL = .5_PR * ( UR - UL )
            !        dUR = .5_PR * ( URR - UR )
            !    End If
            !Case ('downwind')
            !    If ( momentum >= 0._PR ) Then
            !        dUL = .5_PR * ( UR - UL )
            !        dUR = .5_PR * ( URR - UR )
            !    Else
            !        dUL = .5_PR * ( UL - ULL )
            !        dUR = .5_PR * ( UR - UL )
            !    End If
            !Case ('centered')
            !    dUL = .25_PR * ( UR - ULL )
            !    dUR = .25_PR * ( URR - UL )
            !Case Default
            !    Write(STDERR,*) "Unknown slope ", space_scheme_specs%slope_str
            !    Call Exit(1)
            !End Select

            ! Centered slope
            dUL = .25_PR * ( UR - ULL )
            dUR = .25_PR * ( URR - UL )

            Select Case (num_scheme%MUSCL_limiter)
            Case ('None')
                limiterL = 1._PR
                limiterR = 1._PR
            Case ('minmod')
                limiterL = minmod(smoothnessL)
                limiterR = minmod(smoothnessR)
            Case ('vanLeer')
                limiterL = vanLeer(smoothnessL)
                limiterR = vanLeer(smoothnessR)
            Case ('superbee')
                limiterL = SuperBee(smoothnessL)
                limiterR = SuperBee(smoothnessR)
            Case ('UMIST')
                limiterL = UMIST(smoothnessL)
                limiterR = UMIST(smoothnessR)
            Case ('genminmod')
                limiterL = generalisedMinmod(smoothnessL, num_scheme%MUSCL_generalised_minmod_parameter)
                limiterR = generalisedMinmod(smoothnessR, num_scheme%MUSCL_generalised_minmod_parameter)
            Case ('order3')
                limiterL = limiter_order3(smoothnessL)
                limiterR = limiter_order3(1._PR/smoothnessR)
            Case Default
                Write(STDERR,*) "Unknown limiter ", num_scheme%MUSCL_limiter
                Call Exit(1)
            End Select

            ULi = UL + limiterL*dUL
            URi = UR - limiterR*dUR

        Case Default
            Write(STDERR, *) "Order greater than 2 not available (Order", num_scheme%space_scheme_order, "specified)"
            Call Exit(1)
        End Select
    End Subroutine reconstructAtInterface

    Function minmod(smoothness)
        ! --- InOut
        Real(PR), Dimension(4), Intent(In) :: smoothness
        Real(PR), Dimension(4) :: minmod

         minmod = MAX( 0._PR, MIN(1._PR, smoothness) )
    End Function minmod

    Function vanLeer(smoothness)
        ! --- InOut
        Real(PR), Dimension(4), Intent(In) :: smoothness
        Real(PR), Dimension(4) :: vanLeer

        vanLeer = ( smoothness + ABS(smoothness) )/( 1._PR + ABS(smoothness) )
    End Function vanLeer

    Function SuperBee(smoothness)
        ! --- InOut
        Real(PR), Dimension(4), Intent(In) :: smoothness
        Real(PR), Dimension(4) :: SuperBee

        SuperBee = MAX( 0._PR, MIN(1._PR, 2._PR*smoothness), MIN(2._PR, smoothness) )
    End Function SuperBee

    Function UMIST(smoothness)
        ! --- InOut
        Real(PR), Dimension(4), Intent(In) :: smoothness
        Real(PR), Dimension(4) :: UMIST

        UMIST = MAX( 0._PR, MIN(2._PR*smoothness, .25_PR*smoothness + .75_PR, .75_PR*smoothness + .25_PR, 2._PR ) )
    End Function UMIST

    Function generalisedMinmod(smoothness, p)
        ! --- InOut
        Real(PR), Dimension(4), Intent(In) :: smoothness
        Real(PR), Intent(In) :: p
        Real(PR), Dimension(4) :: generalisedMinmod

        generalisedMinmod = MAX( 0._PR, &
            & MIN( p * smoothness, .5_PR * ( smoothness + 1._PR ), &
            & p ) &
            & )
    End Function generalisedMinmod

    Function limiter_order3(smoothness)
        ! --- InOut
        Real(PR), Dimension(4), Intent(In) :: smoothness
        Real(PR), Dimension(4) :: limiter_order3

        limiter_order3 = MAX(0._PR, MIN((2._PR+smoothness)/3._PR, &
            & MAX(-.5_PR*smoothness, MIN(2._PR*smoothness, (2._PR+smoothness)/3._PR, 1.6_PR)) &
            & ))
    End Function limiter_order3

End Module mod_limiters
