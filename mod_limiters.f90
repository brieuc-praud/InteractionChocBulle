Module mod_limiters

    Use mod_parameters
    Implicit None

Contains

    Subroutine reconstructAtInterface(axis, ULi, URi, ULLc, ULc, URc, URRc)
        ! --- InOut
        Character, Intent(In) :: axis
        Real(PR), Dimension(4), Intent(In) :: ULLc, ULc, URc, URRc
        Real(PR), Dimension(4), Intent(InOut) :: ULi, URi ! approximations at the interface
        ! --- Parameters
        Real(PR), Parameter :: eps = 1.e0_PR*EPSILON(Real(PR)) !to avoid division by 0. when computing the smoothness
        Real(PR), Parameter :: radiusO3 = 1.e-1
        ! --- Locals
        Real(PR), Dimension(4) :: dUL, dUR 
        Real(PR), Dimension(4) :: smoothnessL, limiterL
        Real(PR), Dimension(4) :: smoothnessR, limiterR
        Real(PR), Dimension(4) :: ULL, UL, UR, URR, ULip, URip
        Real(PR), Dimension(4) :: etaL, etaR
        Real(PR) :: dd

        Call conservativeToPrimitive(ULLc, ULL(1), ULL(2), ULL(3), ULL(4))
        Call conservativeToPrimitive(ULc, UL(1), UL(2), UL(3), UL(4))
        Call conservativeToPrimitive(URc, UR(1), UR(2), UR(3), UR(4))
        Call conservativeToPrimitive(URRc, URR(1), URR(2), URR(3), URR(4))

        ! local smoothness measure
        smoothnessL = ( UL - ULL + SIGN(eps, UL-ULL)) / ( UR - UL + SIGN(eps, UR-UL))
        smoothnessR = ( UR - UL + SIGN(eps, UR-UL)) / ( URR - UR + SIGN(eps, URR-UR))

        ! Centered slope
        dUL = .25_PR * ( UR - ULL )
        dUR = .25_PR * ( URR - UL )

        Select Case (num_scheme%space_scheme_order)
        Case (1) ! No reconstruction to do for the first order
            ULip = UL
            URip = UR
        Case (2) ! Second order TVD-MUSCL : linear reconstruction
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
            Case Default
                Write(STDERR,*) "Unknown limiter ", num_scheme%MUSCL_limiter
                Call Exit(1)
            End Select

            ULip = UL + limiterL*dUL
            URip = UR - limiterR*dUR
        Case (3) ! LimO3
            Select Case (axis)
            Case ('x')
                dd = deltax
            Case ('y')
                dd = deltay
            Case Default
                Write(STDERR,*) "Unknown limiter ", num_scheme%MUSCL_limiter
                Call Exit(1)
            End Select

            etaL = ( (UR - UL)**2 + (UL - ULL)**2 )/radiusO3/dd**2
            etaR = ( (URR - UR)**2 + (UR - UL)**2 )/radiusO3/dd**2

            limiterL = limiter_O3(etaR, smoothnessL)
            limiterR = limiter_O3(etaL, 1._PR/smoothnessR)

            ULip = UL + .5_PR*limiterL*(UR-UL)
            URip = UR - .5_PR*limiterR*(UR-UL)
        Case Default
            Write(STDERR, *) "Order greater than 2 not available (Order", num_scheme%space_scheme_order, "specified)"
            Call Exit(1)
        End Select

        Call primitiveToConservative(ULip, ULi(1), ULi(2), ULi(3), ULi(4))
        Call primitiveToConservative(URip, URi(1), URi(2), URi(3), URi(4))

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

    Function limiter_O3(eta, smoothness)
        ! --- InOut
        Real(PR), Dimension(4), Intent(In) :: smoothness
        Real(PR), Dimension(4), Intent(In) :: eta
        Real(PR), Dimension(4) :: limiter_O3
        ! --- Parameters
        Real(PR), Parameter :: eps = 1e3*EPSILON(Real(PR))
        Integer :: i

        limiter_O3 = MAX(0._PR, MIN((2._PR+smoothness)/3._PR, &
            & MAX(-.5_PR*smoothness, MIN(2._PR*smoothness, (2._PR+smoothness)/3._PR, 1.6_PR)) &
            & ))

        Do i=1,4
            If (eta(i) <= 1._PR - eps) Then
                limiter_O3(i) = (2._PR + smoothness(i))/3._PR
            Else If (eta(i) < 1._PR + eps) Then
                limiter_O3(i) = (1._PR + (eta(i)-1)/eps)*limiter_O3(i)
                limiter_O3(i) = .5_PR * ( (1._PR - (eta(i)-1)/eps)*(2._PR+smoothness(i))/3._PR + limiter_O3(i) )
            End If
        End Do
    End Function limiter_O3

End Module mod_limiters
