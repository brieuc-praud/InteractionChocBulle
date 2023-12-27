Module mod_test

    Use mod_parameters
    Use mod_schemes
    Implicit None

Contains
    Subroutine PerformTestsAndExit(gammagp)
        ! ------------------- Intent In -----------------------
        Real(PR), Intent(In) :: gammagp
        ! ----------------------------------------------------------
        ! === Check fluxes consistency ===
        Write(STDOUT,*) '=== FLUXES CONSISTENCY ==='

        Write(STDOUT,*) "-- Rusanov"
        Call testNumericalFlux(Rusanov, gammagp, fluxFuncF, fluxFuncG)
        Write(STDOUT,*) " -> SUCCESS" ! The program exits if a test fails

        Write(STDOUT,*) "-- HLL"
        Call testNumericalFlux(HLL, gammagp, fluxFuncF, fluxFuncG)
        Write(STDOUT,*) " -> SUCCESS"

        Write(STDOUT,*) "-- HLLC"
        Call testNumericalFlux(HLLC, gammagp, fluxFuncF, fluxFuncG)
        Write(STDOUT,*) " -> SUCCESS"

        Call Exit()
    End Subroutine PerformTestsAndExit

    Function randomU(gammagp)
        Real(PR), Parameter :: eps = 1.e-3_PR
        Real(PR), Parameter :: hug = 1.e3_PR
        Real(PR), Dimension(4), Parameter :: lowerBounds = &
            & (/ Real(PR) :: eps, -hug, -hug, 0._PR /)
        Real(PR), Dimension(4), Parameter :: upperBounds = &
            & (/ Real(PR) :: hug, hug, hug, hug /)

        Real(PR), Intent(In) :: gammagp
        Real(PR), Dimension(4) :: randomU

        Real(PR), Dimension(4) :: randomUniform

        Real(PR) :: r, u, v, q, p, e

        Call RANDOM_NUMBER(randomUniform)
        randomU = upperBounds * randomUniform + lowerBounds * ( 1._PR - randomUniform )
        ! The energy is recovered from the pressure
        r = randomU(1)
        u = randomU(2)/r
        v = randomU(3)/r
        p = randomU(4)
        q = .5_PR * ( u**2 + v**2 )
        e = p / (gammagp - 1._PR) + r*q
        randomU(4) = e
    End Function randomU

    Subroutine testNumericalFlux(numericalFlux, gammagp, fluxFuncF, fluxFuncG)
        Real(PR), Parameter :: safety_factor = 1.e6_PR
        Integer, Parameter :: number_of_tests = 1000
        ! ------------------- Intent In -----------------------
        Real(PR), Intent(In) :: gammagp
        Interface
            Function numericalFlux(axis, UL, UR, gammagp)
                Import PR
                Real(PR), Dimension(4), Intent(In) :: UL, UR
                Real(PR), Intent(In) :: gammagp
                Character, Intent(In) :: axis
                Real(PR), Dimension(4) :: numericalFlux
            End Function numericalFlux
        End Interface
        Interface
            Function fluxFuncF(Uvect, gammagp)
                Import PR
                Real(PR), Dimension(4), Intent(In) :: Uvect
                Real(PR), Intent(In) :: gammagp
                Real(PR), Dimension(4) :: fluxFuncF
            End Function fluxFuncF
        End Interface
        Interface
            Function fluxFuncG(Uvect, gammagp)
                Import PR
                Real(PR), Dimension(4), Intent(In) :: Uvect
                Real(PR), Intent(In) :: gammagp
                Real(PR), Dimension(4) :: fluxFuncG
            End Function fluxFuncG
        End Interface
        ! ----------------------------------------------------------
        Integer :: i
        Real(PR), Dimension(4) :: U
        Real(PR) :: diffF, diffG
        Logical :: failF, failG

        Do i=1, number_of_tests
            U = randomU(gammagp)
            diffF = MAXVAL( ABS(numericalFlux('x', U, U, gammagp) - fluxFuncF(U, gammagp)) )
            diffG = MAXVAL( ABS(numericalFlux('y', U, U, gammagp) - fluxFuncG(U, gammagp)) )
            failF = diffF > safety_factor*EPSILON(U(1))
            failG = diffG > safety_factor*EPSILON(U(1))
            If (failF .OR. failG) Then
                Write(STDOUT, *) diffF, diffG
                Call Exit(1)
            End If
        End Do
    End Subroutine testNumericalFlux

End Module mod_test
