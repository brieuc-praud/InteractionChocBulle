Module mod_fluxes

    Use mod_parameters
    Use mod_limiters
    Implicit None

Contains

    ! ===== Generical flux interface =====
    Function numericalFlux(axis, ULL, UL, UR, URR)
        ! --- InOut
        Character, Intent(In) :: axis ! x,  y
        Real(PR), Dimension(4), Intent(In) :: ULL, UL, UR, URR
        Real(PR), Dimension(4) :: numericalFlux
        ! --- Locals
        Real(PR), Dimension(4) :: ULi, URi

        Call reconstructAtInterface(axis, ULi, URi, ULL, UL, UR, URR)

        Select Case (TRIM(ADJUSTL(num_scheme%space_scheme_name)))
        Case ('Rusanov')
            numericalFlux = Rusanov(axis, ULi, URi)
        Case ('HLL')
            numericalFlux = HLL(axis, ULi, URi)
        Case ('HLLC')
            numericalFlux = HLLC(axis, ULi, URi)
        Case Default
            Write(STDERR,*) "Unknown flux name ", TRIM(ADJUSTL(num_scheme%space_scheme_name))
            Call Exit(1)
        End Select
    End Function numericalFlux

    ! ===== Numerical fluxes implementations =====
    Function Rusanov(axis, UL, UR)
        ! --- InOut
        Real(PR), Dimension(4), Intent(In) :: UL, UR
        Character, Intent(In) :: axis ! x,  y
        Real(PR), Dimension(4) :: Rusanov
        ! --- Locals
        Real(PR) :: rL, velocity_uL, velocity_vL, velocityL, pressureL, aL, b_minusL, b_plusL, bL
        Real(PR) :: rR, velocity_uR, velocity_vR, velocityR, pressureR, aR, b_minusR, b_plusR, bR
        Real(PR) :: b

        ! == Left
        Call conservativeToPrimitive(UL, rL, velocity_uL, velocity_vL, pressureL)
        aL = SQRT(gammagp * pressureL / rL)
        ! == Right
        Call conservativeToPrimitive(UR, rR, velocity_uR, velocity_vR, pressureR)
        aR = SQRT(gammagp * pressureR / rR)

        Select Case (axis)
        Case ('x')
            velocityL = velocity_uL
            velocityR = velocity_uR
        Case ('y')
            velocityL = velocity_vL
            velocityR = velocity_vR
        Case Default
            Write(STDERR, *) "Unknown axis ", axis
            Call Exit(1)
        End Select

        b_minusL = velocityL - aL
        b_plusL = velocityL + aL
        b_minusR = velocityR - aR
        b_plusR = velocityR + aR

        bL = MAX( ABS(b_minusL), ABS(b_plusL) )
        bR = MAX( ABS(b_minusR), ABS(b_plusR) )
        b = MAX( bL, bR )

        Rusanov = fluxFunc(axis, UL) + fluxFunc(axis, UR)
        Rusanov = .5_PR * ( Rusanov - b*(UR - UL) )
    End Function Rusanov

    Function HLL(axis, UL, UR)
        ! --- InOut
        Real(PR), Dimension(4), Intent(In) :: UL, UR
        Character, Intent(In) :: axis ! x,  y
        Real(PR), Dimension(4) :: HLL
        ! --- Locals
        Real(PR) :: rL, velocity_uL, velocity_vL, velocityL, pressureL, aL, b_minusL, b_plusL
        Real(PR), Dimension(4) :: fluxL
        Real(PR) :: rR, velocity_uR, velocity_vR, velocityR, pressureR, aR, b_minusR, b_plusR
        Real(PR), Dimension(4) :: fluxR
        Real(PR) :: b_minus, b_plus

        ! == Left
        Call conservativeToPrimitive(UL, rL, velocity_uL, velocity_vL, pressureL)
        aL = SQRT(gammagp * pressureL / rL)
        ! == Right
        Call conservativeToPrimitive(UR, rR, velocity_uR, velocity_vR, pressureR)
        aR = SQRT(gammagp * pressureR / rR)

        Select Case (axis)
        Case ('x')
            velocityL = velocity_uL
            velocityR = velocity_uR
        Case ('y')
            velocityL = velocity_vL
            velocityR = velocity_vR
        Case Default
            Write(STDERR, *) "Unknown axis ", axis
            Call Exit(1)
        End Select

        fluxL = fluxFunc(axis, UL)
        fluxR = fluxFunc(axis, UR)

        b_minusL = velocityL - aL
        b_plusL = velocityL + aL
        b_minusR = velocityR - aR
        b_plusR = velocityR + aR
        b_minus = MIN( b_minusL, b_minusR )
        b_plus = MAX( b_plusL, b_plusR )

        If ( b_minus >= 0._PR ) Then
            HLL = fluxL
        Else If ( b_plus >= 0._PR ) Then
            HLL = ( b_plus * fluxL - b_minus * fluxR &
                &   +   b_plus * b_minus * (UR - UL) ) / ( b_plus - b_minus )
        Else
            HLL = fluxR
        End If
    End Function HLL

    Function HLLC(axis, UL, UR)
        ! --- InOut
        Real(PR), Dimension(4), Intent(In) :: UL, UR
        Character, Intent(In) :: axis ! x,  y
        Real(PR), Dimension(4) :: HLLC
        ! --- Locals
        Real(PR) :: rL, velocity_uL, velocity_vL, velocityL, TvelocityL, pressureL, aL, b_minusL, b_plusL
        Real(PR), Dimension(4) :: fluxL, fluxL_star, UL_star
        Real(PR) :: rR, velocity_uR, velocity_vR, velocityR, TvelocityR, pressureR, aR, b_minusR, b_plusR
        Real(PR), Dimension(4) :: fluxR, fluxR_star, UR_star
        Real(PR) :: b_minus, b_star, b_plus

        ! == Left
        Call conservativeToPrimitive(UL, rL, velocity_uL, velocity_vL, pressureL)
        aL = SQRT(gammagp * pressureL / rL)
        ! == Right
        Call conservativeToPrimitive(UR, rR, velocity_uR, velocity_vR, pressureR)
        aR = SQRT(gammagp * pressureR / rR)

        Select Case (axis)
        Case ('x')
            velocityL = velocity_uL
            velocityR = velocity_uR
            TvelocityL = velocity_vL
            TvelocityR = velocity_vR
        Case ('y')
            velocityL = velocity_vL
            velocityR = velocity_vR
            TvelocityL = velocity_uL
            TvelocityR = velocity_uR
        Case Default
            Write(STDERR, *) "Unknown axis ", axis
            Call Exit(1)
        End Select

        fluxL = fluxFunc(axis, UL)
        fluxR = fluxFunc(axis, UR)

        b_minusL = velocityL - aL
        b_plusL = velocityL + aL
        b_minusR = velocityR - aR
        b_plusR = velocityR + aR
        b_minus = MIN( b_minusL, b_minusR )
        b_plus = MAX( b_plusL, b_plusR )

        b_star = pressureR - pressureL &
            & + rL*velocityL*(b_minus - velocityL) &
            & - rR*velocityR*(b_plus - velocityR)
        b_star = b_star / ( rL*(b_minus - velocityL) &
            & - rR*(b_plus - velocityR) )

        ! == * Left
        UL_star(1) = 1._PR
        Select Case (axis)
        Case ('x')
            UL_star(2) = b_star
            UL_star(3) = TvelocityL
        Case ('y')
            UL_star(2) = TvelocityL
            UL_star(3) = b_star
        Case Default
            Write(STDERR, *) "Unknown axis ", axis
            Call Exit(1)
        End Select
        UL_star(4) = UL(4)/rL + ( b_star - velocityL ) &
            & * ( b_star + pressureL / ( rL*(b_minus - velocityL) ) )
        UL_star = UL_star * rL &
            & * ( b_minus - velocityL )/( b_minus - b_star )
        ! == * Right
        UR_star(1) = 1._PR
        Select Case (axis)
        Case ('x')
            UR_star(2) = b_star
            UR_star(3) = TvelocityR
        Case ('y')
            UR_star(2) = TvelocityR
            UR_star(3) = b_star
        Case Default
            Write(STDERR, *) "Unknown axis ", axis
            Call Exit(1)
        End Select
        UR_star(4) = UR(4)/rR + ( b_star - velocityR ) &
            & * ( b_star + pressureR / ( rR*(b_plus - velocityR) ) )
        UR_star = UR_star * rR &
            & * ( b_plus - velocityR )/( b_plus - b_star )

        If ( b_minus >= 0._PR ) Then
            HLLC = fluxL
        Else If ( b_star >= 0._PR ) Then
            fluxL_star = fluxL + b_minus * ( UL_star - UL )
            HLLC = fluxL_star
        Else If ( b_plus >= 0._PR ) Then
            fluxR_star = fluxR + b_plus * ( UR_star - UR )
            HLLC = fluxR_star
        Else
            HLLC = fluxR
        End If

    End Function HLLC


    ! ===== Flux function =====
    Function fluxFunc(axis, Uvect)
        ! --- InOut
        Character, Intent(In) :: axis
        Real(PR), Dimension(4), Intent(In) :: Uvect
        Real(PR), Dimension(4) :: fluxFunc
        ! --- Locals
        Real(PR) :: r, ru, rv, e, u, v, p

        ru = Uvect(2)
        rv = Uvect(3)
        e = Uvect(4)
        Call conservativeToPrimitive(Uvect, r, u, v, p)

        Select Case (axis)
        Case ('x')
            fluxFunc(1) = ru
            fluxFunc(2) = ru*u + p
            fluxFunc(3) = rv*u
            fluxFunc(4) = (e + p)*u
        Case ('y')
            fluxFunc(1) = rv
            fluxFunc(2) = ru*v
            fluxFunc(3) = rv*v + p
            fluxFunc(4) = (e + p)*v
        Case Default
            Write(STDERR, *) "Unknown axis ", axis
            Call Exit(1)
        End Select
    End Function fluxFunc

End Module mod_fluxes
