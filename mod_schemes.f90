Module mod_schemes

   Use mod_parameters
   Implicit None

Contains
    Function Rusanov(axis, UL, UR, gammagp)
        ! --- InOut
        Real(PR), Dimension(4), Intent(In) :: UL, UR
        Real(PR), Intent(In) :: gammagp
        Character, Intent(In) :: axis ! x,  y
        Real(PR), Dimension(4) :: Rusanov
        ! --- Locals
        Real(PR) :: rL, velocity_uL, velocity_vL, eL, pressureL, qL, aL, bL, l1L, l3L
        Real(PR) :: rR, velocity_uR, velocity_vR, eR, pressureR, qR, aR, bR, l1R, l3R
        Real(PR) :: b

        ! == Left
        rL = UL(1)
        velocity_uL = UL(2)/rL
        velocity_vL = UL(3)/rL
        eL = UL(4)
        qL = .5_PR * ( velocity_uL**2 + velocity_vL**2 )
        pressureL = (gammagp - 1._PR)*(eL - rL*qL)
        aL = SQRT(gammagp * pressureL / rL)

        ! == Right
        rR = UR(1)
        velocity_uR = UR(2)/rR
        velocity_vR = UR(3)/rR
        eR = UR(4)
        qR = .5_PR * ( velocity_uR**2 + velocity_vR**2 )
        pressureR = (gammagp - 1._PR)*(eR - rR*qR)
        aR = SQRT(gammagp * pressureR / rR)

        Select Case (axis)
        Case ('y')
            l1L = ABS(velocity_vL - aL)
            l3L = ABS(velocity_vL + aL)
            l1R = ABS(velocity_vR - aR)
            l3R = ABS(velocity_vR + aR)
        Case Default ! 'x'
            l1L = ABS(velocity_uL - aL)
            l3L = ABS(velocity_uL + aL)
            l1R = ABS(velocity_uR - aR)
            l3R = ABS(velocity_uR + aR)
        End Select

        bL = MAX(l1L, l3L)
        bR = MAX(l1R, l3R)
        b = MAX( bL, bR )

        Select Case (axis)
        Case ('y')
            Rusanov = fluxFuncG(UL, gammagp) + fluxFuncG(UR, gammagp)
        Case Default ! 'x'
            Rusanov = fluxFuncF(UL, gammagp) + fluxFuncF(UR, gammagp)
        End Select
        Rusanov = 0.5_PR * ( Rusanov - b*(UR - UL) )
    End Function Rusanov

    Function HLL(axis, UL, UR, gammagp)
        ! --- InOut
        Real(PR), Dimension(4), Intent(In) :: UL, UR
        Real(PR), Intent(In) :: gammagp
        Character, Intent(In) :: axis ! x,  y
        Real(PR), Dimension(4) :: HLL
        ! --- Locals
        Real(PR) :: rL, velocity_uL, velocity_vL, velocityL, eL, pressureL, qL, aL, b_minusL, b_plusL
        Real(PR), Dimension(4) :: fluxL
        Real(PR) :: rR, velocity_uR, velocity_vR, velocityR, eR, pressureR, qR, aR, b_minusR, b_plusR
        Real(PR), Dimension(4) :: fluxR
        Real(PR) :: b_minus, b_plus

        ! == Left
        rL = UL(1)
        velocity_uL = UL(2)/rL
        velocity_vL = UL(3)/rL
        eL = UL(4)
        qL = .5_PR * ( velocity_uL**2 + velocity_vL**2 )
        pressureL = (gammagp - 1._PR)*(eL - rL*qL)
        aL = SQRT(gammagp * pressureL / rL)

        ! == Right
        rR = UR(1)
        velocity_uR = UR(2)/rR
        velocity_vR = UR(3)/rR
        eR = UR(4)
        qR = .5_PR * ( velocity_uR**2 + velocity_vR**2 )
        pressureR = (gammagp - 1._PR)*(eR - rR*qR)
        aR = SQRT(gammagp * pressureR / rR)

        Select Case (axis)
        Case ('y')
            velocityL = velocity_vL
            velocityR = velocity_vR
            fluxL = fluxFuncG(UL, gammagp)
            fluxR = fluxFuncG(UR, gammagp)
        Case Default ! 'x'
            velocityL = velocity_uL
            velocityR = velocity_uR
            fluxL = fluxFuncF(UL, gammagp)
            fluxR = fluxFuncF(UR, gammagp)
        End Select

        b_minusL = velocityL - aL
        b_plusL = velocityL + aL
        b_minusR = velocityR - aR
        b_plusR = velocityR + aR
        b_minus = MIN( b_minusL, b_minusR, 0._PR )
        b_plus = MAX( b_plusL, b_plusR, 0._PR )
        
        HLL = b_plus * fluxL - b_minus * fluxR
        HLL = ( HLL   +   b_plus * b_minus * (UR - UL) ) / ( b_plus - b_minus )
    End Function HLL

    Function HLLC(axis, UL, UR, gammagp)
        ! --- InOut
        Real(PR), Dimension(4), Intent(In) :: UL, UR
        Real(PR), Intent(In) :: gammagp
        Character, Intent(In) :: axis ! x,  y
        Real(PR), Dimension(4) :: HLLC
        ! --- Locals
        Real(PR) :: rL, velocity_uL, velocity_vL, velocityL, TvelocityL, eL, pressureL, qL, aL, b_minusL, b_plusL
        Real(PR), Dimension(4) :: fluxL, fluxL_star, UL_star
        Real(PR) :: rR, velocity_uR, velocity_vR, velocityR, TvelocityR, eR, pressureR, qR, aR, b_minusR, b_plusR
        Real(PR), Dimension(4) :: fluxR, fluxR_star, UR_star
        Real(PR) :: b_minus, b_star, b_plus

        ! == Left
        rL = UL(1)
        velocity_uL = UL(2)/rL
        velocity_vL = UL(3)/rL
        eL = UL(4)
        qL = .5_PR * ( velocity_uL**2 + velocity_vL**2 )
        pressureL = (gammagp - 1._PR)*(eL - rL*qL)
        aL = SQRT(gammagp * pressureL / rL)

        ! == Right
        rR = UR(1)
        velocity_uR = UR(2)/rR
        velocity_vR = UR(3)/rR
        eR = UR(4)
        qR = .5_PR * ( velocity_uR**2 + velocity_vR**2 )
        pressureR = (gammagp - 1._PR)*(eR - rR*qR)
        aR = SQRT(gammagp * pressureR / rR)

        Select Case (axis)
        Case ('y')
            velocityL = velocity_vL
            velocityR = velocity_vR
            TvelocityL = velocity_uL
            TvelocityR = velocity_uR
            fluxL = fluxFuncG(UL, gammagp)
            fluxR = fluxFuncG(UR, gammagp)
        Case Default ! 'x'
            velocityL = velocity_uL
            velocityR = velocity_uR
            TvelocityL = velocity_vL
            TvelocityR = velocity_vR
            fluxL = fluxFuncF(UL, gammagp)
            fluxR = fluxFuncF(UR, gammagp)
        End Select

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
        Case ('y')
            UL_star(2) = TvelocityL
            UL_star(3) = b_star
        Case Default ! 'x'
            UL_star(2) = b_star
            UL_star(3) = TvelocityL
        End Select
        UL_star(4) = eL/rL + ( b_star - velocityL ) &
            & * ( b_star + pressureL / ( rL*(b_minus - velocityL) ) )
        UL_star = UL_star * rL &
            & * ( b_minus - velocityL )/( b_minus - b_star )
        ! == * Right
        UR_star(1) = 1._PR
        Select Case (axis)
        Case ('y')
            UR_star(3) = b_star
            UR_star(2) = TvelocityR
        Case Default ! 'x'
            UR_star(2) = b_star
            UR_star(3) = TvelocityR
        End Select
        UR_star(4) = eR/rR + ( b_star - velocityR ) &
            & * ( b_star + pressureR / ( rR*(b_plus - velocityR) ) )
        UR_star = UR_star * rR &
            & * ( b_plus - velocityR )/( b_plus - b_star )

        fluxL_star = fluxL + b_minus * ( UL_star - UL )
        fluxR_star = fluxR + b_plus * ( UR_star - UR )

        If ( b_minus >= 0._PR ) Then
            HLLC = fluxL
        Else If ( b_star >= 0._PR ) Then
            HLLC = fluxL_star
        Else If ( b_plus >= 0._PR ) Then
            HLLC = fluxR_star
        Else
            HLLC = fluxR
        End If

    End Function HLLC


    ! Fluxes functions
    Function fluxFuncF(Uvect, gammagp)
        Real(PR), Dimension(4), Intent(In) :: Uvect
        Real(PR), Intent(In) :: gammagp
        Real(PR), Dimension(4) :: fluxFuncF
        ! --- Locals
        Real(PR) :: r, ru, rv, e, u, v, p, q

        r = Uvect(1)
        ru = Uvect(2)
        rv = Uvect(3)
        e = Uvect(4)
        u = ru/r
        v = rv/r
        q = 0.5_PR * ( u**2 + v**2 )
        p = (gammagp - 1._PR)*(e - r*q)

        fluxFuncF(1) = ru
        fluxFuncF(2) = ru*u + p
        fluxFuncF(3) = rv*u
        fluxFuncF(4) = (e + p)*u
    End Function fluxFuncF

    Function fluxFuncG(Uvect, gammagp)
        Real(PR), Dimension(4), Intent(In) :: Uvect
        Real(PR), Intent(In) :: gammagp
        Real(PR), Dimension(4) :: fluxFuncG
        ! --- Locals
        Real(PR) :: r, ru, rv, e, u, v, p, q

        r = Uvect(1)
        ru = Uvect(2)
        rv = Uvect(3)
        e = Uvect(4)
        u = ru/r
        v = rv/r
        q = 0.5_PR * ( u**2 + v**2 )
        p = (gammagp - 1._PR)*(e - r*q)

        fluxFuncG(1) = rv
        fluxFuncG(2) = ru*v
        fluxFuncG(3) = rv*v + p
        fluxFuncG(4) = (e + p)*v
    End Function fluxFuncG

End Module mod_schemes
