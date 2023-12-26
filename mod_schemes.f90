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
        Real(PR) :: rL, velocity_uL, velocity_vL, eL, pressureL, qL, aL, b_minusL, b_plusL, l1L, l3L
        Real(PR) :: rR, velocity_uR, velocity_vR, eR, pressureR, qR, aR, b_minusR, b_plusR, l1R, l3R
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
            l1L = velocity_vL - aL
            l3L = velocity_vL + aL
            l1R = velocity_vR - aR
            l3R = velocity_vR + aR
        Case Default ! 'x'
            l1L = velocity_uL - aL
            l3L = velocity_uL + aL
            l1R = velocity_uR - aR
            l3R = velocity_uR + aR
        End Select

        b_minusL = MIN( l1L, l3L )
        b_minusR = MIN( l1R, l3R )
        b_minus = MIN( b_minusL, b_minusR, 0._PR )
        b_plusL = MAX( l1L, l3L )
        b_plusR = MAX( l1R, l3R )
        b_plus = MAX( b_plusL, b_plusR, 0._PR )

        Select Case (axis)
        Case ('y')
            HLL = b_plus * fluxFuncG(UL, gammagp) - b_minus * fluxFuncG(UR, gammagp)
            HLL = ( HLL   +   b_plus * b_minus * (UR - UL) ) / ( b_plus - b_minus )
        Case Default ! 'x'
            HLL = b_plus * fluxFuncF(UL, gammagp) - b_minus * fluxFuncF(UR, gammagp)
            HLL = ( HLL   +   b_plus * b_minus * (UR - UL) ) / ( b_plus - b_minus )
        End Select
    End Function HLL


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
