Module mod_schemes

   Use mod_parameters
   Implicit None

Contains
   Function Rusanov(Ui, Uip1, gamma, F)
      Real(PR), Dimension(4), Intent(In) :: Ui, Uip1
      Real(PR), Intent(In) :: gamma
      Interface
         Function F(Uvect, gamma)
            Import PR
            Real(PR), Dimension(4), Intent(In) :: Uvect
            Real(PR), Intent(In) :: gamma
            Real(PR), Dimension(4) :: F
         End Function F
      End Interface
      Real(PR), Dimension(4) :: Rusanov
      ! --- Locals
      Real(PR) :: r, ru, rv, e, u, v, p, q, a, b

      r = Ui(1)
      ru = Ui(2)
      rv = Ui(3)
      e = Ui(4)
      u = ru/r
      v = rv/r
      q = 0.5_PR * ( u**2 + v**2 )
      p = (gamma - 1._PR)*(e - r*q)
      a = SQRT(gamma*p/r)

      b = MAX(ABS(u-a), ABS(u), ABS(u+a))

      Rusanov = 0.5_PR * ( F(Ui, gamma) + F(Uip1, gamma) - b*(Uip1 - Ui) )
   End Function Rusanov

   Function HLL(Ui, Uip1, gamma, F)
      Real(PR), Dimension(4), Intent(In) :: Ui, Uip1
      Real(PR), Intent(In) :: gamma
      Interface
         Function F(Uvect, gamma)
            Import PR
            Real(PR), Dimension(4), Intent(In) :: Uvect
            Real(PR), Intent(In) :: gamma
            Real(PR), Dimension(4) :: F
         End Function F
      End Interface
      Real(PR), Dimension(4) :: HLL
      ! --- Locals
      Real(PR) :: r, ru, rv, e, u, v, p, q, a
      Real(PR) :: l1, l2, l3, bm, bp

      r = Ui(1)
      ru = Ui(2)
      rv = Ui(3)
      e = Ui(4)
      u = ru/r
      v = rv/r
      q = 0.5_PR * ( u**2 + v**2 )
      p = (gamma - 1._PR)*(e - r*q)
      a = SQRT(gamma*p/r)

      l1 = ABS(u-a)
      l2 = ABS(u)
      l3 = ABS(u+a)
      bm = MIN( l1, l2, l3 )
      bp = MAX( l1, l2, l3 )

      If ( bm > 0 ) Then
         HLL = F(Ui, gamma)
      Else If ( bp < 0 ) Then
         HLL = F(Uip1, gamma)
      Else
         HLL = ( bp * F(Ui, gamma) - bm * F(Uip1, gamma)    +    bp * bm * (Uip1 - Ui) )/( bp - bm )
      End If
   End Function HLL

End Module mod_schemes
