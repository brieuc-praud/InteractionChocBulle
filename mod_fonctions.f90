Module mod_fonctions

   Use mod_params
   Implicit None

Contains

   Function D(x,y)
      Real(PR), Intent(In) :: x, y
      Real(PR) :: D

      D = 1._PR + x + y

   End Function D

   Function S(t,x,y)
      Real(PR), Intent(In) :: t,x, y
      Real(PR) :: S

      S = exp(-t) * ((1._PR + 2._PR*(x+y))*cos(x)*sin(y) - cos(x+y))

   End Function S

   Function Tinit(x,y)
      Real(PR), Intent(In) :: x, y
      Real(PR) :: Tinit

      Tinit = cos(x)*sin(y)

   End Function Tinit

   Function Touest(t,y)
      Real(PR), Intent(In) :: t,y
      Real(PR) :: Touest

      Touest = exp(-t)*sin(y)

   End Function Touest

   Function Test(t,y)
      Real(PR), Intent(In) :: t,y
      Real(PR) :: Test

      Test = exp(-t)*sin(y)

   End Function Test

   Function Tnord(t,x)
      Real(PR), Intent(In) :: t,x
      Real(PR) :: Tnord

      Tnord = 0._PR!exp(-t)*sin(x)

   End Function Tnord

   Function Tsud(t,x)
      Real(PR), Intent(In) :: t,x
      Real(PR) :: Tsud

      Tsud = 0._PR!exp(-t)*sin(x)

   End Function Tsud

End Module mod_fonctions
