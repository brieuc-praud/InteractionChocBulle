Module mod_test

   Use mod_parameters
   Use mod_schemes
   Implicit None

Contains
   Subroutine PerformTest(gamma, fluxFuncF, fluxFuncG)
      ! ------------------- Intent In -----------------------
      Real(PR), Intent(In) :: gamma
      Interface
         Function fluxFuncF(Uvect, gamma)
            Import PR
            Real(PR), Dimension(4), Intent(In) :: Uvect
            Real(PR), Intent(In) :: gamma
            Real(PR), Dimension(4) :: fluxFuncF
         End Function fluxFuncF
      End Interface
      Interface
         Function fluxFuncG(Uvect, gamma)
            Import PR
            Real(PR), Dimension(4), Intent(In) :: Uvect
            Real(PR), Intent(In) :: gamma
            Real(PR), Dimension(4) :: fluxFuncG
         End Function fluxFuncG
      End Interface
      ! ----------------------------------------------------------
      ! === Check fluxes consistency ===
      Write(*,*) '=== FLUXES CONSISTENCY ==='
      Write(*,*) "-- Rusanov"
      Call testNumericalFlux(Rusanov, gamma, fluxFuncF, fluxFuncG)
      Write(*,*) " -> SUCCESS" ! The program exits if a test fails

      Write(*,*) "-- HLL"
      Call testNumericalFlux(HLL, gamma, fluxFuncF, fluxFuncG)
      Write(*,*) " -> SUCCESS"
   End Subroutine PerformTest

   Subroutine PerformTestAndExit(gamma, fluxFuncF, fluxFuncG)
      ! ------------------- Intent In -----------------------
      Real(PR), Intent(In) :: gamma
      Interface
         Function fluxFuncF(Uvect, gamma)
            Import PR
            Real(PR), Dimension(4), Intent(In) :: Uvect
            Real(PR), Intent(In) :: gamma
            Real(PR), Dimension(4) :: fluxFuncF
         End Function fluxFuncF
      End Interface
      Interface
         Function fluxFuncG(Uvect, gamma)
            Import PR
            Real(PR), Dimension(4), Intent(In) :: Uvect
            Real(PR), Intent(In) :: gamma
            Real(PR), Dimension(4) :: fluxFuncG
         End Function fluxFuncG
      End Interface
      ! ----------------------------------------------------------
      Call PerformTest(gamma, fluxFuncF, fluxFuncG)
      Call Exit()
   End Subroutine PerformTestAndExit

   Function randomU()
      Real(PR), Parameter :: safety_factor = 2._PR
      Real(PR), Parameter :: eps = EPSILON(Real(PR))
      Real(PR), Parameter :: hug = 1._PR/eps
      Real(PR), Dimension(4), Parameter :: lowerBounds = &
      & (/ Real(PR) :: eps, -hug, -hug, eps /)
      Real(PR), Dimension(4), Parameter :: upperBounds = &
      & (/ Real(PR) :: hug, hug, hug, hug /)

      Real(PR), Dimension(4) :: randomUniform
      Real(PR), Dimension(4) :: randomU

      Real(PR) :: r, u, v, q

      Call RANDOM_NUMBER(randomUniform)
      randomU = upperBounds * randomUniform + lowerBounds * ( 1._PR - randomUniform )
      ! The lower bound for the energy depends on the three first variables
      r = randomU(1)
      u = randomU(2)/r
      v = randomU(3)/r
      q = 0.5_PR * ( u**2 + v**2 )
      randomU(4) = randomU(4) + r*q*(safety_factor - randomUniform(4))
   End Function randomU

   Subroutine testNumericalFlux(numericalFlux, gamma, fluxFuncF, fluxFuncG)
      Real(PR), Parameter :: safety_factor = 2._PR
      Integer, Parameter :: number_of_tests = 1000
      ! ------------------- Intent In -----------------------
      Real(PR), Intent(In) :: gamma
      Interface
         Function numericalFlux(Ui, Uip1, gamma, F)
            Import PR
            Real(PR), Dimension(4), Intent(In) :: Ui, Uip1
            Real(PR), Intent(In) :: gamma
            Real(PR), Dimension(4) :: numericalFlux
            Interface
               Function F(Uvect, gamma)
                  Import PR
                  Real(PR), Dimension(4), Intent(In) :: Uvect
                  Real(PR), Intent(In) :: gamma
                  Real(PR), Dimension(4) :: F
               End Function F
            End Interface
         End Function numericalFlux
      End Interface
      Interface
         Function fluxFuncF(Uvect, gamma)
            Import PR
            Real(PR), Dimension(4), Intent(In) :: Uvect
            Real(PR), Intent(In) :: gamma
            Real(PR), Dimension(4) :: fluxFuncF
         End Function fluxFuncF
      End Interface
      Interface
         Function fluxFuncG(Uvect, gamma)
            Import PR
            Real(PR), Dimension(4), Intent(In) :: Uvect
            Real(PR), Intent(In) :: gamma
            Real(PR), Dimension(4) :: fluxFuncG
         End Function fluxFuncG
      End Interface
      ! ----------------------------------------------------------
      Integer :: i
      Real(PR), Dimension(4) :: U
      Logical :: failF, failG

      Do i=0, number_of_tests
         U = randomU()
         failF = SUM( ABS(numericalFlux(U, U, gamma, fluxFuncF) - fluxFuncF(U, gamma)) ) > safety_factor*EPSILON(U(1))
         failG = SUM( ABS(numericalFlux(U, U, gamma, fluxFuncG) - fluxFuncG(U, gamma)) ) > safety_factor*EPSILON(U(1))
         If (failF .OR. failG) Then
            Call Exit(1)
         End If
      End Do
   End Subroutine testNumericalFlux

End Module mod_test
