Module mod_functions

   Use mod_parameters
   Implicit None

Contains
   Function Uinit(x,y)
      Real(PR), Intent(In) :: x, y
      Real(PR), Dimension(4) :: Uinit

      Uinit = (/ 1._PR, 0._PR, 0._PR, 1._PR /)
   End Function Uinit
End Module mod_functions
