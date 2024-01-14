Module mod_parameters

   Implicit None

   Integer, Parameter :: PR = 8
   Integer, Parameter :: STDERR = 0, STDIN = 5, STDOUT = 6, FILEINPUT = 500, FILEOUTPUT = 600
   Real(PR), Parameter :: PI = 4._PR * atan(1._PR)

   Integer, Parameter :: quadrature_points_number = 4 ! Achieve a 3rd order approximation

   Type :: numerical_scheme
       ! Space
       Character(len=20) :: space_scheme_name = "Rusanov"
       Character(len=20) :: MUSCL_limiter = "minmod"
       Real(PR) :: MUSCL_generalised_minmod_parameter = 1._PR
       Integer :: space_scheme_order = 1
       ! Time
       Integer :: time_scheme_order = 1
   End Type numerical_scheme

   ! To read in the parameter file
   Real(PR) :: time_max, cfl, gammagp
   Integer :: imax, jmax, output_modulo
   Character(len=20) :: case_name, norm_str
   Type(numerical_scheme) :: num_scheme
   ! Arrays
   Real(PR), Dimension(:,:,:), Allocatable :: Uvect, Uvect_e, K1vect, K2vect, K3vect
   Real(PR), Dimension(:,:,:), Allocatable :: fluxF, fluxG, fluxK1F, fluxK1G, fluxK2F, fluxK2G, fluxK3F, fluxK3G
   Real(PR), Dimension(:,:,:), Allocatable :: source
   Real(PR), Dimension(:), Allocatable :: x, y, xm, ym
   ! Other
   Real(PR) :: xmin, xmax, ymin, ymax
   Real(PR) :: deltax, deltay, deltat, time
   Logical :: exact_solution_available

Contains

End Module mod_parameters
