Program euler
   Use mod_functions
   Use mod_schemes
   Use mod_output
   Use mod_test

   Implicit None

   Character(14), Parameter :: parameters = "parameters.dat"

! Parameters
   Real(PR) xmin, xmax, ymin, ymax, tmax, cfl, gamma
   Integer :: imax, jmax
   Integer :: nplot, nb_fich
! Arrays
   Real(PR), Dimension(:,:,:), Allocatable :: Uvect, fluxF, fluxG
   Real(PR), Dimension(:), Allocatable :: x, y, xm, ym
! Loop indices
   Integer :: i, j, n
! Other
   Real(PR) :: deltax, deltay, deltat, tn
   Integer :: nmax
! Read parameters
   Open(111, File=parameters)
   Read(111, *) xmin, xmax, ymin, ymax, tmax, cfl, imax, jmax, gamma
   Close(111)
! Allocate
   Allocate(x(0:imax), y(0:jmax), xm(imax), ym(jmax))
   Allocate(Uvect(4,imax,jmax), fluxF(4,0:imax, 0:jmax), fluxG(4,0:imax, 0:jmax))
! Compute the grid
   deltax = (xmax - xmin) / imax
   deltay = (ymax - ymin) / jmax
   x = (/ Real(PR) :: (xmin + i*deltax, i=0, imax) /)
   y = (/ Real(PR) :: (ymin + j*deltay, j=0, jmax) /)
   xm = (/ Real(PR) :: (xmin + 0.5_PR*deltax + i*deltax, i=0, imax-1) /)
   ym = (/ Real(PR) :: (ymin + 0.5_PR*deltay + j*deltay, j=0, jmax-1) /)
! Initialise U
   Do i=1, imax
      Do j=1 , jmax
         Uvect(:,i,j) = Uinit(xm(i), ym(j))
      End Do
   End Do
! Compute time quantities
   Call compute_CFL(Uvect, deltax, deltay, deltat)
   nmax = Int(tmax / deltat)
   nb_fich = 1
   nplot = Int(nmax / 10)

   Call PerformTestAndExit(gamma, fluxFuncF, fluxFuncG)

! Boucle temporelle
   Call output(Uvect, x, y, 0)
   tn = 0._PR
   Do n=0, nmax
      tn = n * deltat
      Do j=1, jmax
         Do i=1, imax-1
            fluxF(:,i,j) = Rusanov(Uvect(:,i,j), Uvect(:,i+1,j), gamma, fluxFuncF)
            !fluxF(:,i,j) = HLL(Uvect(:,i,j), Uvect(:,i+1,j), gamma, fluxFuncF)
         End Do
         ! DIRICHLET
         !F(0,j) = -2._PR*D(xmin, ym(j)) * (T(1, j) - Touest(tn, ym(j))) / deltax
         !F(imax,j) = 2._PR*D(xmax, ym(j)) * (T(imax,j) - Test(tn,ym(j))) / deltax
         ! NEUMANN
         fluxF(:,imax,j) = (/ Real(PR) :: 0._PR, 0._PR, 0._PR, 0._PR /)
         fluxF(:,0,j) = (/ Real(PR) :: 0._PR, 0._PR, 0._PR, 0._PR /)
      End Do
      Do i=1, imax
         Do j=1, jmax-1
            fluxG(:,i,j) = Rusanov(Uvect(:,i,j), Uvect(:,i,j+1), gamma, fluxFuncG)
            fluxG(:,i,j) = HLL(Uvect(:,i,j), Uvect(:,i,j+1), gamma, fluxFuncG)
         End Do
         ! DIRICHLET
         !G(i,0) = -2._PR*D(xm(i), ymin) * ( T(i,1) - Tnord(tn,xm(i)) ) / deltay
         !G(i,jmax) = 2._PR*D(xm(i), ymax) * (T(i,jmax) - Tsud(tn,xm(i))) / deltay
         ! NEUMANN
         fluxG(:,i,0) = (/ Real(PR) :: 0._PR, 0._PR, 0._PR, 0._PR /)
         fluxG(:,i,jmax) = (/ Real(PR) :: 0._PR, 0._PR, 0._PR, 0._PR /)
      End Do
      Do i=1, imax
         Do j=1, jmax
            Uvect(:,i,j) = Uvect(:,i,j) &
            & - deltat/deltax * (fluxF(:,i,j) - fluxF(:,i-1,j)) &
            & - deltat/deltay * (fluxG(:,i,j) - fluxG(:,i,j-1))
         End Do
      End Do


      If ( Modulo(n, nplot) == 0) Then
         Call output(Uvect, x, y, nb_fich)
         nb_fich = nb_fich + 1
      End If

   End Do

   Deallocate(x, y, xm, ym)
   Deallocate(Uvect, fluxF, fluxG)

Contains
   Function fluxFuncF(Uvect, gamma)
      Real(PR), Dimension(4), Intent(In) :: Uvect
      Real(PR), Intent(In) :: gamma
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
      p = (gamma - 1._PR)*(e - r*q)

      fluxFuncF(1) = ru
      fluxFuncF(2) = ru*u + p
      fluxFuncF(3) = rv*u
      fluxFuncF(4) = (e + p)*u
   End Function fluxFuncF

   Function fluxFuncG(Uvect, gamma)
      Real(PR), Dimension(4), Intent(In) :: Uvect
      Real(PR), Intent(In) :: gamma
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
      p = (gamma - 1._PR)*(e - r*q)

      fluxFuncG(1) = rv
      fluxFuncG(2) = ru*v
      fluxFuncG(3) = rv*v + p
      fluxFuncG(4) = (e + p)*v
   End Function fluxFuncG

   Subroutine compute_CFL(U, dx, dy, dt)
      ! --- Parameters
      Real(PR), Parameter :: threshold = 1.e-6_PR
      ! --- InOut ---
      Real(PR), Dimension(:,:,:), Intent(In) :: U
      Real(PR), Intent(In) :: dx, dy
      Real(PR), Intent(Out) :: dt
      ! --- Locals ---
      Real(PR) :: rho, u2, u3, umax, vmax

      umax = 0._PR
      vmax = 0._PR
      Do i=1, imax
         Do j=1 , jmax
            rho = U(1,i,j)
            u2 = ABS( U(2,i,j) / rho )
            u3 = ABS( U(3,i,j) / rho )

            If (u2 > umax) Then
               umax = u2
            End If
            If (u3 > vmax) Then
               vmax = u3
            End If
         End Do
      End Do

      ! Avoid dividing by zero
      If (umax < threshold) Then
         umax = threshold
      End If
      If (vmax < threshold) Then
         vmax = threshold
      End If

      dt = 0.5_PR * MIN(dx/umax, dy/vmax)
   End Subroutine compute_CFL

End Program euler
