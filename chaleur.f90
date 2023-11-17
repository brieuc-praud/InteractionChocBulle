Program chaleur

   Use mod_fonctions
   Use mod_sorties

   Implicit None

   Character(14), Parameter :: fichier_parametres = "parametres.dat"

! Paramètres
   Real(PR) xmin, xmax, ymin, ymax, tmax, cfl
   Integer :: imax, jmax
   Integer :: nplot, nb_fich
! Tableaux
   Real(PR), Dimension(:,:), Allocatable :: T, F, G
   Real(PR), Dimension(:), Allocatable :: x, y, xm, ym
! Indices de boucles
   Integer :: i, j, n
! Autre
   Real(PR) :: Dmax
   Real(PR) :: deltax, deltay, deltat, tn
   Integer :: nmax

! Lecture des parametres
   Open(111, File=fichier_parametres)
   Read(111, *) xmin, xmax, ymin, ymax, tmax, cfl, imax, jmax
   Close(111)
! Allocations
   Allocate(x(0:imax), y(0:imax), xm(imax), ym(imax))
   Allocate(T(imax,jmax), F(0:imax, jmax), G(imax, 0:jmax))
! Calcul de deltax et de deltay
   deltax = (xmax - xmin) / imax
   deltay = (ymax - ymin) / jmax
! Calcul de x, y, xm et ym
   x = (/ (xmin + i*deltax, i=0, imax) /)
   y = (/ (ymin + j*deltay, j=0, jmax) /)
   xm = (/ (xmin + 0.5_PR*deltax + i*deltax, i=0, imax-1) /)
   ym = (/ (ymin + 0.5_PR*deltay + j*deltay, j=0, jmax-1) /)
! Calcul de deltat, de nmax et de nplot
   Dmax = 0._PR
   Do i = 0, imax
      Do j = 0, jmax
         If ( D(x(i),y(j)) > Dmax ) Then
            Dmax = D(x(i),y(j))
         End If
      End Do
   End Do
   deltat = cfl / ( 2._PR * Dmax * (1._PR / deltax**2 + 1._PR / deltay**2) )
   nmax = Int(tmax / deltat)
   nb_fich = 1
   nplot = Int(nmax / 10)
! Initialisation de T
   Do i=1, imax
      Do j=1 , jmax
         T(i,j) = Tinit( xm(i), ym(j) )
      End Do
   End Do
! Boucle temporelle
   Call sorties_fichier(T, x, y, 0)
   Do n=0, nmax
      tn = n * deltat
      Do j=1, jmax
         Do i=1, imax-1
            F(i,j) = -D(x(i), ym(j)) * (T(i+1, j) - T(i,j)) / deltax
         End Do
         ! DIRICHLET
         !F(0,j) = -2._PR*D(xmin, ym(j)) * (T(1, j) - Touest(tn, ym(j))) / deltax
         !F(imax,j) = 2._PR*D(xmax, ym(j)) * (T(imax,j) - Test(tn,ym(j))) / deltax
         ! NEUMANN
         F(imax,j) = 0._PR
         F(0,j) = 0._PR
      End Do
      Do i=1, imax
         Do j=1, jmax-1
            G(i,j) = -D(xm(i), y(j)) * (T(i, j+1) - T(i,j)) / deltay
         End Do
         ! DIRICHLET
         !G(i,0) = -2._PR*D(xm(i), ymin) * ( T(i,1) - Tnord(tn,xm(i)) ) / deltay
         !G(i,jmax) = 2._PR*D(xm(i), ymax) * (T(i,jmax) - Tsud(tn,xm(i))) / deltay
         ! NEUMANN
         G(i,0) = 0._PR
         G(i,jmax) = 0._PR
      End Do
      Do i=1, imax
         Do j=1, jmax
            T(i,j) = T(i,j)   - deltat/deltax * (F(i,j) - F(i-1,j)) &
            & - deltat/deltay * (G(i,j) - G(i,j-1)) &
            & + deltat * S( tn, xm(i), ym(j) )
         End Do
      End Do


      If ( Modulo(n, nplot) == 0) Then
         Call sorties_fichier(T, x, y, nb_fich)
         nb_fich = nb_fich + 1
      End If

   End Do

   Deallocate(x, y, xm, ym)
   Deallocate(T, F, G)

Contains

End Program chaleur
