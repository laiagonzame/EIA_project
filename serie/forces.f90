module forces_routines
implicit none

contains

subroutine forces(M,r,F,boxlength,sigma,epsil,epot)

integer, intent(in) :: M
double precision, dimension(3,M), intent(in) :: r
double precision, dimension(3,M), intent(out) :: F
double precision, intent(in) :: boxlength, sigma, epsil
double precision, dimension(3) :: rij 
double precision :: rcut, rijl, dist2, sigmar, rc2, force,epot
integer :: i, j, l

! rcut es el umbral, si la distancia es mayor a rc no se tiene en cuenta la interacción
rcut=3.*sigma 
rc2=rcut*rcut

! inicializacion de la fuerza y la energia potencial
F = 0d0
epot = 0d0

do i=1,M-1
   do j=i+1,M
      ! Calculamos la distancia
      dist2=0d0
      do l=1,3
         rijl = r(l,j) - r(l,i)
         rij(l) = rijl - boxlength*nint(rijl/boxlength)
         dist2 = dist2 + rij(l)*rij(l)
      enddo
      ! Comprovamos que la distancia es menor a rc
      if (dist2<rc2) then
         sigmar = (dist2**(-3d0))*(sigma**6d0)
         ! Calculamos el modulo de la fuerza dividido por pa distancia
         force = 24d0*epsil*(2d0*sigmar**(2d0)-sigmar)/dist2
         ! Actualizamos a las fuerzas debido a esta interaccion
         do l=1,3
            F(l,i) = F(l,i) - force*rij(l)
            F(l,j) = F(l,j) + force*rij(l)
         enddo
         ! Actualizamos a energia potencial debido a esta interaccion
         epot=epot+4d0*epsil*(sigmar**(2.)-sigmar)
      endif
   enddo   
enddo

endsubroutine

end module
