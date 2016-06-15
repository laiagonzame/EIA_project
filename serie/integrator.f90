module integrator
implicit none
public

! Modulo que contiene un integrador tipo velocity verlet
!
! Desarrollado por Pablo M. Blanco, Univerisdad de Barcelona

contains

subroutine Verlet_Coord(coord,frz,vel,N,ts,m,L)
real*8, dimension(3,N),intent(in) :: frz,vel
real*8, dimension(3,N),intent(inout) :: coord
integer, intent(in) :: N
real*8, intent(in) :: ts,m,L
integer :: i
real*8, dimension(3) :: dr
! Subrutina que integra la posicion
! Como input acepta:
! el numero de particulas de la dinamica (N)
! una matriz (3,N) con las posiciones de las N particulas (coord) 
! una matriz (3,N) con las velocidades de las N particulas (vel)
! una matriz (3,N) con las fuerzas de las N particulas (frz)
! el paso de tiempo de la dinamica (ts)
! la massa de las particulas (m)
! la longitud de la caja de simulacion (L)
! Como output saca la matriz coord con las posiciones actualizadas



do i=1,N
   dr=vel(:,i)*ts+(ts**2)*frz(:,i)/(2*m)
   call  Check_properdisplacement(dr,L)
   coord(:,i)=dr+coord(:,i)
end do

end subroutine

subroutine Verlet_Vel(frz,frzt,vel,N,ts,m)
real*8, dimension(3,N), intent(in) :: frz,frzt
real*8, dimension(3,N), intent(inout) :: vel
integer, intent(in):: N
real*8, intent(in)::ts,m
integer:: i

! Subrutina que integra la velocidad
! Como input acepta:
! el numero de particulas de la dinamica (N)
! una matriz (3,N) con las velocidades de las N particulas (vel)
! una matriz (3,N) con las fuerzas de las N particulas en el paso de tiempo anterior(frz)
! una matriz (3,N) con las fuerzas de las N particulas en el paso de tiempo actual(frzt)
! el paso de tiempo de la dinamica (ts)
! la massa de las particulas (m)
! Como output saca  la velocidad actualizada t+ts (vel)

do i=1,N
   vel(:,i)=vel(:,i)+(frz(:,i)+frzt(:,i))*ts/(2*m)
end do

end subroutine

subroutine Check_properdisplacement(dr,L)
real*8,dimension(3), intent (in) :: dr
real*8, intent(in) :: L
! Subrutina que comprueba que la particula no se desplace una longitud mayor que la mitad de la caja de simulacion
! De ser asi, el programa da un mensaje de error y termina
! como input acepta: 
! el desplazamiento de la particula (dr)
! la longitud de la caja de simulacion (L)


if ( (sqrt(dr(1)**2+dr(2)**2+dr(3)**2) > L/2.) ) then
 print *, "ERROR: desplazamiento mayor que la mitad de la longitud de la caja de simulacion, considera reducir el paso de tiempo"
 print *, "Program end: fatal error"
 stop 
end if

end subroutine 


end module
