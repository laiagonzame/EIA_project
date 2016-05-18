module integrator
implicit none
public

! Modulo que contiene un integrador tipo velocity verlet
! Como input acepta:
! el numero de particulas de la dinamica (N)
! una matriz (3,N) con las posiciones de las N particulas (coord) 
! una matriz (3,N) con las velocidades de las N particulas (vel)
! una matriz (3,N) con las posiciones de las N particulas (frz)
! el paso de tiempo de la dinamica (ts)
! la massa de las particulas (m)
! Como output saca la matriz coord con las posiciones actualizadas



contains

subroutine Verlet_Coord(coord,frz,vel,N,ts,m)
real*8, dimension(3,N),intent(in) :: frz,vel
real*8, dimension(3,N),intent(inout) :: coord
integer, intent(in) :: N
real*8, intent(in) :: ts,m
integer :: i

do i=1,N
   coord(:,i)=coord(:,i)+vel(:,i)*ts+(ts**2)*frz(:,i)/(2*m)
end do

end subroutine

subroutine Verlet_Vel(frz,frzt,vel,N,ts,m)
real*8, dimension(3,N),intent(in) :: frz,frzt
real*8, dimension(3,N),intent(inout) :: vel
integer, intent(in):: N
real*8,intent(in)::ts,m
integer:: i
! Calculo de la velocidad a t+ts utilizando el integrador Velocity Verlet


do i=1,N
   vel(:,i)=vel(:,i)+(frz(:,i)+frzt(:,i))*ts/(2*m)
end do

end subroutine

end module
