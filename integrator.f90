module integrator
implicit none
public

! Modulo que contiene un integrador tipo velocity verlet
! Como input acepta:
! el numero de particulas de la dinamica (N)
! una matriz (N,3) con las posiciones de las N particulas (coord) 
! una matriz (N,3) con las velocidades de las N particulas (vel)
! una matriz (N,3) con las posiciones de las N particulas (frz)
! el paso de tiempo de la dinamica (ts)
! la massa de las particulas (m)
! Como output saca la matriz coord con las posiciones actualizadas



contains

subroutine Verlet_integrator(coord,frz,vel,N,ts,m)
 real*8, dimension(N,3),intent(in) :: frz,vel
 real*8, dimension(N,3),intent(inout) :: coord
 integer, intent(in) :: N
 real*8, intent(in) :: ts,m
 integer :: i


 do i=1,N

  coord(i,1)=coord(i,1)+vel(i,1)*ts+(ts**2)*frz(i,1)/(2*m)
  coord(i,2)=coord(i,2)+vel(i,2)*ts+(ts**2)*frz(i,2)/(2*m)
  coord(i,3)=coord(i,3)+vel(i,3)*ts+(ts**2)*frz(i,3)/(2*m)

 end do



end subroutine
end module
