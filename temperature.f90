module temperature
implicit none
public

! Modul que calcula la temperatura del sistema, kbt
! Com a input accepta:
! el número de particulas de la dinamica (N)¡
! una matriu (N,3) amb les velocitats de las N particulas (vel)
! el paso de tiempo de la dinamica (dt)
! la massa de las particulas (m)
! Retorna la velocitat mitja de les partícules i la temperatura  (kbt)



contains

subroutine compute_temperature(kbt,vel,N,dt,m)
 real*8, dimension(3,N),intent(in) :: vel
 integer, intent(in) :: N
 real*8, intent(in) :: dt,m
 integer :: i
 real*8, intent(inout) :: kbt

 kbt=0.;
 
 do i=1,N

  kbt=kbt+m*(vel(1,i)**2+vel(2,i)**2+vel(3,i)**2)/3.

 end do



end subroutine

subroutine kinetic_energy(kine,kbt)
 real*8, intent(in) :: kbt
 real*8, intent(out) :: kine

 kine=3.*kbt/2.;

end subroutine



end module
