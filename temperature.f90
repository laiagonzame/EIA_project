module temperature
implicit none
public

! Modul que calcula la temperatura del sistema, kbt
! Com a input accepta:
! el número de particulas de la dinamica (M)¡
! una matriu (3,M) amb les velocitats de las N particulas (vel)
! el paso de tiempo de la dinamica (dt)
! la massa de las particulas (mass)
! Retorna la velocitat mitja de les partícules i la temperatura  (kbt)

contains

subroutine compute_temperature(kbt, vel, M, nf)
integer, intent(in) :: M, nf
double precision, dimension(3,M),intent(in) :: vel
double precision, intent(out) :: kbt
integer :: i

kbt=0d0
 
do i=1,M
   kbt = kbt +(vel(1,i)**2+vel(2,i)**2+vel(3,i)**2)/3.
end do

kbt = kbt / dfloat(nf)

end subroutine

subroutine kinetic_energy(kine, kbt)
double precision, intent(in) :: kbt
double precision, intent(out) :: kine
 
kine = 3.* kbt/ 2.

end subroutine

end module
