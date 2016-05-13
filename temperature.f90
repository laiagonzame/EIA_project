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

subroutine compute_temperature(kbt ,vel, M, dt)
integer, intent(in) :: M
real*8, dimension(3,M),intent(in) :: vel
real*8, intent(in) :: dt
real*8, intent(out) :: kbt
integer :: i

kbt=0.
 
do i=1,M
! 1. Falta dividir pel numero de graus de llibertat: al main se
! li diu nf ---> nf = 3*M-3 
! 2. Si treballem en unitats de massa no caldria multiplicar-la aqui
 kbt = kbt +(vel(1,i)**2+vel(2,i)**2+vel(3,i)**2)/3.

end do



end subroutine

subroutine kinetic_energy(kine, kbt)
real*8, intent(in) :: kbt
real*8, intent(out) :: kine

kine=3.* kb t/ 2.

end subroutine



end module
