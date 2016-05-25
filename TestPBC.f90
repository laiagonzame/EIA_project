program main

use PeriodicBoundaryConditions

implicit none

!Programa per testejar que la subrutina de PBC esta ben implementada
!Atenció: si es modifica el valor del costat de la caixa aqui també s'ha 
!de modificar a l'script ("script_vmd_TESTPBC.tcl"( si es vol visualitzar.


!Lcaja= costat de la caixa
!posicion= matriu de coordenades 
!M= # de particules

real*8 :: Lcaja, t
real*8, dimension(:,:),allocatable :: posicion
integer :: M, seed, contador


!----Paramatres---!
M=8
Lcaja=4
seed=4321

allocate(posicion(3,M))


t=0.
contador=1

call Ini_test(posicion)

open (unit=20, file= "test_pbc.xyz")

call output


t=1.
contador=2

call PBC(M,posicion,Lcaja)
call output

contains

subroutine Ini_test (posicion)
!Subrutina per inicialitzar les posicions en un radi de dos caixes.
real*8, dimension(3,M), intent(out) :: posicion
integer :: i
real*8 :: aux

!inicialitzar el rand
aux=rand(seed)


	do i=1,M 
 
		aux=rand()
		posicion(1,i)=(2.*aux-1.)*Lcaja
		aux=rand()
		posicion(2,i)=(2.*aux-1.)*Lcaja
		aux=rand()
		posicion(3,i)=(2.*aux-1.)*Lcaja

	enddo


endsubroutine

subroutine output
!Subrutina que me guarda las coordenadas para que las pueda ver con el vmd
integer :: i


write(20,*) M
write(20,*) "trajectoria,iter:",contador,"temps",t,"ns"

   do i=1,M
   Write(20,*) "C ", posicion(1,i), posicion(2,i), posicion(3,i)
   enddo

end subroutine

end program

