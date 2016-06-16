program main
use PeriodicBoundaryConditions
implicit none
include 'mpif.h'

!Programa per testejar que la subrutina de PBC esta ben implementada
!feta per : Cristina Roncero
!Atenció: si es modifica el valor del costat de la caixa aqui també s'ha 
!de modificar a l'script ("script_vmd_TESTPBC.tcl"( si es vol visualitzar.


!Lcaja= costat de la caixa
!posicion= matriu de coordenades 
!M= # de particules

real*8 :: Lcaja, t
real*8, dimension(:,:),allocatable :: posicion
integer :: M, seed(1), contador,rank,numproc,ierror
integer ::stat(MPI_STATUS_SIZE)

call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)


!----Paramatres---!
M=8
Lcaja=4
seed(1)=4321

allocate(posicion(3,M))


t=0.
contador=1
     call Ini_test(posicion)

     open (unit=20, file= "test_pbc.xyz")

if (rank .eq. 0) then
     call output
end if

t=1.
contador=2


call MPI_BARRIER(MPI_COMM_WORLD,ierror)
call PBC(M,posicion,Lcaja,numproc,rank)

if (rank .eq. 0) then
    call output
end if

call MPI_BARRIER (MPI_COMM_WORLD,ierror)
call MPI_FINALIZE(ierror)

contains

subroutine Ini_test (posicion)
!Subrutina per inicialitzar les posicions en un radi de dos caixes.
real*8, dimension(3,M), intent(out) :: posicion
integer :: i,count
real*8 :: aux(1)

!call system_clock(count)
!seed=count

call random_seed(put=seed)
	do i=1,M 
 
		call random_number(aux)
		posicion(1,i)=(2.*aux(1)-1.)*Lcaja
		call random_number(aux)
		posicion(2,i)=(2.*aux(1)-1.)*Lcaja
		call random_number(aux)
		posicion(3,i)=(2.*aux(1)-1.)*Lcaja

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

