module initial
implicit none
public

!Modulo que inicializa el sistema a partir de una densidad dada


contains

!N -- # de particulas
!pos -- matriz de posiciones
!vel -- matriz de velocidades
!box -- Longitud de la caja
!T -- temperatura 

!**************************************************************
subroutine inirandom(N,pos,vel,box,T)
!Subrutina que inicializa las particulas en la caja de manera
! aleatoria en posicion y velocidad

integer :: N,i,j,seed
real*8,dimension(3,N) ::pos,vel
real*8 :: box,aux,T

!Inicializacion de los numeros aleatorios
seed=1221
call random_seed(seed)

do j=1,N
   do i=1,3
      ! Posiciones
      call random_number(aux)
      pos(i,j)=(2.*aux-1.)*box
      ! Velocidades
      call random_number(aux)
      vel(i,j)=(2.*aux-1.)*dsqrt(T) 
   enddo
enddo

end subroutine

!***************************************************************
