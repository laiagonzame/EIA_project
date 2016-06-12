module PeriodicBoundaryConditions
implicit none
public

!Modulo para implementar condiciones periodicas de contorno
!en la dinamica molecular

contains

subroutine PBC(num_particulas,posicion,Lcaja,numproc,rank)
!Necesita de input el numero de particulas, la matriz con las coordenadas,
!y el lado de la caja de simulacion.
!Devuelve como output las matriz con las coordenadas de las particulas 
!o su imagen dentro de la caja

include 'mpif.h'

integer, intent(in) :: num_particulas, numproc, rank
double precision, dimension(3,num_particulas), intent(inout) :: posicion
double precision, intent(in) :: Lcaja
integer :: i,j,partxproc,ini(0:numproc-1),fin(0:numproc-1),size_m
integer :: stat(MPI_STATUS_SIZE),request(0:numproc-1),ierror


partxproc=nint (real(num_particulas)/real(numproc))
do i=0,numproc-2 
      ini(i)=i*partxproc+1
      fin(i)=(i+1)*partxproc
enddo
 ini(numproc-1)=(numproc-1)*partxproc+1
 fin(numproc-1)=num_particulas

   do j= ini(rank),fin(rank)
        do i=1,3
            If (posicion(i,j)>Lcaja/2.) then
                posicion(i,j)=posicion(i,j)-Lcaja
            else if (posicion(i,j)<-Lcaja/2.) then
                  posicion(i,j)=posicion(i,j)+Lcaja
            endif
    
       enddo
    enddo

!Primero juntamos todos los trozos de la matriz en el MASTER (rank=0)

if (rank .ne. 0) then
    request(rank)=0
    size_m=3*(fin(rank)-ini(rank)+1)
    print*,size_m
    print*,posicion(:,ini(rank):fin(rank))
    call MPI_ISEND(posicion(:,ini(rank):fin(rank)),size_m, MPI_REAL, 0, 1, MPI_COMM_WORLD, request(rank), ierror)
endif

call MPI_BARRIER(MPI_COMM_WORLD, ierror) !Esperamos que todos hayan mandado

!Hacemos la recepcion por parte del MASTER (rank=0)
if (rank .eq. 0) then

   do i =1, numproc-1
       size_m=3*(fin(i)-ini(i)+1)
       print*,size_m
       call MPI_RECV(posicion(:,ini(i):fin(i)), size_m, MPI_REAL, i, 1, MPI_COMM_WORLD, stat, ierror)
       print*,posicion(:,ini(i):fin(i))
   enddo


!Copiamos la matriz entera a todos los workers para que asi las siguientes partes del codigo paralelizadas funcionen con toda la informacion actualizada.
   do i=1, numproc-1
        request(0)=i
        call MPI_ISEND(posicion,3*num_particulas,MPI_REAL,i,1, MPI_COMM_WORLD,request(i),ierror)
   enddo
endif

call MPI_BARRIER(MPI_COMM_WORLD, ierror)


if (rank .ne. 0) then
     call MPI_RECV(posicion,3*num_particulas, MPI_REAL, 0, 1, MPI_COMM_WORLD, stat, ierror)
end if

call MPI_BARRIER(MPI_COMM_WORLD, ierror)

endsubroutine


end module
