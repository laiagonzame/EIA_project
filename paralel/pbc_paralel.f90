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
integer :: i,j,partxproc,ini,fin
integer :: stat(MPI_STATUS_SIZE),request,ierror

partxproc=nint (real(num_particulas)/real(numproc))

If (rank .eq.(numproc-1)) then
     ini=rank*partxproc+1
     fin=num_particulas
else
     ini=rank*partxproc+1
     fin=(rank+1)*partxproc
end if

   do j= ini,fin
       do i=1,3
            If (posicion(i,j)>Lcaja/2.) then
                posicion(i,j)=posicion(i,j)-Lcaja
                stop
            else if (posicion(i,j)<-Lcaja/2.) then
                posicion(i,j)=posicion(i,j)+Lcaja
            endif
    
       enddo
    enddo

!Primero juntamos todos los trozos de la matriz en el MASTER (rank=0)

if (rank .ne. 0) then
    call MPI_ISEND(posicion(:,ini:fin),1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, request, ierror)
endif

call MPI_BARRIER(MPI_COMM_WORLD, ierror) !Esperamos que todos hayan mandado

!Hacemos la recepcion por parte del MASTER (rank=0)
if (rank .eq. 0) then
   do i =1, numproc-1
       call MPI_RECV(posicion(:,ini:fin), 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD, stat, ierror)
   enddo


!Copiamos la matriz entera a todos los workers para que asi las siguientes partes del codigo paralelizadas funcionen con toda la informacion actualizada.
   do i=1, numproc-1
        call MPI_ISEND(posicion,1,MPI_INTEGER,i,1, MPI_COMM_WORLD,request,ierror)
   enddo

endif

call MPI_BARRIER(MPI_COMM_WORLD, ierror)


if (rank .ne. 0) then
     call MPI_RECV(posicion,1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, stat, ierror)
end if

call MPI_BARRIER(MPI_COMM_WORLD, ierror)

endsubroutine

end module
