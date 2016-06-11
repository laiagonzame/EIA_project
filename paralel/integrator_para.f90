module integrator_para
implicit none
public

! Modulo que contiene un integrador tipo velocity verlet
! Como input acepta:
! el numero de particulas de la dinamica (N)
! una matriz (3,N) con las posiciones de las N particulas (coord) 
! una matriz (3,N) con las velocidades de las N particulas (vel)
! una matriz (3,N) con las fuerzas de las N particulas (frz)
! el paso de tiempo de la dinamica (ts)
! la massa de las particulas (m)
! Como output saca la matriz coord con las posiciones actualizadas


contains

subroutine Verlet_Coord(coord,frz,vel,N,ts,m,L,rank,numproc)
include 'mpif.h'
real*8, dimension(3,N),intent(in) :: frz,vel
real*8, dimension(3,N),intent(inout) :: coord
integer, intent(in) :: N, rank, numproc
real*8, intent(in) :: ts,m,L
integer :: i, dproc, N_end, MASTER, ierror, request, N_ini
real*8, dimension(3) :: dr
integer :: stat(MPI_STATUS_SIZE)

! Subrutina que integra la posicion
! Como input acepta:
! el numero de particulas de la dinamica (N)
! una matriz (3,N) con las posiciones de las N particulas (coord) 
! una matriz (3,N) con las velocidades de las N particulas (vel)
! una matriz (3,N) con las fuerzas de las N particulas (frz)
! el paso de tiempo de la dinamica (ts)
! la massa de las particulas (m)
! la longitud de la caja de simulacion (L)
! Como output saca la matriz coord con las posiciones actualizadas




dproc=nint(real(N)/real(numproc))
MASTER=0
N_ini=(rank*dproc)+1

if (rank /= (numproc-1)) then

N_end=(rank+1)*dproc

else 

N_end=N

end if

do i=N_ini,N_end
    
 dr=coord(:,i)+vel(:,i)*ts+(ts**2)*frz(:,i)/(2*m)
 call  Check_properdisplacement(dr,L)
 coord(:,i)=dr

end do

! el trabajador rank le envia al master la matriz que ha integrado
if (rank .ne. MASTER) then
   call MPI_ISEND(coord(:,N_ini:N_end), 1, MPI_INTEGER,MASTER, 1, MPI_COMM_WORLD, request, ierror)
end if


call MPI_BARRIER(MPI_COMM_WORLD, ierror)

! el master recive los mensajes de los workers y los junta en una sola matriz
if (rank .eq. MASTER) then
   do i = 1, numproc-1
      
   call MPI_RECV(coord(:,N_ini:N_end), 1, MPI_INTEGER, i, 1,MPI_COMM_WORLD, stat, ierror)

   end do

! el master envia la matriz completa integrada a todos los workers

   do i = 1, numproc-1
      call MPI_ISEND(coord, 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD, request, ierror)
   end do

end if

call MPI_BARRIER(MPI_COMM_WORLD, ierror)

! el trabajador rank recive la matriz completa integrada 

if (rank .ne. MASTER) then
 
  call MPI_RECV(coord, 1, MPI_INTEGER, MASTER, 1, MPI_COMM_WORLD, stat, ierror)
 
end if

call MPI_BARRIER(MPI_COMM_WORLD, ierror)



end subroutine

subroutine Verlet_Vel(frz,frzt,vel,N,ts,m,rank,numproc)
include 'mpif.h'
real*8, dimension(3,N), intent(in) :: frz,frzt
real*8, dimension(3,N), intent(inout) :: vel
integer, intent(in):: N,rank,numproc
real*8, intent(in)::ts,m
integer:: i, dproc, N_end, MASTER, N_ini, request, ierror
integer :: stat(MPI_STATUS_SIZE)

! Subrutina que integra la posicion
! Como input acepta:
! el numero de particulas de la dinamica (N)
! una matriz (3,N) con las velocidades de las N particulas (vel)
! una matriz (3,N) con las fuerzas de las N particulas en el paso de tiempo anterior(frz)
! una matriz (3,N) con las fuerzas de las N particulas en el paso de tiempo actual(frzt)
! el paso de tiempo de la dinamica (ts)
! la massa de las particulas (m)

dproc=nint(real(N)/real(numproc))
MASTER=0
N_ini=(rank*dproc)+1

if (rank /= (numproc-1)) then

N_end=(rank+1)*dproc

else 

N_end=N

end if


do i=N_ini,N_end
   
  vel(:,i)=vel(:,i)+(frz(:,i)+frzt(:,i))*ts/(2*m)

end do


! el trabajador rank le envia al master la matriz que ha integrado
 

if (rank .ne. MASTER) then
   call MPI_ISEND(vel(:,N_ini:N_end), 1, MPI_INTEGER,MASTER, 1, MPI_COMM_WORLD, request, ierror)
end if



call MPI_BARRIER(MPI_COMM_WORLD, ierror)

!  el master recive los mensajes de los workers y los junta en una sola matriz

if (rank .eq. MASTER) then
   do i = 1, numproc-1
      
   call MPI_RECV(vel(:,N_ini:N_end), 1, MPI_INTEGER, i, 1,MPI_COMM_WORLD, stat, ierror)

   end do
  ! el master envia la matriz completa integrada a todos los workers

   do i = 1, numproc-1
      call MPI_ISEND(vel, 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD, request, ierror)
   end do

end if

call MPI_BARRIER(MPI_COMM_WORLD, ierror)

! el trabajador rank recive la matriz completa integrada

if (rank .ne. MASTER) then
 
  call MPI_RECV(vel, 1, MPI_INTEGER, MASTER, 1, MPI_COMM_WORLD, stat, ierror)
 
end if

call MPI_BARRIER(MPI_COMM_WORLD, ierror)


end subroutine

subroutine Check_properdisplacement(dr,L)
real*8,dimension(3), intent (in) :: dr
real*8, intent(in) :: L
! Subrutina que comprueba que la particula no se desplace una longitud mayor que la mitad de la caja de simulacion
! De ser asi, el programa da un mensaje de error y termina
! como input acepta: 
! el desplazamiento de la particula (dr)
! la longitud de la caja de simulacion (L)


if ( (sqrt(dr(1)**2+dr(2)**2+dr(3)**2) > L) ) then
 print *, "ERROR: desplazamiento mayor que la mitad de la longitud de la caja de simulacion, considera reducir el paso de tiempo"
 print *, "Program end: fatal error"
 stop 
end if

end subroutine 


end module
