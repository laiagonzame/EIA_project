module integrator_para
implicit none
public

! Modulo que contiene un integrador tipo velocity verlet
!
! Realizado por Pablo M. Blanco, Universidad de Barcelona

contains

subroutine Verlet_Coord(coord,frz,vel,N,ts,m,L,rank,numproc)
include 'mpif.h'
real*8, dimension(3,N),intent(in) :: frz,vel
real*8, dimension(3,N),intent(inout) :: coord
integer, intent(in) :: N, rank, numproc
real*8, intent(in) :: ts,m,L
integer :: i, dproc,  MASTER, ierror, request
integer, dimension(:), allocatable :: N_ini,N_end, siz
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


! Se calcula el numero de particulas por procesador redondeado al entero 
! mas cercano

dproc=nint(real(N)/real(numproc))
MASTER=0 ! el procesador 0 sera el master

allocate ( N_ini(numproc), N_end(numproc), siz(numproc))

! se asigna el rango de particulas de cada procesador i 
! de ((i-1)*dproc)+1 a i*dproc
! menos el ultimo procesador que se ajusta para tener en cuenta
! el redondeo

do i=1,numproc


N_ini(i)=((i-1)*dproc)+1


if ( i /= numproc) then

N_end(i)=i*dproc

else 

N_end(i)=N

end if

siz(i)=(N_end(i)-N_ini(i)+1)*3

end do

! Bucle de la integracion

do i=N_ini(rank+1),N_end(rank+1)
    
 dr=vel(:,i)*ts+(ts**2)*frz(:,i)/(2*m)
 call  Check_properdisplacement(dr,L)
 coord(:,i)=dr+coord(:,i)

end do

! el trabajador rank le envia al master la matriz que ha integrado
if (rank .ne. MASTER) then
   call MPI_ISEND(coord(:,N_ini(rank+1):N_end(rank+1)), siz(rank+1), MPI_REAL8,MASTER, 1, MPI_COMM_WORLD, request, ierror)
end if


call MPI_BARRIER(MPI_COMM_WORLD, ierror)

! el master recive los mensajes de los workers y los junta en una sola matriz
if (rank .eq. MASTER) then
   do i = 1, numproc-1
      
   call MPI_RECV(coord(:,N_ini(i+1):N_end(i+1)), siz(i+1), MPI_REAL8, i, 1,MPI_COMM_WORLD, stat, ierror)

   end do

! el master envia la matriz completa integrada a todos los workers

   do i = 1, numproc-1
      call MPI_ISEND(coord, N*3, MPI_REAL8, i, 1, MPI_COMM_WORLD, request, ierror)
   end do

end if

call MPI_BARRIER(MPI_COMM_WORLD, ierror)

! el trabajador rank recive la matriz completa integrada 

if (rank .ne. MASTER) then
 
  call MPI_RECV(coord, N*3, MPI_REAL8, MASTER, 1, MPI_COMM_WORLD, stat, ierror)
 
end if

call MPI_BARRIER(MPI_COMM_WORLD, ierror)



end subroutine

subroutine Verlet_Vel(frz,frzt,vel,N,ts,m,rank,numproc)
include 'mpif.h'
real*8, dimension(3,N), intent(in) :: frz,frzt
real*8, dimension(3,N), intent(inout) :: vel
integer, intent(in):: N,rank,numproc
real*8, intent(in)::ts,m
integer:: i, dproc, MASTER, request, ierror
integer :: stat(MPI_STATUS_SIZE)
integer, dimension(:), allocatable :: N_ini,N_end,siz
! Subrutina que integra la posicion
! Como input acepta:
! el numero de particulas de la dinamica (N)
! una matriz (3,N) con las velocidades de las N particulas (vel)
! una matriz (3,N) con las fuerzas de las N particulas en el paso de tiempo anterior(frz)
! una matriz (3,N) con las fuerzas de las N particulas en el paso de tiempo actual(frzt)
! el paso de tiempo de la dinamica (ts)
! la massa de las particulas (m)

! Se calcula el numero de particulas por procesador redondeado al entero 
! mas cercano

dproc=nint(real(N)/real(numproc))
MASTER=0 ! se asigna el procesador 0 como master


allocate ( N_ini(numproc), N_end(numproc),siz(numproc))


! se asigna el rango de particulas de cada procesador i 
! de ((i-1)*dproc)+1 a i*dproc
! menos el ultimo procesador que se ajusta para tener en cuenta
! el redondeo

do i=1,numproc


N_ini(i)=((i-1)*dproc)+1

if (i == numproc) then

N_end(i)=i*dproc

else 

N_end(i)=N

end if

siz(i)=(N_end(i)-N_ini(i)+1)*3
end do

! Bucle de la integracion

do i=N_ini(rank+1),N_end(rank+1)

   
  vel(:,i)=vel(:,i)+(frz(:,i)+frzt(:,i))*ts/(2*m)

end do


! el trabajador rank le envia al master la matriz que ha integrado
 

if (rank .ne. MASTER) then
   call MPI_ISEND(vel(:,N_ini(rank+1):N_end(rank+1)), siz(rank+1), MPI_REAL8,MASTER, 1, MPI_COMM_WORLD, request, ierror)
end if



call MPI_BARRIER(MPI_COMM_WORLD, ierror)

!  el master recive los mensajes de los workers y los junta en una sola matriz

if (rank .eq. MASTER) then
   do i = 1, numproc-1
      
   call MPI_RECV(vel(:,N_ini(i+1):N_end(i+1)), siz(i+1), MPI_REAL8, i, 1,MPI_COMM_WORLD, stat, ierror)

   end do
  ! el master envia la matriz completa integrada a todos los workers

   do i = 1, numproc-1
      call MPI_ISEND(vel, N*3, MPI_REAL8, i, 1, MPI_COMM_WORLD, request, ierror)
   end do

end if

call MPI_BARRIER(MPI_COMM_WORLD, ierror)

! el trabajador rank recive la matriz completa integrada

if (rank .ne. MASTER) then
 
  call MPI_RECV(vel, N*3, MPI_REAL8, MASTER, 1, MPI_COMM_WORLD, stat, ierror)
 
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
