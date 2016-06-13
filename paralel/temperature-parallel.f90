module temperature
implicit none
public

! Modul que calcula la temperatura del sistema, kbt
! Com a input accepta:
! el número de particulas de la dinamica (M)
! una matriu (3,M) amb les velocitats de las N particulas (vel)
! el paso de tiempo de la dinamica (dt)
! la massa de las particulas (mass)
! Retorna la velocitat mitja de les partícules i la temperatura  (kbt)

contains


subroutine compute_temperature(kbt, vel, M, nf)
include 'mpif.h'
integer, intent(in) :: M,nf
double precision, dimension(3,M), intent(in) :: vel
!f2py depend(M) :: vel
double precision, intent(out) :: kbt
integer :: i,Mtask,Mini,MASTER,ierror,taskid,numproc,aerr,partner,sender
 integer, allocatable, dimension(:) :: request
 integer stat(MPI_STATUS_SIZE)
double precision message

!--------- MPI inicialization ------------------
call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierror)
 allocate (request(numproc), STAT=aerr)

 MASTER = 0
 ! dividing work
Mtask=M/numproc
 Mini=taskid*Mtask+1
 if (taskid .eq. numproc-1) then 
   Mtask=Mtask+mod(M,numproc)
 endif
! print *,"Hello from task ", taskid," of ", numproc, " Mtask=",Mtask," of ",M
!---------------------------

!--------- do their part of the job ------------------
message=0d0
do i=1,Mtask
   message = message + (vel(1,Mini+i)**2+vel(2,Mini+i)**2+vel(3,Mini+i)**2)/dfloat(nf)
    if (taskid .eq. numproc-1) then
     endif
end do

!--------- Barrier: waiting for everyone ------------------
 call MPI_BARRIER(MPI_COMM_WORLD, ierror)
 if (taskid .eq. numproc-1) then
 endif
 
 !--------- send information to everyone ------------------
 
  do i=0,numproc-1
    request(i+1) = i 
    partner = i
    call MPI_ISEND(message, 1, MPI_INTEGER, partner, 1, MPI_COMM_WORLD, request(i+1),ierror)
 enddo
 
 !--------- receive kbt contributions and join them ------------------
 kbt=0d0
  do i=0,numproc-1
    sender = i
    call MPI_RECV(message, 1, MPI_INTEGER, sender, 1, MPI_COMM_WORLD, stat, ierror)
    kbt=kbt+message
    ! Print partner info and continue
 enddo

  deallocate (request, STAT=aerr)




end subroutine


subroutine kinetic_energy(kine, kbt, M, nf)
integer, intent(in) :: M, nf
double precision, intent(in) :: kbt
double precision, intent(out) :: kine
 
kine = 3d0 * kbt / 2d0*dfloat(nf)

end subroutine

end module
