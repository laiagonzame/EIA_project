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
integer :: i,Mtask,MASTER,ierror,taskid,numproc,request
 integer stat(MPI_STATUS_SIZE)
double precision message

!--------- MPI inicialization ------------------
call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierror)

 MASTER = 0
 ! dividing work
Mtask=M/numproc
 if (taskid .eq. MASTER) then 
   Mtask=Mtask+mod(M,numproc)
 endif
! print *,"Hello from task ", taskid," of ", numproc, " Mtask=",Mtask," of ",M
!---------------------------

!--------- do their part of the job ------------------
message=0d0
do i=1,Mtask
   message = message + (vel(1,i)**2+vel(2,i)**2+vel(3,i)**2)/dfloat(nf)
end do

!--------- Barrier: waiting for everyone ------------------
 call MPI_BARRIER(MPI_COMM_WORLD, ierror)
 if (taskid .eq. MASTER) then
!   print *, "   ==> END BARRIER TASK <== "
 endif
 
 !--------- send information to MASTEr ------------------
    request=0
    call MPI_ISEND(message, 1, MPI_INTEGER, MASTER, 1, MPI_COMM_WORLD, request,ierror)
!    print *, "Sender",taskid," message  ",message 
    
 if (taskid .eq. MASTER) then 
   kbt=0d0
   do i=0,numproc-1
      call MPI_RECV(message, 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD, stat, ierror)
      kbt=kbt+message
   enddo
!   print *, "MASTER got kbt=",kbt
 endif




end subroutine


subroutine kinetic_energy(kine, kbt, M, nf)
integer, intent(in) :: M, nf
double precision, intent(in) :: kbt
double precision, intent(out) :: kine
 
kine = 3d0 * kbt / 2d0*dfloat(nf)

end subroutine

end module