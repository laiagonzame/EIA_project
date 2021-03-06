module forces_routines
implicit none

contains

subroutine forces(M,r,F,boxlength,sigma,epsil,epot,taskid,numproc,MASTER)
include 'mpif.h'
integer, intent(in) :: M, taskid, numproc, MASTER
double precision, dimension(3,M), intent(in) :: r
double precision, dimension(3,M), intent(out) :: F
double precision, intent(in) :: boxlength, sigma, epsil
double precision, dimension(3) :: rij 
double precision :: rcut, rijl, dist2, sigmar, rc2, force,epot
integer :: i, j, l
! parallelization variables
integer :: ierror, extras, pairs, counter, request
integer :: stat(MPI_STATUS_SIZE)
integer, dimension(:), allocatable :: ppw, indx_ppw ! pairs per worker
integer, dimension(:,:), allocatable :: vec_pairs
double precision, dimension (:,:), allocatable :: Fij
double precision, dimension (:), allocatable :: epot_vec

F = 0d0
epot = 0d0
rc2 = 25d0 * sigma * sigma

allocate(ppw(numproc))

! get number of pairs
pairs = M * (M-1) / 2
allocate(vec_pairs(2,pairs), Fij(3,pairs), indx_ppw(0:numproc), epot_vec(pairs))
epot = 0d0
Fij = 0d0

! pairs per worker
ppw = pairs / numproc
extras = mod(pairs, numproc)
do i = 1, extras
   ppw(i) = ppw(i) + 1
end do

! index pairs per worker
indx_ppw(0) = 1
do i = 1, numproc
   indx_ppw(i) = indx_ppw(i-1) + ppw(i)
end do

! The vector that contains all pairs interactions is computed
counter = 1
do i = 1, M-1
   do j = i+1, M
      vec_pairs(:,counter) = (/ i, j /)
      counter = counter + 1
   end do
end do

! Each worker computes the force and potencial energy due to its pairs
do i = indx_ppw(taskid), indx_ppw(taskid+1)-1
   call lj(r(:,vec_pairs(1,i)), r(:,vec_pairs(2,i)), boxlength, sigma, epsil, rc2, Fij(:,i), epot_vec(i))
end do

! ----- WORKER(i) to MASTER -----
if (taskid .ne. MASTER) then 
   call MPI_ISEND(Fij(:,indx_ppw(taskid):indx_ppw(taskid+1)-1), 3*ppw(taskid), MPI_REAL8, MASTER, 1, MPI_COMM_WORLD, request, ierror)
   call MPI_ISEND(epot_vec(indx_ppw(taskid):indx_ppw(taskid+1)-1), ppw(taskid), MPI_REAL8, MASTER, 1, MPI_COMM_WORLD, request, ierror)
end if


! ---- waiting all WORKERS -----
call MPI_BARRIER(MPI_COMM_WORLD, ierror)

! ---- MASTER recive and merge
if (taskid .eq. MASTER) then
   do i = 1, numproc-1
      call MPI_RECV(Fij(:,indx_ppw(i):indx_ppw(i+1)-1), 3*ppw(i), MPI_REAL8, i, 1, MPI_COMM_WORLD, stat, ierror)
      call MPI_RECV(epot_vec(indx_ppw(i):indx_ppw(i+1)-1), ppw(i), MPI_REAL8, i, 1, MPI_COMM_WORLD, stat, ierror)
   end do
   
   do i = 1, pairs
      F(:, vec_pairs(1,i)) = F(:, vec_pairs(1,i)) + Fij(:,i) 
      F(:, vec_pairs(2,i)) = F(:, vec_pairs(2,i)) - Fij(:,i) 
   end do
   epot = sum(epot_vec)

   ! ----- MASTER to WORKERS -----
   do i = 1, numproc-1
      call MPI_ISEND(F, 3*M, MPI_REAL8, i, 1, MPI_COMM_WORLD, request, ierror)
      call MPI_ISEND(epot, 1, MPI_REAL8, i, 1, MPI_COMM_WORLD, request, ierror)
   end do
end if

call MPI_BARRIER(MPI_COMM_WORLD, ierror)

! ---- all WORKERS recive ---
if (taskid .ne. MASTER) then
   call MPI_RECV(F, 3*M, MPI_REAL8, MASTER, 1, MPI_COMM_WORLD, stat, ierror)
   call MPI_RECV(epot, 1, MPI_REAL8, MASTER, 1, MPI_COMM_WORLD, stat, ierror)
end if

call MPI_BARRIER(MPI_COMM_WORLD, ierror)

endsubroutine

subroutine lj(ri, rj, boxlength, sigma, epsil, rc2, F, epot)
double precision, dimension(3), intent(in) :: ri, rj
double precision, intent(in) :: boxlength, sigma, epsil, rc2
double precision, dimension(3), intent(out) :: F
double precision, intent(out) :: epot
double precision, dimension(3) :: rij 
double precision :: rijl, dist2, sigmar, force
integer :: l

! Initialize the distance, force and potential energy
dist2=0d0
F = 0d0
epot = 0d0

! Compute distance
do l=1,3
   rijl = rj(l) - ri(l)
   rij(l) = rijl - boxlength*nint(rijl/boxlength)
   dist2 = dist2 + rij(l)*rij(l)
enddo
! Ensure that distance is smaller than rc
if (dist2<rc2) then
   sigmar = (dist2**(-3d0))*(sigma**6d0)
   ! Compute force modulus over distance
   force = 24d0*epsil*(2d0*sigmar**(2d0)-sigmar)/dist2
   ! Update force
   do l=1,3
      F(l) = -force*rij(l)
   enddo
   ! Update potential energy
   epot=4d0*epsil*(sigmar**(2.)-sigmar)
endif

endsubroutine
end module
