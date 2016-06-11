program van_der_waals


use temperature


implicit none

include 'mpif.h'

! sigma: sigma parameter of lennard jones potencial
! epsil: epsilon parameter of lennard jones potencial

! parametres del mpi
 integer stat(MPI_STATUS_SIZE)
 integer comm, taskid, numproc, ierror, aerr, partner
 integer MASTER,message,sender
 integer, allocatable, dimension(:) :: request

! system parameters
integer :: M, nf
double precision :: boxL, mass, kBTref
double precision :: sigma, epsil, utime
! Integration parameters
integer :: i, j, N
double precision :: dt
! General variables
double precision, dimension(:,:), allocatable :: r, v, F,F_t
! Ouputs
double precision :: ecin, epot, temp




! defining parameters
M = 300 
N = 100
mass = 1.
sigma=1
epsil=1

! Defining dimensions
allocate(r(3,M), v(3,M), F(3,M),F_t(3,M))

! Init variables
kBTref = 2.
boxL=10d0 / sigma
utime = dsqrt(mass * sigma**2 / epsil) ! unit of time in LJ units
dt = 1. / 400 / utime
nf = 3 * M - 3 ! number of degrees of freedom 

do j=1,M
      v(1,j)=1
      v(2,j)=1
      v(3,j)=1
enddo

   call MPI_INIT(ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierror)
   
   call compute_temperature(temp, v, M, nf)
   call kinetic_energy(ecin, temp, M, nf)
   



    if (taskid .eq. MASTER) then 
      write (*,*) "T, nf, ecin"
      write (*,*) temp,nf, ecin
    endif

   call MPI_FINALIZE(ierror);


end program
