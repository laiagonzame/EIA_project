program van_der_waals

use initial
use PeriodicBoundaryConditions
use integrator_para
use forces_routines
use vmd

implicit none

include 'mpif.h'



! N: number of time steps
! M: number of atoms
! boxL: boxlength
! mass: particle mass
! rho: density
! kBTref: reference temperature
! nf: degrees of freedom
! ecin: kinetic energy
! epot: potential energy
! temp: temperature
! F: Force Matrix
! sigma: sigma parameter of lennard jones potencial
! epsil: epsilon parameter of lennard jones potencial
! nhis: number of edges on g(r) histogram
! tterm: temps que tarda en termalitzar. A partir d'aquí guardem trajs
! stepwrite: cada quants steps de temps escrivim

! system parameters
integer :: M, nf, nhis, stepwrite, stepwrite_count
double precision :: boxL, mass, kBTref, tterm
double precision :: sigma, epsil, utime
! Integration parameters
integer :: i, j, N
double precision :: dt
! General variables
double precision, dimension(:,:), allocatable :: r, real_r, v, F,F_t
! Ouputs
double precision :: ecin, epot, temp
! defining parameters
parameter(M = 300, N = 1000, nhis = 400)
parameter(mass = 1., sigma=1, epsil=1)
parameter(tterm = 0, stepwrite=1)
! variables de paralelizacio
integer :: stat(MPI_STATUS_SIZE)
integer ::  ierror, request, rank, numproc,MASTER

call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierror)
MASTER=0


! Defining dimensions
allocate(r(3,M),real_r(3,M), v(3,M), F(3,M),F_t(3,M))

! Init variables
kBTref = 2.
boxL=10d0 / sigma
utime = dsqrt(mass * sigma**2 / epsil) ! unit of time in LJ units
dt = 1. / 300 / utime
nf = 3 * M - 3 ! number of degrees of freedom 
stepwrite_count = 0 ! number of snapshots

! Open files
if (rank .eq. MASTER) then 
   open(unit=10, file='data/params.data', status='unknown')
   open(unit=11, file='data/posicions.data', status='unknown')
   open(unit=12, file='data/velocitats.data', status='unknown')
   open(unit=13, file='data/posicions_reals.data', status='unknown')
   open(unit=14, file='data/ener_potencial.data', status='unknown')
   open(unit=20, file='data/traj_vmd.data', status='unknown')
end if

! save parameters
write(10,*) boxL, nhis, M, sigma, epsil, mass, dt, kBTref, tterm, stepwrite

! Initial configuration+velocity

call inirandom(M,r,v,boxL,kBTref,numproc,rank)
real_r = r
call PBC(M,r,boxL,numproc,rank)
call forces(M,r,F,boxL,sigma,epsil,epot,rank,numproc,MASTER)
! Temporal loop
do i = 1, N
   
   call Verlet_Coord(r,F,v,M,dt,mass,boxL,rank,numproc)
   call Verlet_Coord(real_r,F,v,M,dt,mass,boxL,rank,numproc)

   ! Apply PBC

   call PBC(M,r,boxL,numproc,rank)
  
   ! Calculate new Forces

   call forces(M,r,F_t,boxL,sigma,epsil,epot,rank,numproc,MASTER)
   call Verlet_Vel(F,F_t,v,M,dt,mass,rank,numproc)
   
   !write trajectories and potencial energy
   if (rank .eq. MASTER) then 
      If (i*dt > tterm .AND. mod(i,stepwrite) ==  0) then
         stepwrite_count = stepwrite_count + 1
         call output(M,i,dt,r)
         do j = 1, M
              write(11,*) r(1,j), r(2,j), r(3,j)
              write(12,*) v(1,j), v(2,j), v(3,j)
              write(13,*) real_r(i,j), real_r(2,j), real_r(3,j)
         end do
         write(14,*) epot
      endif
   endif

   F = F_t
   if (rank == MASTER) then
   print *, i
   end if
end do

! save number of snapshots
if (rank .eq. MASTER) then 
   write(10,*) stepwrite_count
end if

close(10)
close(11)
close(12)
close(13)
close(14)
close(20)

call MPI_FINALIZE(ierror)

end program
