program van_der_waals

use temperature
use initial
use PeriodicBoundaryConditions
use integrator
use forces_routines
use statistics

implicit none
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
parameter(M = 300, N = 500)
parameter(mass = 1., sigma=1, epsil=1)

! Defining dimensions
allocate(r(3,M), v(3,M), F(3,M),F_t(3,M))

! Init variables
kBTref = 2.
boxL=10d0 / sigma
utime = dsqrt(mass * sigma**2 / epsil) ! unit of time in LJ units
dt = 1. / 400 / utime
nf = 3 * M - 3 ! number of degrees of freedom 

! Initial configuration+velocity

call inirandom(M,r,v,boxL,kBTref)

! Open files
open(unit=2, file='data/energy.data', status='unknown')
open(unit=3, file='data/temp.data', status='unknown')
! open(unit=4, file='../data/RDF.data', status='unknown')

call forces(M,r,F,boxL,sigma,epsil,epot)

! Temporal loop
do i = 1, N

   call Verlet_Coord(r,F,v,M,dt,mass)

   ! Apply PBC

   call PBC(M,r,boxL)

   ! Calculate new Forces

   call forces(M,r,F_t,boxL,sigma,epsil,epot)

   call Verlet_Vel(F,F_t,v,M,dt,mass)

   F = F_t

   ! Compute magnitudes

   call compute_temperature(temp, v, M, nf)
   call kinetic_energy(ecin, temp)

   ! Save temporal serie 
   write(2,*) i*dt, temp
   write(3,*) i*dt, ecin, epot 

   ! Compute RDF and average

end do

! Save RDF

close(2)
close(3)
!close(4)

end program
