program van_der_waals
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
double precision :: ecin, epot, temp
! Integration parameters
integer :: i, j, N
double precision :: dt
! General variables
double precision, dimension(:,:), allocatable :: r, v, F,F_t
! Ouputs
double precision :: temp
! defining parameters
parameter(M = 400, N = 100)
parameter(mass = 4., sigma=1,epsil=1)

! Defining dimensions
allocate(r(3,M), v(3,M), F(3,M),F_t(3,M))

! Init variables
kBTref = 2.
dt = 1. / 500
nf = 3 * M - 3 ! number of degrees of freedom 

! Initial configuration+velocity

call inirandom(M,r,v,boxL,kBTref)

! Open files
open(unit=2, file='data/energy-temp.data', status='unknown')
! open(unit=3, file='../data/RDF.data', status='unknown')

call forces(M,r,F,boxL,sigma,epsil)

! Temporal loop
do i = 1, N

   call Verlet_Coord(r,F,v,M,dt,mass)

   ! Apply PBC

   call PBC(M,r,boxL)

   ! Calculate new Forces

   call forces(M,r,F_t,boxL,sigma,epsil)

   call Verlet_Vel(F,F_t,v,M,dt,mass)

   F = F_t

   ! Compute magnitudes

   call compute_temperature(temp, v, M, nf)

   ! Save temporal serie 
   write(*,*) i*dt, temp

   ! Compute RDF and average

end do

! Save RDF

close(2)
close(3)

end program
