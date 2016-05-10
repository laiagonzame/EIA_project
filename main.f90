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

! system parameters
integer :: M, nf
double precision :: boxL, mass, kBTref
double precision :: ecin, epot, temp
! Integration parameters
integer :: i, j, N
double precision :: dt
! General variables
double precision, dimension(:,:), allocatable :: r, v, accel
! defining parameters
parameter(M = 400, N = 100)
parameter(mass = 4.)

! Defining dimensions
allocate(r(3,M), v(3,M), accel(3,M), rho(P))

! Init variables
kBTref = 2.
dt = 1. / 500
nf = 3 * M - 3 ! number of degrees of freedom 

! Initial configuration

! Open files
open(unit=2, file='../data/energy-temp.data', status='unknown')
open(unit=3, file='../data/RDF.data', status='unknown')

! Temporal loop
do i = 1, N
   ! Update forces

   ! Update r and v

   ! Compute magnitdes

   ! Save temporal serie 

   ! Compute RDF and average

end do

! Save RDF

close(2)
close(3)

end program
