program van_der_waals

use temperature
use initial
use PeriodicBoundaryConditions
use integrator
use forces_routines
use statistics

implicit none

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
parameter(M = 300, N = 100)
parameter(mass = 1., sigma=1, epsil=1)

! Defining dimensions
allocate(r(3,M), v(3,M), F(3,M),F_t(3,M))

! Init variables
kBTref = 2.
boxL=10d0 / sigma
utime = dsqrt(mass * sigma**2 / epsil) ! unit of time in LJ units
dt = 1. / 400 / utime
nf = 3 * M - 3 ! number of degrees of freedom 

do j=1,M
      v(1,j)=300
      v(2,j)=0
      v(3,j)=0
      write (*,*) v(1,j),v(2,j),v(3,j)
enddo

   call compute_temperature(temp, v, M, nf)
   
   write (*,*) "T, nf"
   write (*,*) temp,nf


end program
