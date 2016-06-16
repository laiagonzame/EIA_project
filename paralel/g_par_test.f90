!Programa para testear g(r)

program main
use statistics
implicit none
  include 'mpif.h'

  real,dimension(:,:),allocatable :: r
  real :: L, dr, x, densidad
  real,dimension(:),allocatable :: g
  integer :: M, total_esp, i, j, count
  integer,dimension (1) :: seed 
  real,dimension(1) :: rnd

    call system_clock(count)

    seed = count

    call random_seed(put=seed)

  M = 10000
  L = 100
  densidad = M/L**3
  total_esp = 200
  
  allocate(r(3,M))
  allocate(g(total_esp))
  
  do i=1,M
     do j=1,3
        call random_number(rnd)
        x = rnd(1)
        r(j,i) = -L/2. + x*L
     end do
  end do

  open(unit=7,file="./pos.xyz")
  do i=1,M
     write(7,*) r(:,i)
  end do
  close(7)

  call declarate_radial_dist(L,total_esp,dr,g)
  call accumulate_radial_dist(M,L,r,total_esp,dr,g)
  call compute_radial_dist(total_esp, dr, 1, densidad, M, g)
  
  open(unit=8,file="./g_test.out")
  do i=1,total_esp
    write(8,*) (i-1)*dr, g(i)
  end do
  close(8)
end program
