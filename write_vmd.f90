module vmd
implicit none

contains

subroutine output(M,i,dt,r)
!Subrutina que me guarda las coordenadas para que las pueda ver con el vmd
integer, intent(in) :: M
double precision, dimension(3,M), intent(in) :: r
double precision, intent(in) :: dt
integer :: i, j

write(20,*) M
write(20,*) "trajectoria,iter:",i,"temps",i*dt,"ns"

do j=1,M
   write(20,*) "C ", r(1,j), r(2,j), r(3,j)
enddo

end subroutine

end module

