module initial
implicit none
public

!Modulo que inicializa el sistema a partir de una densidad dada


contains

!N -- # de particulas
!pos -- matriz de posiciones
!vel -- matriz de velocidades
!box -- Longitud de la caja
!T -- temperatura 

!**************************************************************
subroutine inirandom(N,pos,vel,box,T)
! Subrutina que inicializa las particulas en la caja de manera
! aleatoria en posicion y velocidad

integer :: N,i,j,seed,cont,x,y,lat,z
double precision, dimension(3,N) :: pos, vel
double precision :: box,aux,T
integer, dimension(:,:,:), allocatable :: lattice
double precision :: valx,valy,valz
double precision, dimension(3) :: pos_cm, vel_cm

lat=floor(N**(1./3.)+1)
allocate(lattice(lat,lat,lat))
lattice = 0
cont = 1

! Inicializacion de los numeros aleatorios
seed=1221
call random_seed(seed)

do while (cont<N+1)

      call random_number(valx)
      call random_number(valy)
      call random_number(valz)
      x=floor(valx*lat+1)
      y=floor(valy*lat+1)
      z=floor(valz*lat+1)
      if (lattice(x,y,z).eq.0) then
      lattice(x,y,z)=1
      pos(:,cont) = (/x, y, z/) * box / dfloat(lat)
      cont=cont+1
      end if
      ! Posiciones
      
      
      ! Velocidades
      call random_number(aux)
      vel(:,cont)= (2.*aux-1.)*dsqrt(T) 
   
enddo
pos_cm = sum(pos, dim=2) / N
vel_cm = sum(vel, dim=2) / N

do i = 1, N
   pos(:,i) = pos(:,i) - pos_cm
   vel(:,i) = vel(:,i) - vel_cm
end do

end subroutine

end module
