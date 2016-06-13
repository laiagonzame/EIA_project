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


subroutine inicubic(N,pos,vel,box,T)
        ! Subrutina que inicializa las particulas en la caja en
        ! una red cubica  
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
                call random_number(aux)
                vel(1,cont)= (2.*aux-1.)*dsqrt(T)
                call random_number(aux)
                vel(2,cont)= (2.*aux-1.)*dsqrt(T)
                call random_number(aux)
                vel(3,cont)= (2.*aux-1.)*dsqrt(T)
                cont=cont+1
        end if
        enddo

        pos_cm = sum(pos, dim=2) / N
        vel_cm = sum(vel, dim=2) / N
        do i = 1, N
                pos(:,i) = pos(:,i) - pos_cm
                vel(:,i) = vel(:,i) - vel_cm
        end do
end subroutine

!**************************************************************
!**************************************************************

subroutine inirandom(N,pos,vel,box,T)
        ! Subrutina que inicializa las particulas en la caja de manera
        ! aleatoria en posicion y velocidad
        integer :: N,i,j,seed,cont,x,y,lat,z
        double precision, dimension(3,N) :: pos, vel
        double precision :: box,aux,T,r,dx,dy,dz
        integer, dimension(:,:,:), allocatable :: lattice
        double precision :: valx,valy,valz
        double precision, dimension(3) :: pos_cm, vel_cm


! Inicializacion de los numeros aleatorios
        seed=1221
        call random_seed(seed)
        cont=1

       do while (cont<N+1)
              
 10           call random_number(aux)
              pos(1,cont)=(2*aux-1.)*box
              call random_number(aux)
              pos(2,cont)=(2*aux-1.)*box
              call random_number(aux)
              pos(3,cont)=(2*aux-1.)*box
              
              vel(1,cont)=0.
              vel(2,cont)=0.
              vel(3,cont)=0.

              
              do i=1,cont-1
                      
                  dx=pos(1,i)-pos(1,cont)
                  dy=pos(2,i)-pos(2,cont)
                  dz=pos(3,i)-pos(3,cont)
                  dx=dx-box*nint(dx/box)
                  dy=dy-box*nint(dy/box)
                  dz=dz-box*nint(dz/box)
                  r=sqrt(dx**2+dy**2+dz**2)
                  if(r<1.4) then
                  go to 10
                  end if
                  
               end do

            
                cont=cont+1
               
         
                
                
        
end do
end subroutine
end module




