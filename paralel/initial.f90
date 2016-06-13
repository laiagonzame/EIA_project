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


subroutine inirandom(N,pos,vel,box,T,numproc,rank)
        ! Subrutina que inicializa las particulas en la caja de manera
        ! aleatoria en posicion y velocidad

        include 'mpif.h'

        integer :: N,i,j,seed,cont,x,y,lat,z
        double precision, dimension(3,N) :: pos, vel
        double precision :: box,aux,T,r,dx,dy,dz
        integer, dimension(:,:,:), allocatable :: lattice
        double precision :: valx,valy,valz
        double precision, dimension(3) :: pos_cm, vel_cm
        integer :: stat(MPI_STATUS_SIZE)
        integer :: request,ierror
        integer :: Np,numproc,rank,nini,nfin,MASTER


        

! Inicializacion de los numeros aleatorios
        seed=1221
        call random_seed(seed)

        Np=floor(float(N)/float(numproc))

        MASTER=0

        if (rank .eq.(numproc-1)) then
          nini=rank*Np+1
          nfin=N
        else
          nini=rank*Np+1
          nfin=(rank+1)*Np
        end if
        
        cont=nini
      do while (cont<nfin+1)
              
 10           call random_number(aux)
              pos(1,cont)=-box/2.+(0.1+0.8*aux+float(rank))*box/float(numproc)
              call random_number(aux)
              pos(2,cont)=(2*aux-1.)*box
              call random_number(aux)
              pos(3,cont)=(2*aux-1.)*box
              
              vel(1,cont)=0.
              vel(2,cont)=0.
              vel(3,cont)=0.

              
              do i=nini,cont-1
                      
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

        if (rank .ne. 0) then
          call MPI_ISEND(pos(:,nini:nfin),1, MPI_REAL, 0, 1, MPI_COMM_WORLD, request, ierror)
        endif
        
        call MPI_BARRIER(MPI_COMM_WORLD, ierror) !Esperamos que todos hayan mandado

        if (rank .eq. 0) then
         do i =1, numproc-1
          call MPI_RECV(pos(:,nini:nfin), 1, MPI_REAL, i, 1, MPI_COMM_WORLD, stat, ierror)
         enddo

         do i=1, numproc-1
          call MPI_ISEND(pos,1,MPI_REAL,i,1, MPI_COMM_WORLD,request,ierror)
         enddo

        
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)


        if (rank .ne. 0) then
           call MPI_RECV(pos,1, MPI_REAL, 0, 1, MPI_COMM_WORLD, stat, ierror)
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, ierror)

        endif
end subroutine
end module
