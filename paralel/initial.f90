module initial
implicit none
public
!Modulo que inicializa el sistema a partir de una densidad dada
!mediante el uso de la paralelizacion de procesadores
contains

!N -- # de particulas
!pos -- matriz de posiciones
!vel -- matriz de velocidades
!box -- Longitud de la caja
!T -- temperatura


!**************************************************************


subroutine inirandom(N,pos,vel,box,T,numproc,rank)
        ! Subrutina que inicializa las particulas en la caja de manera
        ! aleatoria en posicion y velocidad utilizando varios procesadores

        include 'mpif.h'

        integer :: N,i,j,seed,cont,x,y,lat,z,aux2
        double precision, dimension(3,N) :: pos, vel
        double precision :: box,aux,T,r,dx,dy,dz
        integer, dimension(:,:,:), allocatable :: lattice
        double precision :: valx,valy,valz
        double precision, dimension(3) :: pos_cm, vel_cm
        integer :: stat(MPI_STATUS_SIZE)
        integer :: request,ierror
        integer :: Np,numproc,rank,nini,nfin,MASTER


        

! Inicializacion de los numeros aleatorios
        seed=11*rank+3
        call random_seed(seed)
        do i=1,N*rank
        call random_number(aux)
	end do
	
        call random_number(aux)
        Np=floor(float(N)/float(numproc)) !Numero de particulas para cada procesador

	!Definimos que particulas tiene que introducir el procesador
        if (rank .eq.(numproc-1)) then 
          nini=rank*Np+1
          nfin=N
        else
          !El ultimo se queda las particulas sobrantes (problema de dividir dos integers)
          nini=rank*Np+1
          nfin=(rank+1)*Np
        end if
        aux=3
        
        !Introducimos las particulas
        cont=nini
      do while (cont<nfin+1)
              !Tenemos dos dimensiones aleatorioas completamente, y una restringida a la seccion
              !del procesador
 10           call random_number(aux)
              pos(1,cont)=-box/2.+(0.1+0.8*aux+float(rank))*box/float(numproc)
              call random_number(aux)
              pos(2,cont)=(2*aux-1.)*box/2.
              call random_number(aux)
              pos(3,cont)=(2*aux-1.)*box/2.
           
              do i=nini,cont-1
                      
                  dx=pos(1,i)-pos(1,cont)
                  dy=pos(2,i)-pos(2,cont)
                  dz=pos(3,i)-pos(3,cont)
                  dx=dx-box*nint(dx/box)
                  dy=dy-box*nint(dy/box)
                  dz=dz-box*nint(dz/box)
                  r=sqrt(dx**2+dy**2+dz**2)
                  if(r<1.0) then !Condicion de no solapamiento
                  go to 10
                  end if
                  
               end do

                cont=cont+1
      
end do

!Enviamos al master la informacion de cada procesador

        if (rank .ne. 0) then
          call MPI_ISEND(pos(:,nini:nfin),(nfin-nini+1)*3, MPI_REAL8, 0, 1, MPI_COMM_WORLD, request, ierror)
        endif
         
 !Esperamos a que acaben de enviar todos
                
        call MPI_BARRIER(MPI_COMM_WORLD, ierror) !Esperamos que todos hayan mandado
        
        if (rank .eq. 0) then !Si somos el master, tenemos que recibir y reenviar la informacion
          do i =1, numproc-1
             if (i .eq.(numproc-1)) then
                nini=i*Np+1
                nfin=N
             else
          	nini=i*Np+1
          	nfin=(i+1)*Np
             end if
	!Recibimos la informacion correspondiente
           call MPI_RECV(pos(:,nini:nfin), (nfin-nini+1)*3, MPI_REAL8, i, 1, MPI_COMM_WORLD, stat, ierror)
         enddo
         
         !Reenviamos la matriz completa para que todos puedan trabajar con ella
         do i=1, numproc-1
          call MPI_ISEND(pos,N*3,MPI_REAL8,i,1, MPI_COMM_WORLD,request,ierror)
         enddo

        end if
        
        !Esperamos a que se haya enviado todo
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        
             !Si somos diferentes al master, recibimos las matriz completa
        if (rank .ne. 0) then
           call MPI_RECV(pos,N*3, MPI_REAL8, 0, 1, MPI_COMM_WORLD, stat, ierror)
        end if

	!Esperamos a que todos lleguen y ya podemos continuar con las siguientes partes del codigo
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        
         if (rank .eq. 0) then !Si somos el master, tenemos que recibir y reenviar la informacion
            vel=0.0
         
         !Reenviamos la matriz completa para que todos puedan trabajar con ella
         do i=1, numproc-1
          call MPI_ISEND(vel,N*3,MPI_REAL8,i,1, MPI_COMM_WORLD,request,ierror)
         enddo

        end if
        
        !Esperamos a que se haya enviado todo
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        
             !Si somos diferentes al master, recibimos las matriz completa
        if (rank .ne. 0) then
           call MPI_RECV(vel,N*3, MPI_REAL8, 0, 1, MPI_COMM_WORLD, stat, ierror)
        end if

	!Esperamos a que todos lleguen y ya podemos continuar con las siguientes partes del codigo
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        
end subroutine
end module
