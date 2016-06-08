module PeriodicBoundaryConditions
implicit none
public

!Modulo para implementar condiciones periodicas de contorno
!en la dinamica molecular

contains

subroutine PBC(num_particulas,posicion,Lcaja,numproc,rank)
!Necesita de input el numero de particulas, la matriz con las coordenadas,
!y el lado de la caja de simulacion.
!Devuelve como output las matriz con las coordenadas de las particulas 
!o su imagen dentro de la caja

integer, intent(in) :: num_particulas, numproc, rank
double precision, dimension(3,num_particulas), intent(inout) :: posicion
double precision, intent(in) :: Lcaja
integer :: i,j,partxproc


partxproc=nint (nump_particules/numproc)

If (rank=numproc-1) then
    
   do j=rank*partxproc+1,num_particulas
       do i=1,3
            If (posicion(i,j)>Lcaja/2.) then
                posicion(i,j)=posicion(i,j)-Lcaja
            else if (posicion(i,j)<-Lcaja/2.) then
                posicion(i,j)=posicion(i,j)+Lcaja
            endif
    
       enddo
    enddo



else

   do j=rank*partxproc+1,(rank+1)*partxproc
       do i=1,3
            If (posicion(i,j)>Lcaja/2.) then
                posicion(i,j)=posicion(i,j)-Lcaja
            else if (posicion(i,j)<-Lcaja/2.) then
                posicion(i,j)=posicion(i,j)+Lcaja
            endif
    
       enddo
    enddo

end if

call MPI_BARRIER(rank, ierror)
endsubroutine

end module
