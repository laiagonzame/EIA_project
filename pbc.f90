module pbc
implicit none
public

!Modulo para implementar condiciones periodicas de contorno
!en la dinamica molecular



contains

subroutine PBC(num_particulas,posicion,Lcaja)
!Necesita de input el numero de particulas, la matriz con las coordenadas,
!y el lado de la caja de simulacion.
!Devuelve como output las matriz con las coordenadas de las particulas 
!o su imagen dentro de la caja

integer,intent(in)::num_particulas
real(8),dimension(num_particulas,3),intent(inout)::posicion
integer,intent(in)::Lcaja
integer::i,j


do i=1,num_particulas
  do j=1,3
    If (posicion(i,j)>0) then
        If (posicioni(i,j)>Lcaja/2.) then
           posicion(i,j)=posicion(i,j)-Lcaja
        endif
    else
       If (posicion(i,j)<-Lcaja/2.) then
          posicion(i,j)=posicion(i,j)+Lcaja
       endif
    endif
  enddo
enddo

endsubroutine
