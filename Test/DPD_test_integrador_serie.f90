program main
use integrator
INTEGER SEED(4)
REAL*8 :: U1, U2, V1, V2, X1, X2
REAL*8:: RANNYU, mu, w, av, y, var
real*8 :: L,ts, t,temp, tmax,  rx_pbc, ry_pbc, rz_pbc, r_pbc,fx,fy,fz,dens,temp_calc,Fc_a,vx,vy,vz,C_VACF 
real*8, dimension(:,:), allocatable :: coord, vel, vel_i, frz, frzt
integer :: iter=1,iseed,Num,i, cont2=0,N
real*8, parameter :: Kb=1, sigma=1, m=1, pi_c=3.14159265359,sig=5


! Programa que realiza una simulación de Dinamica Molecular de particulas disipativas realitzado y testeado previamente
! Se ha añadido el modulo de la integracion velocity verlet en serie para comprobar que da los mismos resultados
! que en el programa original
! Programa realizado por Pablo M. Blanco Andres, Universidad de Barcelona.




! Inicializacion del sistema, lectura de variables, preparacion de archivos de output
 



call lectura_input()



allocate(coord(Num,3),vel(Num,3),  frz(Num,3), frzt(Num,3))



if (iseed.eq.0) then
         call setseed(iseed)
      end if
      call setr1279(iseed)



call inic_sist(coord,vel)


open (UNIT=15,file="Temp.txt")
write(15,*) "#temps         temperatura"
open (UNIT=18,file="VACF.txt")
write(18,*) "#temps                   VACF"



! Control de que la velocidad del centro de massas sea 0 inicialmente

call check(vel,Num)
vel_i=vel
call  VACF(vel_i,vel,C_VACF)
write(18,*)t,C_VACF

frz=0.

! Bucle principal de la dinamica

  do while (t < tmax) 
  
! Para el primer tiempo, se calculan las fuerzas entre las particulas
! en la posicion inicial  

  if ( t == 0) then 
  frz=0.

! bucle del calculo de fuerzas entre particulas

   do i=1,Num

    do j=i,Num

      if ( i /= j) then

   ! Se calcula la distancia entre particulas teniendo en cuenta las PBC

      call dist_pbc(coord(i,1),coord(i,2),coord(i,3),coord(j,1),coord(j,2),coord(j,3),rx_pbc,ry_pbc,rz_pbc)
      r_pbc=sqrt(rx_pbc**2+ry_pbc**2+rz_pbc**2)
      if (r_pbc < 1) then
          vx=vel(i,1)-vel(j,1)
          vy=vel(i,2)-vel(j,2)
          vz=vel(i,3)-vel(j,3)

          call Force(rx_pbc,ry_pbc,rz_pbc,r_pbc,vx,vy,vz,fx,fy,fz)
 
! Asignacion de las fuerzas entre las particulas i i j (son iguales y de sentido contrario)


          frz(i,1)=frz(i,1)+fx
          frz(i,2)=frz(i,2)+fy
          frz(i,3)=frz(i,3)+fz
           
          frz(j,1)=frz(j,1)-fx
          frz(j,2)=frz(j,2)-fy
          frz(j,3)=frz(j,3)-fz
        end if

      end if
     

    end do
 
  end do


! Calculo de la temperatura inicial

  call  Temperature(vel,temp_calc)
  write(15,*)t,temp_calc 

  end if

! Calculo de la nueva posicion de las particulas usando el modulo integrator en serie y aplicacion de las PBC

  call Verlet_Coord(coord,frz,vel,ts,m,L)

    
  call Coord_PBC(coord)
   
! Escritura de la posicio de les particules al temps t

  call out_coord
  call Vel_Corr(frz,vel)

! Calculo de la nueva fuerza entre particulas en la nueva posicion

  frzt=0.


   do i=1,Num

    do j=i,Num

      if ( i /= j) then

 ! Se calcula la distancia entre particulas teniendo en cuenta las PBC

      call dist_pbc(coord(i,1),coord(i,2),coord(i,3),coord(j,1),coord(j,2),coord(j,3),rx_pbc,ry_pbc,rz_pbc)
      r_pbc=sqrt(rx_pbc**2+ry_pbc**2+rz_pbc**2)
      if (r_pbc < 1) then
          vx=vel(i,1)-vel(j,1)
          vy=vel(i,2)-vel(j,2)
          vz=vel(i,3)-vel(j,3)



          call Force(rx_pbc,ry_pbc,rz_pbc,r_pbc,vx,vy,vz,fx,fy,fz)


! Asignacion de las fuerzas entre las particulas i i j (son iguales y de sentido contrario)

           frzt(i,1)=frzt(i,1)+fx
           frzt(i,2)=frzt(i,2)+fy
           frzt(i,3)=frzt(i,3)+fz
           frzt(j,1)=frzt(j,1)-fx
           frzt(j,2)=frzt(j,2)-fy
           frzt(j,3)=frzt(j,3)-fz

         end if

      end if


    end do

   end do




! Calculo de la velocidad a tiempo t (usando el modulo integrator en serie)

  call Verlet_Vel(frz,frzt,vel,N,ts,m)
  call New_vel(frz,frzt,vel)

! Calculo de la temperatura a tiempo t

 call  Temperature(vel,temp_calc)

 write(15,*)t/60.,temp_calc

   call  VACF(vel_i,vel,C_VACF)
   write(18,*)t,C_VACF

t=t+ts


 frz=frzt

 iter=iter+1 

! Fin del bucle temporal, comprobacion de que la velocidad del centro de massas se conserva
! Esta comprobacion solo se realiza una vez cada 100 pasos de tiempo.

if (cont2 == 100) then

 call check(vel,Num)

cont2=0

end if


  
cont2=cont2+1


  end do
  


 call radial_distribution(coord)



contains 

subroutine lectura_input()
character(20) :: trash

! Subrutina que lee las variables necesarias para la simulacion 
! (numero de particulas, densidad, paso de tiempo, tiempo de simulaion, temperatura, 
! radio de cut-off y semilla para los generadores de numeros aleatorios)
! del archivo de input "Input.txt"


open (unit=12, file='Input.txt')

read (12,*), trash, Num
read (12,*), trash, dens
read (12,*), trash, ts
read (12,*), trash, tmax
read (12,*), trash, temp
read (12,*), trash, iseed
read (12,*), trash, Fc_a


L=(Num*1.0/dens)**(1./3.)
close(12)

call vmd_script


end subroutine

subroutine inic_sist(coord,vel)
real*8, dimension(Num,3), intent(out) :: coord, vel
real*8 ,dimension(3) :: vcm
integer :: i,j
real*8 :: rn,temp_calc,rd

! Subrutina que prepara el sistema, colocando las particulas de manera aleatoria
! dentro de la caja y asignandoles una velocidad inicial aleatoria

! inicializacion del generador de numeros aleatorios con distribucion uniforme
call srand(iseed)

vcm=0.
t=0.
temp_calc=0.

! Distribucion de las particulas de manera aleatoria dentro de la caja

do i=1,Num


  rn=rand()
  coord(i,1)=L*(rn-0.5)
  rn=rand()
  coord(i,2)=L*(rn-0.5)
  rn=rand()
  coord(i,3)=L*(rn-0.5)


! Asignacion de las velocidades iniciales mediante una distrbucion aleatoria uniforme entre -0.5 i 0.5.

  rn=rand()
  vel(i,1)=rn-0.5
  rn=rand()
  vel(i,2)=rn-0.5
  rn=rand()
  vel(i,3)=rn-0.5

! Calculo de la temperatura instantanea del sistema asociada a las velocidades asignadas

  temp_calc=temp_calc+m*(vel(i,1)**2+vel(i,2)**2+vel(i,3)**2)/(Kb*3*(Num-1))


  end do

! Se reescala la velocidad en funcion de la temperatura del sistema y la deseada

  vel=vel*sqrt(temp/temp_calc)


! Calculo de la velocidad del centro de massas 

  vcm(1)=sum(vel(:,1))
  vcm(2)=sum(vel(:,2))
  vcm(3)=sum(vel(:,3))

  vcm=vcm/Num

! Cambio del origen de las velocidades para centrar la velocidad del centro de massas
! del sistema en 0.

   vel(:,1)= vel(:,1)-vcm(1)
   vel(:,2)= vel(:,2)-vcm(2)
   vel(:,3)= vel(:,3)-vcm(3)
 
  


end subroutine

subroutine out_coord
integer:: i
! Subrutina que actualiza el output con las coordenadas de todas las particulas
! del sistema a cada tiempo. Este output se ha escrito en un formato compatible
! con el programa VMD que permite seguir la evolucion del sistema en forma de 
! pelicula.



open (UNIT=9, file="coord.xyz")

write(9,*) Num
write(9,*) "trajectoria, iter:",iter, "temps",t, "ns"

 do i=1,Num
 write(9,*) "C ", coord(i,1), coord(i,2), coord(i,3)
 end do



end subroutine

subroutine check(vel,Num)
real*8, dimension(Num,3), intent(in) :: vel
integer, intent(in) :: Num
integer:: i
real*8, dimension(3) :: vcm_calc

! Subrutina de control que calcula la velocidad del centro de massa
! a tiempo t y muestra un aviso por pantalla en caso de ser diferente 
! de 0 a fin de detectar un mal funcionamiento del programa.


 vcm_calc=0.

  do i=1,Num
  
  vcm_calc(1)=vcm_calc(1)+vel(i,1)
  vcm_calc(2)=vcm_calc(2)+vel(i,2)
  vcm_calc(3)=vcm_calc(3)+vel(i,3)
  
  end do

  vcm_calc=vcm_calc/Num
  
 
  if ((int(vcm_calc(1)) /= 0 ).or.(int(vcm_calc(2)) /= 0).or.(int(vcm_calc(3)) /= 0 )) then

  print *, "WARNING!!! la velocidad del centro de massas es diferente de 0!"

  end if

 
 end subroutine

subroutine dist_pbc(x1,y1,z1,x2,y2,z2,rx_pbc,ry_pbc,rz_pbc)
 real*8, intent(in) :: x1,y1,z1,x2,y2,z2
 real*8, intent(out) :: rx_pbc,ry_pbc, rz_pbc
 real*8 :: rx,ry,rz
 
! Subrutina que dadas las coordenadas de dos particulas devuelve la distancia menor
! entre ellas teniendo en cuenta las PBC 


  rx=x1-x2
  ry=y1-y2
  rz=z1-z2





if (abs(rx) > (L/2.)) then     
      rx_pbc=rx-L*nint(rx/L) 
      else
      rx_pbc=rx
      end if
if (abs(ry) > (L/2.)) then
      ry_pbc=ry-L*nint(ry/L)
      else
      ry_pbc=ry
      end if
if (abs(rz) > (L/2.)) then
      rz_pbc=rz-L*nint(rz/L)
      else
      rz_pbc=rz
      end if




end subroutine

subroutine Force(rx_pbc,ry_pbc,rz_pbc,r_pbc,vx,vy,vz,fx,fy,fz)
 real*8, intent(in) :: rx_pbc,ry_pbc,rz_pbc,r_pbc,vx,vy,vz
 real*8, intent(out) :: fx,fy,fz
 real*8, dimension(3) :: F_c,F_D,F_R,F_t
 real*8 :: w_D
 integer :: k

! Subrutina que calcula la fuerza entre dos particulas utilizando
! una repulsion soft-core, una fuerza aleatoria y una fuerza disipativa

call Conservative_Force(rx_pbc,ry_pbc,rz_pbc,r_pbc,F_c)
call Disipative_Force(rx_pbc,ry_pbc,rz_pbc,r_pbc,vx,vy,vz,F_D,w_D)
call Random_Force(rx_pbc,ry_pbc,rz_pbc,r_pbc,w_D,F_R)


do k=1,3

F_t(k)=F_c(k)+F_D(k)+F_R(k)

end do

fx=F_t(1)
fy=F_t(2)
fz=F_t(3)


end subroutine

subroutine Conservative_Force(rx,ry,rz,r,F_c)
real*8,intent(in):: rx,ry,rz,r
real*8,dimension(3),intent(out):: F_c
real*8,dimension(3):: r_vec
integer :: k

r_vec(1)=rx/r
r_vec(2)=ry/r
r_vec(3)=rz/r


if (r < 1) then

do k=1,3

F_c(k)=Fc_a*(1-r)*r_vec(k)


end do
else

do k=1,3


F_c(k)=0

end do


end if



end subroutine

subroutine Disipative_Force(rx,ry,rz,r,vx,vy,vz,F_D,w_D)
real*8,intent(in):: rx,ry,rz,r,vx,vy,vz
real*8,dimension(3),intent(out):: F_D
real*8,intent(out):: w_D
real*8,dimension(3):: r_vec,v_vec
real*8 :: p_esc
integer :: k

r_vec(1)=rx/r
r_vec(2)=ry/r
r_vec(3)=rz/r
v_vec(1)=vx
v_vec(2)=vy
v_vec(3)=vz

if (r<1) then

w_D=(1-r)**2

else

w_D=0

endif

p_esc=0

do k=1,3

p_esc=p_esc+r_vec(k)*v_vec(k)

end do


do k=1,3

F_D(k)=-1*(sig**2/2.)*w_D*p_esc*r_vec(k)

end do



end subroutine 

subroutine Random_Force(rx,ry,rz,r,w_D,F_R)
real*8,intent(in):: rx,ry,rz,r,w_D
real*8,dimension(3),intent(out):: F_R
real*8,dimension(3):: r_vec
integer :: k
real*8 ::w_R,rn_2

r_vec(1)=rx/r
r_vec(2)=ry/r
r_vec(3)=rz/r


rn_2=gauss()
w_R=sqrt(w_D)

do k=1,3

F_R(k)=sig*w_R*rn_2*r_vec(k)*ts**(-1./2.)
end do


end subroutine
subroutine New_Pos(coord,vel,frz)
 real*8, dimension(Num,3),intent(in) :: frz,vel
 real*8, dimension(Num,3),intent(inout) :: coord

!Calculo de la posicion a t+ts utilizando el integrador Velocity Verlet

 do i=1,Num

  coord(i,1)=coord(i,1)+vel(i,1)*ts+(ts**2)*frz(i,1)/(2*m)
  coord(i,2)=coord(i,2)+vel(i,2)*ts+(ts**2)*frz(i,2)/(2*m)
  coord(i,3)=coord(i,3)+vel(i,3)*ts+(ts**2)*frz(i,3)/(2*m)
 
 end do


end subroutine

subroutine Vel_Corr(frz,vel)
 real*8, dimension(Num,3),intent(in) :: frz
 real*8, dimension(Num,3),intent(inout) :: vel
! Calculo de la velocidad a t+ts utilizando el integrador Velocity Verlet


 do i=1,Num


   vel(i,1)=vel(i,1)+(frz(i,1))*ts/(2*m)

   vel(i,2)=vel(i,2)+(frz(i,2))*ts/(2*m)

   vel(i,3)=vel(i,3)+(frz(i,3))*ts/(2*m)


end do




end subroutine







subroutine New_vel(frz,frzt,vel)
 real*8, dimension(Num,3),intent(in) :: frz,frzt
 real*8, dimension(Num,3),intent(inout) :: vel
! Calculo de la velocidad a t+ts utilizando el integrador Velocity Verlet

 
 do i=1,Num


   vel(i,1)=vel(i,1)+(frz(i,1)+frzt(i,1))*ts/(2*m)

   vel(i,2)=vel(i,2)+(frz(i,2)+frzt(i,2))*ts/(2*m)

   vel(i,3)=vel(i,3)+(frz(i,3)+frzt(i,3))*ts/(2*m)


end do


  

end subroutine

subroutine Coord_PBC(coord)

 real*8, dimension(Num,3),intent(inout) :: coord
 
! Subrutina que dadas las nuevas posiciones de las particulas
! vuelve a colocar dentro de la caja de simulacion todas sus posiciones
! aplicando las PBC


 do i=1,Num


   if ( coord(i,1) > L/2.) then
   coord(i,1)=coord(i,1)-L
   else if (coord(i,1) < -L/2) then
   coord(i,1)=coord(i,1)+L
   end if

   if (coord(i,2) > L/2) then
   coord(i,2)=coord(i,2)-L
   else if (coord(i,2) < -L/2) then
   coord(i,2)=coord(i,2)+L
   end if

   if (coord(i,3) > L/2) then
   coord(i,3)=coord(i,3)-L
   else if (coord(i,3) < -L/2) then
   coord(i,3)=coord(i,3)+L
   end if

 end do
   
  end subroutine 



 subroutine Temperature(vel,temp_calc)
 real*8, dimension(Num,3),intent(in) :: vel
 real*8, intent(out) :: temp_calc
! Subrutina que calcula la temperatura instantanea del sistema
! en un instante t 


 temp_calc=0

 do i=1,Num
 temp_calc=temp_calc+m*(vel(i,1)**2+vel(i,2)**2+vel(i,3)**2)/(Kb*3*(Num-1)) 
 end do

 end subroutine

 subroutine radial_distribution(coord)
 real*8, dimension(Num,3),intent(in) :: coord
 real*8 :: r, dr=0.2, N_r, g_r
! Subrutina que calcula la distribuacion radial de unas particulas
! a partir de una matriz de posiciones xyz 


 open (UNIT=15, file='radial_distribution.txt')

  r=0.

  write(15,*) '#       Radi           g(r) '

  do while (r < L/2.)

  r=r+dr
  g_r=0 
 
  do i=1,Num

  N_r=0.
    do j=1,Num
      if (i /= j) then
       call dist_pbc(coord(i,1),coord(i,2),coord(i,3),coord(j,1),coord(j,2),coord(j,3),rx_pbc,ry_pbc,rz_pbc)
      r_pbc=sqrt(rx_pbc**2+ry_pbc**2+rz_pbc**2)
      
         if ((r_pbc>r).and.(r_pbc < (r+dr))) then 
         N_r=N_r+1
         end if 
   
      end if

    end do
   
   g_r=g_r+N_r/(dens*4*pi_c*((r+dr)**3-r**3)/3)
   

  end do

  g_r=g_r/Num

  write(15,*) r, g_r

 end do

end subroutine

subroutine vmd_script
open (UNIT=16, file='View_vmd.tcl')


write (16,*) " "
write (16,*) "mol delete top"
write (16,*) "mol load xyz coord.xyz"
write (16,*) "mol delrep 0 top"
write (16,*) "display resetview"
write (16,*) " "
write (16,*) "set selC [atomselect top \""name C\""]"
write (16,*) "$selC set radius  0.5 "
write (16,*) "color Display Background white"
write (16,*) "mol representation VDW 0.500000 16.000000"
write (16,*) "mol selection name C"
write (16,*) "mol material Opaque"
write (16,*) "mol color ColorID 2"
write (16,*) "mol addrep top"
write (16,*) " "

close(16)

end subroutine

subroutine VACF(vel_i,vel,C_VACF)
real*8, dimension(Num,3), intent(in)::vel_i,vel
real*8, intent(out):: C_VACF
integer :: k

C_VACF=0

do k=1,Num

C_VACF=C_VACF+vel_i(k,1)*vel(k,1)+vel_i(k,2)*vel(k,2)+vel_i(k,3)*vel(k,3)

end do

C_VACF=C_VACF/Num

end subroutine




      function gauss()
 10   continue
      U1 = r1279()
      U2 = r1279()
      V1 = 2.*U1 - 1.
      V2 = 2.*U2 - 1.
      S = V1**2 + V2**2
      IF (S.GE.1.) GOTO 10
      X1 = V1*SQRT(-2.*LOG(S)/S)
      X2 = V2*SQRT(-2.*LOG(S)/S)
      gauss=x1
      END function

 end program

