module statistics
implicit none
public
  !! Rutinas de estadisticas: g(r)
 
  !!! declarate_radial_dist: declara las variables para la g(r)
  !   Entrada:
  !     Tamaño L de la caja. La g(r) se mide entre 0 y 0.5*L
  !     Numero de espacios (n_esp) en los que discretizar la g(r)
  !   Salida:
  !     histograma (g) con n_esp huecos iniciados a cero
  !     ancho (dr) de cada hueco
  !
  !!! accumulate_radial_dist: actualiza la g(r) cada vez que se la llama
  !   Entrada:
  !     Numero de atomos (M)
  !     Tamaño L de la caja
  !     vector de posiciones
  !     numero de espacios (total_esp) y anchura (dr)
  !   Entrada/Salida:
  !     histograma (g)
  !
  !!! compute_radial_dist: normaliza la g(r), resultado definitivo
  !   Entrada:
  !     numero de espacions (total_esp) y anchura (dr)
  !     numero de acumulaciones (n_acc) Es el numero de llamadas a accumulate_radial_dis
  !     densidad
  !   Entrada/Salida:
  !     histograma (g)

  !! Cómo calcular la g(r): llamar a la declaración, bucle sobre accumulate para distintos snapshots y, finalmente, compute
    
contains

  subroutine declarate_radial_dist(L,total_esp,dr,g)
  real*8,intent(in) :: L
  integer,intent(in) :: total_esp
  real*8,intent(out) :: dr
  real*8,dimension(total_esp),intent(out) :: g
  real*8 :: rmax
  
  rmax = 0.5*L
  dr = rmax/total_esp
  g(:) = 0
    
  end subroutine
  
  subroutine accumulate_radial_dist(M,L,r,total_esp,dr,g)
  ! Datos de entrada y salida
  integer,intent(in) :: M
  real*8,intent(in) :: L
  real*8,dimension(3,M),intent(in) :: r
  integer,intent(in) :: total_esp
  real*8,intent(in) :: dr
  real*8,intent(out) :: g(total_esp)  !salida g(r) (estadistica acumulada)

  ! Variables de calculo
  real*8,dimension(3) :: Ri, Rj, Rij
  real*8 :: distance
  integer :: i,j
  integer :: indice
  
  g = 0

  do i=2,M          !empezar con el par i=2 j=1. Evitar i=1 j=0
    Ri = r(:,i)
    do j=1,(i-1)    !bucle sobre todos los pares de atomos

       Rj(:) = r(:,j)
       
       Rij(:) = Ri(:) - Rj(:)   !vector que los une

       Rij(:) = Rij(:) - L*nint(Rij(:)/L)  !imagen minima

       distance = sqrt(Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3))    !distancia relativa
       indice = floor(distance/dr) + 1  !posicion en histograma
       if (indice <= total_esp) then
          g(indice) = g(indice) + 1  ! Si distance<r_max, la g(r) suma 1
       end if
    end do
  end do
  
  end subroutine
  
  subroutine compute_radial_dist(total_esp, dr, n_acc, density, g)
  integer,intent(in) :: total_esp
  real,intent(in) :: dr
  integer,intent(in) :: n_acc                   !numero de acumulaciones (normalizacion)
  real,intent(in) :: density                    !densidad (normalizacion)
  real,dimension(total_esp),intent(inout) :: g  !salida g(r) normalizada
  real :: vol, pi, R, f
  integer :: i
  
  pi = 3.14159265359
  f = 4*pi*density*n_acc/3.  !factor de correccion y normalizacion
  
  do i=1, total_esp       !corregir con el volumen de los casquetes
    R = (i-1)*dr
    vol = f*( (R + dr)**3 - R**3 )
    g(i) = g(i)/vol
  end do
    
  end subroutine
  
end module
