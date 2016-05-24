module statistics
implicit none
public
  !! Rutinas de estadisticas: g(r)
 
  !!! declarate_radial_dist: declara las variables para la g(r)
  !   Entrada:
  !     Tamaño L de la caja. La g(r) se mide entre 0 y 0.4*L
  !     Numero de espacios (n_esp) en los que discretizar la g(r)
  !   Salida:
  !     histograma (g) con n_esp huecos iniciados a cero
  !     ancho (dr) de cada hueco
  !
  !!! accumulate_radial_dist: actualiza la g(r) cada vez que se la llama
  !   Entrada:
  !     Numero de atomos (M)
  !     vector de posiciones
  !     numero de espacios (total_esp) y anchura (dr)
  !   Entrada/Salida:
  !     histograma (g)
  !
  !!! compute_radial_dist: normaliza la g(r) y la escribe en un fichero
  !   Entrada:
  !     canal de escritura (channel)
  !     numero de espacions (total_esp) y anchura (dr)
  !   Entrada/Salida:
  !     histograma (g)
    
contains

  subroutine declarate_radial_dist(L,total_esp,dr,g)
  real,intent(in) :: L
  integer,intent(in) :: total_esp
  real,intent(out) :: dr
  real,dimension(total_esp),intent(out) :: g
  real :: rmax
  
  rmax = 0.4*L
  dr = rmax/total_esp
  g(:) = 0
    
  end subroutine
  
  subroutine accumulate_radial_dist(M,r,total_esp,dr,g)
  ! Datos de entrada y salida
  integer,intent(in) :: M
  real,dimension(3,M),intent(in) :: r
  integer,intent(in) :: total_esp
  real,intent(in) :: dr
  real,dimension(total_esp),intent(inout) :: g  !salida g(r) (estadistica acumulada)
  ! Variables de calculo
  real,dimension(3) :: Ri, Rj, Rij
  real :: distance
  integer :: i,j
  integer :: indice
  
  do i=2,M          !empezar con el par i=2 j=1. Evitar i=1 j=0
    Ri = r(:,i)
    do j=1,(i-1)    !bucle sobre todos los pares de atomos
       Rj = r(:,j)
       Rij = Ri(:) - Rj(:)      !vector que une los dos atomos
       distance = sqrt(Rij(1)**2 + Rij(2)**2 + Rij(3)**2)    !distancia relativa
       indice = floor(distance/dr) + 1  !posicion en histograma
       if (indice <= total_esp) then
          g(indice) = g(indice) + 1  ! Si distance<r_max, la g(r) suma 1
       end if
    end do
  end do
  
  end subroutine
  
  subroutine compute_radial_dist(channel, total_esp, dr, g)
  integer,intent(in) :: channel
  integer,intent(in) :: total_esp
  real,intent(in) :: dr
  real,dimension(total_esp),intent(inout) :: g  !salida g(r) normalizada
  real :: vol, pi, R, f
  integer :: factor_normalizacion
  integer :: i
  
  pi = 3.14159265359
  f = 4*pi/3.
  
  do i=1, total_esp       !corregir con el volumen de los casquetes
    R = (i-1)*dr
    vol = f*( (R + dr)**3 - R**3 )
    g(i) = g(i)/vol
  end do
  
  do i=1,total_esp
    write(channel,*) dr*i, g(i)  
  end do
  
  end subroutine
  
end module
