module statistics
implicit none
public
  !! Rutinas de estadisticas: g(r)
 
  !!! declarate_radial_dist: declara las variables para la g(r)
  !   Entrada:
  !     Distancia cutoff del potencial. La g(r) se mide entre 0 y 10*cutoff
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

  subroutine declarate_radial_dist(cutoff,total_esp,dr,g)
  real,intent(in) :: cutoff
  integer,intent(in) :: total_esp
  real,intent(out) :: dr
  real,dimension(total_esp),intent(out) :: g
  real :: rmax
  
  rmax = 10*cutoff
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
  real :: vol, dr3, pi
  integer :: factor_normalizacion
  integer :: i
  
  pi = 3.14159265359
  dr3 = dr*dr*dr
  
  do i=1, total_esp       !corregir con el volumen de los casquetes
    vol = 4*pi*dr3*i**2   !4pi*r^2*dr = 4pi*(i*dr)^2*dr = 4pi*dr^3*i^2
    g(i) = g(i)/vol
  end do
  
  factor_normalizacion = sum(g(:)) !DUDA: ¿Hay que normalizar aqui o mejor lo quito?
  g(:) = g(:)/factor_normalizacion
  
  do i=1,total_esp
    write(channel,*) dr*i, g(i)  
  end do
  
  end subroutine
  
end module
