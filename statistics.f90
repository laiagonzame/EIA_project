module statistics
implicit none
public

!! Rutinas de analisis estadistico sobre un archivo con las trayectorias

contains

  !! DISTRIBUCION RADIAL g(r)
 
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
  integer, intent(in) :: total_esp
  integer, intent(in) :: M
  real*8, intent(in) :: L, dr
  !f2py intent(in) :: M, L, total_esp, dr
  real*8,intent(in) :: r(3,M)
  !f2py intent(in) :: r
  !f2py depend(M) :: r
  real*8, intent(out) :: g(total_esp)  !salida g(r) (estadistica acumulada)
  !f2py intent(out) :: g

  ! Variables de calculo
  real*8,dimension(3) :: Ri, Rj, Rij
  real*8 :: distance
  integer :: i,j
  integer :: indice


  g = 0

  do i=2,M          !empezar con el par i=2 j=1. Evitar i=1 j=0
    Ri = r(:,i)
    do j=1,(i-1)    !bucle sobre todos los pares de atomos

       Rj = r(:,j)
       
       Rij = Ri - Rj   !vector que los une

       Rij = Rij - L*nint(Rij/L)  !imagen minima

       distance = sqrt(Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3))    !distancia relativa
       indice = floor(distance/dr) + 1  !posicion en histograma
       if (indice <= total_esp) then
          g(indice) = g(indice) + 2  ! Si distance<r_max, la g(r) suma 1
       end if
    end do
  end do
  
  end subroutine
  
  subroutine compute_radial_dist(total_esp, dr, n_acc, density, M, g, gout)
  integer,intent(in) :: total_esp, M
  real*8,intent(in) :: dr
  integer,intent(in) :: n_acc                   !numero de acumulaciones (normalizacion)
  real*8,intent(in) :: density                    !densidad (normalizacion)
  real*8, intent(in) :: g(total_esp)  !salida g(r) normalizada
  !f2py depend(total_esp) :: g
  real*8,intent(out) :: gout(total_esp)  !salida g(r) normalizada
  !f2py depend(total_esp) :: gout
  real*8 :: vol, pi, R, f
  integer :: i
  
  pi = 3.14159265359
  f = 4*pi*density*n_acc/3.  !factor de correccion y normalizacion
  
  do i=1, total_esp       !corregir con el volumen de los casquetes
    R = (i-1)*dr
    vol = f*( (R + dr)**3 - R**3 )
    gout(i) = g(i)/(M * vol)
  end do
    
  end subroutine
  
  !! DESPLAZAMIENTO CUADRATICO MEDIO
  
  !!! desp_cuad_medio: Calcula el desplazamiento cuadratico medio de una configuracion
  !                    respecto a otra
  !   Entrada:
  !     numero de atomos (M)
  !     configuracion de referencia (pos0). Posiciones en el tiempo de referencia T0
  !     configuracion actual (pos). Posiciones en el tiempo actual T
  !   Salida:
  !     dcm: desplazamiento cuadratico medio
  
  !!  Cómo calcular el DCM: Guardar la config a tiempo T0 en pos0. Leer posiciones en un tiempo, pos, calcular el dcm
  !   y escribir en un archivo el tiempo (T-T0) y el dcm
  
  !Aviso a Miquel y Laia: Para calcular el DCM hay que tener en cuenta la diferencia de tiempos T-T0.
  !Yo creo que lo mas comodo seria incluir en el archivo de trayectorias el instante de tiempo de cada configuracion.
  !En este modulo, estoy suponiendo que ese dato no está en el archivo, y que Miquel lo debe calcular en el main
  !de statistics.py. Avisadme si se cambia el archivo de la trayectoria incluyendo el dato del tiempo,
  !para que yo actualice mi subrutina
  
  !Aviso para Cris: Esta rutina no tiene en cuenta las salidas de la caja.
  !El DCM se tiene que calcular restando las posiciones reales respecto a las de referencia,
  !no las posiciones aplicadas las PBC. ¿Hay ahora mismo algun seguimiento de la posicion real de las particulas?
  
  subroutine desp_cuad_medio(M,pos0,pos,dcm)
  integer,intent(in) :: M
  real,dimension(3,M),intent(in) :: pos0, pos
  real,intent(out) :: dcm
  real,dimension(3) :: r
  integer :: i
  
  dcm = 0
  
  do i=1, M
    r(:) = pos(:,i) - pos0(:,i)
    dcm = dcm + r(1)*r(1)+r(2)*r(2)+r(3)*r(3)
  end do
  
  dcm = dcm / M
  
  end subroutine
  
end module
