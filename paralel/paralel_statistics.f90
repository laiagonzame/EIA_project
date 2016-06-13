module statistics
implicit none
public
  !! Modulo con las rutinas de estadisticas
  !! Paralelizacion de la g(r)
  
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
  real,intent(in) :: L
  integer,intent(in) :: total_esp
  real,intent(out) :: dr
  real,dimension(total_esp),intent(out) :: g
  real*8 :: rmax
  
  rmax = 0.5*L
  dr = rmax/total_esp
  g(:) = 0
    
  end subroutine
  
  subroutine accumulate_radial_dist(M,L,r,total_esp,dr,dist_rad)
  include 'mpif.h'
  ! Datos de entrada y salida
  integer,intent(in) :: M
  real,intent(in) :: L
  real,dimension(3,M),intent(in) :: r
  integer,intent(in) :: total_esp
  real,intent(in) :: dr
  real,dimension(total_esp),intent(out) :: dist_rad  !salida g(r)

  ! Variables de calculo
  real,dimension(3) :: Ri, Rj, Rij
  real*8 :: distance
  integer :: indice
  integer :: i, j
  
  ! parallelization variables
  integer :: MASTER, taskid, ierror, numproc, num_pairs, counter, request
  integer :: stat(MPI_STATUS_SIZE), pairs_per_worker
  integer, dimension(:), allocatable :: list_i, list_j, start, finish
  double precision, dimension (:,:), allocatable :: g
  integer :: k
  
  MASTER = 0
  
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierror)
  
  num_pairs = M*(M-1)/2
  pairs_per_worker = num_pairs/(numproc-1)
  
  allocate(list_i(num_pairs),list_j(num_pairs))
  allocate(start(numproc-1),finish(numproc-1))
  allocate(g(total_esp,numproc-1))
  
  counter = 1
  do i=2,M
    do j=1,(i-1)
      list_i(counter) = i
      list_j(counter) = j
      counter = counter + 1
    end do
  end do
  
  do i=1,numproc-1
    start(i) = (taskid-1)*pairs_per_worker+1
    finish(i) = taskid*pairs_per_worker
  end do
  
  finish(numproc-1) = num_pairs
    
  g(:,:) = 0
  dist_rad(:) = 0

  if(taskid .ne. MASTER) then   
 
    do k=start(taskid),finish(taskid)
       i = list_i(k)
       j = list_j(k)

       Ri = r(:,i)
       Rj(:) = r(:,j)
       
       Rij(:) = Ri(:) - Rj(:)   !vector que los une
  
       Rij(:) = Rij(:) - L*nint(Rij(:)/L)  !imagen minima

       distance = sqrt(Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3))    !distancia relativa
       indice = floor(distance/dr) + 1  !posicion en histograma
       if (indice <= total_esp) then
          g(indice,taskid) = g(indice,taskid) + 1  ! Si distance<r_max, la g(r) suma 1
       end if
    end do
    
    call MPI_ISEND(g(:,taskid), total_esp, MPI_REAL8, MASTER, 1, MPI_COMM_WORLD, request, ierror)
  end if
  
  call MPI_BARRIER(MPI_COMM_WORLD, ierror)
  
  if (taskid .eq. MASTER) then
    do i = 1, numproc-1
       call MPI_RECV(g(:,i), total_esp, MPI_REAL8, i, 1, MPI_COMM_WORLD, stat, ierror)
       dist_rad(:) = dist_rad(:) + g(:,i)
    end do
   
    do i = 1, numproc-1
       call MPI_ISEND(dist_rad(:), total_esp, MPI_REAL8, i, 1, MPI_COMM_WORLD, request, ierror)
    end do
  end if
 
  call MPI_BARRIER(MPI_COMM_WORLD, ierror)

  if (taskid .ne. MASTER) then
    call MPI_RECV(dist_rad(:), total_esp, MPI_REAL8, MASTER, 1, MPI_COMM_WORLD, stat, ierror)
  end if

  call MPI_BARRIER(MPI_COMM_WORLD, ierror)

  call MPI_FINALIZE(ierror)
  
  end subroutine
  
  subroutine compute_radial_dist(total_esp, dr, n_acc, density, g)
  integer,intent(in) :: total_esp
  real,intent(in) :: dr
  integer,intent(in) :: n_acc                   !numero de acumulaciones (normalizacion)
  real,intent(in) :: density                    !densidad (normalizacion)
  real,dimension(total_esp),intent(inout) :: g  !salida g(r) normalizada
  real*8 :: vol, pi, R, f
  integer :: i
  
  pi = 3.14159265359
  f = 4*pi*density*n_acc/3.  !factor de correccion y normalizacion
  
  do i=1, total_esp       !corregir con el volumen de los casquetes
    R = (i-1)*dr
    vol = f*( (R + dr)**3 - R**3 )
    g(i) = g(i)/vol
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
  
  !!  Cómo calcular el DCM: Guardar la config a tiempo T0 en pos0. Leer posiciones en un tiempo, pos, calcular el dcm y escribir en un archivo el tiempo (T-T0) y el dcm
  
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
