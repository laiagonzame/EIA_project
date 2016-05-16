module forces
implicit none

contains

subroutine forces(M,r,F,boxlength,sigma,epsil)

integer, intent(in) :: M
real(8), dimension(3,M), intent(in) :: r
real(8), dimension(3,M), intent(out) :: F
real(8), intent(in) :: boxlength, sigma, epsil
real(8), dimension(3) :: rij 
real(8) :: rcut, rijl, dist2, sigmar, rc2
integer :: i, j, l

rcut=boxlength
rc2=rcut*rcut

do i=1,M
   do l=1,l
      F(l,i) = 0.0
   enddo
enddo

do i=1,M-1
   
   do j=i+1,M

      dist2=0.0
      do l=1,3
         rijl = r(l,j) - r(l,i)
         rij(l) = rijl - boxlength*nint(rijl/boxlength)
         dist2 = dist2 + rij(l)*rij(l)
      enddo

      if (dist2<rc2) then
         sigmar = (dist2**(-3.0))*(sigma**6.0)
         force = 24.0*epsil*(2.0*sigmar**(2.0)-sigmar)/dist2
         do l=1,3
            F(l,i) = F(l,i) - force*rij(l)
            F(l,j) = F(l,j) + force*rij(l)
         enddo
      endif

   enddo   

enddo


endsubroutine
