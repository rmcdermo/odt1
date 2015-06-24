      subroutine
     & BAdv(N,r,tstep,dom,visc,force) 
      implicit none
      integer N, j
      double precision tstep, dom, visc, force, dz
      double precision r(N), f(N)
c
c Advance the equation dr/dt = visc * d^2r/dz^2 + force for 
c a time increment tstep, where z is the spatial coordinate, 
c enforcing r(0) = r(dom) = 0, where z = 0 and dom are the 
c lower and upper boundaries.  visc is kinematic viscosity 
c or a mathematically equivalent quantity and force is 
c (pressure gradient)/density when r is the streamwise 
c velocity component, otherwise force is zero.
c
      dz=dom/dfloat(N)  !physical size of mesh cell
c
c Mesh nodes 1 through N correspond to z = dom/N through 
c z = dom, with uniform spacing.
c
c Compute fluxes by first differencing.
c
      f(1)=-visc*r(1)/dz  !here using r(0)=0.
      do j=2, N
       f(j)=visc*(r(j-1)-r(j))/dz
      enddo
c
c Advance the r array.
c
      do j=1, N-1
       r(j)=r(j)+tstep*((f(j)-f(j+1))/dz+force)
      enddo
c
c There is no process that modifies r(N), 
c so it is not updated.
c
      return
      end
