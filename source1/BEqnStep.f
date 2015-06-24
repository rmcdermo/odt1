      subroutine BEqnStep
     & (N,u,v,w,delt,dom,visc,pgrad,rpars,istat,cstat,
     & NVAL,mVAL,nstatVAL)
      implicit none
      integer N, i, irat, istat, NVAL, mVAL, nstatVAL
      double precision delt, dom, visc, pgrad, zero
      double precision tfrac, tcfl, et
      double precision u(N), v(N), w(N), rpars(100)
      double precision cstat(NVAL,mVAL,nstatVAL)
      zero=0.d0
c
c Subdivide the time advancement interval delt into equal 
c sub-intervals et that are no larger than a specified maximum, 
c which is tfrac times the CFL time for the given domain size 
c dom, number N of mesh cells, and kinematic viscosity (or 
c mathematically equivalent transport coefficient) visc.
c
c If an implicit advancement scheme is introduced in a future 
c version of the code, then the sub-division into time steps 
c of size et will not be needed.  In the present version, the 
c bound td on the time between invocations of BEqnStep is 
c superfluous due to the sub-cycling implemented below, but 
c the bound td would be needed to control the time step in an 
c implicit formulation.
c
      tfrac=rpars(5)
      tcfl=0.5*(dom/dfloat(N))**2/visc  !CFL time
      irat=delt/(tfrac*tcfl)
c
c irat is the integer truncation of delt divided by the 
c maximum time-step magnitude.  For uniform time steps, et, 
c during delt that are as large as possible without 
c exceeding the bound, divide delt by 1 + irat.
c 
      et=delt/(dfloat(irat)+1.d0)
c
c Advance property fields for 1+irat time steps of duration et.
c
      do i=1, irat+1
      call BStats
     & (N,u,v,w,dom,visc,cstat,et,NVAL,mVAL,nstatVAL,istat)
       call BAdv(N,u,et,dom,visc,pgrad)
       call BAdv(N,v,et,dom,visc,zero)
       call BAdv(N,w,et,dom,visc,zero)
      enddo
      return
      end
