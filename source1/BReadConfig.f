      subroutine BReadConfig
     & (niter,nstat,ntseg,tend,dom,visc,pgrad,ratefac,viscpen,N)
c
c Read the inputs needed for the default case .
c
      implicit none
      integer N, niter, nstat, ntseg
      double precision tend, dom, visc, pgrad, C, Z
      double precision ratefac, viscpen
      open(12, file="BConfig.dat", status="old")
c
c N and nstat must not exceed the values of NVAL and nstatVAL, 
c respctively, that are hardwired in parameter statements in the main 
c program (Bodt) .
c The hardwired values can be changed in Bodt to suit the application.
c
      read(12, *) niter  !number of iterations (flow realizations)
      read(12, *) nstat  !number of data-averaging periods per iteration
      read(12, *) ntseg  !number of snapshot intervals per data-avg'g period
      read(12, *) tend   !time duration of a realization (s)
      read(12, *) dom    !domain size (m)
      read(12, *) visc   !viscosity (m^2/s)
      read(12, *) pgrad
c pgrad = (1/density) * imposed mean pressure gradient (m s^-2)
      read(12, *) C	!ODT event rate parameter
      read(12, *) Z	!ODT viscous penalty parameter
      close(12)
c convert the ODT parameters to forms used in function BProb
      ratefac=3*C*dfloat(N)/dom
      viscpen=Z*(visc*dfloat(N)/dom)**2
      return
      end
