      subroutine BInitIter
     & (N,u,v,w,itime,ii,Np,pa,tmark,t,to,dt,td,
     & nrpars,rpars,dom,visc,i1,i2)
      implicit none
      integer N, itime, ii, Np, i1, i2, j, nrpars
      double precision u(N), v(N), w(N),rpars(100), BExp, pav
      double precision dom, visc, pa
      double precision x, tmark, t, to, dt, td, tdfac, test
 100  format(2g16.8)
c
c Initialize counters.
c
      itime=0
      ii=0
      Np=0
c
c Initialize acceptance probability accumulation.
c
      pa=0.d0
c
c Initialize time variables.
c
      tmark=0.d0  !initial elapsed time during subinterval
      to=0.d0  !initial physical time
c
c Compute dt as the product of:
c     viscous time scale                             (dom/N)^2 / visc
c     target average acceptance probability           pav
c     inverse of number of simultaneous small eddies  1/N
c (default dt computation)
c
      pav=rpars(2)
      dt=pav*dom**2/(visc*dfloat(N)**3)
      if(nrpars.ge.6) then
       test=rpars(6)
       if(test.gt.0.d0) then !test is the initial dt value
        dt=test
       else                  !dt = |test| * dt
        dt=abs(test)*dt
       endif
      endif
c
c Set the scheduled time t of the first candidate eddy
c by sampling from the exponential distribution with mean dt.
c
      t=BExp(dt,i1,i2)
      tdfac=rpars(4)  !td/dt (enforce a fixed ratio)
      td=tdfac*dt  !maintain specified proportionality td/dt
c
c Read the initial velocity profiles.
c
      open(1, file="U.dat", status="unknown")
      open(2, file="V.dat", status="unknown")
      open(3, file="W.dat", status="unknown")
      do j=1, N
       read(1, 100) x, u(j)
       read(2, 100) x, v(j)
       read(3, 100) x, w(j)
      enddo
      do j=1, 3
       close(j)
      enddo
      return
      end
