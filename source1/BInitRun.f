      subroutine BInitRun(N,i1,i2,dom)
      implicit none
      integer N, j, k, i1, i2
      double precision c, r, rannum, x, dom
 100  format(2g16.8)
      open(11, file="Brun.dat", status="unknown")
      write(11, *) "Starting the run"
      call flush(11)
c
c Generate initial (u,v,w) profiles now for all realizations 
c so they won't have different random perturbations 
c for different realizations.  The random perturbations 
c prevent ambiguity in the choice of solution of the 
c quadratic equation determining the amplitudes of 
c applied kernels.
c
c The physical domain is [0,dom].  Mesh points are located 
c at j*dom/N for integers j from 1 to N.  Thus there is an 
c asymmetry in the representation of the lower and upper 
c domain boundaries.  The model is coded so there is no 
c statistical asymmetry in model implementation.
c
      c=1.d-8
      open(1, file="U.dat", status="unknown")
      open(2, file="V.dat", status="unknown")
      open(3, file="W.dat", status="unknown")
      do j=1, N-1
       x=dom*j/dfloat(N)
       do k=1, 3
        call Brng(rannum,i1,i2)
        r=c*(rannum-0.5d0)
        write(k,100) x, r
       enddo
      enddo
      c=0.d0
      do k=1, 3 !no perturbation at upper wall
       write(k,100) dom, c
       close(k)
      enddo
      return
      end

