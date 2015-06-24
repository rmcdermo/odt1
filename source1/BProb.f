       subroutine
     &  BProb(N,M,L3,u,v,w,dt,PL,ratefac,viscpen,uK,vK,wK,pp)
       implicit none
       integer N, M, L, L3
       double precision pp, dt, uK, vK, wK, p
       double precision disfac, BsKd, ratefac, viscpen
       double precision u(N), v(N), w(N), PL(N)
c
c Compute the eddy acceptance probability pp.  The model 
c expression for pp is BI(lambda) * dt / (Prob(L) * Prob(M)), 
c where lambda is the eddy rate distribution evaluated for the 
c candidate eddy, dt is the eddy sampling time step, Prob( )
c denotes the probability of sampling the current value of the 
c argument, and BI denotes the integral of lambda over the ranges
c of eddy size and location corresponding to the size and location
c bins used in the sampling of L (discrete size) and M (discrete 
c location).
c
c Prob(L) and Prob(M) are evaluated using the probability
c distributions used to sample L and M, as in subroutines
c BLenProb, BLength, and BSampleEddy.  Instead of performing the
c bin integral (BI) over the continuous variable lambda, lambda is
c approximated as constant within the bin and assigned the value
c corresponding to the discrete values L and M (integer numbers of
c mesh cells).  The approximate evaluation is then
c BI(lambda) = lambda*(3*dom/N)*(dom/N)
c where the second and third factors are the dimensional widths of
c the L and M bins respectively.  dom/N is the physical size of a  
c mesh cell.  (dom is the the physical size of the ODT domain and 
c N is the number of mesh cells.)
c
c The continuum model expression for lambda is 
c     (C/S^3)[u_K^2 + v_K^2 + w_K^2 - Z(visc/S)^2]^(1/2)
c where visc is the kinematic viscosity and S is the physical 
c (dimensional) eddy size L*dom/N.  The quantity in square brackets is 
c evaluated dimensionally.  The viscous penalty Z(visc/S)^2 is 
c evaluated as viscpen/(L*L) where viscpen = Z(visc*N/dom)^2 is an 
c input to the subroutine.  On this basis, 
c             BI(lambda) = (3*dom/N)*(dom/N)*lambda
c    = (ratefac/L^3)[u_K^2 + v_K^2 + w_K^2 - viscpen/(L*L)]^(1/2)
c where ratefac = 3*C*N/dom is an input to the subroutine.
c
       L=3*L3
       uK=BsKd(N,M,L,u)
       vK=BsKd(N,M,L,v)
       wK=BsKd(N,M,L,w)
       p=(uK*uK)+(vK*vK)+(wK*wK)-viscpen/dfloat(L)**2
       if (p .gt. 0.d0) then
c
c The mean-square displacement by a size-L discrete triplet map is 
c (1 - 3/L) times the mean-square displacement by a continuum 
c triplet map (mathematical definition) of the same size.  (This is 
c conveniently derived by applying mathematical induction separately 
c to discrete eddies for which L/3 is even and odd, respectively.)
c The turbulent transport induced by eddies is the event rate times 
c the mean-square displacement.  This suggests compensating for the 
c discrete deviation from continuum displacement by changing the 
c event rate so that the product of rate and mean-square 
c displacement equals the continuum value of the product.  This 
c correction is found to improve convergence of results as mesh 
c resolution increases.
c
        disfac=(1.d0-3.d0/dfloat(L)) !discrete/continuum <displacement^2>
c
c Divide rate by disfac to compensate for discretization error.
c
        pp = (ratefac/(disfac*dfloat(L)**3))*dsqrt(p) *dt 
     &   *(N-L) /(PL(L3)-PL(L3-1))
c
c On the right hand side, the four terms separated by spaces are BI(lambda), 
c eddy-sampling time increment, inverse of eddy-location sampling 
c probability, and inverse of eddy-size sampling probability.
c
       else
        pp=0.d0
       endif
       return
       end
