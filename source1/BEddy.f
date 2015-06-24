       subroutine BEddy(N,M,L,u,v,w,uK,vK,wK)
       implicit none
       integer N, M, L
       double precision uK, vK, wK, dom, root, cu, cv, cw
       double precision u(N), v(N), w(N)
       double precision disfac, cfac
c
c Add c_s*K to triplet-mapped velocity component s (= u, v, or w), 
c where the physical model expression for c_s*K is
c c_s*K = (27/4)*(K/L)*[-sK + sign(sK)*((uK^2+vK^2+wK^2)/3)^(1/2)].
c c_s is evaluated using physical units for quantities in the square 
c brackets but taking the eddy size L to be an integer number of 
c mesh spacings.  In subroutine BAddK, c_s is multiplied by K in 
c the same units as L.  Thus the factor (K/L) above is evaluated 
c consistently.
c 
c In the physical model, 27/4 is the value of 
c l^3/(2 * integral of K^2).  In the numerical implementation, K 
c is an integer rather than a physical length (see subroutine BAddK) 
c and its numerical evaluation on the discrete mesh gives 
c (integral of K^2) = (4/27)*L^3*(1-3/L).  This reflects displacements 
c as implemented numerically, so the physical model expression for 
c c_s*K is multiplied by the factor 1-3/L in the numerical 
c implementation.
c
c compute 27/(4 * L * discrete-mesh modification of K^2 integral)
c
       disfac=(1.d0-3.d0/L)  !discrete-mesh modification of K^2 integral
       cfac=6.75d0/(disfac*dfloat(L))
       root=dsqrt((uK*uK+vK*vK+wK*wK)/3.d0)
       cu=cfac*(-uK+sign(1.d0,uK)*root)
       cv=cfac*(-vK+sign(1.d0,vK)*root)
       cw=cfac*(-wK+sign(1.d0,wK)*root)
       call BTriplet(N,M,L,u)
       call BTriplet(N,M,L,v)
       call BTriplet(N,M,L,w)
c
c After the property profiles are triplet mapped, add c_s*K to 
c velocity component s for s = u, v, and w, respectively.
c
       call BAddK(N,M,L,u,cu)
       call BAddK(N,M,L,v,cv)
       call BAddK(N,M,L,w,cw)
       return
       end
