      subroutine BSampleEddy
     & (N,M,L,u,v,w,dt,td,PL,ratefac,viscpen,uK,vK,wK,
     & pp,ii,pa,Np,rpars,Io,Ip,Cv,Co,i1,i2)
      implicit none
      integer N, M, L, L3, Np, ii, Io, Ip, i1, i2
      double precision dt, td, pp, pa, u(N), v(N), w(N)
      double precision ratefac, viscpen, uK, vK, wK
      double precision PL(N), rpars(100), rannum, Cv, Co
c
c Sample eddy size and location and compute the eddy acceptance probability.
c
      ii=ii+1  !increment eddy-trial counter
c
c When ii reaches a threshold value, an increase of  
c the eddy-trial time increment dt is attempted in BRaisedt.
c
      call BLength(N,PL,L3,Io,Ip,Cv,Co,i1,i2)  !sample L3=L/3
      L=3*L3  !eddy size
c
c The location of the leftmost cell is sampled uniformly within [1,N-L].
c The possibility N-L+1 is excluded because the rightmost cell would 
c then be at location N.  This is disallowed because the continuum 
c spatial range is taken to be (dom/N)*[0,N], where dom is the physical 
c domain size.  Property array values correspond to locations j*(dom/N) 
c where the integer j ranges from 1 to N.  Therefore the leftmost cell 
c of an eddy can be no lower than j=1.  For consistent treatment at the 
c other boundary, the rightmost cell of an eddy is allowed to be no 
c higher than j=N-1.  Therefore size N eddies cannot be implemented.  
c Accordingly, in the run initialization, Iv is bounded from above to  
c enforce L<N.
c
      call Brng(rannum,i1,i2)
      M=1+min(N-L-1,int(rannum*(N-L)))  !sample the leftmost cell
      call BProb(N,M,L3,u,v,w,dt,PL,ratefac,viscpen,uK,vK,wK,pp)
c
c If the acceptance probability pp exceeds a specified maximum value
c pmax, then reduce dt so that pp=pmax.  Also, gather acceptance 
c probability statistics needed in subroutine BRaisedt.
c
      call BLowerdt(dt,td,pp,pa,Np,rpars)
      return
      end
