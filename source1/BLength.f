       subroutine BLength(N,PL,I,Io,Ip,Cv,Co,i1,i2)
       implicit none
       integer N, I, Io, Ip, nn, i1, i2
       double precision Cv, Co, rannum, x, PL(N)
c
c PL is the exact cumulative distribution function from which 
c the candidate eddy-image size (1/3 of the number of mesh cells 
c in the eddy) is sampled.  The size is sampled by choosing a 
c random number from the uniform distribution over [0,1] and 
c finding the smallest eddy-image size whose PL value exceeds 
c the chosen random number.  To make the search efficient, an 
c accurate initial guess is used.  This guess is obtained by 
c treating eddy-image size as a continuous variable, in which 
c case there is a  closed form solution, denoted x, for the 
c eddy-image size corresponding to a given value, rannum, of 
c the cumulative distribution function of eddy-image sizes.  
c To derive x, f(s) = C s^{-2} exp(-2 Ip/s) is integrated over 
c the allowed eddy-image size range [Io,Iv] and this integral 
c is set equal to unity to determine C.  Then f(s) is integrated 
c over [Io,x] and this integral is set equal to rannum.  This 
c equation is then solved for x, yielding the expression given 
c below, where Cv and Co are evaluated in subroutine BLenProb.  
c In that subroutine, the distribution PL is defined as a
c discretization of the continuous distribution f(s).
c
       call Brng(rannum,i1,i2)
       x=-2.d0*Ip/dlog((Cv*rannum)+(Co*(1.d0-rannum)))
       nn=int(x)-1  !initial guess of the eddy-image size
       if (rannum .gt. PL(nn)) then
 10     nn=nn+1  !increase nn if the guess is too small
        if (rannum .gt. PL(nn)) goto 10
       endif
       if (rannum .lt. PL(nn-1)) then
 20     nn=nn-1  !decrease nn if the guess is too large
        if (rannum .lt. PL(nn-1)) goto 20
       endif
c
c prevent a choice nn < Io which could occur if rannum=0.
c
       I=max(Io, nn)  !3*I is the candidate eddy size
       return
       end
