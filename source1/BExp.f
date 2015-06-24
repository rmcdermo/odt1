       function BExp(dt,i1,i2)
c
c Sample from an exponential distribution with mean dt
c by sampling uniformly in [0,1] to obtain a value of 
c the cumulative distribution and solving for the 
c corresponding value of the random variable.
c
       implicit none
       integer i1, i2
       double precision BExp, dt, rannum
       call Brng(rannum,i1,i2)
       BExp=-dt*dlog(1.d0-rannum)
       return
       end
