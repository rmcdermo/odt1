       subroutine BSeries(NVAL,N,itime,u,umoms)
c
c Accumulate centerline values of u and u**2
c (for odd N, not exactly at centerline).
c
       implicit none
       integer NVAL, N, itime, ictr
       double precision u(NVAL), umoms(2,10000)
       ictr=0.5d0*N
       if(itime.le.10000) then
        umoms(1,itime)=umoms(1,itime)+u(ictr)
        umoms(2,itime)=umoms(2,itime)+u(ictr)**2
       endif
       return
       end
