       subroutine BLowerdt(dt,td,pp,pa,Np,rpars)
       implicit none
       integer Np
       double precision dt, td, tdfac,pp, pa, pmax
       double precision rpars(100)
c
c Reduce dt if the acceptance probability pp is too high.
c
       if (pp .eq. 0.d0) return  !no action if eddy is disallowed
       pmax=rpars(1)  !don't let pp exceed this value
       if (pp .gt. pmax) then  !pp is too high; reduce dt
        dt=dt*pmax/pp  !dt value for which pp=pmax
        tdfac=rpars(4)  !td/dt (enforce a fixed ratio)
        td=tdfac*dt  !maintain specified proportionality td/dt
        pp=pmax  !new pp value
       endif
c
c Gather statistics used in BRaisedt to decide whether 
c dt should be increased.
c
        pa=pa+pp  !accumulate non-zero pp values
        Np=Np+1  !number of non-zero pp values summed
       return
       end
