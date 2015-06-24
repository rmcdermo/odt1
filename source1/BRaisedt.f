       subroutine BRaisedt(Np,dt,td,pa,ii,ipars,rpars)
       implicit none
       integer Np, ii, iwait
       integer ipars(100)
       double precision dt, td, pa, pmin, dtfac, tdfac
       double precision rpars(100)
c
c If acceptance probabilities averaged over Np values, pa,  
c fall below a specified min, pmin, then increase dt 
c enough so that the average will be pmin (unless flow 
c statistic change), but do not increase dt by more than 
c the factor dtfac.  Then adjust td to maintain the 
c specified ratio, tdfac, of td to dt.
c
       iwait=ipars(2)  !number of eddy trials before dt adjustment
       if (mod(ii,iwait) .ne. 0) return  !adjust now?
c
c During eddy trials, pa is incremented only if the acceptance 
c probability is non-zero.  There are Np such values, where 
c 0.leq.Np.leq.ii .
c
       if (Np .gt. 0) pa=pa/(1.d0*Np)  !convert sum to average
       pmin=rpars(2)  !aim for this average acceptance probability
       dtfac=rpars(3)  !maximum dt increase factor
       tdfac=rpars(4)  !td/dt (enforce a fixed ratio)
       if (pa .lt. pmin) then
        if (pa .lt. pmin/dtfac) then
         dt=dt*dtfac !don't increase dt by more than dtfac
        else
         dt=dt*pmin/pa  !increase dt to target value
        endif
       endif
       td=tdfac*dt  !maintain specified proportionality td/dt
       pa=0.d0  !reinitialize pa
       Np=0  !reinitialize count of non-zero terms in pa sum
       return
       end
