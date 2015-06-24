       subroutine BSetOld(N,M,L,u,v,w,old,NVAL)
       implicit none
       integer N, M, L, NVAL, j
       double precision u(N), v(N), w(N), old(NVAL,3)
c
c Property values in the range of the eddy about 
c to be implemented are saved in the array old.  
c It is dimensioned using NVAL rather than N so 
c the array structure will conform to its 
c specification in the main program.
c
       do j=M, M+L-1
        old(j,1)=u(j)
       enddo
       do j=M, M+L-1
        old(j,2)=v(j)
       enddo
       do j=M, M+L-1
        old(j,3)=w(j)
       enddo
       return
       end
