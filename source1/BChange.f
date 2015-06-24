       subroutine 
     &  BChange(N,M,L,u,v,w,old,edstat,NVAL,nstatVAL,istat,jj)
       implicit none
       integer M, L, NVAL, nstatVAL, i, j, k, N, istat, idist
       integer jj, LL
       double precision u(N), v(N), w(N)
       double precision edstat(NVAL,4,4,nstatVAL), old(NVAL,3)
c
c Eddy-induced changes of first and second powers of 
c property values in the range of the eddy are saved 
c in the array edstat.  It is dimensioned using NVAL 
c and nstatVAL rather than N and nstat so the array 
c structure will conform to its specification in the 
c main program.
c
c jj is a flag equal to zero for eddy statistics, 
c 2 for BAdv statistics.
c
       do j=1,2
        i=1
        do k=M,M+L-1
         edstat(k,j+jj,i,istat)=edstat(k,j+jj,i,istat)+
     &     u(k)**j-old(k,i)**j
        enddo
        i=2
        do k=M,M+L-1
         edstat(k,j+jj,i,istat)=edstat(k,j+jj,i,istat)+
     &     v(k)**j-old(k,i)**j
        enddo
        i=3
        do k=M,M+L-1
         edstat(k,j+jj,i,istat)=edstat(k,j+jj,i,istat)+
     &     w(k)**j-old(k,i)**j
        enddo
       enddo
       if(jj.eq.2) return !no eddy statistics if applied to BAdv
c
c Accumulate counts of number of accepted eddies as a function 
c of eddy size, separated into eddies whose distance from a 
c wall is less than or equal to its size (i=4) and eddies that 
c do not obey this condition.
c
       idist=min(M,N-M-L+1)  !distance from wall
       j=1
       if(idist .gt. L) j=2
       LL=L/3
       edstat(LL,j,4,istat)=edstat(LL,j,4,istat)+1.d0
       return
       end
