       subroutine BAddTerm(i,N,cstat,term,NVAL,mVAL,nstatVAL,istat)
       implicit none
       integer i, j, N, NVAL, mVAL, nstatVAL, istat
       double precision cstat(NVAL,mVAL,nstatVAL), term(N)
       do j=1, N
        cstat(j,i,istat)=cstat(j,i,istat)+term(j)
       enddo
       return
       end
