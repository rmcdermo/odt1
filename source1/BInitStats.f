       subroutine BInitStats
     &  (N,NVAL,mVAL,nstatVAL,edstat,cstat,umoms)
c
c Initialize the arrays that accumulate statistics over
c all realizations by setting all values equal to zero.
c
       implicit none
       integer N, i, j, k, ii, NVAL, mVAL, nstatVAL
       double precision cstat(NVAL,mVAL,nstatVAL)
       double precision edstat(NVAL,2,4,nstatVAL)
       double precision umoms(2,10000)
       do k=1, nstatVAL
        do j=1, 4
         do ii=1, 2
          do i=1, N
           edstat(i,ii,j,k)=0.d0
          enddo
         enddo
        enddo
       enddo
       do k=1, nstatVAL
        do j=1, mVAL
         do i=1, N
          cstat(i,j,k)=0.d0
         enddo
        enddo
       enddo
       do k=1, 10000
        do i=1, 2
         umoms(i,k)=0.d0
        enddo
       enddo
       return
       end
