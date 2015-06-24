       subroutine BStats
     &  (N,u,v,w,dom,visc,cstat,tstep,NVAL,mVAL,nstatVAL,istat)
       implicit none
       integer N, i, j, k, NVAL, mVAL, nstatVAL, istat
       double precision u(N), v(N), w(N)
       double precision dom, visc, tstep, fluxfac, fluxfac2
       double precision cstat(NVAL,mVAL,nstatVAL), term(N)
c
c Increment cstat array with quantities used to construct
c output statistics.
c
       do j=1, N  !time interval
        term(j)=tstep
       enddo
       i=1
       call BAddTerm(i,N,cstat,term,NVAL,mVAL,nstatVAL,istat)  !increment array cstat
       do j=1, N
        term(j)=u(j)*tstep
       enddo
       i=2
       call BAddTerm(i,N,cstat,term,NVAL,mVAL,nstatVAL,istat)  !increment array cstat
       do j=1, N
        term(j)=v(j)*tstep
       enddo
       i=3
       call BAddTerm(i,N,cstat,term,NVAL,mVAL,nstatVAL,istat)  !increment array cstat
       do j=1, N
        term(j)=w(j)*tstep
       enddo
       i=4
       call BAddTerm(i,N,cstat,term,NVAL,mVAL,nstatVAL,istat)  !increment array cstat
       do j=1, N
        term(j)=u(j)*u(j)*tstep
       enddo
       i=5
       call BAddTerm(i,N,cstat,term,NVAL,mVAL,nstatVAL,istat)  !increment array cstat
       do j=1, N
        term(j)=v(j)*v(j)*tstep
       enddo
       i=6
       call BAddTerm(i,N,cstat,term,NVAL,mVAL,nstatVAL,istat)  !increment array cstat
       do j=1, N
        term(j)=w(j)*w(j)*tstep
       enddo
       i=7
       call BAddTerm(i,N,cstat,term,NVAL,mVAL,nstatVAL,istat)  !increment array cstat
       term(1)=u(1)**2*tstep
       do j=2, N
        term(j)=(u(j)-u(j-1))**2*tstep
       enddo
       i=8
       call BAddTerm(i,N,cstat,term,NVAL,mVAL,nstatVAL,istat)  !increment array cstat
       term(1)=v(1)**2*tstep
       do j=2, N
        term(j)=(v(j)-v(j-1))**2*tstep
       enddo
       i=9
       call BAddTerm(i,N,cstat,term,NVAL,mVAL,nstatVAL,istat)  !increment array cstat
       term(1)=w(1)**2*tstep
       do j=2, N
        term(j)=(w(j)-w(j-1))**2*tstep
       enddo
       i=10
       call BAddTerm(i,N,cstat,term,NVAL,mVAL,nstatVAL,istat)  !increment array cstat
       return
       end
