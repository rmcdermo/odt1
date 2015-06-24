       subroutine XRecord(ifile,N,x,s)
       implicit none
       integer ifile, N, j
       double precision s(N), h(N), x(N)
 100   format(2E15.7)
       do j=1, N
        h(j)=s(j)
        if(abs(h(j)) .lt. 0.1d0**30) h(j)=0.d0	!prevents format problem
        write(ifile,100) x(j),h(j)
       enddo
       call flush(ifile)
       return
       end
