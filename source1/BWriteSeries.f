       subroutine BWriteSeries(niter,itime,tend,ntseg,umoms,ioptions)
       implicit none
       integer niter, itime, i, j, ifile, num, ntseg
       integer ioptions(100)
       double precision umoms(2,1000)
       double precision time(1000), temp(1000), tend, tfac
 100   format(E15.7)
 200   format(2E15.7)
       ifile=40
       tfac=tend/float(itime)
       do j=1, itime  !time for which output is computed
        time(j)=j*tfac
       enddo
       do j=1, itime  !variance of u at centerpoint
        temp(j)=umoms(2,j)/float(niter)-(umoms(1,j)/float(niter))**2
       enddo
       if(ioptions(1).eq.0) then
        write(ifile,*) itime
        write(ifile,100) (time(j),j=1,itime)
        write(ifile,100) (temp(j),j=1,itime)
       endif
       if(ioptions(1).eq.1) then  !xmgrace file format
        do j=1, itime
         write(ifile,200) time(j),temp(j)
        enddo
       endif
       return
       end
