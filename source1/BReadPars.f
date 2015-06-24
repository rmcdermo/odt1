      subroutine BReadPars(ipars,rpars,nipars,nrpars)  
c
c Read algorithm parameters.
c
      implicit none
      integer ipars(100), nipars, nrpars, i
      double precision rpars(100)
      open(13, file="BPars.dat", status="old")
      read(13, *) nipars  !number of integer parameters
      do i=1, 100
       ipars(i)=0
       if(i.eq.3) ipars(i)=2 !smallest eddy default size 6
       if(i.eq.4) ipars(i)=3 !most probable eddy default size 9
c Assign 1/3 of max eddy size the default value N.
c In BLenProb, 1/3 of max eddy size is truncated to the integer part of N/3.
       if(i.eq.5) ipars(i)=ipars(1)
       if(i.eq.6) ipars(i)=1  !default value of isarg
       if(i.le.nipars) read(13, *) ipars(i)  !read the integer parameters
      enddo
      read(13, *) nrpars  !number of real parameters
      do i=1, 100
       rpars(i)=0.d0
       if(i.eq.5) rpars(i)=0.5d0 !default timestep relative to CFL
       if(i.eq.6) rpars(i)=-1.d0 !forces default to computed initial dt
       if(i.le.nrpars) read(13, *) rpars(i)  !read the real parameters
      enddo
      close(13)
      return
      end
