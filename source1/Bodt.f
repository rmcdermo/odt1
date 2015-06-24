      program Bodt
      implicit none
c
c Specify variables for default case.
c
      integer N, NVAL, nstatVAL, mVAL
      parameter(NVAL=10000)  !max number of mesh cells
      parameter(nstatVAL=4)  !max number of data-averaging periods
      parameter(mVAL=10)  !max number of properties gathered in cstat array
      integer ioptions(100), nopt  !from BOptions.dat
      integer ipars(100), nipars, nrpars  !from BPars.dat
      integer niter, nstat, ntseg  !from BConfig.dat
      integer iter, istat, itseg, itime  !time period counters

c variables associated with eddy sampling
      integer i1, i2  !random number seeds (updated during run)
      integer M, L, ii, Np, Io, Ip
      double precision rannum  !random number
      double precision Co, Cv, pa, pp, PL(NVAL)
      double precision uK, vK, wK

      double precision rpars(100)  !from BPars.dat
      double precision tend, dom, visc, pgrad, ratefac, viscpen  !from BConfig.dat
      double precision u(NVAL), v(NVAL), w(NVAL)  !default property fields
      double precision t, to, delt, dt, td, tmark, BExp

c data-gathering arrays
      integer jj
      double precision edstat(NVAL,4,4,nstatVAL), old(NVAL,3)
      double precision cstat(NVAL,mVAL,nstatVAL), umoms(2,10000)

c
c Read flags that specify any non-default options.
c
      call BReadOptions(ioptions,nopt)

      call BReadPars(ipars,rpars,nipars,nrpars)  !read algorithm parameters
      N=ipars(1)  !number of mesh cells for this run
      call BReadConfig  !read inputs needed for the default case
     & (niter,nstat,ntseg,tend,dom,visc,pgrad,ratefac,viscpen,N)
c
c output files
c
      if(ioptions(1).eq.0) then  !output in intercomparison format
       open(40, file="T.dat", status="unknown")
       open(41, file="A1.dat", status="unknown")
       open(42, file="B1.dat", status="unknown")
       open(43, file="C1.dat", status="unknown")
       open(44, file="D1.dat", status="unknown")
       open(45, file="A2.dat", status="unknown")
       open(46, file="B2.dat", status="unknown")
       open(47, file="C2.dat", status="unknown")
       open(48, file="D2.dat", status="unknown")
       open(49, file="A3.dat", status="unknown")
       open(50, file="B3.dat", status="unknown")
       open(51, file="C3.dat", status="unknown")
       open(52, file="D3.dat", status="unknown")
       open(53, file="A4.dat", status="unknown")
       open(54, file="B4.dat", status="unknown")
       open(55, file="C4.dat", status="unknown")
       open(56, file="D4.dat", status="unknown")
      endif
      if(ioptions(1).eq.1) then  !output in xmgrace format
       open(40, file="T.dat", status="unknown")
       open(41, file="A1.dat", status="unknown")
       open(42, file="B1.dat", status="unknown")
       open(43, file="C1.dat", status="unknown")
       open(44, file="D1.dat", status="unknown")
       open(45, file="E1.dat", status="unknown")
       open(46, file="F1.dat", status="unknown")
       open(47, file="G1.dat", status="unknown")
       open(51, file="A2.dat", status="unknown")
       open(52, file="B2.dat", status="unknown")
       open(53, file="C2.dat", status="unknown")
       open(54, file="D2.dat", status="unknown")
       open(55, file="E2.dat", status="unknown")
       open(56, file="F2.dat", status="unknown")
       open(57, file="G2.dat", status="unknown")
       open(61, file="A3.dat", status="unknown")
       open(62, file="B3.dat", status="unknown")
       open(63, file="C3.dat", status="unknown")
       open(64, file="D3.dat", status="unknown")
       open(65, file="E3.dat", status="unknown")
       open(66, file="F3.dat", status="unknown")
       open(67, file="G3.dat", status="unknown")
       open(71, file="A4.dat", status="unknown")
       open(72, file="B4.dat", status="unknown")
       open(73, file="C4.dat", status="unknown")
       open(74, file="D4.dat", status="unknown")
       open(75, file="E4.dat", status="unknown")
       open(76, file="F4.dat", status="unknown")
       open(77, file="G4.dat", status="unknown")
      endif
       open(81, file="H1.dat", status="unknown")
       open(82, file="H2.dat", status="unknown")
       open(83, file="H3.dat", status="unknown")
       open(84, file="H4.dat", status="unknown")
       open(91, file="I1.dat", status="unknown")
       open(92, file="I2.dat", status="unknown")
       open(93, file="I3.dat", status="unknown")
       open(94, file="I4.dat", status="unknown")
c
c Generate initial values of the random number seeds.
c
      call BSeeds(ipars,i1,i2)
c
c Initialize the flow state.
c
      call BInitRun(N,i1,i2,dom)
c
c Initialize the output variables.
c
      call BInitStats
     &  (N,NVAL,mVAL,nstatVAL,edstat,cstat,umoms)
c
c Construct the eddy size PDF.
c
      call BLenProb(N,Io,Ip,PL,Co,Cv,ipars,nipars)
c
c Loop over realizations.
c
      do iter=1,niter
c
c Initialize the default property fields (three velocity components), 
c counters, and time variables.
c
       call BInitIter
     & (N,u,v,w,itime,ii,Np,pa,tmark,t,to,dt,td,
     & nrpars,rpars,dom,visc,i1,i2)
c
c Loop over data-averaging periods within a realization.
c
       do istat=1,nstat  !loop over data-averaging periods
        write(11,*) "iter,istat= ",iter,istat
        call flush(11)
c
c Loop over sub-intervals within data-averaging periods.
c
        do itseg=1,ntseg
        tmark=tmark+tend/dfloat(nstat*ntseg)  !end time of the current sub-interval
        itime=itime+1  !integer index of the current sub-interval
c
c Perform eddy trials scheduled to occur before the end time, 
c tmark, of the current sub-interval.  During this loop, advance 
c other processes as needed.
c
         do while (t .le. tmark)
c
c The physical time to, corresponding to the cumulative time 
c advancement of processes other than eddies, lags (or equals)
c the eddy-trial time.  If the lag exceeds the threshold td,
c then advance the other processes to time t before sampling
c a candidate eddy.
c
          if ((t-to) .ge. td) then
           delt=t-to  !advancement time increment
           M=1
           L=N-1
           call BSetOld(N,M,L,u,v,w,old,NVAL) !gather pre-advancement statistics
           call BEqnStep  !implement the advancement
     &      (N,u,v,w,delt,dom,visc,pgrad,rpars,istat,cstat,
     &      NVAL,mVAL,nstatVAL)
c
c Gather statistics of changes during advancement.
c
           jj=2  !flags advancement statistics
           call BChange
     &      (N,M,L,u,v,w,old,edstat,NVAL,nstatVAL,istat,jj)
c
c update physical time
c
           to=t
          endif
c
c Sample a candidate eddy.
c
          call BSampleEddy
     &     (N,M,L,u,v,w,dt,td,PL,ratefac,viscpen,uK,vK,wK,
     &     pp,ii,pa,Np,rpars,Io,Ip,Cv,Co,i1,i2)
c
c Decide whether to accept the eddy.
c
          call Brng(rannum,i1,i2)
          if (rannum .lt. pp) then  !pp is the acceptance probability
c
c Copy the pre-eddy state.
c
           call BSetOld(N,M,L,u,v,w,old,NVAL)
c
c Implement the accepted eddy. 
c
           call BEddy(N,M,L,u,v,w,uK,vK,wK)
c
c Gather statistics of eddy-induced changes.
c
           jj=0  !flags eddy statistics
           call BChange
     &      (N,M,L,u,v,w,old,edstat,NVAL,nstatVAL,istat,jj)
c
c Time advance other processes to the nominal time 
c at which the eddy was sampled.
c

           delt=t-to  !advancement time increment
           M=1
           L=N-1
           call BSetOld(N,M,L,u,v,w,old,NVAL) !gather pre-advancement statistics
           call BEqnStep  !implement the advancement
     &      (N,u,v,w,delt,dom,visc,pgrad,rpars,istat,cstat,
     &      NVAL,mVAL,nstatVAL)
c
c Gather statistics of changes during advancement.
c
           jj=2  !flags advancement statistics
           call BChange
     &      (N,M,L,u,v,w,old,edstat,NVAL,nstatVAL,istat,jj)
c
c update physical time
c
           to=t
          endif
c
c Update the scheduled time of the next candidate eddy. 
c The thinning algorithm that is used to sample eddies involves 
c sampling at a rate 1/dt, with sampling time increments that are
c sampled from the exponential distribution with mean dt.
c
          t=t+BExp(dt,i1,i2)
c
c Increase the dt value used for future eddy trials if necessary.
c
          call BRaisedt(Np,dt,td,pa,ii,ipars,rpars)
         enddo  !end of loop over eddy trials within a sub-interval
c
c The next eddy trial will occur at t > tmark, which is in a
c later sub-interval.  The current physical time to is at most tmark.  
c Advance the physical time to tmark before gathering statistics 
c for the sub-interval.
c
         delt=tmark-to  !advancement time increment
         if(delt .gt. 0.d0) then
          M=1
          L=N-1
          call BSetOld(N,M,L,u,v,w,old,NVAL) !gather pre-advancement statistics
          call BEqnStep  !implement the advancement
     &     (N,u,v,w,delt,dom,visc,pgrad,rpars,istat,cstat,
     &     NVAL,mVAL,nstatVAL)
c
c Gather statistics of changes during advancement.
c
          jj=2  !flags advancement statistics
          call BChange
     &     (N,M,L,u,v,w,old,edstat,NVAL,nstatVAL,istat,jj)
c
c update physical time
c
          to=tmark
         endif
c
c Gather sub-interval statistics.
c
         call BSeries(NVAL,N,itime,u,umoms)
        enddo  !end of loop over sub-intervals
c
c During the last iteration, reduce the data 
c for the data-gathering period
c
        if(iter.eq.niter) call BSnap
     &  (N,u,v,w,dom,visc,istat,edstat,cstat,ioptions,
     &  NVAL,mVAL,nstatVAL)
       enddo  !end of loop over data-gathering periods
      enddo  !end of loop over iterations
      call BWriteSeries(niter,itime,tend,ntseg,umoms,ioptions)
      if(ioptions(1).eq.0) then
      close(40)
      close(41)
      close(42)
      close(43)
      close(44)
      close(45)
      close(46)
      close(47)
      close(48)
      close(49)
      close(50)
      close(51)
      close(52)
      close(53)
      close(54)
      close(55)
      close(56)
      endif

      if(ioptions(1).eq.1) then
      close(40)
      close(41)
      close(42)
      close(43)
      close(44)
      close(45)
      close(46)
      close(47)
      close(51)
      close(52)
      close(53)
      close(54)
      close(55)
      close(56)
      close(57)
      close(61)
      close(62)
      close(63)
      close(64)
      close(65)
      close(66)
      close(67)
      close(71)
      close(72)
      close(73)
      close(74)
      close(75)
      close(76)
      close(77)
      endif
      close(11)
      close(81)
      close(82)
      close(83)
      close(84)
      close(91)
      close(92)
      close(93)
      close(94)
      stop
      end

       include 'BAddK.f'
       include 'BAddTerm.f'
       include 'BAdv.f'
       include 'BChange.f'
       include 'BEddy.f'
       include 'BEqnStep.f'
       include 'BExp.f'
       include 'BInitIter.f'
       include 'BInitRun.f'
       include 'BInitStats.f'
       include 'BLength.f'
       include 'BLenProb.f'
       include 'BLowerdt.f'
       include 'BProb.f'
       include 'BRaisedt.f'
       include 'BReadConfig.f'
       include 'BReadOptions.f'
       include 'BReadPars.f'
       include 'BRecord.f'
       include 'Brng.f'
       include 'BrngGet.f'
       include 'BrngPut.f'
       include 'BSampleEddy.f'
       include 'BSeeds.f'
       include 'BSeries.f'
       include 'BSetOld.f'
       include 'BsKd.f'
       include 'BSnap.f'
       include 'BStats.f'
       include 'BTriplet.f'
       include 'BWriteSeries.f'
       include 'XRecord.f'
