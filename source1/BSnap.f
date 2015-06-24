       subroutine BSnap
     &  (N,u,v,w,dom,visc,istat,edstat,cstat,ioptions,
     &  NVAL,mVAL,nstatVAL)
       implicit none
       integer NVAL, mVAL, nstatVAL
       integer N, i, j, k, ii, istat, ifile, ioptions(100), N1
       double precision dom, visc, width, fluxfac, fluxfac2, zero(N)
       double precision u(N), v(N), w(N)
       double precision temp(N), tv(N), ta(N), p(N), d(N), bal(N)
       double precision ht(N), empty(N)
       double precision edstat(NVAL,4,4,nstatVAL), eavg(NVAL,4,3)
       double precision cstat(NVAL,mVAL,nstatVAL), cavg(NVAL,mVAL)
       zero(1)=0.d0
       N1=1
       width=dom/dfloat(N) !mesh spacing
       fluxfac=visc/width
       fluxfac2=visc/width**2
c
c cstat arguments: index of mesh cell, index of quantity gathered in  
c subroutine BStats, index of statistics-gathering time interval
c
c edstat arguments: index of mesh cell, velocity moment (1 or 2), index
c of velocity component, index of statistics-gathering time interval
c

c to form averages from sums, divide by the data-collection time
       do j=2, 10		
        do i=1, N
         cavg(i,j)=cstat(i,j,istat)/cstat(i,1,istat) 
        enddo
       enddo

c divide eddy-induced changes by elapsed time to obtain rates of change
       do k=1, 3
        do j=1, 4
         do i=1, N
          eavg(i,j,k)=edstat(i,j,k,istat)/cstat(i,1,istat)
         enddo
        enddo
       enddo

c
c compute terms of the budget of TKE = (<u'^2>+<v'>^2+<w'>^2)/2
c
c begin viscous transport term
       i=2          !form (u')^2
       do j=1, N
        temp(j)=cavg(j,i+3)-cavg(j,i)**2
       enddo
       do i=3, 4    !add (v')^2 and (w')^2
        do j=1, N
         temp(j)=temp(j)+cavg(j,i+3)-cavg(j,i)**2
        enddo
       enddo
       do j=2, N-1
        tv(j)=0.5d0*fluxfac2*(temp(j+1)+temp(j-1)-2.d0*temp(j))
       enddo
       tv(1)=0.5d0*fluxfac2*(temp(2)-2.d0*temp(1))
       tv(N)=tv(N-1)
c end viscous transport term
c begin viscous dissipation term (this is obsolete; see alternate method below)
       i=8
       k=2
       j=1
       temp(j)=cavg(j,i)-cavg(j,k)**2
       do j=2, N !form <du^2>-<du>^2
        temp(j)=cavg(j,i)-(cavg(j,k)-cavg(j-1,k))**2
       enddo
       do ii=1, 2
        i=i+1
        k=k+1
        j=1
        temp(j)=temp(j)+cavg(j,i)-cavg(j,k)**2
        do j=2, N !add <dv^2>-<dv>^2 and <dw^2>-<dw>^2
         temp(j)=temp(j)+cavg(j,i)-(cavg(j,k)-cavg(j-1,k))**2
        enddo
       enddo
       do j=1, N !dimensionalize the viscous dissipation
        d(j)=fluxfac2*temp(j)
       enddo
c end viscous dissipation term
c begin shear production term
       i=1 !u contribution
       temp(N)=.5d0*eavg(N,1,i)
       do j=N-1,1,-1
        temp(j)=temp(j+1)+.5d0*(eavg(j,1,i)+eavg(j+1,1,i))
       enddo
       do j=2,N-1
        p(j)=-.5d0*(cavg(j+1,i+1)-cavg(j-1,i+1))*temp(j)
       enddo
        p(1)=-.5d0*cavg(2,i+1)*temp(1)
       do i=2, 3 !v and w contributions
        temp(N)=.5d0*eavg(N,1,i)
        do j=N-1,1,-1
         temp(j)=temp(j+1)+.5d0*(eavg(j,1,i)+eavg(j+1,1,i))
        enddo
        do j=2,N-1
         p(j)=p(j)-.5d0*(cavg(j+1,i+1)-cavg(j-1,i+1))*temp(j)
        enddo
        p(1)=p(1)-.5d0*cavg(2,i+1)*temp(1)
       enddo
       p(N)=0.d0
c
c Divide by cell width to convert velocity differences to velocity derivatives 
c and multiply by cell width to convert sums of velocity changes to integrals 
c of velocity changes.  These operations cancel so p is correctly 
c dimensioned.
c
c end shear production term
c begin advective transport term
       i=1 !u contribution
       do j=1, N
        ta(j)=0.5d0*eavg(j,2,i)
     &   -cavg(j,i+1)*eavg(j,1,i)-p(j)
       enddo
       do i=2, 3 !v and w contributions
        do j=1, N
         ta(j)=ta(j)+0.5d0*eavg(j,2,i)
     &    -cavg(j,i+1)*eavg(j,1,i)
        enddo
       enddo
c end advective transport term (no factors of cell width required)

c
c alternate evaluation of viscous dissipation:
c Viscous dissipation of TKE is the mean rate of reduction of total KE 
c minus the mean rate of reduction of mean KE plus the viscous transport 
c contribution to the rate of TKE change.  This method is adopted as 
c standard because it gives much better balance of terms than the other 
c method, which is subject to the order-of-accuracy limitation of 
c advective advancement, in particular the coupling of viscous and 
c pressure-gradient effects.
c
       i=1 !u contribution
c       if(i.eq.0) then  !bypass alternate computation
       do j=1, N
        d(j)=-0.5d0*eavg(j,4,i)
     &   +cavg(j,i+1)*eavg(j,3,i)+tv(j)
       enddo
       do i=2, 3 !v and w contributions
        do j=1, N
        d(j)=d(j)-0.5d0*eavg(j,4,i)
     &   +cavg(j,i+1)*eavg(j,3,i)
        enddo
       enddo
c       endif

c
c Compute output statistics and write each record.
c
       do j=1, N
        empty(j)=0.d0 !use to fill records not used
        ht(j)=dom*(1.d0*j)/(1.d0*N) !data locations
       enddo
      if(ioptions(1).eq.0) then !intercomparison file format

       ifile=37+4*istat		!data sets A1 - A4
       N=-N			!flags write header
       call BRecord(ifile,N,ht)	!write file header and heights
       i=2
       do j=1, N			!mean u velocity
        temp(j)=cavg(j,i)
       enddo
        call BRecord(ifile,N,temp)


       ifile=38+4*istat			!data sets B1 - B4
       N=-N				!flags write header
       call BRecord(ifile,N,ht)	!write file header and heights
       do i=2, 4			!velocity variances
        do j=1, N
         temp(j)=(cavg(j,i+3)-cavg(j,i)**2)
        enddo
        call BRecord(ifile,N,temp)
       enddo
 
       ifile=39+4*istat			!data sets C1 - C4
       N=-N				!flags write header
       call BRecord(ifile,N,ht)	!write file header and heights	
       temp(N)=width*eavg(N,1,1)	!u advective flux
       do j=N-1,1,-1
        temp(j)=temp(j+1)+width*eavg(j,1,1)
       enddo
       call BRecord(ifile,N,temp)

       ifile=40+4*istat			!data sets D1 - D4
       N=-N				!flags write header
       call BRecord(ifile,N,ht)	!write file header and heights
       call BRecord(ifile,N,p)  !shear production
       call BRecord(ifile,N,ta) !advective transport
       call BRecord(ifile,N,tv) !viscous transport
       call BRecord(ifile,N,d)  !dissipation

      endif

      if(ioptions(1).eq.1) then !xmgrace file format

       ifile=31+10*istat		!data sets A1 - A4
       call XRecord(ifile,N1,zero,zero)	!write first line
       i=2
       do j=1, N			!mean u velocity
        temp(j)=cavg(j,i)
       enddo
        call XRecord(ifile,N,ht,temp)

       ifile=32+10*istat		!data sets B1 - B4
       call XRecord(ifile,N1,zero,zero)	!write first line
       do i=2, 4			!velocity variances
        do j=1, N
         temp(j)=(cavg(j,i+3)-cavg(j,i)**2)
        enddo
        call XRecord(ifile,N,ht,temp)
       enddo
 
       ifile=33+10*istat		!data sets C1 - C4
       call XRecord(ifile,N1,zero,zero)	!write first line
       temp(N)=width*eavg(N,1,1)	!u advective flux
       do j=N-1,1,-1
        temp(j)=temp(j+1)+width*eavg(j,1,1)
       enddo
       call XRecord(ifile,N,ht,temp)

       ifile=34+10*istat		!data sets D1 - D4
       call XRecord(ifile,N1,zero,zero)	!write first line
       call XRecord(ifile,N,ht,p)  !shear production

       ifile=35+10*istat		!data sets E1 - E4
       call XRecord(ifile,N1,zero,zero)	!write first line
       call XRecord(ifile,N,ht,ta) !advective transport

       ifile=36+10*istat		!data sets F1 - F4
       temp(1)=tv(N)
       call XRecord(ifile,N1,zero,temp)	!wall values forced to be equal 
       call XRecord(ifile,N,ht,tv) !viscous transport

       ifile=37+10*istat		!data sets G1 - G4
       temp(1)=d(N)
       call XRecord(ifile,N1,zero,temp)	!wall values forced to be equal 
       call XRecord(ifile,N,ht,d)  !dissipation

      endif

      ifile=80+istat

      do j=1, N
       bal(j)=ta(j)+tv(j)+p(j)-d(j)  !zero if exact steady-state balance
      enddo

      call XRecord(ifile,N1,zero,zero)	!write first line
      call XRecord(ifile,N,ht,bal)      !balance

c
c write size distributions of eddy counts
c
      ifile=90+istat

      do j=1,N-1
       write(ifile,*) j,edstat(j,1,4,istat),edstat(j,2,4,istat)
      enddo
      return
      end
