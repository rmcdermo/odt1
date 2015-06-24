      subroutine Brng(rannum,i1,i2)
c
c     This subroutine generates pseudo-random numbers uniformly in 
c     the exclusive interval [0,1] and was supplied by Tom Dreeben.
c
c     November 1990
c
c     method: see f. james, computer physics communications, vol.60
c     p.329, 1990.
c
      integer i1, i2
      double precision rannum
c
      ik=i1/53668
      i1=40014*(i1-ik*53668)-ik*12211
      if (i1.lt.0) i1=i1+2147483563
c
      ik=i2/52774
      i2=40692*(i2-ik*52774)-ik*3791
      if (i2.lt.0) i2=i2+2147483399
c
      ix=i1-i2
      if (ix.lt.1) ix=ix+2147483562
c
c
c     change by Alan Kerstein on 9/5/06:
c
c     The line commented out below is the single-precision
c     prescription of f. james.  The line below it gives a 
c     closer approximation to sampling from the interval 
c     [0,1] for a double precision variable.
c
c      rannum=float(ix)*4.656613e-10
c
      rannum=float(ix)*4.65661305956e-10
c
      return
      end
