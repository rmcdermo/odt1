      subroutine BSeeds(ipars,i1,i2)
      implicit none
      integer isarg,i1,i2,i
      integer is(5), ipars(100)
      double precision rannum
      data is/984573405,238400857,757234405,430904431,175728293/
c
c Assign one of five possible pairs of initial random number seeds 
c based on the input isarg.
c
      isarg=ipars(6)
      i1=is(isarg)
      i2=is(1+mod(isarg,5))
c
c To avoid possible artifacts of initial seed choice, 
c warm up the random number generator.
c
      do i=1, 100
       call Brng(rannum,i1,i2)
      enddo
      return
      end
