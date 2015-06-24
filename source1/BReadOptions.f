      subroutine BReadOptions(ioptions,nopt)
c
c Read integers that specify case options .
c File with options list is BOptions.dat .
c Default value of each option is zero .
c To run the default case, set nopt=0 .
c
      implicit none
      integer ioptions(100),k,nopt
c
c ioptions is the array of option flags .
c k is the running argument of ioptions .
c nopt is the number of options listed in BOptions.dat .
c BOptions.dat has one integer entry per line.
c The first entry is nopt, and the remaining nopt 
c entries are the option flags .
c
 200  format(i12)
      open(14, file="BOptions.dat", status="old")
      read(14, *) nopt  !number of options to read
      do k=1,100  !read the case options
       ioptions(k)=0
       if(k.le.nopt) read(14, *) ioptions(k)
      enddo
      close(14)
      return
      end
