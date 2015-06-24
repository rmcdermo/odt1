       subroutine BAddK(N,M,L,r,c)
       implicit none
       integer N, M, L, Lo, j, j1, j2, j3
       double precision c, y1, y2, y3
       double precision r(N)
c
c Add c*K to the array r.  For a given location,
c K is the array index of that location minus the 
c index of the location that is mapped to it by a 
c discrete triplet map applied to the index range 
c [M,M+L-1].  Here, K is a non-dimensional quantity 
c rather than the triplet-map displacement in units 
c of length (as it is defined physically).  The 
c eddy length L in the denominator of c (introduced 
c in the statement cfac=... in the calling routine 
c BEddy) is evaluated in the same units, giving 
c consistent evaluation of the quantity (K/L).
c
       Lo=L/3  !number of cells in a triplet-map image
       do j=1, Lo  !loop over cells in an image
c
c array index after map minus array index before map 
c for the jth cell of map images 1, 2, and 3, respectively
c
        y1=-2.d0*(j-1)
        y2=4.d0*(j+Lo-1)-2.d0*(L-1)
        y3=2.d0*(L-1)-2.d0*(j+Lo+Lo-1)
c
c index of the r array corresponding to the jth cell of 
c each of the three map images
c
        j1=M+j-1
        j2=M+j+Lo-1
        j3=M+j+Lo+Lo-1
c
c Increment u within each of the map images.
c
        r(j1)=r(j1)+(c*y1)
        r(j2)=r(j2)+(c*y2)
        r(j3)=r(j3)+(c*y3)
       enddo
       return
       end
