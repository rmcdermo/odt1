       function BsKd(N,M,L,f)
       implicit none
       integer N, M, L, Lo, j, j1, j2, j3
       double precision BsKd, sum, y1, y2, y3
       double precision s(N), f(N)
c
c Compute sK based on the discrete-map definition of
c triplet-map displacements.  sK is defined as (1/L)^2 
c times the integral of (s' times the displacement 
c distance), where L is the map size and ' denotes the 
c triplet-mapped array.  This is coded in length units 
c for which the size of an array cell is unity.  This 
c differs from the units of the code, but is consistent 
c because the integral and the prefactor (1/L)^2 have 
c units of (length^2 times units of s) and (1/length)^2, 
c so the units reduce to the units of s regardless of 
c the length unit used in the computation.
c
c Work with the mapped portion of the input array f, 
c which is written to array s.
c
       do j=M, M+L-1  !triplet map is applied to cells [M,M+L-1]
        s(j)=f(j)  
       enddo
       call BTriplet(N,M,L,s)  !apply triplet map to s
       Lo=L/3  !number of cells in a mapped image
       sum=0.d0  !initialize integration
       do j=1, Lo  !loop over cells of an image
c
c displacements to respective map images
c
        y1=-2.d0*(j-1)
        y2=4.d0*(j+Lo-1)-2.d0*(L-1)
        y3=2.d0*(L-1)-2.d0*(j+Lo+Lo-1)
c
c indices of locations of cells mapped to the respective images
c
        j1=M+j-1
        j2=M+j+Lo-1
        j3=M+j+Lo+Lo-1
c
c To integrate, add the contribution of each image to the sum.
c
        sum=sum+s(j1)*y1
        sum=sum+s(j2)*y2
        sum=sum+s(j3)*y3
       enddo
       BsKd=sum/dfloat(L*L)  !final result for sK
       return
       end
