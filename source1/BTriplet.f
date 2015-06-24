       subroutine BTriplet(N,M,L,s)
       implicit none
       integer N, M, L, Lo, j, k
       double precision s(N), x(N)
c
c Permute cells M through M+L-1 of the array s 
c as prescribed by the discrete triplet map, where 
c L is an integer multiple of 3.
c
       Lo=L/3  !number of cells in each map image
       do j=1, Lo
        k=M+3*(j-1)
        x(j)=s(k)  !gather the cells going to the 1st image
       enddo
       do j=1, Lo
        k=M+L+1-(3*j)  !minus sign because second image is flipped
        x(j+Lo)=s(k)  !gather the cells going to the 2nd image
       enddo
       do j=1, Lo
        k=M+(3*j)-1
        x(j+Lo+Lo)=s(k)  !gather the cells going to the 3rd image
       enddo
       do j=1, L
        k=M+j-1
        s(k)=x(j)  !write the mapped values to the original array
       enddo
       return
       end
