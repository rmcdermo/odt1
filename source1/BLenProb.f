       subroutine BLenProb(N,Io,Ip,PL,Co,Cv,ipars,nipars)
       implicit none
       integer N, Io, Ip, Iv, I, ipars(100), nipars
       double precision Co, Cv, PL(N), C ,z
c
c Compute the probability distribution from which I is sampled, 
c where 3*I is the candidate eddy size.  Here the eddy size is the 
c integer number of mesh cells in the eddy.  I corresponds to the 
c size of each triplet-map image of the original eddy interval.
c
c The allowed I values are Io through Iv.  The probability of an 
c allowed I value is taken to be the integral from I to I+1 of 
c f(s) = 2 Ip C s^{-2} exp(-2 Ip/s), where C is chosen so that the 
c probabilities of the allowed I values sum to unity.  The 
c quantity z is the value of this integral for C=1.  Therefore 
c the normalization factor C is assigned to be the inverse of 
c the sum of z over the allowed I values.
c
c Set default values of distribution parameters
c
      Io=2
      Ip=4
      Iv=(N-1)/3 !integer truncation, so 3*Iv <= N-1 (see BSampleEddy.f)
c
c Read alternate values if provided in BPars.dat
c
      if(nipars.ge.3) Io=ipars(3)
      if(nipars.ge.4) Ip=ipars(4)
      if(nipars.ge.5) Iv=min(Iv,ipars(5))  !cannot exceed default value
c
c Compute quantities used in BLength to find a good starting point 
c for searching the PL array when sampling an I value.
c
       Co=dexp(-2.d0*Ip/(1.d0*Io))
       Cv=dexp(-2.d0*Ip/(1.d0*Iv))
c
c Compute the normalization factor.
c
       C=0.d0
       do I=Io, Iv
        z=dexp(-2.d0*Ip/(1.d0*I))*(dexp(2.d0*Ip/(I*(I+1.d0)))-1.d0)
        C=C+z
       enddo
       C=1.d0/C
c
c Compute the cumulative probability distribution of I values.
c
       do I=1, Io-1
        PL(I)=0.d0
       enddo
       do I=Io, Iv-1
        z=dexp(-2.d0*Ip/(1.d0*I))*(dexp(2.d0*Ip/(I*(I+1.d0)))-1.d0)
        PL(I)=PL(I-1)+C*z
       enddo
       do I=Iv, N
        PL(I)=1.d0  !use exact value to avoid numerical error
       enddo
       return
       end
