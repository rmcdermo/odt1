Files that are subroutines or functions are named [name].f , 
where [name] is the name of the subroutine or function .  Input 
files are of the form [name].dat .  Names of files used for  
the default version begin with B, except for XRecord.f.  Other 
files might be needed to compile other versions of the code but 
the default version does not use them.

File *[name].dat is read by subroutine *Read[name] , where * 
indicates that a B appears at the beginning of the 
subroutine or file name in the default version.

Inputs specifying the flow configuration are read from BConfig.dat 
by subroutine BReadConfig.  These inputs are defined in the comments 
in BReadConfig.f .

Case options: The option index is the argument k of the array 
ioptions(k) .  The option value 0 is the default.  Subroutine 
BReadOptions reads the option values from file BOptions.dat .

OPTION          OPTION          MEANING
INDEX           VALUE

1		0		output in intercomparison (ABL) format
		1		output in xmgrace format

2		0		vector (3 component) ODT
		1		1 component ODT

3		0		no scalar fields
		n>0		n scalar fields

4		0		bottom velocity B.C.s: Dirichlet
		1		bottom velocity B.C.s: Neumann

5		0		top velocity B.C.s: Dirichlet
		1		top velocity B.C.s: Neumann

6		0		bottom B.C.s of scalars: Dirichlet
		1		bottom scalar B.C.s of scalars: Neumann

7		0		top scalar B.C.s of scalars: Dirichlet
		1		top scalar B.C.s of scalars: Neumann

8		0		all B.C.s = 0. (value or derivative)
		1		read B.C.s from an input file

9		0		all initial profiles are nominally zero
		1		initial profiles linear between endpoints 
				  that are read from a file
		2		initial profiles are read from a file 
				  (e.g., restart file)
		3		initial profiles are computed in a subroutine, 
				  which might read inputs from a file

10		0		no buoyancy
		1		buoyancy; requires n>0, where the first scalar 
				  is taken to be density or a density surrogate

11		0		no Coriolis force
		1		Coriolis force



In the main program, Bodt , the array size of property fields 
is hardwired as 100000 by a parameter statement.  The input 
value of N must not exceed this.  The hardwired value can be 
changed as needed to suit the application.

Structure of input file BPars.dat (algorithm parameters):

First line:	number of integer parameters (one per each line that follows)

	Integer parameters:

INDEX		MEANING

1		N = number of mesh cells for this run (<= NVAL)

2		number of eddy trials before an attempt to increase dt

3		number of mesh cells in one image of the smallest eddy: >=  2, <= N/3

4		number of mesh cells in one image of the most probable eddy: >= 2, <= N/3

5		number of mesh cells in one image of the largest eddy: >= 2, <= N/3
                (in code, defaults to largest integer <= N/3 and truncates any larger 
		value that is specified)

6		isarg = index for initializing random number seeds


Next line:	number of real parameters (one per each line that follows)

	Real parameters:

INDEX		MEANING

1		maximum allowed acceptance probability

2		target average acceptance probability

3		largest allowed multiplicative increase of dt at a given time

4		specified ratio td/dt

5               fraction of CFL time for which advancement is allowed

6		initial value of dt
                (negative means multiply the dt value computed in BInitIter by 
                the absolute value of this parameter)
