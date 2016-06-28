c\BeginDoc
c
c\Name: pslarnv
c
c Message Passing Layer: MPI
c
c\Description:
c
c  Parallel Version of ARPACK utility routine slarnv
c
c  PSLARNV returns a vector of n (nloc) random real numbers from a uniform or
c  normal distribution. It is assumed that X is distributed across a 1-D array 
c  of processors ( nprocs < 1000 )
c
c\Arguments
c  COMM    MPI Communicator for the processor grid
c
c  IDIST   (input) INTEGER
c          Specifies the distribution of the random numbers:
c          = 1:  uniform (0,1)
c          = 2:  uniform (-1,1)
c          = 3:  normal (0,1)
c
c  ISEED   (input/output) INTEGER array, dimension (4)
c          On entry, the seed of the random number generator; the array
c          elements must be between 0 and 4095, and ISEED(4) must be
c          odd.
c          On exit, the seed is updated.
c
c  N       (input) INTEGER
c          The number of random numbers to be generated.
c
c  X       (output) Real array, dimension (N)
c          The generated random numbers.
c
c\Author: Kristi Maschhoff
c
c\Details
c
c  Simple parallel version of LAPACK auxiliary routine slarnv 
c  for X distributed across a 1-D array of processors.
c  This routine calls the auxiliary routine SLARNV to generate random
c  real numbers from a uniform (0,1) distribution. Output is consistent
c  with serial version. 
c
c\SCCS Information: 
c FILE: larnv.F   SID: 1.1   DATE OF SID: 1/23/96   
c
c-----------------------------------------------------------------------
c
      subroutine pslarnv( comm, idist, iseed, n, x )
c
      include  'mpif.h'
c
c     .. MPI VARIABLES AND FUNCTIONS ..
      integer   comm, myid, numprocs, ierr
c     ..
c     .. Scalar Arguments ..
      integer			i, n
c     ..
c     .. Array Arguments ..
      integer			iseed( 4 )
      Real			
     &                  x( * )
c     ..
c     .. Local Array Arguments ..
      integer			dist1(1000), dist2(1000)
c     ..
c     .. External Subroutines ..
      external			slarnv
c     ..
c     .. Executable Statements ..
c
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, numprocs, ierr )
c
	  do 10 i=1,numprocs
         dist1(i) = 0
 10   continue
c
      dist1( myid+1 ) = n
      call MPI_ALLREDUCE( dist1(1), dist2(1), numprocs, 
     &                    MPI_INTEGER, MPI_SUM, comm, ierr )
c
      do 15 i = 1,myid + 1 
         call slarnv( idist, iseed, dist2(i), x )
 15   continue
c
      return
      end
