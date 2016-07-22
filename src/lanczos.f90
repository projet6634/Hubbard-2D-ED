module lanczos


    implicit none

contains

    subroutine lanczos_diagonalize(nstep, a, b, ev, v)
        INCLUDE 'mkl_lapack.fi'
        integer, intent(in) :: nstep
        double precision, intent(in) :: a(nstep), b(nstep)
        double precision, intent(out) :: ev(nstep), v(nstep,nstep)

        double precision :: d(nstep), e(nstep-1), abstol, &
                            work(5*nstep)
        integer :: il, iu, m, iwork(5*nstep), ifail(nstep), &
                   info
        double precision :: w(nstep), z(nstep,nstep)

        integer :: i

        ! copy a
        d = a
        e(1:nstep-1) = b(2:nstep)

        il = 1     ! lowest eigenvalue index
        iu = nstep ! highest eigenvalue index
        abstol = 2*DLAMCH('S') ! recommended tolerance

        call dstevx( 'V', 'A', nstep, d, e, 0.d0, 0.d0, il, iu, abstol, &
                     m, w, z, nstep, work, iwork, ifail, info )

        if (info/=0) then
            write(*,*) "dstevx error. info = ", info
            if (info>0) then
                write(*,*) "ifail = "
                write(*,*), ifail
            endif
            stop
        endif

        ev(1:m) = w(1:m)
        v(1:nstep,1:m) = z(1:nstep,1:m) 
    end subroutine 

    ! calculates matrix elements of matrix M in the lanczos basis.
    ! M  =  ( a1  b2  0   0   ... 0   ) 
    !       ( b2  a2  b3  0   ... 0   )
    !       ( 0   b3  a3  b4  ... 0   )
    !       ( 0   0   b4  a4  ... 0   )
    !       ( ...             ... bn  )
    !       ( ...             bn  an  )
    !
    ! uses an external routine matmult that handles the matrix-vector product.
    ! refer the interface declared inside the subroutine.
    subroutine lanczos_iteration(matmult, nloc, v_init, maxnstep, nstep, a, b)
        use numeric_utils, only: mpi_dot_product, mpi_norm
        integer, intent(in) :: &
            nloc, &       ! dimension of the vector local to the node
            maxnstep         ! maximum number of iteration steps

        double precision, intent(in) :: &
            v_init(nloc)     ! starting vector for lanczos iteration

        integer, intent(out) :: &
            nstep         ! number of calculated steps <= maxnstep
                          ! if b(i) = 0 for some i, nstep = i-1
        double precision, intent(out) :: &
            a(maxnstep), & ! diagonal matrix element in lanczos basis
            b(maxnstep)    ! off-diagonal matrix element in lanczos basis

        interface 
            ! external subroutine for matrix multiplication
            ! Y = M*X
            subroutine matmult(n,x,y)
                integer, intent(in) :: n
                double precision, intent(in) :: x(n)
                double precision, intent(out) :: y(n)
            end subroutine matmult
        end interface

        ! temporary lanczos vectors
        double precision :: v(nloc,2), w(nloc)
        double precision :: norm_v
        
        integer :: j, ierr
        
        ! Lanczos steps
        ! ref: https://en.wikipedia.org/wiki/Lanczos_algorithm#Iteration 

        ! normalize the initial vector
        
        norm_v = mpi_norm( v_init, nloc)
        v(:,2) = v_init/norm_v
        v(:,1) = 0.0D0

        a(maxnstep) = 0.0D0
        b(maxnstep) = 0.0D0

        ! v(:,1) = v_(j-1)
        ! v(:,2) = v_j
        ! w(:)   = w_j
        lanczos_loop: do j=1,maxnstep-1
            nstep = j

            ! w_j = H*v_j
            call matmult( nloc, v(:,2), w(:) )

            ! a_j = dot(w_j,v_j)
            a(j) = mpi_dot_product(w(:), v(:,2), nloc)

            ! w_j = w_j - a_j * v_j - b_j * v_(j-1)
            w(:) = w(:) - a(j)*v(:,2) - b(j)*v(:,1)

            ! b_(j+1) = norm(w_j)
            b(j+1) = mpi_norm(w(:), nloc)

            if (b(j+1).lt.1.d-8) then
                exit lanczos_loop
            endif

            ! v_(j-1) = v_j
            v(:,1) = v(:,2)

            ! v_(j+1) = w_j/b(j+1)
            v(:,2) = w(:)/b(j+1)
        enddo lanczos_loop

        ! handles the last step
        if (nstep.eq.maxnstep-1) then
            call matmult( nloc, v(:,2), w(:) )
            a(maxnstep) = mpi_dot_product(w(:),v(:,2),nloc)
            nstep = maxnstep
        endif
    end subroutine lanczos_iteration
end module lanczos
