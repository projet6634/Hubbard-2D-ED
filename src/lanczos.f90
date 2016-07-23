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
    end subroutine lanczos_diagonalize

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
    ! This routine uses only two n-vector
    subroutine lanczos_iteration(hx, hxpy, nloc, v_init, maxnstep, nstep, a, b)

        use numeric_utils, only: mpi_dot_product, mpi_norm
        use timer
        use mpi

        interface 
            ! Y = H*X
            subroutine hx(n,x,y)
                integer(kind=8), intent(in) :: n
                double precision, intent(in) :: x(n)
                double precision, intent(out) :: y(n)
            end subroutine hx
            ! Y = Y + H*X
            subroutine hxpy(n,x,y)
                integer(kind=8), intent(in) :: n
                double precision, intent(in) :: x(n)
                double precision, intent(out) :: y(n)
            end subroutine hxpy
        end interface

        integer(kind=8), intent(in) :: &
            nloc        ! dimension of the vector local to the node
        integer, intent(in) :: &
            maxnstep         ! maximum number of iteration steps

        double precision, intent(in) :: &
            v_init(nloc)     ! starting vector for lanczos iteration

        integer, intent(out) :: &
            nstep         ! number of calculated steps <= maxnstep
                          ! if b(i) = 0 for some i, nstep = i-1
        double precision, intent(out) :: &
            a(maxnstep), & ! diagonal matrix element in lanczos basis
            b(maxnstep)    ! off-diagonal matrix element in lanczos basis

        ! temporary lanczos vectors
        double precision, allocatable :: v(:), w(:)
        double precision :: norm_v, t  
        
        integer(kind=8) :: i, j, ierr, k
        character(len=100) :: msg

        allocate(v(nloc), w(nloc))
        
        a(maxnstep) = 0.0D0
        b(maxnstep) = 0.0D0

        ! Lanczos steps
        ! ref: G. Golub, Matrix Computations, 4th ed., p.562 (2013)

        if (master) then
            write(msg,"(a,I4)") "Lanczos iteration", 1
            call timestamp(msg)
        endif

        ! normalize the initial vector
        norm_v = mpi_norm( v_init, nloc)

        w = v_init/norm_v

        call hx( nloc, w(:), v(:) )

        a(1) = mpi_dot_product( w, v, nloc )
        b(1) = 0.d0

        if (maxnstep==1) then
            nstep = 1
            deallocate(v,w)
            return
        endif

        call daxpy( nloc, -a(1), w, 1, v, 1 )
        b(2) = mpi_norm( v, nloc )

        k = 2
        do while (1)
            if (master) then
                write(msg,"(a,I4)") "Lanczos iteration", k
                call timestamp(msg)
            endif
            do i=1,nloc
                t = w(i)
                w(i) = v(i)/b(k)
                v(i) = -b(k)*t
            enddo

            call hxpy( nloc, w(:), v(:) )
            a(k) = mpi_dot_product( w, v, nloc ) 

            if (k>=maxnstep.or.b(k)<1.d-16) then
                nstep = k
                exit
            endif

            call daxpy( nloc, -a(k), w, 1, v, 1 )
            b(k+1) = mpi_norm( v, nloc )

            k = k+1
        enddo

        deallocate(v,w)
    end subroutine lanczos_iteration

    subroutine lanczos_ground_state(hx, hxpy, nloc, v_init, nstep, &
                                       a, b, coeff, E0, gs, residual)

        use numeric_utils, only: mpi_dot_product, mpi_norm

        interface 
            ! Y = H*X
            subroutine hx(n,x,y)
                integer(kind=8), intent(in) :: n
                double precision, intent(in) :: x(n)
                double precision, intent(out) :: y(n)
            end subroutine hx
            ! Y = Y + H*X
            subroutine hxpy(n,x,y)
                integer(kind=8), intent(in) :: n
                double precision, intent(in) :: x(n)
                double precision, intent(out) :: y(n)
            end subroutine hxpy
        end interface

        integer(kind=8), intent(in) :: &
            nloc, &       ! dimension of the vector local to the node
            nstep         

        double precision, intent(in) :: &
            v_init(nloc), &   ! starting vector for lanczos iteration
            E0,           &  
            a(nstep),     &   ! diagonal matrix element in lanczos basis
            b(nstep),     &   ! off-diagonal matrix element in lanczos basis
            coeff(nstep)  

        double precision, intent(out) :: &
            gs(nloc), residual

        ! temporary lanczos vectors
        double precision, allocatable :: v(:), w(:)
        double precision :: norm_v, t
        
        integer(kind=8) :: i, j, ierr, k

        gs = 0.d0

        allocate(v(nloc), w(nloc))
        ! normalize the initial vector
        norm_v = mpi_norm( v_init, nloc)

        w = v_init/norm_v ! v1

        call daxpy( nloc, coeff(1), w, 1, gs, 1) 

        call hx( nloc, w(:), v(:) )

        if (nstep==1) then
            deallocate(v,w)
            return
        endif

        call daxpy( nloc, -a(1), w, 1, v, 1 )

        k = 2
        do while (1)
            do i=1,nloc
                t = w(i)
                w(i) = v(i)/b(k)
                v(i) = -b(k)*t
            enddo
            ! w = |v_k>
            call daxpy( nloc, coeff(k), w, 1, gs, 1 ) 

            call hxpy( nloc, w(:), v(:) )

            if (k>=nstep.or.b(k)<1.d-16) then
                exit
            endif
            
            call daxpy( nloc, -a(k), w, 1, v, 1 )
            k = k+1
        enddo

        ! test eigenvector residual | H|gs>-E0|gs> |^2
        call hx( nloc, gs, v )
        call daxpy( nloc, -E0, gs, 1, v, 1)

        residual = mpi_dot_product( v, v, nloc )        

        deallocate(v,w)
    end subroutine lanczos_ground_state
end module lanczos
