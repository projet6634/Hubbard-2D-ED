module lanczos

    use numeric_utils, only: mpi_dot_product, mpi_norm

    implicit none

contains

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
    subroutine lanczos_iteration(matmult, nloc, vec, nstep, a, b)
        integer, intent(in) :: &
            nloc, &       ! dimension of the vector local to the node
            nstep         ! maximum number of iteration steps

        double precision, intent(in) :: &
            vec(nloc)     ! starting vector for lanczos iteration
        double precision, intent(out) :: &
            a(nstep), & ! diagonal matrix element in lanczos basis
            b(nstep)    ! off-diagonal matrix element in lanczos basis

        interface 
            ! external subroutine for matrix multiplication
            ! Y = M*X
            subroutine matmult(n,X,Y)
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
        norm_v = mpi_norm( vec, nloc)
        v(:,2) = vec/norm_v
        v(:,1) = 0.0D0

        a(:) = 0.0D0
        b(:) = 0.0D0

        ! v(:,1) = v_(j-1)
        ! v(:,2) = v_j
        ! w(:)   = w_j
        lanczos_loop: do j=1,nstep-1
            ! w_j = H*v_j
            call matmult( nloc, v(:,2), w(:) )

            ! a_j = dot(w_j,v_j)
            a(j) = mpi_dot_product(w(:), v(:,2), nloc)

            ! w_j = w_j - a_j * v_j - b_j * v_(j-1)
            w(:) = w(:) - a(j)*v(:,2) - b(j)*v(:,1)

            ! b_(j+1) = norm(w_j)
            b(j+1) = mpi_norm(w(:), nloc)

            if (b(j+1).lt.1e-8) then
                exit lanczos_loop
            endif

            ! v_(j-1) = v_j
            v(:,1) = v(:,2)

            ! v_(j+1) = w_j/b(j+1)
            v(:,2) = w(:)/b(j+1)
        enddo lanczos_loop

        ! handles the last step
        call matmult( nloc, v(:,2), w(:) )
        a(nstep) = mpi_dot_product(w(:),v(:,2),nloc)

    end subroutine lanczos_iteration
end module lanczos
