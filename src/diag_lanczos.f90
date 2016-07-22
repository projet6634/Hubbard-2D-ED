subroutine diag_lanczos(isector, basis, nnz, H, row_idx, col_ptr, &
                        maxnev, nev, eigval, gs)
    use utils, only: die
    use ed_params, only: maxnstep
    use ed_basis, only: basis_t
    use lanczos
    use mpi

    type(basis_t), intent(in) :: basis
    integer, intent(in) :: isector, nnz, row_idx(nnz), col_ptr(basis%nloc), maxnev
    double precision, intent(in) :: H(nnz) 
    double precision, intent(out) :: eigval(maxnev), gs(basis%nloc)

    integer :: nstep
    double precision, allocatable :: a(:), b(:), v_init(:), x_all(:), v(:,:), ev(:)

    integer :: i
    double precision :: t1, t2, r


    allocate(a(maxnstep),b(maxnstep))
    allocate(v_init(basis%nloc))
    allocate(x_all(basis%ntot))

    ! initial random v
    do i=1,basis%nloc
        call random_number(r)
        v_init(i) = r
    enddo

    t1 = mpi_wtime(mpierr)
    call lanczos_iteration(matmult, basis%nloc, v_init, maxnstep, nstep, a, b)
    t2 = mpi_wtime(mpierr)
    if (master) then
        print *, "lanczos iteration time = ", (t2-t1), " sec."
    endif

    if (master) then
        allocate(ev(nstep),v(nstep,nstep))
        t1 = mpi_wtime(mpierr)
        call lanczos_diagonalize(nstep, a, b, ev, v)
        t2 = mpi_wtime(mpierr)
        if (master) then
            print *, "lanczos diagonalization time = ", (t2-t1), " sec."
        endif

        nev = maxnev
        do i=1,maxnev
            eigval(i) = ev(i)
        enddo
    endif

    ! @TODO ground state vector
    gs = 0.d0
    call mpi_bcast(eigval, maxnev, mpi_double_precision, &
                    0, comm, mpierr)

    deallocate(v_init,x_all,a,b)
contains

    subroutine matmult(n, x, y)
        integer, intent(in) :: n
        double precision, intent(in) :: x(n)
        double precision, intent(out) :: y(n)

        integer :: icol, i, nnzcol, irow, inz
        double precision :: matel, colsum

        call mpi_allgatherv(x,basis%nloc,mpi_double_precision,x_all,&
            basis%nlocals,basis%offsets,mpi_double_precision,comm,mpierr)

        y = 0.0d0

        do icol=1,basis%nloc 
            colsum = 0.0d0
            nnzcol = col_ptr(icol+1)-col_ptr(icol)

            do i=1,nnzcol
                inz = col_ptr(icol)+i-1
                irow = row_idx(inz)
                matel = H(inz)
                colsum = colsum + x_all(irow)*matel
            enddo
            y(icol) = colsum
        enddo
    end subroutine matmult
end subroutine diag_lanczos
