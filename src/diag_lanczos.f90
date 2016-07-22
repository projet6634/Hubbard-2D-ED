subroutine diag_lanczos(isector, basis, nnz, H, row_idx, col_ptr, &
                        maxnev, nev, eigval, gs)
    use utils, only: die
    use ed_params, only: maxnstep
    use ed_basis, only: basis_t
    use lanczos
    use mpi

    type(basis_t), intent(in) :: basis
    integer, intent(in) :: isector, nnz, row_idx(nnz), &
                           col_ptr(basis%nloc+1), maxnev
    double precision, intent(in) :: H(nnz) 
    double precision, intent(out) :: eigval(maxnev), gs(basis%nloc)

    integer :: nstep
    double precision, allocatable :: a(:), b(:), v_init(:), x_all(:), &
                                     lanczos_v(:,:), ev(:)

    integer :: i
    double precision :: t1, t2, r, residual

    allocate(a(maxnstep),b(maxnstep))
    allocate(v_init(basis%nloc))
    allocate(x_all(basis%ntot))

    ! initial random v
    do i=1,basis%nloc
        call random_number(r)
        v_init(i) = r
    enddo

    t1 = mpi_wtime(mpierr)
    call lanczos_iteration(hx, hxpy, basis%nloc, v_init, maxnstep, nstep, a, b)
    t2 = mpi_wtime(mpierr)
    if (master) then
        print *, "nstep = ", nstep
        print *, "lanczos iteration time = ", (t2-t1), " sec."
    endif

    nev = maxnev

    allocate(ev(nstep),lanczos_v(nstep,nstep))
    if (master) then
        t1 = mpi_wtime(mpierr)
        call lanczos_diagonalize(nstep, a, b, ev, lanczos_v)
        t2 = mpi_wtime(mpierr)
        
        print *, "lanczos diagonalization time = ", (t2-t1), " sec."

        do i=1,nev
            eigval(i) = ev(i)
        enddo
    endif

    call mpi_bcast(lanczos_v, nstep*nstep, mpi_double_precision, 0, comm, mpierr)
    call mpi_bcast(eigval, nev, mpi_double_precision, 0, comm, mpierr)

    t1 = mpi_wtime(mpierr)
    call lanczos_ground_state(hx, hxpy, basis%nloc, v_init, nstep, &
                              a, b, lanczos_v(:,1), eigval(1), gs, residual)
    t2 = mpi_wtime(mpierr)

    if (master) then
        print *, "lanczos ground state time = ", (t2-t1), " sec."
        print *, "residual = ", residual
    endif

    deallocate(v_init,x_all,a,b)
contains

    subroutine hx(n, x, y)
        integer, intent(in) :: n
        double precision, intent(in) :: x(n)
        double precision, intent(out) :: y(n)

        integer :: icol, i, nnzcol, irow, inz
        double precision :: matel, colsum

        call mpi_allgatherv(x,basis%nloc,mpi_double_precision,x_all,&
            basis%nlocals,basis%offsets,mpi_double_precision,comm,mpierr)

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
    end subroutine hx

    subroutine hxpy(n, x, y)
        integer, intent(in) :: n
        double precision, intent(in) :: x(n)
        double precision, intent(out) :: y(n)

        integer :: icol, i, nnzcol, irow, inz
        double precision :: matel, colsum

        call mpi_allgatherv(x,basis%nloc,mpi_double_precision,x_all,&
            basis%nlocals,basis%offsets,mpi_double_precision,comm,mpierr)

        do icol=1,basis%nloc 
            colsum = 0.0d0
            nnzcol = col_ptr(icol+1)-col_ptr(icol)

            do i=1,nnzcol
                inz = col_ptr(icol)+i-1
                irow = row_idx(inz)
                matel = H(inz)
                colsum = colsum + x_all(irow)*matel
            enddo
            y(icol) = y(icol) + colsum 
        enddo

    end subroutine hxpy
end subroutine diag_lanczos
