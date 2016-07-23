module solver_csr
    use mpi
    use ed_params, only: maxnstep, read_hamiltonian, nsite, eta
    use ed_basis, only: generate_basis, basis_t, ed_basis_get
    use ed_hamiltonian, only: make_hamiltonian, load_hamiltonian
    use utils, only: die
    use lanczos
    use ed_grid

    public :: &
        solve_csr

    ! Hamiltonian related
    integer :: nnz
    integer, allocatable :: row_idx(:), col_ptr(:)
    double precision, allocatable :: H(:), x_all(:)

    type(basis_t) :: basis

    private
contains 

    subroutine solve_csr( E0, nloc, gs, G )
        double precision, intent(out) :: &
            E0

        double complex, intent(out) :: &
            G(nwloc)      ! local Green's function
        double precision, allocatable, intent(out) :: &
            gs(:)   ! ground state eigenvector
        integer, intent(out) :: nloc

        double precision :: t1, t2

        call generate_basis(nsite/2, nsite/2, basis)

        nloc = basis%nloc
        allocate(gs(nloc))

        if (.not.read_hamiltonian) then
            if (taskid==0) write(*,*) "Generating Hamiltonian..."
            t1 = mpi_wtime(mpierr)
            call make_hamiltonian(basis)
            t2 = mpi_wtime(mpierr)
            if(taskid==0) write(6,'(1x,a,3x,f10.5,3x,a)') &
                                 "walltime = ",(t2-t1)/60.D0," min."
        endif

        if (taskid==0) write(*,*) "Loading Hamiltonian..."
        t1 = mpi_wtime(mpierr)
        call load_hamiltonian( basis, nnz, H, row_idx, col_ptr )
        t2 = mpi_wtime(mpierr)
        if(taskid==0) write(6,'(1x,a,3x,f10.5,3x,a)') &
                             "walltime = ",(t2-t1)/60.D0," min."

        allocate(x_all(basis%ntot))

        t1 = mpi_wtime(mpierr)
        call diag_lanczos( basis, E0, gs )
        t2 = mpi_wtime(mpierr)
        if(taskid==0) write(6,'(a,3x,f10.5,3x,a)') &
                             "walltime = ",(t2-t1)/60.D0," min."

        ! call green_diag( E0, nloc, gs, G )
        G = 0.d0

        deallocate(H,row_idx,col_ptr,x_all)
    end subroutine solve_csr

    subroutine diag_lanczos(basis, E0, gs)
        type(basis_t), intent(in) :: basis
        double precision, intent(out) :: E0, gs(basis%nloc)

        integer :: nstep
        double precision, allocatable :: a(:), b(:), v_init(:), &
                                         lanczos_v(:,:), ev(:), coeff(:)

        integer :: i
        double precision :: t1, t2, r, residual

        allocate(a(maxnstep),b(maxnstep))
        allocate(v_init(basis%nloc))

        call random_seed
        ! initial random v
        do i=1,basis%nloc
            call random_number(r)
            v_init(i) = r
        enddo

        t1 = mpi_wtime(mpierr)
        call lanczos_iteration(hx, hxpy, basis%nloc, v_init, &
                               maxnstep, nstep, a, b)
        t2 = mpi_wtime(mpierr)

        if (master) then
            print *, "nstep = ", nstep
            print *, "lanczos iteration time = ", (t2-t1), " sec."
        endif

        allocate(ev(nstep),lanczos_v(nstep,nstep))
        allocate(coeff(nstep))

        if (master) then
            t1 = mpi_wtime(mpierr)
            call lanczos_diagonalize(nstep, a, b, ev, lanczos_v)
            t2 = mpi_wtime(mpierr)
            
            print *, "lanczos diagonalization time = ", (t2-t1), " sec."

            E0 = ev(1)
            coeff = lanczos_v(:,1)
        endif

        call mpi_bcast(coeff, nstep, mpi_double_precision, 0, comm, mpierr)
        call mpi_bcast(E0, 1, mpi_double_precision, 0, comm, mpierr)

        t1 = mpi_wtime(mpierr)
        call lanczos_ground_state(hx, hxpy, basis%nloc, v_init, nstep, &
                                  a, b, coeff, E0, gs, residual)
        t2 = mpi_wtime(mpierr)

        if (master) then
            print *, "lanczos ground state time = ", (t2-t1), " sec."
            print *, "residual = ", residual
        endif

        deallocate(v_init,a,b)
    end subroutine diag_lanczos

    subroutine green_diag( E0, nloc, gs, G )
        use numeric_utils
                         
        integer, intent(in) :: nloc
        double precision, intent(in) :: E0, gs(nloc)

        double complex, intent(out) :: G(nwloc)

        integer :: nstep
        double precision, allocatable :: a(:), b(:), v_init(:)
        double precision :: t1, t2
        type(basis_t) :: basis_out

        integer :: isite, ispin, iw
        double complex :: gf, z
        
        G = 0.d0
        isite = 1
        ispin = 1

        allocate( a(maxnstep), b(maxnstep) )

        ! G+
        ! create 1st site, spin up
        call generate_basis( nsite/2+1, nsite/2, basis_out )
        allocate( v_init(basis_out%nloc) )
        call apply_c(basis, gs, x_all, 1, isite, ispin, basis_out, v_init)

        call lanczos_iteration(hx, hxpy, basis_out%nloc, v_init, &
                               maxnstep, nstep, a, b)
        b(1) = mpi_dot_product( v_init, v_init, basis_out%nloc )

        do iw=1,nwloc
            z = cmplx( w(iw)+E0, eta )
            gr = continued_fraction_p( z, nstep, a, b )
            G(iw) = G(iw) + gr
        enddo

        deallocate( v_init )

        deallocate( a, b )
    end subroutine green_diag

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
end module solver_csr
