module solver_otf
    use mpi
    use ed_params, only: maxnstep, nsite, U, mu
    use ed_basis, only: generate_basis, basis_t, ed_basis_get, &
                        ed_basis_idx, kind_basis
    use utils, only: die
    use lanczos
    use ed_grid
    use ed_lattice, only: nnsite, nnn
    use ed_hamiltonian, only: permsgn

    public :: &
        solve_otf

    type(basis_t) :: basis

    double precision, allocatable :: x_all(:)

    private
contains

    subroutine solve_otf( E0, nloc, gs, G )
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

        allocate(x_all(basis%ntot))

        t1 = mpi_wtime(mpierr)
        call diag_lanczos( basis, E0, gs )
        t2 = mpi_wtime(mpierr)
        if(taskid==0) write(6,'(a,3x,f10.5,3x,a)') &
                             "walltime = ",(t2-t1)/60.D0," min."

        ! call green_diag( E0, nloc, gs, G )
        G = 0.d0

        deallocate(x_all)
    end subroutine solve_otf

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

        gs = 0.d0
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

    subroutine hx(n, x, y)
        integer, intent(in) :: n
        double precision, intent(in) :: x(n)
        double precision, intent(out) :: y(n)

        integer :: isite, jsite
        integer :: ket, bra, icol, irow, a
        double precision :: rowsum, val
        logical :: nj(2), ni

        call mpi_allgatherv(x, basis%nloc, mpi_double_precision, x_all,&
            basis%nlocals, basis%offsets, mpi_double_precision, comm, mpierr)

        do irow=1,basis%nloc 
            rowsum = 0.0d0
            ket = ed_basis_get(basis, irow)

            do jsite=1,nsite
                ! diagonal
                nj(1) = BTEST(ket,jsite-1)
                nj(2) = BTEST(ket,nsite+jsite-1)
                if (nj(1).and.nj(2)) then
                    rowsum = rowsum + (U-2*mu)*x(irow)
                else if (nj(1).or.nj(2)) then
                    rowsum = rowsum - mu*x(irow)
                endif

                ! off-diagonal
                do ispin=1,2
                    if (nj(ispin)) then
                        ! site index is due to the spin
                        is = (ispin-1)*nsite

                        ! for each neighbors
                        do i=1,nnn(jsite)
                            ! nearest-neighbor site index
                            isite = nnsite(1,i,jsite)
                            if (.not.BTEST(ket,is+isite-1)) then
                                ! c^+_{isite,ispin} c_{jsite,ispin} 
                                ! t == -1 
                                val = -nnsite(2,i,jsite)&
                                      *permsgn(ket,is,isite,jsite)
                                bra = IBCLR(ket,is+jsite-1)
                                bra = IBSET(bra,is+isite-1)
                                icol = ed_basis_idx(basis, bra)
                                rowsum = rowsum + val*x_all(icol)
                            endif
                        enddo
                    endif
                enddo
            enddo

            y(irow) = rowsum
        enddo
    end subroutine hx

    subroutine hxpy(n, x, y)
        integer, intent(in) :: n
        double precision, intent(in) :: x(n)
        double precision, intent(out) :: y(n)

        integer :: isite, jsite
        integer :: ket, bra, icol, irow
        double precision :: rowsum, val
        logical :: nj(2), ni

        call mpi_allgatherv(x, basis%nloc, mpi_double_precision, x_all,&
            basis%nlocals, basis%offsets, mpi_double_precision, comm, mpierr)

        do irow=1,basis%nloc 
            rowsum = 0.0d0
            ket = ed_basis_get(basis, irow)

            do jsite=1,nsite
                ! diagonal
                nj(1) = BTEST(ket,jsite-1)
                nj(2) = BTEST(ket,nsite+jsite-1)
                if (nj(1).and.nj(2)) then
                    rowsum = rowsum + (U-2*mu)*x(irow)
                else if (nj(1).or.nj(2)) then
                    rowsum = rowsum - mu*x(irow)
                endif

                ! off-diagonal
                do ispin=1,2
                    if (nj(ispin)) then
                        is = (ispin-1)*nsite

                        ! for each neighbors
                        do i=1,nnn(jsite)
                            ! nearest-neighbor site index
                            isite = nnsite(1,i,jsite)
                            if (.not.BTEST(ket,is+isite-1)) then
                                ! c^+_{isite,ispin} c_{jsite,ispin} 
                                ! t == -1 
                                val = -nnsite(2,i,jsite)&
                                      *permsgn(ket,is,isite,jsite)
                                bra = IBCLR(ket,is+jsite-1)
                                bra = IBSET(bra,is+isite-1)
                                icol = ed_basis_idx(basis, bra)
                                rowsum = rowsum + val*x_all(icol)
                            endif
                        enddo
                    endif
                enddo
            enddo

            y(irow) = y(irow) + rowsum
        enddo
    end subroutine hxpy
end module solver_otf
