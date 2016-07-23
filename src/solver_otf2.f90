module solver_otf2
    use mpi
    use ed_params, only: maxnstep, nsite, U, mu
    use ed_basis, only: generate_basis, basis_t, ed_basis_get_g, &
                        ed_basis_idx, kind_basis
    use utils, only: die
    use ed_grid
    use ed_lattice, only: nnsite, nnn
    use ed_hamiltonian, only: permsgn
    use numeric_utils
    use timer


    public :: &
        solve_otf2

    type(basis_t) :: basis

    private
contains

    subroutine solve_otf2( E0, nloc, gs, G )
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

        t1 = mpi_wtime(mpierr)
        call diag_lanczos( basis, E0, gs )
        t2 = mpi_wtime(mpierr)
        if(taskid==0) write(6,'(a,3x,f10.5,3x,a)') &
                             "walltime = ",(t2-t1)/60.D0," min."

        ! call green_diag( E0, nloc, gs, G )
        G = 0.d0
    end subroutine solve_otf2

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
        call lanczos_iteration(basis%nloc, v_init, nstep, a, b)
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
        ! t1 = mpi_wtime(mpierr)
        ! call lanczos_ground_state(hx, hxpy, basis%nloc, v_init, nstep, &
        !                           a, b, coeff, E0, gs, residual)
        ! t2 = mpi_wtime(mpierr)

        ! if (master) then
        !     print *, "lanczos ground state time = ", (t2-t1), " sec."
        !     print *, "residual = ", residual
        ! endif

        deallocate(v_init,a,b)
    end subroutine diag_lanczos

    subroutine hx(n, x, y, work, py)
        integer, intent(in) :: n
        double precision, intent(in) :: x(n)
        logical, intent(in) :: py ! plus y or not
        double precision, intent(out) :: y(n), work(n+1)

        integer :: isite, jsite, ip
        integer :: ket, bra, icol, irow
        integer :: nrows
        logical :: nj(2)

        double precision :: locsum

        do ip=0,nprocs-1
            nrows = basis%nlocals(ip)

            do irow=1,nrows
                ket = ed_basis_get_g(basis, basis%offsets(ip)+irow)

                if (py.and.ip==taskid) then
                    locsum = y(irow)
                else
                    locsum = 0.d0
                endif

                do jsite=1,nsite
                    ! diagonal
                    nj(1) = BTEST(ket,jsite-1)
                    nj(2) = BTEST(ket,nsite+jsite-1)

                    if (ip==taskid) then
                        if (nj(1).and.nj(2)) then
                            locsum = locsum + (U-2*mu)*x(irow)
                        else if (nj(1).or.nj(2)) then
                            locsum = locsum - mu*x(irow)
                        endif
                    endif

                    ! collect all the off-diagonal terms
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
                                    bra = IBSET(IBCLR(ket,is+jsite-1),is+isite-1)
                                    icol = ed_basis_idx(basis, bra)
                                    if (basis%gidx1 <= icol .and. &
                                        icol <= basis%gidx2) then
                                        icol = icol-basis%gidx1+1
                                        ! t == -1 
                                        locsum = locsum - nnsite(2,i,jsite) &
                                            *permsgn(ket,is,isite,jsite) &
                                            *x(icol)
                                    endif
                                endif
                            enddo
                        endif
                    enddo
                enddo

                work(irow) = locsum
            enddo

            call mpi_reduce(work, y, nrows, mpi_double_precision, mpi_sum, &
                            ip, comm, mpierr)
        enddo

    end subroutine hx

    subroutine lanczos_iteration(nloc, v_init, nstep, a, b)
        integer, intent(in) :: nloc
        double precision, intent(in) :: &
            v_init(nloc)     ! starting vector for lanczos iteration

        integer, intent(out) :: &
            nstep         ! number of calculated steps <= maxnstep
                          ! if b(i) = 0 for some i, nstep = i-1
        double precision, intent(out) :: &
            a(maxnstep), & ! diagonal matrix element in lanczos basis
            b(maxnstep)    ! off-diagonal matrix element in lanczos basis

        ! temporary lanczos vectors
        double precision, allocatable :: v(:), w(:), work(:)
        double precision :: norm_v, t  
        
        integer :: i, j, ierr, k
        character(len=100) :: msg

        allocate(v(nloc), w(nloc), work(nloc+1))
        
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

        call hx( nloc, w, v, work, .false. )

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

            call hx( nloc, w, v, work, .true. )
            a(k) = mpi_dot_product( w, v, nloc ) 

            if (k>=maxnstep.or.b(k)<1.d-16) then
                nstep = k
                exit
            endif

            call daxpy( nloc, -a(k), w, 1, v, 1 )
            b(k+1) = mpi_norm( v, nloc )

            k = k+1
        enddo

        deallocate(v,w,work)
    end subroutine lanczos_iteration

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

end module solver_otf2
