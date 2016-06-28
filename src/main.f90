program main
    use mpi
    use ed_params, only: read_input, Lx, Ly, Nsite, U, Nsector, sectors, &
                        read_hamiltonian, nev, mu
    use fdf
    use ed_basis, only: generate_basis, basis_t, ed_basis_get
    use utils, only: die
    use numeric_utils, only: icom
    use ed_hamiltonian, only: make_hamiltonian, load_hamiltonian
    use ed_lattice, only: ed_lattice_init
    implicit none

    ! local variables
    type(basis_t) basis

    double precision :: t1, t2, ti, tf

    integer :: isector, nup, ndown, i,j,k

    ! Hamiltonian related
    integer :: nnz
    integer, allocatable :: row_idx(:), col_ptr(:)
    double precision, allocatable :: H(:)

    double precision, allocatable :: &
        eigval(:), & ! few lowest eigenvalues
        gs(:)        ! ground state eigenvector


    call mpi_setup
    ti = mpi_wtime(mpierr)

    call fdf_init('input.fdf', 'fdf.out')

    call read_input

    call print_header

    allocate(eigval(nev))

    ! ==========================================================================
    ! setup a lattice geometry 
    ! ==========================================================================
    call ed_lattice_init

    ! ==========================================================================
    ! solve for the ground state for each sector
    ! ==========================================================================
    do isector = 1,nsector
        nup = sectors(isector,1)
        ndown = sectors(isector,2)

        if (taskid==0) write(*,*) "Generating basis states..."
        t1 = mpi_wtime(mpierr)
        call generate_basis(nup, ndown, basis)
        t2 = mpi_wtime(mpierr)
        if(taskid==0) write(6,'(a,3x,f10.5,3x,a)') &
                             "walltime = ",(t2-t1)/60.D0," min."

        if (.not.read_hamiltonian) then
            if (taskid==0) write(*,*) "Generating Hamiltonian..."
            t1 = mpi_wtime(mpierr)
            call make_hamiltonian(isector,basis)
            t2 = mpi_wtime(mpierr)
            if(taskid==0) write(6,'(a,3x,f10.5,3x,a)') &
                                 "walltime = ",(t2-t1)/60.D0," min."
        endif

        if (taskid==0) write(*,*) "Loading Hamiltonian..."
        t1 = mpi_wtime(mpierr)
        call load_hamiltonian(isector,basis,nnz,H,row_idx,col_ptr)
        t2 = mpi_wtime(mpierr)
        if(taskid==0) write(6,'(a,3x,f10.5,3x,a)') &
                             "walltime = ",(t2-t1)/60.D0," min."

        allocate(gs(basis%nloc))

        if (taskid==0) write(*,*) "Diagonalizing sector", isector
        t1 = mpi_wtime(mpierr)
        call diag(isector,basis,nnz,H,row_idx,col_ptr,nev,eigval,gs)
        t2 = mpi_wtime(mpierr)
        if(taskid==0) write(6,'(a,3x,f10.5,3x,a)') &
                             "walltime = ",(t2-t1)/60.D0," min."
        deallocate(H,row_idx,col_ptr)

        ! Zero-temperature Green's function.
        ! @TODO needs to be generalized to multiple sector cases
        ! call green(basis,eigval(1),gs)

        if (master) then
            write(*,"(a)") repeat("=",80)
            print *, "eigenvalues for sector ", isector
            write(*,"(a)") repeat("=",80)

            do i=1,nev
                write(*,*) eigval(i)/Nsite
            enddo
            write(*,"(a)") repeat("=",80)
        endif
    enddo
    

    call fdf_shutdown
    
    tf = mpi_wtime(mpierr)
    if(taskid==0) then
        write(*,"(a)") repeat("=",80)
        write(*,*) "End of run."
        write(6,'(a,3x,f10.5,3x,a)') &
                         "Total elapsed time = ",(tf-ti)/60.D0," min."
        write(*,"(a)") repeat("=",80)
    endif

    call mpi_shutdown

contains

    subroutine print_header
        if (master) then
            write(*,"(a)") repeat("=",80)
            write(*,*) "< Hubbard Model Exact Diagonalization >"
            write(*,*) "Half-filled, nearest-neighbor hoppings only."
            write(*,*) "Number of processors = ", nprocs
            write(*,"(a)") repeat("=",80)
            write(*,*) "Input Parameters"
            write(*,"(a)") repeat("=",80)
            write(*,*) "Lattice Size x              = ", Lx
            write(*,*) "Lattice Size y              = ", Ly
            write(*,*) "Number of Sites             = ", Nsite
            write(*,*) "U                           = ", U
            write(*,*) "mu                          = ", mu
            write(*,*) "nev                         = ", nev

            write(*,*) "Number of (up,down) sectors = ", nsector
            do i=1,nsector
                write(*,*) "(",sectors(i,1),",",sectors(i,2),"), dim = ", &
                    sectors(i,3)
            enddo
            write(*,*)
        endif
    end subroutine
end program main
