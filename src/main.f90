program main
    use mpi
    use fdf
    use ed_params, only: read_input, Lx, Ly, Nsite, U, mu, diag_method, &
                         nw, wmin, wmax
    use ed_lattice, only: ed_lattice_init
    use ed_basis, only: generate_basis, basis_t, ed_basis_get
    use solver_csr, only: solve_csr
    use solver_otf, only: solve_otf
    use ed_grid

    implicit none

    ! local variables
    type(basis_t) :: basis

    double precision :: t1, t2, ti, tf

    integer :: nup, ndown, i, j, k

    integer :: nloc
    double precision :: E0, dw
    double precision, allocatable :: &
        gs(:) ! ground state eigenvector
    double complex, allocatable :: &
        G(:)    

    call mpi_setup
    ti = mpi_wtime(mpierr)

    call fdf_init('input.fdf', 'fdf.out')

    call read_input

    call print_header

    call ed_lattice_init

    ! call ed_grid_init

    allocate(G(nwloc))

    select case(diag_method)
        case(1)
            call solve_csr( E0, nloc, gs, G )
        case(2)
            call solve_otf( E0, nloc, gs, G )
        case default
            stop "invalid diag_method"
    end select

    if (master) then
        write(*,"(a)") repeat("=",80)
        write(*,*) "  E                       E/Nsite "
        write(*,*) E0, E0/Nsite
        write(*,"(a)") repeat("=",80)
    endif

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
            write(*,*)
        endif
    end subroutine

end program main
