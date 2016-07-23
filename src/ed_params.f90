module ed_params

    use fdf
    use mpi
    use utils, only: die
    use numeric_utils, only: icom

    integer :: &
        Lx, Ly,         &    ! lattice size x,y
        Nsite,          &    ! number of sites
        maxnstep,       &    ! maximum lanczos iteration steps
        diag_method,    &    ! 1 = csr 2 = otf
        nw

    double precision :: &
        U,       &       ! Coulomb repulsion
        mu,      &       ! chemical potential
        wmin,    &       ! wmin
        wmax,    &       ! wmax
        eta              ! broadening      

    logical :: &
        read_hamiltonian

    character(len=200) :: &
        outdir

contains

    subroutine read_input
        type(block_fdf)            :: bfdf
        type(parsed_line), pointer :: pline
        integer :: i
        logical :: res

        Lx = fdf_get("Lx",2)
        Ly = fdf_get("Ly",2)
        Nsite = Lx*Ly ! even number of sites only

        if (mod(Nsite,2)==1) then
            call die("main","Only even number of sites are allowed.")
        endif

        U = fdf_get("U",1.0d0)
        mu = fdf_get("mu",0.5d0)

        read_hamiltonian = fdf_get("UseSavedHamiltonian",.false.)

        maxnstep = fdf_get("maxnstep",30)
        outdir = fdf_get("outdir","./")
        res = makedirqq(outdir)

        nw = fdf_get("nw", 1000)
        wmin = fdf_get("wmin", -2.d0)
        wmax = fdf_get("wmax", 2.d0)
        eta = fdf_get("broadening", 0.01d0)

        diag_method = fdf_get("Diag", 1)
    end subroutine read_input
end module ed_params
