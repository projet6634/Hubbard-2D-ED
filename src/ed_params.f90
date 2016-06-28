module ed_params

    use fdf
    use mpi
    use utils, only: die
    use numeric_utils, only: icom

    integer :: &
        Lx, Ly,     &    ! lattice size x,y
        Nsite,      &    ! number of sites
        nsector,    &    ! number of (Q,Sz) sectors
        nev

    integer, allocatable :: &
        sectors(:,:) ! (nup,ndown) arrays

    double precision :: &
        U, &         ! Coulomb repulsion
        mu           ! chemical potential

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

        if (mod(Nsite,1)) then
            call die("main","Only even number of sites are allowed.")
        endif

        nsector = fdf_get("Nsector",1)
        allocate(sectors(nsector,3))
        if (fdf_block('Sectors', bfdf)) then
            i = 1
            do while((fdf_bline(bfdf, pline)) .and. (i .le. nsector))
                sectors(i,1) = fdf_bintegers(pline,1) ! nup
                sectors(i,2) = fdf_bintegers(pline,2) ! ndown
                ! dimension of the sector
                sectors(i,3) = icom(sectors(i,1)+sectors(i,2),sectors(i,1))* &
                               icom(sectors(i,1)+sectors(i,2),sectors(i,2)) 
                i = i + 1
            enddo
        else
            call die("main", "please specify the sectors")
        endif

        U = fdf_get("U",1.0d0)
        mu = fdf_get("mu",0.5d0)

        read_hamiltonian = fdf_get("UseSavedHamiltonian",.false.)

        nev = fdf_get("nev",3)
        nstep = fdf_get("continuedfractionsteps",30)
        outdir = fdf_get("outdir","./")
        res = makedirqq(outdir)
    end subroutine read_input
end module ed_params
