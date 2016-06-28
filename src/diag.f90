subroutine diag(isector,basis,nnz,H,row_idx,col_ptr,nev,eigval,gs)
    use utils, only: die
    use ed_basis, only: basis_t
    use ed_hamiltonian, only: multiply_H
    use mpi
    ! include 'debug.h'
    include 'stat.h'
    type(basis_t), intent(in) :: basis
    integer, intent(in) :: isector, nnz, row_idx(nnz), col_ptr(basis%nloc), nev
    double precision, intent(in) :: H(nnz) 
    double precision, intent(out) :: eigval(nev), gs(basis%nloc)
    double precision, allocatable :: workl(:), workd(:), d(:,:), resid(:), &
                                     ax(:), v(:,:)
    logical, allocatable :: select(:)
    character, parameter :: bmat = 'I'
    character(len=2), parameter :: which = 'SA'
    integer :: iparam(11), ipntr(11), lworkl, info, ido, nconv, i, j, &
        maxitr, mode, ishfts, ncv, ldv
    double precision :: sigma, tol, pdnorm2
    external :: pdnorm2, daxpy

    ! ndigit = -3
    ! logfil = 6
    ! msaupd = 3
    ncv = 2*nev+nev/2
    ldv=basis%nloc
    lworkl = ncv*(ncv+8)
    allocate(workl(ncv*(ncv+8)))
    allocate(workd(3*basis%nloc),d(ncv,2))
    allocate(resid(basis%nloc),ax(basis%nloc),v(basis%nloc,ncv))
    allocate(select(ncv))
   
    tol = 0.0

    ishfts = 1
    maxitr = 500
    mode   = 1

    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode

    info = 0
    ido = 0

    do
        call pdsaupd( comm, ido, bmat, basis%nloc, which, nev, tol, resid, &
            ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
        if (ido .eq. -1 .or. ido .eq. 1) then
            call multiply_H(basis,nnz,H,row_idx,col_ptr,&
                            workd(ipntr(1)),workd(ipntr(2)))
        else
            exit
        endif
    enddo

    if ( info .lt. 0 ) then
        if ( node.eq. 0 ) then
            print *, ' Error with pdsaupd, info = ', info
            print *, iparam(5)
            call die("diag","")
        endif
    else
        call pdseupd(comm,.true.,'All',select,d,v,ldv,sigma, &
            bmat, basis%nloc, which, nev, tol, resid, ncv, v, ldv, &
            iparam, ipntr, workd, workl, lworkl, ierr )
        if ( ierr .ne. 0) then
            if ( node .eq. 0 ) then
                print *, ' Error with pdseupd, info = ', ierr
                call die("diag","")
            endif
        else
            nconv =  iparam(5)
            do j=1, nconv
                call multiply_H(basis,nnz,H,row_idx,col_ptr,v(:,j),ax)
                call daxpy(basis%nloc, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = pdnorm2( comm, basis%nloc, ax, 1 )
            enddo
            call pdmout(comm, 6, nconv, 2, d, ncv, -6, &
                'Ritz values and direct residuals')
        end if
    endif

    eigval = d(:,1)
    gs = v(:,1)
end subroutine diag
