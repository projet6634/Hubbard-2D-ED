module ed_hamiltonian
    use ed_params, only: U, Lx, Ly, Nsite, mu, outdir
    use ed_basis, only: basis_t, ed_basis_get, ed_basis_idx, kind_basis
    use utils, only: die
    use ed_lattice, only: nnsite, nnn
    use mpi

    implicit none

    public :: &
        make_hamiltonian, &
        multiply_H, &
        load_hamiltonian, &
        permsgn
    
    private
contains

    subroutine make_hamiltonian(basis)
        type(basis_t), intent(in) :: basis

        character(len=100) :: fn
        integer :: nnz, nnz_in_col

        integer(kind=kind_basis) :: bra,ket
        logical :: nj(2), ni
        integer :: a,b,isite,jsite,is,ispin,i,sgn,iuv,iuh,isector

        double precision :: val, diag

        isector = 1

        write(fn,"(2a,I2.2,a,I4.4,a)") trim(outdir),"/h_val_",isector,"-",taskid,".dat"
        iuv = 101+2*taskid
        open(unit=iuv,file=fn,form="unformatted",status="replace")

        write(fn,"(2a,I2.2,a,I4.4,a)") trim(outdir),"/h_header_",isector,"-",taskid,".dat"
        iuh = 101+2*taskid+1
        open(unit=iuh,file=fn,form="unformatted",status="replace")

        write(iuh) basis%nloc

        nnz = 0

        ! locate nonzero elements for each column
        do a = 1,basis%nloc
            nnz_in_col = 0
            diag = 0.0d0

            ! |a>
            ket = ed_basis_get(basis,a)

            do jsite=1,nsite
                ! n_{j,up} |a> = nj(1) |a>
                nj(1) = BTEST(ket,jsite-1)
                ! n_{j,dn} |a> = nj(2) |a>
                nj(2) = BTEST(ket,nsite+jsite-1)

                ! diagonal term
                if (nj(1).and.nj(2)) then
                    diag = diag + U-2*mu
                else if (nj(1).or.nj(2)) then
                    diag = diag - mu
                endif

                ! off-diagonal hopping element
                ! for spin up,down
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
                                val = -nnsite(2,i,jsite)*permsgn(ket,is,isite,jsite)
                                bra = IBCLR(ket,is+jsite-1)
                                bra = IBSET(bra,is+isite-1)
                                b = ed_basis_idx(basis,bra)

                                write(iuv) b, val
                                nnz_in_col = nnz_in_col + 1
                            endif
                        enddo
                    endif
                enddo
            enddo ! jsite loop

            nnz_in_col = nnz_in_col + 1
            write(iuv) ed_basis_idx(basis,ket), diag

            ! col_ptr(icol)
            write(iuh) (1+nnz)
            nnz = nnz + nnz_in_col
        enddo ! |a> loop

        write(iuh) nnz

        close(iuh)
        close(iuv)
    end subroutine make_hamiltonian

    ! permutation sign for c^+_{isite}c_{jsite} | ... >
    integer function permsgn(ket,offset,isite,jsite) result(s)
        integer :: isite, jsite, offset
        integer :: i,j,k
        integer(kind=kind_basis) :: ket

        if (isite>jsite) then
            s = -1
        else
            s = 1
        endif
        i = min(isite,jsite)
        j = max(isite,jsite)

        do k=i,j-1 
            if (BTEST(ket,offset+k-1)) then
                s = -s
            endif
        enddo
    end function permsgn

    ! read hamiltonian in ccs format
    subroutine load_hamiltonian(basis,nnz,H,row_idx,col_ptr)
        type(basis_t), intent(in) :: basis
        integer, intent(out) :: nnz
        double precision, allocatable, intent(out) :: H(:)
        integer, allocatable, intent(out) :: row_idx(:), col_ptr(:)

        character(len=100) :: fn
        integer :: iuh, iuv, nloc, i, j, isector

        isector = 1

        write(fn,"(2a,I2.2,a,I4.4,a)") trim(outdir),"/h_val_",isector,"-",taskid,".dat"
        iuv = 101+2*taskid
        open(unit=iuv,file=fn,form="unformatted",status="old",action="read")

        write(fn,"(2a,I2.2,a,I4.4,a)") trim(outdir),"/h_header_",isector,"-",taskid,".dat"
        iuh = 101+2*taskid+1
        open(unit=iuh,file=fn,form="unformatted",status="old",action="read")

        read(iuh) nloc
        if (nloc /= basis%nloc) then
            write(*,*) "nloc read = ",nloc, ", supposed to be ",basis%nloc
            call die("load_hamiltonian", "matrix dimension mismatch")
        endif

        allocate(col_ptr(nloc+1))
        do i=1,nloc
            read(iuh) col_ptr(i)
        enddo
        read(iuh) nnz
        col_ptr(nloc+1) = 1+nnz

        allocate(H(nnz),row_idx(nnz))

        do i=1,nnz
            read(iuv) row_idx(i), H(i)
        enddo

        close(iuh)
        close(iuv)
    end subroutine load_hamiltonian

    subroutine multiply_H(basis,nnz,H,row_idx,col_ptr,x,y)
        type(basis_t), intent(in) :: basis
        integer, intent(in) :: nnz, row_idx(nnz), col_ptr(basis%nloc+1)
        double precision, intent(in) :: x(basis%nloc), H(nnz)
        double precision, intent(out) :: y(basis%nloc)

        integer :: icol, i, nnzcol, irow, inz
        double precision :: matel, colsum, x_all(basis%ntot)

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

    end subroutine multiply_H

end module ed_hamiltonian
