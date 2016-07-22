subroutine diag_feast(isector, basis, nnz, H, row_idx, col_ptr, &
                        maxnev, nev, eigval, gs)
    use utils, only: die
    use ed_params, only: maxnstep
    use ed_basis, only: basis_t
    use lanczos
    use mpi

    type(basis_t), intent(in) :: basis
    integer, intent(in) :: isector, nnz, row_idx(nnz), col_ptr(basis%nloc), maxnev
    double precision, intent(in) :: H(nnz) 
    double precision, intent(out) :: eigval(maxnev), gs(basis%nloc)

end subroutine diag_feast
