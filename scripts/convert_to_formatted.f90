program main
    implicit none

    integer :: isector, nnz, taskid
    double precision, allocatable :: H(:)
    integer, allocatable :: row_idx(:), col_ptr(:)

    character(len=100) :: fn
    integer :: iuh, iuv, nloc, i, j

    read(*,*) isector, taskid

    write(fn,"(a,I2.2,a,I4.4,a)") "h_val_",isector,"-",taskid,".dat"
    iuv = 101+2*taskid
    open(unit=iuv,file=fn,form="unformatted",status="old",action="read")

    write(fn,"(a,I2.2,a,I4.4,a)") "h_header_",isector,"-",taskid,".dat"
    iuh = 101+2*taskid+1
    open(unit=iuh,file=fn,form="unformatted",status="old",action="read")

    read(iuh) nloc

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

    ! write formatted
    write(fn,"(a,I2.2,a,I4.4,a)") "h_val_",isector,"-",taskid,".txt"
    iuv = 101+2*taskid
    open(unit=iuv,file=fn,form="formatted",status="replace")

    write(fn,"(a,I2.2,a,I4.4,a)") "h_header_",isector,"-",taskid,".txt"
    iuh = 101+2*taskid+1
    open(unit=iuh,file=fn,form="formatted",status="replace")

    write(iuh,*) nloc
    do i=1,nloc
        write(iuh,*) col_ptr(i)
    enddo
    write(iuh,*) nnz

    do i=1,nnz
        write(iuv,*) row_idx(i), H(i)
    enddo
end program main
