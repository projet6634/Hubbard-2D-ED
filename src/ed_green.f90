module ed_green
    use mpi
    use lanczos
    use ed_basis, only: basis_t, generate_basis, ed_basis_get, &
                        ed_basis_idx, kind_basis
    use ed_params, only: nstep
    implicit none

    integer :: nnz
    integer, allocatable :: row_idx(:), col_ptr(:)
    double precision, allocatable :: H(:), eigval
contains

    subroutine green(basis,eigval,gs)
        type(basis_t), intent(in) :: basis
        double precision, intent(in) :: gs(basis%nloc), eigval

        type(basis_t) :: basis_out
        double precision :: t1, t2
        double precision, allocatable :: vec(:),a(:),b(:)

        ! site 1, spin up
        ! c^+ |gs>
        call apply_c(basis,gs,1,1,1,basis_out,vec)

        if (taskid==0) write(*,*) "Generating Hamiltonian for "
        t1 = mpi_wtime(mpierr)
        call make_hamiltonian(isector,basis)
        t2 = mpi_wtime(mpierr)
        if(taskid==0) write(6,'(a,3x,f10.5,3x,a)') &
                             "walltime = ",(t2-t1)/60.D0," min."
        

        call lanczos_iteration(mult_H,basis_out%nloc,vec,nstep,a,b)

        contains
            subroutine mult_H(n,x,y)
                integer, intent(in) :: n
                double precision, intent(in) :: x(n)
                double precision, intent(out) :: y(n)

                call multiply_H(basis_out,nnz,

            end subroutine mult_H

    end subroutine green

    ! Apply creation/destruction operator to a state vector.
    ! pm = 1 : create
    ! pm = 2 : destroy
    subroutine apply_c(basis,vec,pm,iorb,ispin,basis_out,vec_out)
        type(basis_t), intent(in) :: basis
        double precision, intent(in) :: vec(basis%nloc)
        integer, intent(in) :: pm, iorb, ispin
        type(basis_t), intent(out) :: basis_out
        double precision, allocatable, intent(out) :: vec_out(:)

        double precision, allocatable :: vec_all(:)
        integer(kind=kind_basis) :: basis_i, basis_j
        integer :: i,j,sgntot, isite

        if (ispin.eq.1) then
            if (pm.eq.1) then
                call generate_basis( basis%ne_up+1, basis%ne_down, basis_out)
            else
                call generate_basis( basis%ne_up-1, basis%ne_down, basis_out)
            endif
        else
            if (pm.eq.1) then
                call generate_basis( basis%ne_up, basis%ne_down+1, basis_out)
            else
                call generate_basis( basis%ne_up, basis%ne_down-1, basis_out)
            endif
        endif

        allocate(vec_out(basis_out%nloc))
        allocate(vec_all(basis%ntot))
        vec_out = 0.0D0

        vec_all = vec 
        call mpi_allgatherv(vec,basis%nloc,mpi_double_precision,vec_all,&
            basis%nlocals,basis%offsets,mpi_double_precision,comm,ierr)

        do i=1,basis_out%nloc
            basis_i = ed_basis_get(basis_out,i)
            isite = (ispin-1)*nsite+iorb
            sgntot = sgn(basis_i,1,isite-1)
            ! transposed. 1=creation, 2=destruction
            if (pm.eq.1) then
                if (.not.BTEST(basis_i,isite-1)) then
                    sgntot = 0
                endif
                basis_j = IBCLR(basis_i, isite-1)
            else
                if (BTEST(basis_i,isite-1)) then
                    sgntot = 0
                endif
                basis_j = IBSET(basis_i, isite-1)
            endif

            if (sgntot.eq.0) then
                cycle
            endif

            j = ed_basis_idx(basis, basis_j)
            vec_out(i) = vec_out(i) + vec_all(j)*sgntot
        enddo

        deallocate(vec_all)
    end subroutine apply_c

    integer function sgn(basis_in,i,j)
        integer(kind=kind_basis) :: basis_in
        integer :: i,j, k, sgnsum

        sgnsum = 0
        do k=i,j
            if (BTEST(basis_in,k-1)) then
                sgnsum = sgnsum + 1
            endif
        enddo

        if (mod(sgnsum,2)) then
            sgn = -1
        else
            sgn = +1
        endif
        return
    end function sgn
end module ed_green
