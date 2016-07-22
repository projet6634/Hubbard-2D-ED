module ed_grid

    integer :: nwloc

    double precision, allocatable :: w(:)
    
contains

    subroutine ed_grid_init
        use ed_params, only: nw, wmin, wmax
        double precision :: dw

        allocate(w(nw))

        ! grid setup
        dw = (wmax-wmin)/(nw-1)
        do i=1,nw
            w(i) = wmin+(i-1)*dw
        enddo

        nwloc = nw

    end subroutine ed_grid_init

end module ed_grid
