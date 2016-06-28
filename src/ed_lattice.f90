module ed_lattice

    use ed_params
    implicit none


    public :: ed_lattice_init

    integer, allocatable, public :: &
        nnn(:),       & ! nnn(nsite) number of nearest neighbors for each site
        nnsite(:,:,:)   ! nnsite(1,nnn(isite),nsite) nn site indices
                        ! nnsite(2,nnn(isite),nsite) weight

    private 
contains

    ! square lattice neighbor info
    ! @TODO only square lattice with no tilting is considered (nnn=4)
    subroutine ed_lattice_init
        integer isite, x, y, xm1, ym1, xp1, yp1, i, j, idx
        logical duplicate
        integer, allocatable :: nnsiteall(:,:)

        allocate(nnn(Lx*Ly))
        allocate(nnsiteall(4,Lx*Ly))
        allocate(nnsite(2,4,Lx*Ly))

        do y=1,Ly
            do x=1,Lx
                ! site index used for basis
                isite = get_site_idx(x,y)

                xm1=x-1
                xp1=x+1
                ym1=y-1
                yp1=y+1

                if (x==1) then
                    xm1 = Lx
                else if (x==Lx) then
                    xp1 = 1
                endif
                
                if (y==1) then
                    ym1 = Ly
                else if (y==Ly) then
                    yp1 = 1
                endif

                if (Lx == 1) then
                    xm1=1
                    xp1=1
                endif
                if (Ly == 1) then
                    ym1=1
                    yp1=1
                endif

                nnsiteall(1,isite) = get_site_idx(xm1,y)
                nnsiteall(2,isite) = get_site_idx(xp1,y)
                nnsiteall(3,isite) = get_site_idx(x,ym1)
                nnsiteall(4,isite) = get_site_idx(x,yp1)
            enddo
        enddo

        ! pick unique sites with weights
        do isite=1,nsite
            ! first nn
            idx = 1
            nnsite(1,idx,isite) = nnsiteall(idx,isite)
            nnsite(2,idx,isite) = 1
            nnn(isite) = 1
            do i=2,4
                duplicate = .false.
                jloop: do j=1,i-1
                    if (nnsiteall(j,isite)==nnsiteall(i,isite)) then
                        duplicate = .true.
                        exit jloop
                    endif
                enddo jloop

                if (duplicate) then
                    nnsite(2,idx,isite) = nnsite(2,idx,isite)+1
                else
                    idx = idx + 1
                    nnn(isite) = nnn(isite)+1
                    nnsite(1,idx,isite) = nnsiteall(i,isite)
                    nnsite(2,idx,isite) = 1
                endif
            enddo
            ! DBG OK
            ! print *, "isite = ", isite
            ! print *, "nnn=",nnn(isite)
            ! do i=1,nnn(isite)
            !     print *, "nn,w",nnsite(1,i,isite),nnsite(2,i,isite)
            ! enddo
            ! print *, ""
        enddo
    end subroutine ed_lattice_init

    integer function get_site_idx(x,y) result(idx)
        integer :: x,y
        idx = (y-1)*Lx+x
    end function get_site_idx
end module ed_lattice
