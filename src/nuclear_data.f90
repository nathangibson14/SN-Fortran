module nuclear_data

type xsdata
    real(kind=8), dimension(:), allocatable     :: tot, nufis, chi
    real(kind=8), dimension(:,:), allocatable   :: scat
    real(kind=8), dimension(:,:), allocatable   :: delta
    integer                                     :: groups
end type xsdata

contains

type(xsdata) function get_nuclear_data(groups,folder)
    ! read in cross section data
    ! each cross section contained in a simple text file inside folder
    implicit none
    
    integer, intent(in)             :: groups
    character(len=*), intent(in)    :: folder
    integer                         :: g, gg
    real(kind=8)                    :: temp
    
    allocate(get_nuclear_data%tot(groups), get_nuclear_data%nufis(groups), &
        & get_nuclear_data%chi(groups), get_nuclear_data%scat(groups,groups))
    
    get_nuclear_data%groups = groups
    
    ! read in total cross section
    open(unit=20,file=folder//'/tot',action='read')
    do g=1,groups
        read(20,*) get_nuclear_data%tot(g)
    enddo
    close(unit=20)

    ! read in nu-fission cross section
    open(unit=20,file=folder//'/nufis',action='read')
    do g=1,groups
        read(20,*) get_nuclear_data%nufis(g)
    enddo
    close(unit=20)
    
    ! read in vector fission spectrum
    open(unit=20,file=folder//'/chi',action='read')
    do g=1,groups
        read(20,*) get_nuclear_data%chi(g)
    enddo
    close(unit=20)
    
    ! read in scattering matrix
    ! row major, -1 indicates end of row
    open(unit=20,file=folder//'/scat',action='read')
    do g=1,groups
        do gg=1,groups
            read(20,*) temp
            if (temp==-1) then
                exit
            else
                get_nuclear_data%scat(g,gg) = temp
            endif
        enddo
    enddo
    close(unit=20)


end function get_nuclear_data

subroutine test_pin_cell(sigma, groups, folder)
    implicit none
    
    character(len=*), intent(in)                :: folder
    integer, intent(in)                         :: groups
    type(xsdata), dimension(2), intent(out)     :: sigma
    type(xsdata), dimension(4), target          :: xs
    
    real(kind=8), dimension(4)                  :: N
    real(kind=8), dimension(:,:), pointer       :: temp_scat
    integer                                     :: i
    
    write(*,*) folder//'/h-h2o...'    
    xs(1) = get_nuclear_data(groups, folder//'/h-h2o')

    write(*,*) folder//'/o16...'
    xs(2) = get_nuclear_data(groups, folder//'/o16')
    
    write(*,*) folder//'/u235...'
    xs(3) = get_nuclear_data(groups, folder//'/u235')
    
    write(*,*) folder//'/u238...'
    xs(4) = get_nuclear_data(groups, folder//'/u238')
    
    allocate(sigma(1)%tot(groups),sigma(1)%chi(groups),sigma(1)%nufis(groups), &
        & sigma(1)%scat(groups,groups))
    allocate(sigma(2)%tot(groups),sigma(2)%chi(groups),sigma(2)%nufis(groups), &
        & sigma(2)%scat(groups,groups))
    allocate(temp_scat(groups,groups))
    sigma(1)%groups = groups
    sigma(2)%groups = groups
    
    write(*,*) 'material 1...'
    ! material 1
    N(1) = 0.0          ! H-H2O
    N(2) = 0.0446       ! O-16
    N(3) = 6.691e-4     ! U-235
    N(4) = 0.0216       ! U-238
    sigma(1)%tot   = N(1)*xs(1)%tot + N(2)*xs(2)%tot + N(3)*xs(3)%tot + N(4)*xs(4)%tot
    sigma(1)%nufis = N(1)*xs(1)%nufis + N(2)*xs(2)%nufis + N(3)*xs(3)%nufis + N(4)*xs(4)%nufis
    sigma(1)%chi   = xs(3)%chi
    sigma(1)%scat = 0.0
    do i=1,4
        temp_scat => xs(i)%scat
        sigma(1)%scat = sigma(1)%scat + N(i)*temp_scat
    enddo
!~     sigma(1)%scat  = N(1)*xs(1)%scat + N(2)*xs(2)%scat + N(3)*xs(3)%scat + N(4)*xs(4)%scat
    
    
    write(*,*) 'material 2...'
    ! material 2
    N(1)  = 0.0669      ! H-H2O
    N(2)  = 0.0335      ! O-16
    N(3)  = 0.0         ! U-235
    N(4)  = 0.0         ! U-238
    sigma(2)%tot   = N(1)*xs(1)%tot + N(2)*xs(2)%tot + N(3)*xs(3)%tot + N(4)*xs(4)%tot
    sigma(2)%nufis = N(1)*xs(1)%nufis + N(2)*xs(2)%nufis + N(3)*xs(3)%nufis + N(4)*xs(4)%nufis
    sigma(2)%chi   = xs(3)%chi
    sigma(2)%scat = 0.0
    do i=1,4
        temp_scat => xs(i)%scat
        sigma(2)%scat = sigma(2)%scat + N(i)*temp_scat
    enddo
!~     sigma(2)%scat  = N(1)*xs(1)%scat + N(2)*xs(2)%scat + N(3)*xs(3)%scat + N(4)*xs(4)%scat
    
    deallocate(xs(1)%scat,xs(2)%scat,xs(3)%scat,xs(4)%scat)


end subroutine test_pin_cell

end module nuclear_data
