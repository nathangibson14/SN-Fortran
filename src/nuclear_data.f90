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

subroutine SHEM_pin_cell(sigma)
    implicit none
    
    type(xsdata), dimension(2), intent(out)     :: sigma
    type(xsdata), dimension(4)                  :: xs
    
    real(kind=8)        :: N_H, N_O, N_235, N_238
    
    xs(1) = get_nuclear_data(361,'SHEM361/h1')
    xs(2) = get_nuclear_data(361,'SHEM361/o16')
    xs(3) = get_nuclear_data(361,'SHEM361/u235')
    xs(4) = get_nuclear_data(361,'SHEM361/u238')
    
    allocate(sigma(1)%tot(361),sigma(1)%chi(361),sigma(1)%nufis(361), &
        & sigma(1)%scat(361,361))
    allocate(sigma(2)%tot(361),sigma(2)%chi(361),sigma(2)%nufis(361), &
        & sigma(2)%scat(361,361))
    sigma(1)%groups = 361
    sigma(2)%groups = 361
    
    ! material 1
    N_H = 0.0
    N_O = 0.0446
    N_235 = 6.691e-4 
    N_238 = 0.0216
    sigma(1)%tot   = N_H*xs(1)%tot + N_O*xs(2)%tot + N_235*xs(3)%tot + N_238*xs(4)%tot
    sigma(1)%nufis = N_H*xs(1)%nufis + N_O*xs(2)%nufis + N_235*xs(3)%nufis + N_238*xs(4)%nufis
    sigma(1)%chi   = xs(3)%chi
    sigma(1)%scat  = N_H*xs(1)%scat + N_O*xs(2)%scat + N_235*xs(3)%scat + N_238*xs(4)%scat

    ! material 2
    N_H = 0.0669
    N_O = 0.0335
    N_235 = 0.0
    N_238 = 0.0
    sigma(2)%tot   = N_H*xs(1)%tot + N_O*xs(2)%tot + N_235*xs(3)%tot + N_238*xs(4)%tot
    sigma(2)%nufis = N_H*xs(1)%nufis + N_O*xs(2)%nufis + N_235*xs(3)%nufis + N_238*xs(4)%nufis
    sigma(2)%chi   = xs(3)%chi
    sigma(2)%scat  = N_H*xs(1)%scat + N_O*xs(2)%scat + N_235*xs(3)%scat + N_238*xs(4)%scat
    

end subroutine SHEM_pin_cell

end module nuclear_data
