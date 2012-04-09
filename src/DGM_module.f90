module DGM_module
use DCT_type
use stop_watch

type DGM_structure
    integer                             :: fine_groups
    integer                             :: coarse_groups
    integer, dimension(:), allocatable  :: cg_map
end type DGM_structure

type coarse_group
    integer                                     :: fine_groups
    integer                                     :: min_group
    real(kind=8), dimension(:,:,:), allocatable :: aflux
end type coarse_group

contains

!=======================================================================
! GET_DGM_STRUCTURE
!=======================================================================
type(DGM_structure) function get_DGM_structure(groups)
    implicit none
    
    integer, intent(in)                 :: groups
    integer                             :: coarse_groups
    integer, dimension(:), allocatable  :: cg_map
    
    integer                             :: g,i

    ! read in structure
    if (groups==3) then
        coarse_groups = 1
        allocate(cg_map(coarse_groups))
        cg_map = 3    
            
    elseif (groups==361) then
        open(unit=41,file='DGM_structure/SHEM361',action='read')
        
        read(41,*) coarse_groups
        allocate(cg_map(coarse_groups))
        do i=1,coarse_groups
            read(41,*) cg_map(i)
        enddo
        
        close(unit=41)
        
    elseif (groups==2042) then
        open(unit=41,file='DGM_structure/SHEM361',action='read')
        
        read(41,*) coarse_groups
        allocate(cg_map(coarse_groups))
        do i=1,coarse_groups
            read(41,*) cg_map(i)
        enddo
        
        close(unit=41)
        
    endif
    
    ! ensure mapping is consistent
    if (sum(cg_map) /= groups) then
        write(*,*) 'Invalid coarse group structure.'
        stop
    endif
    
    
    ! return values
    allocate(get_DGM_structure%cg_map(coarse_groups))
    get_DGM_structure%fine_groups = groups
    get_DGM_structure%coarse_groups = coarse_groups
    get_DGM_structure%cg_map = cg_map

end function get_DGM_structure

!=======================================================================
! FLUX_MOMENTS
!=======================================================================
subroutine flux_moments(psi,cg)
    implicit none
    real(kind=8), dimension(:,:,:), intent(in)      :: psi
    type(coarse_group), dimension(:), intent(inout) :: cg
    
    integer         :: xi, K, j, ai, g
    integer         :: NFM, angles, coarse_groups
    
    NFM = size(psi,2)
    angles = size(psi,3)
    coarse_groups = size(cg)
    
    do g=1,coarse_groups
!~         write(*,*) coarse_groups, NFM, angles, cg(g)%fine_groups-1
        
        cg(g)%aflux = 0.0
        do xi=1,NFM
          do ai=1,angles
            do K=0,cg(g)%fine_groups-1
              do j=0,cg(g)%fine_groups-1
                cg(g)%aflux(j,xi,ai) = cg(g)%aflux(j,xi,ai) + & 
                  & DCT(cg(g)%fine_groups)%d(j,K)*psi(K+cg(g)%min_group,xi,ai)
              enddo
            enddo
          enddo
        enddo
        
    enddo
 
end subroutine flux_moments

!=======================================================================
! BUILD_FINE_FLUX
!=======================================================================
subroutine build_fine_flux(psi,cg)
    implicit none
    
    real(kind=8), dimension(:,:,:), intent(inout)   :: psi
    type(coarse_group), dimension(:), intent(in)    :: cg
    
    integer         :: xi, K, j, ai, g
    integer         :: NFM, angles, coarse_groups
    
    NFM = size(psi,2)
    angles = size(psi,3)
    coarse_groups = size(cg)
    
    psi=0.0
    do g=1,coarse_groups
      do xi=1,NFM
        do ai=1,angles
          do K=0,cg(g)%fine_groups-1
            do j=0,cg(g)%fine_groups-1
              psi(K+cg(g)%min_group,xi,ai)=psi(K+cg(g)%min_group,xi,ai) + &
                & 2*DCT(cg(g)%fine_groups)%d(j,K)*cg(g)%aflux(j,xi,ai)/cg(g)%fine_groups
              if (j==0) then
                psi(K+cg(g)%min_group,xi,ai)=0.5_8*psi(K+cg(g)%min_group,xi,ai)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
            

end subroutine build_fine_flux

!=======================================================================
! TEST_EXPANSIONS
!=======================================================================
subroutine test_expansions()
    implicit none
    
    integer, parameter                              :: groups = 361
    integer, parameter                              :: meshes = 1
    integer, parameter                              :: angles = 1
    
    real(kind=8), dimension(groups,meshes,angles)   :: arr1, arr2
    type(DGM_structure)                             :: my_structure
    type(coarse_group), dimension(:), allocatable   :: my_cg
    
    integer                                         :: i,j,g
    ! real(kind=8)                                    :: temp        
    
    call srand(0)

    my_structure=get_DGM_structure(groups)
    allocate(my_cg(my_structure%coarse_groups))
    
    do g=1,my_structure%coarse_groups
        my_cg(g)%fine_groups = my_structure%cg_map(g)
        my_cg(g)%min_group = sum(my_structure%cg_map(1:g-1))+1
        allocate(my_cg(g)%aflux(0:my_cg(g)%fine_groups-1,meshes,angles))
    enddo
    
    
    open(unit=42, file='arrays', action='read')
    read(42,*)
    do g=1,groups
        read(42,*) arr1(g,1,1) ! , temp
    enddo
    close(unit=42)
    call random_number(arr1)
    
    
    call start_timer(1)
    
    allocate(DCT(maxval(my_structure%cg_map)))
    do g=1,my_structure%coarse_groups
        call compute_DCTs(my_structure%cg_map(g))
    enddo
     
    call flux_moments(arr1,my_cg)     
    call build_fine_flux(arr2,my_cg)
     
    write(*,*) 'Inf norm arr1-arr2: ', maxval(abs(arr1-arr2))

end subroutine test_expansions

end module DGM_module
