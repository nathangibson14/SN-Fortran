module DGM_module
use DCT_type
use stop_watch
use nuclear_data
use mesh_type
use state_module
use Quadrature
use solvers

type DGM_structure
    integer                             :: fine_groups
    integer                             :: coarse_groups
    integer, dimension(:), allocatable  :: cg_map
end type DGM_structure

type coarse_group
    integer                                     :: fine_groups
    integer                                     :: min_group
    real(kind=8), dimension(:,:,:), allocatable :: aflux        ! (x, angle, moment)
    
    real(kind=8), dimension(:), allocatable     :: tot, nufis   ! (x)
    real(kind=8), dimension(:,:), allocatable   :: chi          ! (x, moment)
    real(kind=8), dimension(:,:,:), allocatable :: scat         ! (x, group, moment)
    real(kind=8), dimension(:,:,:), allocatable :: delta      ! (x, angle, moment)  
    
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
        open(unit=41,file='DGM_structure/NG2042',action='read')
        
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
! DATA_MOMENTS
!=======================================================================
subroutine data_moments(xs,mesh,phi,psi,cg)
    implicit none
    
    type(xsdata), dimension(:), intent(in)          :: xs
    type(mesh_1d), intent(in)                       :: mesh
    real(kind=8), dimension(:,:), intent(in)        :: phi
    real(kind=8), dimension(:,:,:), intent(in)      :: psi
    type(coarse_group), dimension(:), intent(inout) :: cg

    integer                 :: g, gg, K, L, j, xi, ai    
    integer                 :: coarse_groups, N, angles
    real(kind=8)            :: del, temp
    
    coarse_groups = size(cg)
    angles = size(psi,3)
            
    do g = 1,coarse_groups
        N = cg(g)%fine_groups
        
        cg(g)%tot = 0.0
        cg(g)%nufis = 0.0
        cg(g)%scat = 0.0
        cg(g)%chi = 0.0
        cg(g)%delta = 0.0
    
        do xi=1,mesh%NFM
            
            ! chi
            do K=0,N-1
                do j=0,N-1
                   cg(g)%chi(xi,j) = cg(g)%chi(xi,j) + &
                        & DCT(N)%d(j,K)*xs(mesh%matl(xi))%chi(K+cg(g)%min_group)
                enddo
            enddo
            
            ! tot0g
            do K=0,N-1
                cg(g)%tot(xi) = cg(g)%tot(xi) + & 
                    & xs(mesh%matl(xi))%tot(K+cg(g)%min_group)*phi(K+cg(g)%min_group,xi)
            enddo
            cg(g)%tot(xi) = cg(g)%tot(xi) / &
                & sum(phi(cg(g)%min_group:cg(g)%min_group+N-1,xi))
            
            ! delta
            do K=0,N-1
                del = xs(mesh%matl(xi))%tot(K+cg(g)%min_group) - cg(g)%tot(xi)
                do j=0,N-1
                    do ai=1,angles
                        cg(g)%delta(xi,ai,j) = cg(g)%delta(xi,ai,j) + &
                            & DCT(N)%d(j,K)*del*psi(K+cg(g)%min_group,xi,ai)
                    enddo
                enddo
            enddo
            do j=0,N-1
                do ai=1,angles
                    cg(g)%delta(xi,ai,j) = cg(g)%delta(xi,ai,j) / &
                        & sum(psi(cg(g)%min_group:cg(g)%min_group+N-1,xi,ai))
                enddo
            enddo
            
            ! nufis
            do L=0,N-1
                cg(g)%nufis(xi) = cg(g)%nufis(xi) + & 
                    & xs(mesh%matl(xi))%nufis(cg(g)%min_group+L)*phi(cg(g)%min_group+L,xi)
            enddo
            cg(g)%nufis(xi) = cg(g)%nufis(xi) / &
                & sum(phi(cg(g)%min_group:cg(g)%min_group+N-1,xi))
                
            ! scat
            do gg=1,coarse_groups
                do j=0,N-1
                    do K=0,N-1
                        temp = 0
                        do L=0,cg(gg)%fine_groups-1
                            temp = temp+phi(cg(gg)%min_group+L,xi)* &
                                & xs(mesh%matl(xi))%scat(cg(g)%min_group+K,cg(gg)%min_group+L)
                        enddo
                        cg(g)%scat(xi,gg,j) = cg(g)%scat(xi,gg,j) + temp*DCT(N)%d(j,K)
                    enddo
                    cg(g)%scat(xi,gg,j) = cg(g)%scat(xi,gg,j) / &
                        & sum(phi(cg(gg)%min_group:cg(gg)%min_group+cg(gg)%fine_groups-1,xi))
                enddo
            enddo
        
        enddo ! xi
        
    enddo ! g

end subroutine data_moments

!=======================================================================
! INITIALIZE_COARSE_GROUP
!=======================================================================
subroutine initialize_coarse_group(cg, fine_groups, min_group, &
    & NFM, angles, coarse_groups)
    implicit none
    
    type(coarse_group), intent(inout)   :: cg
    integer, intent(in)                 :: fine_groups, min_group, NFM, &
        &                                  coarse_groups, angles
    
    cg%fine_groups = fine_groups
    cg%min_group = min_group
    allocate(cg%aflux(NFM, angles, 0:fine_groups-1))
    allocate(cg%tot(NFM), cg%nufis(NFM), cg%chi(NFM,0:fine_groups-1), &
        & cg%scat(NFM,coarse_groups,0:fine_groups-1), &
        & cg%delta(NFM, angles, 0:fine_groups-1))

end subroutine initialize_coarse_group

!=======================================================================
! TEST_EXPANSIONS
!=======================================================================
subroutine test_expansions()
    implicit none
    
    integer, parameter                              :: groups = 361
    integer, parameter                              :: meshes = 20
    integer, parameter                              :: angles = 8
    
    real(kind=8), dimension(groups,meshes,angles)   :: arr1, arr2
    type(DGM_structure)                             :: my_structure
    type(coarse_group), dimension(:), allocatable   :: my_cg
    
    integer                                         :: i,j,g     
    
    call srand(0)

    my_structure=get_DGM_structure(groups)
    allocate(my_cg(my_structure%coarse_groups))
    
    do g=1,my_structure%coarse_groups
        my_cg(g)%fine_groups = my_structure%cg_map(g)
        my_cg(g)%min_group = sum(my_structure%cg_map(1:g-1))+1
        allocate(my_cg(g)%aflux(0:my_cg(g)%fine_groups-1,meshes,angles))
    enddo
    
    
!~     open(unit=42, file='arrays', action='read')
!~     read(42,*)
!~     do g=1,groups
!~         read(42,*) arr1(g,1,1)
!~     enddo
!~     close(unit=42)
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

!=======================================================================
! DGM_SN
!=======================================================================
subroutine dgm_sn(mesh,quad,bc,xs,state,structure)
    implicit none
    
    type(mesh_1d), intent(in)                   :: mesh
    type(xsdata), dimension(:), intent(in)      :: xs
    type(QuadGL), intent(in)                    :: quad
    real(kind=8), dimension(2), intent(in)      :: bc
    type(state_type), intent(inout)             :: state
    type(DGM_structure), intent(in)             :: structure
    
    type(state_type)                                :: exact, dgmstate
    type(coarse_group), dimension(:), allocatable   :: cg
    type(xsdata), dimension(:), allocatable         :: dgmxs
    type(mesh_1d)                                   :: dgmmesh    
    
    integer                                     :: g, i, ai
    
    allocate(cg(structure%coarse_groups))
    do g=1,structure%coarse_groups
        call initialize_coarse_group(cg(g), structure%cg_map(g), &
            & sum(structure%cg_map(1:g-1))+1, mesh%NFM, quad%order, &
            structure%coarse_groups)
    enddo
    
    allocate(dgmxs(mesh%NFM))
    do i=1,mesh%NFM
        allocate(dgmxs(i)%tot(structure%coarse_groups))
        allocate(dgmxs(i)%nufis(structure%coarse_groups))
        allocate(dgmxs(i)%scat(structure%coarse_groups,structure%coarse_groups))
        allocate(dgmxs(i)%chi(structure%coarse_groups))
        allocate(dgmxs(i)%delta(structure%coarse_groups,quad%order))
    enddo
    
    call initialize_state(exact,xs(1)%groups,mesh%NFM,quad%order)
    call initialize_state(dgmstate,structure%coarse_groups,mesh%NFM,quad%order)
    
    allocate(DCT(maxval(structure%cg_map)))
    do g=1,structure%coarse_groups
        call compute_DCTs(structure%cg_map(g))
    enddo
    
    ! guess flux
    ! get exact solution to use
    call power_iteration(mesh,quad,bc,xs,exact,.false.)
    
    write(*,*) 'STARTING DGM CALCULATION...',read_timer(1)    
    
    call data_moments(xs,mesh,exact%phi,exact%psi_mesh,cg)
    do i=1,mesh%NFM
        dgmxs(i)%groups = structure%coarse_groups
        do g=1,structure%coarse_groups
            dgmxs(i)%tot(g) = cg(g)%tot(i)
            dgmxs(i)%nufis(g) = cg(g)%nufis(i)
            dgmxs(i)%scat(g,:) = cg(g)%scat(i,:,0)
            dgmxs(i)%chi(g) = cg(g)%chi(i,0)
            dgmxs(i)%delta(g,:) = cg(g)%delta(i,:,0)
        enddo
    enddo
    
    dgmmesh = make_mesh_1d(mesh%x_edge,(/(1,i=1,mesh%NFM)/),(/(i,i=1,mesh%NFM)/))
    
!~     open(unit=41,file='dgmxs',action='write',status='replace')
!~     do i=1,mesh%NFM
!~     write(41,*) 'tot nufis chi delta1 delta2 delta3 delta4 delta5 delta6'
!~     do g=1,structure%coarse_groups
!~         write(41,101) dgmxs(i)%tot(g), dgmxs(i)%nufis(g), dgmxs(i)%chi(g), dgmxs(i)%delta(g,1), & 
!~             & dgmxs(i)%delta(g,2), dgmxs(i)%delta(g,3), dgmxs(i)%delta(g,4) , dgmxs(i)%delta(g,5), &
!~             & dgmxs(i)%delta(g,6), dgmxs(i)%delta(g,7), dgmxs(i)%delta(g,8)
!~         write(41,101) dgmxs(2)%scat(g,:)
!~     enddo
!~     enddo
!~     101 format (1X, 12E12.4)
!~     close(unit=41)
    
    
    call power_iteration(dgmmesh,quad,bc,dgmxs,dgmstate,.true.)
    write(*,*) ' DGM KEIG: ', dgmstate%keig
    
    
    state = exact

end subroutine dgm_sn

end module DGM_module
