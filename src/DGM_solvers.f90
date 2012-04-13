module DGM_solvers
use DCT_type
use stop_watch
use nuclear_data
use mesh_type
use state_module
use Quadrature
use DGM_module
use solvers

contains
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
    
    integer                                     :: g, i, ai, gg, j
    real(kind=8), dimension(:), allocatable     :: Q_iso
    real(kind=8), dimension(:,:), allocatable   :: Q_disc
    real(kind=8), dimension(:), allocatable     :: Q_f, tot
    real(kind=8)                                :: Q_s, bc_err, fission_rate
    real(kind=8), dimension(:,:), allocatable   :: psi_edge
    real(kind=8)                                :: k_old
    real(kind=8), dimension(:,:), allocatable   :: phi_old        
    
    integer                                     :: bciter, dgmiter, update
    integer, parameter                          :: bciter_max = 100
    integer, parameter                          :: dgmiter_max = 15
    integer, parameter                          :: update_max = 4               
    real(kind=8), parameter                     :: dgm_tol = 1e-5
    real(kind=8), parameter                     :: bc_tol = 1e-9
    
    
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
    
    allocate(Q_iso(mesh%NFM), Q_disc(mesh%NFM,quad%order), Q_f(mesh%NFM))
    allocate(tot(mesh%NFM))
    allocate(psi_edge(mesh%NFM+1,quad%order))
    allocate(phi_old(structure%fine_groups, mesh%NFM))
    
    ! guess flux
    ! get exact solution to use
!~     call power_iteration(mesh,quad,bc,xs,exact,.false.)
!~     state = exact
    exact%keig = 0.0
    exact%psi_mesh = 0.0

    state%psi_mesh = 1.0_8
    state%phi = 1.0_8
    state%keig = 1.0_8
    
    write(*,*) ''
    write(*,*) ''
    write(*,*) 'STARTING DGM CALCULATION...',read_timer(1)    
    
    dgmmesh = make_mesh_1d(mesh%x_edge,(/(1,i=1,mesh%NFM)/),(/(i,i=1,mesh%NFM)/))
    
    dgm_iteration: do dgmiter=1,dgmiter_max
        write(*,*) 'DGM ITERATION ', dgmiter, ': ', read_timer(1)
    
        k_old = state%keig
        phi_old = state%phi
        
        ! flux moments, data moments
        call flux_moments(state%psi_mesh,cg)
        call data_moments(xs,mesh,state%phi,state%psi_mesh,cg)
                
        ! 0th order solution
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
        call power_iteration(dgmmesh,quad,bc,dgmxs,dgmstate,.true.)
        state%keig = dgmstate%keig
        
        write(*,*) ' DGM KEIG: ', state%keig
        write(*,*) ' FG KEIG:  ', exact%keig
        
        
        
        do g=1,structure%coarse_groups
            do i=1,mesh%NFM
                cg(g)%aflux(i,:,0) = dgmstate%psi_mesh(g,i,:)
                cg(g)%flux(i) = dgmstate%phi(g,i)
            enddo
        enddo
        
        write(*,*) 'Starting higher order moments:', read_timer(1)
        
        ! higher order moments
        
        ! calculate fission source
        Q_f = 0.0
        do g=1,structure%coarse_groups
            do i=1,mesh%NFM
                Q_f(i) = Q_f(i) + cg(g)%nufis(i)*cg(g)%flux(i)/dgmstate%keig
            enddo
        enddo
        
        psi_edge = 1.0_8
        
        do g=1,structure%coarse_groups
        
        
            do j=1,cg(g)%fine_groups-1
                
                ! build RHS
                do i=1,mesh%NFM
                    Q_s = 0.0
                    do gg=1,structure%coarse_groups
                        Q_s = Q_s + cg(g)%scat(i,gg,j)*cg(gg)%flux(i)
                    enddo
                
                    Q_iso(i) = 0.5_8*(Q_f(i)*cg(g)%chi(i,j)+Q_s)
                    do ai=1,quad%order
                        Q_disc(i,ai) = -cg(g)%delta(i,ai,j)*cg(g)%aflux(i,ai,0)
                    enddo
                enddo
                
                
                
                ! BC iterations
                bc_iteration: do bciter=1,bciter_max
                
                    call sweep_sd(psi_edge,cg(g)%aflux(:,:,j),bc,mesh, &
                        & cg(g)%tot,quad,Q_iso,Q_disc=Q_disc)
                    
                    ! check convergence
                    bc_err = 0
                    do ai=1,quad%order/2
                        bc_err = bc_err + abs(psi_edge(mesh%NFM+1,ai) - &
                            & bc(2)*psi_edge(mesh%NFM+1,quad%order-ai+1))
                    enddo
                    do ai=quad%order/2+1,quad%order
                        bc_err = bc_err + abs(psi_edge(1,ai) - &
                            & bc(1)*psi_edge(1,quad%order-ai+1))
                    enddo
                    
                    
                    if (bc_err < bc_tol) then
                        exit bc_iteration
                    elseif (bciter == bciter_max) then
                        write(*,*) 'BC iteration not converged!', bc_err, g,j
                        stop
                    endif
                
                enddo bc_iteration
              
        
            enddo
        enddo
        
        call build_fine_flux(state%psi_mesh, cg)    
        
        ! calculate scalar flux
        do g=1,structure%fine_groups
            state%phi(g,:) = 0.0
            do ai=1,quad%order
                state%phi(g,:) = state%phi(g,:)+quad%weight(ai)*state%psi_mesh(g,:,ai)
            enddo
        enddo
!~ 
!~         fission_rate = 0
!~         do g=1,structure%fine_groups
!~             do i=1,mesh%NFM
!~                 fission_rate = fission_rate+state%phi(g,i)*xs(mesh%matl(i))%nufis(g)
!~             enddo
!~         enddo        
!~         state%phi=state%phi/fission_rate

        write(*,*) 'Starting flux updates:', read_timer(1)
        
        ! flux updates go here
        flux_updates: do update=1,update_max
            
            call Lei_update(state, mesh, xs, quad, structure, bc)
!~             call Rahnema_update(state, mesh, xs, quad, structure, bc, phi_old, cg)
!~             call kill_negative_fluxes(state, mesh, xs, quad, structure, bc)
        
        enddo flux_updates
        
        write(*,*) 'Infinity norm of psi_mesh error: ', maxval(abs(state%psi_mesh-exact%psi_mesh))
        write(*,*) 'Minimum psi_mesh: ', minval(state%psi_mesh)
        write(*,*) ''
        write(*,*) ''
        
        ! convergence check
        if (abs(state%keig-k_old) < dgm_tol) then
            write(*,*) 'DGM iteration converged at ', dgmiter
            exit dgm_iteration
        elseif (dgmiter==dgmiter_max) then
            write(*,*) 'DGM ITERATION NOT CONVERGED!'
            stop 
        endif
    
    enddo dgm_iteration
    
    
    
!~     open(unit=41,file='dgm_aflux',action='write',status='replace')
!~     open(unit=42,file='exact_aflux',action='write',status='replace')
!~     do i=1,mesh%NFM
!~         write(41,102) state%psi_mesh(50,i,:)
!~         write(42,102) exact%psi_mesh(50,i,:)
!~     enddo
!~     102 format (8E12.4)
!~     close(unit=41)
    

end subroutine dgm_sn

!=======================================================================
! NO_HIGH_ORDER
!=======================================================================
subroutine no_high_order(mesh,quad,bc,xs,state,structure)
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
    
    integer                                     :: g, i, ai, gg, j
    real(kind=8), dimension(:), allocatable     :: Q_iso
    real(kind=8), dimension(:,:), allocatable   :: Q_disc
    real(kind=8), dimension(:), allocatable     :: Q_f, tot
    real(kind=8)                                :: Q_s, bc_err, fission_rate
    real(kind=8), dimension(:,:), allocatable   :: psi_edge
    real(kind=8)                                :: k_old
    real(kind=8), dimension(:,:), allocatable   :: phi_old        
    
    integer                                     :: bciter, dgmiter, update
    integer, parameter                          :: bciter_max = 100
    integer, parameter                          :: dgmiter_max = 30
    integer, parameter                          :: update_max = 1                
    real(kind=8), parameter                     :: dgm_tol = 1e-5
    real(kind=8), parameter                     :: bc_tol = 1e-9
    
    
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
    
    allocate(Q_iso(mesh%NFM), Q_disc(mesh%NFM,quad%order), Q_f(mesh%NFM))
    allocate(tot(mesh%NFM))
    allocate(psi_edge(mesh%NFM+1,quad%order))
    allocate(phi_old(structure%fine_groups, mesh%NFM))
    
    ! guess flux
    ! get exact solution to use
    call power_iteration(mesh,quad,bc,xs,exact,.false.)
!~     state = exact

    state%psi_mesh = 1.0_8
    state%phi = 1.0_8
    state%keig = 1.0_8
    
    write(*,*) ''
    write(*,*) ''
    write(*,*) 'STARTING DGM CALCULATION...',read_timer(1)    
    
    dgmmesh = make_mesh_1d(mesh%x_edge,(/(1,i=1,mesh%NFM)/),(/(i,i=1,mesh%NFM)/))
    
    dgm_iteration: do dgmiter=1,dgmiter_max
        write(*,*) 'DGM ITERATION ', dgmiter, ': ', read_timer(1)
    
        k_old = state%keig
        phi_old = state%phi
        
        ! flux moments, data moments
        call flux_moments(state%psi_mesh,cg)
        call data_moments0(xs,mesh,state%phi,state%psi_mesh,cg)
                
        ! 0th order solution
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
        call power_iteration(dgmmesh,quad,bc,dgmxs,dgmstate,.true.)
        state%keig = dgmstate%keig
        
        write(*,*) ' DGM KEIG: ', state%keig
        write(*,*) ' FG KEIG:  ', exact%keig
        
        
        
        do g=1,structure%coarse_groups
            do i=1,mesh%NFM
                cg(g)%aflux(i,:,0) = dgmstate%psi_mesh(g,i,:)
                cg(g)%flux(i) = dgmstate%phi(g,i)
            enddo
        enddo
        
        ! flux updates go here
        flux_updates: do update=1,update_max
            
!~             call Lei_update(state, mesh, xs, quad, structure, bc)
            call Rahnema_update(state, mesh, xs, quad, structure, bc, phi_old, cg)
!~             call kill_negative_fluxes(state, mesh, xs, quad, structure, bc)
        
        enddo flux_updates
        
        write(*,*) 'Infinity norm of psi_mesh error: ', maxval(abs(state%psi_mesh-exact%psi_mesh))
        write(*,*) 'Minimum psi_mesh: ', minval(state%psi_mesh)
        write(*,*) ''
        write(*,*) ''
        
        ! convergence check
        if (abs(state%keig-k_old) < dgm_tol) then
            write(*,*) 'DGM iteration converged at ', dgmiter
            exit dgm_iteration
        elseif (dgmiter==dgmiter_max) then
            write(*,*) 'DGM ITERATION NOT CONVERGED!'
            stop 
        endif
    
    enddo dgm_iteration

end subroutine no_high_order

!=======================================================================
! LEI_UPDATE
!=======================================================================
subroutine Lei_update(state, mesh, xs, quad, structure, bc)
    implicit none
    type(state_type), intent(inout)             :: state
    type(mesh_1d), intent(in)                   :: mesh    
    type(xsdata), dimension(:), intent(in)      :: xs
    type(QuadGL), intent(in)                    :: quad
    type(DGM_structure), intent(in)             :: structure
    real(kind=8), dimension(2), intent(in)      :: bc    
    
    real(kind=8), dimension(mesh%NFM)               :: Q_f, Q_iso, tot
    real(kind=8), dimension(mesh%NFM+1,quad%order)  :: psi_edge    
    real(kind=8)                                    :: Q_s
    
    integer                             :: i, g, ai
    
    integer                             :: bciter
    integer, parameter                  :: bciter_max = 100
    real(kind=8), parameter             :: bc_tol = 1.e-9
    real(kind=8)                        :: bc_err
    
    Q_f = 0.0
    do g=1,structure%fine_groups
        do i=1,mesh%NFM
            Q_f(i) = Q_f(i) + xs(mesh%matl(i))%nufis(g)*state%phi(g,i)/state%keig
        enddo
    enddo

    psi_edge = 0.0

    do g=1,structure%fine_groups
        
        ! build tot
        do i=1,mesh%NFM
            tot(i) = xs(mesh%matl(i))%tot(g)
        enddo
        
        ! build RHS
        do i=1,mesh%NFM
            Q_s = dot_product(xs(mesh%matl(i))%scat(g,:),state%phi(:,i))
        
            Q_iso(i) = 0.5_8*(Q_f(i)*xs(mesh%matl(i))%chi(g)+Q_s)
        enddo
        
       
        
        ! BC iterations
        bc_iteration2: do bciter=1,bciter_max
        
            call sweep_sd(psi_edge,state%psi_mesh(g,:,:),bc,mesh, &
                & tot,quad,Q_iso)
            
            ! check convergence
            bc_err = 0
            do ai=1,quad%order/2
                bc_err = bc_err + abs(psi_edge(mesh%NFM+1,ai) - &
                    & bc(2)*psi_edge(mesh%NFM+1,quad%order-ai+1))
            enddo
            do ai=quad%order/2+1,quad%order
                bc_err = bc_err + abs(psi_edge(1,ai) - &
                    & bc(1)*psi_edge(1,quad%order-ai+1))
            enddo
            
            
            if (bc_err < bc_tol) then
                exit bc_iteration2
            elseif (bciter == bciter_max) then
                write(*,*) 'BC iteration in update not converged!', bc_err
                stop
            endif
        
        enddo bc_iteration2

        ! calculate scalar flux
        state%phi(g,:) = 0.0
        do ai=1,quad%order
            state%phi(g,:) = state%phi(g,:)+quad%weight(ai)*state%psi_mesh(g,:,ai)
        enddo

    enddo 


end subroutine Lei_update

!=======================================================================
! KILL_NEGATIVE_FLUXES
!=======================================================================
subroutine kill_negative_fluxes(state,mesh,xs,quad, structure, bc)
    implicit none
    
    type(state_type), intent(inout)             :: state
    type(mesh_1d), intent(in)                   :: mesh    
    type(xsdata), dimension(:), intent(in)      :: xs
    type(QuadGL), intent(in)                    :: quad
    type(DGM_structure), intent(in)             :: structure
    real(kind=8), dimension(2), intent(in)      :: bc   
    
    integer                                     :: i, g, ai
    real(kind=8)                                :: fission_rate
    
!~     if (minval(state%psi_mesh)<0.0) then
!~         state%psi_mesh = state%psi_mesh - 2._8*minval(state%psi_mesh)
!~     endif
    where (state%psi_mesh < 0)
        state%psi_mesh = -state%psi_mesh
    endwhere
    
    ! calculate scalar flux
    do g=1,structure%fine_groups
        state%phi(g,:) = 0.0
        do ai=1,quad%order
            state%phi(g,:) = state%phi(g,:)+quad%weight(ai)*state%psi_mesh(g,:,ai)
        enddo
    enddo
    
    fission_rate = 0
    do g=1,structure%fine_groups
        do i=1,mesh%NFM
            fission_rate = fission_rate+state%phi(g,i)*xs(mesh%matl(i))%nufis(g)
        enddo
    enddo        
    state%phi=state%phi/fission_rate

end subroutine kill_negative_fluxes

!=======================================================================
! RAHNEMA_UPDATE
!=======================================================================
subroutine Rahnema_update(state, mesh, xs, quad, structure, bc, phi_old, cg)
    implicit none
    
    type(state_type), intent(inout)             :: state
    type(mesh_1d), intent(in)                   :: mesh    
    type(xsdata), dimension(:), intent(in)      :: xs
    type(QuadGL), intent(in)                    :: quad
    type(DGM_structure), intent(in)             :: structure
    real(kind=8), dimension(2), intent(in)      :: bc    
    real(kind=8), dimension(:,:), intent(in)    :: phi_old
    type(coarse_group), dimension(:), intent(in):: cg        
    
    real(kind=8), dimension(mesh%NFM)               :: Q_f, Q_iso, tot
    real(kind=8), dimension(mesh%NFM+1,quad%order)  :: psi_edge    
    real(kind=8)                                    :: Q_s
    
    integer                             :: i, g, ai, gg, K, L, ggg
    
    integer                             :: bciter
    integer, parameter                  :: bciter_max = 100
    real(kind=8), parameter             :: bc_tol = 1.e-9
    real(kind=8)                        :: bc_err, temp
    
    real(kind=8), dimension(:,:,:), allocatable   :: scat_rec
    
    ! scat_rec
    allocate(scat_rec(structure%fine_groups, structure%coarse_groups, mesh%NFM))
!~     allocate(scat_rec(1, structure%coarse_groups, 1))
    do i=1,mesh%NFM
        do g=1,structure%fine_groups
            do gg=1,structure%coarse_groups
                temp = 0
                do L=0,cg(gg)%fine_groups-1
                    temp = temp+phi_old(cg(gg)%min_group+L,i)* &
                        & xs(mesh%matl(i))%scat(g,cg(gg)%min_group+L)
                enddo
                scat_rec(g,gg,i) = temp / sum(phi_old(cg(gg)%min_group: &
                    cg(gg)%min_group+cg(gg)%fine_groups-1,i))
            enddo
        enddo
    enddo
    
    
    
    Q_f = 0.0
    do g=1,structure%coarse_groups
        do i=1,mesh%NFM
            Q_f(i) = Q_f(i) + cg(g)%nufis(i)*cg(g)%flux(i)/state%keig
        enddo
    enddo

    psi_edge = 0.0

    do g=1,structure%fine_groups
        
        ! build tot
        do i=1,mesh%NFM
            tot(i) = xs(mesh%matl(i))%tot(g)
        enddo
        
        ! build RHS
        do i=1,mesh%NFM
            Q_s = 0
            do gg=1,structure%coarse_groups
!~                 do ggg=1,structure%coarse_groups
!~                     temp = 0
!~                     do L=0,cg(gg)%fine_groups-1
!~                         temp = temp+phi_old(cg(gg)%min_group+L,i)* &
!~                             & xs(mesh%matl(i))%scat(g,cg(gg)%min_group+L)
!~                     enddo
!~                     scat_rec(1,gg,1) = temp / sum(phi_old(cg(gg)%min_group: &
!~                         cg(gg)%min_group+cg(gg)%fine_groups-1,i))
!~                 enddo
                Q_s = Q_s +scat_rec(g,gg,i)*cg(gg)%flux(i)
!~                 Q_s = Q_s +scat_rec(1,gg,1)*cg(gg)%flux(i)
            enddo
        
            Q_iso(i) = 0.5_8*(Q_f(i)*xs(mesh%matl(i))%chi(g)+Q_s)
        enddo
        
       
        
        ! BC iterations
        bc_iteration2: do bciter=1,bciter_max
        
            call sweep_sd(psi_edge,state%psi_mesh(g,:,:),bc,mesh, &
                & tot,quad,Q_iso)
            
            ! check convergence
            bc_err = 0
            do ai=1,quad%order/2
                bc_err = bc_err + abs(psi_edge(mesh%NFM+1,ai) - &
                    & bc(2)*psi_edge(mesh%NFM+1,quad%order-ai+1))
            enddo
            do ai=quad%order/2+1,quad%order
                bc_err = bc_err + abs(psi_edge(1,ai) - &
                    & bc(1)*psi_edge(1,quad%order-ai+1))
            enddo
            
            
            if (bc_err < bc_tol) then
                exit bc_iteration2
            elseif (bciter == bciter_max) then
                write(*,*) 'BC iteration in update not converged!', bc_err
                stop
            endif
        
        enddo bc_iteration2

        ! calculate scalar flux
        state%phi(g,:) = 0.0
        do ai=1,quad%order
            state%phi(g,:) = state%phi(g,:)+quad%weight(ai)*state%psi_mesh(g,:,ai)
        enddo

    enddo 
    
end subroutine Rahnema_update

end module DGM_solvers
