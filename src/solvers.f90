module solvers
use Quadrature
use mesh_type
use nuclear_data
use stop_watch
use state_module

contains

!=======================================================================
! SN_1D
!=======================================================================
subroutine sn_1d(mesh,quad,bc,xs,Q_iso,state,delta_flag,Q_disc)
    implicit none

    type(mesh_1d), intent(in)                           :: mesh
    type(xsdata), dimension(:), intent(in)              :: xs
    type(QuadGL), intent(in)                            :: quad
    real(kind=8), dimension(2), intent(in)              :: bc
    real(kind=8), dimension(:,:,:), intent(in), optional:: Q_disc ! group, mesh, angle
    real(kind=8), dimension(:,:), intent(in)            :: Q_iso  ! group, mesh
    type(state_type), intent(inout)                     :: state
    logical, intent(in)                                 :: delta_flag
       
        
    integer                             :: gsiter, siter
    integer                             :: gsiter_max = 500
    integer                             :: siter_max  = 500
    
    real(kind=8)                        :: gs_tol = 1.e-8 ! in phi 1 norm
    real(kind=8)                        :: s_tol  = 1.e-10 ! in psi infinity norm
    
    integer                             :: g, g_min, i, angle   
    integer                             :: groups, NFM, angles
    
    real(kind=8), dimension(:,:,:), allocatable :: psi_old_s
    real(kind=8), dimension(:,:), allocatable   :: phi_old
    real(kind=8), dimension(:), allocatable     :: Q_s, Q
    real(kind=8)                                :: Q_d
    
    real(kind=8), dimension(:), allocatable     :: tot    
    real(kind=8), dimension(:,:), allocatable   :: delta
   
    groups = xs(1)%groups
    NFM = mesh%NFM
    angles = quad%order
    
    allocate(psi_old_s(groups,NFM+1,angles))
    allocate(phi_old(groups,NFM))
    allocate(Q(NFM),Q_s(NFM))

    allocate(tot(NFM))
    allocate(delta(NFM,angles))
    delta = 0.0
    
    state%psi_edge = 0.0
    state%psi_mesh = 0.0
    phi_old = 0.0
    g_min = 1

    gauss_seidel: do gsiter=1,gsiter_max
             
        if (gsiter==1) then
            g_min = 1
        else
            if (groups < 100) then
                g_min = 1
            elseif (groups < 1000) then
                g_min = groups-100
            else
                g_min = groups-1000
            endif
        endif
        
!~         if (gsiter==1) then
!~             write(*,*) '    GS iter ', gsiter
!~         else
!~             write(*,*) '    GS iter ', gsiter,  &
!~                 & maxval(abs(state%phi(g_min:groups,:)-phi_old(g_min:groups,:)))
!~         endif
        
        phi_old = state%phi
        
        group_sweep: do g=g_min,groups
            ! build scattering source
            Q_s = 0.0
            do i=1,NFM
                Q_s(i) = dot_product(xs(mesh%matl(i))%scat(g,:),state%phi(:,i)) - &
                    & xs(mesh%matl(i))%scat(g,g)*state%phi(g,i)
            enddo
                
            do i=1,NFM
                tot(i)=xs(mesh%matl(i))%tot(g)
                if (delta_flag) delta(i,:) = xs(mesh%matl(i))%delta(g,:)
            enddo    
                
            source_iteration: do siter=1,siter_max
                psi_old_s(g,:,:) = state%psi_edge(g,:,:)
                
                ! build total iso source
                do i=1,NFM
                    Q(i) = 0.5*(Q_s(i)+xs(mesh%matl(i))%scat(g,g)*state%phi(g,i) + Q_iso(g,i))
                enddo
                
                ! SWEEP                
                call sweep_sd(state%psi_edge(g,:,:), state%psi_mesh(g,:,:), bc, &
                    & mesh, tot, quad, Q, delta=delta)
                
                ! calculate scalar flux
                state%phi(g,:) = 0.0
                do angle=1,angles
                    state%phi(g,:) = state%phi(g,:)+quad%weight(angle)*state%psi_mesh(g,:,angle)
                enddo
                
                ! check convergence 
                if (maxval(abs(psi_old_s(g,:,:)-state%psi_edge(g,:,:))) < s_tol) then
                    ! write(*,*) '       Source converged at ', siter
                    exit source_iteration
                elseif (siter == siter_max) then
                    write(*,*) 'Source iteration not converged!'
                    ! stop
                endif
                
            enddo source_iteration
            
        enddo group_sweep
        
        ! check convergence
!~         if (sum(abs(phi(g_min:groups,:)-phi_old(g_min:groups,:))) &
!~             & /(groups-g_min+1+NFM) < gs_tol) then
        if (sum(abs(state%phi-phi_old))/(groups+NFM) < gs_tol) then
            ! write(*,*) '       GS converged at ', gsiter
            exit gauss_seidel
        elseif (gsiter == gsiter_max) then
            write(*,*) 'Gauss-Seidel iteration not converged!'
            ! stop
        endif
    enddo gauss_seidel
    
    
!~     open(unit=51,file='/home/nathan/Documents/MATLAB/DGM_tools/SN/fdata',status='replace',action='write')
!~     write(51,*) 'x_mesh phi Q_f'
!~     do i=1,NFM
!~         write(51,*) mesh%x_mesh(i), phi(1,i), Q_iso(1,i)
!~     enddo
!~     close(unit=51)
!~     stop    
    
end subroutine sn_1d

!=======================================================================
! POWER_ITERATION
!=======================================================================
subroutine power_iteration(mesh,quad,bc,xs,state,delta_flag)
    implicit none
    
    type(mesh_1d), intent(in)                   :: mesh
    type(xsdata), dimension(:), intent(in)      :: xs
    type(QuadGL), intent(in)                    :: quad
    real(kind=8), dimension(2), intent(in)      :: bc
    type(state_type), intent(inout)             :: state
    logical, intent(in)                         :: delta_flag
    
    
    real(kind=8)                                :: k_old
    real(kind=8)                                :: fission_rate
    real(kind=8), dimension(:,:), allocatable   :: phi_old
    real(kind=8), dimension(:,:), allocatable   :: Q_f
    
    integer                                     :: kiter
    integer                                     :: kiter_max = 100
    real(kind=8)                                :: ktol = 1.e-6
    
    integer                                     :: groups, NFM
    integer                                     :: i, g   
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - -
    
    groups = xs(1)%groups
    NFM = mesh%NFM
    
    allocate(phi_old(groups,NFM))
    allocate(Q_f(groups,NFM))
    
    ! initialize scalar flux
    ! normalize to total fission rate
    state%phi = 1.0
    fission_rate = 0
    do g=1,groups
        do i=1,NFM
            fission_rate = fission_rate+state%phi(g,i)*xs(mesh%matl(i))%nufis(g)
        enddo
    enddo
    
!~     open(unit=41,file='nufis_check',action='write',status='replace')
!~     write(41,*) 'matl phi nufis'
!~     do g=1,groups
!~         do i=1,NFM
!~             write(41,*) mesh%matl(i), state%phi(g,i), xs(mesh%matl(i))%nufis(g)
!~         enddo
!~     enddo
!~     write(41,*) fission_rate
!~     close(unit=41)
    
    state%phi = state%phi/fission_rate
    
    ! initial guess for keig
    state%keig = 1.0
    
    
    ! begin power iteration loop
    kiter_loop: do kiter=1,kiter_max
        write(*,*) 'Power iteration ', kiter, ': ', read_timer(1)
    
        ! store old flux, keig
        phi_old = state%phi
        k_old = state%keig
        
        ! build fission source
        do g=1,groups
            do i=1,NFM
                Q_f(g,i) = dot_product(xs(mesh%matl(i))%nufis,state%phi(:,i)) * &
                    & xs(mesh%matl(i))%chi(g)/state%keig
            enddo
        enddo
        
        ! CALL SN_1D TO GET NEW FLUX
        call sn_1d(mesh,quad,bc,xs,Q_f,state,delta_flag)
        
        ! update flux, keig
        fission_rate = 0
        do g=1,groups
            do i=1,NFM
                fission_rate = fission_rate+state%phi(g,i)*xs(mesh%matl(i))%nufis(g)
            enddo
        enddo
        state%phi = state%phi/fission_rate
        state%keig = state%keig*fission_rate
        
        write(*,*) '    keig = ', state%keig
        
        ! check for convergence
        if ( abs(state%keig-k_old) < ktol ) then
            write(*,*) 'Power iteration converged at ', kiter
            exit kiter_loop
        elseif (kiter == kiter_max) then
            write(*,*) 'Power iteration not converged!'
            ! stop
        endif        
    
    enddo kiter_loop
    

end subroutine power_iteration

!=======================================================================
! SWEEP_SD
!=======================================================================
subroutine sweep_sd(psi_edge,psi_mesh,bc,mesh,tot,quad,Q_iso,Q_disc,delta)
    implicit none
    
    real(kind=8), dimension(:,:), intent(inout)         :: psi_edge,psi_mesh
    real(kind=8), dimension(2), intent(in)              :: bc
    type(mesh_1d), intent(in)                           :: mesh
    real(kind=8), dimension(:), intent(in)              :: tot
    real(kind=8), dimension(:,:), intent(in), optional  :: delta
    type(QuadGL), intent(in)                            :: quad
    real(kind=8), dimension(:), intent(in)              :: Q_iso
    real(kind=8), dimension(:,:), intent(in), optional  :: Q_disc
    
    integer         :: angles, NFM
    integer         :: angle, i
    real(kind=8)    :: Q_d, del    
    
    angles = quad%order
    NFM = mesh%NFM        

    Q_d = 0
    del = 0
    
    ! sweep over negative angles
    do angle = 1,angles/2
        psi_edge(NFM+1,angle) = bc(2)*psi_edge(NFM+1,angles-angle+1)
        do i=NFM,1,-1
            if (present(Q_disc)) Q_d = Q_disc(i,angle)
            if (present(delta)) del = delta(i,angle)
            psi_edge(i,angle) = &
                & psi_edge(i+1,angle)*(quad%point(angle)/ & 
                & (quad%point(angle)-mesh%dx(i)*(tot(i)+del))) + &
                & (Q_iso(i)+Q_d)*(-mesh%dx(i)/(quad%point(angle)- &
                & mesh%dx(i)*(tot(i)+del)))
        enddo
        psi_mesh(:,angle) = psi_edge(1:NFM,angle)
    enddo
    
    ! sweep over positive angles
    do angle = angles/2+1,angles
        psi_edge(1,angle) = bc(1)*psi_edge(1,angles-angle+1)
        do i=2,NFM+1
            if (present(Q_disc)) Q_d = Q_disc(i-1,angle)
            if (present(delta)) del = delta(i-1,angle)
            psi_edge(i,angle) = &
                & psi_edge(i-1,angle)*(quad%point(angle)/ &
                & (quad%point(angle)+mesh%dx(i-1)*(tot(i-1)+del))) + &
                & (Q_iso(i-1)+Q_d)*(mesh%dx(i-1)/(quad%point(angle)+ &
                & mesh%dx(i-1)*(tot(i-1)+del)))
        enddo
        psi_mesh(:,angle) = psi_edge(2:NFM+1,angle)
    enddo
    
end subroutine sweep_sd


end module solvers
