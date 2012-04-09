module solvers
use Quadrature
use mesh_type
use nuclear_data
use stop_watch

contains

!=======================================================================
! SN_1D
!=======================================================================
subroutine sn_1d(mesh,quad,bc,xs,Q_iso,phi,Q_disc)
    implicit none

    type(mesh_1d), intent(in)                           :: mesh
    type(xsdata), dimension(:), intent(in)              :: xs
    type(QuadGL), intent(in)                            :: quad
    real(kind=8), dimension(2), intent(in)              :: bc
    real(kind=8), dimension(:,:,:), intent(in), optional:: Q_disc ! group, mesh, angle
    real(kind=8), dimension(:,:), intent(in)            :: Q_iso  ! group, mesh
    real(kind=8), dimension(:,:), intent(inout)         :: phi
       
        
    integer                             :: gsiter, siter
    integer                             :: gsiter_max = 500
    integer                             :: siter_max  = 500
    
    real(kind=8)                        :: gs_tol = 1.e-8 ! in phi 1 norm
    real(kind=8)                        :: s_tol  = 1.e-10 ! in psi infinity norm
    
    integer                             :: g, g_min, i, angle   
    integer                             :: groups, NFM, angles
    
    real(kind=8), dimension(:,:,:), allocatable :: psi_edge, psi_mesh, psi_old_s
    real(kind=8), dimension(:,:), allocatable   :: phi_old
    real(kind=8), dimension(:), allocatable     :: Q_s, Q
    
    real(kind=8), dimension(10)                 :: t
    
    t = 0.0
    
    groups = xs(1)%groups
    NFM = mesh%NFM
    angles = quad%order
    
    allocate(psi_edge(groups,NFM+1,angles), psi_mesh(groups,NFM,angles), &
        & psi_old_s(groups,NFM+1,angles))
    allocate(phi_old(groups,NFM))
    allocate(Q(NFM),Q_s(NFM))
    
    psi_edge = 0.0
    psi_mesh = 0.0
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
        
        write(*,*) '    GS iter ', gsiter, maxval(abs(phi(g_min:groups,:)-phi_old(g_min:groups,:)))
        phi_old = phi
        
        group_sweep: do g=g_min,groups
            ! build scattering source
            Q_s = 0.0
            do i=1,NFM
                Q_s(i) = dot_product(xs(mesh%matl(i))%scat(g,:),phi(:,i)) - &
                    & xs(mesh%matl(i))%scat(g,g)*phi(g,i)
            enddo
                
            source_iteration: do siter=1,siter_max
                psi_old_s(g,:,:) = psi_edge(g,:,:)
                
                ! build total iso source
                do i=1,NFM
                    Q(i) = 0.5*(Q_s(i)+xs(mesh%matl(i))%scat(g,g)*phi(g,i) + Q_iso(g,i))
                enddo
                
                ! sweep over negative angles
                do angle = 1,angles/2
                    psi_edge(g,NFM+1,angle) = bc(2)*psi_edge(g,NFM+1,angles-angle+1)
                    do i=NFM,1,-1
                        psi_edge(g,i,angle) = &
                            & psi_edge(g,i+1,angle)*(quad%point(angle)/ & 
                            & (quad%point(angle)-mesh%dx(i)*xs(mesh%matl(i))%tot(g))) + &
                            & Q(i)*(-mesh%dx(i)/(quad%point(angle)-mesh%dx(i)*xs(mesh%matl(i))%tot(g)))
                    enddo
                    psi_mesh(g,:,angle) = psi_edge(g,1:NFM,angle)
                enddo
                
                ! sweep over positive angles
                do angle = angles/2+1,angles
                    psi_edge(g,1,angle) = bc(1)*psi_edge(g,1,angles-angle+1)
                    do i=2,NFM+1
                        psi_edge(g,i,angle) = &
                            & psi_edge(g,i-1,angle)*(quad%point(angle)/ &
                            & (quad%point(angle)+mesh%dx(i-1)*xs(mesh%matl(i-1))%tot(g))) + &
                            & Q(i-1)*(mesh%dx(i-1)/(quad%point(angle)+mesh%dx(i-1)*xs(mesh%matl(i-1))%tot(g)))
                    enddo
                    psi_mesh(g,:,angle) = psi_edge(g,2:NFM+1,angle)
                enddo
                
                ! calculate scalar flux
                phi(g,:) = 0.0
                do angle=1,angles
                    phi(g,:) = phi(g,:)+quad%weight(angle)*psi_mesh(g,:,angle)
                enddo
                
                ! check convergence 
                if (maxval(abs(psi_old_s(g,:,:)-psi_edge(g,:,:))) < s_tol) then
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
        if (sum(abs(phi-phi_old))/(groups+NFM) < gs_tol) then
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
subroutine power_iteration(mesh,quad,bc,xs,keig)
    implicit none
    
    type(mesh_1d), intent(in)                   :: mesh
    type(xsdata), dimension(:), intent(in)      :: xs
    type(QuadGL), intent(in)                    :: quad
    real(kind=8), dimension(2), intent(in)      :: bc
    real(kind=8), intent(out)                   :: keig
    
    
    real(kind=8)                                :: k_old
    real(kind=8)                                :: fission_rate
    real(kind=8), dimension(:,:), allocatable   :: phi, phi_old
    real(kind=8), dimension(:,:), allocatable   :: Q_f
    
    integer                                     :: kiter
    integer                                     :: kiter_max = 100
    real(kind=8)                                :: ktol = 1.e-6
    
    integer                                     :: groups, NFM
    integer                                     :: i, g   
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - -
    
    groups = xs(1)%groups
    NFM = mesh%NFM
    
    allocate(phi(groups,NFM), phi_old(groups,NFM))
    allocate(Q_f(groups,NFM))
    
    
    
    ! initialize scalar flux
    ! normalize to total fission rate
    phi = 1.0
    fission_rate = 0
    do g=1,groups
        do i=1,NFM
            fission_rate = fission_rate+phi(g,i)*xs(mesh%matl(i))%nufis(g)
        enddo
    enddo
    phi = phi/fission_rate
    
    ! initial guess for keig
    keig = 1.0
    
    
    ! begin power iteration loop
    kiter_loop: do kiter=1,kiter_max
        write(*,*) 'Power iteration ', kiter, ': ', read_timer(1)
    
        ! store old flux, keig
        phi_old = phi
        k_old = keig
        
        ! build fission source
        do g=1,groups
            do i=1,NFM
                Q_f(g,i) = dot_product(xs(mesh%matl(i))%nufis,phi(:,i)) * &
                    & xs(mesh%matl(i))%chi(g)/keig
            enddo
        enddo
        
        ! CALL SN_1D TO GET NEW FLUX
        call sn_1d(mesh,quad,bc,xs,Q_f,phi)
        
        ! update flux, keig
        fission_rate = 0
        do g=1,groups
            do i=1,NFM
                fission_rate = fission_rate+phi(g,i)*xs(mesh%matl(i))%nufis(g)
            enddo
        enddo
        phi = phi/fission_rate
        keig = keig*fission_rate
        
        write(*,*) '    keig = ', keig
        
        ! check for convergence
        if ( abs(keig-k_old) < ktol ) then
            write(*,*) 'Power iteration converged at ', kiter
            exit kiter_loop
        elseif (kiter == kiter_max) then
            write(*,*) 'Power iteration not converged!'
            ! stop
        endif        
    
    enddo kiter_loop
    

end subroutine power_iteration


end module solvers
