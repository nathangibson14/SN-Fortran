module state_module

type state_type
    real(kind=8), dimension(:,:), allocatable   :: phi
    real(kind=8), dimension(:,:,:), allocatable :: psi_mesh, psi_edge ! group, mesh, angle
    real(kind=8)                                :: keig
end type state_type

contains

subroutine initialize_state(this,groups,NFM,angles)
    implicit none
    
    type(state_type), intent(inout)     :: this
    integer, intent(in)                 :: groups, NFM, angles
    
    allocate(this%phi(groups,NFM), this%psi_mesh(groups,NFM,angles), &
        & this%psi_edge(groups,NFM+1,angles))

end subroutine initialize_state

end module state_module
