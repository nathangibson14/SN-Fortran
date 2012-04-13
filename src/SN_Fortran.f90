program SN_Fortran
use DGM_module
use DGM_solvers
use nuclear_data
use solvers
use mesh_type
use stop_watch
use state_module

implicit none

type(xsdata), dimension(2)  :: sigma
type(mesh_1d)               :: mesh
type(QuadGL)                :: quad  
real(kind=8), dimension(2)  :: bc
type(state_type)            :: state
type(DGM_structure)         :: structure 

call test_pin_cell(sigma,2042,'NG2042')
!~ sigma(1) = get_nuclear_data(3,'test_3g')

!~ sigma(1) = get_nuclear_data(14767,'UF14767/u235')
!~ sigma(2) = get_nuclear_data(14767,'UF14767/h-h2o')

mesh = make_mesh_1d((/0._8,0.393_8,0.631_8/),(/30,10/),(/1,2/))
quad = get_QuadGL(8)
bc = 1.0
call initialize_state(state,sigma(1)%groups,mesh%NFM,quad%order)
structure = get_DGM_structure(sigma(1)%groups)

call start_timer(1)

write(*,*) ''
write(*,*) ''
write(*,*) 'Starting power iteration'
!~ call power_iteration(mesh,quad,bc,sigma,state,.false.)
!~ call no_high_order(mesh,quad,bc,sigma,state,structure)
call dgm_sn(mesh,quad,bc,sigma,state,structure)

!~ write(*,*) '    keig = ', state%keig
write(*,*) 'Program executed in ', read_timer(1), ' seconds'

!~ call start_timer(1)
!~ call test_expansions()
!~ write(*,*) 'Program executed in ', read_timer(1), ' seconds'

end program SN_Fortran
