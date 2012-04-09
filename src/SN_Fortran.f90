program SN_Fortran
use DGM_module
use nuclear_data
use solvers
use mesh_type
use stop_watch


implicit none

type(xsdata), dimension(2)  :: sigma
type(mesh_1d)               :: mesh
type(QuadGL)                :: quad
real(kind=8)                :: keig    
real(kind=8), dimension(2)  :: bc

!~ call test_pin_cell(sigma,2042,'NG2042')
!~ sigma(1) = get_nuclear_data(3,'test_3g')

!~ sigma(1) = get_nuclear_data(14767,'UF14767/u235')
!~ sigma(2) = get_nuclear_data(14767,'UF14767/h-h2o')

!~ mesh = make_mesh_1d((/0._8,2._8,3._8/),(/5,5/),(/1,2/))
!~ quad = get_QuadGL(8)
!~ bc = 1.0

!~ call start_timer(1)
!~ write(*,*) 'Starting power iteration'
!~ call power_iteration(mesh,quad,bc,sigma,keig)
!~ 
!~ write(*,*) '    keig = ', keig
!~ write(*,*) 'Program executed in ', read_timer(1), ' seconds'

call start_timer(1)
call test_expansions()
write(*,*) 'Program executed in ', read_timer(1), ' seconds'

end program SN_Fortran
