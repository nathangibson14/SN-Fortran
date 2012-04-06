program SN_Fortran
use DCT_type
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

call SHEM_pin_cell(sigma)
! sigma(1) = get_nuclear_data(3,'test_3g')
mesh = make_mesh_1d((/0._8,2._8,3._8/),(/10,10/),(/1,2/))
quad = get_QuadGL(8)
bc = 1.0

call start_timer(1)
write(*,*) 'Starting power iteration'
call power_iteration(mesh,quad,bc,sigma,keig)

write(*,*) '    keig = ', keig
write(*,*) 'Program executed in ', read_timer(1), ' seconds'

end program SN_Fortran
