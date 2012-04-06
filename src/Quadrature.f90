module Quadrature

type QuadGL
    real(kind=8), dimension(:), allocatable         :: point
    real(kind=8), dimension(:), allocatable         :: weight 
    integer                                         :: order   
end type QuadGL

contains

type(QuadGL) function get_QuadGL(N)
    implicit none
    integer, intent(in)         :: N

    allocate(get_QuadGL%point(N), get_QuadGL%weight(N))
    get_QuadGL%order = N
    
    if (N == 8) then
        get_QuadGL%point = (/-0.9602898564,-0.7966664774, &
            & -0.5255324099,-0.1834346424, 0.1834346424, &
            & 0.5255324099,0.7966664774,0.9602898564 /)
        get_QuadGL%weight = (/0.1012285363,0.2223810344, &
            & 0.3137066459,0.3626837834, 0.3626837834, &
            & 0.3137066459,0.2223810344,0.1012285363 /)
    else
        write(*,*) 'Quadrature not supported.'
        stop
    endif
end function get_QuadGL

end module Quadrature
