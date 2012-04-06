module DCT_type

real(kind=8), parameter :: pi = 3.14159265358979323846

contains

real(kind=8) function get_DCT(i,K,N)
    implicit none
    integer, intent(in) :: i,K,N
    
    if (i .GE. N .OR. K .GE. N .OR. i < 0 .OR. K < 0 .OR. N < 0) then
        write(*,*) 'Bad inputs to get_DCT'
        stop
    endif
        
    get_DCT =  cos((i*pi/N)*(K+.5))

end function get_DCT

end module DCT_type
