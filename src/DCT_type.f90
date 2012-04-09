module DCT_type

real(kind=8), parameter :: pi = 3.141592653589793238462643383279_8

type discreteCosineTransform
    real(kind=8), dimension(:,:), allocatable :: d
end type discreteCosineTransform

type(discreteCosineTransform), dimension(:), allocatable :: DCT

contains

real(kind=8) function get_DCT(i,K,N)
    implicit none
    integer, intent(in) :: i,K,N
    
    if (i .GE. N .OR. K .GE. N .OR. i < 0 .OR. K < 0 .OR. N < 0) then
        write(*,*) 'Bad inputs to get_DCT'
        stop
    endif
        
    get_DCT =  dcos((real(i,8)*pi/real(N,8))*(real(K,8)+.5_8))
!~     write(*,*) (real(i,8)*pi/real(N,8))*(real(K,8)+.5_8)

end function get_DCT

subroutine compute_DCTs(N)
    implicit none
    
    integer, intent(in)         :: N
    integer                     :: i, K
    
    if (.NOT. allocated(DCT(N)%d)) then
        allocate(DCT(N)%d(0:N-1,0:N-1))
        do K=0,N-1
            do i=0,N-1
                DCT(N)%d(i,K) = get_DCT(i,K,N)
            enddo
        enddo
    endif
    
end subroutine compute_DCTs

end module DCT_type
