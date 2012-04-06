module stop_watch

real(kind=8), dimension(2,10) :: time

contains

subroutine start_timer(i)
    integer, intent(in) :: i
    call cpu_time(time(1,i))
end subroutine start_timer

real(kind=8) function read_timer(i)
    integer, intent(in) :: i

    call cpu_time(time(2,i))
    read_timer = time(2,i)-time(1,i)
    return 
end function read_timer

end module stop_watch
