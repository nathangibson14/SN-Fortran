module mesh_type

type mesh_1d
    real(kind=8), dimension(:), allocatable     :: x_edge
    real(kind=8), dimension(:), allocatable     :: x_mesh
    real(kind=8), dimension(:), allocatable     :: dx
    integer, dimension(:), allocatable          :: matl 
    integer                                     :: NFM  
end type mesh_1d


contains

type(mesh_1d) function make_mesh_1d(edges, NFM, matl)
    implicit none
    
    real(kind=8), dimension(:), intent(in)      :: edges
    integer, dimension(:), intent(in)           :: NFM, matl
    integer                                     :: i
    
    make_mesh_1d%NFM = sum(NFM)
    
    allocate(make_mesh_1d%x_edge(sum(NFM)+1), &
        make_mesh_1d%x_mesh(sum(NFM)), make_mesh_1d%dx(sum(NFM)), &
        make_mesh_1d%matl(sum(NFM)))
    
    do i=1,size(NFM)
        if (i==1) then
            make_mesh_1d%dx(1:NFM(1)) = (edges(2)-edges(1))/NFM(1)
            make_mesh_1d%matl(1:NFM(1)) = matl(i)
        else
            make_mesh_1d%dx(sum(NFM(1:i-1))+1:sum(NFM(1:i))) = &
                & (edges(i+1)-edges(i))/NFM(i)
            make_mesh_1d%matl(sum(NFM(1:i-1))+1:sum(NFM(1:i))) = matl(i)
        endif
    enddo
    
    make_mesh_1d%x_edge(1) = edges(1)
    do i=1,sum(NFM)
        make_mesh_1d%x_edge(i+1) = make_mesh_1d%x_edge(i) + make_mesh_1d%dx(i)
        make_mesh_1d%x_mesh(i) = 0.5*(make_mesh_1d%x_edge(i)+make_mesh_1d%x_edge(i+1))
    enddo

end function make_mesh_1d


end module mesh_type
