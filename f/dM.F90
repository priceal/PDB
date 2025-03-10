
!! program to calculate the distance matrix from an input
!! list of coordinates. first line of input file must have number
!! of coordinates. must be 3 dimensional, each coordinate triple
!! on one line.
!! output will have number of lines equal to number of coordinates
!! each line with number of coordinates values (a square array)

Program dM
    implicit none
    
    ! declare and initialize default values
    integer :: i, j, numCoordinates=0, dims=3
    
    ! first index is dimension and second is coordinate number
    real(8), dimension(:,:), allocatable :: x
    
    ! square matrix of distances
    real(8), dimension(:,:), allocatable :: distanceMatrix
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! parameter input from standard input
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    read *, numCoordinates
    allocate ( x(dims,numCoordinates) )
    allocate ( distanceMatrix(numCoordinates,numCoordinates) )
    do i = 1, numCoordinates
        read *, x(1,i), x(2,i), x(3,i)
    end do    
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    ! main algorithm - enforces symmetry and zero diagonal
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, numCoordinates
        do j = i+1, numCoordinates
            distanceMatrix(i,j)=magnitude(x(:,i)-x(:,j))
            distanceMatrix(i,i)=0.0
            distanceMatrix(j,i)=distanceMatrix(i,j)
        end do
    end do
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    ! output matrix
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, numCoordinates
        print *, distanceMatrix(i,:)
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
! procedures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
contains
    function magnitude(xx) result (mag)
        implicit none

    	! declarations, only done first time function is called
        real(8), dimension(:) :: xx
        real(8) :: mag, norm2
        integer :: i
        
        ! distance calculation
        norm2 = 0.0
    	do i = 1, size(xx)
    	    norm2 = norm2 + xx(i)*xx(i)
    	end do
        mag = sqrt(norm2)
        
        return 
        
    end function magnitude

End Program dM


