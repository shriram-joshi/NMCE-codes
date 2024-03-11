Module Data_

    use kind
    use omp_lib
    implicit none

    ! Declare Global variables that will be frequently used during simulation
    ! Examples are:
    real(kind=rk) :: C1, C2
    real(kind=rk) :: Oh, Gbond
    real(kind=rk) :: ph1(3,6),dph1(3,6)
    real(kind=rk), dimension(6) :: xiarr
    real(kind=rk), dimension(3) :: wei

    real(kind=rk) :: Pe
    real(kind=rk), dimension(:), allocatable :: R, xMesh
    real(kind=rk), dimension(:,:), allocatable :: J, jl
    real(kind=rk), dimension(2) :: ubc, xSpan

    integer :: n, nl, c

contains

    ! You can make subroutines here

    subroutine Init_problem()
        n=0
        do while (n <= 0)
            ! Total number of ELEMENTS
            print*, "Number of elements (n > 0) - "
            read*, n
            print*
        end do

        Pe = -1
        do while (Pe < 0)
            ! Setup Pe number and Neumann boundary conditions
            print*, "Peclet Number (Pe >= 0) - "
            read*, Pe
            print*
        end do

        ! Size of the local matrix, depends on type of elements
        ! 2 for linear, 3 for quadratic, so on
        nl = 2 

        ! BC for u
        ubc = (/0.0_rk, 1.0_rk/)
        ! Range of x
        xSpan = (/0.0_rk, 1.0_rk/)

        ! Initiallize variables
        allocate(jl(nl,nl))
        allocate(R(n+1))
        allocate(xMesh(n+1))
        ! We allocate n+2 columns since the FullGaussSolverp function requires a matrix with dimesnions (n,n+1)
        ! Check the FullGaussSolverp.f file for details
        allocate(J(n+1,n+2)) 
        
        R(1:n+1) = 0.0_rk
        jl(1:nl,1:nl) = 0.0_rk
        J(1:n+1,1:n+2) = 0.0_rk

    end subroutine Init_problem

    subroutine Gauss_points()
        implicit none

        xiarr(1)=(1.0_rk-0.774596669241483_rk)*0.5_rk
        xiarr(2)=0.5_rk
        xiarr(3)=(1.0_rk+0.774596669241483_rk)*0.5_rk
        xiarr(4)=0.0_rk
        xiarr(5)=0.5_rk
        xiarr(6)=1.0_rk
        wei(1)=(0.555555555555556_rk)*0.5_rk
        wei(2)=(0.888888888888889_rk)*0.5_rk
        wei(3)=(0.555555555555556_rk)*0.5_rk

    end subroutine Gauss_points
end module Data_
