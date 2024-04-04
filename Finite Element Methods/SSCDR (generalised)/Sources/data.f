Module Data_

    use kind
    use omp_lib
    implicit none

    ! Declare Global variables that will be frequently used during simulation
    ! Examples are:
    real(kind=rk) :: Pe, Da
    real(kind=rk), allocatable, dimension(:,:) :: ph, dph
    real(kind=rk), dimension(3) :: xiarr
    real(kind=rk), dimension(3) :: wei
    real(kind=rk), dimension(:), allocatable :: R, xMesh
    real(kind=rk), dimension(:,:), allocatable :: JG, jL
    real(kind=rk), dimension(2) :: uBC, xSpan
    integer :: n, nl, nVar

contains

    ! You can make subroutines here

    subroutine Init_problem()
        
        call param_input()
        ! BC for u
        uBC = (/0.0_rk, 1.0_rk/)
        ! Range of x
        xSpan = (/0.0_rk, 1.0_rk/)

        ! Initiallize variables
        allocate(jL(nl,nl))
        allocate(R(nVar))
        ! We allocate nVar+1 columns since the FullGaussSolverp function requires a matrix with dimesnions (j,j+1)
        ! Check the FullGaussSolverp.f file for details
        allocate(JG(nVar,nVar+1))
        allocate(xMesh(n+1))
        allocate(ph(nl,3))
        allocate(dph(nl,3))

        R = 0.0_rk
        jL = 0.0_rk
        JG = 0.0_rk

    end subroutine Init_problem

    subroutine Gauss_points()
        implicit none

        xiarr(1)=(1.0_rk-0.774596669241483_rk)*0.5_rk
        xiarr(2)=0.5_rk
        xiarr(3)=(1.0_rk+0.774596669241483_rk)*0.5_rk
        wei(1)=(0.555555555555556_rk)*0.5_rk
        wei(2)=(0.888888888888889_rk)*0.5_rk
        wei(3)=(0.555555555555556_rk)*0.5_rk

    end subroutine Gauss_points

    subroutine param_input()
        n=0
        do while (n <= 0)
            ! Total number of ELEMENTS
            print*, "Number of elements (n > 0) - "
            read*, n
            print*
        end do

        nVar = n+1

        Pe = -1
        do while (Pe < 0)
            ! Setup Pe number and Neumann boundary conditions
            print*, "Peclet Number (Pe >= 0) - "
            read*, Pe
            print*
        end do

        Da = -1
        do while (Da < 0)
            ! Setup Pe number and Neumann boundary conditions
            print*, "Damkohler Number (Da >= 0) - "
            read*, Da
            print*
        end do

        ! Size of the local matrix, depends on type of elements
        ! 2 for linear, 3 for quadratic, so on
        nl = 0
        do while (.not.((nl == 2) .or. (nl == 3)))
            ! Total number of ELEMENTS
            print*, "Type of linear basis function."
            print*, "2 for Linear Basis Functions"
            print*, "3 for Quadratric Basis Functions"
            read*, nl
            print*
        end do

    end subroutine param_input 

end module Data_
