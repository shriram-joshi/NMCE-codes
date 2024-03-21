Module Data_

    use kind
    use omp_lib
    implicit none

    ! Declare Global variables that will be frequently used during simulation
    ! Examples are:

    real(kind=rk):: Pe
    real(kind=rk), allocatable, dimension(:,:):: ph1, dph1
    real(kind=rk), dimension(3):: xiarr
    real(kind=rk), dimension(3):: wei
    real(kind=rk), dimension(:), allocatable:: RG, rL, xMesh
    real(kind=rk), dimension(:,:), allocatable:: JG, jL
    real(kind=rk), dimension(2):: uBC, xSpan
    integer :: nE, bft, nMesh, nVar, nLVar, nLP
    integer, dimension(:,:), allocatable:: nop

contains
    ! You can make subroutines here

    subroutine init_problem()
        implicit none
        
        call param_input()

        allocate(nop(nLP, nE))

        call generate_nop()

        ! Range of independant variable
        xSpan = (/0.0_rk, 1.0_rk/)
        ! BC for dependant variable
        uBC = (/0.0_rk, 1.0_rk/)

        ! Initiallize variables
        allocate(jL(nLVar,nLVar))
        allocate(rL(nLVar))
        allocate(RG(nVar))
        ! We allocate nE+2 columns since the FullGaussSolverp function requires a matrix with dimesnions (j,j+1)
        ! Check the FullGaussSolverp.f file for details
        allocate(JG(nVar,nVar+1))
        allocate(xMesh(nVar))
        allocate(ph1(nLP,3))
        allocate(dph1(nLP,3))

        ! Initialize all matrices
        rL = 0.0_rk
        RG = 0.0_rk
        jL = 0.0_rk
        JG = 0.0_rk

    end subroutine init_problem

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
        implicit none

        nE=0
        do while (nE <= 0)
            ! Total number of ELEMENTS
            print*, "Number of elements (n > 0) - "
            read*, nE
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
        bft = 0
        do while ((bft <= 0) .or. (bft > 2))
            ! Total number of ELEMENTS
            print*, "Type - "
            print*, "1 for Linear Basis Functions"
            print*, "2 for Quadratric Basis Functions"
            read*, bft
            print*
        end do

        ! Set number of local points nLP based on type of basis function
        if (bft == 1) then
            nLP = 2;
        else if (bft == 2) then
            nLP = 3;
        end if

        nLVar = nLP;
        nMesh = (nLP-1)*(nE-1) + nLP
        nVar = nMesh ! Change this if there are multiple unkowns points

    end subroutine param_input 

    subroutine generate_nop()
        implicit none
        
        integer:: i, j

        do i = 1,nLP
            do j = 1,nE
                nop(i,j) = (nLP-1)*(j-1)+i
            end do
        end do

        print*, "nop ="
        do i = 1, nLP
            do j = 1, nE
                write(*, fmt='(I8)', advance='no') nop(i, j)
            end do
            print *  ! Move to the next line after printing each row
        end do
    end subroutine generate_nop

end module Data_
