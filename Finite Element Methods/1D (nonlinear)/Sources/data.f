module data

    use kind
    use io
    use omp_lib
    implicit none

    ! Declare Global variables that will be frequently used during simulation
    ! Examples are:
    
    ! Constant K from the problem
    real(kind=rk):: KK, err = 1.0_rk, tol = 1e-6
    ! Gauss quadrature points xiarr and weights wei
    real(kind=rk), dimension(3):: xiarr, wei
    
    ! The basis function phi and d(phi)/d(xi) matrices for different gauss points 
    real(kind=rk), allocatable, dimension(:,:):: ph, dph
    ! Global and local assemble vectors for the JJu = RR problem
    real(kind=rk), dimension(:), allocatable:: RG, rL, xMesh, UG
    ! Global and local asseblmy matrix for the JJ.u = RR problem
    real(kind=rk), dimension(:,:), allocatable:: JG, jL
    
    ! uBC - the boundary condition for the unkown U
    ! xSpan - range of domain of solution
    real(kind=rk), dimension(2):: uBC = (/0.0_rk, 1.0_rk/), xSpan = (/0.0_rk, 1.0_rk/)

    ! Count keeping variables
    ! nUnk = # of unkowns solving for in problem
    ! nE = # of elements
    ! bft = Integer that determines basis function type (bft)
    !       Linear Basis Function bft = 1
    !       Quadratic Basis Function bft = 2
    ! nLP = # of points per element. For linear BF nLP = 2, quadratic BF nLP = 3, etc
    ! nMesh = # of mesh nodes.  = (nLP-1)*(nE-1) + nLP
    ! nVar = # of variables in the problem. nVar = nUnk * nMesh
    ! nLVar = # of variables in local problem. nLVar = nUnk * nLP
    integer :: nUnk = 1, nE, bft, nLP, nMesh, nVar,  nLVar, itr = 0

    ! Global mesh point number matrix 
    integer, dimension(:,:), allocatable:: nop

contains
    ! You can make subroutines here

    subroutine init_problem()
        implicit none
        
        call gauss_points()

        call param_input()

        allocate(nop(nLP, nE))

        call generate_nop()

        ! Initiallize variables
        allocate(jL(nLVar,nLVar))
        allocate(rL(nLVar))
        allocate(RG(nVar))
        allocate(UG(nVar))
        ! We allocate nE+2 columns since the FullGaussSolverp function requires a matrix with dimesnions (j,j+1)
        ! Check the FullGaussSolverp.f file for details
        allocate(JG(nVar,nVar+1))
        allocate(xMesh(nVar))
        allocate(ph(nLP,3))
        allocate(dph(nLP,3))

        ! Initialize all matrices
        rL = 0.0_rk
        RG = 0.0_rk
        jL = 0.0_rk
        JG = 0.0_rk
        UG = 0.0_rk

    end subroutine init_problem

    subroutine gauss_points()
        implicit none

        xiarr(1)=(1.0_rk-0.774596669241483_rk)*0.5_rk
        xiarr(2)=0.5_rk
        xiarr(3)=(1.0_rk+0.774596669241483_rk)*0.5_rk
        wei(1)=(0.555555555555556_rk)*0.5_rk
        wei(2)=(0.888888888888889_rk)*0.5_rk
        wei(3)=(0.555555555555556_rk)*0.5_rk

    end subroutine gauss_points

    subroutine generate_nop()
        implicit none
        
        integer:: i, j

        do i = 1,nLP
            do j = 1,nE
                nop(i,j) = (nLP-1)*(j-1)+i
            end do
        end do

    end subroutine generate_nop

    subroutine param_input()
        implicit none

        nE=0
        do while (nE <= 0)
            ! Total number of ELEMENTS
            print*, "Number of elements (n > 0) - "
            read*, nE
            print*
        end do

        KK = -1
        do while (KK < 0)
            print*, "K (K >= 0) - "
            read*, KK
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
            nLP = 2
        else if (bft == 2) then
            nLP = 3
        end if

        nLVar = nUnk * nLP
        nMesh = (nLP-1)*(nE-1) + nLP
        nVar = nUnk * nMesh ! Change this if there are multiple unkowns points

    end subroutine param_input 

end module data
