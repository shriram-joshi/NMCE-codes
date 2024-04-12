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

    ! Book keeping variables
    ! nUnk = # of unkowns solving for in problem
    ! nE = # of elements
    ! bft = Integer that determines basis function type (bft)
    !       Linear Basis Function bft = 1
    !       Quadratic Basis Function bft = 2
    ! nLP = # of points per element. For linear BF nLP = 2, quadratic BF nLP = 3, etc
    ! nMesh = # of mesh nodes.  = (nLP-1)*(nE-1) + nLP
    ! nVar = # of variables in the problem. nVar = nUnk * nMesh
    ! nLVar = # of variables in local problem. nLVar = nUnk * nLP
    ! nop(local point #, element #) = Global mesh point number matrix
    ! nopp(global node #) = gives first global variable number for the argument global node number
    ! MDF(global node #) = number of unkowns for the argument global node number
    integer :: nUnk = 1, nE, bft, nLP, nMesh, nVar,  nLVar, itr = 0
    integer, dimension(:,:), allocatable:: nop
    integer, dimension(:), allocatable:: MDF, nopp
contains

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

    subroutine generate_MDF()
        implicit none
        integer :: i
        
        do i = 1, nMesh
            MDF(i) = 1
        end do

    end subroutine generate_MDF

    subroutine generate_nopp()
        implicit none
            integer :: i
            nopp(1)=1
            do i=2,nMesh
                nopp(i)=nopp(i-1) + MDF(i-1)
            end do
    end subroutine generate_nopp

    subroutine initialize()
        implicit none

        call gauss_points()

        call inputs()

        ! Set number of local points nLP based on type of basis function
        if (bft == 1) then
            nLP = 2
        else if (bft == 2) then
            nLP = 3
        end if

        nLVar = nUnk * nLP
        nMesh = (nLP-1)*(nE-1) + nLP
        nVar = nUnk * nMesh ! Change this if there are multiple unkowns points

        allocate(nop(nLP, nE))

        call generate_nop()

        allocate(nopp(nMesh))

        call generate_nopp()

        allocate(MDF(nMesh))

        call generate_MDF()

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

    end subroutine initialize

    subroutine inputs()
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
        
    end subroutine inputs

end module data
