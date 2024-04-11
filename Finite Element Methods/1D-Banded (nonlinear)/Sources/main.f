    ! --------------------------------------------------
    ! Non-linear problem u"(x)-Ku(x)^2=0, u(0)=0, u(1)=1
    ! --------------------------------------------------
program main

    use kind
    use data
    use library
    use meshing
    use assembly
    use io
    use omp_lib 
    implicit none

    integer :: i
    real(kind=rk), dimension(:,:), allocatable:: outputData
    character(len=50):: fileName

    call initialize()

    if ( bft == 1 ) then
        call lin_basis_phi1d()
    else if ( bft == 2 ) then
        call quad_basis_phi1d()
    end if
    
    call generate_unifmesh()

    UG = 1.0_rk
    UG(nVar) = uBC(1)
    UG(nVar) = uBC(2)

    do while (err > tol)

        itr = itr + 1
    
        ! Initialize all global matrices
        RG = 0.0_rk
        JG = 0.0_rk

        call global_assembly()

        call FullGaussSolverp(JG, RG, nVar)

        UG = UG + RG

        err = 0.0_rk
        do i = 1, nVar
            err = err + RG(i)*RG(i)
        end do
        err = sqrt(err)
        
        print*, "Iteration- ", itr, "| Error = ", err

    end do

    ! allocate(outputData(nVar, 2))
    ! outputData(:,1) = xMesh
    ! outputData(:,2) = UG

    ! fileName = 'solution'
    ! call output_data(outputData, fileName)

contains

end program
